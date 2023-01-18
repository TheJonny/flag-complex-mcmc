use std::fs::File;
use std::io::prelude::*;

use flag_complex::prelude::*;

use crate::{MCMCSampler, Parameters, Precomputed};
use rand::Rng;
use serde::Serialize;
use serde::de::DeserializeOwned;

use hdf5;
use ndarray;

use crate::Edge;

// Flag File
pub fn read_flag_file<G: DirectedGraphNew>(fname:&str) -> G {
    let mut file = File::open(fname).expect("could not find .flag input file");
    let mut fcontents = String::new();
    file.read_to_string(&mut fcontents).unwrap();

    let mut lines = fcontents.lines();
    lines.next(); // skip dim 0
    let nnodes = lines.next().unwrap().split(' ').filter(|s| *s != "").count(); // todo: mehr beleidigungen // todo: better parsing
    let mut graph = G::new_disconnected(nnodes);
    lines.next(); // skip dim 1
    for line in lines {
        let mut ijw = line.split(' ').filter(|s| *s != "");
        if let (Some(i), Some(j)) = (ijw.next(), ijw.next()) {
            graph.add_edge(i.parse().unwrap(), j.parse().unwrap());
        }
    }
    return graph;
}

pub fn save_flag_file<G: DirectedGraph>(fname:&str, graph:&G) -> std::io::Result<()> {
    let mut content = "dim 0:\n".to_string();
    content += &("1 ".repeat(graph.nnodes()).trim_end().to_owned() + "\n"); // add vertices
    content += "dim 1:\n";
    let mut edges = graph.edges();
    edges.sort_unstable();
    for [i,j] in edges {
        content += &(format!("{} {} 1\n", i, j));
    }
    std::fs::write(fname, content).expect("Unable to write file");
    return Ok(());
}


// Save/Restore state (serde/bincode)
pub fn save_state<R: Rng+Serialize>(fname:&str, sampler:&MCMCSampler<R>) -> anyhow::Result<()> {
    let mut f = std::io::BufWriter::new(File::create(format!("{fname}.tmp")).unwrap());
    sampler.save(&mut f)?;
    f.flush()?;
    std::fs::rename(format!("{fname}.tmp"), fname).expect("moving temp state file to correct location failed");
    return Ok(());
}

pub fn load_state<'pa, 'pc, R: Rng+DeserializeOwned>(fname:&str, parameters: &'pa Parameters, precomputed: &'pc Precomputed,) -> anyhow::Result<MCMCSampler<'pa, 'pc, R>> {
    let mut f = std::io::BufReader::new(File::open(fname).unwrap());
    let sampler = MCMCSampler::load(&mut f, parameters, precomputed)?;
    return Ok(sampler);
}


// HDF5
pub fn save_to_hdf<G: DirectedGraph>(state_store_dir:&str, label:&str, seed:u64, sample_number:usize, graph:&G, flag_count:&Vec<usize>) -> hdf5::Result<()> {
    let file = hdf5::File::open_rw(format!("{state_store_dir}/{label}-{seed:03}.hdf5"))?;
    let groupname = format!("/{seed:03}/{sample_number:06}");
    if file.link_exists(&groupname) {
        file.unlink(&groupname)?;
    }
    let group = file.create_group(&groupname)?;
    let builder = group.new_dataset_builder();
    let mut edges = graph.edges();
    edges.sort_unstable();
    let ds = builder.deflate(4).with_data(&ndarray::arr2(&edges)).create("edgelist")?;
    ds.new_attr::<usize>().shape(flag_count.len()).create("flag_count")?.write(flag_count)?;
    ds.new_attr::<usize>().shape(1).create("number_of_vertices")?.write(&[graph.nnodes()])?;
    Ok(())
}

pub fn new_hdf_file(state_store_dir:&str, label:&str, seed:u64) -> hdf5::Result<()> {
    hdf5::File::create(format!("{state_store_dir}/{label}-{seed:03}.hdf5"))?;
    //TODO: Add metadata/comment to this file
    return Ok(());
}

pub fn save_dot<G: DirectedGraph, W: std::io::Write>(writer: &mut W, graph:&G) -> std::io::Result<()> {
    writeln!(writer, "digraph x {{")?;
    for [a,b] in graph.edges() {
        writeln!(writer, "{} -> {};", a, b)?;
    }
    writeln!(writer, "}}")?;
    Ok(())
}

pub fn load_graph_hdf5<G: DirectedGraphNew>(filename: &str, seed: u64, sample_number: u64) -> Result<G,Box<dyn std::error::Error>> {
    let groupname = format!("/{seed:03}/{sample_number:06}/edgelist");
    let h = hdf5::File::open(filename)?;
    let ds = h.dataset(&groupname)?;
    let arr = ds.read_2d()?;
    let nedges = arr.nrows();
    assert!(arr.ncols() == 2);
    let nnodes = ds.attr("number_of_vertices")?.read_1d()?[0];
    let mut g = G::new_disconnected(nnodes);
    for i in 0 .. nedges {
        let a = arr[(i,0)];
        let b = arr[(i,1)];
        g.add_edge(a, b);
    }
    return Ok(g);
}

pub fn load_edgelist_hdf5(filename: &str, groupname: &str) -> Result<Vec<Edge>, Box<dyn std::error::Error>> {
    let h = hdf5::File::open(filename)?;
    let arr = h.dataset(groupname)?.read_2d()?;
    let n = arr.nrows();
    assert!(arr.ncols() == 2);
    let mut res = Vec::with_capacity(n);
    for i in 0..n {
        let a = arr[(i,0)];
        let b = arr[(i,1)];
        res.push([a, b]);
    }
    Ok(res)
}

pub struct BitOutput {
    edges: Vec<Edge>,
    chunk_size: usize,
    index_in_file: usize,
    index_in_dir: usize,
    current_file: Option<std::io::BufWriter<std::fs::File>>,
    dir: String,
}

impl BitOutput {
    pub fn new<G: DirectedGraph>(graph: &G, dir: &str) -> std::io::Result<Self> {
        let res = std::fs::create_dir(dir);
        match res {
            Ok(_) => {
                // ok
            }
            Err(e) if e.kind() == std::io::ErrorKind::AlreadyExists => {
                // ok
            }
            Err(e) => {
                return Err(e);
            }
        }
        save_flag_file(&format!("{}/{}", dir, "graph.flag"), graph)?;
        let mut edges = graph.edges();
        for i in 0..edges.len() {
            let [a,b] = edges[i];
            edges.push([b,a]);
        }
        use std::cmp::{min, max};
        edges.sort_unstable_by_key(|&[a,b]| (max(a,b), min(a,b), a < b));
        edges.dedup();
        Ok(BitOutput {
            chunk_size: max(2_000_000_000 / (edges.len()/8), 1),
            edges,
            index_in_dir: 0,
            index_in_file: 0,
            current_file: None,
            dir: dir.to_owned(),
        })
    }
    pub fn save<G: DirectedGraph>(&mut self, graph: &G) -> std::io::Result<()> {
        // open file
        if self.index_in_file == 0 {
            assert!(self.current_file.is_none());
            let f = std::fs::File::create(&format!("{}/{}.edgebits", self.dir, self.index_in_dir))?;
            let b = std::io::BufWriter::new(f);
            self.current_file = Some(b);
        }
        // write graph
        {
            let f = self.current_file.as_mut().unwrap();
            let mut bit: u8 = 1;
            let mut byte = 0;
            for &[a,b] in &self.edges {
                if graph.has_edge(a, b) {
                    byte |= bit;
                }
                bit = bit.rotate_left(1);
                if bit == 1 {
                    f.write(&[byte])?;
                    byte = 0;
                }
            }
            if bit != 1 {
                f.write(&[byte])?;
            }
            self.index_in_file += 1;
        }
        // close file
        if self.index_in_file == self.chunk_size {
            let mut f = self.current_file.take().unwrap();
            f.flush()?;
            self.index_in_file = 0;
            self.index_in_dir += 1;
        }
        Ok(())
    }
    pub fn flush(&mut self) -> std::io::Result<()> {
        if let Some(f) = self.current_file.as_mut() {
            f.flush()?;
        }
        Ok(())
    }
}
