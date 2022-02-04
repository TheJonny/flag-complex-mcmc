use std::fs::File;
use std::io::prelude::*;
use crate::graph::*;

use crate::MCMCSampler;
use rand::Rng;
use serde::Serialize;
use serde::de::DeserializeOwned;
use bincode;

use hdf5;
use ndarray;

// Flag File
pub fn read_flag_file<G: DirectedGraphNew>(fname:&str) -> G {
    let mut file = File::open(fname).unwrap();
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
    for [i,j] in graph.edges() {
        content += &(format!("{} {} 1\n", i, j));
    }
    std::fs::write(fname, content).expect("Unable to write file");
    return Ok(());
}

// Save/Restore state (serde/bincode)
pub fn save_state<R: Rng+Serialize>(fname:&str, sample_number:usize, sampler:MCMCSampler<R>) -> std::io::Result<()> {
    let mut f = std::io::BufWriter::new(File::create(format!("{fname}.tmp")).unwrap());
    bincode::serialize_into(&mut f, &(sample_number, sampler)).unwrap();
    std::fs::rename(format!("{fname}.tmp"), fname).expect("moving temp state file to correct location failed");
    return Ok(());
}

pub fn load_state<R: Rng+DeserializeOwned>(fname:&str) -> std::io::Result<(usize, MCMCSampler<R>)> {
    let mut f = std::io::BufReader::new(File::open(fname).unwrap());
    let (sample_index_end, sampler) = bincode::deserialize_from(&mut f).unwrap();
    return Ok((sample_index_end, sampler))
}


// HDF5
pub fn save_to_hdf<G: DirectedGraph>(label:&str, seed:u64, sample_number:usize, graph:&G, flag_count:&Vec<usize>) -> hdf5::Result<()> {
    let file = hdf5::File::open_rw(format!("{label}-{seed:03}.hdf5"))?;
    let group = file.create_group(&format!("/{seed:03}/{sample_number:05}"))?;
    let builder = group.new_dataset_builder();
    let ds = builder.deflate(4).with_data(&ndarray::arr2(&graph.edges())).create("edgelist")?;
    ds.new_attr::<usize>().shape(flag_count.len()).create("flag_count")?.write(flag_count)?;
    ds.new_attr::<usize>().shape(1).create("number_of_vertices")?.write(&[graph.nnodes()])?;
    Ok(())
}

pub fn new_hdf_file(label:&str, seed:u64) -> hdf5::Result<()> {
    hdf5::File::create(format!("{label}-{seed:03}.hdf5"))?;
    //TODO: Add metadata/comment to this file
    return Ok(());
}
