use std::fs::File;
use std::io::prelude::*;
use crate::DirectedGraph;

pub fn read_flag_file(fname:&str) -> DirectedGraph {
    let mut file = File::open(fname).unwrap();
    let mut fcontents = String::new();
    file.read_to_string(&mut fcontents).unwrap();

    let mut lines = fcontents.lines();
    lines.next(); // skip dim 0
    let nnodes = lines.next().unwrap().split(' ').filter(|s| *s != "").count(); // todo: mehr beleidigungen // todo: better parsing
    let mut graph = DirectedGraph::new_disconnected(nnodes);
    lines.next(); // skip dim 1
    for line in lines {
        let mut ijw = line.split(' ').filter(|s| *s != "");
        if let (Some(i), Some(j)) = (ijw.next(), ijw.next()) {
            graph.add_edge(i.parse().unwrap(), j.parse().unwrap());
        }
    }
    return graph;
}


pub fn save_flag_file(fname:&str, graph:&DirectedGraph) -> std::io::Result<()> {
    let mut content = "dim 0:\n".to_string();
    content += &("1 ".repeat(graph.nnodes).trim_end().to_owned() + "\n"); // add vertices
    content += "dim 1:\n";
    for [i,j] in graph.edges() {
        content += &(format!("{} {} 1\n", i, j));
    }
    std::fs::write(fname, content).expect("Unable to write file");
    return Ok(());
}
