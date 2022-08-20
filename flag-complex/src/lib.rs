mod graph;

pub mod prelude {
    pub use super::graph::{DirectedGraph, DirectedGraphExt, DirectedGraphNew};
}

pub use graph::{Node, Edge};

pub type Graph = graph::EdgeMapGraph;

mod complex;
pub use complex::*;

