mod graph;

pub mod prelude {
    pub use super::graph::{DirectedGraph, DirectedGraphExt, DirectedGraphNew};
}

pub use graph::{Node, Edge};

pub type Graph = graph::EdgeMapGraph;

pub mod flag_complex;
pub use flag_complex::*;

