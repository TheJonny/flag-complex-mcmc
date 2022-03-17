pub type Node = u32;
pub type Edge = [Node; 2];

//use rand::Rng;
//use rand::{Rng, prelude::SliceRandom};
use rand::prelude::*;
use serde::{Serialize, Deserialize};

use indexmap::set::IndexSet;

use std::cmp::{min, max};

pub trait DirectedGraph: Sync {
    fn has_edge(&self, from: Node, to: Node) -> bool;
    fn add_edge(&mut self, from: Node, to: Node);
    fn remove_edge(&mut self, from: Node, to: Node);
    fn nnodes(&self) -> usize;
    fn iter_nodes(&self) -> std::ops::Range<Node> {
        (0.. (self.nnodes() as Node)).into_iter()
    }

    fn edges(&self) -> Vec<[Node; 2]> {
        let mut edges = vec![];
        for from in self.iter_nodes() {
            for to in self.iter_nodes() {
                if self.has_edge(from, to) {
                    edges.push([from, to]);
                }
            }
        }
        return edges
    }
    fn set_edge(&mut self, from: Node, to: Node, create: bool) {
        if create {
            self.add_edge(from, to);
        }
        else {
            self.remove_edge(from, to);
        }
    }

    fn edges_from(&self, from: Node) -> Vec<Node>{
        let mut res = vec![];
        for to in self.iter_nodes() {
            if self.has_edge(from, to) {
                res.push(to);
            }
        }
        return res;
    }

}

pub trait DirectedGraphNew: DirectedGraph + Sized {
    fn new_disconnected(nnodes: usize) -> Self;
    fn subgraph<G: DirectedGraph>(ori: &G, vertices: &[Node]) -> Self {
        // macht komische dinge wenn Vertices in dem Slice doppelt vorkommen
        // FIXME?
        let mut sub = Self::new_disconnected(vertices.len());
        for (new_from, &ori_from) in (0..).zip(vertices.iter()) {
            for (new_to, &ori_to) in (0..).zip(vertices.iter()) {
                if ori.has_edge(ori_from, ori_to) {
                    sub.add_edge(new_from, new_to);
                }
            }
        }
        return sub;
    }
    fn copy<G: DirectedGraph>(g: &G) -> Self {
        let mut n = Self::new_disconnected(g.nnodes());
        for [a,b] in g.edges() {
            n.add_edge(a, b);
        }
        return n
    }
}


pub trait DirectedGraphExt: DirectedGraph {
    fn compute_maximal_cliques(&self) -> Vec<Vec<Node>> {
        // TODO: Kommentierung
        // undirected flagser statt eigenbau?
        let mut r = vec![];
        let x = vec![];
        let p = (0..(self.nnodes() as Node)).collect();
        let mut res = vec![];
        self.bron_kerbosch(&mut r, p, x, &mut res);
        return res;
    }
    fn compute_cliques(&self) -> Vec<Vec<Node>> {
        let mut edges = self.edges();
        for e in &mut edges {
            e.sort();
        }
        let mut normalized_graph = EdgeMapGraph::new_disconnected(self.nnodes());
        for [a,b] in edges {
            normalized_graph.add_edge(a, b);
        }
        let mut result = vec![];
        let mut add = |cell: &[Node]| {
            result.push(cell.into())
        };
        crate::flag_complex::for_each_cell(&normalized_graph, &mut add, 0, self.nnodes());
        return result
    }
    /// recursion for compute_maximal_cliques
    fn bron_kerbosch(&self, r: &mut Vec<Node>, p: Vec<Node>, mut x: Vec<Node>, res: &mut Vec<Vec<Node>>) {
        // from wikipedia: https://en.wikipedia.org/wiki/Bron%E2%80%93Kerbosch_algorithm
        // algorithm BronKerbosch1(R, P, X) is
        // if P and X are both empty then
        //     report R as a maximal clique
        // for each vertex v in P do
        //     BronKerbosch1(R ⋃ {v}, P ⋂ N(v), X ⋂ N(v))
        //     P := P \ {v}
        //     X := X ⋃ {v}
        
        if p.is_empty() && x.is_empty() {
            res.push(r.clone());
        }
        for (i,&v) in (0..).zip(&p) {
            let newp = p[i..].iter().cloned().filter(|&u| self.has_edge(u,v) || self.has_edge(v,u)).collect();
            let newx = x.iter().cloned().filter(|&u| self.has_edge(u,v) || self.has_edge(v,u)).collect();
            let addv = !r.contains(&v);
            if addv {
                r.push(v);
            }
            self.bron_kerbosch(r, newp, newx, res);
            if addv {
                r.pop();
            }
            x.push(v);
        }
    }
    
    fn flagser_count(&self) -> Vec<usize> {
        //crate::flagser::count_unweighted(self.nnodes(), &self.edges())
        crate::flag_complex::count_cells(self)
    }

}

impl<G: DirectedGraph> DirectedGraphExt for G {}

type Chunk = u64;
const CHUNK_SIZE: usize = Chunk::BITS as usize;

#[derive(Debug, Serialize, Deserialize)]
pub struct EdgeMapGraph{
    nnodes: usize,
    edges: IndexSet<Edge>,
    double_edges: IndexSet<Edge>,

    out_matrix: Vec<Chunk>,
    row_len: usize,
}

impl DirectedGraphNew for EdgeMapGraph {
    fn new_disconnected(nnodes: usize) -> Self {
        // nnodes / Chunk::BITS, rounded up.
        let row_len = (nnodes + Chunk::BITS as usize - 1) / (Chunk::BITS as usize);

        let matsize = row_len.checked_mul(nnodes).expect("size of adjacency matrix overflows");
        let out_matrix = vec![0; matsize];
        //let out_degrees = vec![0; nnodes];
        EdgeMapGraph { nnodes, edges: Default::default(), double_edges: Default::default() , out_matrix, row_len}
    }
}

impl DirectedGraph for EdgeMapGraph {
    fn nnodes(&self) -> usize {
        self.nnodes
    }
    fn has_edge(&self, from: Node, to: Node) -> bool {
        if from as usize >= self.nnodes {
            panic!("from out of bounds: {} >= {}", from, self.nnodes);
        }
        if to as usize >= self.nnodes {
            panic!("to out of bounds: {} >= {}", to, self.nnodes);
        }
        let toh = to as usize / CHUNK_SIZE;
        let tol = to as usize % CHUNK_SIZE;
        let row = from as usize;
        return self.out_matrix[row * self.row_len + toh] & (1<<tol) != 0;
    }
    fn add_edge(&mut self, from: Node, to: Node) {
        // update hashmaps
        self.edges.insert([from, to]);
        if self.has_edge(to, from) {
            let big = max(from, to);
            let small = min(from, to);
            self.double_edges.insert([big, small]);
        }

        // update bits
        let toh = to as usize / CHUNK_SIZE;
        let tol = to as usize % CHUNK_SIZE;
        let row = from as usize;
        self.out_matrix[row * self.row_len + toh] |= 1<<tol;
    }
    fn remove_edge(&mut self, from: Node, to: Node){
        // update hashmaps
        if self.has_edge(from, to) && self.has_edge(to, from) {
            let big = max(from, to);
            let small = min(from, to);
            self.double_edges.remove(&[big, small]); // does nothing if it does not exist
        }
        self.edges.remove(&[from, to]);

        //  update bits
        let toh = to as usize / CHUNK_SIZE;
        let tol = to as usize % CHUNK_SIZE;
        let row = from as usize;
        self.out_matrix[row * self.row_len + toh] &= !(1<<tol);
    }

    // FIXME: API ändern, dass nicht kopiert werden muss
    fn edges(&self) -> Vec<[Node; 2]> {
        self.edges.iter().cloned().collect()
    }
}

impl EdgeMapGraph {
    pub fn sample_edge<R: Rng>(&self, rng: &mut R) -> Option<Edge> {
        if self.edges.len() == 0 {
            return None;
        }
        let i = rng.gen_range(0 .. self.edges.len());
        return Some(self.edges[i]);
    }
    pub fn sample_double_edge<R: Rng>(&self, rng: &mut R) -> Option<Edge> {
        if self.double_edges.len() == 0 {
            return None;
        }
        let i = rng.gen_range(0 .. self.double_edges.len());
        return Some(self.double_edges[i]);
    }

}
