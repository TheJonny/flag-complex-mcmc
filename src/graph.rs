pub type Node = u32;
pub type CliqueId = usize;

pub type Edge = [Node; 2];

//use rand::Rng;
//use rand::{Rng, prelude::SliceRandom};
use rand::prelude::*;
use serde::{Serialize, Deserialize};

use indexmap::set::IndexSet;

use std::cmp::{min, max};

use seahash;

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

#[derive(Clone)]
pub struct BoolMatrixGraph {
    pub adjmat: Vec<bool>,
    pub nnodes: usize,
}
impl std::fmt::Debug for BoolMatrixGraph{
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        return write!(f, "A graph with {} nodes and the following edges: {:?}", self.nnodes, self.edges());
    }
}

impl BoolMatrixGraph {
    pub fn edge_mut (&mut self, from: Node, to: Node) -> &mut bool{
        if from as usize >= self.nnodes {
            panic!("from out of bounds");
        }
        if to as usize >= self.nnodes {
            panic!("to out of bounds");
        }
        return &mut self.adjmat[(from as usize) * self.nnodes + (to as usize)];
    }
}
impl DirectedGraphNew for BoolMatrixGraph {
    fn new_disconnected(nnodes: usize) -> Self{
        let adjsize = nnodes.checked_mul(nnodes).expect("adjacency matrix size overflows");
        let adjmat = vec![false; adjsize];
        return BoolMatrixGraph{adjmat, nnodes};
    }
}

impl DirectedGraph for BoolMatrixGraph {
    fn nnodes(&self) -> usize {
        return self.nnodes;
    }
    fn has_edge(&self, from: Node, to: Node) -> bool{
        if from as usize >= self.nnodes {
            panic!("from out of bounds");
        }
        if to as usize >= self.nnodes {
            panic!("to out of bounds");
        }
        return self.adjmat[(from as usize) * self.nnodes + (to as usize)];
    }

    fn add_edge(&mut self, from: Node, to: Node) {
        *self.edge_mut(from, to) = true;
    }
    fn remove_edge(&mut self, from: Node, to: Node) {
        *self.edge_mut(from, to) = false;
    }
    fn set_edge(&mut self, from: Node, to: Node, create: bool) {
        *self.edge_mut(from, to) = create;
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
    fn compute_maximal_cliques2(&self) -> Vec<Vec<Node>> {
        // TODO: Kommentierung
        // undirected flagser statt eigenbau?
        let mut r = vec![];
        let x = vec![];
        let p = self.iter_nodes().collect();
        let mut res = vec![];
        let mut rng = rand::thread_rng();
        self.bron_kerbosch_rand(&mut r, p, x, &mut res, &mut rng);
        return res;
    }
    /// recursion for compute_maximal_cliques
    fn bron_kerbosch_rand<R: rand::Rng>(&self, r: &mut Vec<Node>, mut p: Vec<Node>, mut x: Vec<Node>, res: &mut Vec<Vec<Node>>, rng: &mut R) {
        // from wikipedia: https://en.wikipedia.org/wiki/Bron%E2%80%93Kerbosch_algorithm
        // algorithm BronKerbosch2(R, P, X) is
        // if P and X are both empty then
        //     report R as a maximal clique
        // choose a pivot vertex u in P ⋃ X
        // for each vertex v in P \ N(u) do
        //     BronKerbosch2(R ⋃ {v}, P ⋂ N(v), X ⋂ N(v))
        //     P := P \ {v}
        //     X := X ⋃ {v}
        
        if p.is_empty() && x.is_empty() {
            res.push(r.clone());
            return;
        }
        let pivot_i = rng.gen_range(0..p.len() + x.len());
        let pivot = if pivot_i < p.len() { p[pivot_i] } else {x[pivot_i-p.len()]};

        for i in 0..p.len() {
            let v = p[i];
            if self.has_edge(v, pivot) || self.has_edge(pivot, v) {
                // skipping will exclude v=p[i] from p[i..].
                // so add it at the back to keep it.
                p.push(v);
                continue;
            }
            let newp = p[i..].iter().cloned().filter(|&u| self.has_edge(u,v) || self.has_edge(v,u)).collect();
            let newx = x.iter().cloned().filter(|&u| self.has_edge(u,v) || self.has_edge(v,u)).collect();
            let addv = !r.contains(&v);
            if addv {
                r.push(v);
            }
            self.bron_kerbosch_rand(r, newp, newx, res, rng);
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

/// undirected edge pair to index in triangular adjacency matrix
/// ```
/// use directed_scm::graph::*;
/// assert_eq!(edge_id(1, 0), 0);
/// assert_eq!(edge_id(2, 0), 1);
/// assert_eq!(edge_id(2, 1), 2);
/// assert_eq!(edge_id(3, 0), 3);
/// assert_eq!(edge_id(3, 1), 4);
/// assert_eq!(edge_id(3, 2), 5);
/// ```
pub fn edge_id(a: Node, b: Node) -> usize{
    assert!(a > b);
    return (a as usize) * ((a as usize)-1) / 2 + (b as usize);
    
}

type Chunk = u64;
const CHUNK_SIZE: usize = Chunk::BITS as usize;

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct CompactMatrixGraph {
    nnodes: usize,
    row_len: usize,
    out_matrix: Vec<Chunk>,
    //out_degrees: Vec<usize>,
}

impl DirectedGraphNew for CompactMatrixGraph {
    fn new_disconnected(nnodes: usize) -> Self {
        // nnodes / Chunk::BITS, rounded up.
        let row_len = (nnodes + Chunk::BITS as usize - 1) / (Chunk::BITS as usize);

        let matsize = row_len.checked_mul(nnodes).expect("size of adjacency matrix overflows");
        let out_matrix = vec![0; matsize];
        //let out_degrees = vec![0; nnodes];
        CompactMatrixGraph { out_matrix, row_len, nnodes}
    }
}

impl DirectedGraph for CompactMatrixGraph {
    fn nnodes(&self) -> usize {
        return self.nnodes;
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
        let toh = to as usize / CHUNK_SIZE;
        let tol = to as usize % CHUNK_SIZE;
        let row = from as usize;
        self.out_matrix[row * self.row_len + toh] |= 1<<tol;
    }
    fn remove_edge(&mut self, from: Node, to: Node) {
        let toh = to as usize / CHUNK_SIZE;
        let tol = to as usize % CHUNK_SIZE;
        let row = from as usize;
        self.out_matrix[row * self.row_len + toh] &= !(1<<tol);
    }

    fn edges(&self) -> Vec<[Node; 2]> {
        let mut result = vec![];
        for from in 0 .. self.nnodes {
            for toh in 0 .. self.row_len {
                let mut bits = self.out_matrix[from * self.row_len + toh];
                while bits != 0 {
                    let b = bits.trailing_zeros();
                    bits &= !(1<<b);
                    let to = toh as u32 * Chunk::BITS | b;
                    result.push([from as Node, to]);
                }
            }
        }
        result
    }
    /// ```
    /// use directed_scm::graph::*;
    /// let mut g = CompactMatrixGraph::new_disconnected(5);
    /// g.add_edge(2, 3);
    /// assert_eq!(g.edges_from(2), vec![3]);
    /// ```
    fn edges_from(&self, from: Node) -> Vec<Node> {
        let offset = self.row_len * from as usize;
        let mut result = vec![];
        for toh in 0 .. self.row_len {
            let mut bits = self.out_matrix[offset + toh];
            while bits != 0 {
                let b = bits.trailing_zeros();
                bits &= !(1<<b);
                let to = toh as u32 * Chunk::BITS | b;
                result.push(to);
            }
        }
        result
    }
}

#[derive(Debug, Serialize, Deserialize)]
pub struct EdgeMapGraph{
    nnodes: usize,
    edges: IndexSet<Edge, std::hash::BuildHasherDefault<seahash::SeaHasher>>,
    double_edges: IndexSet<Edge, std::hash::BuildHasherDefault<seahash::SeaHasher>>,

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

    // fn choose_double_edge ->
}

