pub type Node = u32;
pub type CliqueId = usize;

pub type Edge = [Node; 2];

use rayon::prelude::*;

#[derive(Clone)]
pub struct DirectedGraph {
    pub adjmat: Vec<bool>,
    pub nnodes: usize,
}
impl std::fmt::Debug for DirectedGraph{
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        return write!(f, "A graph with {} nodes and the following edges: {:?}", self.nnodes, self.edges());
    }
}

impl DirectedGraph {
    pub fn new_disconnected(nnodes: usize) -> Self{
        let adjsize = nnodes.checked_mul(nnodes).expect("adjacency matrix size overflows");
        let adjmat = vec![false; adjsize];
        return DirectedGraph{adjmat, nnodes};
    }
    pub fn edge(&self, from: Node, to: Node) -> bool{
        if from as usize >= self.nnodes {
            panic!("from out of bounds");
        }
        if to as usize >= self.nnodes {
            panic!("to out of bounds");
        }
        return self.adjmat[(from as usize) * self.nnodes + (to as usize)];
    }
    pub fn edge_mut (&mut self, from: Node, to: Node) -> &mut bool{
        if from as usize >= self.nnodes {
            panic!("from out of bounds");
        }
        if to as usize >= self.nnodes {
            panic!("to out of bounds");
        }
        return &mut self.adjmat[(from as usize) * self.nnodes + (to as usize)];
    }

    pub fn add_edge(&mut self, from: Node, to: Node) {
        *self.edge_mut(from, to) = true;
    }

    pub fn subgraph(&self, vertices: &[Node]) -> DirectedGraph {
        // macht komische dinge wenn Vertices in dem Slice doppelt vorkommen
        // FIXME?
        let mut sub = DirectedGraph::new_disconnected(vertices.len());
        for (new_from, &ori_from) in (0..).zip(vertices.iter()) {
            for (new_to, &ori_to) in (0..).zip(vertices.iter()) {
                *sub.edge_mut(new_from, new_to) = self.edge(ori_from, ori_to);
            }
        }
        return sub;
    }

    pub fn compute_maximal_cliques(&self) -> Vec<Vec<Node>> {
        // TODO: Kommentierung
        // undirected flagser statt eigenbau?
        let mut r = vec![];
        let x = vec![];
        let p = (0..(self.nnodes as Node)).collect();
        let mut res = vec![];
        self.bron_kerbosch(&mut r, p, x, &mut res);
        return res;
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
            let newp = p[i..].iter().cloned().filter(|&u| self.edge(u,v) || self.edge(v,u)).collect();
            let newx = x.iter().cloned().filter(|&u| self.edge(u,v) || self.edge(v,u)).collect();
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
    pub fn compute_maximal_cliques2(&self) -> Vec<Vec<Node>> {
        // TODO: Kommentierung
        // undirected flagser statt eigenbau?
        let mut r = vec![];
        let x = vec![];
        let p = (0..(self.nnodes as Node)).collect();
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
            if self.edge(v, pivot) || self.edge(pivot, v) {
                // skipping will exclude v=p[i] from p[i..].
                // so add it at the back to keep it.
                p.push(v);
                continue;
            }
            let newp = p[i..].iter().cloned().filter(|&u| self.edge(u,v) || self.edge(v,u)).collect();
            let newx = x.iter().cloned().filter(|&u| self.edge(u,v) || self.edge(v,u)).collect();
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

    pub fn edges(&self) -> Vec<[Node; 2]> {
        let mut edges = vec![];
        for from in 0..(self.nnodes as Node) {
            for to in 0..(self.nnodes as Node) {
                if self.edge(from, to) {
                    edges.push([from, to]);
                }
            }
        }
        return edges
    }

    pub fn flagser_count(&self) -> Vec<usize> {
        crate::flagser::count_unweighted(self.nnodes, &self.edges())
    }

    pub fn iter_nodes(&self) -> impl Iterator<Item=Node> {
        (0.. (self.nnodes as Node)).into_iter()
    }

    /// for every edge, this gathers the nodes that are connected to both ends.
    pub fn compute_edge_neighborhoods(&self) -> Vec<Vec<Node>>{
        let mut undirected_adj_lists = vec![vec![]; self.nnodes];
        let mut undirected_edges = self.edges();
        use std::cmp::{min, max};
        for e in &mut undirected_edges {
            let a = max(e[0], e[1]);
            let b = min(e[0], e[1]);
            *e = [a, b];
        }
        undirected_edges.sort_unstable();
        undirected_edges.dedup();

        for [a,b] in self.edges() {
            undirected_adj_lists[a as usize].push(b);
            undirected_adj_lists[b as usize].push(a);
        }
        for v in &mut undirected_adj_lists {
            v.sort_unstable();
        }
        
        let mut respairs = vec![];
        undirected_edges.par_iter().map(|&[a, b]| {
            assert!(a > b);

            let mut l = crate::util::intersect_sorted(&undirected_adj_lists[a as usize], &undirected_adj_lists[b as usize]);
            l.shrink_to_fit();
            (edge_id(a, b), l)
        }).collect_into_vec(&mut respairs);
        let mut res = vec![vec![]; self.nnodes * (self.nnodes-1) / 2];
        for (eid, l) in respairs {
            res[eid] = l;
        }
        
        return res;
    }
}

pub fn edge_id(a: Node, b: Node) -> usize{
    assert!(a > b);
    return (a as usize) * ((a as usize)-1) / 2 + (b as usize);
    
}
