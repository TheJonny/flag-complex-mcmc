use std::cmp::{max, min};
use rand;
use rand::prelude::*;

pub mod io;

mod util;

pub mod graph;
use graph::*;

mod flagser;

//pub mod flag_complex;

type Graph = CompactMatrixGraph;
//type Graph = BoolMatrixGraph;

#[derive(Clone, Debug)]
pub struct State {
    pub cliques: Vec<Vec<Node>>,
    pub which_cliques: Vec<Vec<CliqueId>>,
    pub edge_neighborhood: Vec<Vec<Node>>,
    pub graph: Graph,
    pub flag_count: Vec<usize>,
    pub flag_count_min: Vec<usize>,
    pub flag_count_max: Vec<usize>,
}

impl State {
    pub fn new(graph: Graph) -> Self {

        println!("undirected maximal cliques");
        let cliques = graph.compute_maximal_cliques();
        let mut which_cliques = vec![Vec::new(); graph.nnodes()];
        for (c,cid) in cliques.iter().zip(0..) {
            for &v in c  {
                which_cliques[v as usize].push(cid);
            }
        }
        println!("initial flagser");
        let flag_count = graph.flagser_count();
        //let flag_count = flag_complex::count_cells(&graph);
        //assert_eq!(flag_count, flag_count2);
        let bandwidth = 1.10;
        let mut flag_count_max: Vec<_> = flag_count.iter().map(|x| (*x as f64 * bandwidth + 10.) as usize).collect();
        flag_count_max.push(2); // some higher dimensional simplices should be ok
        let flag_count_min: Vec<_> = flag_count.iter().map(|x| (*x as f64 / bandwidth - 10.) as usize).collect();
        println!("We have {:?},\n lower limit {:?},\n upper limit {:?}\n", &flag_count, &flag_count_min, &flag_count_max);

        println!("computing edge neighborhoods");
        let edge_neighborhood = graph.compute_edge_neighborhoods();

        State { graph, cliques, which_cliques, flag_count, flag_count_min, flag_count_max, edge_neighborhood}
    }

    /// applies transition, returns the change in simplex counts
    pub fn apply_transition(&mut self, t: &Transition) -> (Vec<usize>, Vec<usize>) {
        let nei = self.edgeset_neighborhood(&t.change_edges.iter().map(|&([a,b], _)| [max(a,b), min(a,b)]).collect::<Vec<Edge>>());
        let pre = Graph::subgraph(&self.graph, &nei).flagser_count();
        for (p,s) in pre.iter().zip(self.flag_count.iter_mut()) {
            assert!(*s >= *p);
            *s -= *p;
        }
        for &([a,b],add) in &t.change_edges {
            self.graph.set_edge(a, b, add);
        }
        let post = Graph::subgraph(&self.graph, &nei).flagser_count();
        if post.len() > self.flag_count.len() {
            self.flag_count.resize(post.len(), 0);
        }
        for (p,s) in post.iter().zip(self.flag_count.iter_mut()) {
            *s += *p;
        }

        return (pre, post);
    }

    pub fn revert_transition(&mut self, t: &Transition, &(ref pre, ref post): &(Vec<usize>, Vec<usize>)) {
        for &([a,b],add) in &t.change_edges {
            self.graph.set_edge(a, b, !add);
        }

        for (p,s) in post.iter().zip(self.flag_count.iter_mut()) {
            assert!(*s >= *p);
            *s -= *p;
        }

        if pre.len() > self.flag_count.len() {
            self.flag_count.resize(pre.len(), 0);
        }
        for (p,s) in pre.iter().zip(self.flag_count.iter_mut()) {
            *s += *p;
        }
    }

    fn valid(&self) -> bool {
        return all_le(&self.flag_count_min, &self.flag_count, &0)
            && all_le(&self.flag_count, &self.flag_count_max, &0);
    }

    /// returns the subgraph affected by shuffling vertices/edges in clique cid.
    /// it consists of all cliques, that share at least one edge with cid.
    pub fn clique_neighborhood(&self, cid: CliqueId) -> Vec<Node> {

        let mut affected_vertices = self.cliques[cid].clone();
        for &a in &self.cliques[cid as usize] {
            for &b in &self.cliques[cid as usize] {
                if a > b {
                    affected_vertices.extend_from_slice(&self.edge_neighborhood[edge_id(a, b)]);
                }
            }
        }
        affected_vertices.sort_unstable();
        affected_vertices.dedup();
        return affected_vertices;
    }

    pub fn edgeset_neighborhood(&self, edges: &[Edge]) -> Vec<Node>{
        let mut affected_vertices = vec![];
        for &[a,b] in edges {
            affected_vertices.extend_from_slice(&self.edge_neighborhood[edge_id(a, b)]);
            affected_vertices.push(a);
            affected_vertices.push(b);
        }
        affected_vertices.sort_unstable();
        affected_vertices.dedup();
        return affected_vertices;
    }
}

fn all_le<T: PartialOrd> (a: &[T], b: &[T], z: &T) -> bool{
    let maxlen = std::cmp::max(a.len(), b.len());
    let left = a.iter().chain(std::iter::repeat(z));
    let right = b.iter().chain(std::iter::repeat(z));
    for (l,r) in left.zip(right).take(maxlen) {
        if l > r {
            return false;
        }
    }
    return true;
}

pub struct MCMCSampler<R: Rng> {
    pub state: State,
    pub burn_in: usize,
    pub sample_distance: usize,
    pub sampled: usize,
    pub accepted: usize,
    pub rng: R,
}

impl<R: Rng> MCMCSampler<R> {
    pub fn burn_in(&mut self) {
        while self.burn_in > 0 {
            self.burn_in -= 1;
            let t = Transition::random_clique_shuffling(&self.state, &mut self.rng);
            let counters = self.state.apply_transition(&t);
            if ! self.state.valid() {
                self.state.revert_transition(&t, &counters);
            }
        }
    }
    pub fn next(&mut self) -> &State{
        for _ in 0..self.sample_distance {
            let t = Transition::random_clique_shuffling(&self.state, &mut self.rng);
            let counters = self.state.apply_transition(&t);
            self.sampled += 1;
            if self.state.valid() {
                self.accepted += 1;
            }
            else {
                self.state.revert_transition(&t, &counters);
            }
        }
        return &self.state;
    }
    pub fn acceptance_ratio(&self) -> f64 {
        return self.accepted as f64 / self.sampled as f64;
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct Transition {
    /// true: add this edge; false: remove this edge
    pub change_edges: Vec<(Edge, bool)>,
}

impl Transition {
    /// new.edge(s[i], s[j]) <-> old.edge(i,j)
    pub fn new_clique_shuffling(state: &State, cid: CliqueId, perm: &[usize]) -> Self{
        let cl = &state.cliques[cid as usize];

        // set new edge directions:
        // create a new subgraph adjacency matrix as hashmap and then apply it to self.
        let mut change_edges = vec![];

        // prepare new adjacency matrix...
        for i in 0 .. cl.len() {
            for j in 0 .. cl.len() {
                let pre = state.graph.edge(cl[perm[i]], cl[perm[j]]);
                let post = state.graph.edge(cl[i], cl[j]);
                if pre != post {
                    change_edges.push(([cl[perm[i]], cl[perm[j]]], post));
                }
            }
        }
        Transition {change_edges}
    }
    pub fn random_clique_shuffling<R: Rng>(state: &State, rng: &mut R) -> Self {
        let cid = rng.gen_range(0..state.cliques.len() as CliqueId);
        let n = state.cliques[cid].len();
        let mut perm = (0..n).collect::<Vec<usize>>();
        perm.shuffle(rng);
        return Transition::new_clique_shuffling(state, cid, &perm);
    }

    fn edges_from_clique(state: &State, cid: usize) -> (Vec<Edge>, Vec<Edge>) {
        // gets a clique id and returns a vector of its single and a vector of its double edges.
        // TODO: refactor, should not be part of Transition
        let vertices = &state.cliques[cid];
        let mut single_edges = vec![];
        let mut double_edges = vec![];
        for (i,&from) in vertices.iter().enumerate() {
            for &to in &vertices[..i]{
                if state.graph.edge(from,to) {
                    if state.graph.edge(to,from) {
                        double_edges.push([from,to]);
                    } else {
                        single_edges.push([from,to]);
                    }
                } else {
                    single_edges.push([to,from]);
                }
            }
        }
        return (single_edges, double_edges);
    }

    pub fn random_edge_swap<R: Rng>(state: &State, rng: &mut R) -> Self {
        let cid = rng.gen_range(0..state.cliques.len() as CliqueId);
        let (single_edges, _) = Transition::edges_from_clique(state, cid);
        if let Some(&[from,to]) = single_edges.choose(rng) {
            return Transition{change_edges: vec![([from,to], false), ([to,from], true)]};
        } else {
            return Transition{change_edges: vec![]};
        }

    }
    
    /* TODO: Jonathan machen lassen
    pub fn random_double_edge_move<R: Rng>(state: &State, rng: &mut R) -> Self {
        let [cid1,cid2] = rand::seq::index::sample(rng, state.cliques.len(), 2).into_vec()[..];
        let (single_edges1, double_edges1) = Transition::edges_from_clique(state, cid1);
        let (single_edges2, double_edges2) = Transition::edges_from_clique(state, cid2);
        if let Some(&[from1,to1]) = single_edges_1.choose(rng) && let Some(&[from2,to2]) = double_edges_2.choose(rng) 
                || let Some(&[from1,to1]) = single_edges_2.choose(rng) && let Some(&[from2,to2]) = double_edges_1.choose(rng) {
            [from2, to2] = [[from2,to2], [to2,from2]].choose(rng);
            return Transition{change_edges: vec![([from2,to2], false), ([to1,from1], true)]};
        } else {
            return Transition{change_edges: vec![]};
        }
    }
    */
}
