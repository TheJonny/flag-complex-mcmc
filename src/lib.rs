use std::cmp::{max, min};
use std::collections::HashMap;

use rand;
use rand::distributions::WeightedIndex;
use rand::prelude::*;
use serde::{Serialize, Deserialize};

pub mod io;

mod util;

pub mod graph;
use graph::*;

mod flagser;

//pub mod flag_complex;

type Graph = EdgeMapGraph;
//type Graph = BoolMatrixGraph;

#[derive(Debug, Serialize, Deserialize)]
pub struct State {
    pub cliques: Vec<Vec<Node>>,
    pub which_cliques: Vec<Vec<CliqueId>>,
    pub edge_neighborhood: HashMap<Edge, Vec<Node>>,
    pub graph: Graph,
    pub flag_count: Vec<usize>,
    pub flag_count_min: Vec<usize>,
    pub flag_count_max: Vec<usize>,
}

impl State {
    pub fn new(graph: CompactMatrixGraph) -> Self {

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

        let graph = Graph::copy(&graph);

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
                    affected_vertices.extend_from_slice(&self.edge_neighborhood[&[a, b]]);
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
            let big = max(a, b);
            let small = min(a,b);
            affected_vertices.extend_from_slice(&self.edge_neighborhood[&[big, small]]);
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

#[derive(Debug, Serialize, Deserialize)]
pub struct MCMCSampler<R: Rng> {
    pub state: State,
    pub move_distribution: WeightedIndex<f64>,
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
    pub fn random_next_move<R: Rng>(state: &State, rng: &mut R, move_distribution: &WeightedIndex<f64>) -> Self {
        let choices = [Transition::single_edge_flip, Transition::double_edge_move, Transition::new_clique_shuffling];
        let random_move = choices[move_distribution.sample(&mut rng)];
        return random_move(state, rng);
    }

    /// new.has_edge(s[i], s[j]) <-> old.has_edge(i,j)
    pub fn new_clique_shuffling(state: &State, cid: CliqueId, perm: &[usize]) -> Self{
        let cl = &state.cliques[cid as usize];

        // set new edge directions:
        // create a new subgraph adjacency matrix as hashmap and then apply it to self.
        let mut change_edges = vec![];

        // prepare new adjacency matrix...
        for i in 0 .. cl.len() {
            for j in 0 .. cl.len() {
                let pre = state.graph.has_edge(cl[perm[i]], cl[perm[j]]);
                let post = state.graph.has_edge(cl[i], cl[j]);
                if pre != post {
                    change_edges.push(([cl[perm[i]], cl[perm[j]]], post));
                }
            }
        }
        return Transition{change_edges};
    }
    pub fn random_clique_shuffling<R: Rng>(state: &State, rng: &mut R) -> Self {
        let cid = rng.gen_range(0..state.cliques.len() as CliqueId);
        let n = state.cliques[cid].len();
        let mut perm = (0..n).collect::<Vec<usize>>();
        perm.shuffle(rng);
        return Transition::new_clique_shuffling(state, cid, &perm);
    }

    pub fn single_edge_flip<R: Rng>(state: &State, rng: &mut R) -> Self {
        if let Some([from, to]) = state.graph.sample_edge(rng) {
            if !state.graph.has_edge(to, from) { // its a single edge
                return Transition{change_edges: vec![([from,to], false), ([to,from], true)]};
            }
        }
        return Transition{change_edges: vec![]};
    }
    
    pub fn double_edge_move<R: Rng>(state: &State, rng: &mut R) -> Self {
        // if there is no edge, return an empty transition below
        if let Some(double_edge) = state.graph.sample_double_edge(rng) {
            // FIXME: assert somewhere, that there are single edges.
            let [a,b] = loop {
                let [a,b] = state.graph.sample_edge(rng).expect("there was a double edge, so there are edges");
                if !state.graph.has_edge(b, a) {
                    break [a,b];
                }
            };
            // a,b is a single edge -> make it double
            // then remove a random side of double_edge
            let delme = if rng.gen_bool(0.5) {
                double_edge
            } else {
                [double_edge[1], double_edge[0]]
            };

            return Transition { change_edges: vec![([b,a], true), (delme, false)] }
        }
        return Transition{change_edges: vec![]};
    }
}
