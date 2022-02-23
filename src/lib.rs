use std::cmp::{max, min};
use std::collections::HashMap;

use rand;
use rand::distributions::WeightedIndex;
use rand::prelude::*;
use serde::{Serialize, Deserialize};

use rayon::prelude::*;

use crate::util::{all_le, intersect_sorted, random_perm, vec_intersect, vec_setminus};

pub mod io;

mod util;

pub mod graph;
use graph::*;

pub mod flag_complex;

type Graph = EdgeMapGraph;
//type Graph = BoolMatrixGraph;

#[derive(Debug, Serialize, Deserialize)]
pub struct EdgeInfo {
    nbhd: Vec<Node>,
    ncliques_in_nbhd: usize,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct State {
    pub max_cliques_by_order: Vec<Vec<Vec<Node>>>,
    pub max_clique_weights_by_order: Vec<Vec<f64>>,
    pub edge_neighborhood: HashMap<Edge, Vec<Node>>,
    pub graph: Graph,
    pub flag_count: Vec<usize>,
    pub flag_count_min: Vec<usize>,
    pub flag_count_max: Vec<usize>,
}

impl State {
    pub fn new(graph: Graph) -> Self {

        println!("undirected maximal max_cliques");
        let max_cliques = graph.compute_maximal_cliques();
        let mut max_cliques_by_order = vec![];
        for c in max_cliques.clone() {
            let max_clique_order = c.len();
            if max_clique_order > max_cliques_by_order.len() {
                max_cliques_by_order.resize(max_clique_order, vec![]);
            }
            max_cliques_by_order[max_clique_order-1].push(c);
        }

        println!("calculate max_clique distribution");
        let undirected_edges = compute_undirected_edges(&graph);
        let mut max_cliques_per_edge = std::collections::HashMap::<Edge, Vec<&Vec<Node>>>::with_capacity(undirected_edges.len());
        for c in max_cliques.iter() { 
           for (vid, &i) in c.iter().enumerate() {
               for &j in c[..vid].iter() {
                   max_cliques_per_edge.entry([i,j]).or_default().push(c);
               }
           }
        }
        let mut max_clique_weights = std::collections::HashMap::<&Vec<Node>, f64>::with_capacity(max_cliques.len());
        for e in undirected_edges {
            let weight_per_max_clique = 1./(max_cliques_per_edge[&e].len() as f64);
            for c in &max_cliques_per_edge[&e] {
                *max_clique_weights.entry(c).or_default() += weight_per_max_clique;
            }
        }
        let max_clique_weights_by_order:Vec<Vec<f64>> = (0..max_cliques_by_order.len()).map(|d| max_cliques_by_order[d].iter().map(|c| max_clique_weights[c]).collect()).collect();

        println!("initial flagser");
        let flag_count = graph.flagser_count();

        println!("computing edge neighborhoods");
        let (edge_neighborhood, max_by_dim) = compute_edge_infos(&graph);
        dbg!(&max_by_dim);
        
        println!("compute number of (non-maximal) cliques");
        let clique_count = compute_normalized_undirected_graph(&graph).flagser_count();
        dbg!(clique_count);

        //let clique_count = compute_normalized_directed_graph(&graph).flagser_count();
        //dbg!(clique_count);

        let target_relax = 1.02; 
        let relax_de_upper = crate::util::calc_relax_de(&flag_count.iter().map(|&scd| ((scd as f64) * target_relax).round() as usize).collect());
        let relax_de_lower = crate::util::calc_relax_de(&flag_count.iter().map(|&scd| ((scd as f64) / target_relax).round() as usize).collect());
        let mut flag_count_max: Vec<usize> = vec![];
        let mut flag_count_min: Vec<usize> = vec![];
        for d in 0..flag_count.len() {
            let relax_upper = relax_de_upper[d]/2;
            let relax_lower = relax_de_upper[d]/2;
            println!("absolute relaxation (upper/lower) in dimension {d} is: {relax_upper} {relax_lower}");
            flag_count_max.push(((flag_count[d] as f64) * target_relax + relax_upper as f64).round() as usize);
            flag_count_min.push(((flag_count[d] as f64) / target_relax - relax_lower as f64).round() as usize);
        }
        //flag_count_max.push(clique_count10); TODO: ADD SOMETHING LIKE THIS
        println!("We have {:?},\n lower limit {:?},\n upper limit {:?}\n", &flag_count, &flag_count_min, &flag_count_max);

        let nchange_dims = max_by_dim.len().checked_sub(2).expect("there should be at least one edge!");

        //for (m, f) in zip_longest(max_by_dim.iter(), 
        // TODO: flag_count_max/min abhängig von max_by_dim

        State { graph, max_cliques_by_order, max_clique_weights_by_order, flag_count, flag_count_min, flag_count_max, edge_neighborhood}
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

#[derive(Debug, Serialize, Deserialize)]
pub struct Distributions {
    pub moves: WeightedIndex<f64>,
    pub max_clique_orders: WeightedIndex<f64>,
    pub max_cliques_by_order: Vec<WeightedIndex<f64>>,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct MCMCSampler<R: Rng> {
    pub state: State,
    pub dists: Distributions,
    pub sample_distance: usize,
    pub sampled: usize,
    pub accepted: usize,
    pub rng: R,
}

impl<R: Rng> MCMCSampler<R> {
    pub fn next(&mut self) -> &State{
        for _ in 0..self.sample_distance {
            let t = Transition::random_move(&self.state, &mut self.rng, &self.dists);
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
    pub fn random_move<R: Rng>(state: &State, rng: &mut R, dists: &Distributions) -> Self {
        let potential_moves = [Transition::single_edge_flip, Transition::double_edge_move,
                                Transition::clique_permute, Transition::clique_swap];
        let random_move = potential_moves[dists.moves.sample(rng)];
        return random_move(state, rng, dists);
    }

    pub fn clique_permute<R: Rng>(state: &State, rng: &mut R, dists: &Distributions) -> Self {
        let dim = dists.max_clique_orders.sample(rng);
        let cl = &state.max_cliques_by_order[dim][dists.max_cliques_by_order[dim].sample(rng)];
        let perm = random_perm(0,cl.len(), rng);

        let mut change_edges = vec![];

        for i in 0 .. cl.len() {
            for j in 0 .. cl.len() {
                let pre = state.graph.has_edge(cl[perm[i]], cl[perm[j]]);
                let post = state.graph.has_edge(cl[i], cl[j]);
                if pre != post {
                    change_edges.push(([cl[perm[i]], cl[perm[j]]], post));
                }
            }
        }
        return Transition {change_edges};
    }

    pub fn clique_swap<R: Rng>(state: &State, rng: &mut R, dists: &Distributions) -> Self {
        let dim = dists.max_clique_orders.sample(rng);
        let m1 = &state.max_cliques_by_order[dim][dists.max_cliques_by_order[dim].sample(rng)];
        let m2 = &state.max_cliques_by_order[dim][dists.max_cliques_by_order[dim].sample(rng)];
        
        let c = vec_intersect(&m1, &m2);
        let d = {let mut x = c.clone(); x.extend(&vec_setminus(&m1,&c)); x.extend(&vec_setminus(&m2,&c)); x};

        let n_c = c.len();
        let n_d = d.len();
        let n_a = m1.len() - n_c;
        
        //let perm_c = (0..n_c).collect::<Vec<usize>>();    //uncomment to not permute common vertices
        let perm_c = random_perm(0, n_c, rng);   //with perm on c
        let perm_a = random_perm(n_c, n_c+n_a, rng);
        let perm_b = random_perm(n_c+n_a, n_d, rng);

        let perm_d = {let mut x = perm_c.clone(); x.extend(&perm_b); x.extend(&perm_a); x};

        //TODO VIEEEEEL SCHÖNER
        let mut new_edges = Vec::<Edge>::new();
        let mut old_edges = Vec::<Edge>::new();
        for i in 0..n_c + n_a {
            for j in 0..n_c + n_a {
                if state.graph.has_edge(d[i], d[j]) {
                    new_edges.push([d[perm_d[i]], d[perm_d[j]]]);
                    old_edges.push([d[i], d[j]]);
                }
            }
        }
        for i in (0..n_c).chain(n_c+n_a..n_d) {
            for j in (0..n_c).chain(n_c+n_a..n_d) {
                if state.graph.has_edge(d[i], d[j]) {
                    new_edges.push([d[perm_d[i]], d[perm_d[j]]]);
                    old_edges.push([d[i], d[j]]);
                }
            }
        }
        new_edges.sort();
        new_edges.dedup();
        old_edges.sort();
        old_edges.dedup();

        let mut change_edges = vec![];
        for ne in new_edges {
            if old_edges.contains(&ne) {
                old_edges.retain(|&e| e != ne);
            } else {
                change_edges.push((ne, true));
            }
        }
        for oe in old_edges {
            change_edges.push((oe, false));
        }
        //println!("{:?}", &change_edges);
        return Transition{change_edges};
    }

    pub fn single_edge_flip<R: Rng>(state: &State, rng: &mut R, _: &Distributions) -> Self {
        if let Some([from, to]) = state.graph.sample_edge(rng) {
            if !state.graph.has_edge(to, from) { // its a single edge
                return Transition{change_edges: vec![([from,to], false), ([to,from], true)]};
            }
        }
        return Transition{change_edges: vec![]};
    }
    
    pub fn double_edge_move<R: Rng>(state: &State, rng: &mut R, _: &Distributions) -> Self {
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
/// for every edge, this gathers the nodes that are connected to both ends.
fn compute_undirected_edges(graph: &Graph) -> Vec<Edge> {
    let mut undirected_edges = graph.edges();
    for e in &mut undirected_edges {
        let a = max(e[0], e[1]);
        let b = min(e[0], e[1]);
        *e = [a, b];
    }
    undirected_edges.sort_unstable();
    undirected_edges.dedup();
    return undirected_edges
}

fn compute_normalized_undirected_graph(graph: &Graph) -> Graph {
    let undirected_edges = compute_undirected_edges(&graph);
    let mut normalized_graph = EdgeMapGraph::new_disconnected(graph.nnodes());
    for &[a,b] in &undirected_edges {
        normalized_graph.add_edge(a, b);
    }
    return normalized_graph;
}

fn compute_normalized_directed_graph(graph: &Graph) -> Graph {
    let edges = {let mut x = graph.edges(); x.sort(); x};
    let mut normalized_graph = EdgeMapGraph::new_disconnected(graph.nnodes());
    for &[a,b] in edges.iter() {
        if a < b {normalized_graph.add_edge(a,b)}
        else if !normalized_graph.has_edge(a,b) {normalized_graph.add_edge(a,b)}
        else {normalized_graph.add_edge(b,a);}
    }
    return normalized_graph;
}

fn compute_edge_infos(graph: &Graph)-> (HashMap<Edge, Vec<Node>>, Vec<usize>){
    let mut undirected_adj_lists = vec![vec![]; graph.nnodes()];
    let undirected_edges = compute_undirected_edges(&graph);

    for [a,b] in graph.edges() {
        undirected_adj_lists[a as usize].push(b);
        undirected_adj_lists[b as usize].push(a);
    }
    for v in &mut undirected_adj_lists {
        v.sort_unstable();
        v.dedup();
    }

    let mut respairs = vec![];
    undirected_edges.par_iter().map(|&[a, b]| {
        assert!(a > b);

        let mut l = intersect_sorted(&undirected_adj_lists[a as usize], &undirected_adj_lists[b as usize]);
        l.shrink_to_fit();
        ([a,b], l)
    }).collect_into_vec(&mut respairs);
    let mut edge_neighborhood = HashMap::with_capacity(respairs.len());
    for (e, l) in respairs {
        edge_neighborhood.insert(e, l);
    }
    
    // to determine the acceptance limits, we need to know, how much damage flipping any single edge
    // can do (in each dimension).
    // therefore we calculate the number of adjacent (not necessarly maximal) cliques...
    let normalized_graph = compute_normalized_undirected_graph(&graph);
    let mut ncliques_by_edge_and_dim = HashMap::<Edge, Vec<usize>>::with_capacity(undirected_edges.len());

    let mut count_for_edges = |simplex: &[Node]| {
        let dim = simplex.len() - 1;
        // iterate over simplex's edges
        for (i,&b) in simplex.iter().enumerate() {
            for &a in &simplex[0..i] {
                assert!(a>b);
                let by_dim = ncliques_by_edge_and_dim.entry([a,b]).or_insert_with(Vec::new);
                if by_dim.len() <= dim {
                    by_dim.resize(dim+1, 0);
                }
                by_dim[dim] += 1;
            }
        }
    };
    flag_complex::for_each_cell(&normalized_graph, &mut count_for_edges,0, normalized_graph.nnodes());

    //  ... and compute the maximums by dimension
    let mut max_by_dim = vec![];
    for by_dim in ncliques_by_edge_and_dim.values() {
        if max_by_dim.len() < by_dim.len() {
            max_by_dim.resize(by_dim.len(), 0);
        }
        for (m, x) in max_by_dim.iter_mut().zip(by_dim) {
            *m = max(*m, *x);
        }
    }

    return (edge_neighborhood, max_by_dim);
}
