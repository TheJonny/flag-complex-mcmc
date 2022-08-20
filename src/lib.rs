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

use ::flag_complex::*;
use ::flag_complex::prelude::*;


//type Graph = BoolMatrixGraph;

#[derive(Debug, Serialize, Deserialize)]
pub struct EdgeInfo {
    nbhd: Vec<Node>,
    ncliques_in_nbhd: usize,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct State {
    pub cliques_by_order: Vec<Vec<Vec<Node>>>,
    pub edge_neighborhood: HashMap<Edge, Vec<Node>>,
    pub graph: Graph,
    pub flag_count: Vec<usize>,
}

impl State {
    pub fn new(graph: Graph) -> Self {

        //println!("undirected maximal cliques");
        let cliques = graph.compute_maximal_cliques();
        let mut cliques_by_order = vec![];
        for c in cliques.clone() {
            let clique_order = c.len();
            if clique_order > cliques_by_order.len() {
                cliques_by_order.resize(clique_order, vec![]);
            }
            cliques_by_order[clique_order-1].push(c);
        }
        //println!("initial flagser");
        let flag_count = graph.flagser_count();

        println!("computing edge neighborhoods");
        let edge_neighborhood = compute_edge_neighborhoods(&graph);
        
        //flag_count_max.push(10); TODO: ADD SOMETHING LIKE THIS

        State { graph, cliques_by_order, flag_count, edge_neighborhood}
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
#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct Bounds {
    pub flag_count_min: Vec<usize>,
    pub flag_count_max: Vec<usize>,
}
impl Bounds {
    pub fn calculate(initial: &State) -> Self {
        // to determine the acceptance limits, we need to know, how much damage flipping any single edge
        // can do (in each dimension).
        // therefore we calculate the number of adjacent (not necessarly maximal) cliques...

        let undirected_edges = initial.graph.undirected_edges();
        let mut normalized_graph = Graph::new_disconnected(initial.graph.nnodes());
        for &[a,b] in &undirected_edges {
            normalized_graph.add_edge(a, b);
        }

        let mut nqliques_by_edge_and_dim = HashMap::<Edge, Vec<usize>>::with_capacity(undirected_edges.len());

        let mut count_for_edges = |simplex: &[Node]| {
            let dim = simplex.len() - 1;
            // iterate over simplex's edges
            for (i,&b) in simplex.iter().enumerate() {
                for &a in &simplex[0..i] {
                    assert!(a>b);
                    let by_dim = nqliques_by_edge_and_dim.entry([a,b]).or_insert_with(Vec::new);
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
        for by_dim in nqliques_by_edge_and_dim.values() {
            if max_by_dim.len() < by_dim.len() {
                max_by_dim.resize(by_dim.len(), 0);
            }
            for (m, x) in max_by_dim.iter_mut().zip(by_dim) {
                *m = max(*m, *x);
            }
        }
        dbg!(&max_by_dim);

        let additional_relax = 1.05; 
        let mut flag_count_max: Vec<usize> = vec![];
        let mut flag_count_min: Vec<usize> = vec![];
        let relax_de = crate::util::calc_relax_de(&initial.flag_count);
        dbg!(&relax_de);
        for d in 0..initial.flag_count.len() {
            let relax : usize = (std::cmp::max(max_by_dim[d]*2, relax_de[d]))/2;
            flag_count_max.push(((initial.flag_count[d] + 2*relax) as f64 * additional_relax) as usize);
            flag_count_min.push(((initial.flag_count[d] - relax) as f64 / additional_relax) as usize);
        }
        //flag_count_max.push(10); TODO: ADD SOMETHING LIKE THIS
        println!("We have {:?},\n lower limit {:?},\n upper limit {:?}\n", &initial.flag_count, &flag_count_min, &flag_count_max);


        let nchange_dims = max_by_dim.len().checked_sub(2).expect("there should be at least one edge!");

        //for (m, f) in zip_longest(max_by_dim.iter(), 
        // TODO: flag_count_max/min abhängig von max_by_dim
        Bounds { flag_count_min, flag_count_max}
    }
    pub fn check(&self, state: &State) -> bool {
        return all_le(&self.flag_count_min, &state.flag_count, &0)
            && all_le(&state.flag_count, &self.flag_count_max, &0);
    }
}

#[derive(Debug, Serialize, Deserialize)]
pub struct MCMCSampler<R: Rng> {
    // Variable State
    pub rng: R,
    pub state: State,

    // Settings
    pub bounds: Bounds,
    pub move_distribution: WeightedIndex<f64>,
    pub clique_order_distribution: WeightedIndex<f64>,
    pub sample_distance: usize,

    // Metrics
    pub sampled: usize,
    pub accepted: usize,
}

impl<R: Rng> MCMCSampler<R> {
    pub fn next(&mut self) -> &State{
        for _ in 0..self.sample_distance {
            let t = Transition::random_move(&self.state, &mut self.rng, &self.move_distribution, &self.clique_order_distribution);
            let counters = self.state.apply_transition(&t);
            self.sampled += 1;
            if self.bounds.check(&self.state) {
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
    pub fn random_move<R: Rng>(state: &State, rng: &mut R, move_distribution: &WeightedIndex<f64>, clique_order_distribution: &WeightedIndex<f64>) -> Self {
        let potential_moves = [Transition::single_edge_flip_wrap, Transition::double_edge_move_wrap,
                                Transition::clique_permute, Transition::clique_swap];
        let random_move = potential_moves[move_distribution.sample(rng)];
        return random_move(state, rng, clique_order_distribution);
    }

    pub fn clique_permute<R: Rng>(state: &State, rng: &mut R, clique_order_distribution: &WeightedIndex<f64>) -> Self {
        let cliques_of_fixed_order = &state.cliques_by_order[clique_order_distribution.sample(rng)];
        let cl = cliques_of_fixed_order.choose(rng).unwrap();

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

    pub fn clique_swap<R: Rng>(state: &State, rng: &mut R, clique_order_distribution: &WeightedIndex<f64>) -> Self {
        let cliques_of_fixed_order = &state.cliques_by_order[clique_order_distribution.sample(rng)];
        let m1 = cliques_of_fixed_order.choose(rng).unwrap();
        let m2 = cliques_of_fixed_order.choose(rng).unwrap();
        
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

    pub fn single_edge_flip<R: Rng>(state: &State, rng: &mut R) -> Self {
        if let Some([from, to]) = state.graph.sample_edge(rng) {
            if !state.graph.has_edge(to, from) { // its a single edge
                return Transition{change_edges: vec![([from,to], false), ([to,from], true)]};
            }
        }
        return Transition{change_edges: vec![]};
    }
    fn single_edge_flip_wrap<R: Rng>(state: &State, rng: &mut R, _: &WeightedIndex<f64>) -> Self {
        Self::single_edge_flip(state, rng)
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
    fn double_edge_move_wrap<R: Rng>(state: &State, rng: &mut R, _: &WeightedIndex<f64>) -> Self {
        Self::double_edge_move(state, rng)
    }
}
/// for every edge, this gathers the nodes that are connected to both ends.
fn compute_edge_neighborhoods(graph: &Graph)-> HashMap<Edge, Vec<Node>>{
    let mut undirected_adj_lists = vec![vec![]; graph.nnodes()];
    for [a,b] in graph.edges() {
        undirected_adj_lists[a as usize].push(b);
        undirected_adj_lists[b as usize].push(a);
    }
    for v in &mut undirected_adj_lists {
        v.sort_unstable();
    }

    let undirected_edges = graph.undirected_edges();
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
    
    return edge_neighborhood;
}
