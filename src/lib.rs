use std::marker::PhantomData;


use std::cmp::{max, min};
use rand;
use rand::prelude::*;

pub mod io;

mod util;

pub mod graph;
use graph::*;

mod flagser;

/// undirected edge pair to index in triangular adjacency matrix
/// ```
/// use directed_scm::*;
/// assert_eq!(edge_id(1, 0), 0);
/// assert_eq!(edge_id(2, 0), 1);
/// assert_eq!(edge_id(2, 1), 2);
/// assert_eq!(edge_id(3, 0), 3);
/// assert_eq!(edge_id(3, 1), 4);
/// assert_eq!(edge_id(3, 2), 5);
/// ```

#[derive(Clone, Debug)]
pub struct State {
    pub cliques: Vec<Vec<Node>>,
    pub which_cliques: Vec<Vec<CliqueId>>,
    pub edge_neighborhood: Vec<Vec<Node>>,
    pub graph: DirectedGraph,
    pub flag_count: Vec<usize>,
    pub flag_count_min: Vec<usize>,
    pub flag_count_max: Vec<usize>,
}

#[derive(Clone, Debug)]
pub struct Shuffling(CliqueId, Vec<usize>);

impl State {
    pub fn new(graph: DirectedGraph) -> Self {

        println!("undirected maximal cliques");
        let cliques = graph.compute_maximal_cliques();
        let mut which_cliques = vec![Vec::new(); graph.nnodes];
        for (c,cid) in cliques.iter().zip(0..) {
            for &v in c  {
                which_cliques[v as usize].push(cid);
            }
        }
        println!("initial flagser");
        let flag_count = graph.flagser_count();
        let bandwidth = 1.10;
        let mut flag_count_max: Vec<_> = flag_count.iter().map(|x| (*x as f64 * bandwidth + 10.) as usize).collect();
        flag_count_max.push(2); // some higher dimensional simplices should be ok
        let flag_count_min: Vec<_> = flag_count.iter().map(|x| (*x as f64 / bandwidth - 10.) as usize).collect();
        println!("We have {:?},\n lower limit {:?},\n upper limit {:?}\n", &flag_count, &flag_count_min, &flag_count_max);

        println!("computing edge neighborhoods");
        let edge_neighborhood = graph.compute_edge_neighborhoods();

        State { graph, cliques, which_cliques, flag_count, flag_count_min, flag_count_max, edge_neighborhood}
    }
    pub fn sample_shuffling<R: Rng>(&self, rng: &mut R) -> Shuffling {
        let cid = rng.gen_range(0..self.cliques.len() as CliqueId);
        let n = self.cliques[cid].len();
        let mut perm = (0..n).collect::<Vec<usize>>();
        perm.shuffle(rng);
        Shuffling(cid, perm)
    }
    /// new.edge(s[i], s[j]) <-> old.edge(i,j)
    pub fn apply_shuffling(&mut self, s: &Shuffling) {
        if true { // toggle between new and old algorithm
            self.apply_transition(&Transition::new_clique_shuffling(self, s.0, &s.1));
            return;
        }
        let Shuffling(cid, perm) = s;
        let cid = *cid;
        let cn = self.clique_neighborhood(cid);
        let pre = self.graph.subgraph(&cn).flagser_count();
        //println!("{:?}", &self.clique_neighborhood(cid));
        //dbg!(cn.nnodes);
        for (p,s) in pre.iter().zip(self.flag_count.iter_mut()) {
            assert!(*s >= *p);
            *s -= *p;
        }

        let cl = &mut self.cliques[cid as usize];

        // set new edge directions:
        // create a new subgraph adjacency matrix as hashmap and then apply it to self.
        {
            use std::collections::HashMap;
            let mut adjmat_new = HashMap::<(Node, Node),  bool>::new();

            // prepare new adjacency matrix...
            for i in 0 .. cl.len() {
                for j in 0 .. cl.len() {
                    adjmat_new.insert((cl[perm[i]], cl[perm[j]]), self.graph.edge(cl[i], cl[j]));
                }
            }
            // ..apply it
            for &from in cl.iter() {
                for &to in cl.iter() {
                    *self.graph.edge_mut(from, to) = adjmat_new[&(from, to)];
                }
            }
        }

        // TODO: vertices der neighborhood müssen nicht neu gesucht werden, kann man speichern
        let post = self.graph.subgraph(&cn).flagser_count();
        if post.len() > self.flag_count.len() {
            self.flag_count.resize(post.len(), 0);
        }
        for (p,s) in post.iter().zip(self.flag_count.iter_mut()) {
            *s += *p;
        }

    }

    /// applies transition, returns the change in simplex counts
    pub fn apply_transition(&mut self, t: &Transition) -> (Vec<usize>, Vec<usize>) {
        let nei = self.edgeset_neighborhood(&t.change_edges.iter().map(|&([a,b], _)| [max(a,b), min(a,b)]).collect::<Vec<Edge>>());
        let pre = DirectedGraph::subgraph(&self.graph, &nei).flagser_count();
        for (p,s) in pre.iter().zip(self.flag_count.iter_mut()) {
            assert!(*s >= *p);
            *s -= *p;
        }
        for &([a,b],add) in &t.change_edges {
            *self.graph.edge_mut(a, b) = add;
        }
        let post = self.graph.subgraph(&nei).flagser_count();
        if post.len() > self.flag_count.len() {
            self.flag_count.resize(post.len(), 0);
        }
        for (p,s) in post.iter().zip(self.flag_count.iter_mut()) {
            *s += *p;
        }

        return (pre, post);
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

impl MarcovState for State {
    fn valid(&self) -> bool {
        return all_le(&self.flag_count_min, &self.flag_count, &0)
            && all_le(&self.flag_count, &self.flag_count_max, &0);
    }
}

#[cfg(test)]
mod tests {
    use crate::*;
    #[test]
    fn it_works() {
        let result = 2 + 2;
        assert_eq!(result, 4);
    }


    #[test]
    fn find_cliques() {
        let g = examples::gengraph();
        // normalize
        let mut cc = g.compute_maximal_cliques();
        for c in &mut cc {
            c.sort();
        }
        cc.sort();
        assert_eq!(cc, vec![
                   vec![0,1,2],
                   vec![0,2,5,6],
                   vec![1,2,3],
                   vec![4],
            ]);

    }

    #[test]
    fn neighborhood() {
        let s = State::new(examples::gengraph());
        assert_eq!(s.edge_neighborhood.len(), 7 * (7-1) / 2);
        assert_eq!(s.edge_neighborhood[edge_id(1,0)], vec![2]);
        assert_eq!(s.edge_neighborhood[edge_id(2,0)],vec![1, 5,6]);
        assert_eq!(s.edge_neighborhood[edge_id(2,1)], vec![0,3]);
        assert_eq!(s.edge_neighborhood[edge_id(3,0)], vec![]); // does not exist
        assert_eq!(s.edge_neighborhood[edge_id(3,1)], vec![2]);
        assert_eq!(s.edge_neighborhood[edge_id(3,2)], vec![1]);
        assert_eq!(s.edge_neighborhood[edge_id(4,0)], vec![]); // 4 is separated
        assert_eq!(s.edge_neighborhood[edge_id(4,1)], vec![]); // 4 is separated
        assert_eq!(s.edge_neighborhood[edge_id(4,2)], vec![]); // 4 is separated
        assert_eq!(s.edge_neighborhood[edge_id(4,3)], vec![]); // 4 is separated
        assert_eq!(s.edge_neighborhood[edge_id(5,0)], vec![2,6]); // 
        assert_eq!(s.edge_neighborhood[edge_id(5,1)], vec![]); // does not exist
        assert_eq!(s.edge_neighborhood[edge_id(5,2)], vec![0,6]);
        assert_eq!(s.edge_neighborhood[edge_id(5,3)], vec![]); // does not exist
        assert_eq!(s.edge_neighborhood[edge_id(5,4)], vec![]); // does not exist
        assert_eq!(s.edge_neighborhood[edge_id(6,0)], vec![2,5]);
        assert_eq!(s.edge_neighborhood[edge_id(6,1)], vec![]); // does not exist
        assert_eq!(s.edge_neighborhood[edge_id(6,2)], vec![0,5]);
        assert_eq!(s.edge_neighborhood[edge_id(6,3)], vec![]); // does not exist
        assert_eq!(s.edge_neighborhood[edge_id(6,4)], vec![]); // 4 i separated
        assert_eq!(s.edge_neighborhood[edge_id(6,5)], vec![0,2]);
        dbg!(&s);
        dbg!(&s.cliques[0]);
        dbg!(s.clique_neighborhood(0));
        dbg!(&s.cliques[1]);
        dbg!(s.clique_neighborhood(1));
    }

    #[test]
    fn flagser_count(){
        let g = examples::gengraph();
        let mut edges = vec![];
        for from in 0..(g.nnodes as Node) {
            for to in 0..(g.nnodes as Node) {
                if g.edge(from, to) {
                    edges.push([from, to]);
                }
            }
        }
        let ncells_by_dim = flagser::count_unweighted(g.nnodes, &edges);
        assert_eq!(&ncells_by_dim, &[7, 12, 6, 1]);
    }
}

pub mod examples {
    use crate::DirectedGraph;
    pub fn gengraph() -> DirectedGraph {
        let mut g = DirectedGraph::new_disconnected(7);

        // 0,1
        g.add_edge(0,1);
        g.add_edge(1,0);

        // 1,2,3
        g.add_edge(1,2);
        g.add_edge(2,3);
        g.add_edge(3,1);

        // 0,2,5,6
        g.add_edge(0,2);
        g.add_edge(0,5);
        g.add_edge(0,6);
        g.add_edge(2,5);
        g.add_edge(2,6);
        g.add_edge(5,6);
        g.add_edge(6, 0);
        
        // this implicitely extended 0,1 to 0,1,2 
        // 0,1,2

        // 4 is left over and maximal
        return g
    }

    pub fn gengraph2() -> DirectedGraph {
        // simplest of all simplices of dimension 2
        let mut g = DirectedGraph::new_disconnected(3);
        g.add_edge(0,1);
        g.add_edge(0,2);
        g.add_edge(1,2);
        return g
    }
}

pub trait Action<State : ?Sized> {
    //type State;
    fn reverse(&self) -> Self;
    fn apply(&self, state: &mut State);
    fn gen<R: Rng>(state: &State, rng: &mut R) -> Self;
}

pub trait MarcovState {
    fn valid(&self) -> bool;
    fn apply<A: Action<Self>>(&mut self, action: &A) {
        action.apply(self);
    }
}

impl Action<State> for Shuffling {
    fn reverse(&self) -> Self{
        let Shuffling(cid,perm) = self;
        let mut inverse = vec![0; perm.len()];
        for i in 0 .. perm.len() {
            inverse[perm[i]] = i;
        }
        return Shuffling(*cid, inverse);
    }

    fn apply(&self, state: &mut State) {
         state.apply_shuffling(self);
    }
    fn gen<R: Rng>(state: &State, rng: &mut R) -> Shuffling {
        state.sample_shuffling(rng)
    }
}

pub struct MCMCSampler<S: MarcovState, R: Rng, A: Action<S>> {
    pub state: S,
    pub burn_in: usize,
    pub sample_distance: usize,
    pub sampled: usize,
    pub accepted: usize,
    pub rng: R,
    pub _a: PhantomData<A>
}

impl<State: MarcovState+Clone, R, A> MCMCSampler<State, R, A>
    where A: Action<State> + std::fmt::Debug, R: Rng
{
    pub fn burn_in(&mut self) {
        while self.burn_in > 0 {
            self.burn_in -= 1;
            let a = A::gen(&self.state, &mut self.rng);
            self.state.apply(&a);
            if ! self.state.valid() {
                let a = a.reverse();
                self.state.apply(&a);
            }
        }
    }
    pub fn next(&mut self) -> &State{
        for _ in 0..self.sample_distance {
            let a = A::gen(&self.state, &mut self.rng);
            //println!("{:?}", &a);
            self.state.apply(&a);
            self.sampled += 1;
            if self.state.valid() {
                self.accepted += 1;
            }
            else {
                let a = a.reverse();
                self.state.apply(&a);
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
}
