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
        let ncells_by_dim = g.flagser_count();
        assert_eq!(&ncells_by_dim, &[7, 12, 6, 1]);
    }
}

pub mod examples {
    use crate::graph::*;
    use crate::Graph;
    pub fn gengraph() -> Graph {
        let mut g = Graph::new_disconnected(7);

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

    pub fn gengraph2() -> Graph {
        // simplest of all simplices of dimension 2
        let mut g = Graph::new_disconnected(3);
        g.add_edge(0,1);
        g.add_edge(0,2);
        g.add_edge(1,2);
        return g
    }
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
}
