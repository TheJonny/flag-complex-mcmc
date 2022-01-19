use std::marker::PhantomData;


use rand;
use rand::prelude::*;

pub mod flagser;
pub mod io;

type Node = u32;
type CliqueId = usize;

#[derive(Clone)]
pub struct DirectedGraph {
    adjmat: Vec<bool>,
    nnodes: usize,
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
        flagser::count_unweighted(self.nnodes, &self.edges())
    }
}

#[derive(Clone, Debug)]
pub struct State {
    pub cliques: Vec<Vec<Node>>,
    pub which_cliques: Vec<Vec<CliqueId>>,
    pub graph: DirectedGraph,
    pub flag_count: Vec<usize>,
    pub flag_count_min: Vec<usize>,
    pub flag_count_max: Vec<usize>,
}

#[derive(Clone, Debug)]
pub struct Shuffling(CliqueId, Vec<usize>);

impl State {
    pub fn new(graph: DirectedGraph) -> Self {
        let cliques = graph.compute_maximal_cliques();
        let mut which_cliques = vec![Vec::new(); graph.nnodes];
        for (c,cid) in cliques.iter().zip(0..) {
            for &v in c  {
                which_cliques[v as usize].push(cid);
            }
        }
        let flag_count = graph.flagser_count();
        let bandwidth = 1.10;
        let mut flag_count_max: Vec<_> = flag_count.iter().map(|x| (*x as f64 * bandwidth + 10.) as usize).collect();
        flag_count_max.push(2); // some higher dimensional simplices should be ok
        let flag_count_min: Vec<_> = flag_count.iter().map(|x| (*x as f64 / bandwidth - 10.) as usize).collect();
        println!("We have {:?},\n lower limit {:?},\n upper limit {:?}\n", &flag_count, &flag_count_min, &flag_count_max);

        State { graph, cliques, which_cliques, flag_count, flag_count_min, flag_count_max}
    }
    pub fn sample_shuffling<R: Rng>(&self, rng: &mut R) -> Shuffling {
        let cid = rng.gen_range(0..self.cliques.len() as CliqueId);
        let n = self.cliques[cid].len();
        let mut perm = (0..n).collect::<Vec<usize>>();
        perm.shuffle(rng);
        Shuffling(cid, perm)
    }
    pub fn apply_shuffling(&mut self, s: &Shuffling) {
        let Shuffling(cid, perm) = s;
        let cid = *cid;
        let cn = self.clique_neighborhood(cid);
        let pre = cn.flagser_count();
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

        // TODO: vertices der neighbourhood müssen nicht neu gesucht werden, kann man speichern
        let post = self.clique_neighborhood(cid).flagser_count();
        if post.len() > self.flag_count.len() {
            self.flag_count.resize(post.len(), 0);
        }
        for (p,s) in post.iter().zip(self.flag_count.iter_mut()) {
            *s += *p;
        }

    }

    /// returns the subgraph affected by shuffling vertices/edges in clique cid.
    /// it consists of all cliques, that share at least one edge with cid.
    pub fn clique_neighborhood(&self, cid: CliqueId) -> DirectedGraph {

        // which cliques are affected, when the edges in the clique "cid" are shuffled?
        // therefore count how many vertices they have in common with cid
        let mut affected_cids = vec![cid];
        //let mut neigh_cliq_count = vec![0; self.cliques.len()];
        let mut neigh_cliqs = std::collections::HashSet::new();
        for &v in &self.cliques[cid as usize] {
            for &vcid in &self.which_cliques[v as usize] {
                if vcid != cid {
                    if neigh_cliqs.contains(&vcid) {
                        affected_cids.push(vcid);
                    }
                    neigh_cliqs.insert(vcid);
                }
            }
        }

        affected_cids.sort_unstable();
        affected_cids.dedup();
        let mut affected_vertices = Vec::new();
        for acid in affected_cids {
            affected_vertices.extend_from_slice(&self.cliques[acid as usize]);
        }
        affected_vertices.sort_unstable();
        affected_vertices.dedup();

        return self.graph.subgraph(&affected_vertices)
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
    pub fn next(&mut self) -> State{
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
        return self.state.clone();
    }
    pub fn acceptance_ratio(&self) -> f64 {
        return self.accepted as f64 / self.sampled as f64;
    }
}
