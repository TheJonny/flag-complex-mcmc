// Authors: "Jonathan Krebs and Florian Unger"
use directed_scm::*;
use rand::SeedableRng;
use rand_xoshiro::Xoshiro256StarStar;

use std::cmp::{min, max};

use flag_complex::{prelude::*, Graph};
use flag_complex::{Node, Edge};

use clap::Parser;


/// MCMC sampler for flag complexes of a directed graph
#[derive(Parser, Debug, Clone)]
#[clap(version, about, long_about = None)]
struct Args{
    /// edge probability
    #[clap(short)]
    p: f64,

    #[clap(short='n', long)]
    nnodes: u32,

    #[clap(short='l', long)]
    iteration_limit: usize,
}

fn main() {
    let args = Args::parse();
 
    let cxs = Box::new(std::sync::Mutex::<std::collections::HashSet::<_>>::default());
    let cxs: &'static _ = Box::leak(cxs);

    let mut threads = Vec::new();
    for i in 0..8 {
        let args = args.clone();
        let t = std::thread::spawn(move || doit(args, cxs));
        threads.push(t);
    }
    for t in threads {
        t.join().unwrap();
    }
}

fn  doit(args: Args, cxs: &std::sync::Mutex::<std::collections::HashSet::<Vec<Edge>>>) {
    loop {
        let mut rng = rand::thread_rng();

        let cx = loop {
            let g = flag_complex::Graph::gen_seo_er(args.nnodes, args.p, &mut rng);
            let all_clique_count = count_all_cliques(&g);
            let mut st = State::new(g.clone());
            let mut success = false;
            for i in 0..args.iteration_limit {
                let t = Transition::single_edge_flip(&st, &mut rng);
                // for exact mode:
                let mut cells_before = Vec::new();
                let mut cells_after = std::collections::HashSet::new();
                let nei = st.edgeset_neighborhood(&t.change_edges.iter().map(|&([a,b], _)| [max(a,b), min(a,b)]).collect::<Vec<Edge>>());
                let localg = flag_complex::Graph::subgraph(&st.graph, &nei);
                flag_complex::for_each_cell(&localg, &mut |cell: &[Node]|{
                    let mut cv = Vec::from(cell);
                    cv.sort();
                    cells_before.push(cv);
                }, 2,2);
                let (pre, post) = st.apply_transition(&t);
                let mut accept = true;
                if post.len() < pre.len() { accept = false; }
                for(&i,&j) in pre.iter().zip(post.iter()) { if i > j { accept = false; break; } }
                if accept {
                    //dbg!(&cells_before);
                    let localg = flag_complex::Graph::subgraph(&st.graph, &nei);
                    flag_complex::for_each_cell(&localg, &mut |cell: &[Node]|{
                        let mut cv = Vec::from(cell);
                        cv.sort();
                        cells_after.insert(cv);
                    }, 2,2);
                    for c in cells_before {
                        if ! cells_after.contains(&c) {
                            accept = false;
                            //          dbg!(&cv);
                        }
                    }
                }
                if !accept {
                    st.revert_transition(&t, &(pre, post));
                }
                if st.flag_count == all_clique_count {
                    success = true;
                    break;
                }
            }
            if ! success {
                break st.graph.clone();
            }
        };

        let mut e = cx.edges();
        e.sort();
        let mut guard = cxs.lock().unwrap();
        let is_new = guard.insert(e);
        if is_new {
            println!("count: {}", guard.len());
        }
        drop(guard);
    }
}

fn count_all_cliques(graph: &flag_complex::Graph) -> Vec<usize>{

    let mut undirected_edges = graph.edges();
    for e in &mut undirected_edges {
        let a = max(e[0], e[1]);
        let b = min(e[0], e[1]);
        *e = [a, b];
    }
    undirected_edges.sort_unstable();
    undirected_edges.dedup();
    let mut normalized_graph = flag_complex::Graph::new_disconnected(graph.nnodes());
    for &[a,b] in &undirected_edges {
        normalized_graph.add_edge(a, b);
    }

    flag_complex::count_cells(&normalized_graph)
}
