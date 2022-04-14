// Authors: "Jonathan Krebs and Florian Unger"
use directed_scm::*;
use rand::SeedableRng;
use rand_xoshiro::Xoshiro256StarStar;

use std::cmp::{min, max};

use flag_complex::prelude::*;

use clap::Parser;


/// MCMC sampler for flag complexes of a directed graph
#[derive(Parser, Debug)]
#[clap(version, about, long_about = None)]
struct Args{
    
    /// random seed 
    #[clap(short, long, default_value_t = 0)]
    seed: u64,


    /// edge probability
    #[clap(short)]
    p: f64,

    #[clap(short='n', long)]
    nnodes: usize,

    #[clap(short='l', long)]
    iteration_limit: usize
}

fn main() {
    let args = Args::parse();
    
    let nnodes = 100;
    let nedges = 1000;
    let mut rng = Xoshiro256StarStar::seed_from_u64(args.seed);
    loop {
        let g = flag_complex::Graph::gen_seo_er(nnodes, nedges, &mut rng);
        let all_clique_count = count_all_cliques(&g);
        let st = State::new(g.clone());
        let mut success = false;
        for _ in 0..args.iteration_limit {
            let t = Transition::single_edge_flip(&st, &mut rng);
            let (pre, post) = st.apply_transition(&t);
            let mut accept = true;
            for(&i,&j) in pre.iter().zip(post.iter()) { if i > j { accept = false; break; } }
            if !accept {
                st.revert_transition(&t, &(pre, post));
            }
            if st.flag_count == all_clique_count {
                success = true;
                break;
            }
        }
        if ! success {
            io::save_flag_file(&format!("counterexample_seo_greedy_{seed}_start.flag", seed=args.seed), &g);
            io::save_flag_file(&format!("counterexample_seo_greedy_{seed}_bad.flag", seed=args.seed), &st.graph);
            break;
        }
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