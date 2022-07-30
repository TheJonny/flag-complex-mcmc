// Authors: "Jonathan Krebs and Florian Unger"
use directed_scm::*;
use rand::SeedableRng;
use rand_xoshiro::Xoshiro256StarStar;

use std::cmp::{min, max};

use flag_complex::prelude::*;
use flag_complex::{Node, Edge};

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
    nnodes: u32,

    #[clap(short='l', long)]
    iteration_limit: usize,

    /// if each simplex should be tracked and not just the numbers
    #[clap(short='x', long)]
    exact: bool,
}

fn main() {
    let args = Args::parse();
    
    let mut rng = Xoshiro256StarStar::seed_from_u64(args.seed);
    loop {
        let g = flag_complex::Graph::gen_seo_er(args.nnodes, args.p, &mut rng);
        let all_clique_count = count_all_cliques(&g);
        println!("all clique count: {:?}", all_clique_count);
        let mut st = State::new(g.clone());
        let mut success = false;
        for i in 0..args.iteration_limit {
            let t = Transition::single_edge_flip(&st, &mut rng);
            // for exact mode:
            let mut cells_before = Vec::new();
            let mut cells_after = std::collections::HashSet::new();
            let nei = st.edgeset_neighborhood(&t.change_edges.iter().map(|&([a,b], _)| [max(a,b), min(a,b)]).collect::<Vec<Edge>>());
            if args.exact {
                let localg = flag_complex::Graph::subgraph(&st.graph, &nei);
                flag_complex::for_each_cell(&localg, &mut |cell: &[Node]|{
                    let mut cv = Vec::from(cell);
                    cv.sort();
                    cells_before.push(cv);
                }, 2,2);
            }
            let (pre, post) = st.apply_transition(&t);
            let mut accept = true;
            if post.len() < pre.len() { accept = false; }
            if post[2] <= pre[2] {accept = false;}
            //for(&i,&j) in pre.iter().zip(post.iter()) { if i > j { accept = false; break; } }

            if accept {
                //dbg!(&cells_before);
                if args.exact {
                    let localg = flag_complex::Graph::subgraph(&st.graph, &nei);
                    flag_complex::for_each_cell(&localg, &mut |cell: &[Node]|{
                        let mut cv = Vec::from(cell);
                        cv.sort();
                        cells_after.insert(cv);
                    }, 2,2);
                }
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
            if i % 1000 == 999 {
                println!("{} steps, flag_count: {:?} of {:?}", i+1, st.flag_count, all_clique_count);
            }
            if st.flag_count == all_clique_count {
                success = true;
                println!("reached maximum after {} steps", i+1);
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
