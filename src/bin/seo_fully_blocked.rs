// Authors: "Jonathan Krebs and Florian Unger"
use directed_scm::*;
use rand::SeedableRng;
use rand_xoshiro::Xoshiro256StarStar;

use std::cmp::{min, max};
use ::flag_complex::*;
use ::flag_complex::Graph;
use ::flag_complex::prelude::*;

use clap::Parser;


#[derive(Parser, Debug)]
#[clap(version, about, long_about = None)]
struct Args{
    /// flag input file location
    #[clap(short, long, default_value = "")]
    input: String,
    

    /// random seed 
    #[clap(short, long, default_value_t = 0)]
    seed: u64,

    /// edge probability
    #[clap(short, default_value_t = 1.0)]
    p: f64,

    #[clap(short='n', long, default_value_t = 0)]
    nnodes: u32,

    #[clap(short='d', long)]
    debug: bool,
}


fn check_blockedness(g: &mut Graph) -> bool {
    g.edges().sort_unstable();
    for [i,j] in g.edges() {
        let mut blocked = false;
        for v in 0..g.nnodes() as u32 {
            if g.has_edge(i,v) & g.has_edge(v,j) {blocked = true; break}
        }
        if !blocked {return false}
    }
    return true

}

fn main() {
    let args = Args::parse();
    let mut seed = args.seed;
    loop {
        let mut rng = Xoshiro256StarStar::seed_from_u64(seed);
        let mut g = Graph::gen_seo_er(args.nnodes, args.p, &mut rng);

        if args.debug {
            io::save_flag_file(&format!("debug_seo_fully_blocked_N{N:02}_s{seed:09}_tart.flag", seed=seed, N=args.nnodes), &g);
        }

        let success = check_blockedness(&mut g);

        if success {
            //println!("Found one with seed {seed}");
            let n_cliques = count_all_cliques(&g)[2];
            let n_simplices = count_cells(&g)[2];
            let cycle_ratio = (n_cliques-n_simplices) as f32 / n_cliques as f32;
            if cycle_ratio < 0.25 {
                println!("With {n_cliques} 3-cliques and {n_simplices} 2-simplices, we have {cycle_ratio}% 3-cycles ");
                io::save_flag_file(&format!("example_seo_fully_blocked_only_once_N{N:02}_s{seed:09}_start.flag", seed=seed, N=args.nnodes), &g);
            }
            //std::process::exit(1);
        } else {
            //println!("no example for seed {seed}");
        }
        if !args.input.is_empty() {
            std::process::exit(0);
        }
        seed += 1000;
    }
}

fn count_all_cliques(graph: &Graph) -> Vec<usize>{
    let mut undirected_edges = graph.edges();
    for e in &mut undirected_edges {
        let a = max(e[0], e[1]);
        let b = min(e[0], e[1]);
        *e = [a, b];
    }
    undirected_edges.sort_unstable();
    undirected_edges.dedup();
    let mut normalized_graph = Graph::new_disconnected(graph.nnodes());
    for &[a,b] in &undirected_edges {
        normalized_graph.add_edge(a, b);
    }
    flag_complex::count_cells(&normalized_graph)
}
