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
}

fn flip(e: Edge) -> Transition {
    Transition { change_edges: vec![(e,false), ([e[1], e[0]], true)] }
}

fn rec(state: &mut State, remaining_edges: &mut Vec<Edge>, target: usize) -> bool {
    if state.flag_count.get(2).copied().unwrap_or(0) == target {
        return true;
    }
    if remaining_edges.is_empty() {
        //println!("no edges remaining");
        return false;
    }

    let mut edges_with_change = remaining_edges.iter().enumerate().map(|(i, &e)| {
        let t = flip(e);
        let (pre, post) = state.apply_transition(&t);
        let balance = (post.get(2).copied().unwrap_or(0) as i32) - 
                      (pre.get(2).copied().unwrap_or(0) as i32);
        state.revert_transition(&t, &(pre, post));

        (balance, e, i)
    }).collect::<Vec<_>>();

    // sort descending by balance:
    //  first element will be the maximum
    edges_with_change.sort_by(|a,b| (b.0).cmp(&a.0));

    if edges_with_change[0].0 < 0 { // no increasing move -> bail out
        //println!("no good moves");
        return false;
    }
    
    let end_index =
        if edges_with_change[0].0 > 0 {
            1
        } else {
            let mut i = 0;
            while i < edges_with_change.len() && edges_with_change[i].0 == 0 {
                i += 1;
            }
            i
        };

    for i in 0 .. end_index {
        let e = edges_with_change[i].1;
        let ei = edges_with_change[i].2;
        let t = flip(e);
        //println!("Flip: {e:?}");
        let (pre, post) = state.apply_transition(&t);
        
        assert!(remaining_edges[ei] == e);
        remaining_edges.swap_remove(edges_with_change[i].2);

        if rec(state, remaining_edges, target) {
            return true;
        }
        //println!("rollback");

        state.revert_transition(&t, &(pre, post));
        remaining_edges.push(e);
        let n = remaining_edges.len();
        remaining_edges.swap(ei, n-1);
        assert!(remaining_edges[ei] == e);
    }

    println!("possibilities exhausted");
    return false;
}

fn main() {
    let args = Args::parse();
    let mut rng = Xoshiro256StarStar::seed_from_u64(args.seed);
    loop {
        let g = if args.input.is_empty() {
            let mut g = Graph::gen_seo_er(args.nnodes, args.p, &mut rng);
            while (count_all_cliques(&g).len() < 3) {                   // prevents graphs without 3-cliques
                g = Graph::gen_seo_er(args.nnodes, args.p, &mut rng);
            }
            g
        } else {
            io::read_flag_file(&args.input)
        };
        let mut edges = g.edges();
        let target = count_all_cliques(&g)[2];

        let mut st = State::new(g.clone());

        let success = rec(&mut st, &mut edges, target);
        if !success {
            println!("FAIL");
            io::save_flag_file(&format!("counterexample_seo_flip_only_once_{seed:02}_N_{N:02}_start.flag", seed=args.seed, N=args.nnodes), &g);
            //io::save_flag_file(&format!("counterexample_seo_flip_only_once_{seed:02}_bad.flag", seed=args.seed), &st.graph);
            std::process::exit(1);
        } else {
            println!("worked this time");
        }
        if !args.input.is_empty() {
            std::process::exit(1);
        }
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
