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

fn flip(e: Edge) -> Transition {
    Transition { change_edges: vec![(e,false), ([e[1], e[0]], true)] }
}

fn calc_degree(nnodes: usize, edges: &Vec<Edge>) -> (Vec<u32>, Vec<u32>) {
    let mut outdeg = vec![0; nnodes];
    let mut indeg = vec![0; nnodes];
    for &[i,j] in edges {
        outdeg[i as usize] += 1;
        indeg[j as usize] += 1;
    }
    return (outdeg, indeg)
}

fn calc_order_and_sc_increase(state: &mut State, edges: &Vec<Edge>) -> Vec<(isize, isize, Edge, usize)> {
    let (outdeg, indeg) = calc_degree(state.graph.nnodes(), &state.graph.edges());
    println!("out {outdeg:?} in {indeg:?}");
    let edges_with_change = edges.iter().enumerate().map(|(ind, &[i,j])| {
        //order change
        let mut order_increase = 0;
        if indeg[i as usize] > outdeg[i as usize] {order_increase += 2;}
        else if indeg[i as usize] == outdeg[i as usize] {order_increase += 1;}
        else {order_increase += -2;}
        if indeg[j as usize] < outdeg[j as usize] {order_increase += 2;}
        else if indeg[j as usize] == outdeg[j as usize] {order_increase += 1;}
        else {order_increase += -2;}
        
        //simplex count change
        let t = flip([i,j]);
        let (pre, post) = state.apply_transition(&t);
        let sc_increase = (post.get(2).copied().unwrap_or(0)) as isize - 
                          (pre.get(2).copied().unwrap_or(0)) as isize;
        state.revert_transition(&t, &(pre, post));

        (order_increase, sc_increase, [i,j], ind)
    }).collect::<Vec<_>>();
    return edges_with_change
}

fn identify_dropable_vertices(state: &State) -> Vec<Node> {
    let edges = state.graph.edges();
    let nnodes = state.graph.nnodes();
    let (outdeg, indeg) = calc_degree(nnodes, &edges);
    // dropable: indeg < 3 OR outdeg < 3. If neither is, not dropable.
    // If both are 0, it was already dropped.
    let dropable_vertices = (0..nnodes as Node).filter(|&i|
                    ((outdeg[i as usize] < 3) | (indeg[i as usize] < 3))
                    & (indeg[i as usize] > 0 | outdeg[i as usize]) )
        .collect::<Vec<Node>>();
    return dropable_vertices
}

fn rec(state: &mut State, remaining_edges: &mut Vec<Edge>, target: usize) -> bool {
    println!("current number of 2-simps is {} while target is {target}", state.flag_count.get(2).copied().unwrap_or(0));
    if state.flag_count.get(2).copied().unwrap_or(0) == target {
        return true;
    }

    //  DROPPING VERTICES
    let dropable_vertices = identify_dropable_vertices(&state);
    println!("vertices to drop: {dropable_vertices:?}");
    if dropable_vertices.len() > 0 {
        let mut new_graph = state.graph.clone();
        for [i,j] in new_graph.edges() {
            if dropable_vertices.binary_search(&i).is_ok() || dropable_vertices.binary_search(&j).is_ok() {
                new_graph.remove_edge(i,j);
            }
        }
        // update desired simplex count and actual simplices
        let mut new_state = State::new(new_graph.clone());
        let new_target = count_all_cliques(&new_state.graph).get(2).copied().unwrap_or(0);

        // update remaining edges
        let mut new_remaining_edges = remaining_edges.iter().cloned().filter(|[i,j]| !(dropable_vertices.binary_search(&i).is_ok() || dropable_vertices.binary_search(&j).is_ok())).collect();

        if rec(&mut new_state, &mut new_remaining_edges, new_target) {
            return true;
        }
    }

    // Looking which edges to flip 
    let mut edges_with_change = calc_order_and_sc_increase(state, remaining_edges);
    
    // sort descending by order balance: first element will be the maximum
    edges_with_change.sort_by(|a,b| (b.0).cmp(&a.0));
    println!("moves: {:?}", edges_with_change);

    /*
    for &(oi, sci, e, ei) in &edges_with_change {
        if (oi > 0) & (sci < 0) {
            println!("edge {e:?} has positive order increase {oi}, but positive simplex count increase {sci}");
            println!("graph: {:?}", state.graph);
        } else if (sci > 0) & (oi < 0) {
            println!("edge {e:?} has negative order increase {oi}, but positive simplex count increase {sci}");
            println!("graph: {:?}", state.graph);
        }
    }
    */
    for (oi, sci, e, ei) in edges_with_change {
        if (oi >= 0) & (sci >= 0) {
            let t = flip(e);
            println!("Flip: {e:?}");
            let (pre, post) = state.apply_transition(&t);
            
            assert!(remaining_edges[ei] == e);
            remaining_edges.swap_remove(ei);


            if rec(state, remaining_edges, target) {
                return true;
            }
            println!("rollback");

            state.revert_transition(&t, &(pre, post));
            remaining_edges.push(e);
            let n = remaining_edges.len();
            remaining_edges.swap(ei, n-1);
            assert!(remaining_edges[ei] == e);
        }
    }
    println!("possibilities exhausted");
    return false;
}

fn main() {
    let args = Args::parse();
    let mut seed = args.seed;
    loop {
        let mut rng = Xoshiro256StarStar::seed_from_u64(seed);
        let g = if args.input.is_empty() {
            let mut g = Graph::gen_seo_er(args.nnodes, args.p, &mut rng);
            while count_all_cliques(&g).len() < 3 {                   // prevents graphs without 3-cliques
                g = Graph::gen_seo_er(args.nnodes, args.p, &mut rng);
            }
            g
        } else {
            io::read_flag_file(&args.input)
        };
        if args.debug {
            io::save_flag_file(&format!("debug_seo_flip_only_once_{seed:02}_N_{N:02}_start.flag", seed=seed, N=args.nnodes), &g);
        }
        let mut edges = g.edges();
        let target = count_all_cliques(&g)[2];

        let mut st = State::new(g.clone());

        let success = rec(&mut st, &mut edges, target);
        if !success {
            println!("FAIL on seed {seed}");
            io::save_flag_file(&format!("counterexample_seo_flip_only_once_{seed:02}_N_{N:02}_start.flag", seed=seed, N=args.nnodes), &g);
            std::process::exit(1);
        } else {
            println!("worked for seed {seed}");
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
