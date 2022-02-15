// Authors: "Jonathan Krebs and Florian Unger"
use directed_scm::*;
use rand::SeedableRng;
use rand_xoshiro::Xoshiro256StarStar;
use std::fs;

use clap::Parser;


/// MCMC sampler for flag complexes of a directed graph
#[derive(Parser, Debug)]
#[clap(version, about, long_about = None)]
struct Args{
    /// label
    #[clap(short, long)]
    label: String,

    /// flag input file location
    #[clap(short, long)]
    input: String,
    
    /// random seed 
    #[clap(short, long, default_value_t = 0)]
    seed: u64,

    /// number of samples to draw
    #[clap(short, long, default_value_t = 1000)]
    number_of_samples: usize,

    /// burn_in: number of steps to burn in MCMC sampler
    #[clap(long, default_value_t = 5000)]
    burn_in: usize,

    /// sample_distance: number of steps between too samples
    #[clap(long, default_value_t = 4000)]
    sample_distance: usize,

    /// continue: look for a serde file and continue from there
    #[clap(short, long, default_value = "")]
    continue_from: String,

    /// samples_store_dir: directory where to store hdf5 files with samples
    #[clap(long, default_value = "./samples/")]
    samples_store_dir: String,

    /// state_store_dir: directory where to store hdf5 files with samples
    #[clap(long, default_value = "./state/")]
    state_store_dir: String,
}

fn main() {
    let args = Args::parse();
    fs::create_dir_all(&args.state_store_dir).expect("could not create state storage directory");
    fs::create_dir_all(&args.samples_store_dir).expect("could not create samples storage directory");
    
    let (sample_index_start, mut sampler) = if !args.continue_from.is_empty() {
        io::load_state(&args.continue_from).expect("unable to load state")
    } else {
        let g = io::read_flag_file(&args.input);
        io::new_hdf_file(&args.samples_store_dir, &args.label, args.seed).unwrap();
        let st = State::new(g);
        println!("we have the following number of maximal k-cliques {:?}", st.cliques_by_order.iter().map(|cs| cs.len()).collect::<Vec<usize>>());
        let rng = Xoshiro256StarStar::seed_from_u64(args.seed);
        let move_distribution = rand::distributions::WeightedIndex::new([0.1, 0.1, 0.06, 0.2]).unwrap();
        let mut sampler = MCMCSampler{state: st, move_distribution, burn_in: args.burn_in, sample_distance: args.sample_distance, accepted: 0, sampled: 0, rng};
        sampler.burn_in();
        (0, sampler)
    };
    let sample_index_end = sample_index_start + args.number_of_samples;
    for i in sample_index_start..sample_index_end {
        let s = sampler.next();
        io::save_to_hdf(&args.samples_store_dir, &args.label, args.seed, i, &s.graph, &s.flag_count).unwrap();

        println!("flag count: {:?}", s.flag_count);
        drop(s);
        dbg!(sampler.acceptance_ratio());
        if (i % 1000) == 0 {
            println!("saving state in step {i}");
            io::save_state(&format!("{state_store_dir}/sampler-{l}-{s:03}.state", state_store_dir = args.state_store_dir, l=args.label, s=args.seed), i, &sampler).unwrap();
        }
    }
    //println!("whole graph flag count: {:?}", &sampler.state.graph.flagser_count());

    io::save_state(&format!("{state_store_dir}/sampler-{l}-{s:03}.state", state_store_dir = args.state_store_dir, l=args.label, s=args.seed), sample_index_end, &sampler).unwrap();
}

