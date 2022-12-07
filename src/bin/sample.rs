// Authors: "Jonathan Krebs and Florian Unger"
use directed_scm::*;
use rand::SeedableRng;
use rand_xoshiro::Xoshiro256StarStar;
use rand::distributions::WeightedIndex;
use std::fs;

use clap::Parser;

// Advanced users may set hardcoded target/relaxed simplex counts boundaries here.
// They overwrite target_relaxation command line option.
const TARGET_BOUNDS:Bounds =  Bounds{flag_count_min:vec![], flag_count_max:vec![]};
const RELAXED_BOUNDS:Bounds =  Bounds{flag_count_min:vec![], flag_count_max:vec![]};

// other stuff to tinker with: Obscure enough to not bother with command line options
const MOVE_DISTRIBUTION_SIMPLE:[f64; 4] = [0.5, 0.5, 0.0, 0.0];        //simple moves only
const MOVE_DISTRIBUTION:[f64; 4] = [0.1, 0.1, 0.6, 0.2];


/// MCMC sampler for flag complexes of a directed graph
#[derive(Parser, Debug)]
#[clap(version, about, long_about = None)]
struct Args{
    // Sampler configuration
    /// flag input file location
    #[clap(short, long)]
    input: String,

    /// target relaxation (percentage)
    #[clap(short, long, default_value_t = 0.01)]
    target_relaxation: f64,
       
    /// number of samples to draw
    #[clap(short, long, default_value_t = 1000)]
    number_of_samples: usize,


    // Sampler configuration: Technical
    /// short human-readable label for reference
    #[clap(short, long)]
    label: String,

    /// random seed 
    #[clap(short, long, default_value_t = 0)]
    seed: u64,

    /// sample_distance: number of steps between too samples
    #[clap(long, default_value_t = 0)]
    sample_distance: usize,


    // Caching, saving, restoring
    /// continue: load given serde file and continue from there.
    /// Takes precedence over almost all other options if present.
    #[clap(short, long, default_value = "")]
    continue_from: String,

    /// samples_store_dir: directory where to store hdf5 files with samples
    #[clap(long, default_value = "./samples/")]
    samples_store_dir: String,

    /// state_store_dir: directory where to store saved states
    #[clap(long, default_value = "./state/")]
    state_store_dir: String,

    /// state_save_interval: Interval in which to save states
    #[clap(long, default_value_t = 100)]
    state_save_interval: usize,

    /// save output as edge list bit array instead of hdf5
    #[clap(long)]
    save_bits: bool,

    /// only use simple single edge flips and double edge moves
    #[clap(long)]
    simple: bool,

}

fn initialize_new_sampler(args: &Args) -> MCMCSampler<Xoshiro256StarStar> {
    let g = io::read_flag_file(&args.input);

    let st = State::new(g);

    println!("we have the following number of maximal k-cliques {:?}", st.cliques_by_order.iter().map(|cs| cs.len()).collect::<Vec<usize>>());
    let rng = Xoshiro256StarStar::seed_from_u64(args.seed);
    let adjusted_clique_order = st.cliques_by_order.iter().map(|cs| (cs.len() as f64).powf(0.2)).collect::<Vec<f64>>();
    let clique_order_distribution = WeightedIndex::new(adjusted_clique_order).unwrap();
    let target_bounds = if TARGET_BOUNDS.flag_count_min.len() == 0 || TARGET_BOUNDS.flag_count_max.len() == 0 {
        let flag_count_min = st.flag_count.iter().enumerate().map(|(d, &scd)| if d < 2 {scd} else {(scd as f64 * (1. - args.target_relaxation)).floor() as usize}).collect();
        let flag_count_max = st.flag_count.iter().enumerate().map(|(d, &scd)| if d < 2 {scd} else {(scd as f64 * (1. + args.target_relaxation)).floor() as usize}).collect();
        Bounds{flag_count_min, flag_count_max}
    } else {
        TARGET_BOUNDS
    };
    let bounds = if RELAXED_BOUNDS.flag_count_min.len() == 0 || RELAXED_BOUNDS.flag_count_max.len() == 0 {
        Bounds::calculate(&st, target_bounds.clone())
    } else {
        RELAXED_BOUNDS
    };
    let move_distribution: WeightedIndex<f64> = WeightedIndex::new(if args.simple { MOVE_DISTRIBUTION_SIMPLE} else { MOVE_DISTRIBUTION }).unwrap();
    let sample_distance = if args.sample_distance == 0 {(2. * st.flag_count[1] as f64 * (st.flag_count[1] as f64).log2()).ceil() as usize} else {args.sample_distance};
    println!("The sampling distance was set to {sample_distance}.");
    return MCMCSampler{state: st, move_distribution, clique_order_distribution, sample_distance: sample_distance, accepted: 0, sampled: 0, rng, bounds};
}

fn main() {
    let args = Args::parse();
    fs::create_dir_all(&args.state_store_dir).expect("could not create state storage directory");
    fs::create_dir_all(&args.samples_store_dir).expect("could not create samples storage directory");

    let save_hdf5: bool = ! args.save_bits;
    
    let (sample_index_start, mut sampler) = if !args.continue_from.is_empty() {
        io::load_state(&args.continue_from).expect("unable to load state")
    } else {
        if save_hdf5 {
            io::new_hdf_file(&args.samples_store_dir, &args.label, args.seed).unwrap();
        }
        (0, initialize_new_sampler(&args))
    };

    let mut bit_output: Option<io::BitOutput> = if args.save_bits {
        Some(io::BitOutput::new(&sampler.state.graph, &format!("{}/{}-{:03}", args.samples_store_dir, args.label, args.seed)).unwrap())
    } else { None };

    let sample_index_end = sample_index_start + args.number_of_samples;
    for i in sample_index_start..sample_index_end {
        if (i % args.state_save_interval) == 0 {
            println!("saving state in step {i}");
            io::save_state(&format!("{state_store_dir}/sampler-{l}-{s:03}.state", state_store_dir = args.state_store_dir, l=args.label, s=args.seed), i, &sampler).unwrap();
        }
        let s = sampler.next();
        if save_hdf5 {
            io::save_to_hdf(&args.samples_store_dir, &args.label, args.seed, i, &s.graph, &s.flag_count).unwrap();
        }
        if let Some(out) = bit_output.as_mut() {
            out.save(&s.graph).unwrap();
        }

        println!("flag count: {:?}", s.flag_count);
        drop(s);
        dbg!(sampler.acceptance_ratio());
    }

    io::save_state(&format!("{state_store_dir}/sampler-{l}-{s:03}.state", state_store_dir = args.state_store_dir, l=args.label, s=args.seed), sample_index_end, &sampler).unwrap();
}

