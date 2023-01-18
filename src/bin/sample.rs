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

    /// only use simple single edge flips and double edge moves
    #[clap(long)]
    simple: bool,

    // Sampler configuration: Technical
    /// short human-readable label for reference
    #[clap(short, long)]
    label: String,

    /// Give seeds to seed the RNG in a comma-separated list.
    /// Each seed starts a separate chain in a seperate thread with a separate RNG.
    /// Usage example: --seeds=0,1,2,3,5 or --seeds=$(seq -s, 1 5) .
    /// Defaults to one seed which is 0.
    #[clap(short, long, use_value_delimiter = true, value_delimiter = ',', default_value = "0")]
    seeds: Vec<u64>,

    /// sample_distance: number of steps between too samples
    #[clap(long, default_value_t = 0)]
    sample_distance: usize,


    // Caching, saving, restoring
    /// continue: load given serde file and continue from there.
    /// Takes precedence over almost all other options if present.
    #[clap(short, long, default_value_t = false, conflicts_with_all(&[
                                                    "input",
                                                    "target_relaxation",
                                                    "simple",
                                                    //"save_bits",
                                                    "sample_distance",
                                                    "state_save_interval",
                                                    ]))]
    resume: bool,

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

}

fn initialize_new(args: &Args) -> (Parameters, Precomputed, MarkovState) {
    let g = io::read_flag_file(&args.input);

    let st = MarkovState::new(g);
    let precomputed = Precomputed::new(&st.graph);

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

    let parameters = Parameters {
        bounds, move_distribution, sample_distance, state_save_distance: args.state_save_distance
    };
    return (parameters, precomputed, st);
}

fn main() {
    let args = Args::parse();
    println!("{:?}", &args.seeds);
    fs::create_dir_all(&args.state_store_dir).expect("could not create state storage directory");
    fs::create_dir_all(&args.samples_store_dir).expect("could not create samples storage directory");

    //let save_hdf5: bool = ! args.save_bits;
    

    let precomputed: Precomputed;
    let parameters: Parameters;
    let samplers : Vec<MCMCSampler<Xoshiro256StarStar>>;
    if args.resume {
        (precomputed, parameters) = io::load_shared(&args.state_dir, &args.label);
        samplers = args.seeds.iter()
            .map(|&seed| io::load_state(&args.state_dir, &args.label, seed).expect("unable to load state"))
            .collect();
    } else {
        let start_state;
        (parameters, precomputed, start_state) = initialize_new(&args);
        io::save_shared(&args.state_dir, &args.label, &parameters, &precomputed);
        for &seed in &args.seeds {
                io::new_hdf_file(&args.samples_store_dir, &args.label, seed, &start_state.graph).unwrap();
        }
        samplers = args.seeds.iter()
            .map(|&seed| {
                let rng = Xoshiro256StarStar::seed_from_u64(args.seed);
                MCMCSampler::new(start_state.clone(), rng, &parameters, &precomputed)
            })
            .collect();
    }

    //let mut bit_output: Option<io::BitOutput> = if args.save_bits {
    //    Some(io::BitOutput::new(&sampler.state.graph, &format!("{}/{}-{:03}", args.samples_store_dir, args.label, args.seed)).unwrap())
    //} else { None };
    std::thread::scope(|scope| {
        for (sampler, &seed) in samplers.iter().zip(args.seeds.iter()) {
            scope.spawn(||{
                // run the sampler for args.number_of_stamples steps.
                // save state
                loop {
                    let s = sampler.step();
                    if sampler.step_number() % parameters.state_save_distance == 0 {
                        println!("saving state in step {}", sampler.step_number());
                        io::save_state(&format!("{state_store_dir}/sampler-{l}-{s:03}.state", state_store_dir = args.state_store_dir, l=args.label, s=seed), &sampler).unwrap();
                        // save state
                    }
                    if sampler.step_number % parameters.sample_distance == 0 {
                        let output_index = sampler.step_number() / parameters.sample_distance;
                        // write output
                        println!("flag count: {:?}", s.flag_count);
                        dbg!(sampler.acceptance_ratio());
                        io::save_to_hdf(&args.samples_store_dir, &args.label, seed, output_index, &s.graph, &s.flag_count).unwrap();
                    }
                }
            });
        }
    });
}
