// Authors: "Jonathan Krebs and Florian Unger"
use directed_scm::*;
use std::marker::PhantomData;
use rand::prelude::*;
use clap::Parser;

/// MCMC sampler for flag complexes of a directed graph
#[derive(Parser, Debug)]
#[clap(version, about, long_about = None)]
struct Args{
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

    // /// outdir: location where to save sampled flags
    // #[clap(short, long, default_value_t = "./")]
    // outdir: String,
}

fn main() {
    let args = Args::parse();

    let g = examples::gengraph();
    let g = io::read_flag_file(&args.input);

    let st = State::new(g);
    println!("we have the following number of maximal cliques {:?}", &st.cliques.iter().map(|c| c.len()).sum::<usize>());
    //let mut sampler = MCMCSampler {state: st, burn_in: 2000, sample_distance: 2000, accepted: 0, sampled: 0, rng: rand::thread_rng(), _a: PhantomData::<Shuffling>::default() };
    let rng = rand::rngs::StdRng::seed_from_u64(args.seed);
    let mut sampler = MCMCSampler {state: st, burn_in: args.burn_in, sample_distance: args.sample_distance, accepted: 0, sampled: 0, rng, _a: PhantomData::<Shuffling>::default() };
    sampler.burn_in();
    for i in 0..args.number_of_samples {
        let s = sampler.next();
        io::save_flag_file(&format!("{}-{:03}-{:03}", args.input, args.seed, i), &s.graph);
        dbg!(sampler.acceptance_ratio());
        println!("flag count: {:?}", s.flag_count);
    }
}

