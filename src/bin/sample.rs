use directed_scm::*;
use std::marker::PhantomData;
use rand::prelude::*;

fn main() {
    let g = examples::gengraph();
    let g = io::read_flag_file("/tmp/c.elegans.flag");
    let st = State::new(g);
    println!("we have the following number of maximal cliques {:?}", &st.cliques.iter().map(|c| c.len()).sum::<usize>());
    //let mut sampler = MCMCSampler {state: st, burn_in: 2000, sample_distance: 2000, accepted: 0, sampled: 0, rng: rand::thread_rng(), _a: PhantomData::<Shuffling>::default() };
    let rng = rand::rngs::StdRng::seed_from_u64(2);
    let mut sampler = MCMCSampler {state: st, burn_in: 000, sample_distance: 2000, accepted: 0, sampled: 0, rng, _a: PhantomData::<Shuffling>::default() };
    sampler.burn_in();
    for i in 0..100 {
        let s = sampler.next();
        dbg!(sampler.acceptance_ratio());
        println!("flag count: {:?}", s.flag_count);
    }
}

