use directed_scm::*;
use std::marker::PhantomData;

fn main() {
    let g = examples::gengraph();
    let st = State::new(g);
    let mut sampler = MCMCSampler {state: st, burn_in: 2000, sample_distance: 2000, accepted: 0, sampled: 0, rng: rand::thread_rng(), _a: PhantomData::<Shuffling>::default() };
    sampler.burn_in();
    for i in 0..1000 {
        let s = sampler.next();
        dbg!(sampler.acceptance_ratio());
        dbg!(s.flag_count);
    }
}
