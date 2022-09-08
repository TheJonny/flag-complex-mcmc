use flag_complex::prelude::*;
use flag_complex;
use directed_scm;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    for seed in 0..50 {
        for sample in 0..100 {
            analyze(seed, sample)?;
        }
    }
    Ok(())
}

fn analyze(seed: u64, sample: u64) -> Result<(), Box<dyn std::error::Error>> {
    println!("== seed: {seed}, sample: {sample} : ==");
    let g: flag_complex::Graph = directed_scm::io::load_graph_hdf5(&format!("/scratch/samples/bbp0_l13-{seed:03}.hdf5"), seed, sample)?;
    let edges = directed_scm::io::load_edgelist_hdf5("/scratch/samples/intersections.hdf5", &format!("{:03}", seed))?;
    let s = directed_scm::State::new(g.clone());
    let nei = s.edgeset_neighborhood(&edges);
    let gg = flag_complex::Graph::subgraph(&g, &nei);
    let nei_count = gg.flagser_count();
    println!("whole graph: {:?}", s.flag_count);
    println!("nbhd(intersection): {:?}", nei_count);

    let mut reverse_nei = vec![0; g.nnodes()];
    for (i,j) in nei.iter().cloned().enumerate() {
        reverse_nei[j as usize] = i as u32;
    }
    let mut cut = gg;
    for &[a,b] in &edges {
        assert!(cut.has_edge(reverse_nei[a as usize], reverse_nei[b as usize]));
        cut.remove_edge(reverse_nei[a as usize], reverse_nei[b as usize]);
    }
    let cut_count = cut.flagser_count();
    println!("nbhd(intersection)-intersection: {:?}", cut_count);
    let mut difference = nei_count.clone();
    for (i,&c) in cut_count.iter().enumerate() {
        difference[i] -= c;
    }
    println!("difference: intersection: {:?}", difference);
    Ok(())
}
