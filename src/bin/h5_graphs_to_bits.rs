use clap::Parser;
use directed_scm::io;
use flag_complex::prelude::*;

#[derive(Parser)]
struct Args {
    #[clap(short, long)]
    input: String,
    #[clap(short, long)]
    output: String,
}

pub fn load_graph_hdf5<G: DirectedGraphNew>(ds: &hdf5::Dataset) -> Result<G, Box<dyn std::error::Error>> {
    let arr = ds.read_2d()?;
    let nedges = arr.nrows();
    assert!(arr.ncols() == 2);
    let nnodes = ds.attr("number_of_vertices")?.read_1d()?[0];
    let mut g = G::new_disconnected(nnodes);
    for i in 0..nedges {
        let a = arr[(i, 0)];
        let b = arr[(i, 1)];
        g.add_edge(a, b);
    }
    return Ok(g);
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args = Args::parse();
    let f = hdf5::File::open(&args.input)?;
    

    let seed_group_vec = f.groups()?;
    match &seed_group_vec[..] {
        [seed_group] => {
            let name = seed_group.name();
            eprintln!("processing seed {}", name);
            let mut sample_iter = seed_group.groups()?.into_iter().map(|g| g.dataset("edgelist").expect("no .../edge_list"));
            let first : flag_complex::Graph = load_graph_hdf5(&sample_iter.next().unwrap())?;

            let mut writer = directed_scm::io::BitOutput::new(&first, &args.output)?;
            writer.save(&first)?;
            
            for sample_ds in sample_iter {
                let graph : flag_complex::Graph = load_graph_hdf5(&sample_ds)?;
                writer.save(&graph)?;
            }
        }
        _ => {
            eprintln!("{} has not exactly one group", args.input);
            std::process::exit(1);
        }
    }
    Ok(())
}
