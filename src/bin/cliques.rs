use directed_scm::io;
use directed_scm::graph::*;

fn main(){
    if let [_, infile] = &std::env::args().collect::<Vec<_>>()[..] {
        let g: directed_scm::graph::BoolMatrixGraph = io::read_flag_file(&infile);

        let cs = g.compute_maximal_cliques();
        println!("{}", cs.len());
    }
}
