use directed_scm::io;

fn main(){
    if let [_, infile] = &std::env::args().collect::<Vec<_>>()[..] {
        let g = io::read_flag_file(&infile);

        let cs = g.compute_maximal_cliques();
        println!("{}", cs.len());
    }
}
