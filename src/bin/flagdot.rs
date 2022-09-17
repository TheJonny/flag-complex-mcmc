use directed_scm::io;
// Authors: "Jonathan Krebs and Florian Unger"

fn main() -> std::io::Result<()>{
    let inname = std::env::args().nth(1);
    if let Some(inname) = inname {
        let g : flag_complex::Graph = io::read_flag_file(&inname);
        return io::save_dot(&mut std::io::stdout(), &g);
    }
    Ok(())
}
