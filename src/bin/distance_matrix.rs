use hdf5;
use ndarray::prelude::*;

fn sorted_distance(a: &Array2<u32>, b: &Array2<u32>) -> usize{
    let mut ai = 0;
    let mut bi = 0;
    let mut d = 0;
    while ai < a.nrows() && bi < b.nrows() {
        use core::cmp::Ordering::*;
        let aa = [a[(ai, 0)], a[(ai, 1)]];
        let bb = [b[(bi, 0)], b[(bi, 1)]];
        match Ord::cmp(&aa, &bb) {
            Less => {
                ai += 1;
                d += 1;
            }
            Equal => {
                ai += 1;
                bi += 1;
            }
            Greater => {
                bi += 1;
                d += 1;
            }
        }
    }
    d += (a.nrows() - ai) + (b.nrows() - bi);

    return d;
}

fn main() -> Result<(), Box<dyn std::error::Error>> {

    let irange = (0..10000).step_by(10);
    let jrange = (0..10000).step_by(10);

    let inname = "samples/c_elegans_singlestep-000.hdf5";
    let f = hdf5::File::open(inname)?;


    let mut mat = ndarray::Array::<usize,_>::zeros((irange.len(), jrange.len()));

    for (out_j, seq_j) in jrange.enumerate() {
        let edgesj = f.dataset(&format!("000/{:06}/edgelist", seq_j))?.read_2d::<u32>()?;
        assert!(edgesj.ncols() == 2);

        for (out_i, seq_i) in irange.clone().enumerate() {
            let ds = f.dataset(&format!("000/{:06}/edgelist", seq_i))?;
            let edgesi = ds.read_2d::<u32>()?;
            let d = sorted_distance(&edgesi, &edgesj);
            mat[(out_j,out_i)] = d;
        }
    }

    let of = hdf5::File::create("distances.hdf5")?;
    of.new_dataset_builder().with_data(&mat).create("d")?;

    Ok(())
}
