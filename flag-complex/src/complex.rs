use crate::graph::*;
use rayon::prelude::*;

// This is a translation of parts from the "flagser" C++ code of Daniel Lütgehetmann.
// Published as "Computing Persistent Homology of Directed Flag Complexes" https://doi.org/10.3390/a13010019
// https://github.com/luetge/flagser/blob/master/include/complex/directed_flag_complex.h

pub fn for_each_cell_par<'a, G: DirectedGraph+Sync+?Sized, R, F>
    (graph: &'a G, f: F, min_dimension: usize, max_dimension: usize) -> impl ParallelIterator<Item=R> + 'a
where
    F: Fn(&mut R, &[Node]) + Send + Sync + 'a,
    R: Default+Send+Sync + 'a,
{
    (0 .. graph.nnodes() as Node).into_par_iter().fold(R::default, move |mut state, v| {
        let possible_next_vertices = graph.edges_from(v);
        let mut prefix = vec![v];
        do_for_each_cell(graph, &mut |simplex| f(&mut state, simplex), min_dimension, max_dimension, &mut prefix, &possible_next_vertices);
        state
    })
}

pub fn for_each_cell<G,F>
    (graph: &G, f: &mut F, min_dimension: usize, max_dimension: usize)
where
    G: DirectedGraph+?Sized,
    F: FnMut(&[Node]),
{
    let mut prefix = vec![];
    for v in graph.iter_nodes() {
        prefix.resize(1, v);
        prefix[0] = v;
        let possible_next_vertices = graph.edges_from(v);
        do_for_each_cell(graph, f, min_dimension, max_dimension, &mut prefix, &possible_next_vertices);
    }
}

fn do_for_each_cell<F,G>(graph: &G, f: &mut F, min_dimension: usize,
                         max_dimension: usize,
                         prefix: &mut Vec<Node>,
                         possible_next_vertices: &[Node]
                         )
    where G: DirectedGraph+?Sized, F: FnMut(&[Node])
{
    if prefix.len() >= min_dimension + 1 {
        f(&prefix);
    }
    if prefix.len() == max_dimension.wrapping_add(1) {
        return;
    }

    for &vertex in possible_next_vertices {
        let new_possible_next_vertices = if prefix.len() == 0 {
            graph.edges_from(vertex)
        } else {
            possible_next_vertices.iter().cloned().filter(|&v| graph.has_edge(vertex, v)).collect()
        };
        prefix.push(vertex);
        do_for_each_cell(graph, f, min_dimension, max_dimension, prefix, &new_possible_next_vertices);
        prefix.pop();
    }
}

pub fn count_cells_par<G: DirectedGraph+Sync+?Sized>(graph: &G) -> Vec<usize>{
    let count = |result: &mut Vec<usize>, cell: &[Node]| {
        if result.len() < cell.len() {
            result.resize(cell.len(), 0);
        }
        result[cell.len()-1] += 1;
    };
    let combine = |mut left: Vec<usize>, mut right: Vec<usize>| {
        if left.len() < right.len() {
            std::mem::swap(&mut left, &mut right)
        }
        for i in 0..right.len() {
            left[i] += right[i];
        }
        left
    };
    let res = for_each_cell_par(graph, &count, 0, graph.nnodes()).reduce(Vec::new, &combine);
    return res;
}

pub fn count_cells<G: DirectedGraph+Sync+?Sized>(graph: &G) -> Vec<usize>{
    let mut result = vec![];
    let mut count = |cell: &[Node]| {
        if result.len() < cell.len() {
            result.resize(cell.len(), 0);
        }
        result[cell.len()-1] += 1;
    };
    for_each_cell(graph, &mut count, 0, graph.nnodes());
    return result;
}
