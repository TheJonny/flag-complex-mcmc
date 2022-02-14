use crate::Node;
use rand::Rng;
use rand::prelude::SliceRandom;

pub fn intersect_sorted<T: std::cmp::Ord + Clone>(a: &[T], b: &[T]) -> Vec<T> {
    let mut out = vec![];
    let mut ai = 0;
    let mut bi = 0;
    while ai < a.len() && bi < b.len() {
        use std::cmp::Ordering;
        match Ord::cmp(&a[ai], &b[bi]) {
            Ordering::Equal =>  {
                out.push(a[ai].clone());
                ai += 1;
                bi += 1;
            }
            Ordering::Less => {
                ai += 1;
            }
            Ordering::Greater => {
                bi += 1;
            }
        }
    }
    return out;
}

pub fn random_perm<R: Rng>(l: usize, h:usize, rng: &mut R) -> Vec<usize> {
    let mut perm : Vec<usize> = (l..h).collect();
    perm.shuffle(rng);
    return perm;
}

pub fn vec_intersect(xs: &Vec<Node>, ys: &Vec<Node>) -> Vec<Node> {
    // assumes uniqueness
    let mut r = vec![];
    for x in xs {
        if ys.contains(&x) {
            r.push(x.clone());
        }
    }
    return r;
}

pub fn vec_setminus(xs: &Vec<Node>, ys: &Vec<Node>) -> Vec<Node> {
    // xs - ys: return vec contains xs minus the elements in ys
    let mut r = xs.clone();
    r.retain(|x| !ys.contains(x));
    return r;
}

// pub fn vec_max() {}

pub fn all_le<T: PartialOrd> (a: &[T], b: &[T], z: &T) -> bool{
    let maxlen = std::cmp::max(a.len(), b.len());
    let left = a.iter().chain(std::iter::repeat(z));
    let right = b.iter().chain(std::iter::repeat(z));
    for (l,r) in left.zip(right).take(maxlen) {
        if l > r {
            return false;
        }
    }
    return true;
}

pub fn factorial(x: usize) -> usize {
    let mut res = 1 as usize;
    for i in 1..x {
        res *= i;
    }
    return res;
}

pub fn calc_relax_de(sc: &Vec<usize>) -> Vec<usize> {
    let mut relax_de = vec![];
    for d in 0..sc.len() {
        let mut ind = 1;
        let mut simplices_lost = vec![];
        while OEIS_A058298[ind] < sc[d] {
            simplices_lost.push(OEIS_A058298[ind] - OEIS_A058298[ind-1]);
            ind += 1;
        }
        let relax_de_a = simplices_lost.iter().max().unwrap();
        let relax_de_b = factorial(d+1);
        relax_de.push(std::cmp::min(*relax_de_a, relax_de_b));
    }
    return relax_de;
}

// OEIS: A058298, see https://oeis.org/A058298
// This describes the maximum number of simplices achieved by having 1.. double edges in a clique.
// Other description: 	Triangle n!/(n-k), 1 <= k < n, read by rows. 
const OEIS_A058298 : [usize;40] = [2,3,6,8,12,24,30,40,60,120,144,180,240,360,720,
 840,1008,1260,1680,2520,5040,5760,6720,8064,10080,
 13440,20160,40320,45360,51840,60480,72576,90720,
 120960,181440,362880,403200,453600,518400,604800];


#[test]
fn test_intersect() {
    let a = [1,2,        5,8,9];
    let b = [0,2,2,2,3,4,5,  9,9,9];
    let i = intersect_sorted(&a, &b);
    assert_eq!(i, vec![2,5,9]);

    let a = [0,2,2,2,3,4,5,  9,9,9];
    let b = [1,2,        5,8,9];
    let i = intersect_sorted(&a, &b);
    assert_eq!(i, vec![2,5,9]);

    let a = [];
    let b = [0,2,2,2,3,4,5,  9,9,9];
    let i = intersect_sorted(&a, &b);
    assert_eq!(i, vec![]);

    let a = [2];
    let b = [3];
    let i = intersect_sorted(&a, &b);
    assert_eq!(i, vec![]);

    let a = [2];
    let b = [0,1,2,3,4,5];
    let i = intersect_sorted(&a, &b);
    let j = intersect_sorted(&a, &b);
    assert_eq!(i,j);
    assert_eq!(i,vec![2]);

    let a = [0];
    let b = [0,1,2,3,4,5];
    let i = intersect_sorted(&a, &b);
    let j = intersect_sorted(&a, &b);
    assert_eq!(i,j);
    assert_eq!(i,vec![0]);

    let a = [0];
    let b = [0,1,2,3,4,5];
    let i = intersect_sorted(&a, &b);
    let j = intersect_sorted(&a, &b);
    assert_eq!(i,j);
    assert_eq!(i,vec![0]);

    let a = [5];
    let b = [0,1,2,3,4,5];
    let i = intersect_sorted(&a, &b);
    let j = intersect_sorted(&a, &b);
    assert_eq!(i,j);
    assert_eq!(i,vec![5]);
}
