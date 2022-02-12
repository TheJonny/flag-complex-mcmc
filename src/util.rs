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
