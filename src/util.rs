
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
