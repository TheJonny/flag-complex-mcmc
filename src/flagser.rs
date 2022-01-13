 //include!(concat!(env!("OUT_DIR"), "/flagser_bindgen.rs"));

use libc::size_t;

use crate::Node;

#[link(name = "flagser")]
extern {
    fn flagser_count_unweighted(nvertices: size_t, nedges: size_t, edges: *const [Node; 2], res_size: *mut size_t) -> *mut size_t;
}


pub fn count_unweighted(nvertices: size_t, edges: &[[Node; 2]]) -> Vec<usize> {
    let mut nres: usize = 0;
    let res_data = unsafe{flagser_count_unweighted(nvertices, edges.len(), edges.as_ptr(), &mut nres)};
    let res_slice = unsafe{std::slice::from_raw_parts(res_data, nres)};
    let v: Vec<usize> = res_slice.into();
    drop(res_slice);
    unsafe{libc::free(res_data as *mut _)};
    return v;
}
