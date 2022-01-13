
use bindgen::builder;

use std::path::Path;

fn main() {

    // Configure and generate bindings.
    let bindings = builder().header("flagser.hpp")
        .allowlist_type("directed_graph_t")
        .allowlist_type("filtered_directed_graph_t")
        .opaque_type("std::.*")
        .allowlist_function("directed_graph_t::add_edge")
        .allowlist_function("filtered_directed_graph_t::add_edge")
        .allowlist_function("count_cells")

        .generate().unwrap();

    // Write the generated bindings to an output file.
    let out_dir = std::env::var("OUT_DIR").unwrap();
    let dest = Path::new(&out_dir).join("flagser_bindgen.rs");
    bindings.write_to_file(&dest).unwrap();
    println!("cargo:rerun-if-changed=build.rs");
    println!("cargo:rerun-if-changed=flagser.hpp");
    println!("cargo:rerun-if-changed={}", dest.to_str().unwrap());


     cc::Build::new()
        .file("libflagser.cpp")
        .compile("flagser");
    println!("cargo:rustc-link-lib=flagser");
    println!("cargo:rustc-link-lib=stdc++");
}
