
fn main() {
    println!("cargo:rerun-if-changed=build.rs");
    println!("cargo:rerun-if-changed=flagser.hpp");
    println!("cargo:rerun-if-changed=libflagser.cpp");

     cc::Build::new()
        .file("libflagser.cpp")
        .compile("flagser");
    println!("cargo:rustc-link-lib=flagser");
    println!("cargo:rustc-link-lib=stdc++");
}

