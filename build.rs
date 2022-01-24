fn main() {
    println!("cargo:rerun-if-changed=build.rs");

    println!("cargo:rerun-if-changed=flagser.hpp");
    println!("cargo:rerun-if-changed=libflagser.cpp");

    cc::Build::new()
        .file("libflagser.cpp")
//        .flag("-march=native")
        .cpp(true)
        .compile("flagser");
}
