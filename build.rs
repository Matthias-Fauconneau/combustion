fn main() { println!("cargo:rustc-link-search=../Cantera/build/lib\ncargo:rustc-link-lib=c++\ncargo:rustc-link-lib=fmt\ncargo:rustc-link-lib=yaml-cpp\ncargo:rustc-link-lib=sundials_cvodes"); }
