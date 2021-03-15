fn main() -> std::result::Result<(),std::io::Error> {
	bindgen::Builder::default().header("ct.h").generate().unwrap().write_to_file(std::path::PathBuf::from(std::env::var("OUT_DIR").unwrap()).join("cantera.rs"))?;
	println!("cargo:rustc-link-lib=dylib=/usr/local/lib/libcantera_shared.so");
	Ok(())
}
