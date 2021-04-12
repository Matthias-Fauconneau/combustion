fn main() -> std::result::Result<(),std::io::Error> {
	//bindgen::Builder::default().header("ct.h").generate().unwrap().write_to_file(std::path::PathBuf::from(std::env::var("OUT_DIR").unwrap()).join("cantera.rs"))?;
	//println!("cargo:rustc-flags=-l dylib=stdc++");
	println!("cargo:rustc-link-lib=stdc++");
	println!("cargo:rustc-link-lib=yaml-cpp");
	println!("cargo:rustc-link-lib=sundials_cvodes");
	println!("cargo:rustc-link-lib=cantera");
	Ok(())
}
