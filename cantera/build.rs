fn main() {
	println!("cargo:rustc-link-lib=cantera");
	println!("cargo:rerun-if-changed=cantera.cpp");
	let o = std::path::PathBuf::from(std::env::var("OUT_DIR").unwrap()).join("cantera.o");
	std::process::Command::new("clang++").args(&["-fPIC","-O2","-c","cantera.cpp","-o",o.to_str().unwrap()]).spawn().unwrap().wait().unwrap();
	let a = std::path::PathBuf::from(std::env::var("OUT_DIR").unwrap()).join("libcantera-rust.a");
	std::process::Command::new("ar").args(&["r", a.to_str().unwrap(), o.to_str().unwrap()]).spawn().unwrap().wait().unwrap();
	println!("cargo:rustc-link-search={}",std::env::var("OUT_DIR").unwrap());
	println!("cargo:rustc-link-lib=static=cantera-rust");
	println!("cargo:rustc-link-lib=c++");
}
