use {fehler::throws, error::Error};
#[throws] fn main() {
	println!("cargo:rustc-link-lib=cantera");
	println!("cargo:rerun-if-changed=cantera.cpp");
	let out = std::env::var("OUT_DIR")?;
	println!("cargo:rustc-link-search={}", out);
	let out = std::path::PathBuf::from(out);
	let o = out.join("cantera.o");
	std::process::Command::new("c++").args(&["-fPIC","-O2","-c","cantera.cpp","-o",o.to_str().unwrap()]).spawn()?.wait()?;
	let a = out.join("libcantera-rust.a");
	std::process::Command::new("ar").args(&["r", a.to_str().unwrap(), o.to_str().unwrap()]).spawn()?.wait()?;
	println!("cargo:rustc-link-lib=static=cantera-rust");
	println!("cargo:rustc-link-lib=c++");
}
