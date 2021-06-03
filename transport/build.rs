#![feature(bool_to_option)]
fn main() {
	let file = "/usr/share/cantera/data/LiDryer.yaml";
	let model = yaml::Loader::load_from_str(std::str::from_utf8(&std::fs::read(file).unwrap()).unwrap()).unwrap();
	let species = yaml::parse(&model).unwrap().species;
	let _ = std::process::Command::new("nvcc").args(&[&format!("-DSPECIES={}", species.len()),"--ptx","main.cu","-o",&format!("{}/main.ptx",std::env::var("OUT_DIR").unwrap())])
		.spawn().map(|mut c| c.wait().unwrap().success().then_some(()).unwrap());
	println!("cargo:rerun-if-changed={}", file);
	println!("cargo:rerun-if-changed=main.cu");
	println!("cargo:rustc-env=SPECIES={}", species.len());
}
