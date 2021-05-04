#![allow(incomplete_features, non_snake_case)]#![feature(const_generics, const_evaluatable_checked, bool_to_option, array_methods, slice_pattern)]

fn main() -> Result<(), Box<dyn std::error::Error>> {
	let model = &std::fs::read("CH4+O2.ron")?;
	use combustion::*;
	let model = model::Model::new(&model)?;
	let (_species_names, species) = combustion::Species::new(&model.species);
	let OUT_DIR = std::env::var("OUT_DIR").unwrap();
	std::process::Command::new("nvcc").args(&["-I",&OUT_DIR,&format!("-DSPECIES={}", species.len()),"--ptx","main.cu","-o",&format!("{}/main.ptx",OUT_DIR)]).spawn()?.wait()?
		.success().then_some(()).unwrap();
	println!("cargo:rerun-if-changed=CH4+O2.ron");
	println!("cargo:rerun-if-changed=main.cu");
	println!("cargo:rustc-env=SPECIES={}", species.len());
	Ok(())
}
