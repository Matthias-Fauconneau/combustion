#![allow(incomplete_features, non_snake_case)]#![feature(const_generics, const_evaluatable_checked, bool_to_option, non_ascii_idents, array_methods, slice_pattern)]
//use core::slice::SlicePattern;

fn main() -> Result<(), Box<dyn std::error::Error>> {
	let model = &std::fs::read("CH4+O2.ron")?;
	use combustion::*;
	let model = model::Model::new(&model)?;
	let (_species_names, species) = combustion::Species::new(&model.species);
	/*let transport_polynomials = species.transport_polynomials();
	let combustion::Species{molar_mass, thermodynamics, ..} = species;
	trait CUDA {
		fn r#type(&self) -> String { std::any::type_name::<Self>().rsplit("::").next().unwrap().to_string() }
		fn cuda(&self) -> String;
	}
	impl CUDA for f64 { fn r#type(&self) -> String { "double".to_string() } fn cuda(&self) -> String { format!("{}", self) } }
	impl<T: CUDA> CUDA for &[T] {
		fn r#type(&self) -> String { Default::default() /*format!("{}[]", self.first().unwrap().r#type())*/ }
		fn cuda(&self) -> String { use itertools::Itertools; format!("{}{{{}}}\n", self.r#type(), self.iter().format_with(",", |e, f| f(&e.cuda()))) }
	}
	impl<T: CUDA, const N: usize> CUDA for [T; N] {
		fn r#type(&self) -> String { Default::default() /*format!("{}[]", self.first().unwrap().r#type())*/ }
		fn cuda(&self) -> String { use itertools::Itertools; format!("{}{{{}}}\n", self.r#type(), self.iter().format_with(",", |e, f| f(&e.cuda()))) }
	}
	impl<T: CUDA> CUDA for Box<[T]> {
		fn r#type(&self) -> String { self.as_slice().r#type() }
		fn cuda(&self) -> String { self.as_slice().cuda() }
	}*/
	let OUT_DIR = std::env::var("OUT_DIR").unwrap();
	//std::fs::write(std::path::PathBuf::from(&OUT_DIR).join("species.h"), species.cuda())?;
	std::process::Command::new("nvcc").args(&["-I", &OUT_DIR,
		&format!("-DSPECIES={}", species.thermodynamics.len()),
		"--ptx","main.cu","-o",&format!("{}/main.ptx",OUT_DIR)
	]).spawn()?.wait()?.success().then_some(()).unwrap();
	println!("cargo:rerun-if-changed=CH4+O2.ron");
	println!("cargo:rerun-if-changed=main.cu");
	println!("cargo:rustc-env=SPECIES={}", species.len());
	Ok(())
}
