#![allow(incomplete_features, non_snake_case)]#![feature(const_generics, const_evaluatable_checked, bool_to_option, non_ascii_idents, array_methods, slice_pattern)]
use core::slice::SlicePattern;

fn main() -> Result<(), Box<dyn std::error::Error>> {
	let model = &std::fs::read("CH4+O2.ron")?;
	use combustion::{*, transport::*};
	let model = model::Model::new(&model)?;
	let (_species_names, species) = combustion::Species::new(model.species);
	let transport_polynomials = species.transport_polynomials();
	let combustion::Species{molar_mass, thermodynamics, ..} = species;
	/*#[derive(Debug)] struct Falloff<const S: usize> where [(); S-1]: {
		reaction: Reaction<S>,
		efficiency: usize,
	}
		elementary: Box<[Reaction<S>]>,
		three_body: Box<[Reaction<S>]>,
		pressure_modification: Box<[Reaction<S>]>,
		efficiencies: Box<[[f64; S]]>,
		falloff: Box<[Falloff<S>]>,*/
	#[derive(Debug)] struct Species {
		molar_mass: Box<[f64]>,
		thermodynamics: Box<[NASA7]>,
		transport_polynomials: TransportPolynomials,
	}
	/*let mut efficiencies_set = Vec::new();
	let falloff = reactions.iter().filter_map(|reaction|
		if let Model::Falloff{efficiencies, ..} = reaction.model {
			if !efficiencies_set.contains(&efficiencies) { efficiencies_set.push(efficiencies); }
			Some(Falloff{reaction: *reaction, efficiency: efficiencies_set.iter().position(|&e| e == efficiencies).unwrap()})
		} else { None }
	).collect();*/

	let species = Species{
		molar_mass,
		thermodynamics,
		/*elementary: reactions.iter().filter(|r| matches!(r.model, Model::Elementary{})).copied().collect(),
		three_body: reactions.iter().filter(|r| matches!(r.model, Model::ThreeBody{..})).copied().collect(),
		pressure_modification: reactions.iter().filter(|r| matches!(r.model, Model::PressureModification{..})).copied().collect(),
		efficiencies: efficiencies_set.into_boxed_slice(),
		falloff,*/
		transport_polynomials,
	};
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
	}
	/*impl<const S: usize> CUDA for Reaction<S> where [(); S-1]: {
		fn r#type(&self) -> String { use Model::*; match self.model {Elementary => "Elementary", ThreeBody{..} => "ThreeBody", PressureModification{..} => "PressureModification", Falloff{..} => "Falloff"}.to_string()}
		fn cuda(&self) -> String {
			impl CUDA for RateConstant { fn cuda(&self) -> String { format!("{}{{{}, {}, {}}}", self.r#type(), self.log_preexponential_factor, self.temperature_exponent, self.activation_temperature) } }
			format!("{}{{Reaction{{{}, {}, {}, {}, {}, {}, {}}}{}", self.r#type(), self.reactants.cuda(), self.products.cuda(), self.net.cuda(), self.Σreactants, self.Σproducts, self.Σnet, self.rate_constant.cuda(),
				{use Model::*; match self.model {
					Elementary => format!("}}"),
					ThreeBody{efficiencies} => format!(", {}}}", efficiencies.as_slice().cuda()),
					PressureModification{efficiencies, k0} => format!(", {}, {}}}", efficiencies.as_slice().cuda(), k0.cuda()),
					Falloff{..} => format!(""),
				}}
			)
		}
	}
	impl<const S: usize> CUDA for Falloff<S> where [(); S-1]: {
		fn r#type(&self) -> String { "Falloff".to_string() }
 		fn cuda(&self) -> String {
			if let Model::Falloff{efficiencies: _, k0, troe: ron::Troe{A, T3, T1, T2}} = self.reaction.model {
				format!("{}, {}, {}, {}, {}, {}, {}}}", self.reaction.cuda(), self.efficiency, k0.cuda(), A, T3, T1, T2)
			} else { unreachable!() }
		}
	}*/
	impl CUDA for TransportPolynomials {
		//fn r#type(&self) -> String { "TransportPolynomials".to_string() }
 		fn cuda(&self) -> String { format!("{}{{{},{},{}}}", self.r#type(), self.sqrt_viscosity_T14.as_slice().cuda(), self.thermal_conductivity_T12.cuda(), self.binary_thermal_diffusion_coefficients_T32.cuda()) }
	}
	impl CUDA for NASA7 { fn r#type(&self) -> String { self.0.r#type() } fn cuda(&self) -> String { self.0.cuda() } }
	impl CUDA for Species {
		fn r#type(&self) -> String { Default::default() /*"Species".to_string()*/ }
		fn cuda(&self) -> String {
			format!("{},\n{},\n{}", //self.r#type(),
				self.molar_mass.cuda(),
				self.thermodynamics.cuda(),
				/*self.elementary.as_ref().cuda(),
				self.three_body.as_ref().cuda(),
				self.pressure_modification.as_ref().cuda(),
				self.efficiencies.as_ref().cuda(),
				self.falloff.as_ref().cuda(),*/
				self.transport_polynomials.cuda()
			)
		}
	}
	println!("SPECIES {}", species.thermodynamics.len());
	/*println!("ELEMENTARY {}", system.elementary.len());
	println!("THREE_BODY {}", system.three_body.len());
	println!("PRESSURE_MODIFICATION {}", system.pressure_modification.len());
	println!("EFFICIENCIES {}", system.efficiencies.len());
	println!("FALLOFF {}", system.falloff.len());*/

	let OUT_DIR = std::env::var("OUT_DIR").unwrap();
	std::fs::write(std::path::PathBuf::from(&OUT_DIR).join("species.h"), species.cuda())?;
	std::process::Command::new("nvcc").args(&["-I", &OUT_DIR,
		&format!("-DSPECIES={}", species.thermodynamics.len()),
		/*&format!("-DELEMENTARY={}", system.elementary.len()),
		&format!("-DTHREE_BODY={}", system.three_body.len()),
		&format!("-DPRESSURE_MODIFICATION={}", system.pressure_modification.len()),
		&format!("-DEFFICIENCIES={}", system.efficiencies.len()),
		&format!("-DFALLOFF={}", system.falloff.len()),*/
		"--ptx","main.cu","-o",&format!("{}/main.ptx",OUT_DIR)
	]).spawn()?.wait()?.success().then_some(()).unwrap();
	println!("cargo:rerun-if-changed=CH4+O2.ron");
	println!("cargo:rerun-if-changed=main.cu");
	Ok(())
}
