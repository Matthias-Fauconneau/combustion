#![allow(incomplete_features, non_snake_case)]#![feature(const_generics, const_evaluatable_checked, type_ascription, non_ascii_idents, array_methods, slice_pattern, bindings_after_at)]
use core::slice::SlicePattern;

fn main() -> Result<(), Box<dyn std::error::Error>> {
	println!("cargo:rerun-if-changed=CH4+O2.ron");
	println!("cargo:rerun-if-changed=main.comp");
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
		if let ReactionModel::Falloff{efficiencies, ..} = reaction.model {
			if !efficiencies_set.contains(&efficiencies) { efficiencies_set.push(efficiencies); }
			Some(Falloff{reaction: *reaction, efficiency: efficiencies_set.iter().position(|&e| e == efficiencies).unwrap()})
		} else { None }
	).collect();*/

	let species = Species{
		molar_mass,
		thermodynamics,
		/*elementary: reactions.iter().filter(|r| matches!(r.model, ReactionModel::Elementary{})).copied().collect(),
		three_body: reactions.iter().filter(|r| matches!(r.model, ReactionModel::ThreeBody{..})).copied().collect(),
		pressure_modification: reactions.iter().filter(|r| matches!(r.model, ReactionModel::PressureModification{..})).copied().collect(),
		efficiencies: efficiencies_set.into_boxed_slice(),
		falloff,*/
		transport_polynomials,
	};
	trait GLSL {
		fn r#type(&self) -> String { std::any::type_name::<Self>().rsplit("::").next().unwrap().to_string() }
		fn glsl(&self) -> String;
	}
	impl GLSL for f64 { fn r#type(&self) -> String { "float".to_string() } fn glsl(&self) -> String { format!("{}", self) } }
	impl<T: GLSL> GLSL for &[T] {
		fn r#type(&self) -> String { format!("{}[]", self.first().unwrap().r#type()) }
		fn glsl(&self) -> String { use itertools::Itertools; format!("{}({})\n", self.r#type(), self.iter().format_with(", ", |e, f| f(&e.glsl()))) }
	}
	impl<T: GLSL, const N: usize> GLSL for [T; N] {
		fn r#type(&self) -> String { format!("{}[]", self.first().unwrap().r#type()) }
		fn glsl(&self) -> String { use itertools::Itertools; format!("{}({})\n", self.r#type(), self.iter().format_with(", ", |e, f| f(&e.glsl()))) }
	}
	impl<T: GLSL> GLSL for Box<[T]> {
		fn r#type(&self) -> String { self.as_slice().r#type() }
		fn glsl(&self) -> String { self.as_slice().glsl() }
	}
	/*impl<const S: usize> GLSL for Reaction<S> where [(); S-1]: {
		fn r#type(&self) -> String { use ReactionModel::*; match self.model {Elementary => "Elementary", ThreeBody{..} => "ThreeBody", PressureModification{..} => "PressureModification", Falloff{..} => "Falloff"}.to_string()}
		fn glsl(&self) -> String {
			impl GLSL for RateConstant { fn glsl(&self) -> String { format!("{}({}, {}, {})", self.r#type(), self.log_preexponential_factor, self.temperature_exponent, self.activation_temperature) } }
			format!("{}(Reaction({}, {}, {}, {}, {}, {}, {}){}", self.r#type(), self.reactants.glsl(), self.products.glsl(), self.net.glsl(), self.Σreactants, self.Σproducts, self.Σnet, self.rate_constant.glsl(),
				{use ReactionModel::*; match self.model {
					Elementary => format!(")"),
					ThreeBody{efficiencies} => format!(", {})", efficiencies.as_slice().glsl()),
					PressureModification{efficiencies, k0} => format!(", {}, {})", efficiencies.as_slice().glsl(), k0.glsl()),
					Falloff{..} => format!(""),
				}}
			)
		}
	}
	impl<const S: usize> GLSL for Falloff<S> where [(); S-1]: {
		fn r#type(&self) -> String { "Falloff".to_string() }
 		fn glsl(&self) -> String {
			if let ReactionModel::Falloff{efficiencies: _, k0, troe: model::Troe{A, T3, T1, T2}} = self.reaction.model {
				format!("{}, {}, {}, {}, {}, {}, {})", self.reaction.glsl(), self.efficiency, k0.glsl(), A, T3, T1, T2)
			} else { unreachable!() }
		}
	}*/
	impl GLSL for TransportPolynomials {
		//fn r#type(&self) -> String { "TransportPolynomials".to_string() }
 		fn glsl(&self) -> String { format!("{}({},{},{})", self.r#type(), self.sqrt_viscosity_T14.as_slice().glsl(), self.thermal_conductivity_T12.as_slice().glsl(), self.binary_thermal_diffusion_coefficients_T32.as_slice().glsl()) }
	}
	impl GLSL for NASA7 { fn r#type(&self) -> String { self.0.r#type() } fn glsl(&self) -> String { self.0.glsl() } }
	impl GLSL for Species {
		//fn r#type(&self) -> String { "Species".to_string() }
		fn glsl(&self) -> String {
			format!("{}({},\n{},\n{})", self.r#type(),
				self.molar_mass.glsl(),
				self.thermodynamics.glsl(),
				/*self.elementary.as_ref().glsl(),
				self.three_body.as_ref().glsl(),
				self.pressure_modification.as_ref().glsl(),
				self.efficiencies.as_ref().glsl(),
				self.falloff.as_ref().glsl(),*/
				self.transport_polynomials.glsl()
			)
		}
	}
	println!("SPECIES {:?}", Some(&species.thermodynamics.len().to_string()));
	/*println!("ELEMENTARY {:?}", Some(&system.elementary.len().to_string()));
	println!("THREE_BODY {:?}", Some(&system.three_body.len().to_string()));
	println!("PRESSURE_MODIFICATION {:?}", Some(&system.pressure_modification.len().to_string()));
	println!("EFFICIENCIES {:?}", Some(&system.efficiencies.len().to_string()));
	println!("FALLOFF {:?}", Some(&system.falloff.len().to_string()));*/

	let mut options = shaderc::CompileOptions::new().unwrap();
	options.add_macro_definition("SPECIES", Some(&species.thermodynamics.len().to_string()));
	/*options.add_macro_definition("ELEMENTARY", Some(&system.elementary.len().to_string()));
	options.add_macro_definition("THREE_BODY", Some(&system.three_body.len().to_string()));
	options.add_macro_definition("PRESSURE_MODIFICATION", Some(&system.pressure_modification.len().to_string()));
	options.add_macro_definition("EFFICIENCIES", Some(&system.efficiencies.len().to_string()));
	options.add_macro_definition("FALLOFF", Some(&system.falloff.len().to_string()));*/
	let OUT_DIR = std::env::var("OUT_DIR").unwrap();
	std::fs::write(std::path::PathBuf::from(&OUT_DIR).join("species.h"), species.glsl())?;
	options.set_include_callback(|_,_,_,_| Ok(shaderc::ResolvedInclude{resolved_name:"species.h".to_string(), content: species.glsl()}));
	std::fs::write(std::path::PathBuf::from(OUT_DIR).join("main.spv"),
											shaderc::Compiler::new().unwrap().compile_into_spirv(include_str!("main.comp"), shaderc::ShaderKind::Compute, "main.comp", "main", Some(&options))?.as_binary_u8())?;
	Ok(())
}
