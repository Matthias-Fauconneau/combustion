#![allow(incomplete_features)]#![feature(const_generics, const_evaluatable_checked, type_ascription, non_ascii_idents)]#![allow(non_snake_case)]

fn main() -> Result<(), Box<dyn std::error::Error>> {
	#![allow(unused_variables,dead_code)]
	println!("cargo:rerun-if-changed=H2+O2.ron");
	println!("cargo:rerun-if-changed=main.comp");
	let system = std::fs::read("H2+O2.ron")?;
	type Simulation<'t> = combustion::Simulation::<'t, 9>;
	let combustion::System{molar_masses, reduced_molar_masses, thermodynamics, reactions, ..} = Simulation::new(&system)?.system;
	use combustion::*;
	struct System<const S: usize> where [(); S-1]: {
		molar_masses: [f64; S],
		reduced_molar_masses: [f64; S-1],
		thermodynamics: [NASA7; S],
		elementary: Box<[Reaction<S>]>,
		three_body: Box<[Reaction<S>]>,
		falloff: Box<[Reaction<S>]>,
	}
	let system = System{
		molar_masses,
		reduced_molar_masses,
		thermodynamics,
		elementary: reactions.iter().filter(|r| matches!(r.model, Model::Elementary{})).copied().collect(),
		three_body: reactions.iter().filter(|r| matches!(r.model, Model::ThreeBody{..})).copied().collect(),
		falloff: reactions.iter().filter(|r| matches!(r.model, Model::Falloff{..})).copied().collect()
	};
	let mut options = shaderc::CompileOptions::new().unwrap();
	options.add_macro_definition("SPECIES", Some(&system.thermodynamics.len().to_string()));
	options.add_macro_definition("ELEMENTARY", Some(&system.elementary.len().to_string()));
	options.add_macro_definition("THREE_BODY", Some(&system.three_body.len().to_string()));
	options.add_macro_definition("FALLOFF", Some(&system.falloff.len().to_string()));
	impl<const S: usize> std::fmt::Display for System<S> where [(); S-1]: {
		fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
			use itertools::Itertools;
			let fmt = |reactions:&[Reaction<S>]| {
				let variant = {use Model::*; match reactions.first().unwrap().model {Elementary => "Elementary", ThreeBody{..} => "ThreeBody", Falloff{..} => "Falloff"}};
				format!("{}[]({})", variant, reactions.iter().format_with(",\\\n", |r, f| {
					let fmt = |r:RateConstant| format!("RateConstant({}, {}, {})", r.log_preexponential_factor, r.temperature_exponent, r.activation_temperature);
					f(&format_args!("{}(Reaction(double[]({}), double[]({}), double[]({}), {}, {}, {}, {})",
						variant, r.reactants.iter().format(", "), r.products.iter().format(", "), r.net.iter().format(", "), r.Σreactants, r.Σproducts, r.Σnet, fmt(r.rate_constant)))?;
					{use Model::*; match r.model {Elementary => f(&")"),
						ThreeBody{efficiencies} => f(&format_args!(", double[]({}))", efficiencies.iter().format(", "))),
						Falloff{efficiencies, k0, troe: ron::Troe{A, T3, T1, T2}} => f(&format_args!(", double[]({}), {}, {}, {}, {}, {})", efficiencies.iter().format(", "), fmt(k0), A, T3, T1, T2)),
					}}
				}))
			};
			write!(f, "double[]({}),\\\ndouble[]({}),\\\ndouble[][][]({}),\\\n{},\\\n{},\\\n{}",
									self.molar_masses.iter().format(", "),
									self.reduced_molar_masses.iter().format(", "),
									self.thermodynamics.iter().format_with(", ", |t, f|
										f(&format_args!("double[][]({})", t.0.iter().format_with(", ", |t, f| f(&format_args!("double[]({})", t.iter().format(", "))))))),
									fmt(&self.elementary), fmt(&self.three_body), fmt(&self.falloff))
    }
	}
	options.add_macro_definition("SYSTEM", Some(&system.to_string()));
	std::fs::write(std::path::PathBuf::from(std::env::var("OUT_DIR").unwrap()).join("main.spv"),
											shaderc::Compiler::new().unwrap().compile_into_spirv(include_str!("main.comp"), shaderc::ShaderKind::Compute, "main.comp", "main", Some(&options))?.as_binary_u8())?;
	Ok(())
}
