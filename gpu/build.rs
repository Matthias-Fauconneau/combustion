#![allow(incomplete_features)]#![feature(const_generics, const_evaluatable_checked, type_ascription, non_ascii_idents)]#![allow(non_snake_case)]

fn main() -> Result<(), Box<dyn std::error::Error>> {
	println!("cargo:rerun-if-changed=CH4+O2.ron");
	println!("cargo:rerun-if-changed=main.comp");
	let system = std::fs::read("CH4+O2.ron")?;
	use combustion::*;
	let System{molar_masses, thermodynamics, reactions, ..} = Simulation::<35>::new(&system)?.system;
	struct Falloff<const S: usize> where [(); S-1]: {
		reaction: Reaction<S>,
		efficiency: usize,
	};
	struct System<const S: usize> where [(); S-1]: {
		molar_masses: [f64; S],
		thermodynamics: [NASA7; S],
		elementary: Box<[Reaction<S>]>,
		three_body: Box<[Reaction<S>]>,
		pressure_modification: Box<[Reaction<S>]>,
		efficiencies: Box<[[f64; S]]>,
		falloff: Box<[Falloff<S>]>,
		transport_polynomials: TransportPolynomials<S>,
	}
	let mut efficiencies_set = Vec::new();
	let falloff = reactions.iter().filter_map(|reaction|
		if let Model::Falloff{efficiencies, ..} = reaction.model {
			if !efficiencies_set.contains(&efficiencies) { efficiencies_set.push(efficiencies); }
			Some(Falloff{reaction: *reaction, efficiency: efficiencies_set.iter().position(|&e| e == efficiencies).unwrap()})
		} else { None }
	).collect();

	let system = System{
		molar_masses,
		thermodynamics,
		elementary: reactions.iter().filter(|r| matches!(r.model, Model::Elementary{})).copied().collect(),
		three_body: reactions.iter().filter(|r| matches!(r.model, Model::ThreeBody{..})).copied().collect(),
		pressure_modification: reactions.iter().filter(|r| matches!(r.model, Model::PressureModification{..})).copied().collect(),
		efficiencies: efficiencies_set.into_boxed_slice(),
		falloff,
		transport_polyonmials
	};
	impl<const S: usize> std::fmt::Display for System<S> where [(); S-1]: {
		fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
			use itertools::Itertools;
			let fmt = |reactions:&[Reaction<S>]| {
				let variant = {use Model::*; match reactions.first().unwrap().model {Elementary => "Elementary", ThreeBody{..} => "ThreeBody", PressureModification{..} => "PressureModification", Falloff{..} => "Falloff"}};
				format!("{}[]({})", variant, reactions.iter().format_with(",\\\n", |r, f| {
					let fmt = |r:RateConstant| format!("RateConstant({}, {}, {})", r.log_preexponential_factor, r.temperature_exponent, r.activation_temperature);
					f(&format_args!("{}(Reaction(double[]({}), double[]({}), double[]({}), {}, {}, {}, {})",
						variant, r.reactants.iter().format(", "), r.products.iter().format(", "), r.net.iter().format(", "), r.Σreactants, r.Σproducts, r.Σnet, fmt(r.rate_constant)))?;
					{use Model::*; match r.model {
						Elementary => f(&")"),
						ThreeBody{efficiencies} => f(&format_args!(", double[]({}))", efficiencies.iter().format(", "))),
						PressureModification{efficiencies, k0} => f(&format_args!(", double[]({}), {})", efficiencies.iter().format(", "), fmt(k0))),
						Falloff{..} => unreachable!(), //Falloff{efficiencies, k0, troe: ron::Troe{A, T3, T1, T2}} => f(&format_args!(", double[]({}), {}, {}, {}, {}, {})", efficiencies.iter().format(", "), fmt(k0), A, T3, T1, T2)),
					}}
				}))
			};
			let fmt_falloff = |falloff:&[Falloff<S>]| {
				format!("Falloff[]({})", falloff.iter().format_with(",\\\n", |falloff, f| {
					let r = falloff.reaction;
					let fmt = |r:RateConstant| format!("RateConstant({}, {}, {})", r.log_preexponential_factor, r.temperature_exponent, r.activation_temperature);
					f(&format_args!("Falloff(Reaction(double[]({}), double[]({}), double[]({}), {}, {}, {}, {})",
						r.reactants.iter().format(", "), r.products.iter().format(", "), r.net.iter().format(", "), r.Σreactants, r.Σproducts, r.Σnet, fmt(r.rate_constant)))?;
					{use Model::*; match r.model {
						Falloff{efficiencies: _, k0, troe: ron::Troe{A, T3, T1, T2}} => f(&format_args!(", {}, {}, {}, {}, {}, {})", falloff.efficiency, fmt(k0), A, T3, T1, T2)),
						_ => unreachable!(),
					}}
				}))
			};
			write!(f, "double[]({}),\\\ndouble[][][]({}),\\\n{},\\\n{},\\\n{},\\\ndouble[][]({}),\\\n{}",
									self.molar_masses.iter().format(", "),
									self.thermodynamics.iter().format_with(", ", |t, f|
										f(&format_args!("double[][]({})", t.0.iter().format_with(", ", |t, f| f(&format_args!("double[]({})", t.iter().format(", "))))))),
									fmt(&self.elementary), fmt(&self.three_body), fmt(&self.pressure_modification), self.efficiencies.iter().format_with(", ",|e,f| f(&format_args!("double[]({})", e.iter().format(", ")))), fmt_falloff(&self.falloff))
    }
	}
	let mut options = shaderc::CompileOptions::new().unwrap();
	options.add_macro_definition("SPECIES", Some(&system.thermodynamics.len().to_string()));
	options.add_macro_definition("ELEMENTARY", Some(&system.elementary.len().to_string()));
	options.add_macro_definition("THREE_BODY", Some(&system.three_body.len().to_string()));
	options.add_macro_definition("PRESSURE_MODIFICATION", Some(&system.pressure_modification.len().to_string()));
	options.add_macro_definition("EFFICIENCIES", Some(&system.efficiencies.len().to_string()));
	options.add_macro_definition("FALLOFF", Some(&system.falloff.len().to_string()));
	let OUT_DIR = std::env::var("OUT_DIR").unwrap();
	std::fs::write(std::path::PathBuf::from(&OUT_DIR).join("model.h"), system.to_string())?;
	options.add_macro_definition("SYSTEM", Some(&system.to_string()));
	std::fs::write(std::path::PathBuf::from(OUT_DIR).join("main.spv"),
											shaderc::Compiler::new().unwrap().compile_into_spirv(include_str!("main.comp"), shaderc::ShaderKind::Compute, "main.comp", "main", Some(&options))?.as_binary_u8())?;
	Ok(())
}
