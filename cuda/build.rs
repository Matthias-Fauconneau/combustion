#![allow(incomplete_features, non_snake_case)]#![feature(const_generics, const_evaluatable_checked, bool_to_option, non_ascii_idents)]

use combustion::*;
struct System<const S: usize> where [(); S-1]: {
	molar_masses: [f64; S],
	thermodynamics: [NASA7; S],
	elementary: Box<[Reaction<S>]>,
	three_body: Box<[Reaction<S>]>,
	falloff: Box<[Reaction<S>]>,
}
impl<const S: usize> std::fmt::Display for System<S> where [(); S-1]: {
	fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
		use itertools::Itertools;
		let fmt = |reactions:&[Reaction<S>]| {
			let variant = {use Model::*; match reactions.first().unwrap().model {Elementary => "Elementary", ThreeBody{..} => "ThreeBody", Falloff{..} => "Falloff"}};
			format!("{{{}}}", reactions.iter().format_with(",\n", |r, f| {
				let fmt = |r:RateConstant| format!("RateConstant{{{}, {}, {}}}", r.log_preexponential_factor, r.temperature_exponent, r.activation_temperature);
				f(&format_args!("{}{{Reaction{{{{{}}}, {{{}}}, {{{}}}, {}, {}, {}, {}}}",
					variant, r.reactants.iter().format(", "), r.products.iter().format(", "), r.net.iter().format(", "), r.Σreactants, r.Σproducts, r.Σnet, fmt(r.rate_constant)))?;
				{use Model::*; match r.model {Elementary => f(&"}"),
					ThreeBody{efficiencies} => f(&format_args!(", {{{}}}}}", efficiencies.iter().format(", "))),
					Falloff{efficiencies, k0, troe: ron::Troe{A, T3, T1, T2}} => f(&format_args!(", {{{}}}, {}, {}, {}, {}, {}}}", efficiencies.iter().format(", "), fmt(k0), A, T3, T1, T2)),
				}}
			}))
		};
		write!(f, "{{{}}},\n{{{}}},\n{},\n{},\n{}",
								self.molar_masses.iter().format(", "),
								self.thermodynamics.iter().format_with(", ", |t, f|
									f(&format_args!("{{{}}}", t.0.iter().format_with(", ", |t, f| f(&format_args!("{{{}}}\n", t.iter().format(", "))))))),
								fmt(&self.elementary), fmt(&self.three_body), fmt(&self.falloff))
	}
}

#[fehler::throws(anyhow::Error)] fn main() {
	println!("cargo:rerun-if-changed=H2+O2.ron");
	println!("cargo:rerun-if-changed=main.cu");
	let system = std::fs::read("H2+O2.ron")?;
	type Simulation<'t> = combustion::Simulation::<'t, 9>;
	let combustion::System{molar_masses, thermodynamics, reactions, ..} = Simulation::new(&system)?.system;
	let system = System{
		molar_masses,
		thermodynamics,
		elementary: reactions.iter().filter(|r| matches!(r.model, Model::Elementary{})).copied().collect(),
		three_body: reactions.iter().filter(|r| matches!(r.model, Model::ThreeBody{..})).copied().collect(),
		falloff: reactions.iter().filter(|r| matches!(r.model, Model::Falloff{..})).copied().collect()
	};
	let OUT_DIR = std::env::var("OUT_DIR").unwrap();
	std::fs::write(std::path::PathBuf::from(&OUT_DIR).join("system.h"), system.to_string())?;
	std::process::Command::new("/opt/cuda/bin/nvcc").args(&["-I", &OUT_DIR,
		&format!("-DSPECIES={}",system.thermodynamics.len()), &format!("-DELEMENTARY={}", system.elementary.len()), &format!("-DTHREE_BODY={}", system.three_body.len()), &format!("-DFALLOFF={}", system.falloff.len()),
		"--ptx","main.cu","-o",&format!("{}/main.ptx",OUT_DIR)
	]).spawn()?.wait()?.success().then_some(()).ok_or(anyhow::anyhow!(""))?;
}
