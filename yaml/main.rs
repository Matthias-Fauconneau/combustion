#![feature(array_methods, non_ascii_idents)]#![allow(non_snake_case)]//#![recursion_limit="8"]

struct Pretty<T>(T);
impl std::fmt::Display for Pretty<&f64> {
	fn fmt(&self, fmt: &mut std::fmt::Formatter) -> std::fmt::Result {
		let value = *self.0;
		if value == f64::floor(value) {
			if value >= 1e4 { write!(fmt, "{:e}", self.0) }
			else { (value as i64).fmt(fmt) }
		}
		else { float_pretty_print::PrettyPrintFloat(value).fmt(fmt) }
	}
}
impl std::fmt::Display for Pretty<&u8> {
	fn fmt(&self, fmt: &mut std::fmt::Formatter) -> std::fmt::Result { self.0.fmt(fmt) }
}

use {fehler::throws, anyhow::Error};

use combustion::ron::*;
#[throws(std::fmt::Error)] fn to_string(system: &System) -> String {
	let mut o = String::new();
	use std::fmt::Write;
	writeln!(o, "#![enable(unwrap_newtypes)]")?;
	writeln!(o, "(")?;
	writeln!(o, "time_step: {},", Pretty(&system.time_step))?;
	//let inline = || ::ron::ser::PrettyConfig::new().with_new_line(" ".to_owned()).with_indentor("".to_owned());
	//writeln!(o, "state: {}", to_string_pretty(&system.state, inline())?);
	/*impl<V> std::fmt::Display for Pretty<&Map<&str, V>> where Pretty<&V>: std::fmt::Display {
		fn fmt(&self, fmt: &mut std::fmt::Formatter) -> std::fmt::Result { write!(fmt, "{{{}}}", self.0.iter().format_with(", ", |(k,v), f| f(&format_args!("\"{}\": {}", k, Pretty(v))))) }
	}*/
	use itertools::Itertools;
	impl std::fmt::Display for Pretty<&Map<&str, f64>> {
		fn fmt(&self, fmt: &mut std::fmt::Formatter) -> std::fmt::Result { write!(fmt, "{{{}}}", self.0.iter().format_with(", ", |(k,v), f| f(&format_args!("\"{}\": {}", k, Pretty(v))))) }
	}
	impl std::fmt::Display for Pretty<&Map<&str, u8>> {
		fn fmt(&self, fmt: &mut std::fmt::Formatter) -> std::fmt::Result { write!(fmt, "{{{}}}", self.0.iter().format_with(", ", |(k,v), f| f(&format_args!("\"{}\": {}", k, Pretty(v))))) }
	}
	writeln!(o, "state: (temperature: {}, pressure: {}, mole_proportions: {}),", system.state.temperature, system.state.pressure, Pretty(&system.state.mole_proportions))?;
	writeln!(o, "species: {{")?;
	for (name, Specie{composition, thermodynamic:NASA7{temperature_ranges, pieces}, transport}) in &system.species {
		writeln!(o, "\"{}\": (", name)?;
		writeln!(o, "\tcomposition: {:?},", composition)?;
		writeln!(o, "\tthermodynamic: (")?;
		/*impl<'t, T: 't> std::fmt::Display for Pretty<&[T]> where Pretty<&'t T>: std::fmt::Display {
			fn fmt(&self, fmt: &mut std::fmt::Formatter) -> std::fmt::Result { write!(fmt, "[{}]", self.0.iter().format_with(", ", |v, f| f(&format_args!("{}", Pretty(v))))) }
		}*/
		impl std::fmt::Display for Pretty<&[f64]> {
			fn fmt(&self, fmt: &mut std::fmt::Formatter) -> std::fmt::Result { write!(fmt, "[{}]", self.0.iter().format_with(", ", |v, f| f(&format_args!("{}", Pretty(v))))) }
		}
		impl<T, const N: usize> std::fmt::Display for Pretty<&[T; N]> where for<'t> Pretty<&'t T>: std::fmt::Display {
			fn fmt(&self, fmt: &mut std::fmt::Formatter) -> std::fmt::Result { write!(fmt, "({})", self.0.iter().format_with(", ", |v, f| f(&format_args!("{}", Pretty(v))))) }
		}
		/*impl<'t, T: 't, const N: usize> std::fmt::Display for Pretty<&[T; N]> where Pretty<&'t T>: std::fmt::Display {
			fn fmt(&self, fmt: &mut std::fmt::Formatter) -> std::fmt::Result { write!(fmt, "({})", self.0.iter().format_with(", ", |v, f| f(&format_args!("{}", Pretty(v))))) }
		}*/
		writeln!(o, "\t\ttemperature_ranges: {},", Pretty(temperature_ranges.as_ref()))?;
		writeln!(o, "\t\t\tpieces: [")?;
		for piece in pieces.iter() { writeln!(o, "\t\t\t{},", Pretty(piece))?; }
		writeln!(o, "\t\t],")?;
		writeln!(o, "\t),")?;
		//writeln!(o, "\ttransport: {}", to_string_pretty(&transport, inline())?);
		//impl<'t, T: IntoIterator<Item=(&'t str, f64)>> std::fmt::Display for Pretty<T> {
		struct Fields<'t, const N: usize>(&'t [(&'t str, f64); N]);
		impl<const N: usize> std::fmt::Display for Fields<'_, N> {
			#[throws(std::fmt::Error)] fn fmt(&self, fmt: &mut std::fmt::Formatter) {
				let nonzero = self.0.into_iter().filter(|(_,v)| *v != 0.).collect::<Box<_>>();
				//if nonzero.len() == 0 { return; }
				write!(fmt, "({})", nonzero.iter().format_with(", ", |(k,v), f| f(&format_args!("{}: {}", k, Pretty(v)))))?
			}
		}
		writeln!(o, "\ttransport: (well_depth_K: {}, diameter_A: {}, geometry: {})", transport.well_depth_K, transport.diameter_Å, {use Geometry::*; match transport.geometry {
			Atom => format!("Atom"),
			Linear{polarizability_Å3, rotational_relaxation} => format!("Linear{}", Fields(&[("polarizability_A3",polarizability_Å3),("rotational_relaxation",rotational_relaxation)])),
			Nonlinear{polarizability_Å3, rotational_relaxation, permanent_dipole_moment_Debye} =>
				format!("Nonlinear{}", Fields(&[("polarizability_A3",polarizability_Å3),("rotational_relaxation",rotational_relaxation),("permanent_dipole_moment_Debye",permanent_dipole_moment_Debye)])),
		}})?;
		writeln!(o, "),")?;
	}
	writeln!(o, "}},")?;
	writeln!(o, "reactions: [")?;
	for Reaction{equation, rate_constant, model} in system.reactions.iter() {
		impl std::fmt::Display for Pretty<&RateConstant> {
			fn fmt(&self, fmt: &mut std::fmt::Formatter) -> std::fmt::Result { write!(fmt, "(A: {}, beta: {}, Ea: {})", Pretty(&self.0.preexponential_factor), Pretty(&self.0.temperature_exponent), Pretty(&self.0.activation_energy)) }
		}
		write!(o, "(equation: {}, rate_constant: {}, model: ", Pretty(equation), Pretty(rate_constant))?;
		if let Model::Falloff{..} = model { write!(o, "\n\t ")?; }
		use Model::*; match model {
			Elementary => write!(o, "Elementary"),
			ThreeBody{efficiencies} =>
				write!(o, "ThreeBody(efficiencies: {})", Pretty(efficiencies)),
			PressureModification{k0, efficiencies} =>
				write!(o, "PressureModification(k0: {}, efficiencies: {})", Pretty(k0), Pretty(efficiencies)),
			Falloff{k0, troe: Troe{A,T3,T1,T2}, efficiencies} =>
				write!(o, "Falloff(k0: {}, troe: {}, efficiencies: {})", Pretty(k0), format!("(A: {}, T3: {}, T1: {}, T2: {})", Pretty(A), Pretty(T3), Pretty(T1), Pretty(T2)), Pretty(efficiencies)),
		}?;
		writeln!(o, "),")?;
	}
	writeln!(o, "]")?;
	write!(o, ")")?;
	o
}

#[throws(Error)] fn main() {
	let system = std::fs::read("CH4+O2.ron")?;
	let mut system: System  = ::ron::de::from_bytes(&system)?;
	let yaml = yaml_rust::YamlLoader::load_from_str(std::str::from_utf8(&std::fs::read("/usr/share/cantera/data/gri30.yaml")?)?).unwrap();
	for yaml in yaml[0]["species"].as_vec().unwrap() {
		if let Some(specie) = system.species.get_mut(yaml["name"].as_str().unwrap()) {
			use std::convert::TryInto;
			specie.thermodynamic.pieces = yaml["thermo"]["data"].as_vec().unwrap().iter().map(|piece| piece.as_vec().unwrap().iter().map(|v| v.as_f64().unwrap()).collect::<Box<_>>().as_ref().try_into().unwrap()).collect();
		}
	}
	//writeln!(o, "{}", ::ron::ser::to_string_pretty(system, ::ron::ser::PrettyConfig{scientific_notation: true, ..default()})?)
	//pretty_assertions::assert_eq!(::ron::de::from_str::<System>(&to_string(&system)?)?, system);
	println!("{}", to_string(&system)?);
}
