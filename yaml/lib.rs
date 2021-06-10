#![feature(default_free_fn,array_map)]#![allow(non_snake_case)]
use std::default::default;
fn box_<T>(t: T) -> Box<T> { Box::new(t) }

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

use combustion::model::*;

extern crate pest;
#[macro_use]
extern crate pest_derive;
use pest::Parser;
#[derive(Parser)]#[grammar_inline = r#"
	WHITESPACE = _{ " " }
	atom = { "H" | "O" | "C" | "AR" | "N" }
	count = @{ '2'..'9' | '1'..'9' ~ ('0'..'9')+ }
	specie = @{ (atom ~ count?)+ ~ "(S)"? }
	term = { count? ~ specie }
	side = { term ~ ("+" ~ term)* }
	equation = { side ~ ("=>"|"<=>") ~ side | side ~ ("+"~"M" | "(+M)") ~ "<=>" ~ side ~ ("+"~"M" | "(+M)") }
"#] struct EquationParser;
fn equation(equation: &str) -> [Map<&str, u8>; 2] {
	use {std::convert::TryInto, std::str::FromStr};
	EquationParser::parse(Rule::equation, equation).unwrap().next().unwrap().into_inner().map(|side|
		side.into_inner().map(|term|
			match term.into_inner().collect::<Box<_>>().as_ref() {
				[specie] => (specie.as_str(), 1),
				[count, specie] => (specie.as_str(), u8::from_str(count.as_str()).unwrap()),
				_ => unreachable!(),
			}
		)
		.fold(Map::new(), |mut map, (k, v)| { *map.entry(k).or_insert(0) += v; map })
	).collect::<Vec<_>>().try_into().unwrap()
}

#[test] fn test() {
	use std::iter::FromIterator;
	assert_eq!(equation("CH2 + CH2 <=> H2 + C2H2"), [Map::from_iter([("CH2",2)].iter().copied()), Map::from_iter([("H2",1),("C2H2",1)].iter().copied())]);
}

pub use yaml_rust::YamlLoader as Loader;

#[throws] pub fn parse(yaml: &[yaml_rust::Yaml]) -> Model {
	use {std::convert::TryInto, std::str::FromStr};
	let data = &yaml[0];
	let species: Map<&str, Specie> = data["species"].as_vec().unwrap().iter().map(|specie| (specie["name"].as_str().unwrap(), Specie{
		composition: specie["composition"].as_hash().unwrap().iter().map(|(k,v)| (Element::from_str(k.as_str().unwrap()).unwrap(), v.as_i64().unwrap().try_into().unwrap())).collect(),
		thermodynamic: (|thermo:&yaml_rust::Yaml| NASA7{
			temperature_ranges: thermo["temperature-ranges"].as_vec().unwrap().iter().map(|limit| limit.as_f64().unwrap()).collect(),
			pieces: thermo["data"].as_vec().unwrap().iter().map(|piece| piece.as_vec().unwrap().iter().map(|v| v.as_f64().unwrap()).collect::<Box<_>>().as_ref().try_into().unwrap()).collect()
		})(&specie["thermo"]),
		transport: (|transport:&yaml_rust::Yaml| Transport{
			well_depth_K: transport["well-depth"].as_f64().unwrap(),
			diameter_Å: transport["diameter"].as_f64().unwrap(),
			geometry: {use Geometry::*; match transport["geometry"].as_str().unwrap() {
				"atom" => Atom,
				"linear" => Linear{
					polarizability_Å3: transport["polarizability"].as_f64().unwrap_or(0.),
					rotational_relaxation: transport["rotational-relaxation"].as_f64().unwrap_or(0.),
				},
				"nonlinear" => Nonlinear{
					polarizability_Å3: transport["polarizability"].as_f64().unwrap_or(0.),
					rotational_relaxation: transport["rotational-relaxation"].as_f64().unwrap_or(0.),
					permanent_dipole_moment_Debye: transport["dipole"].as_f64().unwrap_or(0.),
				},
				_ => unimplemented!()
		}}
		})(&specie["transport"])
	})).collect();
	//let species = {let mut species = species; species.insert("AR", species.remove("AR").unwrap()); species};
	let reactions : Box<_> = data["reactions"].as_vec().unwrap().iter().map(|reaction| {
		let equation: [Map<&str, u8>; 2] = equation(reaction["equation"].as_str().unwrap());
		let reactants: u8 = equation[0].iter().map(|(_,n)| n).sum();
		let reactants = reactants + if let Some("three-body") = reaction["type"].as_str() {1} else {0};
		let rate_constant = |rate_constant:&yaml_rust::Yaml, concentration_cm3_unit_conversion_factor_exponent: u8| RateConstant{
			preexponential_factor: rate_constant["A"].as_f64().unwrap()*f64::powf(1e-6, concentration_cm3_unit_conversion_factor_exponent as f64),
			temperature_exponent: rate_constant["b"].as_f64().unwrap(),
			activation_energy: rate_constant["Ea"].as_f64().unwrap(),
		};
		Reaction{
			equation,
			rate_constant: rate_constant(Some(&reaction["rate-constant"]).filter(|v| !v.is_badvalue()).unwrap_or(&reaction["high-P-rate-constant"]), reactants-1),
			model: {
				let efficiencies = reaction["efficiencies"].as_hash().map(|h| h.iter().map(|(k,v)| (k.as_str().unwrap(), v.as_f64().unwrap())).collect()).unwrap_or_default();
				use ReactionModel::*;
				match reaction["type"].as_str() {
					None|Some("elementary") => {
						if reaction["equation"].as_str().unwrap().contains(" <=> ") { Elementary }
						else if reaction["equation"].as_str().unwrap().contains(" => ") { Irreversible }
						else { unimplemented!() }
					},
					Some("three-body") => ThreeBody{efficiencies},
					// Pr = c x k0 / k∞ [1 = mol/m³ x [k0] / [k∞]] => [k0] = [k∞]/(mol/m³) = ((mol/m³)^(-r))/s => k0[m] = k0[cm] / (1e-2³)^(-r)
					Some("falloff") if reaction["Troe"].is_badvalue() => PressureModification{efficiencies, k0: rate_constant(&reaction["low-P-rate-constant"], reactants)},
					Some("falloff") => Falloff{efficiencies, k0: rate_constant(&reaction["low-P-rate-constant"], reactants), troe: Troe{
						A: reaction["Troe"]["A"].as_f64().unwrap(),
						T1: reaction["Troe"]["T1"].as_f64().unwrap(),
						T3: reaction["Troe"]["T3"].as_f64().unwrap(),
						T2: reaction["Troe"]["T2"].as_f64().unwrap_or(0.)
					}},
					_ => unimplemented!()
				}
			}
		}
	}).collect();
	let species_names : Box<_> = species.iter().map(|(name,_)| *name).collect();
	let species = if species_names.iter().map(|s| reactions.iter().all(|Reaction{equation:[a,b],..}| *a.get(s).unwrap_or(&0) as i8 == *b.get(s).unwrap_or(&0) as i8))
													.position(|inert| inert == true).filter(|&s| s == species.len()-1).is_some() { species } else {
		let mut species = species;
		species.insert("INERT", Specie{composition: default(), thermodynamic: NASA7{temperature_ranges: box_([]), pieces: box_([])},
		                                                     transport: Transport{well_depth_K: 0., diameter_Å: 0.,geometry: Geometry::Atom}});
		species
	};
	Model{
		state: State{volume: 1., temperature: 1000., pressure: 101325., amount_proportions: species.iter().map(|(&name,_)| (name, 1.)).collect()},
		species,
		reactions,
		time_step: 1e-8,
	}
}

#[throws(std::fmt::Error)] pub fn to_string(model: &Model) -> String {
	let mut o = String::new();
	use std::fmt::Write;
	writeln!(o, "#![enable(unwrap_newtypes)]")?;
	writeln!(o, "(")?;
	writeln!(o, "time_step: {},", Pretty(&model.time_step))?;
	use itertools::Itertools;
	impl std::fmt::Display for Pretty<&Map<&str, f64>> {
		fn fmt(&self, fmt: &mut std::fmt::Formatter) -> std::fmt::Result { write!(fmt, "{{{}}}", self.0.iter().format_with(", ", |(k,v), f| f(&format_args!("\"{}\": {}", k, Pretty(v))))) }
	}
	impl std::fmt::Display for Pretty<&Map<&str, u8>> {
		fn fmt(&self, fmt: &mut std::fmt::Formatter) -> std::fmt::Result { write!(fmt, "{{{}}}", self.0.iter().format_with(", ", |(k,v), f| f(&format_args!("\"{}\": {}", k, Pretty(v))))) }
	}
	writeln!(o, "state: (temperature: {}, pressure: {}, volume: {}, amount_proportions: {}),",
		model.state.temperature, model.state.pressure, model.state.volume, Pretty(&model.state.amount_proportions))?;
	writeln!(o, "species: {{")?;
	for (name, Specie{composition, thermodynamic:NASA7{temperature_ranges, pieces}, transport}) in &model.species {
		writeln!(o, "\"{}\": (", name)?;
		writeln!(o, "\tcomposition: {:?},", composition)?;
		writeln!(o, "\tthermodynamic: (")?;
		impl std::fmt::Display for Pretty<&[f64]> {
			fn fmt(&self, fmt: &mut std::fmt::Formatter) -> std::fmt::Result { write!(fmt, "[{}]", self.0.iter().format_with(", ", |v, f| f(&format_args!("{}", Pretty(v))))) }
		}
		impl<T, const N: usize> std::fmt::Display for Pretty<&[T; N]> where for<'t> Pretty<&'t T>: std::fmt::Display {
			fn fmt(&self, fmt: &mut std::fmt::Formatter) -> std::fmt::Result { write!(fmt, "({})", self.0.iter().format_with(", ", |v, f| f(&format_args!("{}", Pretty(v))))) }
		}
		writeln!(o, "\t\ttemperature_ranges: {},", Pretty(temperature_ranges.as_ref()))?;
		writeln!(o, "\t\t\tpieces: [")?;
		for piece in pieces.iter() { writeln!(o, "\t\t\t{},", Pretty(piece))?; }
		writeln!(o, "\t\t],")?;
		writeln!(o, "\t),")?;
		struct Fields<'t, const N: usize>(&'t [(&'t str, f64); N]);
		impl<const N: usize> std::fmt::Display for Fields<'_, N> {
			#[throws(std::fmt::Error)] fn fmt(&self, fmt: &mut std::fmt::Formatter) {
				let nonzero = self.0.into_iter().filter(|(_,v)| *v != 0.).collect::<Box<_>>();
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
	for Reaction{equation, rate_constant, model} in model.reactions.iter() {
		impl std::fmt::Display for Pretty<&RateConstant> {
			fn fmt(&self, fmt: &mut std::fmt::Formatter) -> std::fmt::Result { write!(fmt, "(A: {}, beta: {}, Ea: {})", Pretty(&self.0.preexponential_factor), Pretty(&self.0.temperature_exponent), Pretty(&self.0.activation_energy)) }
		}
		write!(o, "(equation: {}, rate_constant: {}, model: ", Pretty(equation), Pretty(rate_constant))?;
		if let ReactionModel::Falloff{..} = model { write!(o, "\n\t ")?; }
		use ReactionModel::*; match model {
			Elementary => write!(o, "Elementary"),
			Irreversible => write!(o, "Irreversible"),
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
