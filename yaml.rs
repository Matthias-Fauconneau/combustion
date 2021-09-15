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

use iter::{list, map};
/*use pest::Parser;
//('A'..'Z' | '0'..'9' | "-" | "*" | "(" | ")" | "," | "#")*
#[derive(pest_derive::Parser)]#[grammar_inline = r##"
	WHITESPACE = _{ " " }
	count = @{ '2'..'9' | '1'..'9' ~ ('0'..'9')+ }
	atom = { "H" | "O" | "C" | "AR" | "N" }
	specie = @{ ("A"|"I"|"P"|"S")? ~ (atom ~ count?)+ ~ ("V")?  }
	term = { count? ~ specie }
	side = { term ~ ("+" ~ term)* }
	equation = { side ~ ("=>"|"<=>") ~ side | side ~ ("+"~"M" | "(+M)") ~ "<=>" ~ side ~ ("+"~"M" | "(+M)") }
"##] struct EquationParser;
fn equation(equation: &str) -> [Map<&str, u8>; 2] {
	map(EquationParser::parse(Rule::equation, equation).unwrap().next().unwrap().into_inner(), |side|
		side.into_inner().map(|term|
			match &*list(term.into_inner()) {
				[specie] => (specie.as_str(), 1),
				[count, specie] => (specie.as_str(), count.as_str().parse().unwrap()),
				_ => unreachable!(),
			}
		)
		.fold(Map::new(), |mut map, (k, v)| { *map.entry(k).or_insert(0) += v; map })
	).into_vec().try_into().unwrap()
}*/
fn equation(equation: &str) -> [Map<&str, u8>; 2] {
	let sides: Box<[&str; 2]> = (if equation.contains("<=>") { list(equation.split("<=>")) } else { list(equation.split("=>")) }).try_into().unwrap();
	sides.map(|side|
		side.trim().trim_end_matches("(+ M)").trim_end_matches("+ M").trim_end_matches("(+M)").split(" + ").map(str::trim).map(|term| {
			if term.contains(' ') { let t:Box<[&str; 2]> = list(term.split_whitespace()).try_into().unwrap(); let [count, specie] = *t; (specie, count.parse().expect(term)) }
			else { (term, 1) }
		})
		.fold(Map::new(), |mut map, (k, v)| { *map.entry(k).or_insert(0) += v; map })
	)
}

pub use yaml::{Yaml, YamlLoader as Loader};
use combustion::model::*;

pub fn parse(yaml: &[Yaml]) -> Model {
	use std::str::FromStr;
	let data = &yaml[0];
	let units = &data["units"];
	assert!(units["length"].as_str()==Some("cm") && units["time"].as_str()==Some("s") && units["quantity"].as_str()==Some("mol"));
	let species = map(data["species"].as_vec().unwrap(), |specie| (specie["name"].as_str().unwrap(), Specie{
		composition: specie["composition"].as_hash().unwrap().iter().map(|(k,v)| (Element::from_str(k.as_str().unwrap()).unwrap(), v.as_i64().unwrap().try_into().unwrap())).collect(),
		thermodynamic: (|thermo:&Yaml| {
			let temperature_ranges = map(thermo["temperature-ranges"].as_vec().unwrap(), |limit| limit.as_f64().unwrap());
			let pieces = map(thermo["data"].as_vec().unwrap(), |piece| map(piece.as_vec().unwrap().iter(), |v| v.as_f64().unwrap()).into_vec().try_into().unwrap());
			fn has_duplicates<T:PartialEq>(slice:&[T]) -> bool { (1..slice.len()).any(|i| slice[i..].contains(&slice[i - 1])) }
			if {
				struct NASA7([f64; 7]);
				impl PartialEq for NASA7 { fn eq(&self, b: &Self) -> bool { self.0[0..6]==b.0[0..6] && f64::abs(self.0[6]-b.0[6])<3e-8 } }
				has_duplicates(&*map(&*pieces, |&x| NASA7(x)))
			} {
				eprintln!("{}: Merged duplicate splines", specie["name"].as_str().unwrap());
				assert!(pieces.len() == 2);
				let temperature_ranges: Box<_> = temperature_ranges.try_into().unwrap();
				let [min, _, max]: [f64; 3] = *temperature_ranges;
				NASA7{temperature_ranges: list([min, max]), pieces: list([pieces[0]])}
			} else {
				NASA7{temperature_ranges, pieces}
			}
		})(&specie["thermo"]),
		transport: (|transport:&Yaml| Transport{
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
	}));
	let reactions = map(data["reactions"].as_vec().unwrap().iter(), |reaction| {
		let equation: [Map<&str, u8>; 2] = equation(reaction["equation"].as_str().unwrap());
		let reactants: u8 = equation[0].iter().map(|(_,n)| n).sum();
		let reactants = reactants + if let Some("three-body") = reaction["type"].as_str() {1} else {0};
		let rate_constant = |rate_constant:&Yaml, concentration_cm3_unit_conversion_factor_exponent: u8| RateConstant{
			preexponential_factor: rate_constant["A"].as_f64().unwrap()*f64::powf(1e-6, concentration_cm3_unit_conversion_factor_exponent as f64),
			temperature_exponent: rate_constant["b"].as_f64().unwrap(),
			activation_temperature: {
				let Ea = rate_constant["Ea"].as_f64().unwrap();
				match units["activation-energy"].as_str().unwrap() {
					"K" => Ea,
					"cal/mol" => {
						const J_per_cal: f64 = 4.184;
						Ea*J_per_cal/R
					},
					_ => unimplemented!(),
				}
			}
		};
		Reaction{
			equation,
			rate_constant: rate_constant(Some(&reaction["rate-constant"]).filter(|v| !v.is_badvalue()).unwrap_or_else(|| Some(&reaction["high-P-rate-constant"]).filter(|v| !v.is_badvalue()).unwrap_or_else(|| {
				assert!(reaction["type"].as_str() == Some("pressure-dependent-Arrhenius"));
				assert!(reaction["rate-constants"].as_vec().unwrap().len() == 1);
				&reaction["rate-constants"].as_vec().unwrap()[0]
			})), reactants-1),
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
						T2: reaction["Troe"]["T2"].as_f64().unwrap_or(f64::INFINITY)
					}},
					Some("pressure-dependent-Arrhenius") => {
						assert!(reaction["rate-constants"].as_vec().unwrap().len() == 1);
						if reaction["equation"].as_str().unwrap().contains(" <=> ") { Elementary }
						else if reaction["equation"].as_str().unwrap().contains(" => ") { Irreversible }
						else { unimplemented!() }
					}
					_ => unimplemented!("{reaction:?}")
				}
			}
		}
	});
	let species = if false { species } else { // Easier Cantera comparison /!\ Breaks rates /!\
		let (active, inert) : (Vec<_>, _) = species.to_vec().into_iter().partition(|(specie,_)| reactions.iter().any(|Reaction{equation,..}| equation[0].get(specie).unwrap_or(&0) != equation[1].get(specie).unwrap_or(&0)));
		[&*active, &*inert].concat().into_boxed_slice()
	};
	let initial = list("H2 O2 CO CH2O CH4 C2H6 C2H4 C2H2 C3H6 pC3H4 aC3H4 1-C4H8 2-C4H8 i-C4H8 1,3-C4H6 C4H4 C4H2 C6H6 C7H8 H2O2 CH3OH CH3CHO C2H5OH CH3OCH3 C3H8 n-C4H10 i-C4H10".split(' ')); // CH3CO
	Model{
		state: State{volume: 1., temperature: 1000., pressure: 101325., amount_proportions: map(&*species, |(name,_)| (*name,
			if initial.contains(name) { 1. } else { 0. }))},
		species,
		reactions,
		time_step: 1e-7,
	}
}
