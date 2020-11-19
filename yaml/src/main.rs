#[fehler::throws(anyhow::Error)] fn main() {
	let system = std::fs::read("H2+O2.ron")?;
	use combustion::*;
	let ron::System{mut species, ..} = ::ron::de::from_bytes(&system)?;
	let yaml = yaml_rust::YamlLoader::load_from_str(std::str::from_utf8(&std::fs::read("/usr/share/cantera/data/gri30.yaml")?)?).unwrap();
	for yaml in yaml[0]["species"].as_vec().unwrap() {
		if let Some(specie) = species.get_mut(yaml["name"].as_str().unwrap()) {
			use std::convert::TryInto;
			specie.thermodynamic.pieces = yaml["thermo"]["data"].as_vec().unwrap().iter().map(|piece| piece.as_vec().unwrap().iter().map(|v| v.as_f64().unwrap()).collect::<Box<_>>().as_ref().try_into().unwrap()).collect();
		}
	}
	for (name, Specie{composition, thermodynamic:ron::NASA7{temperature_ranges, pieces}, ..}) in species {
		println!("\"{}\": (", name);
		println!("  composition: {:?},", composition);
		println!("  thermodynamic: (");
		println!("    temperature_ranges: {:?},", temperature_ranges);
		println!("    pieces: [");
		use itertools::Itertools;
		for piece in pieces.iter() { println!("      ({}),", piece.iter().map(|&x| float_pretty_print::PrettyPrintFloat(x)).format(", ")); }
		println!("    ],");
		println!("  ),");
		println!("),");
	}
}
