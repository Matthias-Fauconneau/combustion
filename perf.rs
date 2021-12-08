#![feature(trait_alias,destructuring_assignment,default_free_fn,let_else,iter_zip,unboxed_closures,fn_traits,box_patterns)]
#![allow(non_snake_case,non_upper_case_globals)]
mod yaml; mod device;
use {anyhow::Result, combustion::*, device::*};
fn main() -> Result<()> {
	let path = std::env::args().skip(1).next().unwrap_or("gri30".to_string());
	let path = if std::path::Path::new(&path).exists() { path } else { format!("/usr/local/share/cantera/data/{path}.yaml") };
	let model = yaml::Loader::load_from_str(std::str::from_utf8(&std::fs::read(&path).expect(&path))?)?;
	let model = yaml::parse(&model);
	let (species_names, ref species, _active, reactions, State{pressure_R, volume, temperature, ref amounts}) = new(&model);
	let rates = reaction::rates(&species.molar_mass, &species.thermodynamics, &reactions, &species_names);
	let block_size = 1;
	let rates = with_repetitive_input(assemble(rates, block_size), 1*block_size);
    let nonbulk_amounts = &amounts[0..amounts.len()-1];
    rates(&[pressure_R, 1./pressure_R], &([[temperature, volume].as_slice(), nonbulk_amounts].concat()))?;
    Ok(())
}