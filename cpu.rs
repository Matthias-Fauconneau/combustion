#![feature(format_args_capture,iter_is_partitioned)]#![allow(non_snake_case)]
pub fn compile(f: &ast::Function) -> impl Fn(&[f32]) -> Box<[f32]> {
	let output = f.output.len();
	let f = ir::assemble({let clif = ir::compile(&f); std::fs::write("/var/tmp/clif", &clif.to_string()).unwrap(); clif});
	move |input| { let mut output = vec![0.; output].into_boxed_slice(); f(input, &mut output); output }
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
	let model = yaml_model::Loader::load_from_str(std::str::from_utf8(&std::fs::read(std::env::args().skip(1).next().unwrap())?)?)?;
	let model = yaml_model::parse(&model);
	use chemistry::*;
	let (ref species_names, ref species) = Species::new(&model.species);
	let reactions = map(&*model.reactions, |r| Reaction::new(species_names, r));
	let ref state = initial_state(&model);
	use {iter::map, itertools::Itertools, ast::let_};
	if true {
		let rates = compile(&reaction::rates(&species.thermodynamics, &reactions));
		assert!(state.volume == 1.);
		let State{temperature: T, pressure_R, amounts, ..} = state;
		let total_amount = amounts.iter().sum::<f64>();
		let active_amounts = &amounts[0..amounts.len()-1];
		let_!{ [dtT_T, rates @ ..] = &*rates(&iter::box_([*pressure_R as f32, total_amount as f32, *T as f32].into_iter().chain(active_amounts.iter().map(|&n| n as f32)))) => {
		println!("{}, dtT_T: {dtT_T:.3e}", rates.iter().zip(&**species_names).format_with(", ", |(rate, name), f| f(&format!("{name}: {rate:.0}").to_string())));
		}}
	}
	#[cfg(feature="transport")] {
		let transport = wrap(transport::properties::<4>(&species));
		let State{temperature: T, pressure_R, amounts, ..} = state;
		let total_amount = amounts.iter().sum::<f64>();
		let active_amounts = &amounts[0..amounts.len()-1];
		let_!{[viscosity, thermal_conductivity, mixture_diffusion_coefficients @ ..] = &*transport(&[&[*pressure_R,total_amount, *T], active_amounts].concat()) => {
		eprintln!("μ: {viscosity:.4e}, λ: {thermal_conductivity:.4}, D: {:.4e}", mixture_diffusion_coefficients.iter().format(" "));
		}}
	}
	Ok(())
}
