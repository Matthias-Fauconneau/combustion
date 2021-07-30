#![feature(destructuring_assignment,format_args_capture,trait_alias,default_free_fn,in_band_lifetimes,unboxed_closures,fn_traits)]#![allow(non_snake_case,non_upper_case_globals)]
mod yaml; mod device;
use {anyhow::Result, iter::{list, map}, combustion::*, device::*};
fn main() -> Result<()> {
	let model = yaml::Loader::load_from_str(std::str::from_utf8(&std::fs::read(std::env::args().skip(1).next().unwrap())?)?)?;
	let model = yaml::parse(&model);
	let (ref species_names, ref species, _, reactions, ref _state) = new(&model);
	let rates = reaction::rates(&species.thermodynamics, &reactions);
	#[cfg(not(feature="f32"))] type T = f64;
	#[cfg(feature="f32")] type T = f32;
	let rates = assemble::<T>(rates, 1);

	let State{pressure_R, temperature, ..} = initial_state(&model);
	fn parse(s:&str) -> std::collections::HashMap<&str,f64> {
		s.split(",").map(|e| { let [key, value] = {let b:Box<[_;2]> = e.split(":").collect::<Box<_>>().try_into().unwrap(); *b}; (key, value.parse().unwrap()) }).collect()
	}
	let amounts = map(&**species_names, |s| *parse("H2:2,O2:1,N2:2").get(s).unwrap_or(&0.));
	let total_amount = amounts.iter().sum::<f64>();
	let state = [&[temperature], &amounts[..amounts.len()-1]].concat();
	let mut input = list([list([total_amount as _])].into_iter().chain(state.iter().map(|&v| list([v as _]))));
	let mut evaluations = 0;
	let f = |u: &[f64], f_u: &mut [f64]| {
		//use itertools::Itertools;
		//println!("{:3} {:.2e}", "u", u.iter().format(", "));
		assert!(u[0]>200.);
		assert!(u.iter().all(|u| u.is_finite()));
		for (input, &u) in input[1..].iter_mut().zip(u) { input[0] = u as _; }
		let output = rates(&[pressure_R as _], &map(&*input, |input| &**input)).unwrap();
		assert!(output.iter().all(|u| u[0].is_finite()), "{input:?} {output:?}");
		for (f_u, output) in f_u.iter_mut().zip(&*output) { *f_u = all_same(output, 1) as _; }
		evaluations += 1;
		//println!("{:3} {:.2e}", "f_u", f_u.iter().format(", "));
		true
	};

	let mut cvode = cvode::CVODE::new(/*relative_tolerance:*/ 1e-4, /*absolute_tolerance:*/ 1e-7, &state);
	let plot_min_time = 0.1;
	let (mut time, mut state) = (0., &*state);
	let start = std::time::Instant::now();
	while time < plot_min_time {
		let next_time = time + model.time_step;
		while time < next_time { (time, state) = cvode.step(&f, next_time, state); }
		dbg!(time/plot_min_time, time, state[0]);
		//assert!(state[0]<1500.);
	}
	//println!("T {}", state[0]);
	let mut min = 0f64;
	let values = map(0..(0.2/model.time_step) as usize, |_| if let [temperature, active_amounts@..] = state {
		let value = ((time-plot_min_time)*1e3, vec![vec![*temperature].into_boxed_slice(), map(active_amounts, |v| v/active_amounts.iter().sum::<f64>())].into_boxed_slice());
		let next_time = time + model.time_step;
		while time < next_time { (time, state) = cvode.step(&f, next_time, state); }
		min = min.min(state[1..].iter().copied().reduce(f64::min).unwrap());
		value
	} else { unreachable!() });
	let end = std::time::Instant::now();
	let time = (end-start).as_secs_f64();
	println!("T {}", state[0]);
	println!("{evaluations} / {time}s = {}", (evaluations as f64)/time);
	//println!("{:.0e}", min);
	let active = map(0..species.len()-1, |k| reactions.iter().any(|Reaction{net,..}| net[k] != 0)).iter().position(|active| !active).unwrap_or(species.len()-1);
	let key = map((0..active).map(|k| values.iter().map(move |(_, sets)| sets[1][k])), |Y| Y.reduce(f64::max).unwrap());
	let mut s = map(0..active, |i| i);
	let (_, _, select) = s.select_nth_unstable_by(species.len()-6, |&a,&b| key[a].partial_cmp(&key[b]).unwrap());
	let species_names = map(&*select, |&k| species_names[k]);
	let values = map(Vec::from(values).into_iter(), |(t, sets)| {
		let [sets_0, sets_1] = *std::convert::TryInto::<Box<[_;2]>>::try_into(sets).unwrap();
		(t, vec![sets_0, map(&*select, |&k| sets_1[k])].into_boxed_slice())
	});
	let mut plot = ui::plot::Plot::new(vec![&["T"] as &[_], &species_names].into_boxed_slice(), values.to_vec());
	if true { Ok(ui::app::run(plot)?) }
	else {
		let mut target = image::Image::zero((3840, 2160).into());
		ui::widget::Widget::paint(&mut plot, &mut target.as_mut())?;
		pub fn as_bytes<T>(slice: &[T]) -> &[u8] { unsafe{std::slice::from_raw_parts(slice.as_ptr() as *const u8, slice.len() * std::mem::size_of::<T>())} }
		Ok(png::save_buffer("/var/tmp/image.png", as_bytes(&target.data), target.size.x, target.size.y, png::ColorType::Rgba8)?)
	}
}
