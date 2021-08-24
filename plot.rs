#![feature(format_args_capture,trait_alias,destructuring_assignment,default_free_fn)]#![allow(non_snake_case,non_upper_case_globals)]
mod yaml; mod device;
use {anyhow::Result, iter::map, combustion::*, device::*};
fn main() -> Result<()> {
	let path = std::env::args().skip(1).next().unwrap_or("LiDryer".to_string());
	let path = if std::path::Path::new(&path).exists() { path } else { format!("/usr/share/cantera/data/{path}.yaml") };
	let model = yaml::Loader::load_from_str(std::str::from_utf8(&std::fs::read(&path).expect(&path))?)?;
	let model = yaml::parse(&model);
	let (species_names, ref species, active, reactions, ref _state) = new(&model);
	let rates = reaction::rates(&species.molar_mass, &species.thermodynamics, &reactions, &species_names);
	#[cfg(not(feature="f32"))] type T = f64;
	#[cfg(feature="f32")] type T = f32;
	let rates = with_repetitive_input(assemble::<T>(rates, 1), 1);

	let State{pressure_R, volume, temperature, ..} = initial_state(&model);
	fn parse(s:&str) -> std::collections::HashMap<&str,f64> {
		s.split(",").map(|e| { let [key, value] = {let b:Box<[_;2]> = e.split(":").collect::<Box<_>>().try_into().unwrap(); *b}; (key, value.parse().unwrap()) }).collect()
	}
	let amounts = map(&*species_names, |s| *parse("H2:2,O2:1,N2:2").get(s).unwrap_or(&0.));
	let state = [&[temperature, volume], &amounts[..active]].concat();
	let mut evaluations = 0;
	let f = |u: &[f64], f_u: &mut [f64]| {
		//use itertools::Itertools;
		//println!("{:3} {:.2e}", "u", u.iter().format(", "));
		assert!(u[0]>200.);
		assert!(u.iter().all(|u| u.is_finite()));
		f_u.copy_from_slice(&rates(&[pressure_R as _, (1./pressure_R) as _], u).unwrap());
		assert!(f_u.iter().all(|u| u.is_finite()), "{u:?} {f_u:?}");
		evaluations += 1;
		//println!("{:3} {:.2e}", "f_u", f_u.iter().format(", "));
		true
	};

	let mut cvode = cvode::CVODE::new(/*relative_tolerance:*/ 1e-4, /*absolute_tolerance:*/ 1e-7, &state);
	let plot_min_time = 0.00005;
	let (mut time, mut state) = (0., &*state);
	let start = std::time::Instant::now();
	while time < plot_min_time {
		let next_time = time + model.time_step;
		while time < next_time { (time, state) = cvode.step(&f, next_time, state); }
		//dbg!(time/plot_min_time, time, state[0]);
		//assert!(state[0]<1500.);
	}
	println!("T {}", state[0]);
	let mut min = 0f64;
	let values = map(0..(0.0002/model.time_step) as usize, |_| if let [temperature, volume, active_amounts@..] = state {
		let value = ((time-plot_min_time)*1e3, vec![vec![*temperature, *volume].into_boxed_slice(), map(active_amounts, |v| v/active_amounts.iter().sum::<f64>())].into_boxed_slice());
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
	let key = map((0..active).map(|k| values.iter().map(move |(_, sets)| sets[1][k])), |Y| Y.reduce(f64::max).unwrap());
	let mut s = map(0..active, |i| i);
	let (_, _, select) = s.select_nth_unstable_by(species.len()-6, |&a,&b| key[a].partial_cmp(&key[b]).unwrap());
	let species_names = map(&*select, |&k| species_names[k]);
	let values = map(Vec::from(values).into_iter(), |(t, sets)| {
		let [sets_0, sets_1] = *std::convert::TryInto::<Box<[_;2]>>::try_into(sets).unwrap();
		(t, vec![sets_0, map(&*select, |&k| sets_1[k])].into_boxed_slice())
	});
	let mut plot = ui::plot::Plot::new(vec![&["T","V"] as &[_], &species_names].into_boxed_slice(), values.to_vec());
	if true { Ok(ui::app::run(plot)?) }
	else {
		let mut target = image::Image::zero((3840, 2160).into());
		ui::widget::Widget::paint(&mut plot, &mut target.as_mut())?;
		pub fn as_bytes<T>(slice: &[T]) -> &[u8] { unsafe{std::slice::from_raw_parts(slice.as_ptr() as *const u8, slice.len() * std::mem::size_of::<T>())} }
		Ok(png::save_buffer("/var/tmp/image.png", as_bytes(&target.data), target.size.x, target.size.y, png::ColorType::Rgba8)?)
	}
}
