#![feature(trait_alias,destructuring_assignment,default_free_fn,let_else,iter_zip,unboxed_closures,fn_traits,box_patterns)]
#![allow(non_snake_case,non_upper_case_globals)]
mod yaml; mod device;
use {std::iter::zip, anyhow::Result, iter::map, itertools::Itertools, combustion::*, device::*};
fn main() -> Result<()> {
	let path = std::env::args().skip(1).next().unwrap_or("LiDryer".to_string());
	let path = if std::path::Path::new(&path).exists() { path } else { format!("/usr/local/share/cantera/data/{path}.yaml") };
	let model = yaml::Loader::load_from_str(std::str::from_utf8(&std::fs::read(&path).expect(&path))?)?;
	let model = yaml::parse(&model);
	let (species_names, ref species, active, reactions, ref _state) = new(&model);
	let rates = reaction::rates(&species.molar_mass, &species.thermodynamics, &reactions, &species_names);
	let block_size = 1;
	let rates = assemble(rates, block_size);

	let ref state@State{pressure_R, volume, temperature, ref amounts} = {
		let State{pressure_R, ..} = initial_state(&model);
		let equivalence_ratio = 1.;
		//let amounts = std::collections::HashMap::<_, _>::from_iter([("H2", equivalence_ratio),("O2", 1./2.),("N2",79./21./2.)]);
		let amounts = std::collections::HashMap::<_, _>::from_iter([("CH4", equivalence_ratio),("O2", 2.),("N2",79./21.)]);
		let amounts = map(&*species_names, |s| *amounts.get(s).unwrap_or(&0.));
		State{pressure_R, volume: 1., temperature: 1100., amounts}
	};
	//let amounts = &state.amounts;
	#[derive(serde::Serialize)] struct StandardState { pressure: f64, temperature: f64, X: Box<[f64]> }
	impl From<&State> for StandardState { 
		fn from(State{pressure_R, temperature, amounts, ..}: &State) -> Self {
			Self{pressure: pressure_R*R, temperature: *temperature, X: amounts.clone()}
		}
	}
	if true { println!("{}", serde_yaml::to_string(&StandardState::from(state))?); }
	assert!(active < amounts.len()); // Assume bulk specie is inert
	let state : Box<[f32]> = map([&[temperature, volume], &amounts[..active]].concat(), |v| v as _);
	let input = [&*state, &*map(&amounts[active..amounts.len()-1], |&v| v as _)].concat();
	let states_len = 1;
	let input = map(input.into_iter(), |x| vec![x as _; states_len]);

	//let mut evaluations = 0;
	#[cfg(not(feature="vpu"))] let mut f = {
		let mut input = input;
		move |u: &[f32], f_u: &mut [f32]| {
			assert!(u[0]>200.);
			assert!(u.iter().all(|u| u.is_finite()));
			for (input, &u) in zip(input[..2+active].iter_mut(), u) { input[0] = u; } // T, V, active, non bulk inert (non bulk inerts are input parameters but not reaction state)
			f_u.copy_from_slice(&map(&*rates(&[pressure_R as _, (1./pressure_R) as _], &map(&*input, |x| &**x)).unwrap(), |y| all_same(y, states_len)));
			assert!(f_u.iter().all(|u| u.is_finite()), "{u:?} {f_u:?}");
			//evaluations += 1;
			true
		}
	};
	#[cfg(feature="vpu")] use {iter::list, vulkan::*};
	#[cfg(feature="vpu")] let Function{ref device, output_len, block_size, pipeline,..} = rates;
	#[cfg(feature="vpu")] 	let mut input = map(&*input, |array| Buffer::new(device, array).unwrap());
	#[cfg(feature="vpu")] 	let output = map(0..output_len, |_| Buffer::new(device, &vec![0.; states_len]).unwrap());
	#[cfg(feature="vpu")] let mut f = {
		let buffers = list(input.iter().chain(&*output));
		device.bind(pipeline.descriptor_set, &buffers)?;
		let mut input = map(input.iter_mut(), |array| array.map_mut(device).unwrap());
		let output = map(&*output, |array| array.map(device).unwrap());
		let constants : [f32; 2] = [pressure_R as _, (1./pressure_R) as _];
		let command_buffer = device.command_buffer(&pipeline, as_u8(&constants), (states_len as u32)/(block_size as u32))?;
		move |u: &[f32], f_u: &mut [f32]| {
			//println!("{:3} {:.2e}", "u", u.iter().format(", "));
			assert!(u[0]>200. && u[0] < 3000., "{u:?}");
			assert!(u.iter().all(|u| u.is_finite()));
			//input[..2+active].copy_from_slice(u); // T, V, active, non bulk inert (non bulk inerts are input parameters but not reaction state)
			//f_u.copy_from_slice(&rates(&[pressure_R as _, (1./pressure_R) as _], &input).unwrap());
			for (input, &u) in zip(input[..2+active].iter_mut(), u) { input[0] = u; } // T, V, active, non bulk inert (non bulk inerts are input parameters but not reaction state)
			//println!("{u:?} {f_u:?}");
			device.submit_and_wait(command_buffer).unwrap(); // constant constants
			for (output, f_u) in zip(&*output, f_u.iter_mut()) { *f_u = output[0]; }
			assert!(f_u.iter().all(|u| u.is_finite()), "{u:?} {f_u:?}");
			//assert!(num::abs(f_u[0]) < 1., "{u:?} {f_u:?}");
			//evaluations += 1;
			//println!("{:3} {:.2e}", "f_u", f_u.iter().format(", "));
			true
		}
	};

	//let mut integrator = cvode::CVODE::new(/*relative_tolerance:*/ 1e-4, /*absolute_tolerance:*/ 1e-7, &state)
	let mut integrator = {
		struct Explicit<F: FnMut(&[f32], &mut [f32])->bool> { f_u: Box<[f32]>, _marker: std::marker::PhantomData<F> }
		impl<F: FnMut(&[f32], &mut [f32])->bool> Explicit<F> {
			pub fn step(&mut self, f: &mut F, dt: f32, u: &mut [f32]) -> f32 {
				f(&map(&*u,|&u| u as _), &mut self.f_u);
				for (u, &f_u) in zip(u.iter_mut(), &*self.f_u) { *u += dt*f_u as f32; }
				dt
			}
		}
		Explicit{f_u: vec![0.; state.len()].into(), _marker: std::marker::PhantomData}
	};
	let plot_min_time = 0.;
	let (mut time, mut state) = (0., state);
	//let start = std::time::Instant::now();
	println!("{plot_min_time}");
	while time < plot_min_time {
		time += integrator.step(&mut f, model.time_step as f32, &mut state) as f64;
		//dbg!(time/plot_min_time, time, state[0]);
		//assert!(state[0]<1500.);
	}
	//eprintln!("T {}", state[0]);
	eprintln!("T {}", state[0]);
	//let mut min:f32 = 0.;
	let duration = 0.3;
	//let steps = duration/model.time_step;
	let points = 3840;
	let point_duration = duration/(points as f64);
	//let points_steps = ((points as f64)/steps).ceil();
	//assert!(points_steps == 1);
	let points = map(0..points, |_| {
		let [temperature, volume, active_amounts@..] = &state[..] else { unreachable!() };
		let point = ((time-plot_min_time)*1e3, vec![vec![*temperature as f64, (*volume*1e3) as f64].into_boxed_slice(), map(active_amounts, |&v| (v as f32/active_amounts.iter().sum::<f32>()) as f64)].into_boxed_slice());
		let point_end_time = time+point_duration;
		while time < point_end_time {
			time += integrator.step(&mut f, model.time_step as f32, &mut state) as f64;
		}
		//min = min.min(state[1..].iter().copied().reduce(f32::min).unwrap());
		point
	});
	//let end = std::time::Instant::now();
	//let time = (end-start).as_secs_f64();
	//println!("T {}", state[0]);
	//println!("{evaluations} / {time}s = {}", (evaluations as f64)/time);
	//println!("{:.0e}", min);
	if true {
		let OH = species_names.iter().position(|&name| name=="OH").unwrap();
		let max_OH = points.iter().map(|point|{
			let (_, box [ box [_, _], ref active_amounts]) = point else { unreachable!() };
			ordered_float::OrderedFloat(active_amounts[OH])
		}).position_max().unwrap();
		println!("{max_OH} {}", (max_OH as f64)*model.time_step);
		let from = /*amounts*/|point: &(_, Box<[Box<[f64]>]>)| -> State {
			let (_, box [ box [temperature, volume], ref active_amounts]) = point else { unreachable!() };
			let amounts = [&active_amounts, &amounts[active..]].concat().into_boxed_slice(); // Inerts stays during reaction
			State{pressure_R, volume: *volume, temperature: *temperature, amounts}
		};
		println!("{}", serde_yaml::to_string(&StandardState::from(&from(&points[max_OH])))?);
		println!("{}", serde_yaml::to_string(&StandardState::from(&from(points.last().unwrap())))?);
	}
	let key = map((0..active).map(|k| points.iter().map(move |(_, sets)| sets[1][k])), |Y| Y.reduce(f64::max).unwrap());
	let mut s = map(0..active, |i| i);
	let (_, _, select) = s.select_nth_unstable_by(species.len()-5, |&a,&b| key[a].partial_cmp(&key[b]).unwrap());
	let species_names = map(&*select, |&k| species_names[k]);
	let points = map(Vec::from(points).into_iter(), |(t, sets)| {
		let [sets_0, sets_1] = *std::convert::TryInto::<Box<[_;2]>>::try_into(sets).unwrap();
		(t as f64, vec![sets_0, map(&*select, |&k| sets_1[k])].into_boxed_slice())
	});
	let keys = vec![&["T","V"] as &[_], &species_names];
	let mut plot = ui::plot::Plot::new(&keys, &points);
	if false { Ok(ui::app::run(plot)?) }
	else {
		let mut target = image::Image::zero((3840, 2160).into());
		ui::widget::Widget::paint(&mut plot, &mut target.as_mut())?;
		pub fn as_bytes<T>(slice: &[T]) -> &[u8] { unsafe{std::slice::from_raw_parts(slice.as_ptr() as *const u8, slice.len() * std::mem::size_of::<T>())} }
		Ok(png::save_buffer("/var/tmp/image.png", as_bytes(&target.data), target.size.x, target.size.y, png::ColorType::Rgba8)?)
	}
}
