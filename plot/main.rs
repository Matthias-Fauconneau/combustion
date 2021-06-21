#![feature(destructuring_assignment)]
#![allow(non_snake_case)]

fn main() -> anyhow::Result<()> {
	use {num::sq, iter::{map, dot}, combustion::{*, reaction::*}};
	let model = &std::fs::read("../LiDryer.ron").unwrap();
	let model = model::Model::new(&model).unwrap();
	let (ref species_names, ref species) = Species::new(&model.species);
	let exp_Gibbs_RT = exp_Gibbs_RT(&species.thermodynamics[0..species.len()-1]);
	let rates = rates(&map(&*model.reactions, |r| Reaction::new(species_names, r)));

	let State{pressure_R, temperature, amounts, ..} = combustion::initial_state(&model);
	let total_amount = amounts.iter().sum::<f64>();
	let f = |u: &[f64], f_u: &mut [f64]| if let &[T, ref active_amounts@..] = u {
		let total_concentration = pressure_R / T; // Constant pressure
		let density = total_concentration / total_amount;
		let active_concentrations = map(active_amounts, |&n| density*f64::max(0., n));
		let inert_concentration = total_concentration - active_concentrations.iter().sum::<f64>();
		let concentrations = [&active_concentrations as &[_],&[inert_concentration]].concat();
		let log_T = f64::log2(T);
		let T2 = T*T;
		let T3 = T*T2;
		let T4 = T*T3;
		let rcp_T = 1./T;
		let exp_Gibbs0_RT = {let mut output = vec![0.; active_concentrations.len()]; exp_Gibbs_RT(&[log_T,T,T2,T3,T4,rcp_T],&[], &mut output); output};
		let P0_RT = NASA7::reference_pressure / T;
		let dtω = &mut f_u[1..];
		rates(&[log_T,T,T2,T4,rcp_T,sq(rcp_T),P0_RT,1./P0_RT], &[&exp_Gibbs0_RT, &concentrations], dtω);
		let dtT_T = dot(dtω.iter().copied().zip(species.thermodynamics.iter().map(|specie| specie.enthalpy_RT(T))))
										/ dot(concentrations.into_iter().zip(species.thermodynamics.iter().map(|specie| specie.molar_heat_capacity_at_constant_pressure_R(T))));
		f_u[0] = T * dtT_T;
		use itertools::Itertools;
		println!("{:3} {:.2e}", "u", u.iter().format(", "));
		println!("{:3} {:.2e}", "f_u", f_u.iter().format(", "));
		true
	} else { unreachable!() };

	let state = [&[temperature], &amounts[..amounts.len()-1]].concat();

	let mut cvode = cvode::CVODE::new(/*relative_tolerance:*/ 1e-4, /*absolute_tolerance:*/ 1e-7, &state);
	let plot_min_time = 0.4;
	let (mut time, mut state) = (0., &*state);
	while time < plot_min_time {
		let next_time = time + model.time_step;
		while time < next_time { (time, state) = cvode.step(&f, next_time, state); }
	}
	let mut min = 0f64;
	let values = map(0..1000, |_| if let [temperature, active_amounts@..] = state {
		let value = ((time-plot_min_time)*1e3, vec![vec![*temperature].into_boxed_slice(), map(active_amounts, |v| v/active_amounts.iter().sum::<f64>())].into_boxed_slice());
		let next_time = time + model.time_step;
		while time < next_time { (time, state) = cvode.step(&f, next_time, state); }
		min = min.min(state[1..].iter().copied().reduce(f64::min).unwrap());
		value
	} else { unreachable!() });
	println!("{:.0e}", min);
	let key = map((0..species.len()).map(|k| values.iter().map(move |(_, sets)| sets[1][k])), |Y| Y.reduce(f64::max).unwrap());
	let mut s = map(0..species.len(), |i| i);
	let (_, _, select) = s.select_nth_unstable_by(species.len()-6, |&a,&b| key[a].partial_cmp(&key[b]).unwrap());
	let species_names = map(&*select, |&k| species_names[k]);
	let values = map(Vec::from(values).into_iter(), |(t, sets)| {
		let [sets_0, sets_1] = *std::convert::TryInto::<Box<[_;2]>>::try_into(sets).unwrap();
		(t, vec![sets_0, map(&*select, |&k| sets_1[k])].into_boxed_slice())
	});
	let mut plot = ui::plot::Plot::new(vec![&["T"] as &[_], &species_names].into_boxed_slice(), values.to_vec());
	if true { ui::app::run(plot) }
	else {
		let mut target = image::Image::zero((3840, 2160).into());
		ui::widget::Widget::paint(&mut plot, &mut target.as_mut())?;
		pub fn as_bytes<T>(slice: &[T]) -> &[u8] { unsafe{std::slice::from_raw_parts(slice.as_ptr() as *const u8, slice.len() * std::mem::size_of::<T>())} }
		Ok(png::save_buffer("/var/tmp/image.png", as_bytes(&target.data), target.size.x, target.size.y, png::ColorType::Rgba8)?)
	}
}
