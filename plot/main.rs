#![allow(incomplete_features, non_snake_case)]
#![feature(const_generics, const_evaluatable_checked, destructuring_assignment, array_map, unboxed_closures, fn_traits, type_ascription, trait_alias, box_patterns)]
fn get_mut<T: Default, U>(cell: &std::cell::Cell<T>, f: impl FnOnce(&mut T) -> U) -> U {
	let mut value = cell.take();
	let result = f(&mut value);
	cell.set(value);
	result
}

/*struct Row<W, const N: usize>([W; N]);
use {xy::xy, ui::widget::{size, Target, Widget, Event, EventContext}};
impl<W:Widget, const N: usize> Widget for Row<W,N> {
	fn size(&mut self, size: size) -> size { self.0[0].size(size) }
	#[fehler::throws(anyhow::Error)] fn paint(&mut self, target: &mut Target) {
			for (column, widget) in self.0.iter_mut().enumerate() {
				widget.paint(&mut target.slice_mut(xy{x: (column as u32)*target.size.x/(N as u32), y: 0}, target.size/xy{x: N as u32, y: 1}))?
			}
	}
	#[fehler::throws(anyhow::Error)] fn event(&mut self, size: size, event_context: &EventContext, event: &Event) -> bool {
		for (_column, widget) in self.0.iter_mut().enumerate() { if widget.event(size, event_context, event)? { return true; } } // todo: translate position
		false
	}
}*/

#[fehler::throws(anyhow::Error)] fn main() {
	let model = &std::fs::read("../LiDryer.ron").unwrap();
	use combustion::{*, reaction::{*, Property::*}};
	let model = model::Model::new(&model).unwrap();
	let (ref species_names, ref species) = Species::new(&model.species);
	let len = species.len();
	let reactions = iter::map(&*model.reactions, |r| Reaction::new(species_names, r));

	use reaction::*;
	let exp_Gibbs_RT = exp_Gibbs_RT(&species.thermodynamics[0..species.len()-1]);
	let rates = rates(&iter::map(&*model.reactions, |r| Reaction::new(species_names, r)));

	let f = |u, f_u| {
		let [temperature, active_amounts..] = u;
		let total_concentration = pressure_R / T;
		let density = total_concentration / total_amount;
		let active_concentrations = map(active_amounts, |&n| density*max(0., n)));
		let inert_concentration = total_concentration - active_concentrations.sum();
		let ref concentrations = [&active_concentrations as &[_],&[inert_concentration]].concat();
		let log_T = f64::log2(T);
		let T2 = T*T;
		let T3 = T*T2;
		let T4 = T*T3;
		let rcp_T = 1./T;
		let exp_Gibbs0_RT = vec![0.; u.len()];
		exp_Gibbs_RT(&[log_T,T,T2,T3,T4,rcp_T],&[], &mut exp_Gibbs0_RT);
		let P0_RT = NASA7::reference_pressure / T;
		rates(&[log_T,T,T2,T4,rcp_T,num::sq(rcp_T),P0_RT,1./P0_RT], &[&exp_Gibbs0_RT, &concentrations], du)
	}

	let ref state = combustion::initial_state(&model);
	let total_amount = state.amounts.iter().sum();
	let state = [&[state.temperature], amounts[..amounts.len()-1]].concat();

	let mut cvode = cvode::CVODE::new(/*relative_tolerance:*/ 1e-4, /*absolute_tolerance:*/ 1e-7, &state);
	let plot_min_time = 0.4;
	let (mut time, mut state) = (0., &*state);
	while time < plot_min_time {
		let next_time = time + model.time_step;
		while time < next_time { (time, state) = cvode.step(f, next_time, state); }
	}
	let mut min = 0f64;
	let values = iter::eval(1000, |_| {
		fn from([temperature, active_amounts..] : &[f64]) -> Box<[Box<[f64]>]> {
        vec![vec![*temperature].into_boxed_slice(), iter::map(&**active_amounts, |v| v/active_amounts.iter().sum::<f64>())].into_boxed_slice()
		}
		let value = ((time-plot_min_time)*1e3))));
		let next_time = time + model.time_step;
		while time < next_time { (time, state) = cvode.step(f, next_time, state); }
		min = min.min(state[1..].iter().copied().reduce(f64::min).unwrap());
		value
	});
	println!("{:.0e}", min);
	let key = iter::map((0..len).map(|k| values.iter().map(move |(_, sets)| sets[1][k])), |Y| Y.reduce(f64::max).unwrap());
	let mut s = iter::eval(len, |i| i);
	let (_, _, select) = s.select_nth_unstable_by(len-6, |&a,&b| key[a].partial_cmp(&key[b]).unwrap());
	let species_names = iter::map(&*select, |&k| species_names[k]);
	let values = iter::map(Vec::from(values).into_iter(), |(t, sets)| {
		use std::convert::TryInto;
		let sets:Box<[_;2]> = sets.try_into().unwrap();
		let box [sets_0, sets_1] = sets;
		(t, vec![sets_0, iter::map(&*select, |&k| sets_1[k])].into_boxed_slice())
	});
	let mut plot = ui::plot::Plot::new(vec![&["T",/*"P",*/"V"] as &[_], &species_names].into_boxed_slice(), values.to_vec());
	if true { ui::app::run(plot)? }
	else {
		let mut target = image::Image::zero((3840, 2160).into());
		plot.paint(&mut target.as_mut())?;
		pub fn as_bytes<T>(slice: &[T]) -> &[u8] { unsafe{std::slice::from_raw_parts(slice.as_ptr() as *const u8, slice.len() * std::mem::size_of::<T>())} }
		png::save_buffer("/var/tmp/image.png", as_bytes(&target.data), target.size.x, target.size.y, png::ColorType::Rgba8).unwrap();
	}
}
