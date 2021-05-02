#![allow(incomplete_features, non_snake_case)]
#![feature(const_generics, const_evaluatable_checked, destructuring_assignment, array_map, unboxed_closures, fn_traits, type_ascription, trait_alias)]
fn get_mut<T: Default, U>(cell: &std::cell::Cell<T>, f: impl FnOnce(&mut T) -> U) -> U {
	let mut value = cell.take();
	let result = f(&mut value);
	cell.set(value);
	result
}

pub trait Rate<const CONSTANT: Property> = Fn(Constant<CONSTANT>, &StateVector<CONSTANT>, &mut Derivative<CONSTANT>, &mut [f64]);
pub fn rate<'t, const CONSTANT: Property>(species: &Species, reactions: impl IntoIterator<Item/*:AsRef<Reaction>*/=&'t Reaction>) -> impl Rate<CONSTANT> {
	let mut module = cranelift_jit::JITModule::new({
		//let mut flag_builder = cranelift_codegen::settings::builder();
		//use cranelift_codegen::settings::Configurable;
		//flag_builder.enable("is_pic").unwrap();
		//flag_builder.set("enable_probestack", "false").unwrap();
		let isa_builder = cranelift_native::builder().unwrap();
		let isa = isa_builder.finish(cranelift_codegen::settings::Flags::new(/*flag_builder*/cranelift_codegen::settings::builder()));
		cranelift_jit::JITBuilder::with_isa(isa, cranelift_module::default_libcall_names())
	});
	use cranelift_module::Module;
	let mut context = module.make_context();
	context.func = combustion::reaction::rate::<_, CONSTANT>(species, reactions, 1);
  let id = module.declare_function(&"", cranelift_module::Linkage::Export, &context.func.signature).unwrap();
  module.define_function(id, &mut context, &mut cranelift_codegen::binemit::NullTrapSink{}, &mut cranelift_codegen::binemit::NullStackMapSink{}).unwrap();
	module.finalize_definitions();
	let function = module.get_finalized_function(id);
	let function = unsafe{std::mem::transmute::<_,extern fn(f64, *const f64, *mut f64, *mut f64)>(function)};
	move |constant:Constant<CONSTANT>, state:&StateVector<CONSTANT>, derivative:&mut Derivative<CONSTANT>, debug: &mut [f64]| {
		let constant = constant.0;
		function(constant, state.0.as_ptr(), derivative.0.as_mut_ptr(), debug.as_mut_ptr());
	}
}

use {fehler::throws, error::Error, combustion::{*, reaction::{*, Property::*}}};

struct Row<W, const N: usize>([W; N]);
use {xy::xy, ui::widget::{size, Target, Widget, Event, EventContext}};
impl<W:Widget, const N: usize> Widget for Row<W,N> {
	fn size(&mut self, size: size) -> size { self.0[0].size(size) }
	#[throws] fn paint(&mut self, target: &mut Target) {
			for (column, widget) in self.0.iter_mut().enumerate() {
				widget.paint(&mut target.slice_mut(xy{x: (column as u32)*target.size.x/(N as u32), y: 0}, target.size/xy{x: N as u32, y: 1}))?
			}
	}
	#[throws] fn event(&mut self, size: size, event_context: &EventContext, event: &Event) -> bool {
		for (_column, widget) in self.0.iter_mut().enumerate() { if widget.event(size, event_context, event)? { return true; } } // todo: translate position
		false
	}
}

fn implicit(u: &[f64]) -> Box<[f64]> { [u[0]].iter().chain(&u[2..]).copied().collect() } // one of P or V imply the other using ideal gas law
fn explicit(total_amount: f64, pressure_R: f64, u: &[f64]) -> Box<[f64]> { // Reconstructs P|V using ideal gas law
	let temperature = u[0];
	let volume = total_amount / pressure_R * temperature;
	[temperature, volume].iter().chain(&u[1..]).copied().collect()
}

#[throws] fn main() {
	let model = &std::fs::read("CH4+O2.ron")?;
	let model = model::Model::new(&model)?;
	let (ref species_names, ref species) = Species::new(&model.species);
	let len = species.len();
	let reactions = iter::map(&*model.reactions, |r| Reaction::new(species_names, r));
	let rate = rate(species, &*reactions);
	let ref state = combustion::initial_state(&model);
	let total_amount = state.amounts.iter().sum();

	struct Derivative<Rate: self::Rate<CONSTANT>, const CONSTANT: Property> {
		rate: Rate,
		constant: Constant<CONSTANT>,
		total_amount: f64,
		derivative: std::cell::Cell<StateVector::<CONSTANT>>
	}
	impl<Rate: self::Rate<CONSTANT>, const CONSTANT: Property> FnOnce<(&[f64],)> for Derivative<Rate, CONSTANT> {
		type Output = Option<Box<[f64]>>;
		extern "rust-call" fn call_once(mut self, args: (&[f64],)) -> Self::Output { self.call_mut(args) }
	}
	impl<Rate: self::Rate<CONSTANT>, const CONSTANT: Property> FnMut<(&[f64],)> for Derivative<Rate, CONSTANT> {
		extern "rust-call" fn call_mut(&mut self, args: (&[f64],)) -> Self::Output { self.call(args) }
	}
	impl<Rate: self::Rate<CONSTANT>, const CONSTANT: Property> Fn<(&[f64],)> for Derivative<Rate, CONSTANT> {
		extern "rust-call" fn call(&self, (u,): (&[f64],)) -> Self::Output {
			let Self{rate, constant, total_amount, derivative} = self;
			get_mut(derivative, |mut derivative| {
				let pressure_R = constant.0 as f64;
				let mut debug = vec![f64::NAN; /*model.reactions.len()*/325*2].into_boxed_slice();
				rate(*constant, &StateVector(explicit(*total_amount, pressure_R, u)), &mut derivative, &mut debug);
				Some(implicit(&derivative.0))
			})
		}
	}
	let constant = state.constant();
	let derivative = /*Derivative*/StateVector(std::iter::repeat(0.).take(2+len-1).collect());
	let derivative = std::cell::Cell::new(derivative);
	let ref derivative = Derivative{
		rate,
		constant,
		total_amount,
		derivative
	};

	fn from(State{/*temperature, pressure, volume,*/ amounts, ..}: &State) -> Box<[Box<[f64]>]> {
		//eprintln!("{} {} {}", temperature, *pressure*NA, *volume);
		Box::new( [/*Box::new([*temperature, *pressure*NA, *volume]) as Box<[_]>,*/ amounts.clone()] ) as Box<[_]>
	}
	let mut time = 0.;
	//let plot = || ui::plot::Plot::new(Box::new( [&["T","P","V"] as &[_], species_names] ), vec![(time, from(&state))]);
	let plot = || ui::plot::Plot::new(Box::new( [species_names] ), vec![(time, from(&state))]);
	let app = ui::app::App::new(Row([plot()/*, plot()*/]))?;
	let state : StateVector<{Pressure}> = state.into();
	let mut state = implicit(&state);
	let mut cvode = cvode::CVODE::new(&state);
	app.run(move |plots| { //: &mut Row<[ui::plot::Plot; 2]>
		let next_time = time + model.time_step;
		let mut steps = 0;
		while time < next_time {
			(time, state) = {
					let (time, state) = cvode.step(derivative, next_time, &state);
					(time, state.to_vec().into_boxed_slice())
			};
			steps += 1;
		}
		println!("{:.0} {:.0} {}", time*1e3, state[0], steps);
		assert_eq!(time, next_time);
		plots.0[0].values.push((time*1e3/*ms*/, from(&State::new(total_amount, constant, &StateVector::<{Pressure}>(explicit(total_amount, constant.0 as f64, &state))))));
		//plots.0[1].values.push((time*1e9/*ns*/, from(&reaction(simulation, time)?)));
		Ok(true)
		//Err(error::anyhow!(""))
	})?
}
