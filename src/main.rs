#![allow(mixed_script_confusables, non_snake_case, incomplete_features)]#![feature(type_ascription, array_map, non_ascii_idents, const_generics, const_evaluatable_checked)]
use combustion::*;
#[fehler::throws(Box<dyn std::error::Error>)] fn main() {
	let system = std::fs::read("CH4+O2.ron")?;
	let simulation = Simulation::<35>::new(&system)?;
	dbg!(cvode(&simulation));
	//test_transport_cantera(&simulation);
}

use {sundials_sys::*, std::ffi::c_void as void};

use std::convert::TryInto;
fn n_vector(v: &[f64]) -> N_Vector { unsafe{N_VMake_Serial(v.len() as i64, v.as_ptr() as *mut _)} }
fn r#ref<'t, const N: usize>(v: N_Vector) -> &'t [f64; N] { unsafe{std::slice::from_raw_parts(N_VGetArrayPointer_Serial(v), N_VGetLength_Serial(v) as usize)}.try_into().unwrap()	}
fn r#mut<'t, const N: usize>(v: N_Vector) -> &'t mut [f64; N] { unsafe{std::slice::from_raw_parts_mut(N_VGetArrayPointer_Serial(v), N_VGetLength_Serial(v) as usize)}.try_into().unwrap()	}

fn cvode<const S: usize>(Simulation{system, time_step, pressure_R, state: combustion::State{temperature, amounts}, ..}: &Simulation<S>) -> (f64, [f64; 1+S-1]) where [(); S-1]:, [(); 1+S-1]: {
	let cvode = unsafe{CVodeCreate(CV_BDF)};
	#[derive(Clone, Copy)] struct UserData<'t, const S: usize> where [(); S-1]: {system: &'t System<S>, pressure_R: f64}
	extern "C" fn rhs<const S: usize>(_t: f64, y: N_Vector, /*mut*/ ydot: N_Vector, user_data: *mut void) -> i32 where [(); S-1]:, [(); 1+S-1]: {
		let UserData{system, pressure_R} = unsafe{*(user_data as *const UserData<S>)};
		let (derivative, /*jacobian*/) = system.derivative/*and_jacobian*/(pressure_R, r#ref(y));
		r#mut::<{1+S-1}>(ydot).copy_from_slice(&derivative);
		0
	}
	let y0: [_; 1+S-1] = {use iter::{array_from_iter as from_iter, into::IntoChain}; from_iter([*temperature].chain(amounts[..S-1].try_into().unwrap():[_;S-1]))};
	fn to_str<'t>(s: *const i8) -> &'t str { unsafe{std::ffi::CStr::from_ptr(s)}.to_str().unwrap() }
	extern "C" fn err(_error_code: i32, _module: *const i8, function: *const i8, msg: *mut i8, _user_data: *mut void) { panic!("{}: {}", to_str(function), to_str(msg)); }
	unsafe{CVodeSetErrHandlerFn(cvode, Some(err), std::ptr::null_mut())};
	assert_eq!(unsafe{CVodeInit(cvode, Some(rhs::<S>), /*t0:*/ 0., n_vector(&y0))}, CV_SUCCESS);
	assert_eq!(unsafe{CVodeSStolerances(cvode, /*relative_tolerance:*/ 1e-8, /*absolute_tolerance:*/ 1e-14)}, CV_SUCCESS);
	assert_eq!(unsafe{CVodeSetUserData(cvode, ((&UserData{system, pressure_R: *pressure_R} as *const UserData::<S>) as *const void) as *mut void)}, CV_SUCCESS);
	let N = 1+S-1;
	let A = unsafe{SUNDenseMatrix(N as i64, N as i64)};
	let /*mut*/ y = n_vector(&y0);
	assert_eq!(unsafe{CVodeSetLinearSolver(cvode, SUNDenseLinearSolver(y, A), A)}, CV_SUCCESS);
	let mut t = f64::NAN; unsafe{CVode(cvode, /*t:*/ *time_step, y, &mut t, CV_ONE_STEP)};
	(t, *r#ref(y))
}

#[allow(dead_code)] fn test_transport_cantera<const S: usize>(Simulation{system, state, pressure_R, species_names, ..}: &Simulation<S>) where [(); S-1]: {
	let transport = system.transport(*pressure_R, &state);
	let cantera = {
		use itertools::Itertools;
		let mole_proportions = format!("{}", species_names.iter().zip(&state.amounts).filter(|(_,&n)| n > 0.).map(|(s,n)| format!("{}:{}", s, n)).format(", "));
		let mole_proportions = std::ffi::CString::new(mole_proportions).unwrap();
		use std::ptr::null;
		let ([mut viscosity, mut thermal_conductivity], mut species_len, mut specie_names, mut mixture_averaged_thermal_diffusion_coefficients) = ([0.; 2], 0, null(), null());
		unsafe {
			let pressure = pressure_R * (combustion::kB*combustion::NA);
			extern "C" { fn cantera(pressure: f64, temperature: f64, mole_proportions: *const std::os::raw::c_char,
				viscosity: &mut f64, thermal_conductivity: &mut f64,
				species_len: &mut usize, species: &mut *const *const std::os::raw::c_char, mixture_averaged_thermal_diffusion_coefficients: &mut *const f64); }
			cantera(pressure, state.temperature, mole_proportions.as_ptr(), &mut viscosity, &mut thermal_conductivity, &mut species_len, &mut specie_names, &mut mixture_averaged_thermal_diffusion_coefficients);
			let specie_names = iter::box_collect(std::slice::from_raw_parts(specie_names, species_len).iter().map(|&s| std::ffi::CStr::from_ptr(s).to_str().unwrap()));
			let order = |o:&[_]| iter::vec::eval(species_names, |s| o[specie_names.iter().position(|&k| k==s.to_uppercase()).expect(&format!("{} {:?}", s, species_names))]);
			let mixture_averaged_thermal_diffusion_coefficients = order(std::slice::from_raw_parts(mixture_averaged_thermal_diffusion_coefficients, species_len));
			Transport{viscosity, thermal_conductivity, mixture_averaged_thermal_diffusion_coefficients}
		}
	};
	dbg!(&transport, &cantera);
	dbg!((transport.viscosity-cantera.viscosity)/cantera.viscosity, (transport.thermal_conductivity-cantera.thermal_conductivity)/cantera.thermal_conductivity);
	assert!(f64::abs(transport.viscosity-cantera.viscosity)/cantera.viscosity < 0.03, "{}", transport.viscosity/cantera.viscosity);
	assert!(f64::abs(transport.thermal_conductivity-cantera.thermal_conductivity)/cantera.thermal_conductivity < 0.05, "{:?}", (transport.thermal_conductivity, cantera.thermal_conductivity));
}
