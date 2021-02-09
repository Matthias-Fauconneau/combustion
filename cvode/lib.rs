#![allow(incomplete_features, non_snake_case)] #![feature(const_generics, const_evaluatable_checked, type_ascription, array_methods)]
use {sundials_sys::*, std::ffi::c_void as void};

use std::convert::TryInto;
fn n_vector<const N: usize>(v: &[f64; N]) -> N_Vector { unsafe{N_VMake_Serial(v.len() as i64, v.as_ptr() as *mut _)} }
fn r#ref<'t, const N: usize>(v: N_Vector) -> &'t [f64; N] {
	unsafe{std::slice::from_raw_parts(N_VGetArrayPointer_Serial(v), N_VGetLength_Serial(v) as usize)}.try_into().unwrap()
}
fn r#mut<'t, const N: usize>(v: N_Vector) -> &'t mut [f64; N] {
	unsafe{std::slice::from_raw_parts_mut(N_VGetArrayPointer_Serial(v), N_VGetLength_Serial(v) as usize)}.try_into().unwrap()
}

pub struct CVODE<F: Fn(&[f64; N])->Option<[f64; N]>, const N: usize>(*mut void, std::marker::PhantomData<F>);

impl<F: Fn(&[f64; N])->Option<[f64; N]>, const N: usize> CVODE<F, N> {
	pub fn new(u: &[f64; N]) -> Self {
		let cvode = unsafe{CVodeCreate(CV_BDF)};
		fn to_str<'t>(s: *const i8) -> &'t str { unsafe{std::ffi::CStr::from_ptr(s)}.to_str().unwrap() }
		extern "C" fn err(_error_code: i32, _module: *const i8, function: *const i8, msg: *mut i8, _user_data: *mut void) { panic!("{}: {}", to_str(function), to_str(msg)); }
		unsafe{CVodeSetErrHandlerFn(cvode, Some(err), std::ptr::null_mut())};
		extern "C" fn shim<F: Fn(&[f64; N])->Option<[f64; N]>, const N: usize>(_t: f64, u: N_Vector, /*mut*/ fu: N_Vector, f: *mut void) -> i32 {
			let u = r#ref(u);
			for u in u { assert!(u.is_finite()); }
			if let Some(f_u) = unsafe{(f as *const F).as_ref()}.unwrap()(u) {
				for f_u in &f_u { assert!(f_u.is_finite()); }
				r#mut::<N>(fu).copy_from_slice(&f_u);
				0
			} else { 1 }
		}
		let /*mut*/ u = n_vector(u);
		assert_eq!(unsafe{CVodeInit(cvode, Some(shim::<F, N>), /*t0:*/ 0., u)}, CV_SUCCESS);
		assert_eq!(unsafe{CVodeSStolerances(cvode, /*relative_tolerance:*/ 1e-8, /*absolute_tolerance:*/ 1e-14)}, CV_SUCCESS);
		let A = unsafe{SUNDenseMatrix(N as i64, N as i64)};
		assert_eq!(unsafe{CVodeSetLinearSolver(cvode, SUNDenseLinearSolver(u, A), A)}, CV_SUCCESS);
		CVODE(cvode, std::marker::PhantomData)
	}
	pub fn step(&mut self, f: F, target_t: f64, u: &[f64; N]) -> (f64, [f64; N]) {
		assert_eq!(unsafe{CVodeSetUserData(self.0, &f as *const F as *mut void)}, CV_SUCCESS);
		unsafe{CVodeSetStopTime(self.0, target_t)};
		let /*mut*/ u = n_vector(u);
		let mut reached_t = f64::NAN; unsafe{CVode(self.0, target_t, u, &mut reached_t, CV_ONE_STEP)};
		(reached_t, *r#ref(u))
	}
}
