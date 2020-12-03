use {std::convert::TryInto, num::{norm, error},
				iter::{array_from_iter as from_iter, into::{Collect, Enumerate, IntoChain, IntoMap}, zip, map, eval,
				vec::{eval, Dot, Suffix, Sub}}};

// Estimate principal eigenvector/value of dyF|y
fn power_iteration<const N: usize>(f: impl Fn(&[f64; N])->[f64; N], tmax: f64, y: &[f64; N], fy: &[f64; N], v: &[f64; N]) -> ([f64; N], f64) {
	let [norm_y, norm_v] = [y,v].map(norm);
	assert!(norm_y > 0.);
	let ε = norm_y * f64::EPSILON.sqrt();
	assert!(norm_v > 0.);
	let ref mut yεv = eval!(y, v; |y, v| y + ε * v / norm_v);
	let mut ρ = 0.;
	for i in 1..=50 {
		let ref fεv = f(yεv).sub(fy);
		let norm_fεv = norm(fεv);
		assert!(norm_fεv > 0.);
		let previous_ρ = ρ;
		ρ = norm_fεv / ε;
		if i >= 2 && f64::abs(ρ - previous_ρ) <= 0.01*ρ.max(1./tmax) { break; }
		*yεv = eval!(y, fεv; |y, fεv| y + (ε / norm_fεv) * fεv);
	}
	(yεv.sub(y), ρ * 1.2)
}

pub fn RKC<const N: usize>(f: impl Fn(&[f64; N])->[f64; N], relative_tolerance: f64, absolute_tolerance: f64, tmax: f64, mut u: [f64; N]) -> [f64; N] {
	let mut fu = f(&u);
	let (mut v, mut jacobian_spectral_radius) = power_iteration(f, tmax, &u, &fu, &fu);
	let max_stages = ((relative_tolerance / (10. * f64::EPSILON)).sqrt().round() as usize).max(2);
	let mut h = {
		let h = (1./jacobian_spectral_radius).min(tmax);
		let ref fu1 = self.h(P, &eval!(&u, &fu; |u, fu| u + h * fu));
		(h/(h*error(zip!(&fu, fu1, &u).map(|(fu, fu1, u):(&f64,&f64,&f64)| (fu1 - fu) / (absolute_tolerance + relative_tolerance * u.abs())))) / 10.).min(tmax)
	};
	let (mut previous_error, mut previous_h) = (0., 0.);
	let mut nstep = 0;
	let mut t = 0.;
	loop {
		let stages = 1 + (1. + 1.54 * h * jacobian_spectral_radius).sqrt().floor() as usize;
		let stages = if stages > max_stages {
			h = (max_stages*max_stages - 1) as f64 / (1.54 * jacobian_spectral_radius);
			max_stages
		} else { stages };
		let ref u1 = {
			let w0 = 1. + 2. / (13.0 * (stages * stages) as f64);
			let sqw01 = w0*w0 - 1.;
			let arg = stages as f64 * (w0 + sqw01.sqrt()).ln();
			let w1 = arg.sinh() * sqw01 / (arg.cosh() * stages as f64 * sqw01.sqrt() - w0 * arg.sinh());
			let mut B = [1. / (4.*w0*w0); 2];
			let mu_t = w1 * B[0];
			let [ref mut u0, mut u1] = [u, eval!(&u, &fu; |u, fu| u + mu_t * h * fu)];
			let mut Z = [w0, 1.];
			let mut dZ = [1., 0.];
			let mut ddZ = [0., 0.];
			for _ in 1..stages {
				let z = 2. * w0 * Z[0] - Z[1];
				let dz = 2. * w0 * dZ[0] - dZ[1] + 2. * Z[0];
				let ddz = 2. * w0 * ddZ[0] - ddZ[1] + 4. * dZ[0];
				let b = ddz / (dz * dz);
				let gamma_t = 1. - (Z[0] * B[0]);
				let nu = - b / B[1];
				let mu = 2. * b * w0 / B[0];
				let mu_t = mu * w1 / w0;
				let ref fu1 = self.h(P, &u1);
				for (u0, u1, fu1, u, fu) in zip!(u0, &mut u1, fu1, &u, &fu) {
					let u0_ = *u0;
					*u0 = *u1;
					*u1 = (1.-mu-nu)*u + nu*u0_ + mu**u1 + h*mu_t*(fu1-(gamma_t*fu));
				}
				B = [b, B[0]];
				Z = [z, Z[0]];
				dZ = [dz, dZ[0]];
				ddZ = [ddz, ddZ[0]];
			}
			u1
		};
		let ref fu1 = self.h(P, &u1);
		let ũ = map!(u1, &u, &fu, fu1; |u1, u, fu, fu1| u1-h*(fu+fu1)/2.-u);
		let error = 0.8 * error(zip!(ũ, &u, u1).map(|(ũ, u, u1):(f64,&f64,&f64)| ũ / (absolute_tolerance + relative_tolerance*u.abs().max(u1.abs()))));
		if error > 1. {
			h *= 0.8 / error.powf(1./3.);
			(v, jacobian_spectral_radius) = self.power_iteration(P, tmax, &u, &fu, &v);
			//println!("{:e} {:e} {:e} {:e} {} {} {}", jacobian_spectral_radius, t, h, error, nstep, stages, previous_error/error);
		} else {
			t += h;
			if t >= tmax { break *u1; }
			u = *u1;
			fu = *fu1;
			nstep += 1;
			if nstep%25 == 0 { (v, jacobian_spectral_radius) = self.power_iteration(P, tmax, &u, &fu, &v); }
			let factor = (0.8 * if previous_error > f64::EPSILON { h/previous_h*(previous_error/error).powf(1./3.) } else { 1./error.powf(1./3.) } ).clamp(0.1, 10.);
			if nstep%10000==0 { println!("{:e} {:e} {:e} {:e} {} {} {} {}", jacobian_spectral_radius, t, h, error, nstep, stages, factor, previous_error/error); }
			previous_error = error;
			previous_h = h;
			h = (factor*h).min(tmax);
		}
	}
}
