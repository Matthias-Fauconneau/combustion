pub fn bucket<I:IntoIterator<Item:Eq>>(iter: I) -> impl std::iter::IntoIterator<Item=(I::Item, Vec<usize>)> {
	let mut map = linear_map::LinearMap::<_, Vec<_>>::new();
	for (index, key) in iter.into_iter().enumerate() { map.entry(key).or_insert(Default::default()).push(index) }
	map
}

use {std::{ops::Deref,default::default,iter::zip}, iter::{Prefix, eval, list, map, DotN}, std::collections::{HashSet as Set, hash_map::HashMap as Map}, ast::*};

#[derive(derive_more::Deref,derive_more::DerefMut)] pub struct Block<'t> {
	#[deref]#[deref_mut] pub block: ast::Block<'t>,
	values: Map<Expression, Value>,
	expressions: Vec<(Expr, String)>,
	after_CSE: Set<Expr>,
}
impl<'t> Block<'t> {
	pub fn new(values: &'t mut Vec<String>) -> Self { Self{block: ast::Block::new(values), values: default(), expressions: vec![], after_CSE: Set::new()} }
	#[track_caller] pub fn def(&mut self, e: impl Into<Expression>, name: impl ToString) -> Value {
		let e = e.into();
		assert!(!e.is_leaf(), "{e:?}");
		*self.values.entry(e.into()).or_insert_with_key(|e| self.block.def(e.clone(), name))
	}
}
#[track_caller] pub fn cdef(e: impl Into<Expression>, f: &mut Block, name: impl ToString) -> Expression {
	let e = e.into(); if let Some(x) = e.f64() { x.into() } else { f.def(e, name).into() }
}

#[track_caller] fn check(e: &Expression, f: &mut Block) -> Result<(),String> {
	if !e.is_leaf() {
		let new = !f.values.contains_key(e) && f.after_CSE.insert(e.deref().clone());
		if !new { return Err(format!("{} already in {} [{}]", e.to_string(f.names), f.names[f.values[e].0], {use itertools::Itertools; f.expressions.iter().filter_map(|(k,v)| (k==e.deref()).then(||v)).format(" ")})); }
	}
	match e.visit(|e| check(e, f)) {
		[Some(Err(e)), _]|[_, Some(Err(e))] => Err(e),
		_  => Ok(()),
	}
}

#[track_caller] fn chk(e: impl Into<Expression>, f: &mut Block) -> Expression { let e = e.into(); check(&e, f).expect("chk"); e }

struct Ratio(Expression, Expression);
impl<E:Into<Expression>> From<E> for Ratio { fn from(e: E) -> Ratio { Ratio(e.into(), (1.).into()) } }
impl std::ops::Mul<Ratio> for Expr { type Output = Expression; fn mul(self, Ratio(n, d): Ratio) -> Expression {
	if let Some(d) = d.f64() { if d == 1. { return self*n } }
	(self*n)/d
}}
impl std::ops::Mul<Expr> for Ratio { type Output = Expression; fn mul(self, e: Expr) -> Expression { e*self } }
impl std::ops::Mul<Ratio> for &Value { type Output = Expression; fn mul(self, r: Ratio) -> Expression { Expr::from(self)*r } }
impl std::ops::Mul<Value> for Ratio { type Output = Expression; fn mul(self, v: Value) -> Expression { Expr::from(v)*self } }

// To reuse partial efficiencies sum
//fn sum(iter: impl IntoIterator<Item:Into<Value>>, f: &mut Block) -> Option<Value> { iter.into_iter().map(|e| e.into()).reduce(|s,t| f.def(s+t, "Σ")) }

fn product_of_exponentiations_<'t, N:Into<i16>>(iter: impl IntoIterator<Item=(&'t Value, N)>, f: &mut Block, name: &str) -> Value {
	let (num, div) : (Vec::<_>,Vec::<_>) = iter.into_iter().map(|(v,n)| (v,n.into())).filter(|&(_,n)| n!=0).partition(|&(_,n)| n>0);
	let num = num.into_iter().fold(None, |mut p:Option<Value>, (v,n)|{ for _ in 0..n { p = Some(if let Some(p) = p { f.def(p*v, name) } else { v.into() }); } p });
	let div = div.into_iter().fold(None, |mut p:Option<Value>, (v,n)|{ for _ in 0..-n { p = Some(if let Some(p) = p { f.def(p*v, name) } else { v.into() }); } p });
	match (num, div) {
		(None, None) => None,
		(Some(num), None) => Some(num),
		(None, Some(div)) => Some(f.def(1./div, name)),
		(Some(num), Some(div)) => Some(f.def(num/div, name))
	}.unwrap()
}
fn product_of_exponentiations__<N:iter::IntoExactSizeIterator>(v: &[Value], n: N, f: &mut Block, name: &str) -> Value where <N as iter::IntoIterator>::Item:Into<i16> {
	use iter::Zip;
	product_of_exponentiations_(v.zip(n), f, name)
}
fn product_of_exponentiations<I:iter::IntoExactSizeIterator,N:Copy+Into<i16>>(v: &[Value], n: I, f: &mut Block, name: &str) -> Value where <I as iter::IntoIterator>::Item:Deref<Target=N> {
	use iter::Map;
	product_of_exponentiations__(v, n.map(|x:<I as iter::IntoIterator>::Item| *x), f, name)
}

#[derive(Clone,Copy)] pub struct T { pub lnT: Value, pub T: Value, pub T2: Value, pub T3: Value, pub T4: Value, pub rcpT: Value, pub rcpT2: Value}

fn molar_heat_capacity_at_constant_pressure_R(a: &[f64; 7], T{T,T2,T3,T4,..}: T) -> Expression { (*a.prefix()).dot([Expr::from(1.),T.into(),T2.into(),T3.into(),T4.into()]) }
fn enthalpy_RT(a: &[f64; 7], T{T,T2,T3,T4,rcpT,..}: T) -> Expression { a[0] + a[1]/2.*T + a[2]/3.*T2 + a[3]/4.*T3 + a[4]/5.*T4 + a[5]*rcpT }

fn Gibbs0_RT(a: &[f64; 7]) -> [f64; 7] { [a[0]-a[6], -a[0], -a[1]/2., (1./3.-1./2.)*a[2], (1./4.-1./3.)*a[3], (1./5.-1./4.)*a[4], a[5]] }
fn thermodynamic<'a, A: 'a, R: std::iter::Sum<<&'a A as std::ops::Mul<Expr>>::Output>>(a: &'a [A; 7], T{lnT,T,T2,T3,T4,rcpT,..}: T) -> R where &'a A: std::ops::Mul<Expr> {
	a.dot([Expr::from(1.), lnT.into(), T.into(), T2.into(), T3.into(), T4.into(), rcpT.into()])
}

fn thermodynamics<const N: usize>(thermodynamics: &[NASA7], expressions: [impl Fn(&[f64; 7], T)->Expression; N], Ts@T{T,..}: T, f: &mut Block, debug: [&str; N]) -> [Box<[Expression]>; N] {
	let mut specie_results = eval(|_| map(thermodynamics, |_| None));
	for (temperature_split, ref species) in bucket(thermodynamics.iter().map(|s| ordered_float::OrderedFloat(s.temperature_split))) {
		if temperature_split.is_nan() {
			for (expression, specie_results) in expressions.iter().zip(specie_results.iter_mut()) { for &specie in species {
				let e = expression(&thermodynamics[specie].pieces[0], Ts);
				let e = if e.is_leaf() { e } else { f.def(e, "").into() };
				assert!(specie_results[specie].replace(e).is_none());
			}}
		} else {
			let (exprs, results):(Vec<_>,Vec<_>) = expressions.iter().zip(specie_results.iter_mut()).zip(debug).map(|((expression, specie_results), debug)| {
				let results = map(species, |specie| f.value(format!("{debug}[{specie}]")));
				use iter::Zip;
				for (&specie, result) in species.zip(&*results) { assert!(specie_results[specie].replace(result.into()).is_none()) }
				let /*mut*/ true_exprs = map(species, |&specie| expression(&thermodynamics[specie].pieces[0], Ts));
				let /*mut*/ false_exprs = map(species, |&specie| expression(&thermodynamics[specie].pieces[1], Ts));
				//let defs = eliminate_common_subexpressions(&mut true_exprs, &mut false_exprs, f);
				//for def in defs { let Statement::Value{value,id} = &def else { panic!() }; assert!(f.values.insert(value.clone(), *id).is_none()); f.statements.push(def); }
				//for e in true_exprs.iter().chain(&*false_exprs) { check(e, f).unwrap(); }
				([true_exprs, false_exprs], results)
			}).unzip();
			let (true_exprs, false_exprs):(Vec<_>,Vec<_>) = exprs.into_iter().map(|[a,b]| (a,b)).unzip();
			let (true_exprs, false_exprs, results) = (true_exprs.concat().into(), false_exprs.concat().into(), results.concat().into());
			push(Statement::Select{condition: /*le!(f,*/(less_or_equal(T, f64::from(temperature_split))), true_exprs, false_exprs, results}, f);
		}
	}
	specie_results.map(|specie_results| map(specie_results.into_vec().into_iter(), Option::unwrap))
}

// A.T^β.exp(-Ea/kT) => A.exp(βlnT-Ea/kT)
fn arrhenius(&RateConstant{preexponential_factor: A, temperature_exponent, activation_temperature}: &RateConstant, T{lnT,T,T2,T3,T4,rcpT,rcpT2}: T, f: &mut Block) -> Expr {
	let (temperature_factor, temperature_exponent_remainder) =
		/*if temperature_exponent <= -8. { unimplemented!("{temperature_exponent}") }
		// T~1000: T^-n << 1 so no need to factorize exponent out of exp to reduce input domain
		else*/ if temperature_exponent <= -1.5 { (Some(rcpT2), temperature_exponent+2.) }
		else if temperature_exponent <= -0.5 { (Some(rcpT), temperature_exponent+1.) }
		//else if temperature_exponent >= 8. { unimplemented!("{temperature_exponent}") }
		// TODO: T~1000: factorize exponent out of exp to reduce input domain
		else if temperature_exponent >= 3.5 { (Some(T4), temperature_exponent-4.) }
		else if temperature_exponent >= 2.5 { (Some(T3), temperature_exponent-3.) }
		else if temperature_exponent >= 1.5 { (Some(T2), temperature_exponent-2.) }
		else if temperature_exponent >= 0.5 { (Some(T), temperature_exponent-1.) }
		else { (None, temperature_exponent) };
	let temperature_factor = temperature_factor.map(|x| x.into());
	//let temperature_exponent_remainder = temperature_exponent_remainder;// as f32 as f64;
	//const T0: f64 = 1024.; // preexponential_factor overflow can be compensated when temperature_exponent is quite negative (TODO: add to exp argument instead of multiplying)
	//eprintln!("{temperature_exponent}, {temperature_exponent_remainder} {activation_temperature}");
	let x = Σ([
		(temperature_exponent_remainder as f32 != 0.).then(|| if temperature_exponent_remainder as f32==1. { lnT } else { f.def(temperature_exponent_remainder * (lnT/*-f64::ln(T0)*/), "β'lnT") }),
		(activation_temperature != 0.).then(|| f.def((-activation_temperature) * rcpT, "-Ea/kT"))
	].into_iter().filter_map(|x| x));
	let mul = |a:Option<Value>,b| if let Some(a) = a { a*b } else { b };
	let k = mul(temperature_factor, if let Some(x) = x {let e=exp(x+f64::ln(A), f); f.def(e, "exp(β'lnT-Ea/kT)").into()} else { f64(A).unwrap().into() });
	//.map(|x| {let e=min(f32::MAX as f64, exp(x, f)); f.def(e, "exp(β'lnT-Ea/kT)").into()}); // Saturates to MAX to resolve any inf-inf in dT_T to an arbitrary finite value
	//assert!(A!=1. || temperature_factor.is_some() || x.is_some(), "{A} {temperature_exponent} {activation_temperature} {temperature_factor:?} {exp:?}");
	min(f32::MAX as f64, k).expr() // Saturates to MAX to resolve any k*0 to 0
}

fn efficiency(efficiencies: &[f64], concentration: Value, concentrations: &[Value], f: &mut Block) -> Expression {
	use iter::Zip;
	let non_one_efficiencies:Box<[Value]> = list(efficiencies.zip(concentrations).into_iter().filter_map(|(&e, C)| if e==1. { None } else if e-1. == 1. { Some(*C) } else { Some(f.def((e-1.)*C,"fC")) }));
	if non_one_efficiencies.is_empty() { concentration.into() } else { concentration + sum(non_one_efficiencies.into_vec()) }
}

use super::*;

fn forward_rate_constant(model: &ReactionModel, k_inf: &RateConstant, T: T, concentration: Value, concentrations: &[Value], f: &mut Block) -> Ratio {
	use ReactionModel::*; match model {
		Elementary|Irreversible => edef(arrhenius(k_inf, T, f), f, "k").into(),
		ThreeBody{efficiencies} => (chk(arrhenius(k_inf, T, f), f) * efficiency(efficiencies, concentration, concentrations, f)).into(),
		PressureModification{efficiencies, k0} => {
			let efficiency = efficiency(efficiencies, concentration, concentrations, f);
			let k0 = arrhenius(k0, T, f);
			let C_k0 = f.def(efficiency * k0, "C_k0");
			let k_inf = edef(arrhenius(k_inf, T, f), f, "k_inf");
			Ratio(C_k0 * k_inf.shallow(), C_k0 + k_inf + (f32::MIN_POSITIVE as f64)) // Resolves undetermined form 0/0 as 0
		},
		Troe{efficiencies, k0, troe} => {
			let efficiency = efficiency(efficiencies, concentration, concentrations, f);
			let (k_inf, Pr) = {
				let (factor, ref k0, ref k_inf) = {
					let (mut k0, mut k_inf) = (*k0, *k_inf);
					let factor = k_inf.preexponential_factor;
					k0.preexponential_factor /= factor;
					//assert!((k0.preexponential_factor as f32).is_finite(), "{k0:?} {k_inf:?}");
					k_inf.preexponential_factor = 1.;
					(factor, k0, k_inf)
				};
				let k0 = {let e = arrhenius(k0, T, f); f.def(e, "k_0")};
				if k_inf.temperature_exponent == 0. && k_inf.activation_temperature == 0. {
					let Pr = f.def(efficiency * k0, "Pr");
					(f64(factor).unwrap().into(), Pr)
				} else {
					let k_inf = edef(arrhenius(k_inf, T, f), f, "k_inf");
					let Pr = f.def(efficiency * k0 / (k_inf.shallow() + f32::MIN_POSITIVE as f64), "Pr"); // Resolves undetermined form 0/0 as 0
					(factor * k_inf, Pr)
				}
			};
			let Fcent = {let model::Troe{A, T3, T1, T2} = *troe; let T{T,rcpT,..}=T; ast::sum([
				/*(T3 > 1e-30)*/true.then(|| { let y = 1.-A; if /*T3<1e30*/true { y * def(exp(T/(-T3), f),  f, "exp(-T/T3)") } else { y.into() }}), // skipping exp(-T/T3~1e-30) increases difference to Cantera from e-8 to e-3 on synthetic test with all mole fractions equals (including radicals)*/
				/*(T1 > 1e-30)*/true.then(|| { let y = A; if T1<1e30 { y * def(exp(T/(-T1), f), f, "exp(-T/T1)") } else { y.into() }}),
				(T2.is_finite()).then(|| def(exp((-T2)*rcpT, f), f, "exp(-T2/T)").into())
			].into_iter().filter_map(|x| x))};
			let lnFcent:Expr = if let Some(x) = Fcent.f64() { x.into() } else { def(ln(1./2., Fcent, f), f, "lnFcent").into() }; // 0.1-0.7 => e-3
			let C = -0.67*lnFcent.shallow() - 0.4*f64::ln(10.);
			let N = -1.27*lnFcent.shallow() + 0.75*f64::ln(10.);
			let lnPr𐊛C = def(ln(1., Pr + f32::MIN_POSITIVE as f64, f) + C, f, "lnPr+C"); // 2m - 2K // Resolves undetermined form -inf/-inf as 1
			let f1 = f.def(lnPr𐊛C / (-0.14*lnPr𐊛C+N), "f1");
			let F = exp(lnFcent/(f1*f1+1.), f);
			Ratio(F * k_inf * Pr, Pr + 1.)
		}
		SRI{/*efficiencies, k0, sri*/..} => {
			/*let efficiency = efficiency(efficiencies, concentration, concentrations, f);
			let k0 = {let e = arrhenius(k0, T, f); f.def(e, "k_0")};
			let k_inf = edef(arrhenius(k_inf, T, f), f, "k_inf");
			let Pr = f.def(efficiency * k0 / (k_inf.shallow() + f32::MIN_POSITIVE as f64), "Pr"); // Resolves undetermined form 0/0 as 0
			let logPr = ln(1., Pr + f32::MIN_POSITIVE as f64, f) ;
			let model::SRI{A, B, C, D, E} = *sri;
			let F = D*pow(A*exp(-B*rcpT)+exp((-1./C)*T), 1./(1.+logPr*logPr))*pow(T, E);
			Ratio(F * k_inf * Pr, Pr + 1.)*/unimplemented!()
		}
	}
}

fn reactions(species: &'t [NASA7], active: usize, Ts@T{T,rcpT,..}: T, concentration: Value, concentrations: &'t [Value], f: &mut Block) -> impl 't+FnMut(&Reaction, &mut Block)->Value {
	/*let mut notone = vec![false; species.len()];
	//fn efficiency(&Reaction{model, ..}) ->
	for Reaction{model, ..} in reactions {
		let efficiencies = {use ReactionModel::*; match model {
			Elementary|Irreversible => { continue; },
			ThreeBody{efficiencies} => efficiencies,
			PressureModification{efficiencies, ..} => efficiencies,
			Falloff{efficiencies, ..} => efficiencies,
		}};
		for (specie, &efficiency) in efficiencies.iter().enumerate() { if efficiency != 1. { notone[specie] = true; } }
	}
	eprintln!("{}", notone.iter().filter(|&&x| x).count());
	let mut add = 0; let mut sub = 0;
	for Reaction{model, ..} in reactions {
		let efficiencies = {use ReactionModel::*; match model {
			Elementary|Irreversible => { continue; },
			ThreeBody{efficiencies} => efficiencies,
			PressureModification{efficiencies, ..} => efficiencies,
			Falloff{efficiencies, ..} => efficiencies,
		}};
		for (specie, &efficiency) in efficiencies.iter().enumerate() {
			if notone[specie] { add+=1; }
			if efficiency != 1. { sub+=1; }
		}
	}
	eprintln!("{add} {}", sub+notone.iter().filter(|&&x| x).count());*/
	let rcp_C0 = f.def((1./NASA7::reference_pressure) * T, "RT/P0");
	let C0 = f.def(NASA7::reference_pressure * rcpT, "P0/RT");
	struct Cache {
		x: Box<[Option<Value>]>,
		exp: Box<[Option<Value>]>,
		exp_neg: Box<[Option<Value>]>
	}
	fn get(cache: &mut Box<[Option<Value>]>, k: usize, x: (&mut impl FnMut(usize,&mut Block)->Expression, &str), f: &mut Block) -> Value {
		*cache[k].get_or_insert_with(|| { let e = x.0(k,f); f.def(e, format!("{}{k}",x.1))})
	}
	fn max() -> f32 {
		let mut max = f32::cbrt(f32::MAX);
		while max*max*max > f32::MAX {
			use float_next_after::NextAfter;
			max = max.next_after(f32::NEG_INFINITY);
		}
		max
	}
	impl Cache {
		fn exp(&mut self, k: usize, (mut x,name): (impl Fn(usize,&mut Block)->Expression, &str), f: &mut Block) -> Value {
			get(&mut self.exp, k, (&mut |k,f|
				min(max() as f64, exp(get(&mut self.x,k,(&mut x, name),f), f)) // Saturates to MAX^⅓ to resolve (max^3).0 as 0
				, &format!("exp{}",name)), f)
		}
		fn exp_neg(&mut self, k: usize, (mut x,name): (impl Fn(usize,&mut Block)->Expression, &str), f: &mut Block) -> Value {
			get(&mut self.exp_neg, k, (&mut |k,f|
				min(f32::cbrt(f32::MAX) as f64, exp(-get(&mut self.x,k,(&mut x, name),f), f)) // Saturates to MAX^⅓ to resolve (max^3).0 as 0
				, &format!("exp_neg{}",name)), f)
		}
	}
	let mut cache = Cache{
		x: vec![None; active].into_boxed_slice(),
		exp: vec![None; active].into_boxed_slice(),
		exp_neg: vec![None; active].into_boxed_slice(),
	};
	move |Reaction{reactants, products, net, Σnet, rate_constant, model, ..}, f| {
		let forward_rate_constant = forward_rate_constant(&model, &rate_constant, Ts, concentration, concentrations, f);
		let Rforward = product_of_exponentiations(concentrations, reactants, f, "Rf");
		let Rnet = if let ReactionModel::Irreversible = model { Rforward } else {
			let Gibbs0_RT = |k, f: &mut Block| -> Expression {
				let NASA7{temperature_split, pieces} = &species[k];
				if temperature_split.is_nan() {
					self::thermodynamic(&self::Gibbs0_RT(&pieces[0]), Ts)
				} else {
					let piece = eval(|i| f.value(format!("g{k}_{i}")));
					f.statements.push(Statement::Select{
						condition: /*le!(f,*/(less_or_equal(T, temperature_split)),
						true_exprs: Box::new(self::Gibbs0_RT(&pieces[0]).map(|x| x.into())),
						false_exprs: Box::new(self::Gibbs0_RT(&pieces[1]).map(|x| x.into())),
						results: Box::new(piece)
					});
					self::thermodynamic(&piece, Ts)
				}
			};
			let rcp_equilibrium_constant_0 = if true {
				let Gibbs0_RT = (Gibbs0_RT, "Gibbs0_RT");
				net.iter().enumerate().fold(None, |mut p:Option</*Value*/Expression>, (k,&n)| {
					let t = if n>0 { cache.exp(k, Gibbs0_RT, f) } else { cache.exp_neg(k, Gibbs0_RT, f) };
					for _ in 0..i8::abs(n) { p = Some(if let Some(p) = p { p*t/*ast::def(p*t, "K0")*/ } else { t.into() }); }
					p
				}).unwrap()
			} else {
				min(f32::MAX as f64, exp(sum(net.iter().enumerate().filter(|(_,&n)| n!=0).map(|(k,&n)| (n as f64)*Gibbs0_RT(k, f))), f)) // Changes inf to MAX to resolve undetermined forms 0.inf as 0
			};
			let rcp_equilibrium_constant = match -Σnet { // reverse_rate_constant / forward_rate_constant
				0 => rcp_equilibrium_constant_0,
				1 => f.def(C0 * rcp_equilibrium_constant_0, "K").into(),
				-1 => f.def(rcp_C0 * rcp_equilibrium_constant_0, "K").into(),
				-2 => f.def(rcp_C0 * rcp_C0 * rcp_equilibrium_constant_0, "K").into(),
				_ => unreachable!(Σnet)
			};
			let Rreverse = {let e = rcp_equilibrium_constant * product_of_exponentiations(concentrations, products, f, "Rr"); f.def(e, "Rr/K")};
			f.def(Rforward - Rreverse, "R").into()
		};
		f.def(forward_rate_constant * Rnet, "cR")
	}
}

pub fn species_rates(species: &[NASA7], reactions: &[Reaction], Ts: T, concentration: Value, concentrations: &[Value], f: &mut Block, species_names: &[&str]) -> Box<[Value]> {
	let active = reactions[0].net.len();
	let mut species_rates: Box<[Option<Value>]> = vec![None; active].into_boxed_slice();
	let mut cR = self::reactions(species, active, Ts, concentration, concentrations, f);
	for reaction@Reaction{net, ..} in reactions {
		let cR = cR(reaction, f);
		for ((sum, &ν), specie) in zip(zip(&mut *species_rates, &**net), species_names).filter(|((_,&ν),_)| ν != 0) {
			let name = format!("rates_{specie}"); //ω̇
			*sum = Some(
				/**/   if ν == 1 { if let Some(sum) = sum { f.def(*sum+cR, name) } else { cR } }
        else if ν == -1 { if let Some(sum) = sum { f.def(*sum-cR, name) } else { f.def(-cR, name) } }
        else { let νcR = f.def((ν as f64)*cR, format!("νcR")); if let Some(sum) = sum { f.def(*sum+νcR, name).into() } else { νcR.into() } })
		}
	}
	map(&*species_rates, |x| x.unwrap())
}

pub fn reactions_rates(species: &[NASA7], reactions: &[Reaction]) -> Function {
	let active = reactions[0].net.len();
	let input@[pressure_R, T, volume, nonbulk_amounts @ ..] = &*map(0..(3+species.len()-1), Value) else { panic!() };
	let mut values = ["pressure_","rcp_pressure_", "T", "volume"].iter().map(|s| s.to_string()).chain((0..species.len()-1).map(|i| format!("active_amounts[{i}]"))).collect();
	let mut function = Block::new(&mut values);
	let ref mut f = function;
	let molar_density = f.def(1./volume, "molar_density");
	let rcpT = f.def(1./T, "rcpT");
	let lnT = def(ln(1024., T, f), f, "lnT");
	let T2 = f.def(T*T, "T2");
	let nonbulk_concentrations = map(nonbulk_amounts, |nonbulk_amount| f.def(molar_density*max(0., nonbulk_amount), "nonbulk_concentrations"));
	let T3 = f.def(T*T2, "T3");
	let concentration = f.def(pressure_R * rcpT, "concentration");
	let rcpT2 = f.def(rcpT*rcpT, "rcpT2");
	let T4 = f.def(T2*T2, "T4");
	let Ts = T{lnT,T: *T,T2,T3,T4,rcpT,rcpT2};
	let bulk_concentration = f.def(concentration - sum(&*nonbulk_concentrations), "bulk_concentration"); // Constant pressure
	let concentrations = list(nonbulk_concentrations.into_vec().into_iter().chain([bulk_concentration].into_iter()));
	let mut cR = self::reactions(species, active, Ts, concentration, &concentrations, f);
	Function{
		output: map(reactions, |r| cR(r, f).into()),
		statements: function.block.statements.into(),
		input: vec![Type::F32; input.len()].into(),
		values: values.into()
	}
}

pub fn rates(molar_mass: &[f64], species: &[NASA7], reactions: &[Reaction], species_names: &[&str]) -> Function {
	let active = reactions[0].net.len();
	let input@[pressure_R, rcp_pressure_R, T, volume, nonbulk_amounts @ ..] = &*map(0..(4+species.len()-1), Value) else { panic!() };
	let mut values = ["pressure_","rcp_pressure_", "T", "volume"].iter().map(|s| s.to_string()).chain((0..species.len()-1).map(|i| format!("active_amounts[{i}]"))).collect();
	let mut function = Block::new(&mut values);
	let ref mut f = function;
	let molar_density = f.def(1./volume, "molar_density");
	let rcpT = f.def(1./T, "rcpT");
	let lnT = def(ln(1024., T, f), f, "lnT");
	let T2 = f.def(T*T, "T2");
	let nonbulk_concentrations = map(nonbulk_amounts, |nonbulk_amount| f.def(molar_density*max(0., nonbulk_amount), "nonbulk_concentrations"));
	let T3 = f.def(T*T2, "T3");
	let concentration = f.def(pressure_R * rcpT, "concentration");
	let rcpT2 = f.def(rcpT*rcpT, "rcpT2");
	let T4 = f.def(T2*T2, "T4");
	let Ts = T{lnT,T: *T,T2,T3,T4,rcpT,rcpT2};
	let bulk_concentration = f.def(concentration - sum(&*nonbulk_concentrations), "bulk_concentration"); // Constant pressure
	let concentrations = list(nonbulk_concentrations.into_vec().into_iter().chain([bulk_concentration].into_iter()));
	let species_rates = species_rates(species, reactions, Ts, concentration, &concentrations, f, species_names);
	let [enthalpy_RT] = thermodynamics(&species[0..active], [enthalpy_RT], Ts, f, ["enthalpy_RT"]);
	use iter::Dot;
	let energy_rate_RT : Expression = (&species_rates).dot(enthalpy_RT.into_vec());
	let [molar_heat_capacity_at_CP_R] = thermodynamics(species, [molar_heat_capacity_at_constant_pressure_R], Ts, f, ["molar_heat_capacity_at_CP_R"]);
	let Cp: Expression = molar_heat_capacity_at_CP_R.dot(concentrations);
	let dtT_T = f.def(-energy_rate_RT/Cp, "dT_T");
	let bulk_molar_mass = molar_mass.last().unwrap();
	let reduced_molar_masses = map(&*molar_mass, |w| 1.-w/bulk_molar_mass); // Ensures conservation of mass (transmutes with bulk specie instead)
	let dtE: Expression = reduced_molar_masses.dot(&species_rates);
	let dtV = volume * (dtT_T + rcp_pressure_R*T*dtE);
	Function{
		output: list([T * dtT_T, dtV/*f64(0.).unwrap().into(), f64(0.).unwrap().into()*/].into_iter().chain(species_rates.into_vec().into_iter().map(|e| e.into()))),
		statements: function.block.statements.into(),
		input: vec![Type::F32; input.len()].into(),
		values: values.into()
	}
}
