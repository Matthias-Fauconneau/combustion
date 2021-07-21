//#![feature(associated_type_bounds,bindings_after_at,array_map,format_args_capture,trait_alias,array_zip)]#![allow(uncommon_codepoints,confusable_idents,non_snake_case)]
fn bucket<I:IntoIterator<Item:Eq>>(iter: I) -> impl std::iter::IntoIterator<Item=(I::Item, Vec<usize>)> {
	let mut map = linear_map::LinearMap::<_, Vec<_>>::new();
	for (index, key) in iter.into_iter().enumerate() { map.entry(key).or_insert(Default::default()).push(index) }
	map
}

use {fehler::throws, std::default::default, std::ops::Deref, iter::{Prefix, generate, list, map, DotN}, std::collections::hash_map::{HashMap,Entry::*}, ast::*};

#[derive(derive_more::Deref,derive_more::DerefMut)] pub struct Block<'t> {
	#[deref]#[deref_mut] block: ast::Block<'t>,
	values: HashMap<Expression, Value>,
	expressions: Vec<(Expr, String)>,
	after_CSE: std::collections::HashSet<Expr>,
}
#[throws(f64)] pub fn def(e: impl Into<Expression>, f: &mut Block, name: String) -> Value {
	match f.values.entry(e.into()) {
		Occupied(entry) => *entry.get(),
		Vacant(entry) => {
			let e = ast::def(entry.key().clone(), &mut f.block, name)?;
			*entry.insert(e)
		},
	}
}
#[macro_export] macro_rules! le {
	($f:ident $e:expr) => (def($e, $f, format!("{}:{}: {}", file!(), line!(), stringify!($e))).unwrap());
	($f:ident; $e:expr) => (def($e, $f, format!("{}:{}: {}", file!(), line!(), stringify!($e))).map(|v| v.into()).unwrap_or_else(|e| e.into()))
}

#[track_caller] fn check(e: &Expression, f: &mut Block) -> Result<(),String> {
	{
		let e = e.deref();
		if !e.is_leaf() {
			let new = f.after_CSE.insert(e.clone());
			if !new { return Err(format!("{}: {}", e.to_string(f.names), {use itertools::Itertools; f.expressions.iter().filter_map(|(k,v)| (k==e).then(||v)).format(" ")})); }
		}
	}
	match e.visit(|e| check(e, f)) {
		[Some(Err(e)), _]|[_, Some(Err(e))] => Err(e),
		_  => Ok(()),
	}
}

#[track_caller] fn chk(e: Expression, f: &mut Block) -> Expression { check(&e, f).unwrap(); e }

struct Ratio(Expression, Expression);
impl<E:Into<Expression>> From<E> for Ratio { fn from(e: E) -> Ratio { Ratio(e.into(), (1.).into()) } }
impl std::ops::Mul<Ratio> for Expr { type Output = Expression; fn mul(self, Ratio(n, d): Ratio) -> Expression {
	if let Some(d) = d.f64() { if d == 1. { return self*n } }
	(self*n)/d
}}
impl std::ops::Mul<Ratio> for &Value { type Output = Expression; fn mul(self, r: Ratio) -> Expression { self.into():Expr*r } }

pub fn sum(iter: impl IntoIterator<Item:Into<Expression>>, f: &mut Block) -> Option<Expression> {
	iter.into_iter().map(|e| e.into()).filter(|e| if let Some(x) = e.f32() {x!=0.} else {true}).reduce(|a,b| le!(f; a+b))
}

fn product_of_exponentiations_<N:Into<i16>>(iter: impl IntoIterator<Item=(&'t Value, N)>, f: &mut Block) -> Expr {
	let (num, div) : (Vec::<_>,Vec::<_>) = iter.into_iter().map(|(v,n)| (v,n.into())).filter(|&(_,n)| n!=0).partition(|&(_,n)| n>0);
	let num = num.into_iter().fold(None, |mut p:Option<Expr>, (v,n)|{ for _ in 0..n { p = Some(match p { Some(p) => le!(f; p*v), None => v.into() }); } p });
	let div = div.into_iter().fold(None, |mut p:Option<Expr>, (v,n)|{ for _ in 0..-n { p = Some(match p { Some(p) => le!(f; p*v), None => v.into() }); } p });
	match (num, div) {
		(None, None) => None,
		(Some(num), None) => Some(num),
		(None, Some(div)) => Some(l!(f; 1./div)),
		(Some(num), Some(div)) => Some(le!(f; num/div))
	}.unwrap()
}
fn product_of_exponentiations__<N:iter::IntoExactSizeIterator>(v: &[Value], n: N, f: &mut Block) -> Expr where <N as iter::IntoIterator>::Item:Into<i16> {
	use iter::Zip;
	product_of_exponentiations_(v.zip(n), f)
}
fn product_of_exponentiations<I:iter::IntoExactSizeIterator,N:Copy+Into<i16>>(v: &[Value], n: I, f: &mut Block) -> Expr where <I as iter::IntoIterator>::Item:Deref<Target=N> {
	use iter::Map;
	product_of_exponentiations__(v, n.map(|x:<I as iter::IntoIterator>::Item| *x), f)
}

//pub fn dotv(iter: impl IntoIterator<Item=(f64, impl Into<Expression>)>) -> Option<Expression> { iter.into_iter().map(|(c,e)| c*e.into()).sum() }
//#[track_caller] pub fn dotv(c: &[f64], v: impl IntoIterator<Item:Into<Expression>>) -> Option<Expression> { zdotv(c.iter().copied().zip(v)) }

#[derive(Clone,Copy)] struct T { ln_T: Value, T: Value, T2: Value, T3: Value, T4: Value, rcp_T: Value, rcp_T2: Value}

fn molar_heat_capacity_at_constant_pressure_R(a: &[f64; 7], T{T,T2,T3,T4,..}: T) -> Expression { (*a.prefix()).dot([Expr::from(1.),T.into(),T2.into(),T3.into(),T4.into()]) }
/*fn enthalpy_RT(a: &[f64; 7], T{T,T2,T3,T4,rcp_T,..}: T) -> Expression { a[0] + a[1]/2.*T + a[2]/3.*T2 + a[3]/4.*T3 + a[4]/5.*T4 + a[5]*rcp_T }
fn Gibbs_RT(a: &[f64; 7], T{ln_T,T,T2,T3,T4,rcp_T,..}: T) -> Expression { -a[0]*ln_T +a[0]-a[6] -a[1]/2.*T +(1./3.-1./2.)*a[2]*T2 +(1./4.-1./3.)*a[3]*T3 +(1./5.-1./4.)*a[4]*T4 +a[5]*rcp_T }*/
fn a1T(a: &[f64; 7], T{T,..}: T) -> Expression { a[1]/2.*T }
fn a5rcpT(a: &[f64; 7], T{rcp_T,..}: T) -> Expression { a[5]*rcp_T }
fn enthalpy_RT(a: &[f64; 7], T{T2,T3,T4,..}: T) -> Expression { a[0] + a[2]/3.*T2 + a[3]/4.*T3 + a[4]/5.*T4 }
fn Gibbs0_RT(a: &[f64; 7], T{ln_T,T2,T3,T4,..}: T) -> Expression { -a[0]*ln_T +(a[0]-a[6]) +(1./3.-1./2.)*a[2]*T2 +(1./4.-1./3.)*a[3]*T3 +(1./5.-1./4.)*a[4]*T4 }

fn thermodynamics<const N: usize>(thermodynamics: &[NASA7], expressions: [impl Fn(&[f64; 7], T)->Expression; N], Ts@T{T,..}: T, f: &mut Block, debug: [&str; N]) -> [Box<[Expression]>; N] {
	use iter::IntoConstSizeIterator;
	let mut specie_results = generate(|_| map(thermodynamics, |_| None)).collect();
	for (temperature_split, ref species) in bucket(thermodynamics.iter().map(|s| ordered_float::OrderedFloat(s.temperature_split))) {
		if temperature_split.is_nan() {
			for (expression, specie_results) in expressions.iter().zip(specie_results.iter_mut()) { for &specie in species {
				let e = expression(&thermodynamics[specie].pieces[0], Ts);
				check(&e, f).unwrap();
				let e = if e.is_leaf() { e } else { l!(f; e) };
				assert!(specie_results[specie].replace(e).is_none());
			}}
		} else {
			let (exprs, results):(Vec<_>,Vec<_>) = expressions.iter().zip(specie_results.iter_mut()).zip(debug).map(|((expression, specie_results), debug)| {
				let results = map(species, |specie| f.value(format!("{debug}[{specie}]")));
				use iter::Zip;
				for (&specie, result) in species.zip(&*results) { assert!(specie_results[specie].replace(result.into()).is_none()) }
				let mut true_exprs = map(species, |&specie| expression(&thermodynamics[specie].pieces[0], Ts));
				let mut false_exprs = map(species, |&specie| expression(&thermodynamics[specie].pieces[1], Ts));
				eliminate_common_subexpressions(&mut true_exprs, &mut false_exprs, f);
				for e in true_exprs.iter().chain(&*false_exprs) { check(e, f).unwrap(); }
				([true_exprs, false_exprs], results)
			}).unzip();
			let (true_exprs, false_exprs):(Vec<_>,Vec<_>) = exprs.into_iter().map(|[a,b]| (a,b)).unzip();
			let (true_exprs, false_exprs, results) = (true_exprs.concat().into(), false_exprs.concat().into(), results.concat().into());
			push(Statement::Select{condition: le!(f; less_or_equal(T, f64::from(temperature_split))), true_exprs, false_exprs, results}, f);
		}
	}
	specie_results.map(|specie_results| map(specie_results.into_vec().into_iter(), Option::unwrap))
}

type Error = String;

// A.T^β.exp(-Ea/kT)
#[throws] fn arrhenius(&RateConstant{preexponential_factor: A, temperature_exponent, activation_temperature}: &RateConstant, T{ln_T,T,T2,T3,T4,rcp_T,rcp_T2}: T, f: &mut Block) -> Expression {
	let (temperature_factor, ln_T_coefficient) =
		/*if temperature_exponent <= -8. { unimplemented!("{temperature_exponent}") }
		// T~1000: T^-n << 1 so no need to factorize exponent out of exp to reduce input domain
		else*/ if temperature_exponent <= -1.5 { (Some(rcp_T2), temperature_exponent+2.) }
		else if temperature_exponent <= -0.5 { (Some(rcp_T), temperature_exponent+1.) }
		else if temperature_exponent >= 8. { unimplemented!("{temperature_exponent}") }
		// TODO: T~1000: factorize exponent out of exp to reduce input domain
		else if temperature_exponent >= 3.5 { (Some(T4), temperature_exponent-4.) }
		else if temperature_exponent >= 2.5 { (Some(T3), temperature_exponent-3.) }
		else if temperature_exponent >= 1.5 { (Some(T2), temperature_exponent-2.) }
		else if temperature_exponent >= 0.5 { (Some(T), temperature_exponent-1.) }
		else { (None, temperature_exponent) };
	let ln_T_coefficient = ln_T_coefficient as f32 as f64;
	//const T0: f64 = 1024.;
	//eprintln!("{temperature_exponent}, {ln_T_coefficient}");
	let exp = [
			(ln_T_coefficient != 0.).then(|| le!(f; ln_T_coefficient * (ln_T/*-f64::ln(T0)*/))),
			(activation_temperature != 0.).then(|| le!(f; (-activation_temperature) * rcp_T))
		].into_iter().filter_map(|x:Option<Expression>| x).sum::<Option<Expression>>()
		.map(|x| le!(f; exp(x, f)));
	Π([
		Some(f64(A).ok_or_else(|| format!("{A:e}"))?/* *f64::powf(T0,ln_T_coefficient)*/),
		temperature_factor.map(|x| x.into()),
		exp
	]).unwrap()
}

fn efficiency(efficiencies: &[f64], concentrations: &[Value], f: &mut Block) -> Expression {
	use iter::Zip;
	let efficiencies:Box<[Expr]> = list(efficiencies.zip(concentrations).into_iter().filter_map(|(&e, C)| if e==0. { None } else if e==1. { Some(C.into()) } else { Some(le!(f; e*C)) }));
	sum(efficiencies.into_vec(), f).unwrap()
}

use super::*;

fn forward_rate_constant(model: &ReactionModel, k_inf: &RateConstant, T: T, concentrations: &[Value], f: &mut Block) -> Ratio {
	use ReactionModel::*; match model {
		Elementary|Irreversible => le!(f; arrhenius(k_inf, T, f).unwrap()),
		ThreeBody{efficiencies} => (chk(arrhenius(k_inf, T, f).unwrap(), f) * efficiency(efficiencies, concentrations, f)).into(),
		PressureModification{efficiencies, k0} => {
			let efficiency = efficiency(efficiencies, concentrations, f);
			let C_k0 = l!(f efficiency * chk(arrhenius(k0, T, f).unwrap(), f));
			let k_inf = l!(f chk(arrhenius(k_inf, T, f).unwrap(), f));
			Ratio(C_k0 * k_inf, C_k0 + k_inf)
		},
		Falloff{efficiencies, k0, troe} => {
			let efficiency = efficiency(efficiencies, concentrations, f);
			let k0 = chk(arrhenius(k0, T, f).expect(&format!("{k0:?}/{k_inf:?}")), f);
			let k_inf = (|e:Expression,f|->Expression{ if let Some(x) = e.f64() { x.into() } else { l!(f; e) }})(chk(arrhenius(k_inf, T, f).unwrap(), f), f);
			//f.block(|f|{
				let Pr = l!(f efficiency * k0 / k_inf.shallow());
				let Fcent = {let Troe{A, T3, T1, T2} = *troe; let T{T,rcp_T,..}=T; ast::sum([
					(T3 > 1e-30).then(|| { let y = 1.-A; if T3<1e30 { y * le!(f exp(T/(-T3), f)) } else { y.into() }}),
					(T1 > 1e-30).then(|| { let y = A; if T1<1e30 { y * le!(f exp(T/(-T1), f)) } else { y.into() }}),
					(T2.is_finite()).then(|| exp((-T2)*rcp_T, f))
				].into_iter().filter_map(|x| x))};
				let lnFcent: Expr = if let Some(x) = Fcent.f64() { x.into() } else { l!(f; ln(1./2., Fcent, f)) }; // 0.1-0.7 => e-3
				let C = -0.67*lnFcent.shallow() - 0.4*f64::ln(10.);
				let N = -1.27*lnFcent.shallow() + 0.75*f64::ln(10.);
				let lnPr𐊛C = l!(f ln(1., Pr, f) + C); // 2m - 2K
				let f1 = l!(f lnPr𐊛C / (-0.14*lnPr𐊛C+N));
				let F = exp(lnFcent/(f1*f1+1.), f);
				Ratio(F * k_inf * Pr, Pr + 1.)
			//})
		}
	}
}

fn reaction_rates(reactions: &[Reaction], T: T, C0: &Value, rcp_C0: &Value, exp_Gibbs0_RT: &[Value], concentrations: &[Value], f: &mut Block) -> Box<[Value]> {
	map(reactions.iter().enumerate(), |(_i, Reaction{reactants, products, net, Σnet, rate_constant, model, ..})| {
		let forward_rate_constant = forward_rate_constant(model, rate_constant, T, concentrations, f);
		let forward = product_of_exponentiations(concentrations, reactants, f);
		let coefficient = if let ReactionModel::Irreversible = model { forward } else {
			let rcp_equilibrium_constant_0 = product_of_exponentiations(exp_Gibbs0_RT, net, f);
			//let rcp_equilibrium_constant_0 = exp(dot(net.iter().map(|&net| net as f64).zip(Gibbs0_RT)), f);
			let rcp_equilibrium_constant = match -Σnet { // reverse_rate_constant / forward_rate_constant
				0 => rcp_equilibrium_constant_0,
				1 => le!(f; C0 * rcp_equilibrium_constant_0),
				-1 => le!(f; rcp_C0 * rcp_equilibrium_constant_0),
				_ => unreachable!()
			};
			let reverse = le!(f rcp_equilibrium_constant * product_of_exponentiations(concentrations, products, f));
			le!(f; forward - reverse)
		};
		le!(f coefficient * forward_rate_constant)
	})
}

pub fn rates(species: &[NASA7], reactions: &[Reaction]) -> Function {
	let active = reactions[0].net.len();
	let_!{ input@[ref pressure_R, ref total_amount, ref T, ref nonbulk_amounts @ ..] = &*map(0..(3+species.len()-1), Value) => {
	let mut values = ["pressure_","total_amount","T"].iter().map(|s| s.to_string()).chain((0..species.len()-1).map(|i| format!("active_amounts[{i}]"))).collect();
	let mut function = Block{block: ast::Block::new(&mut values), values: default(), expressions: vec![], after_CSE: default()};
	let ref mut f = function;
	let ln_T = l!(f ln(1024., T, f));
	let T2 = l!(f T*T);
	let T3 = l!(f T*T2);
	let T4 = l!(f T2*T2);
	let rcp_T = l!(f 1./T);
	let rcp_T2 = l!(f rcp_T*rcp_T);
	let rcp_C0 = l!(f chk((1./NASA7::reference_pressure) * T, f));
	let C0 = l!(f NASA7::reference_pressure * rcp_T);
	let total_concentration = l!(f pressure_R / T); // Constant pressure
	let Ts = T{ln_T, T: *T,T2,T3,T4,rcp_T,rcp_T2};
	let [a1T, a5rcpT, Gibbs0_RT] = thermodynamics(&species[0..active], [a1T as fn(&[f64;7],T)->Expression, a5rcpT, Gibbs0_RT], Ts, f, ["a1T","a5/T","Gibbs0/RT"]);
	let Gibbs0_RT = map(a1T.iter().zip(&*a5rcpT).zip(Gibbs0_RT.into_vec()), |((a1T,a5rcpT),g)| -a1T.shallow()+g+a5rcpT.shallow());
	let exp_Gibbs0_RT = map(Gibbs0_RT.into_vec(), |g| l!(f exp(g, f)));
	let density = l!(f total_concentration / total_amount);
	let nonbulk_concentrations = map(0..active, |k| def(density*max(0., &nonbulk_amounts[k]),f, format!("C[{k}]")).unwrap());
	let bulk_concentration = l!(f total_concentration - sum(&*nonbulk_concentrations, f).unwrap());
	let concentrations = list(nonbulk_concentrations.into_vec().into_iter().chain([bulk_concentration].into_iter()));
	let rates = reaction_rates(reactions, Ts, &C0, &rcp_C0, &exp_Gibbs0_RT, &concentrations, f);
	pub fn dot(iter: impl IntoIterator<Item=(f64, impl Into<Expression>)>, f: &mut Block) -> Option<Expression> {
		iter.into_iter().fold(None, |sum, (c, e)| {
			let e = e.into();
			Some(if c == 1. { sum+e }
			else if c == -1. { if let Some(sum) = sum { sum-e } else { le!(f; -e) } } // fixme: reorder -a+b -> b-a to avoid some neg
			else { le!(f; sum+le!(f; c*e)) })
		})
	}
	let ref rates = map(0..active, |specie| l!(f dot(reactions.iter().map(|Reaction{net, ..}| net[specie] as f64).zip(&*rates), f).unwrap()));
	let [enthalpy_RT] = thermodynamics(&species[0..active], [enthalpy_RT], Ts, f, ["enthalpy_RT"]);
	let enthalpy_RT = map(a1T.into_vec().into_iter().zip(a5rcpT.into_vec().into_iter()).zip(enthalpy_RT.into_vec()), |((a1T,a5rcpT),h)| a1T+h+a5rcpT);
	use iter::Dot;
	let energy_rate_RT : Expression = rates.dot(enthalpy_RT.into_vec());
	let [molar_heat_capacity_at_CP_R] = thermodynamics(species, [molar_heat_capacity_at_constant_pressure_R], Ts, f, ["molar_heat_capacity_at_CP_R"]);
	let Cp : Expression = molar_heat_capacity_at_CP_R.dot(concentrations);
	let dtT_T = Ratio(- energy_rate_RT, Cp);
	Function{output: list([T * dtT_T].into_iter().chain(rates.iter().map(|v:&Value| v.into()))), statements: function.block.statements.into(), input: input.len(), values: values.into()}
}}}
