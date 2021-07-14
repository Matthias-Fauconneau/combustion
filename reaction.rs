//#![feature(associated_type_bounds,bindings_after_at,array_map,format_args_capture,trait_alias)]#![allow(uncommon_codepoints,confusable_idents,non_snake_case)]
fn bucket<I:IntoIterator<Item:Eq>>(iter: I) -> impl IntoIterator<Item=(I::Item, Vec<usize>)> {
	let mut map = linear_map::LinearMap::<_, Vec<_>>::new();
	for (index, key) in iter.into_iter().enumerate() { map.entry(key).or_insert(Default::default()).push(index) }
	map
}

use ast::{*/*, dbg*/};

fn product_of_exponentiations(b: &[Value], n: &[impl Copy+Into<i16>]) -> Expression {
	let (num, div) : (Vec::<_>,Vec::<_>) = n.iter().map(|&n| n.into()).zip(b).filter(|&(n,_)| n!=0).partition(|&(n,_)| n>0);
	let num : Option<Expression> = num.into_iter().fold(None, |mut p, (n,e)|{ for _ in 0..n { p = Some(match p { Some(p) => p*e, None => e.into() }); } p });
	let div : Option<Expression> = div.into_iter().fold(None, |mut p, (n,e)|{ for _ in 0..-n { p = Some(match p { Some(p) => p*e, None => e.into() }); } p });
	match (num, div) {
		(None, None) => None,
		(Some(num), None) => Some(num),
		(None, Some(div)) => Some(1./div),
		(Some(num), Some(div)) => Some(num/div)
	}.unwrap()
}

type Error = String;
use {fehler::throws, iter::{list, map}, super::*};

#[derive(Clone, Copy)] struct T<'t> { ln_T: &'t Value, T: &'t Value, T2: &'t Value, T3: &'t Value, T4: &'t Value, rcp_T: &'t Value, rcp_T2: &'t Value}
fn thermodynamics(thermodynamics: &[NASA7], expression: impl Fn(&[f64], T<'_>, &mut Block)->Expression, T: T<'_>, f: &mut Block, debug: &str) -> Box<[Value]> {
	let mut specie_results = map(thermodynamics, |_| None);
	for (temperature_split, ref species) in bucket(thermodynamics.iter().map(|s| s.temperature_split.to_bits())) {
		let results = map(species, |specie| f.value(format!("{debug}[{specie}]")));
		for (&specie, result) in species.iter().zip(&*results) { assert!(specie_results[specie].replace(result.clone()).is_none()) }
		push(Statement::Select{
			condition: less_or_equal(T.T, f64::from_bits(temperature_split)),
			true_exprs: map(species, |&specie| expression(&thermodynamics[specie].pieces[0], T, f)),
			false_exprs: map(species, |&specie| expression(&thermodynamics[specie].pieces[1], T, f)),
			results
		}, f);
	}
	map(specie_results.into_vec().into_iter(), Option::unwrap)
}

fn molar_heat_capacity_at_constant_pressure_R(a: &[f64], T{T,T2,T3,T4,..}: T, _: &mut Block) -> Expression { a[0]+a[1]*T+a[2]*T2+a[3]*T3+a[4]*T4 }
fn enthalpy_RT(a: &[f64], T{T,T2,T3,T4,rcp_T,..}: T<'_>, _: &mut Block) -> Expression { a[0]+a[1]/2.*T+a[2]/3.*T2+a[3]/4.*T3+a[4]/5.*T4+a[5]*rcp_T }
fn Gibbs_RT(a: &[f64], T{ln_T,T,T2,T3,T4,rcp_T,..}: T<'_>, _: &mut Block) -> Expression {
	(a[0]-a[6])-a[0]*ln_T-a[1]/2.*T+(1./3.-1./2.)*a[2]*T2+(1./4.-1./3.)*a[3]*T3+(1./5.-1./4.)*a[4]*T4+a[5]*rcp_T
}
fn exp_Gibbs_RT(a: &[f64], T: T<'_>, f: &mut Block) -> Expression { exp(Gibbs_RT(a, T, f), f) }

// A.T^β.exp(-Ea/kT)
#[track_caller]#[throws] fn arrhenius(&RateConstant{preexponential_factor: A, temperature_exponent, activation_temperature}: &RateConstant, T{ln_T,T,T2,T3,T4,rcp_T,rcp_T2}: T<'_>, f: &mut Block) -> Expression {
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
	const T0: f64 = 1024.;
	[Some(float(A*f64::powf(T0,ln_T_coefficient)).ok_or(format!("A:{A:e}"))?), temperature_factor.map(|x| x.into()), [
	 (ln_T_coefficient != 0.).then(|| ln_T_coefficient * (ln_T-f64::ln(T0))),
	 (activation_temperature != 0.).then(|| -activation_temperature * rcp_T)
	].into_iter().filter_map(|x| x).sum::<Option<_>>().map(|x| exp(x, f))].into_iter().filter_map(|x| x).product::<Option<_>>().unwrap()
}

fn forward_rate_constant(model: &ReactionModel, k_inf: &RateConstant, T: T, concentrations: &[Value], f: &mut Block) -> Expression {
	use ReactionModel::*; match model {
		Elementary|Irreversible => arrhenius(k_inf, T, f).unwrap(),
		ThreeBody{efficiencies} => arrhenius(k_inf, T, f).unwrap() * dot(efficiencies, concentrations),
		PressureModification{efficiencies, k0} => f.block(|f|{
			let ref C_k0 = l!(f dot(efficiencies, concentrations) * arrhenius(k0, T, f).unwrap());
			let ref k_inf = l!(f arrhenius(k_inf, T, f).unwrap());
			(C_k0 * k_inf) / (C_k0 + k_inf)
		}),
		Falloff{efficiencies, k0, troe} => {let Troe{A, T3, T1, T2} = *troe;/*ICE inside*/ f.block(|f|{
			let k0 = arrhenius(k0, T, f).expect(&format!("{k0:?}/{k_inf:?}"));
			let ref k_inf = l!(f arrhenius(k_inf, T, f).unwrap());
			let ref Pr = l!(f dot(efficiencies, concentrations) * k0 / k_inf);
			let Fcent = {let T{T,rcp_T,..}=T; sum([
				(T3 > 1e-30).then(|| { let y = 1.-A; if T3<1e30 { y * exp(T/(-T3), f) } else { y.into() }}),
				(T1 > 1e-30).then(|| { let y = A; if T1<1e30 { y * exp(T/(-T1), f) } else { y.into() }}),
				(T2 != 0.).then(|| exp((-T2)*rcp_T, f))
			].into_iter().filter_map(|x| x))};
			let ref lnFcent = l!(f ln(1./2., Fcent, f)); // 0.1-0.7 => e-3
			let C =-0.67*lnFcent - 0.4*f64::ln(10.);
			let N = -1.27*lnFcent + 0.75*f64::ln(10.);
			let ref lnPr𐊛C = l!(f ln(1., Pr, f) + C); // 2m - 2K
			let ref f1 = l!(f lnPr𐊛C / (-0.14*lnPr𐊛C+N));
			let F = exp(lnFcent/(f1*f1+1.), f);
			k_inf * Pr / (Pr + 1.) * F
		})}
	}
}

fn reaction_rates(reactions: &[Reaction], T: T, C0: &Value, rcp_C0: &Value, exp_Gibbs0_RT: &[Value], concentrations: &[Value], f: &mut Block) -> Box<[Value]> {
	map(reactions.iter().enumerate(), |(_i, Reaction{reactants, products, net, Σnet, rate_constant, model, ..})| {
		let forward_rate_constant = forward_rate_constant(model, rate_constant, T, concentrations, f); // todo: CSE
		let forward = product_of_exponentiations(concentrations, reactants);
		let coefficient = if let ReactionModel::Irreversible = model { forward } else {
			let rcp_equilibrium_constant_0 = product_of_exponentiations(exp_Gibbs0_RT, net);
			//dbg!(f; rcp_equilibrium_constant_0);
			//let rcp_equilibrium_constant_0 = exp2(idot(net.iter().map(|&net| net as f64).zip(Gibbs0_RT)), f);
			let rcp_equilibrium_constant = match -Σnet { // reverse_rate_constant / forward_rate_constant
				0 => rcp_equilibrium_constant_0.into(),
				1 => C0 * rcp_equilibrium_constant_0,
				-1 => rcp_C0 * rcp_equilibrium_constant_0,
				_ => unreachable!()
			};
			//dbg!(f; rcp_equilibrium_constant);
			let reverse = rcp_equilibrium_constant * product_of_exponentiations(concentrations, products);
			//dbg!(f; forward, reverse);
			forward - reverse
		};
		//dbg!(f; forward_rate_constant, coefficient);
		l!(f forward_rate_constant * coefficient)
	})
}

pub fn rates(species: &[NASA7], reactions: &[Reaction]) -> Function {
	let active = {
		let ref mut iter = (0..species.len()-1).map(|k| reactions.iter().any(|Reaction{net,..}| net[k] != 0));
		let active = iter.take_while(|&is_active| is_active).count();
		assert!(iter.all(|is_active| !is_active));
		active
	};
	let_!{ input@[ref pressure_R, ref total_amount, ref T, ref nonbulk_amounts @ ..] = &*map(0..(3+species.len()-1), Value) => {
	let mut values = ["pressure_","total_amount","T"].iter().map(|s| s.to_string()).chain((0..species.len()-1).map(|i| format!("active_amounts[{i}]"))).collect();
	let mut function = Block::new(&mut values);
	let ref mut f = function;
	let ref ln_T = l!(f ln(1024., T, f));
	let ref T2 = l!(f T*T);
	let ref T3 = l!(f T*T2);
	let ref T4 = l!(f T*T3);
	let ref rcp_T = l!(f 1./T);
	let ref rcp_T2 = l!(f num::sq(rcp_T));
	let ref rcp_C0 = l!(f (1./NASA7::reference_pressure) * T);
	let ref C0 = l!(f NASA7::reference_pressure * rcp_T);
	let ref total_concentration = l!(f pressure_R / T); // Constant pressure
	let T = T{ln_T,T,T2,T3,T4,rcp_T,rcp_T2};
	let ref exp_Gibbs0_RT = thermodynamics(&species[0..active], exp_Gibbs_RT, T, f, "exp_Gibbs0_RT");
	let ref density = l!(f total_concentration / total_amount);
	let nonbulk_concentrations = map(0..active, |k| l!(f density*max(0., &nonbulk_amounts[k])));
	let bulk_concentration = l!(f total_concentration - sum(&*nonbulk_concentrations));
	let concentrations = list(nonbulk_concentrations.into_vec().into_iter().chain([bulk_concentration].into_iter()));
	let rates = reaction_rates(reactions, T, C0, rcp_C0, exp_Gibbs0_RT, &concentrations, f);
	let rates = map(0..active, |specie| l!(f idot(reactions.iter().map(|Reaction{net, ..}| net[specie] as f64).zip(&*rates))));
	let enthalpy_RT = thermodynamics(&species[0..active], enthalpy_RT, T, f, "enthalpy_RT");
	let dot = |a:&[Value], b:&[Value]| iter::dot(a.iter().zip(b.iter()));
	let energy_rate_RT : Expression = dot(&rates, &enthalpy_RT);
	let Cp : Expression = dot(&concentrations, &thermodynamics(species, molar_heat_capacity_at_constant_pressure_R, T, f, "molar_heat_capacity_at_CP_R"));
	let dtT_T = - energy_rate_RT / Cp;
	Function{output: list([T.T * dtT_T].into_iter().chain(rates.iter().map(|v| v.into()))), statements: function.statements.into(), input: input.len(), values: values.into()}
}}}
