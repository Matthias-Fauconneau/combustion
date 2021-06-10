fn bucket<I:IntoIterator<Item:Eq>>(iter: I) -> impl IntoIterator<Item=(I::Item, Vec<usize>)> {
	let mut map = linear_map::LinearMap::<_, Vec<_>>::new();
	for (index, key) in iter.into_iter().enumerate() { map.entry(key).or_insert(std::default::default()).push(index) }
	map
}

use std::f64::consts::LN_2;
use super::{*, program::*};

struct T<T> { log_T: T, T: T, T2: T, T3: T, T4: T, rcp_T: T, rcp_T2: T }
macro_rules! T{ {$($field:ident),*} => (T{$($field: r#use($field)),*}) }
impl From<&T<Value>> for T<Expression> { fn from(&T{log_T,T,T2,T3,T4,rcp_T,rcp_T2}:&T<Value>) -> Self { T!{log_T,T,T2,T3,T4,rcp_T,rcp_T2} } }

fn thermodynamic_function(thermodynamics: &[NASA7], expression: impl Fn(&[f64], T<Expression>)->Expression) -> Subroutine {
	let (parameters, [log_T,T,T2,T3,T4,rcp_T]) = parameters(self::stringify![log_T,T,T2,T3,T4,rcp_T]);
	let ref Ts = T{log_T,T,T2,T3,T4,rcp_T,rcp_T2:Value::NONE};
	Subroutine{
		parameters: parameters.to_vec(),
		output: thermodynamics.len(),
		statements: bucket(thermodynamics.iter().map(|s| s.temperature_split.to_bits())).into_iter().map(|(temperature_split, species)| Statement::ConditionalStatement{
			condition: less(r#use(T), f64::from_bits(temperature_split)),
			consequent: species.iter().map(|&specie| output(specie, expression(&thermodynamics[specie].pieces[0], Ts.into()))).collect(),
			alternative: species.iter().map(|&specie| output(specie, expression(&thermodynamics[specie].pieces[1], Ts.into()))).collect()
		}).collect()
	}
}

pub fn molar_heat_capacity_at_constant_pressure_R(thermodynamics: &[NASA7]) -> Subroutine {
		thermodynamic_function(&thermodynamics, |a, T{T,T2,T3,T4,..}| a[0]+a[1]*T+a[2]*T2+a[3]*T3+a[4]*T4 )
}
pub fn enthalpy_RT(thermodynamics: &[NASA7]) -> Subroutine {
	thermodynamic_function(&thermodynamics, |a, T{T,T2,T3,T4,rcp_T,..}| a[0]+a[1]/2.*T+a[2]/3.*T2+a[3]/4.*T3+a[4]/5.*T4+a[5]*rcp_T )
}
pub fn exp_Gibbs_RT(thermodynamics: &[NASA7]) -> Subroutine {
	thermodynamic_function(&thermodynamics, |a, T{log_T,T,T2,T3,T4,rcp_T,..}|
			exp2((a[0]-a[6])/LN_2-a[0]*log_T-a[1]/2./LN_2*T+(1./3.-1./2.)*a[2]/LN_2*T2+(1./4.-1./3.)*a[3]/LN_2*T3+(1./5.-1./4.)*a[4]/LN_2*T4+a[5]/LN_2*rcp_T) )
}

// A.T^β.exp(-Ea/kT)
fn arrhenius(&RateConstant{preexponential_factor: A, temperature_exponent, activation_temperature}: &RateConstant,
											T{log_T,T,T2,T4,rcp_T,rcp_T2,..}: T<Expression>) -> Expression {
	if [0.,-1.,1.,2.,4.,-2.].contains(&temperature_exponent) && activation_temperature == 0. {
		if temperature_exponent == 0. { A.into() }
		else if temperature_exponent == -1. { A * rcp_T }
		else if temperature_exponent == 1. { A * T }
		else if temperature_exponent == 2. { A * T2 }
		else if temperature_exponent == 4. { A * T4 }
		else if temperature_exponent == -2. { A * rcp_T2 }
		else { unreachable!() }
	} else {
		let βlogT𐊛logA = if temperature_exponent == 0. { f64::log2(A).into() } else { temperature_exponent * log_T + f64::log2(A) };
		let log_arrhenius = if activation_temperature == 0. { βlogT𐊛logA } else { -activation_temperature/LN_2 * rcp_T + βlogT𐊛logA };
		exp2(log_arrhenius)
	}
}

fn rcp_arrhenius(&RateConstant{preexponential_factor: A, temperature_exponent, activation_temperature}: &RateConstant,
														T{log_T,T,T2,rcp_T,rcp_T2,..}: T<Expression>) -> Expression {
	if [0.,-1.,1.,2.,4.,-2.].contains(&temperature_exponent) && activation_temperature == 0. {
		if temperature_exponent == 0. { (1./A).into() }
		else if temperature_exponent == -1. { T / A }
		else if temperature_exponent == 1. { rcp_T / A }
		else if temperature_exponent == 2. { rcp_T2 / A }
		else if temperature_exponent == -2. { T2 / A }
		else { unreachable!() }
	} else {
		let m_βlogT𐊛logA = if temperature_exponent == 0. { (-f64::log2(A)).into() } else { - temperature_exponent * log_T - f64::log2(A) };
		let m_log_arrhenius = if activation_temperature == 0. { m_βlogT𐊛logA } else { activation_temperature/LN_2 * rcp_T + m_βlogT𐊛logA };
		exp2(m_log_arrhenius)
	}
}

impl ReactionModel {
fn efficiency(&self, k_inf: &RateConstant, T: &T<Value>, concentrations: Value, f: &Block) -> Expression {
	use ReactionModel::*; match self {
		Elementary|Irreversible => arrhenius(k_inf, T.into()),
		ThreeBody{efficiencies} => arrhenius(k_inf, T.into()) * dot(efficiencies, concentrations),
		PressureModification{efficiencies, k0} => f.block(|def|{
			let ref Pr = def(dot(efficiencies, concentrations) * arrhenius(k0, T.into()));
			Pr / (rcp_arrhenius(k_inf, T.into()) * Pr + 1.)
		}),
		Falloff{efficiencies, k0, troe} => f.block(|def|{
			let ref Pr = def(dot(efficiencies, concentrations) * arrhenius(k0, T.into()));
			let model::Troe{A, T3, T1, T2} = *troe;
			let Fcent = {let &T{T,rcp_T,..}=T; (1.-A) * exp2(r#use(T)/(-LN_2*T3)) + A * exp2(r#use(T)/(-LN_2*T1)) + exp2((-T2/LN_2)*r#use(rcp_T))};
			let ref logFcent = def(log2(Fcent));
			let c =-0.67*logFcent - 0.4*f64::log2(10.);
			let N = -1.27*logFcent + 0.75*f64::log2(10.);
			let ref logPr𐊛c = def(log2(Pr) + c);
			let ref f1 = def(logPr𐊛c / (-0.14*logPr𐊛c+N));
			let F = exp2(logFcent/(f1*f1+1.));
			Pr / (rcp_arrhenius(k_inf, T.into()) * Pr + 1.) * F
		})
	}
}
}

pub fn rates(reactions: &[Reaction]) -> Subroutine {
	let species_len = reactions[0].reactants.len();
	let (parameters,                   [log_T,T,T2,T4,rcp_T,rcp_T2,P0_RT,rcp_P0_RT, exp_Gibbs0_RT,concentrations])
	= parameters(self::stringify![log_T,T,T2,T4,rcp_T,rcp_T2,P0_RT,rcp_P0_RT, exp_Gibbs0_RT,concentrations]);
	let ref T = T{log_T,rcp_T,T,T2,T3:Value::NONE,T4,rcp_T2};
	let mut f = Block::new(&parameters);
	let rates = iter::map(reactions.iter().enumerate(), |(_i, Reaction{reactants, products, net, Σnet, rate_constant, model, ..})| {
		let efficiency = model.efficiency(rate_constant, T, concentrations, &f); // todo: CSE
		let forward_rate = product_of_exponentiations(reactants, concentrations);
		let rate = if let ReactionModel::Irreversible = model { forward_rate } else {
			let rcp_equilibrium_constant = product_of_exponentiations(net, exp_Gibbs0_RT)
																											* match -Σnet { 0=>(1.).into(), 1=>r#use(P0_RT), -1=>r#use(rcp_P0_RT), _ => unreachable!()};
			let reverse_rate = rcp_equilibrium_constant * product_of_exponentiations(products, concentrations);
			forward_rate - reverse_rate
		};
		let rate = f.def(efficiency * rate);
		//f.statements.push(output(species_len-1+_i, r#use(rate)));
		rate
	});
	f.extend((0..species_len-1).map(|specie|
		output(specie, rates.iter().zip(reactions.iter()).fold(None, |dtω, (&rate, Reaction{net, ..})| {
				let rate = r#use(rate);
				match net[specie] {
					0 => dtω,
					1 => match dtω { None => Some(rate), Some(dtω) => Some(dtω + rate) }
					-1 => match dtω { None => Some(-rate), Some(dtω) => Some(dtω - rate) }
					ν => match dtω { None => Some((ν as f64) * rate), Some(dtω) => Some((ν as f64) * rate + dtω) }
				}
		}).unwrap())
	));
	Subroutine {parameters: parameters.into(), output: species_len-1/*+rates.len()*/, statements: f.statements}
}
