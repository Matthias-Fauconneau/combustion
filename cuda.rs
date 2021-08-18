#![feature(format_args_capture,in_band_lifetimes,default_free_fn,associated_type_bounds,unboxed_closures,fn_traits,trait_alias,iter_zip,bool_to_option)]
#![allow(non_snake_case,non_upper_case_globals)]

use {std::iter::zip, ast::*, iter::{list, map}, itertools::Itertools};

type Value = String;

struct Builder<'t> {
	builder: Vec<String>,
	values: Box<[Option<(ast::Type, Value)>]>,
	names: &'t [String],
}
impl std::ops::Deref for Builder<'_> { type Target=Vec<String>; fn deref(&self) -> &Self::Target { &self.builder } }
impl std::ops::DerefMut for Builder<'_> { fn deref_mut(&mut self) -> &mut Self::Target { &mut self.builder } }

impl Builder<'_> {
fn rtype(&self, e: &Expression) -> ast::Type { e.rtype(&|value| self.values[value.0].as_ref().unwrap().0) }
fn expr(&mut self, expr: &Expression, parent: Option<&Expr>) -> Value {
	match expr {
		Expression::Expr(e) => {
			use Expr::*;
			match e {
				&F32(value) => value.to_string(),
				&F64(value) => if *value==0. { "0.".to_string() } else if *value==1. { "1.".to_string() } else { format!("{:<26}",*value) },
				Value(v) => self.names[v.0].clone(),
				Neg(x) => { let x = self.expr(x, Some(&e)); format!("-{x}") }
				Max(a, b) => { let [a,b] = [a,b].map(|x| self.expr(x, Some(&e))); format!("max({a},{b})") }
				Add(a, b) => { let [a,b] = [a,b].map(|x| self.expr(x, Some(&e))); let y = format!("{a}+{b}"); let y = if y.len() > 139 { format!("{a}\n\t+{b}") } else { y };
					if let None|Some(Add(_,_)) = parent { y } else { format!("({y})") } }
				Sub(a, b) => { let [a,b] = [a,b].map(|x| self.expr(x, Some(&e))); let y = format!("{a}-{b}");
					if let None|Some(Add(_,_)) = parent { y } else { format!("({y})") } }
				LessOrEqual(a, b) => { let [a,b] = [a,b].map(|x| self.expr(x, Some(&e))); assert!(parent.is_none()); format!("{a}<{b}") }
				Mul(a, b) => { let [a,b] = [a,b].map(|x| self.expr(x, Some(&e))); let y = format!("{a}*{b}");
					if let Some(Div(_,_)) = parent { format!("({y})") } else { y } }
				Div(a, b) => { let [a,b] = [a,b].map(|x| self.expr(x, Some(&e))); let y = format!("{a}/{b}");
					if let Some(Div(_,_)) = parent { format!("({y})") } else { y } }
				Sqrt(x) => { let x = self.expr(x, None); format!("sqrt({x})") }
				Exp(x) => { let x = self.expr(x, None); format!("exp({x})") }
				Ln{x,..} => { let x = self.expr(x, None); format!("log({x})") }
				Sq(x) => { let x = self.expr(x, None); format!("sq({x})") }
			}
		}
		Expression::Block { statements, result } => {
			for s in &**statements { self.extend(s) }
			self.expr(result, None)
		}
	}
}
fn extend(&mut self, s: &Statement) {
	use Statement::*;
	match s {
		Value { id, value } => {
			let result = self.expr(value, None);
			let rtype = self.rtype(value);
			self.builder.push(format!("{rtype} {} = {result};", self.names[id.0]));
			assert!(self.values[id.0].replace((rtype,self.names[id.0].clone())).is_none());
		},
		Select { condition, true_exprs, false_exprs, results } => {
			let types = map(&**true_exprs, |e| self.rtype(e));
			for (t,e) in zip(&*types, &**false_exprs) { assert!(self.rtype(e) == *t); }
			for (rtype, id) in zip(&*types, &**results) { self.builder.push(format!("{rtype} {};", self.names[id.0])); }

			let condition = self.expr(condition, None);
			let scope = self.values.clone();
			let true_values = map(&**true_exprs, |e| self.expr(e, None));
			self.values = scope;
			assert!(results.len() == true_values.len());
			let scope = self.values.clone();
			let false_values = map(&**false_exprs, |e| self.expr(e, None));
			self.values = scope;
			assert!(results.len() == false_values.len());

			let true_values = zip(&**results, &*true_values).map(|(id, value)| format!("{} = {value};",self.names[id.0])).format("\n\t");
			let false_values = zip(&**results, &*false_values).map(|(id, value)| format!("{} = {value};",self.names[id.0])).format("\n\t");

			self.push(format!("if {condition} {{\n\t{true_values}\n}} else {{\n\t{false_values}\n}}"));
			for (rtype, id) in zip(&*types, &**results) {
				assert!(self.values[id.0].replace((*rtype, format!("phi"))).is_none());
			}
		}
	}
}
}

pub fn compile(constants_len: usize, ast: &ast::Function, name: &str) -> String {
	let output_types = {
		let mut types = Types(ast.input.iter().copied().map(Some).chain((ast.input.len()..ast.values.len()).map(|_| None)).collect());
		for s in &*ast.statements { types.push(s); }
		map(&*ast.output, |output| {
			types.expr(output);
			types.rtype(output)
		})
	};
	let parameters =
		ast.input[0..constants_len].iter().enumerate().map(|(i, atype)| format!(/*const*/"{atype} {}",ast.values[i]))
		.chain( ast.input.iter().enumerate().skip(constants_len).map(|(i, atype)| format!(/*const*/"{atype} in{i}[]")) )
		.chain( output_types.iter().enumerate().map(|(i, rtype)| format!("{rtype} out{i}[]")) ).format(", ");
	let input_values = ast.input.iter().enumerate().skip(constants_len).map(|(i, atype)| format!(/*const*/"{atype} {} = in{i}[id];",ast.values[i])).format("\n");
	let mut b = Builder{
		builder: vec![],
		values: list(ast.input.iter().enumerate().map(|(i, atype)| Some((*atype, ast.values[i].clone()))).chain((ast.input.len()..ast.values.len()).map(|_| None))),
		names: &ast.values
	};
	for s in &*ast.statements { b.extend(s); }
	let instructions = b.builder.iter().format("\n").to_string();
	let store = ast.output.iter().enumerate().map(|(i, expr)| {
		let value = b.expr(expr, None);
		format!("out{i}[id] = {value};")
	}).format("\n");
	format!(r#"__global__ void {name}({parameters}) {{
const unsigned int id = blockIdx.x * blockDim.x + threadIdx.x;
{input_values}
{instructions}
{store}
}}"#).replace("f64","double")
}

mod yaml;
mod device;
use {anyhow::{Error, Context}, iter::Dot, std::env::*, combustion::*};

#[fehler::throws] fn main() {
	let path = args().skip(1).next().unwrap_or("gri30".to_string());
	let path = if std::path::Path::new(&path).exists() { path } else { format!("/usr/share/cantera/data/{path}.yaml") };
	let model = yaml::Loader::load_from_str(std::str::from_utf8(&std::fs::read(&path).context(path)?)?)?;
	let model = yaml::parse(&model);
	let (ref _species_names, ref species@Species{ref molar_mass, ref thermodynamics, ..}, _, _, State{ref amounts, temperature, pressure_R, ..}) = new(&model);
	let K = species.len();
	let total_amount = amounts.iter().sum::<f64>();
	let mole_fractions = map(&**amounts, |n| n/total_amount);
	let diffusivity = 1.;
	let concentration = pressure_R / temperature;
	let mean_molar_mass = zip(&**molar_mass, &*mole_fractions).map(|(m,x)| m*x).sum::<f64>();
	let density = concentration * mean_molar_mass;
	let Vviscosity = f64::sqrt(density * diffusivity);
	let mean_molar_heat_capacity_at_CP_R:f64 = thermodynamics.iter().map(|a| a.molar_heat_capacity_at_constant_pressure_R(temperature)).dot(mole_fractions);
	let R = kB*NA;
	let thermal_conductivity = mean_molar_heat_capacity_at_CP_R * R / mean_molar_mass * diffusivity;
	let transport::Polynomials{thermal_conductivityIVT, VviscosityIVVT, binary_thermal_diffusionITVT} = transport::Polynomials::<4>::new(&species, temperature);
	let VviscosityIVVT = map(&*VviscosityIVVT, |P| P.map(|p| (f64::sqrt(f64::sqrt(temperature))/Vviscosity)*p));
	let thermal_conductivityIVT = map(&*thermal_conductivityIVT, |P| P.map(|p| (f64::sqrt(temperature)/(2.*thermal_conductivity))*p));
	let binary_thermal_diffusionITVT = map(&*binary_thermal_diffusionITVT, |P| P.map(|p| (f64::sqrt(temperature)/(R*density*diffusivity))*p));

	let thermal_conductivityNIVT = {
		let_!{ input@[ref lnT, ref lnT2, ref lnT3, ref mole_proportions @ ..] = &*map(0..(3+K), Value) => {
		let mut values = ["lnT","lnT2","lnT3"].iter().map(|s| s.to_string()).chain((0..K).map(|i| format!("X{i}"))).collect::<Vec<_>>();
		assert!(input.len() == values.len());
		let mut function = Block::new(&mut values);
		Function{
			output: list([transport::thermal_conductivityNIVT(&thermal_conductivityIVT, &[(1.).into(), lnT.into(), lnT2.into(), lnT3.into()], mole_proportions, &mut function)]),
			statements: function.statements.into(),
			input: vec![Type::F64; input.len()].into(),
			values: values.into()
		}
	}}};

	let viscosityIVT = {
		let_!{ input@[ref lnT, ref lnT2, ref lnT3, ref mole_proportions @ ..] = &*map(0..(3+K), Value) => {
		let mut values = ["lnT","lnT2","lnT3"].iter().map(|s| s.to_string()).chain((0..K).map(|i| format!("X{i}"))).collect::<Vec<_>>();
		assert!(input.len() == values.len());
		let mut function = Block::new(&mut values);
		Function{
			output: list([transport::viscosityIVT(molar_mass, &VviscosityIVVT, &[(1.).into(), lnT.into(), lnT2.into(), lnT3.into()], mole_proportions, &mut function)]),
			statements: function.statements.into(),
			input: vec![Type::F64; input.len()].into(),
			values: values.into()
		}
	}}};

	let density_diffusion = {
		let_!{ input@[ref mean_molar_mass_VT, ref VT, ref lnT, ref lnT2, ref lnT3, ref mole_proportions @ ..] = &*map(0..(5+K), Value) => {
		let mut values = ["mean_molar_mass_VT","VT","lnT","lnT2","lnT3"].iter().map(|s| s.to_string()).chain((0..K).map(|i| format!("X{i}"))).collect::<Vec<_>>();
		assert_eq!(input.len(), values.len());
		let mut function = Block::new(&mut values);
		Function{
			output: list(transport::density_diffusion(molar_mass, &binary_thermal_diffusionITVT, mean_molar_mass_VT, VT, &[(1.).into(), lnT.into(), lnT2.into(), lnT3.into()], mole_proportions, &mut function)),
			statements: function.statements.into(),
			input: vec![Type::F64; input.len()].into(),
			values: values.into()
		}
	}}};

	let compile = |f:&Function,name| {
		let mut s = self::compile(0, f, name);
		s = s.replace("__global__","__device__").replace("const unsigned int id = blockIdx.x * blockDim.x + threadIdx.x;\n","");
		for i in 0..f.input.len() { s = s.replace(&format!("in{i}[]"),&f.values[i]).replace(&format!(/*const*/"double {} = in{i}[id];\n",&f.values[i]),""); }
		if f.output.len() == 1 { s = s.replace(", double out0[]","").replace("out0[id] = ","return ").replace("__device__ void","__device__ double"); }
		//else { for i in 0..output { let i = format!("out{}", i); s = s.replace(&format!("{i}[]"),&format!("&{i}")).replace(&format!("{i}[id]"), &i); } }
		eprintln!("{}", s.lines().map(|l| l.len()).max().unwrap());
		s//.replace("+-","-")
	};
	println!("{}", compile(&thermal_conductivityNIVT,"thermal_conductivityNIVT"));
	println!("{}", compile(&viscosityIVT,"viscosityIVT").replace("rcp_VviscosityIVVT","r").replace("VviscosityIVVT","v"));
	println!("{}", compile(&density_diffusion,"density_diffusion").replace("density_diffusion(","density_diffusion(int id, "));

	/*let transport = transport::properties::<4>(&species, temperature, Vviscosity, thermal_conductivity, density*diffusivity);
	use device::*;
	let transport = with_repetitive_input(assemble::<f64>(transport, 1), 1);
	let thermal_conductivityNIVT = with_repetitive_input(assemble::<f64>(thermal_conductivityNIVT, 1), 1);
	let viscosityIVT = with_repetitive_input(assemble::<f64>(viscosityIVT, 1), 1);
	let density_diffusion = with_repetitive_input(assemble::<f64>(density_diffusion, 1), 1);
	let temperature0 = temperature;
	let total_amount = amounts.iter().sum::<f64>();
	let nonbulk_amounts = map(&amounts[0..amounts.len()-1], |&n| n);
	let_!{ [thermal_conductivity, viscosity, ref_density_diffusion @ ..] = &*transport(&[], &([&[total_amount, temperature/temperature0], &*nonbulk_amounts].concat())).unwrap() => {
		let mole_fractions = map(&**amounts, |n| n/total_amount);
		let mean_molar_mass : f64 = molar_mass.dot(&mole_fractions);
		let mass_fractions = map(zip(&*mole_fractions,&**molar_mass), |(x,&m)| m / mean_molar_mass * x);
		let T = temperature/temperature0;
		let VT =	f64::sqrt(T);
		let mean_molar_mass_VTN = VT; // M * (VT * sum mole proportions) = M * VT / M
		let lnT = f64::ln(T);
		let lnT2 = lnT*lnT;
		let lnT3 = lnT2*lnT;
		let mut mole_proportions = vec![0.; K];
		let mut sum_nonbulk_mass_fractions = 0.;
		for k in 0..K-1 {
			let mass_fraction = f64::max(0., mass_fractions[k]);
			mole_proportions[k] = mass_fraction/molar_mass[k];
			sum_nonbulk_mass_fractions += mass_fraction;
		}
		mole_proportions[K-1] = (1. - sum_nonbulk_mass_fractions) / molar_mass[K-1];
		let_!{ [thermal_conductivityNIVT] = &*thermal_conductivityNIVT(&[], &([&[lnT, lnT2, lnT3], &*mole_proportions].concat())).unwrap() => {
		{let e = num::relative_error(*thermal_conductivity, thermal_conductivityNIVT*mean_molar_mass*VT); assert!(e<4e-3, "{e:e}");}
		let_!{ [viscosityIVT] = &*viscosityIVT(&[], &([&[lnT, lnT2, lnT3], &*mole_proportions].concat())).unwrap() => {
		{let e = num::relative_error(*viscosity, viscosityIVT*VT); assert!(e<=0., "{e:e}");}
		let density_diffusion = density_diffusion(&[], &([&[mean_molar_mass_VTN, VT, lnT, lnT2, lnT3], &*mole_proportions].concat())).unwrap();
		for (&a,&b) in zip(ref_density_diffusion, &*density_diffusion) {let e = num::relative_error(a,b); assert!(e<1e-15, "{e:e}");}
	}}}}}}*/
}
