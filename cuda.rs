#![feature(format_args_capture,in_band_lifetimes,default_free_fn,associated_type_bounds,unboxed_closures,fn_traits,trait_alias,iter_zip,bool_to_option)]
#![allow(non_snake_case,non_upper_case_globals)]

use {std::iter::zip, ast::*, iter::{list, map}, itertools::Itertools};

type Value = String;

struct Builder<'t> {
	builder: Vec<String>,
	values: Box<[Option<(ast::Type, Value)>]>,
	#[allow(dead_code)] names: &'t [String],
}
impl std::ops::Deref for Builder<'_> { type Target=Vec<String>; fn deref(&self) -> &Self::Target { &self.builder } }
impl std::ops::DerefMut for Builder<'_> { fn deref_mut(&mut self) -> &mut Self::Target { &mut self.builder } }

impl Builder<'_> {
fn rtype(&self, e: &Expression) -> ast::Type { e.rtype(&|value| self.values[value.0].as_ref().unwrap().0) }
fn expr(&mut self, expr: &Expression) -> Value {
	match expr {
		Expression::Expr(e) => {
			use Expr::*;
			match e {
				&F32(value) => value.to_string(),
				&F64(value) => if *value==0. { "0.".to_string() } else { value.to_string() },
				Value(v) => format!("v{}", v.0),
				Neg(x) => { let x = self.expr(x); format!("-{x}") }
				Max(a, b) => { let [a,b] = [a,b].map(|x| self.expr(x)); format!("max({a},{b})") }
				Add(a, b) => { let [a,b] = [a,b].map(|x| self.expr(x)); format!("({a}+{b})") }
				Sub(a, b) => { let [a,b] = [a,b].map(|x| self.expr(x)); format!("({a}-{b})") }
				LessOrEqual(a, b) => { let [a,b] = [a,b].map(|x| self.expr(x)); format!("({a}<{b})") }
				Mul(a, b) => { let [a,b] = [a,b].map(|x| self.expr(x)); format!("({a}*{b})") }
				Div(a, b) => { let [a,b] = [a,b].map(|x| self.expr(x)); format!("({a}/{b})") }
				Sqrt(x) => { let x = self.expr(x); format!("sqrt({x})") }
				Exp(x) => { let x = self.expr(x); format!("exp({x})") }
				Ln{x,..} => { let x = self.expr(x); format!("log({x})") }
			}
		}
		Expression::Block { statements, result } => {
			for s in &**statements { self.extend(s) }
			self.expr(result)
		}
	}
}
fn extend(&mut self, s: &Statement) {
	use Statement::*;
	match s {
		Value { id, value } => {
			let result = self.expr(value);
			let rtype = self.rtype(value);
			self.push(format!("{rtype} v{} = {result};", id.0));
			assert!(self.values[id.0].replace((rtype,result)).is_none());
		},
		Select { condition, true_exprs, false_exprs, results } => {
			let types = map(&**true_exprs, |e| self.rtype(e));
			for (t,e) in zip(&*types, &**false_exprs) { assert!(self.rtype(e) == *t); }
			for (rtype, id) in zip(&*types, &**results) { self.push(format!("{rtype} v{};", id.0)); }

			let condition = self.expr(condition);
			let scope = self.values.clone();
			let true_values = map(&**true_exprs, |e| self.expr(e));
			self.values = scope;
			assert!(results.len() == true_values.len());
			let scope = self.values.clone();
			let false_values = map(&**false_exprs, |e| self.expr(e));
			self.values = scope;
			assert!(results.len() == false_values.len());

			let true_values = zip(&**results, &*true_values).map(|(id, value)| format!("{} = {value};",id.0)).format("\n\t");
			let false_values = zip(&**results, &*false_values).map(|(id, value)| format!("{} = {value};",id.0)).format("\n\t");

			self.push(format!("if {condition} {{\n\t{true_values}\n}} else {{\n\t{false_values}\n}}"));
			for (rtype, id) in zip(&*types, &**results) {
				assert!(self.values[id.0].replace((*rtype, format!("phi"))).is_none());
			}
		}
	}
}
}

pub fn compile(constants_len: usize, ast: ast::Function, name: &str) -> String {
	let output_types = {
		let mut types = Types(ast.input.iter().copied().map(Some).chain((ast.input.len()..ast.values.len()).map(|_| None)).collect());
		for s in &*ast.statements { types.push(s); }
		map(&*ast.output, |output| {
			types.expr(output);
			types.rtype(output)
		})
	};
	let parameters =
		ast.input[0..constants_len].iter().enumerate().map(|(i, atype)| format!("const {atype} v{i}"))
		.chain( ast.input.iter().enumerate().skip(constants_len).map(|(i, atype)| format!("const {atype} in{i}[]")) )
		.chain( output_types.iter().enumerate().map(|(i, rtype)| format!("{rtype} out{i}[]")) ).format(", ");
	let input_values = ast.input.iter().enumerate().skip(constants_len).map(|(i, atype)| format!("const {atype} v{i} = in{i}[id];")).format("\n");
	let mut b = Builder{
		builder: vec![],
		values: list(ast.input.iter().enumerate().map(|(i, atype)| Some((*atype, format!("v{i}")))).chain((ast.input.len()..ast.values.len()).map(|_| None))),
		names: &ast.values
	};
	for s in &*ast.statements { b.extend(s); }
	let instructions = b.builder.iter().format("\n").to_string();
	let store = ast.output.iter().enumerate().map(|(i, expr)| {
		let value = b.expr(expr);
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

	let thermal_conductivityIVT = {
		let_!{ input@[ref lnT, ref lnT2, ref lnT3, ref mole_proportions @ ..] = &*map(0..(3+K), Value) => {
		let mut values = ["lnT","lnT2","lnT3"].iter().map(|s| s.to_string()).chain((0..K).map(|i| format!("X{i}"))).collect::<Vec<_>>();
		assert!(input.len() == values.len());
		let mut function = Block::new(&mut values);
		Function{
			output: list([transport::thermal_conductivityIVT(&thermal_conductivityIVT, &[(1.).into(), lnT.into(), lnT2.into(), lnT3.into()], mole_proportions, &mut function)]),
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
		let_!{ input@[_, ref mean_molar_mass, ref VT, ref lnT, ref lnT2, ref lnT3, ref mole_proportions @ ..] = &*map(0..(6+K), Value) => {
		let mut values = ["id","mean_molar_mass","VT","lnT","lnT2","lnT3"].iter().map(|s| s.to_string()).chain((0..K).map(|i| format!("X{i}"))).collect::<Vec<_>>();
		assert!(input.len() == values.len());
		let mut function = Block::new(&mut values);
		Function{
			output: list(transport::density_diffusion(molar_mass, &binary_thermal_diffusionITVT, mean_molar_mass, VT, &[(1.).into(), lnT.into(), lnT2.into(), lnT3.into()], mole_proportions, &mut function)),
			statements: function.statements.into(),
			input: vec![Type::F64; input.len()].into(),
			values: values.into()
		}
	}}};

	let compile = |f:Function,name| {
		let input = f.input.len();
		let output = f.output.len();
		let mut s = self::compile(0, f, name);
		s = s.replace("__global__","__device__").replace("const unsigned int id = blockIdx.x * blockDim.x + threadIdx.x;","");
		for i in 0..input { let i = format!("in{}", i); s = s.replace(&format!("{i}[]"),&i).replace(&format!("{i}[id]"), &i); }
		if output == 1 { s = s.replace(", double out0[]","").replace("out0[id] = ","return ").replace("__device__ void","__device__ double"); }
		//else { for i in 0..output { let i = format!("out{}", i); s = s.replace(&format!("{i}[]"),&format!("&{i}")).replace(&format!("{i}[id]"), &i); } }
		eprintln!("{}", s.lines().map(|l| l.len()).max().unwrap());
		s
	};
	println!("{}", compile(thermal_conductivityIVT,"thermal_conductivityIVT"));
	println!("{}", compile(viscosityIVT,"viscosityIVT"));
	println!("{}", compile(density_diffusion,"density_diffusion"));
}
