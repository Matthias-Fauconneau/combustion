#![feature(format_args_capture,in_band_lifetimes,default_free_fn,associated_type_bounds,unboxed_closures,fn_traits,trait_alias,iter_zip,bool_to_option,type_ascription)]
#![allow(non_snake_case,non_upper_case_globals)]

use {std::iter::zip, ast::*, iter::{list, map}, itertools::Itertools};

struct Builder<'t> {
	builder: Vec<String>,
	values: Box<[Option<(ast::Type, String)>]>,
	registers: &'t [Option<u8>],
	names: &'t [String],
}
impl std::ops::Deref for Builder<'_> { type Target=Vec<String>; fn deref(&self) -> &Self::Target { &self.builder } }
impl std::ops::DerefMut for Builder<'_> { fn deref_mut(&mut self) -> &mut Self::Target { &mut self.builder } }

impl Builder<'_> {
fn rtype(&self, e: &Expression) -> ast::Type { e.rtype(&|value| self.values[value.0].as_ref().unwrap().0) }
fn expr(&mut self, expr: &Expression, parent: Option<&Expr>) -> String {
	match expr {
		Expression::Expr(e) => {
			use Expr::*;
			match e {
				&F32(value) => value.to_string(),
				&F64(value) => if *value==0. { "0.".to_string() } else if *value==1. { "1.".to_string() } else { format!("{}",*value) },
				Value(v) => if let Some(r) = self.registers[v.0] { format!("r[{}]/*{}*/", r, self.names[v.0]) } else { self.names[v.0].clone() },
				Neg(x) => { let x = self.expr(x, Some(&e)); format!("-{x}") }
				Max(a, b) => { let [a,b] = [a,b].map(|x| self.expr(x, Some(&e))); format!("max({a},{b})") }
				Add(a, b) => { let [a,b] = [a,b].map(|x| self.expr(x, Some(&e))); let y = format!("{a}+{b}"); let y = if y.len() > 139 { format!("{a}\n +{b}") } else { y };
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
			self.builder.push(format!("r[{}]/*{}*/ = {result};", self.registers[id.0].unwrap(), self.names[id.0]));
			assert!(self.values[id.0].replace((rtype,format!("r[{}]",self.registers[id.0].unwrap()))).is_none());
		},
		Select { condition, true_exprs, false_exprs, results } => {
			let types = map(&**true_exprs, |e| self.rtype(e));
			for (t,e) in zip(&*types, &**false_exprs) { assert!(self.rtype(e) == *t); }
			//for (rtype, id) in zip(&*types, &**results) { self.builder.push(format!("{rtype} r{};", self.registers[id.0])); }

			let condition = self.expr(condition, None);
			let scope = self.values.clone();
			let true_values = map(&**true_exprs, |e| self.expr(e, None));
			self.values = scope;
			assert!(results.len() == true_values.len());
			let scope = self.values.clone();
			let false_values = map(&**false_exprs, |e| self.expr(e, None));
			self.values = scope;
			assert!(results.len() == false_values.len());

			let true_values = zip(&**results, &*true_values).map(|(id, value)| format!("r[{}]/*{}*/ = {value};", self.registers[id.0].unwrap(), self.names[id.0])).format("\n\t");
			let false_values = zip(&**results, &*false_values).map(|(id, value)| format!("r[{}]/*{}*/ = {value};", self.registers[id.0].unwrap(), self.names[id.0])).format("\n\t");

			self.push(format!("if {condition} {{\n\t{true_values}\n}} else {{\n\t{false_values}\n}}"));
			for (rtype, id) in zip(&*types, &**results) {
				assert!(self.values[id.0].replace((*rtype, format!("r[{}]",self.registers[id.0].unwrap()))).is_none());
			}
		}
	}
}
}

pub fn compile(constants_len: usize, ast: &ast::Function, registers: &[Option<u8>], name: impl std::fmt::Display) -> String {
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
	let register_count = registers.iter().filter_map(|&x| x).max().unwrap()+1;
	let mut b = Builder{
		builder: vec![],
		values: list(ast.input.iter().enumerate().map(|(i, atype)| Some((*atype, ast.values[i].clone()))).chain((ast.input.len()..ast.values.len()).map(|_| None))),
		registers,
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
double r[{register_count}];
{instructions}
{store}
}}"#).replace("f64","double")
}

mod yaml;
#[cfg(feature = "check")] mod device;
use {iter::Dot, num::sqrt, std::{path::Path,fs::read,env::args}, self::yaml::{Loader,parse}, combustion::*};

#[fehler::throws(anyhow::Error)] fn main() {
	let path = args().skip(1).next().unwrap_or("gri30".to_string());
	let path = if Path::new(&path).exists() { path } else { format!("/usr/share/cantera/data/{path}.yaml") };
	let model = Loader::load_from_str(std::str::from_utf8(&read(&path)?)?)?;
	let model = parse(&model);
	let (ref _species_names, ref species@Species{ref molar_mass, ref thermodynamics, ..}, _, _, State{ref amounts, temperature, pressure_R, ..}) = new(&model);
	let K = species.len();
	let total_amount = amounts.iter().sum::<f64>();
	let mole_fractions = map(&**amounts, |n| n/total_amount);
	let diffusivity = 1.;
	let concentration = pressure_R / temperature;
	let mean_molar_mass = zip(&**molar_mass, &*mole_fractions).map(|(m,x)| m*x).sum::<f64>();
	let density = concentration * mean_molar_mass;
	let viscosity = density * diffusivity;
	let specific_heat_capacity = R / mean_molar_mass * thermodynamics.iter().map(|a| a.molar_heat_capacity_at_constant_pressure_R(temperature)).dot(mole_fractions):f64;
	let conductivity = specific_heat_capacity * diffusivity;
	let transport::Polynomials{conductivityIVT, VviscosityIVVT, diffusivityITVT} = transport::Polynomials::<5>::new(&species, temperature);
	let VviscosityIVVT = map(&*VviscosityIVVT, |P| P.map(|p| (sqrt(sqrt(temperature)/viscosity))*p));
	let conductivityIVT = map(&*conductivityIVT, |P| P.map(|p| (sqrt(temperature)/(2.*conductivity))*p));
	let diffusivityITVT = map(&*diffusivityITVT, |P| P.map(|p| (sqrt(temperature)/(R*viscosity))*p));

	let conductivityNIVT = {
		let_!{ input@[ref sum_mole_proportions, ref lnT, ref lnT2, ref lnT3, ref lnT4, ref mole_proportions @ ..] = &*map(0..(5+K), Value) => {
		let mut values = ["sum_mole_proportions", "lnT","lnT2","lnT3","lnT4"].iter().map(|s| s.to_string()).chain((0..K).map(|i| format!("mole_proportions{i}"))).collect::<Vec<_>>();
		assert!(input.len() == values.len());
		let mut function = Block::new(&mut values);
		Function{
			output: list([transport::conductivityNIVT(&conductivityIVT, sum_mole_proportions, &[(1.).into(), lnT.into(), lnT2.into(), lnT3.into(), lnT4.into()], mole_proportions, &mut function)]),
			statements: function.statements.into(),
			input: vec![Type::F64; input.len()].into(),
			values: values.into()
		}
	}}};

	let viscosityIVT = {
		let_!{ input@[ref lnT, ref lnT2, ref lnT3, ref lnT4, ref mole_proportions @ ..] = &*map(0..(4+K), Value) => {
		let mut values = ["lnT","lnT2","lnT3","lnT4"].iter().map(|s| s.to_string()).chain((0..K).map(|i| format!("mole_proportions{i}"))).collect::<Vec<_>>();
		assert!(input.len() == values.len());
		let mut function = Block::new(&mut values);
		Function{
			output: list([transport::viscosityIVT(molar_mass, &VviscosityIVVT, &[(1.).into(), lnT.into(), lnT2.into(), lnT3.into(), lnT4.into()], mole_proportions, &mut function)]),
			statements: function.statements.into(),
			input: vec![Type::F64; input.len()].into(),
			values: values.into()
		}
	}}};

	let density_diffusivity = {
		let_!{ input@[ref mean_molar_mass_VTN, ref VT, ref lnT, ref lnT2, ref lnT3, ref lnT4, ref mole_proportions @ ..] = &*map(0..(6+K), Value) => {
		let mut values = ["mean_molar_mass_VTN","VT","lnT","lnT2","lnT3","lnT4"].iter().map(|s| s.to_string()).chain((0..K).map(|i| format!("mole_proportions{i}"))).collect::<Vec<_>>();
		assert_eq!(input.len(), values.len());
		let mut function = Block::new(&mut values);
		Function{
			output: list(transport::density_diffusivity(molar_mass, &diffusivityITVT, mean_molar_mass_VTN, VT, &[(1.).into(), lnT.into(), lnT2.into(), lnT3.into(), lnT4.into()], mole_proportions, &mut function)),
			statements: function.statements.into(),
			input: vec![Type::F64; input.len()].into(),
			values: values.into()
		}
	}}};

	let compile = |mut f: Function, name| {
		let mut last_use = vec![0; f.values.len()];
		fn visitor(last_use: &mut Vec<usize>, program_counter: usize, e: &Expression) {
			if let Expression::Expr(Expr::Value(id)) = e { last_use[id.0] = program_counter; } else { e.visit(|e| visitor(last_use, program_counter, e)); }
		}
		for (program_counter, statement) in f.statements.iter().enumerate() {
			use Statement::*;
			match statement {
				Value{value,..} => visitor(&mut last_use, program_counter, value),
				_/*Select { condition, true_exprs, false_exprs, results }*/ => unimplemented!() /*{
					condition.visit(|e| if Expr(Value(id))=e { last_use[i] = i;});
					for e in true_exprs.chain(false_exprs) { e.visit(|e| if Expr(Value(id))=e { last_use[id.0] = program_counter;}); }
				}*/
			}
		}
		for (output_counter, value) in f.output.iter().enumerate() { visitor(&mut last_use, f.statements.len()+output_counter, value); }

		/*let ref names = f.values;
		for (program_counter, statement) in f.statements.iter().enumerate() {
			eprintln!("{program_counter}: {:?} {:?}", statement.to_string(names), last_use.iter().enumerate().filter(|&(_,&pc)| pc==program_counter).map(|(id,_)| &names[id]).format(" "));
		}*/
		let mut registers: Vec<Option<ast::Value>> = vec![];
		let mut map = vec![None; f.values.len()];
		for (program_counter, statement) in f.statements.iter_mut().enumerate() {
			for (_r, register) in registers.iter_mut().enumerate() { if let Some(id) = register { if program_counter>=last_use[id.0] {
				//eprintln!("-{_r}: {}", names[id.0]);
				*register= None;
				}
			} }
			use Statement::*;
			match statement {
				Value{id,..} => /*{ value.visit(|e| if let Expression::Expr(Expr::Value(id)) = e {*/ if id.0 >= f.input.len() {
					if !registers.contains(&Some(*id)) {
						let r =
							if let Some((r, register)) = registers.iter_mut().enumerate().filter(|(_,register)| register.is_none()).next() { *register = Some(*id); r }
							else { registers.push(Some(*id)); registers.len()-1 };
						//eprintln!("+{r}: {}", f.values[id.0]);
						map[id.0] = Some(r.try_into().unwrap());
					}
				}//}); }
				_/*Select { condition, true_exprs, false_exprs, results }*/ => unimplemented!()
			}
		}
		let mut s = self::compile(0, &f, &map, name);
		s = s.replace("__global__","__NEKRK_DEVICE__").replace("const unsigned int id = blockIdx.x * blockDim.x + threadIdx.x;\n","");
		for i in 0..f.input.len() { s = s.replace(&format!("in{i}[]"),&f.values[i]).replace(&format!(/*const*/"double {} = in{i}[id];\n",&f.values[i]),""); }
		s = s.replace(", double mole_proportions0",", double mole_proportions[]");
		for k in (0..K).rev() { s = s.replace(&format!(", double mole_proportions{k},"),",").replace(&format!("mole_proportions{k}"), &format!("mole_proportions[{k}]")); }
		if f.output.len() == 1 { s = s.replace(", double out0[]","").replace("out0[id] = ","return ").replace("__ void","__ double"); }
		else {
			s = s.replace(", double out0[]",", double* out[]");
			for k in 0..K { s = s.replace(&format!(", double out{k}[]"),"").replace(&format!("out{k}["), &format!("out[{k}][")); }
		}
		s.replace("+-","-").replace("double","dfloat")
	};
	let target = std::path::PathBuf::from(args().skip(2).next().unwrap_or(std::env::var("HOME")?+"/nekRK/share/mechanisms"));
	let name = std::path::Path::new(&path).file_stem().unwrap();
	let name = name.to_str().unwrap();
	let prefix = args().skip(3).next().unwrap_or("nekrk_".to_string());
	std::fs::write(target.join(format!("{name}/conductivity.c")), compile(conductivityNIVT,format!("{prefix}conductivityNIVT")))?;
	std::fs::write(target.join(format!("{name}/viscosity.c")), [
		"__NEKRK_DEVICE__ double sq(double x) { return x*x; }",
		&compile(viscosityIVT,format!("{prefix}viscosityIVT")).replace("rcp_VviscosityIVVT","r").replace("VviscosityIVVT","v")
	].join("\n"))?;
	std::fs::write(target.join(format!("{name}/diffusivity.c")), compile(density_diffusivity,format!("{prefix}density_diffusivity")).replace("density_diffusivity(","density_diffusivity(unsigned int id, "))?;

	#[cfg(feature="check")] {
		let transport = transport::properties::<4>(&species, temperature, viscosity, conductivity);
		use device::*;
		let transport = with_repetitive_input(assemble::<f64>(transport, 1), 1);
		let conductivityNIVT = with_repetitive_input(assemble::<f64>(conductivityNIVT, 1), 1);
		let viscosityIVT = with_repetitive_input(assemble::<f64>(viscosityIVT, 1), 1);
		let density_diffusivity = with_repetitive_input(assemble::<f64>(density_diffusivity, 1), 1);
		let temperature0 = temperature;
		let total_amount = amounts.iter().sum::<f64>();
		let nonbulk_amounts = map(&amounts[0..amounts.len()-1], |&n| n);
		let_!{ [conductivity, viscosity, ref_density_diffusivity @ ..] = &*transport(&[], &([&[total_amount, temperature/temperature0], &*nonbulk_amounts].concat())).unwrap() => {
		let mole_fractions = map(&**amounts, |n| n/total_amount);
		let mean_molar_mass : f64 = molar_mass.dot(&mole_fractions);
		let mass_fractions = map(zip(&*mole_fractions,&**molar_mass), |(x,&m)| m / mean_molar_mass * x);
		let T = temperature/temperature0;
		let VT =	sqrt(T);
		let mean_molar_mass_VTN = VT; // M * (VT * sum mole proportions) = M * VT / M
		let lnT = num::ln(T);
		let lnT2 = lnT*lnT;
		let lnT3 = lnT2*lnT;
		let mut mole_proportions = vec![0.; K];
		let mut sum_nonbulk_mass_fractions = 0.;
		let mut mean_rcp_molar_mass = 0.;
		for k in 0..K-1 {
			let mass_fraction = f64::max(0., mass_fractions[k]);
			mole_proportions[k] = mass_fraction/molar_mass[k];
			sum_nonbulk_mass_fractions += mass_fraction;
			mean_rcp_molar_mass += mole_proportions[k];
		}
		mole_proportions[K-1] = (1. - sum_nonbulk_mass_fractions) / molar_mass[K-1];
		mean_rcp_molar_mass += mole_proportions[K-1];
		let_!{ [conductivityNIVT] = &*conductivityNIVT(&[], &([&[lnT, lnT2, lnT3], &*mole_proportions].concat())).unwrap() => {
		{let e = num::relative_error(*conductivity, conductivityNIVT*VT/mean_rcp_molar_mass); assert!(e<4e-3, "{e:e}");}
		dbg!(conductivityNIVT, VT, mean_rcp_molar_mass);
		let_!{ [viscosityIVT] = &*viscosityIVT(&[], &([&[lnT, lnT2, lnT3], &*mole_proportions].concat())).unwrap() => {
		{let e = num::relative_error(*viscosity, viscosityIVT*VT); assert!(e<=0., "{e:e}");}
		let density_diffusivity = density_diffusivity(&[], &([&[mean_molar_mass_VTN, VT, lnT, lnT2, lnT3], &*mole_proportions].concat())).unwrap();
		for (&a,&b) in zip(ref_density_diffusivity, &*density_diffusivity) {let e = num::relative_error(a,b); assert!(e<1e-15, "{e:e}");}
	}}}}}}}
}
