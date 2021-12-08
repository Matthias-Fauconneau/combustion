#![feature(iter_zip,array_zip,destructuring_assignment,let_else)]#![allow(non_snake_case,non_upper_case_globals)]
use {std::iter::zip, ast::*, iter::{list, map}, itertools::Itertools};

struct Builder<'t> {
	builder: Vec<String>,
	values: Box<[Option<(ast::Type, String)>]>,
	registers: &'t [Option<u16>],
	names: &'t [String],
}
impl std::ops::Deref for Builder<'_> { type Target=Vec<String>; fn deref(&self) -> &Self::Target { &self.builder } }
impl std::ops::DerefMut for Builder<'_> { fn deref_mut(&mut self) -> &mut Self::Target { &mut self.builder } }

pub const WRAP: usize = 139;

impl Builder<'_> {
fn rtype(&self, e: &Expression) -> ast::Type { e.rtype(&|value| self.values[value.0].as_ref().unwrap().0) }
fn value(&self, &Value(id): &Value) -> String {
	if let Some(r) = self.registers[id] { if !self.names[id].is_empty() { format!("r[{r}]/*{}*/", self.names[id]) } else { format!("r[{r}]") } } else {
		self.names[id].replace(",","comma").replace("/","div").replace("-","minus").replace("*","star").replace("#","sharp").replace("+","plus").replace("(","leftpar").replace(")","rightpar").replace("'","prime").replace("β","beta").replace("ν","nu") //.replace("ω̇","rates_")
	}
}
fn expr(&mut self, expr: &Expression, parent: Option<&Expr>, wrap: usize) -> String {
	match expr {
		Expression::Expr(e) => {
			use Expr::*;
			match e {
				&F32(value) => value.to_string(),
				&F64(value) => if *value==0. { "0.".to_string() } else if *value==1. { "1.".to_string() } else { let f = format!("{}",*value); if f.contains('.') { f } else { f+"." } },
				Value(value) => self.value(value),
				Neg(x) => { let x = self.expr(x, Some(&e), wrap); format!("-{x}") }
				Min(a, b) => { let [a,b] = [a,b].map(|x| self.expr(x, Some(&e), wrap)); format!("min({a},{b})") }
				Max(a, b) => { let [a,b] = [a,b].map(|x| self.expr(x, Some(&e), wrap)); format!("max({a},{b})") }
				Add(A, b) => {
					let [a,b] = [A,b].map(|x| self.expr(x, Some(&e), wrap));
					let y = format!("{a}+{b}");
					let y = if wrap > 0 && y.len() > wrap { let a = self.expr(A, Some(&e), 0); format!("{a}\n +{b}") } else { y };
					if let None|Some(Add(_,_)) = parent { y } else { format!("({y})") }
				}
				Sub(a, b) => { let [a,b] = [a,b].map(|x| self.expr(x, Some(&e), wrap)); let y = format!("{a}-{b}");
					if let None|Some(Add(_,_)) = parent { y } else { format!("({y})") } }
				LessOrEqual(a, b) => { let [a,b] = [a,b].map(|x| self.expr(x, Some(&e), wrap)); assert!(parent.is_none()); format!("{a}<{b}") }
				Mul(a, b) => { let [a,b] = [a,b].map(|x| self.expr(x, Some(&e), wrap)); let y = format!("{a}*{b}");
					if let Some(Div(_,_)) = parent { format!("({y})") } else { y } }
				Div(a, b) => { let [a,b] = [a,b].map(|x| self.expr(x, Some(&e), wrap)); let y = format!("{a}/{b}");
					if let Some(Div(_,_)) = parent { format!("({y})") } else { y } }
				Sqrt(x) => { let x = self.expr(x, None, wrap); format!("sqrt({x})") }
				Exp(x) => { let x = self.expr(x, None, wrap); format!("exp({x})") }
				Ln{x,..} => { let x = self.expr(x, None, wrap); format!("log({x})") }
				Sq(x) => { let x = self.expr(x, None, wrap); format!("sq({x})") }
				_ => unreachable!()
			}
		}
		Expression::Block { statements, result } => {
			for s in &**statements { self.extend(s) }
			self.expr(result, None, wrap)
		}
	}
}
fn extend(&mut self, s: &Statement) {
	use Statement::*;
	match s {
		Value { id, value } => {
			let result = self.expr(value, None, WRAP);
			let rtype = self.rtype(value);
			let value = self.value(id);
			self.builder.push(format!("{} = {result};",
				if self.registers[id.0].is_none() && !self.values.iter().filter_map(|v| v.as_ref()).any(|v| v.1 == value) { format!("{rtype} {value}") } else { value.clone() }
			));
			assert!(self.values[id.0].replace((rtype,value)).is_none());
		},
		Select { condition, true_exprs, false_exprs, results } => {
			let types = map(&**true_exprs, |e| self.rtype(e));
			for (t,e) in zip(&*types, &**false_exprs) { assert!(self.rtype(e) == *t); }
			for (rtype, id) in zip(&*types, &**results) { if self.registers[id.0].is_none() {
				let value = self.value(id);
				if !self.values.iter().filter_map(|v| v.as_ref()).any(|v| v.1 == value) { self.builder.push(format!("{rtype} {value};")) }
			}}

			let condition = self.expr(condition, None, WRAP);
			let scope = self.values.clone();
			let true_values = map(&**true_exprs, |e| self.expr(e, None, WRAP));
			self.values = scope;
			assert!(results.len() == true_values.len());
			let scope = self.values.clone();
			let false_values = map(&**false_exprs, |e| self.expr(e, None, WRAP));
			self.values = scope;
			assert!(results.len() == false_values.len());

			let true_values = zip(&**results, &*true_values).map(|(id, value)| format!("{} = {value};", self.value(id))).format(" ").to_string();
			let false_values = zip(&**results, &*false_values).map(|(id, value)| format!("{} = {value};", self.value(id))).format(" ").to_string();

			self.push(format!("if ({condition}) {{\n\t{true_values}\n}} else {{\n\t{false_values}\n}}"));
			for (rtype, id) in zip(&*types, &**results) {
				assert!(self.values[id.0].replace((*rtype, self.value(id))).is_none());
			}
		}
	}
}
}

pub fn compile(constants_len: usize, ast: &ast::Function, registers: &[Option<u16>], name: impl std::fmt::Display) -> String {
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
	let register_count = registers.iter().filter_map(|&x| x).max().map(|r| r+1).unwrap_or(0);
	let mut b = Builder{
		builder: vec![],
		values: list(ast.input.iter().enumerate().map(|(i, atype)| Some((*atype, ast.values[i].clone()))).chain((ast.input.len()..ast.values.len()).map(|_| None))),
		registers,
		names: &ast.values
	};
	for s in &*ast.statements { b.extend(s); }
	let instructions = b.builder.iter().format("\n").to_string();
	let store = ast.output.iter().enumerate().map(|(i, expr)| {
		let value = b.expr(expr, None, WRAP);
		format!("out{i}[id] = {value};")
	}).format("\n");
	format!(r#"__global__ void {name}({parameters}) {{
const unsigned int id = blockIdx.x * blockDim.x + threadIdx.x;
{input_values}
{}{instructions}
{store}
}}"#, if register_count > 0 { format!("double r[{register_count}];\n") } else { "".into() }).replace("f64","double")
}
mod yaml;
#[cfg(feature = "check")] mod device;
use {anyhow::Context, std::{path::Path,fs::read,env::args}, self::yaml::{Loader,parse}, combustion::*};

#[fehler::throws(anyhow::Error)] fn main() {
	let path = args().skip(1).next().unwrap_or("gri30".to_string());
	let path = if Path::new(&path).exists() { path } else { format!("/usr/local/share/cantera/data/{path}.yaml") };
	let model = Loader::load_from_str(std::str::from_utf8(&read(&path).context(path.clone())?)?)?;
	let model = parse(&model);
	let (ref species_names, ref species@Species{ref molar_mass, ref thermodynamics, ..}, _, reactions, State{ref amounts, temperature, pressure_R, ..}) = new(&model);
	let K = species.len();
	let _ = (molar_mass, amounts, temperature, pressure_R);
	#[cfg(feature="transport")] {
	use {iter::Dot, num::sqrt};
	let diffusivity = 1.;
	let concentration = pressure_R / temperature;
	let density = concentration * mean_molar_mass;
	let viscosity = density * diffusivity;
	let total_amount = amounts.iter().sum::<f64>();
	let mole_fractions = map(&**amounts, |n| n/total_amount);
	let mean_molar_mass = zip(&**molar_mass, &*mole_fractions).map(|(m,x)| m*x).sum::<f64>();
	let molar_heat_capacity_R: f64 = thermodynamics.iter().map(|a| a.molar_heat_capacity_at_constant_pressure_R(temperature)).dot(mole_fractions);
	let specific_heat_capacity = R / mean_molar_mass * molar_heat_capacity_R;
	let conductivity = specific_heat_capacity * diffusivity;
	let transport::Polynomials{conductivityIVT, VviscosityIVVT, diffusivityITVT} = transport::Polynomials::<5>::new(&species, temperature);
	let VviscosityIVVT = map(&*VviscosityIVVT, |P| P.map(|p| (sqrt(sqrt(temperature)/viscosity))*p));
	let conductivityIVT = map(&*conductivityIVT, |P| P.map(|p| (sqrt(temperature)/(2.*conductivity))*p));
	let diffusivityITVT = map(&*diffusivityITVT, |P| P.map(|p| (sqrt(temperature)/(R*viscosity))*p));

	let conductivityNIVT = {
		let_!{ input@[sum_mole_proportions, lnT, lnT2, lnT3, lnT4, ref mole_proportions @ ..] = &*map(0..(5+K), Value) => {
		let values = list(["sum_mole_proportions", "lnT","lnT2","lnT3","lnT4"].iter().map(|s| s.to_string()).chain((0..K).map(|i| format!("mole_proportions{i}"))));
		assert!(input.len() == values.len());
		let lnT = [(1.).into(), lnT.into(), lnT2.into(), lnT3.into(), lnT4.into()];
		let mut values = values.into();
		let mut function = Block::new(&mut values);
		Function{
			output: list([transport::conductivityNIVT(&conductivityIVT, sum_mole_proportions, &lnT, mole_proportions, &mut function)]),
			statements: function.statements.into(),
			input: vec![Type::F64; input.len()].into(),
			values: values.into()
		}
	}}};

	let viscosityIVT = {
		let_!{ input@[lnT, lnT2, lnT3, lnT4, mole_proportions @ ..] = &*map(0..(4+K), Value) => {
		let values = ["lnT","lnT2","lnT3","lnT4"].iter().map(|s| s.to_string()).chain((0..K).map(|i| format!("mole_proportions{i}"))).collect::<Vec<_>>();
		assert!(input.len() == values.len());
		let lnT = [(1.).into(), lnT.into(), lnT2.into(), lnT3.into(), lnT4.into()];
		let mut values = values.into();
		let mut function = Block::new(&mut values);
		Function{
			output: list([transport::viscosityIVT(molar_mass, &VviscosityIVVT, &lnT, mole_proportions, &mut function)]),
			statements: function.statements.into(),
			input: vec![Type::F64; input.len()].into(),
			values: values.into()
		}
	}}};

	let density_diffusivity = {
		let_!{ input@[mean_molar_mass_VTN, VT, lnT, lnT2, lnT3, lnT4, mole_proportions @ ..] = &*map(0..(6+K), Value) => {
		let values = list(["mean_molar_mass_VTN","VT","lnT","lnT2","lnT3","lnT4"].iter().map(|s| s.to_string()).chain((0..K).map(|i| format!("mole_proportions{i}"))));
		assert_eq!(input.len(), values.len());
		let lnT = [(1.).into(), lnT.into(), lnT2.into(), lnT3.into(), lnT4.into()];
		let mut values = values.into();
		let mut function = Block::new(&mut values);
		Function{
			output: list(transport::density_diffusivity(molar_mass, &diffusivityITVT, mean_molar_mass_VTN, VT, &lnT, mole_proportions, &mut function)),
			statements: function.statements.into(),
			input: vec![Type::F64; input.len()].into(),
			values: values.into()
		}
	}}};
	}

	let rates = {
		let input@[pressure_R, lnT, T, T2, T3, T4, rcpT, concentrations @ ..] = &*map(0..(7+K), Value) else { panic!() };
		let values = list(["pressure_R","lnT","T","T2","T3","T4","rcpT"].iter().map(|s| s.to_string()).chain((0..K).map(|i| format!("concentrations{i}"))));
		assert_eq!(input.len(), values.len());
		let mut values = values.into();
		let mut function = Block::new(&mut values);
		let ref mut f = function;
		use reaction::*;
		let T = T{lnT:*lnT,T:*T,T2:*T2,T3:*T3,T4:*T4,rcpT:*rcpT,rcpT2: f.def(rcpT*rcpT,"rcpT2")};
		let concentration = f.def(pressure_R * rcpT, "concentration");
		Function{
			output: map(&*species_rates(thermodynamics, &reactions, T, concentration, concentrations, f, species_names), |x| x.into()),
			statements: function.block.statements.into(),
			input: vec![Type::F64; input.len()].into(),
			values: values.into()
		}
	};

	let compile = |mut f: Function, name, array_input, direct| {
		let [constants, SFU, ALU] = {
			fn add<T: std::ops::Add, const N: usize>(a: [T;N], b: [T;N]) -> [<T as std::ops::Add>::Output; N] { a.zip(b).map(|(a,b)| a+b) }
			fn visitor(e: &Expression) -> [usize; 3] {
				let [constants, SFU, ALU] = e.visit(visitor).into_iter().filter_map(|x| x).reduce(add).unwrap_or([0,0,0]);
				use Expr::*;
				match e { Expression::Expr(e) => match e {
					F32{..}|F64{..} => [constants+1, SFU, ALU],
					Value(_)|Neg(_/*TODO: fmsub*/) => [constants, SFU, ALU],
					Div{..}|Sqrt{..}|Exp{..}|Ln{..} => [constants, SFU+1, ALU], //TODO: 1/sqrt
					Add{..}|Sub{..}|Mul{..}|Sq{..} |Min{..}|Max{..}|LessOrEqual{..} => [constants, SFU, ALU+1], //TODO: FMA
					_ => unreachable!()
					},
					Expression::Block{..} => [constants, SFU, ALU]
				}
			}
			f.statements.iter().map(|statement| {
				use Statement::*;
				match statement {
					Value{value,..} => visitor(value),
					Select{condition, true_exprs, false_exprs, ..} =>
						[visitor(condition)].into_iter().chain(true_exprs.iter().map(|e| visitor(e) )).chain(false_exprs.iter().map(|e| visitor(e) )).reduce(add).unwrap()
				}
			}).reduce(add).unwrap()
		};
		//fn gcd(mut a: usize, mut b: usize) -> usize { while b != 0 { (a,b) = (b, a % b) } a } //ALU/gcd(SFU,ALU), SFU/gcd(SFU,ALU)
		println!("{name} {SFU} {ALU} {:.0} {constants}", ALU as f64/SFU as f64);//{}+{}/{}, (ALU%SFU)/gcd(SFU,ALU), SFU/gcd(SFU,ALU));

		let mut last_use = vec![0; f.values.len()];
		fn visitor(last_use: &mut Vec<usize>, program_counter: usize, e: &Expression) {
			if let Expression::Expr(Expr::Value(id)) = e { last_use[id.0] = program_counter; } else { e.visit(|e| visitor(last_use, program_counter, e)); }
		}
		for (program_counter, statement) in f.statements.iter().enumerate() {
			use Statement::*;
			match statement {
				Value{value,..} => visitor(&mut last_use, program_counter, value),
				Select{condition, true_exprs, false_exprs, ..} => {
					visitor(&mut last_use, program_counter, condition);
					for e in &**true_exprs { visitor(&mut last_use, program_counter, e); }
					for e in &**false_exprs { visitor(&mut last_use, program_counter, e); }
				}
			}
		}
		for (output_counter, value) in f.output.iter().enumerate() { visitor(&mut last_use, f.statements.len()+output_counter, value); }

		let mut registers: Vec<Option<ast::Value>> = vec![];
		let mut map = vec![None; f.values.len()];
		for (program_counter, statement) in f.statements.iter_mut().enumerate() {
			for (_r, register) in registers.iter_mut().enumerate() { if let Some(id) = register { if program_counter>=last_use[id.0] {
				//eprintln!("-{_r}: {}", names[id.0]);
				*register= None;
				}
			} }
			use Statement::*;
			let mut visitor = |id:ast::Value| if id.0 >= f.input.len() {
				if !registers.contains(&Some(id)) {
					let r =
						if let Some((r, register)) = registers.iter_mut().enumerate().filter(|(_,register)| register.is_none()).next() { *register = Some(id); r }
						else { registers.push(Some(id)); registers.len()-1 };
					map[id.0] = Some(r.try_into().unwrap());
				}
			};
			match statement {
				Value{id,..} => visitor(*id),
				Select{results,..} => for id in &**results { visitor(*id) },
			}
		}
		map = vec![None; f.values.len()]; // Clear register allocation
		let mut s = self::compile(0, &f, &map, name);
		s = s.replace("__global__","__NEKRK_DEVICE__").replace("const unsigned int id = blockIdx.x * blockDim.x + threadIdx.x;\n","");
		for i in 0..f.input.len() { s = s.replace(&format!("in{i}[]"),&f.values[i]).replace(&format!(/*const*/"double {} = in{i}[id];\n",&f.values[i]),""); }
		s = s.replace(&format!(", double {array_input}0"),&format!(", double {array_input}[]"));
		for k in (0..K).rev() { s = s.replace(&format!(", double {array_input}{k},"),",").replace(&format!("{array_input}{k}"), &format!("{array_input}[{k}]")); }
		if f.output.len() == 1 { s = s.replace(", double out0[]","").replace("out0[id] = ","return ").replace("__ void","__ double"); }
		else if direct {
			//s = s.replace(", double out0[]",", double* out[]");
			//for k in 0..K { s = s.replace(&format!(", double out{k}[]"),"").replace(&format!("out{k}[id]"), &format!("out[{k}][id]")); }
			s = s.replace(", double out0[]",", double out[], unsigned int n_states");
			for k in 0..K { s = s.replace(&format!(", double out{k}[]"),"").replace(&format!("out{k}[id]"), &format!("out[{k}*n_states+id]")); }
		} else {
			s = s.replace(", double out0[]",", double out[]");
			for k in 0..K { s = s.replace(&format!(", double out{k}[]"),"").replace(&format!("out{k}[id]"), &format!("out[{k}]")); }
		}
		s.replace("+-","-").replace("double","dfloat")
	};
	let target = std::path::PathBuf::from(args().skip(2).next().unwrap_or(std::env::var("HOME")?+"/nekRK/share/mechanisms"));
	eprintln!("{target:?}");
	let name = std::path::Path::new(&path).file_stem().unwrap();
	let name = name.to_str().unwrap();
	let prefix = args().skip(3).next().unwrap_or("nekrk_".to_string());
	let write = |module, content| std::fs::write(target.join(format!("{name}/{module}.c")), content);
	#[cfg(feature="transport")] {
	write("conductivity",
		compile(conductivityNIVT, format!("{prefix}conductivityNIVT"), "mole_proportions", false)
	)?;
	write("viscosity",
		[
			"__NEKRK_DEVICE__ double sq(double x) { return x*x; }",
			&compile(viscosityIVT, format!("{prefix}viscosityIVT"), "mole_proportions", false)
				.replace("rcp_VviscosityIVVT","r").replace("VviscosityIVVT","v")
		].join("\n")
	)?;
	write("diffusivity",
		compile(density_diffusivity, format!("{prefix}density_diffusivity"), "mole_proportions", true)
		.replace("density_diffusivity(","density_diffusivity(unsigned int id, ")
	)?;
	}
	write("rates",
		compile(rates, format!("{prefix}species_rates"), "concentrations", false)
	)?;

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
