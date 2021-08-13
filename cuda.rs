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

pub fn compile(constants_len: usize, ast: ast::Function) -> String {
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
	format!(r#"__global__ void function({parameters}) {{
const unsigned int id = blockIdx.x * blockDim.x + threadIdx.x;
{input_values}
{instructions}
{store}
}}"#).replace("f64","double")
}

mod yaml;
use {anyhow::{Result, Context, anyhow as err}, num::sq, iter::Dot, std::env::*, combustion::*};

fn main() -> Result<()> {
	let path = args().skip(1).next().unwrap_or("gri30".to_string());
	let path = if std::path::Path::new(&path).exists() { path } else { format!("/usr/share/cantera/data/{path}.yaml") };
	let model = yaml::Loader::load_from_str(std::str::from_utf8(&std::fs::read(&path).context(path)?)?)?;
	let model = yaml::parse(&model);
	let (ref _species_names, ref species@Species{ref molar_mass, ref thermodynamics, ..}, _, _, State{ref amounts, temperature, pressure_R, ..}) = new(&model);
	let total_amount = amounts.iter().sum::<f64>();
	let mole_fractions = map(&**amounts, |n| n/total_amount);
	let length = 1.;
	let velocity = 1.;
	let time = length / velocity;
	let concentration = pressure_R / temperature;
	let mean_molar_mass = zip(&**molar_mass, &*mole_fractions).map(|(m,x)| m*x).sum::<f64>();
	let density = concentration * mean_molar_mass;
	let Vviscosity = f64::sqrt(density * time) * velocity;
	let mean_molar_heat_capacity_at_CP_R:f64 = thermodynamics.iter().map(|a| a.molar_heat_capacity_at_constant_pressure_R(temperature)).dot(mole_fractions);
	let R = kB*NA;
	let thermal_conductivity = mean_molar_heat_capacity_at_CP_R * R / mean_molar_mass * density * length * velocity;
	let mixture_diffusion = sq(length) / time;
	let function = transport::properties::<4>(&species, temperature, Vviscosity, thermal_conductivity, density*mixture_diffusion);
	let function = compile(0, function);
	if std::fs::metadata("/var/tmp/main.cu").map_or(true, |cu|
		std::fs::read("/var/tmp/main.cu").unwrap() != function.as_bytes() ||
		std::fs::metadata("/var/tmp/main.ptx").map_or(true, |bin| bin.modified().unwrap() < cu.modified().unwrap())
	) {
		std::fs::write("/var/tmp/main.cu", &function)?;
		std::process::Command::new("nvcc").args(&["--ptx","/var/tmp/main.cu","-o","/var/tmp/main.ptx"]).spawn().unwrap_or_else(|_| panic!("{function}")).wait()?.success().then_some(()).ok_or(err!(""))?;
	}
	use cuda::{init, prelude::*};
	init(CudaFlags::empty())?;
	let device = Device::get_device(0)?;
	let _context = Context::create_and_push(ContextFlags::SCHED_BLOCKING_SYNC, device)?;
	let module = Module::load_from_file(&std::ffi::CString::new("/var/tmp/main.ptx")?)?;
	let stream = Stream::new(/*irrelevant*/StreamFlags::NON_BLOCKING, None)?;
	let_!{ input/*@[total_amount, temperature, nonbulk_amounts @ ..]*/ =
		&*map([&[total_amount, temperature], &amounts[0..amounts.len()-1]].concat(), |x| DeviceBuffer::from_slice(&[x]).unwrap()) => {
	let_!{ output/*@[conductivity, viscosity, diffusion @ ..]*/ = &*map(0..(2+species.len()), |_| unsafe{DeviceBuffer::zeroed(1)}.unwrap()) => {
	let function = module.get_function(&std::ffi::CString::new("function")?)?;
	//let start = std::time::Instant::now();
	unsafe{stream.launch(&function, 1, 1, 0, &*map(input.iter().chain(&*output), |x| x as *const _ as *mut ::std::ffi::c_void))}?;
	unsafe{libc::_exit(0)} // Exit process without dropping DeviceBuffers (Failed to deallocate CUDA Device memory.: InvalidValue)
	Ok(())
}}}}}
