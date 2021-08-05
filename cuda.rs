#![feature(format_args_capture,in_band_lifetimes,default_free_fn,associated_type_bounds,unboxed_closures,fn_traits,trait_alias,iter_zip)]
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
				&F64(value) => value.to_string(),
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
				Ln{x,..} => { let x = self.expr(x); format!("ln({x})") }
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
	format!(r#"__global__ void kernel({parameters}) {{"
const unsigned int id = blockIdx.x * blockDim.x + threadIdx.x;
{input_values}
{instructions}
{store}
}}"#)
}

mod yaml;
use {anyhow::{Result, Context}, std::env::*, combustion::*};

fn main() -> Result<()> {
	let path = args().skip(1).next().unwrap();
	let path = if std::path::Path::new(&path).exists() { path } else { format!("/usr/share/cantera/data/{path}.yaml") };
	let model = yaml::Loader::load_from_str(std::str::from_utf8(&std::fs::read(&path).context(path)?)?)?;
	let model = yaml::parse(&model);
	let (ref _species_names, ref species, _, _, ref _state) = new(&model);
	let transport = transport::properties::<5>(&species);
	println!("{}", compile(1, transport));
	Ok(())
}
