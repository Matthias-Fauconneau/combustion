use super::*;

#[derive(PartialEq, Debug, Clone)] enum DataValue { None, Bool(bool), I32(u32), F32(f32), F64(f64) }
impl DataValue {
	fn bool(self) -> bool { if let Self::Bool(v) = self { v } else { panic!("{self:?}") } }
	fn i32(self) -> u32 { if let Self::I32(v) = self { v } else { panic!("{self:?}") } }
	fn f32(self) -> f32 { if let Self::F32(v) = self { v } else { panic!("{self:?}") } }
	//#[track_caller] fn f64(self) -> f64 { if let Self::F64(v) = self { v } else { panic!("{self:?}") } }
}
impl std::fmt::Display for DataValue { fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result { match self { Self::F32(v) => write!(f,"{v:e}"), _ => unimplemented!() } } }

type State = Vec<DataValue>;

fn to_string(state: &State, expression: &Expression) -> String {
	use Expression::*;
	match expression {
		Float(x) => format!("{:e}", x),
		Value(id) => format!("{}", state[id.0].clone()),
		Add(a, b) => format!("{} + {}", eval(state, a), eval(state, b)),
		Mul(a, b) => format!("{} * {}", eval(state, a), eval(state, b)),
		Exp(x) => format!("exp({})", eval(state, x)),
		e => panic!("{:?}", e),
	}
}

fn eval(state: &State, expression: &Expression) -> DataValue {
	use {Expression::*, DataValue::{Bool, I32, F64, F32}};
	let result = match expression {
		&Expression::I32(value) => I32(value),
		&Expression::F32(value) => F32(value),
		&Expression::F64(value) => F64(value),
		&Expression::Float(value) => F32(value as f32),
		Cast(Type::I32, from) => I32(eval(state, from).f32().to_bits()),
		Cast(Type::F32, from) => F32(f32::from_bits(eval(state, from).i32())),
		And(a,b) => I32(eval(state, a).i32()&eval(state, b).i32()),
		Or(a,b) => I32(eval(state, a).i32()|eval(state, b).i32()),
		//IShLImm(_,_)|UShRImm(_,_)|IAdd(_,_)|ISub(_,_)|FPromote(_)|FCvtToSInt(_)|FCvtFromSInt(_) => panic!("{e:?}"),
		Value(id) => state[id.0].clone(),
		//Load(variable) => { self.variables[variable.0].unwrap() }
		Neg(x) if let F32(x) = eval(state, x) => F32(-x),
		Neg(x) if let F64(x) = eval(state, x) => F64(-x),
		Max(a,b) if let [F32(a), F32(b)] = [a,b].map(|x| eval(state, x)) => F32(f32::max(a,b)),
		Max(a,b) if let [F64(a), F64(b)] = [a,b].map(|x| eval(state, x)) => F64(f64::max(a,b)),
		Add(a,b) if let [F32(a), F32(b)] = [a,b].map(|x| eval(state, x)) => F32(a+b),
		Add(a,b) if let [F64(a), F64(b)] = [a,b].map(|x| eval(state, x)) => F64(a+b),
		Sub(a,b) if let [F32(a), F32(b)] = [a,b].map(|x| eval(state, x)) => F32(a-b),
		Sub(a,b) if let [F64(a), F64(b)] = [a,b].map(|x| eval(state, x)) => F64(a-b),
		LessOrEqual(a,b)  if let [F32(a), F32(b)] = [a,b].map(|x| eval(state, x)) => Bool(a<=b),
		LessOrEqual(a,b)  if let [F64(a), F64(b)] = [a,b].map(|x| eval(state, x)) => Bool(a<=b),
		Mul(a,b) if let [F32(a), F32(b)] = [a,b].map(|x| eval(state, x)) => F32(a*b),
		Mul(a,b) if let [F64(a), F64(b)] = [a,b].map(|x| eval(state, x)) => F64(a*b),
		//MulAdd(a,b,c) => { let [a,b,c] = [a,b,c].map(|x| eval(state, x).f64()); F64(f64::mul_add(a,b,c)) },
		Div(a,b) if let [F32(a), F32(b)] = [a,b].map(|x| eval(state, x)) => F32(a/b),
		Div(a,b) if let [F64(a), F64(b)] = [a,b].map(|x| eval(state, x)) => F64(a/b),
		//FDemote(x) => F32(eval(state, x).f64() as f32),
		//FPromote(x) => F64(eval(state, x).f32() as f64),
		Sqrt(x) if let F32(x) = eval(state, x) => F32(f32::sqrt(x)),
		Sqrt(x) if let F64(x) = eval(state, x) => F64(f64::sqrt(x)),
		Exp(x) if let F32(x) = eval(state, x) => F32(f32::exp(x)),
		Exp(x) if let F64(x) = eval(state, x) => F64(f64::exp(x)),
		Ln{x,..} if let F32(x) = eval(state, x) => F32(f32::ln(x)),
		Ln{x,..} if let F64(x) = eval(state, x) => F64(f64::ln(x)),
		Block { statements, result } => {
			let ref mut state = state.clone();
			run(state, statements);
			let result = eval(state, result);
			/*for statement in statements.iter() {
				match statement {
					Statement::Value{ id, .. } => { state[id.0] = DataValue::None; }
					//Statement::Trace{..} => {},
					_ => { unreachable!() }
				}
			}*/
			result
		},
		expression => panic!("{expression:?}"),
	};
	if let F32(x) = result { assert!(false || f32::abs(x)<1e28, "{x}: {expression:?} {}", to_string(state, expression)); }
	result
}

fn run(state: &mut State, statements: &[Statement]) {
	for statement in statements {
		use Statement::*;
		match statement {
			Value { id, value } => {
				let value = eval(state, value);
				if state.len() <= id.0 { state.resize_with(id.0+1, || DataValue::None); }
				assert!(state[id.0] == DataValue::None);
				state[id.0] = value;
			},
			/*Store { variable, value } => {
				let value = eval(state, value);
				assert!(self.variables[variable.0].replace(value).is_none());
			}*/
			Select { condition, true_exprs, false_exprs, results } => {
				let values = if eval(state, condition).bool() { true_exprs } else { false_exprs };
				for (id, value) in results.iter().zip(&**values) {
					let value = eval(state, value);
					if state.len() <= id.0 { state.resize_with(id.0+1, || DataValue::None); }
					assert!(state[id.0] == DataValue::None);
					state[id.0] = value;
				}
			},
			/*Trace { id } => {
				state.trace.push(state[id.0].clone().f64());
			}*/
		}
	}
}

pub fn call(f: &Function, input: &[f32], output: &mut [f32]) {
	assert!(input.len() == f.input);
	let mut state = iter::map(input, |&v| DataValue::F32(v)).into_vec();
	run(&mut state, &f.statements);
	for (slot, e) in output.iter_mut().zip(&*f.output) { if let DataValue::F32(v) = eval(&mut state, e) { *slot = v; } else { panic!("{e:?}"); } }
}

impl FnOnce<(&[f32], &mut [f32])> for Function { type Output = (); extern "rust-call" fn call_once(mut self, args: (&[f32], &mut [f32])) -> Self::Output { self.call_mut(args) } }
impl FnMut<(&[f32], &mut [f32])> for Function { extern "rust-call" fn call_mut(&mut self, args: (&[f32], &mut [f32])) -> Self::Output { self.call(args) } }
impl Fn<(&[f32], &mut [f32])> for Function { extern "rust-call" fn call(&self, (input, output): (&[f32], &mut [f32])) -> Self::Output { call(&self, input, output); } }
