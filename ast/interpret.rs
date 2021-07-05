use super::*;

#[derive(PartialEq, Debug, Clone)] enum DataValue { None, Bool(bool), I32(u32), F32(f32), F64(f64) }
impl DataValue {
	fn bool(self) -> bool { if let Self::Bool(v) = self { v } else { panic!("{self:?}") } }
	fn i32(self) -> u32 { if let Self::I32(v) = self { v } else { panic!("{self:?}") } }
	fn f32(self) -> f32 { if let Self::F32(v) = self { v } else { panic!("{self:?}") } }
	#[track_caller] fn f64(self) -> f64 { if let Self::F64(v) = self { v } else { panic!("{self:?}") } }
}

#[derive(derive_more::Deref,derive_more::DerefMut)] struct State {
	#[deref]#[deref_mut] values: Vec<DataValue>,
	trace: Vec<f64>,
}

fn eval(state: &mut State, e: &Expression) -> DataValue {
	use {Expression::*, DataValue::{Bool, I32, F64, F32}};
	let result = match e {
		&Expression::I32(value) => I32(value),
		&Expression::F64(value) => F64(value),
		&Expression::F32(value) => F32(value),
		Cast(Type::I32, from) => I32(eval(state, from).f32().to_bits()),
		Cast(Type::F32, from) => F32(f32::from_bits(eval(state, from).i32())),
		And(a,b) => I32(eval(state, a).i32()&eval(state, b).i32()),
		Or(a,b) => I32(eval(state, a).i32()|eval(state, b).i32()),
		//IShLImm(_,_)|UShRImm(_,_)|IAdd(_,_)|ISub(_,_)|FPromote(_)|FCvtToSInt(_)|FCvtFromSInt(_) => panic!("{e:?}"),
		Value(id) => state[id.0].clone(),
		//Load(variable) => { self.variables[variable.0].unwrap() }
		Neg(x) => F64(-eval(state, x).f64()),
		Max(a,b) => F64(f64::max(eval(state, a).f64(), eval(state, b).f64())),
		Add(a,b) => F64(eval(state, a).f64() + eval(state, b).f64()),
		Sub(a,b) => F64(eval(state, a).f64() - eval(state, b).f64()),
		LessOrEqual(a,b) => Bool(eval(state, a).f64() <= eval(state, b).f64()),
		Mul(a,b) => F64(eval(state, a).f64() * eval(state, b).f64()),
		MulAdd(a,b,c) => { let [a,b,c] = [a,b,c].map(|x| eval(state, x).f64()); F64(f64::mul_add(a,b,c)) },
		Div(a,b) => F64(eval(state, a).f64() / eval(state, b).f64()),
		FDemote(x) => F32(eval(state, x).f64() as f32),
		FPromote(x) => F64(eval(state, x).f32() as f64),
		Sqrt(x) => F64(f64::sqrt(eval(state, x).f64())),
		Block { statements, result } => {
			run(state, statements);
			let result = eval(state, result);
			for statement in statements.iter() {
				match statement {
					Statement::Value{ id, .. } => { state[id.0] = DataValue::None; }
					//Statement::Trace{..} => {},
					_ => { unreachable!() }
				}
			}
			result
		},
		e => panic!("{e:?}"),
	};
	//assert!(result.is_finite());//, "{result}: {expression:?} {}", self.debug(expression));
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
			Branch { condition, consequent, alternative, results } => {
				let values = if eval(state, condition).bool() { consequent } else { alternative };
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

pub fn call(f: &Function, input: &[f64], output: &mut [f64]) -> Vec<f64> {
	assert!(input.len() == f.input);
	let mut state = State{ values: iter::map(input, |&v| DataValue::F64(v)).into_vec(), trace: vec![] };
	run(&mut state, &f.statements);
	for (slot, e) in output.iter_mut().zip(&*f.output) { if let DataValue::F64(v) = eval(&mut state, e) { *slot = v; } else { panic!("{e:?}"); } }
	state.trace
}

impl FnOnce<(&[f64], &mut [f64])> for Function {
	type Output = ();
	extern "rust-call" fn call_once(mut self, args: (&[f64], &mut [f64])) -> Self::Output { self.call_mut(args) }
}
impl FnMut<(&[f64], &mut [f64])> for Function {
	extern "rust-call" fn call_mut(&mut self, args: (&[f64], &mut [f64])) -> Self::Output { self.call(args) }
}
impl Fn<(&[f64], &mut [f64])> for Function {
	extern "rust-call" fn call(&self, (input, output): (&[f64], &mut [f64])) -> Self::Output { call(&self, input, output); }
}
