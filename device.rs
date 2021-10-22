use {anyhow::Result, iter::map};
type Output<T> = Result<Box<[Box<[T]>]>>;
fn convert(mut f: ast::Function) -> ast::Function {
	use ast::*;
	fn visit(_input_len: usize, e: &mut Expression) { match e {
		Expression::Expr(Expr::F64(x)) => { *e = f32(f64::from(*x) as _).expect(&format!("{:e} overflows f32, retry with f64",f64::from(*x))).into(); } // Demote constants
		_ => { e.visit_mut(|e| visit(_input_len, e)); }
	}}
	for s in &mut *f.statements { use Statement::*; match s {
		Value{value,..} => visit(f.input.len(), value),
		Select { condition, true_exprs, false_exprs, .. } => {
			visit(f.input.len(), condition);
			for e in &mut **true_exprs { visit(f.input.len(), e); }
			for e in &mut **false_exprs { visit(f.input.len(), e); }
		}
	}}
	for e in &mut *f.output { visit(f.input.len(), e); }
	f
}

#[cfg(not(feature="vpu"))] mod device {
	use {std::default::default, iter::{list, map}, ast::*, super::*};
	#[cfg(not(feature="interpret"))] pub trait Interpret = Sized;
	#[cfg(feature="interpret")] pub trait Interpret = Into<interpret::DataValue>+From<interpret::DataValue>;
	pub fn assemble<'t, T:'t+Copy+Default+Interpret>(function: Function, _block_size: usize) -> impl 't+Fn(&[T], &[&[T]]) -> Output<T> {
		let function = convert(function);
		let input_len = function.input.len();
		let output_len = function.output.len();
		#[cfg(feature="ir")] let function = ir::assemble(ir::compile(&function));
		move |constants:&[T], inputs:&[&[T]]| {
			assert!(constants.len() <= 2);
			assert_eq!(constants.len()+inputs.len(), input_len);
			let states_len = inputs[0].len();
			let mut outputs = map(0..output_len, |_| vec![default(); states_len].into_boxed_slice());
			let time = std::time::Instant::now();
			let mut input = list(constants.iter().copied().chain(inputs.iter().map(|_| default())));
			let mut output = vec![default(); output_len].into_boxed_slice();
			for state_id in 0..states_len {
				for (input, array) in input[constants.len()..].iter_mut().zip(inputs) { *input = array[state_id]; }
				#[cfg(any(feature="interpret",feature="ir"))] function(&input, &mut output);
				#[cfg(not(any(feature="interpret",feature="ir")))] compile_error!("Either feature \"interpret\" or \"ir\" must be enabled for this crate.");
				for (array, output) in outputs.iter_mut().zip(&*output) { array[state_id] = *output; }
			}
			let time= time.elapsed().as_secs_f64();
			if states_len > 1 { println!("{:.0}K in {:.0}ms = {:.0}ns, {:.2}M/s", states_len as f64/1e3, time*1e3, time/(states_len as f64)*1e9, (states_len as f64)/1e6/time); }
			Ok(outputs)
		}
	}
}
#[cfg(feature="vpu")] mod device {
use {std::{default::default, mem::size_of, borrow::Borrow}, iter::{list, map}, vulkan::*, super::*};
pub struct Function<D: Borrow<Device>> {
	pub device: D,
	input_len: usize,
	pub output_len: usize,
	pub block_size: usize,
	pub pipeline: Pipeline,
}
pub fn call<D: Borrow<Device>, T:Plain+Default>(Function{input_len, output_len, device, block_size, pipeline,..}: &Function<D>, constants: &[T], input: &[&[T]]) -> Output<T> {
	assert!(size_of::<T>() == 4 && constants.len() <= 2 && constants.len()+input.len() == *input_len);
	let states_len = input[0].len();
	let device = device.borrow();
	let input = map(&*input, |array| Buffer::new(device, array).unwrap());
	let output = map(0..*output_len, |_| Buffer::new(device, &vec![default(); states_len]).unwrap());
	let buffers = list(input.iter().chain(&*output));
	device.bind(pipeline.descriptor_set, &buffers)?;
	let command_buffer = device.command_buffer(&pipeline, as_u8(&constants), (states_len as u32)/(*block_size as u32))?;
	for _ in 0..1 { // First iteration loads from host visible to device local memory
		let time = device.submit_and_wait(command_buffer)?;
		if states_len > 1 { println!("{block_size}: {:.0}K in {:.0}ms = {:.0}ns, {:.2}M/s", states_len as f32/1e3, time*1e3, time/(states_len as f32)*1e9, (states_len as f32)/1e6/time); }
	}
	Ok(map(&*output, |array| (*array.map(device).unwrap()).into()))
}

pub fn assemble(function: ast::Function, block_size: usize) -> Function<Device> {
	let function = convert(function);
	assert!(function.input.iter().all(|&t| t==ast::Type::F32));
	let output_types = {
		let mut types = ast::Types(function.input.iter().copied().map(Some).chain((function.input.len()..function.values.len()).map(|_| None)).collect());
		for s in &*function.statements { types.push(s); }
		map(&*function.output, |output| {
			types.expr(output);
			types.rtype(output)
		})
	};
	assert!(output_types.iter().all(|&t| t==ast::Type::F32));
	let device = vulkan::Device::new().unwrap();
	let constants_len = 2;
	let pipeline = device.pipeline(
		&spirv::compile(constants_len, &function).unwrap(),
		block_size as u32,
		constants_len*size_of::<f32>(),
		function.input.len()-constants_len+function.output.len()
	).unwrap();
	Function{device, input_len: function.input.len(), output_len: function.output.len(), block_size, pipeline}
}

impl<D: Borrow<Device>, T:Plain+Default> FnOnce<(&[T], &[&[T]])> for Function<D> { type Output = self::Output<T>; extern "rust-call" fn call_once(mut self, args: (&[T], &[&[T]])) -> Self::Output { self.call_mut(args) } }
impl<D: Borrow<Device>, T:Plain+Default> FnMut<(&[T], &[&[T]])> for Function<D> { extern "rust-call" fn call_mut(&mut self, args: (&[T], &[&[T]])) -> Self::Output { self.call(args) } }
impl<D: Borrow<Device>, T:Plain+Default> Fn<(&[T], &[&[T]])> for Function<D> { extern "rust-call" fn call(&self, (constants, input): (&[T], &[&[T]])) -> Self::Output { call(&self, constants, input) } }
}
pub use device::*;
#[allow(dead_code)] pub fn all_same<T:PartialEq+Copy+std::fmt::Debug>(array:&[T], times: usize) -> T { assert!(array.len() == times); for &v in array { assert_eq!(v, array[0]); } array[0] }
#[allow(dead_code)] pub fn with_repetitive_input<T:Copy+PartialEq+std::fmt::Debug>(f: impl Fn(&[T],&[&[T]])->Output<T>, times: usize) -> impl Fn(&[T],&[T])->Result<Box<[T]>> {
	move |constants, inputs| Ok(map(&*f(&map(constants, |x| *x as _), &map(&*map(inputs, |x| vec![*x as _; times]), |x| &**x))?, |y| all_same(y, times)))
}
