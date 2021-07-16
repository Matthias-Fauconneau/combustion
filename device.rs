use {anyhow::Result, iter::map};
#[allow(non_camel_case_types)] pub type float = f32;
type Output = Result<Box<[Box<[float]>]>>;
use ast::*;
fn demote(mut f: ast::Function) -> ast::Function {
	use Expr::*;
	fn demote(e: &mut Expression) { if let Expression::Expr(F64(ref x)) = e { *e = F32(R32::new(f64::from(*x) as _).unwrap()).into() } else { e.visit_mut(demote) } }
	for s in &mut *f.statements { use Statement::*; match s {
		Value{value,..} => demote(value),
		Select { condition, true_exprs, false_exprs, .. } => {
			demote(condition);
			for e in &mut **true_exprs { demote(e); }
			for e in &mut **false_exprs { demote(e); }
		}
	}}
	f
}
#[cfg(not(feature="vpu"))] mod device {
	use {iter::{list, map}, ast::*, super::*};
	pub fn assemble<'t>(function: &'t Function) -> impl 't+Fn(&[float], &[&[float]]) -> super::Output {
		let input_len = function.input;
		let output_len = function.output.len();
		#[cfg(feature="ir")] let function = ir::assemble(ir::compile(function));
		move |constants:&[float], inputs:&[&[float]]| {
			assert!(constants.len() == 1 && constants.len()+inputs.len() == input_len);
			let states_len = inputs[0].len();
			let mut outputs = map(0..output_len, |_| vec![0.; states_len].into_boxed_slice());
			let time = std::time::Instant::now();
			let mut input = list(constants.iter().copied().chain(inputs.iter().map(|_| 0.)));
			let mut output = vec![0.; output_len].into_boxed_slice();
			for state_id in 0..states_len {
				for (input, array) in input[constants.len()..].iter_mut().zip(inputs) { *input = array[state_id]; }
				function(&input, &mut output);
				for (array, output) in outputs.iter_mut().zip(&*output) { array[state_id] = *output; }
			}
			let time= time.elapsed().as_secs_f64();
			if states_len > 1 { println!("{:.0}K in {:.0}ms = {:.0}ns, {:.2}M/s", states_len as f64/1e3, time*1e3, time/(states_len as f64)*1e9, (states_len as f64)/1e6/time); }
			Ok(outputs)
		}
	}
}
#[cfg(feature="vpu")] mod device {
use {std::borrow::Borrow, iter::{list, map}, vulkan::*, super::*};
pub struct Function<D: Borrow<Device>> {
	device: D,
	input_len: usize,
	output_len: usize,
	function: Box<[u32]>,
}
pub fn call<D: Borrow<Device>>(Function{input_len, output_len, device, function}: &Function<D>, constants: &[float], input: &[&[float]]) -> super::Output {
	assert!(std::mem::size_of::<float>() == 4 && constants.len() == 1 && constants.len()+input.len() == *input_len);
	let states_len = input[0].len();
	let local_size = std::cmp::min(512, states_len as u32);
	let device = device.borrow();
	let input = map(&*input, |array| Buffer::new(device, array).unwrap());
	let output = map(0..*output_len, |_| Buffer::new(device, &vec![0.; states_len]).unwrap());
	let buffers = list(input.iter().chain(&*output));
	let pipeline = device.pipeline(&function, local_size, as_u8(&constants), &buffers)?;
	device.bind(pipeline.descriptor_set, &buffers)?;
	let command_buffer = device.command_buffer(&pipeline, as_u8(&constants), (states_len as u32)/local_size)?;
	for _ in 0..1 { // First iteration loads from host visible to device local memory
		let time = device.submit_and_wait(command_buffer)?;
		println!("{local_size}: {:.0}K in {:.0}ms = {:.0}ns, {:.2}M/s", states_len as f32/1e3, time*1e3, time/(states_len as f32)*1e9, (states_len as f32)/1e6/time);
	}
	Ok(map(&*output, |array| (*array.map(device).unwrap()).into()))
}

/*#[throws] pub fn assemble<'t>(device: &'t Device, function: &ast::Function) -> Function<&'t Device> {
	Function{device, input_len: function.input, output_len: function.output.len(), function: spirv::compile(1, function)?}
}*/
pub fn assemble(function: ast::Function) -> Function<Device> {
	Function{device: vulkan::Device::new().unwrap(), input_len: function.input, output_len: function.output.len(), function: spirv::compile(1, &demote(function)).unwrap()}
}

impl<D: Borrow<Device>> FnOnce<(&[float], &[&[float]])> for Function<D> { type Output = super::Output; extern "rust-call" fn call_once(mut self, args: (&[float], &[&[float]])) -> Self::Output { self.call_mut(args) } }
impl<D: Borrow<Device>> FnMut<(&[float], &[&[float]])> for Function<D> { extern "rust-call" fn call_mut(&mut self, args: (&[float], &[&[float]])) -> Self::Output { self.call(args) } }
impl<D: Borrow<Device>> Fn<(&[float], &[&[float]])> for Function<D> { extern "rust-call" fn call(&self, (constants, input): (&[float], &[&[float]])) -> Self::Output { call(&self, constants, input) } }
}
pub use {device::*, ast::let_};
pub fn all_same(array:&[float], times: usize) -> float { assert!(array.len() == times); for &v in array { assert_eq!(v, array[0]); } array[0] }
pub fn with_repetitive_input(f: impl Fn(&[float],&[&[float]])->Output, times: usize) -> impl Fn(&[float],&[float])->Result<Box<[float]>> {
	move |constants, inputs| Ok(map(&*f(&map(constants, |x| *x as _), &map(&*map(inputs, |x| vec![*x as _; times]), |x| &**x))?, |y| all_same(y, times)))
}
