type Output = anyhow::Result<Box<[Box<[f32]>]>>;
#[cfg(not(feature="gpu"))] mod device {
	use {iter::{list, map}, ast::*};
	pub fn assemble<'t>(function: &'t Function) -> impl 't+Fn(&[f32], &[&[f32]]) -> super::Output {
		let input_len = function.input;
		let output_len = function.output.len();
		#[cfg(feature="ir")] let function = ir::assemble(ir::compile(function));
		move |uniforms:&[f32], inputs:&[&[f32]]| {
			assert!(uniforms.len() == 1 && uniforms.len()+inputs.len() == input_len);
			let states_len = inputs[0].len();
			let mut outputs = map(0..output_len, |_| vec![0.; states_len].into_boxed_slice());
			let time = std::time::Instant::now();
			let mut input = list(uniforms.iter().copied().chain(inputs.iter().map(|_| 0.)));
			let mut output = vec![0.; output_len].into_boxed_slice();
			for state_id in 0..states_len {
				for (input, array) in input[uniforms.len()..].iter_mut().zip(inputs) { *input = array[state_id]; }
				function(&input, &mut output);
				for (array, output) in outputs.iter_mut().zip(&*output) { array[state_id] = *output; }
			}
			let time= time.elapsed().as_secs_f64();
			if states_len > 1 { println!("{:.0}K in {:.0}ms = {:.0}ns, {:.2}M/s", states_len as f64/1e3, time*1e3, time/(states_len as f64)*1e9, (states_len as f64)/1e6/time); }
			Ok(outputs)
		}
	}
}
#[cfg(feature="gpu")] mod device {
use {iter::{list, map}, ast::*, vulkan::*, super::Result};
pub struct Function<Device: AsRef<Device>> {
	device: Device,
	input_len: usize,
	output_len: usize,
	function: Box<[u32]>,
}
pub fn call(Function{input_len, output_len, device, function}: Function, uniforms: &[f32], input: &[&[f32]]) -> super::Output {
	assert!(uniforms.len() == 1 && uniforms.len()+input.len() == input_len);
	let local_size = 512;
	let states_len = input[0].len();
	let device = device.as_ref();
	let input = map(&*input, |array| Buffer::new(device, array.iter().copied()).unwrap());
	let output = map(0..output_len, |_| Buffer::new(device, vec![0.; states_len]).unwrap());
	let buffers = list(input.iter().chain(&*output));
	pub fn cast<T>(slice: &[T]) -> &[u8] { unsafe{std::slice::from_raw_parts(slice.as_ptr() as *const u8, slice.len() * std::mem::size_of::<T>())} }
	let pipeline = device.pipeline(&function, local_size, cast(&uniforms), &buffers)?;
	device.bind(pipeline.descriptor_set, &buffers)?;
	let command_buffer = device.command_buffer(&pipeline, cast(&uniforms), (states_len as u32)/local_size)?;
	let time = device.submit_and_wait(command_buffer)?;
	println!("{local_size}: {:.0}K in {:.0}ms = {:.0}ns, {:.2}M/s", states_len as f32/1e3, time*1e3, time/(states_len as f32)*1e9, (states_len as f32)/1e6/time);
	Ok(map(&*output, |array| (*array.map(device).unwrap()).into()))
}

/*#[throws] pub fn assemble<'t>(device: &'t Device, function: &ast::Function) -> Function<&'t Device> {
	Function{device, input_len: function.input, output_len: function.output.len(), function: spirv::compile(1, function)?}
}*/
#[throws] pub fn assemble(function: &ast::Function) -> Function<Device> {
	Function{device: vulkan::Device::new()?, input_len: function.input, output_len: function.output.len(), function: spirv::compile(1, function).unwrap()}
}

impl FnOnce<(&[f32], &mut [f32])> for Function { type Output = super::Output; extern "rust-call" fn call_once(mut self, args: (&[f32], &mut [f32])) -> Self::Output { self.call_mut(args) } }
impl FnMut<(&[f32], &mut [f32])> for Function { extern "rust-call" fn call_mut(&mut self, args: (&[f32], &mut [f32])) -> Self::Output { self.call(args) } }
impl Fn<(&[f32], &mut [f32])> for Function { extern "rust-call" fn call(&self, (input, output): (&[f32], &mut [f32])) -> Self::Output { call(&self, input, output); } }
}
pub use device::*;
