#![feature(format_args_capture,in_band_lifetimes,default_free_fn,associated_type_bounds,iter_is_partitioned)]#![allow(non_snake_case)]
mod vulkan;

use {iter::{box_, map}, ast::*, vulkan::*, fehler::throws, anyhow::{Error,Result}};
#[throws] fn compile(device: &'t Device, function: &Function) -> impl 't+Fn(&[f32], &[&[f32]]) -> Result<Box<[Box<[f32]>]>> {
	let input_len = function.input;
	let output_len = function.output.len();
	let function = spirv::compile(1, function)?;
	move |uniforms:&[f32], input:&[&[f32]]| {
		assert!(uniforms.len() == 1 && uniforms.len()+input.len() == input_len);
		let local_size = 512;
		let states_len = input[0].len();
		let input = map(&*input, |array| Buffer::new(device, array.iter().copied()).unwrap());
		let output = map(0..output_len, |_| Buffer::new(device, vec![0.; states_len]).unwrap());
		let buffers = box_(input.iter().chain(&*output));
		pub fn cast<T>(slice: &[T]) -> &[u8] { unsafe{std::slice::from_raw_parts(slice.as_ptr() as *const u8, slice.len() * std::mem::size_of::<T>())} }
		let pipeline = device.pipeline(&function, local_size, cast(&uniforms), &buffers)?;
		device.bind(pipeline.descriptor_set, &buffers)?;
		let command_buffer = device.command_buffer(&pipeline, cast(&uniforms), (states_len as u32)/local_size)?;
		let time = device.submit_and_wait(command_buffer)?;
		println!("{}: {:.0}K in {:.0}ms = {:.0}ns, {:.2}M/s", local_size, states_len as f32/1e3, time*1e3, time/(states_len as f32)*1e9, (states_len as f32)/1e6/time);
		Ok(map(&*output, |array| (*array.map(device).unwrap()).into()))
	}
}
