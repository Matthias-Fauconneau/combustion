#![feature(default_free_fn)]
#[fehler::throws(Box<dyn std::error::Error>)] pub fn main() {
	use {std::default::default, wgpu::*};
	let (device, queue) = futures::executor::block_on(async{Instance::new(BackendBit::PRIMARY).request_adapter(&default()).await.unwrap().request_device(&DeviceDescriptor{shader_validation: true, ..default()}, None).await})?;
	let ref module = device.create_shader_module(include_spirv!(env!("f.spv")));
	let bind_group_layout = device.create_bind_group_layout(&BindGroupLayoutDescriptor{label: None, entries: &[]});
	let pipeline_layout = device.create_pipeline_layout(&PipelineLayoutDescriptor{label: None, bind_group_layouts: &[&bind_group_layout], push_constant_ranges: &[]});
	let compute_pipeline = device.create_compute_pipeline(&ComputePipelineDescriptor{label: None, layout: Some(&pipeline_layout), compute_stage: ProgrammableStageDescriptor{module, entry_point: "main"}});
	let bind_group = device.create_bind_group(&BindGroupDescriptor{label: None, layout: &bind_group_layout, entries: &[]});
	let mut encoder = device.create_command_encoder(&CommandEncoderDescriptor{label: None});
	{let mut pass = encoder.begin_compute_pass();
		pass.set_bind_group(0, &bind_group, &[]);
		pass.set_pipeline(&compute_pipeline);
		pass.dispatch(1, 1, 1);}
	queue.submit(Some(encoder.finish()));
  println!("OK");
}
