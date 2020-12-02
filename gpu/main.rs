fn main() -> Result<(), Box<dyn std::error::Error>> {
	use {std::{sync::Arc, ffi::CStr}, vulkano::{
		instance::*, device::*, pipeline::{shader::ShaderModule, ComputePipeline}, buffer::*, descriptor::{pipeline_layout::*, descriptor::*, descriptor_set::*}, command_buffer::*, sync::{self, *}}};
	let instance = Instance::new(None, &InstanceExtensions::none(), None)?;
	let physical = PhysicalDevice::enumerate(&instance).next().unwrap();
	let queue_family = physical.queue_families().find(|&q| q.supports_compute()).unwrap();
	let extensions = &DeviceExtensions{khr_storage_buffer_storage_class: true, ..DeviceExtensions::none()};
	let (device, mut queues) = Device::new(physical, physical.supported_features(), extensions, [(queue_family, 0.5)].iter().cloned())?;
	let queue = queues.next().unwrap();
	let pipeline = Arc::new({
		#[derive(Clone)] struct Main;
		unsafe impl PipelineLayoutDesc for Main {
			fn num_sets(&self) -> usize { 1 }
			fn num_bindings_in_set(&self, set: usize) -> Option<usize> { match set { 0 => Some(1), _ => None } }
			fn descriptor(&self, set: usize, binding: usize) -> Option<DescriptorDesc> {
				match (set, binding) {
					(0, 0) => Some(
										DescriptorDesc{ty: DescriptorDescTy::Buffer(DescriptorBufferDesc{dynamic: Some(false), storage: true}), array_count: 1, stages: ShaderStages::compute(), readonly: true}
					),
					_ => None,
				}
			}
			fn num_push_constants_ranges(&self) -> usize { 0 }
			fn push_constants_range(&self, _: usize) -> Option<PipelineLayoutDescPcRange> { None }
		}
		let shader = unsafe{ ShaderModule::from_words(device.clone(), vk_shader_macros::include_glsl!("main.comp"))? };
		ComputePipeline::new(device.clone(), &unsafe{ shader.compute_entry_point(CStr::from_bytes_with_nul(b"main\0") ?, Main) }, &())?
	});
	let layout = pipeline.layout().descriptor_set_layout(0).unwrap();
	let buffer = {
		let data = (0..65536).map(|n| n);
		CpuAccessibleBuffer::from_iter(device.clone(), BufferUsage::all(), false, data)?
	};
	let set = Arc::new(PersistentDescriptorSet::start(layout.clone()).add_buffer(buffer.clone())?.build()?);
	let mut builder = AutoCommandBufferBuilder::primary_one_time_submit(device.clone(), queue.family())?;
	builder.dispatch([1024, 1, 1], pipeline.clone(), set.clone(), ())?;
	let command_buffer = builder.build()?;
	let future = sync::now(device.clone()).then_execute(queue.clone(), command_buffer)?.then_signal_fence_and_flush()?;
	future.wait(None).unwrap();
	let buffer = buffer.read().unwrap();
	for n in 0..65536 { assert_eq!(buffer[n as usize], n * 12); }
	println!("OK");
	Ok(())
}
