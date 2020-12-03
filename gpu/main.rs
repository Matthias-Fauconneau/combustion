#![feature(default_free_fn)]
#[fehler::throws(anyhow::Error)] fn main() {
	let system = std::fs::read("H2+O2.ron")?;
	type Simulation<'t> = combustion::Simulation::<'t, 9>;
	let Simulation{pressure, volume, state: combustion::State{temperature, amounts}, ..} = Simulation::new(&system)?;

	use {std::default::default, ash::{*, vk::*, version::*, Device}};
	unsafe {
		let instance = Entry::new()?.create_instance(&InstanceCreateInfo{p_application_info: &ApplicationInfo{api_version: make_version(1, 0, 0), ..default()}, ..default()}, None)?;
		let device = *instance.enumerate_physical_devices()?.first().unwrap();
		let queue_family_index = instance.get_physical_device_queue_family_properties(device).iter().position(|p| p.queue_flags.contains(QueueFlags::COMPUTE)).unwrap() as u32;
		let ref device = instance.create_device(device,
									&DeviceCreateInfo::builder().queue_create_infos(&[DeviceQueueCreateInfo::builder().queue_family_index(queue_family_index).queue_priorities(&[1.]).build()]), None)?;
		let queue = device.get_device_queue(queue_family_index, 0);

		let stride = 64;
		let len = 1*stride;
		#[allow(non_snake_case)] #[fehler::throws(Result)] unsafe fn Buffer_new<I: ExactSizeIterator>(device: &Device, iter: I) -> Buffer {
			use std::mem::size_of;
			let byte_len = iter.len() * size_of::<I::Item>();
			let buffer = device.create_buffer(&BufferCreateInfo{size: byte_len as u64, usage: {type F = BufferUsageFlags; F::STORAGE_BUFFER|F::TRANSFER_SRC|F::TRANSFER_DST},
				sharing_mode: SharingMode::EXCLUSIVE, ..default()}, None)?;
			let memory_requirements = device.get_buffer_memory_requirements(buffer);
			let memory = device.allocate_memory(&MemoryAllocateInfo{allocation_size: memory_requirements.size, memory_type_index: 0, ..default()}, None)?;
			device.bind_buffer_memory(buffer, memory, 0)?;
			let map = std::slice::from_raw_parts_mut(device.map_memory(memory, 0, byte_len as u64, default())? as *mut I::Item, byte_len/size_of::<I::Item>());
			for (item, cell) in iter.zip(map) { *cell = item; }
			device.unmap_memory(memory);
			buffer
		}

		let pressure = Buffer_new(device, (0..len).map(|_| pressure))?;
		let temperature = Buffer_new(device, (0..len).map(|_| temperature))?;
		let volume = Buffer_new(device, (0..len).map(|_| volume))?;
		let amounts = amounts.iter().map(|&n| Buffer_new(device, (0..len).map(|_| n)).unwrap()).collect::<Box<_>>();
		type T = DescriptorType; type S = ShaderStageFlags;
		use std::iter::once;
		let bindings = (0..3).map(|i| DescriptorSetLayoutBinding{binding: i, descriptor_type: T::STORAGE_BUFFER, descriptor_count: 1, stage_flags: S::COMPUTE, ..default()})
								.chain(once(DescriptorSetLayoutBinding{binding: 3, descriptor_type: T::STORAGE_BUFFER, descriptor_count: amounts.len() as u32, stage_flags: S::COMPUTE, ..default()}));
		let descriptor_set_layouts = [device.create_descriptor_set_layout(&DescriptorSetLayoutCreateInfo::builder().bindings(&bindings.collect::<Box<_>>()), None)?];
		let spv = include_bytes!(concat!(env!("OUT_DIR"), "/main.spv"));
		let code = {let (h, code, t) = spv.align_to(); assert!(h.is_empty() && t.is_empty()); code};
		let module = device.create_shader_module(&ShaderModuleCreateInfo::builder().code(code), None)?;
		let layout = device.create_pipeline_layout(&PipelineLayoutCreateInfo::builder().set_layouts(&descriptor_set_layouts), None)?;
		let pipeline = device.create_compute_pipelines(default(), &[ComputePipelineCreateInfo{
			stage: PipelineShaderStageCreateInfo::builder().stage(S::COMPUTE).module(module).name(&std::ffi::CStr::from_bytes_with_nul(b"main\0")?).build(), layout, ..default()}],
			None).map_err(|(_,e)| e)?[0];

		let descriptor_pool = device.create_descriptor_pool(
			&DescriptorPoolCreateInfo::builder().pool_sizes(&[DescriptorPoolSize{ty: T::STORAGE_BUFFER, descriptor_count: 3+amounts.len() as u32}]).max_sets(1), None)?;
		let descriptor_set = device.allocate_descriptor_sets(&DescriptorSetAllocateInfo::builder().descriptor_pool(descriptor_pool).set_layouts(&descriptor_set_layouts))?[0];
		for (i, buffer) in [pressure, temperature, volume].iter().enumerate() {
			device.update_descriptor_sets(&[WriteDescriptorSet::builder().dst_set(descriptor_set).dst_binding(i as u32).descriptor_type(T::STORAGE_BUFFER)
											.buffer_info(&[DescriptorBufferInfo{buffer: *buffer, offset: 0, range: WHOLE_SIZE}]).build()], &[]);
		}
		device.update_descriptor_sets(&[WriteDescriptorSet::builder().dst_set(descriptor_set).dst_binding(3).descriptor_type(T::STORAGE_BUFFER)
											.buffer_info(&amounts.iter().map(|buffer| DescriptorBufferInfo{buffer: *buffer, offset: 0, range: WHOLE_SIZE}).collect::<Box<_>>()).build()], &[]);

		let command_pool =
			device.create_command_pool(&CommandPoolCreateInfo{flags: CommandPoolCreateFlags::RESET_COMMAND_BUFFER, queue_family_index, ..default()}, None)?;
		let command_buffer =
			device.allocate_command_buffers(&CommandBufferAllocateInfo{command_pool, level: CommandBufferLevel::PRIMARY, command_buffer_count: 1, ..default()})?[0];
		device.begin_command_buffer(command_buffer, &CommandBufferBeginInfo{flags: CommandBufferUsageFlags::ONE_TIME_SUBMIT, ..default()})?;
		device.cmd_bind_pipeline(command_buffer, PipelineBindPoint::COMPUTE, pipeline);
		device.cmd_bind_descriptor_sets(command_buffer, PipelineBindPoint::COMPUTE, layout, 0, &[descriptor_set], &[]);
		device.cmd_dispatch(command_buffer, len/stride, 1, 1);
		device.cmd_pipeline_barrier(command_buffer, PipelineStageFlags::ALL_COMMANDS, PipelineStageFlags::HOST, default(),
			&[MemoryBarrier{src_access_mask: AccessFlags::MEMORY_WRITE, dst_access_mask: AccessFlags::HOST_READ, ..default()}], &[], &[]);
		device.end_command_buffer(command_buffer)?;
		let fence = device.create_fence(&default(), None)?;
		device.queue_submit(queue, &[SubmitInfo::builder().command_buffers(&[command_buffer]).build()], fence)?;
		device.wait_for_fences(&[fence], true, !0)?;
		device.reset_fences(&[fence])?;
	}
	println!("OK");
}
