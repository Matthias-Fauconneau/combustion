#![feature(default_free_fn,in_band_lifetimes, array_map)]
#[macro_use] mod include_bytes_align_as;
use {std::{default::default, mem::size_of, ffi::CStr}, ash::{*, vk::*, version::*, Device, extensions::ext::DebugUtils}};

pub struct Buffer {len: usize, buffer: vk::Buffer, memory: DeviceMemory}
struct Map<'t, T> {device: &'t Device, memory: &'t DeviceMemory, map: &'t [T]}
impl<T> Drop for Map<'_, T> { fn drop(&mut self) { unsafe { self.device.unmap_memory(*self.memory); } } }
impl<T> core::ops::Deref for Map<'t, T> { type Target = &'t [T]; fn deref(&self) -> &Self::Target { &self.map } }
struct MapMut<'t, T> {device: &'t Device, memory: &'t DeviceMemory, map: &'t mut [T]}
impl<T> Drop for MapMut<'_, T> { fn drop(&mut self) { unsafe { self.device.unmap_memory(*self.memory); } } }
impl<T> core::ops::Deref for MapMut<'t, T> { type Target = &'t mut [T]; fn deref(&self) -> &Self::Target { &self.map } }
impl<T> core::ops::DerefMut for MapMut<'t, T> { fn deref_mut(&mut self) -> &mut Self::Target { &mut self.map } }

impl Buffer {
	#[fehler::throws(Result)] unsafe fn map<T>(&'t self, device: &'t Device) -> Map<T> {
		Map{device, memory: &self.memory,
			map: std::slice::from_raw_parts(device.map_memory(self.memory, 0, (self.len * size_of::<T>()) as u64, default())? as *const T, self.len)}
	}

	#[fehler::throws(Result)] unsafe fn map_mut<T>(&'t mut self, device: &'t Device) -> MapMut<T> {
		MapMut{device, memory: &self.memory,
			map: std::slice::from_raw_parts_mut(device.map_memory(self.memory, 0, (self.len * size_of::<T>())as u64, default())? as *mut T, self.len)}
	}

	#[fehler::throws(Result)] unsafe fn new<I: ExactSizeIterator>(device: &Device, iter: I) -> Self {
		let len = iter.len();
		let buffer = device.create_buffer(&BufferCreateInfo{size: (len * size_of::<I::Item>()) as u64,
			usage: {type F = BufferUsageFlags; F::STORAGE_BUFFER|F::TRANSFER_SRC|F::TRANSFER_DST},
			sharing_mode: SharingMode::EXCLUSIVE, ..default()}, None)?;
		let memory_requirements = device.get_buffer_memory_requirements(buffer);
		let memory = device.allocate_memory(&MemoryAllocateInfo{allocation_size: memory_requirements.size, memory_type_index: 0, ..default()}, None)?;
		device.bind_buffer_memory(buffer, memory, 0)?;
		let mut buffer = Self{len, buffer, memory};
		for (item, cell) in iter.zip(buffer.map_mut(device)?.into_iter()) { *cell = item; }
		buffer
	}
}

#[fehler::throws(anyhow::Error)] fn main() {
	let system = std::fs::read("H2+O2.ron")?;
	type Simulation<'t> = combustion::Simulation::<'t, 9>;
	let Simulation{system, pressure, volume, state: combustion::State{temperature, amounts}, ..} = Simulation::new(&system)?;
	use iter::{array_from_iter as from_iter, vec::{Prefix, eval}};
	let state = {use iter::into::IntoChain; from_iter([temperature,volume].chain(*amounts.prefix::<{Simulation::species_len-1}>()))};
	let f = system.dt(pressure, &state);

	unsafe {
		let ref main = CStr::from_bytes_with_nul(b"main\0")?;
		let entry = Entry::new()?;
		let instance = entry.create_instance(&InstanceCreateInfo::builder().application_info(
			&ApplicationInfo::builder().api_version(make_version(1, 0, 0)).application_name(main).application_version(0).engine_name(main))
			.enabled_layer_names(&[CStr::from_bytes_with_nul(b"VK_LAYER_KHRONOS_validation\0")?.as_ptr()]).enabled_extension_names(&[DebugUtils::name().as_ptr()]),
			None)?;
		let debug_utils = DebugUtils::new(&entry, &instance);
		unsafe extern "system" fn vulkan_debug_callback(severity: DebugUtilsMessageSeverityFlagsEXT,
																																											r#type: DebugUtilsMessageTypeFlagsEXT,
																																											data: *const DebugUtilsMessengerCallbackDataEXT, _user_data: *mut std::os::raw::c_void) -> Bool32 {
			let data = *data;
			println!("{:?}:\n{:?} [{:?} ({})] : {:?}\n", severity, r#type, Some(data.p_message_id_name).filter(|p| !p.is_null()).map(|p| CStr::from_ptr(p)), data.message_id_number,
																																									Some(data.p_message).filter(|p| !p.is_null()).map(|p| CStr::from_ptr(p)) );
			FALSE
		}

		let _debug_utils_messenger = debug_utils.create_debug_utils_messenger(&DebugUtilsMessengerCreateInfoEXT{
			message_severity: DebugUtilsMessageSeverityFlagsEXT::ERROR|DebugUtilsMessageSeverityFlagsEXT::WARNING,
			message_type: DebugUtilsMessageTypeFlagsEXT::all(), pfn_user_callback: Some(vulkan_debug_callback), ..default()}, None)?;

		let device = *instance.enumerate_physical_devices()?.first().unwrap();
		let queue_family_index = instance.get_physical_device_queue_family_properties(device).iter().position(|p| p.queue_flags.contains(QueueFlags::COMPUTE)).unwrap() as u32;
		let ref device = instance.create_device(device, &DeviceCreateInfo::builder()
			.queue_create_infos(&[DeviceQueueCreateInfo::builder().queue_family_index(queue_family_index).queue_priorities(&[1.]).build()])
			.enabled_features(&PhysicalDeviceFeatures{shader_float64: TRUE, ..default()})
			, None)?;
		let queue = device.get_device_queue(queue_family_index, 0);

		let stride = 64;
		let len = 1*stride;

		let [pressure, temperature, volume] = [pressure, temperature, volume].map(|initial| Buffer::new(device, (0..len).map(|_| initial)).unwrap());
		let amounts = eval(amounts, |n| Buffer::new(device, (0..len).map(|_| n)).unwrap());
		type T = DescriptorType; type S = ShaderStageFlags;
		use std::iter::once;
		let bindings = (0..3).map(|i| DescriptorSetLayoutBinding{binding: i, descriptor_type: T::STORAGE_BUFFER, descriptor_count: 1, stage_flags: S::COMPUTE, ..default()})
								.chain(once(DescriptorSetLayoutBinding{binding: 3, descriptor_type: T::STORAGE_BUFFER, descriptor_count: amounts.len() as u32, stage_flags: S::COMPUTE, ..default()}));
		let descriptor_set_layouts = [device.create_descriptor_set_layout(&DescriptorSetLayoutCreateInfo::builder().bindings(&bindings.collect::<Box<_>>()), None)?];
		let spv = include_bytes_align_as!(u32, concat!(env!("OUT_DIR"), "/main.spv"));
		let code = {let (h, code, t) = spv.align_to(); assert!(h.is_empty() && t.is_empty(), "{} {} {}", h.len(), code.len(), t.len()); code};
		let module = device.create_shader_module(&ShaderModuleCreateInfo::builder().code(code), None)?;
		let layout = device.create_pipeline_layout(&PipelineLayoutCreateInfo::builder().set_layouts(&descriptor_set_layouts), None)?;
		let pipeline = device.create_compute_pipelines(default(), &[ComputePipelineCreateInfo{
			stage: PipelineShaderStageCreateInfo::builder().stage(S::COMPUTE).module(module).name(&CStr::from_bytes_with_nul(b"main\0")?).build(), layout, ..default()}],
			None).map_err(|(_,e)| e)?[0];

		let descriptor_pool = device.create_descriptor_pool(
			&DescriptorPoolCreateInfo::builder().pool_sizes(&[DescriptorPoolSize{ty: T::STORAGE_BUFFER, descriptor_count: 3+amounts.len() as u32}]).max_sets(1), None)?;
		let descriptor_set = device.allocate_descriptor_sets(&DescriptorSetAllocateInfo::builder().descriptor_pool(descriptor_pool).set_layouts(&descriptor_set_layouts))?[0];
		for (i, buffer) in [&pressure, &temperature, &volume].iter().enumerate() {
			device.update_descriptor_sets(&[WriteDescriptorSet::builder().dst_set(descriptor_set).dst_binding(i as u32).descriptor_type(T::STORAGE_BUFFER)
											.buffer_info(&[DescriptorBufferInfo{buffer: buffer.buffer, offset: 0, range: WHOLE_SIZE}]).build()], &[]);
		}
		device.update_descriptor_sets(&[WriteDescriptorSet::builder().dst_set(descriptor_set).dst_binding(3).descriptor_type(T::STORAGE_BUFFER)
											.buffer_info(&amounts.iter().map(|buffer| DescriptorBufferInfo{buffer: buffer.buffer, offset: 0, range: WHOLE_SIZE}).collect::<Box<_>>()).build()], &[]);

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

		let gpu_f =
			from_iter({use iter::into::{IntoChain, IntoMap}; (&[temperature, volume]).chain(amounts.prefix::<{Simulation::species_len-1}>()).map(|buffer:&Buffer| buffer.map(device).unwrap()[0])});
		use itertools::Itertools;
		assert!(f == gpu_f, "\n{:e}\n{:e}", f.iter().format(" "), gpu_f.iter().format(" "));
	};
	println!("OK");
}
