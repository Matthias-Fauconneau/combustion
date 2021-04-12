#[macro_use] mod include_bytes_align_as;
use {std::{default::default, mem::size_of, ffi::CStr}, ash::{*, vk::*, version::*, extensions::ext::DebugUtils}};

pub struct Device {
	_entry: Entry,
	_instance: ash::Instance,
	_debug_utils: DebugUtils,
	_debug_utils_messenger: DebugUtilsMessengerEXT,
	device: ash::Device,
	queue: Queue,
	timestamp_period: f32,
	command_pool: CommandPool,
	query_pool: QueryPool,
	fence: Fence,
}
impl std::ops::Deref for Device { type Target = ash::Device; fn deref(&self) -> &Self::Target { &self.device } }

impl Device {
	#[fehler::throws(Box<dyn std::error::Error>)] pub fn new() -> Self {
		let ref main = CStr::from_bytes_with_nul(b"main\0")?;
		unsafe {
			let entry = Entry::new()?;
			let instance = entry.create_instance(
				&InstanceCreateInfo::builder().application_info(
					&ApplicationInfo::builder().api_version(make_version(1, 2, 0)).application_name(main).application_version(0).engine_name(main)
				)
				//.enabled_layer_names(&[CStr::from_bytes_with_nul(b"VK_LAYER_KHRONOS_validation\0")?.as_ptr()])
				.enabled_extension_names(&[DebugUtils::name().as_ptr()])
				//.push_next(&mut ValidationFeaturesEXT::builder().enabled_validation_features(&[ValidationFeatureEnableEXT::DEBUG_PRINTF])
				, None
			)?;
			let debug_utils = DebugUtils::new(&entry, &instance);
			unsafe extern "system" fn vulkan_debug_callback(severity: DebugUtilsMessageSeverityFlagsEXT, r#type: DebugUtilsMessageTypeFlagsEXT, data: *const DebugUtilsMessengerCallbackDataEXT, _user_data: *mut std::os::raw::c_void) -> Bool32 {
				let data = *data;
				println!("{:?}:\n{:?} [{:?} ({})] : {:?}\n", severity, r#type, Some(data.p_message_id_name).filter(|p| !p.is_null()).map(|p| CStr::from_ptr(p)), data.message_id_number, Some(data.p_message).filter(|p| !p.is_null()).map(|p| CStr::from_ptr(p)) );
				FALSE
			}

			let _debug_utils_messenger = debug_utils.create_debug_utils_messenger(&DebugUtilsMessengerCreateInfoEXT{
				message_severity: DebugUtilsMessageSeverityFlagsEXT::ERROR|DebugUtilsMessageSeverityFlagsEXT::WARNING,
				message_type: DebugUtilsMessageTypeFlagsEXT::all(), pfn_user_callback: Some(vulkan_debug_callback), ..default()}, None)?;

			let device = *instance.enumerate_physical_devices()?.first().unwrap();
			let timestamp_period = instance.get_physical_device_properties(device).limits.timestamp_period;
			//dbg!(instance.get_physical_device_properties(device).limits.max_compute_work_group_size);
			//dbg!(instance.get_physical_device_properties(device).limits.max_compute_work_group_count);
			let queue_family_index = instance.get_physical_device_queue_family_properties(device).iter().position(|p| p.queue_flags.contains(QueueFlags::COMPUTE)).unwrap() as u32;
			let device = instance.create_device(
				device,
				&DeviceCreateInfo::builder()
					.queue_create_infos(&[DeviceQueueCreateInfo::builder().queue_family_index(queue_family_index).queue_priorities(&[1.]).build()])
					.enabled_features(&PhysicalDeviceFeatures{shader_float64: TRUE, ..default()})
					.enabled_extension_names(&[CStr::from_bytes_with_nul(b"VK_KHR_shader_non_semantic_info\0").unwrap().as_ptr()]),
				None
			)?;
			let queue = device.get_device_queue(queue_family_index, 0);
			let command_pool = device.create_command_pool(&CommandPoolCreateInfo{flags: CommandPoolCreateFlags::RESET_COMMAND_BUFFER, queue_family_index, ..default()}, None)?;
			let query_pool = device.create_query_pool(&vk::QueryPoolCreateInfo{query_type: vk::QueryType::TIMESTAMP, query_count: 2, ..default()}, None)?;
			let fence = device.create_fence(&default(), None)?;
			Self{_entry: entry, _instance: instance, _debug_utils: debug_utils, _debug_utils_messenger, device, queue, timestamp_period, command_pool, query_pool, fence}
		}
	}
}

pub struct Buffer {len: usize, buffer: vk::Buffer, memory: DeviceMemory}
pub struct Map<'t, T> {device: &'t ash::Device, memory: &'t DeviceMemory, map: &'t [T]}
impl<T> Drop for Map<'_, T> { fn drop(&mut self) { unsafe { self.device.unmap_memory(*self.memory); } } }
impl<T> core::ops::Deref for Map<'t, T> { type Target = &'t [T]; fn deref(&self) -> &Self::Target { &self.map } }
pub struct MapMut<'t, T> {device: &'t ash::Device, memory: &'t DeviceMemory, map: &'t mut [T]}
impl<T> Drop for MapMut<'_, T> { fn drop(&mut self) { unsafe { self.device.unmap_memory(*self.memory); } } }
impl<T> core::ops::Deref for MapMut<'t, T> { type Target = &'t mut [T]; fn deref(&self) -> &Self::Target { &self.map } }
impl<T> core::ops::DerefMut for MapMut<'t, T> { fn deref_mut(&mut self) -> &mut Self::Target { &mut self.map } }

impl Buffer {
	#[fehler::throws(Result)] pub fn map<T>(&'t self, device: &'t Device) -> Map<T> {
		Map{device, memory: &self.memory, map: unsafe { std::slice::from_raw_parts(device.map_memory(self.memory, 0, (self.len * size_of::<T>()) as u64, default())? as *const T, self.len)} }
	}

	#[fehler::throws(Result)] pub fn map_mut<T>(&'t mut self, device: &'t Device) -> MapMut<T> {
		MapMut{device, memory: &self.memory, map: unsafe { std::slice::from_raw_parts_mut(device.map_memory(self.memory, 0, (self.len * size_of::<T>())as u64, default())? as *mut T, self.len)} }
	}

	#[fehler::throws(Result)] pub fn new<I: ExactSizeIterator+Iterator<Item=f64/*Prevent &f64 by mistake*/>>(device: &Device, iter: I) -> Self {
		let len = iter.len();
		let mut buffer = unsafe {
			let buffer = device.create_buffer(&BufferCreateInfo{size: (len * size_of::<I::Item>()) as u64,
				usage: BufferUsageFlags::STORAGE_BUFFER|BufferUsageFlags::TRANSFER_SRC|BufferUsageFlags::TRANSFER_DST,
				sharing_mode: SharingMode::EXCLUSIVE, ..default()}, None)?;
			let memory_requirements = device.get_buffer_memory_requirements(buffer);
			let memory = device.allocate_memory(&MemoryAllocateInfo{allocation_size: memory_requirements.size, memory_type_index: 0, ..default()}, None)?;
			device.bind_buffer_memory(buffer, memory, 0)?;
			Self{len, buffer, memory}
		};
		for (item, cell) in iter.zip(buffer.map_mut(device)?.into_iter()) { *cell = item; }
		buffer
	}
}

pub struct Pipeline {
	_module: ShaderModule,
	_descriptor_pool: DescriptorPool,
	layout: PipelineLayout,
	pipeline: ash::vk::Pipeline,
	pub descriptor_set: DescriptorSet,
}
//impl std::ops::Deref for Pipeline { type Target = ash::vk::Pipeline; fn deref(&self) -> &Self::Target { &self.pipeline } }

impl Device {
	#[fehler::throws(Result)] pub fn pipeline/*<T>*/(&self, constants: &[/*T*/f64], buffers: &[&[&Buffer]]) -> Pipeline {
		let ty = DescriptorType::STORAGE_BUFFER;
		let stage_flags = ShaderStageFlags::COMPUTE;
		let bindings = buffers.iter().enumerate().map(|(binding, buffers)| DescriptorSetLayoutBinding{binding: binding as u32, descriptor_type: ty, descriptor_count: buffers.len() as u32, stage_flags, ..default()}).collect::<Box<_>>();
		let Self{device, ..} = self;
		unsafe {
			let spv = include_bytes_align_as!(u32, concat!(env!("OUT_DIR"), "/main.spv"));
			let code = {let (h, code, t) = spv.align_to(); assert!(h.is_empty() && t.is_empty(), "{} {} {}", h.len(), code.len(), t.len()); code};
			let module = device.create_shader_module(&ShaderModuleCreateInfo::builder().code(code), None)?;
			let descriptor_pool = device.create_descriptor_pool(&DescriptorPoolCreateInfo::builder().pool_sizes(&[DescriptorPoolSize{ty, descriptor_count: buffers.iter().map(|b| b.len()).sum::<usize>() as u32}]).max_sets(1), None)?;
			let descriptor_set_layouts = [device.create_descriptor_set_layout(&DescriptorSetLayoutCreateInfo::builder().bindings(&bindings), None)?];
			let layout = device.create_pipeline_layout(&PipelineLayoutCreateInfo::builder().set_layouts(&descriptor_set_layouts)
																																																																	  .push_constant_ranges(&[PushConstantRange{stage_flags, offset: 0, size: (constants.len() * std::mem::size_of::<f64/*T*/>()) as u32}]), None)?;
			let pipeline = [ComputePipelineCreateInfo{stage: PipelineShaderStageCreateInfo::builder().stage(stage_flags).module(module).name(&CStr::from_bytes_with_nul(b"main\0").unwrap()).build(), layout, ..default()}];
			let pipeline_cache = device.create_pipeline_cache(&default(), None)?;
			let pipeline = device.create_compute_pipelines(pipeline_cache, &pipeline, None).map_err(|(_,e)| e)?[0];
			std::fs::write("/var/tmp/pipeline", device.get_pipeline_cache_data(pipeline_cache)?).unwrap();
			let descriptor_set = device.allocate_descriptor_sets(&DescriptorSetAllocateInfo::builder().descriptor_pool(descriptor_pool).set_layouts(&descriptor_set_layouts))?[0];
			Pipeline{_module: module, _descriptor_pool: descriptor_pool, descriptor_set, layout, pipeline}
		}
	}
	#[fehler::throws(Result)] pub fn bind(&self, descriptor_set: DescriptorSet, buffers: &[&[&Buffer]]) {
		let descriptor_type = DescriptorType::STORAGE_BUFFER;
		let Self{device, ..} = self;
		unsafe {
			for (binding, buffers) in buffers.iter().enumerate() {
				device.update_descriptor_sets(&[WriteDescriptorSet::builder()
					.dst_set(descriptor_set).dst_binding(binding as u32).descriptor_type(descriptor_type).buffer_info(&buffers.iter().map(|buffer| DescriptorBufferInfo{buffer: buffer.buffer, offset: 0, range: WHOLE_SIZE}).collect::<Box<_>>())
					.build()], &[]);
			}
		}
	}
	#[fehler::throws(Result)] pub fn command_buffer/*<T>*/(&self, pipeline: &Pipeline, constants: &[/*T*/f64], stride: usize, len: usize) -> CommandBuffer {
		let Self{device, command_pool, query_pool, ..} = self;
		let Pipeline{descriptor_set, layout, pipeline, ..} = pipeline;
		unsafe {
			let command_buffer = device.allocate_command_buffers(&CommandBufferAllocateInfo{command_pool: *command_pool, level: CommandBufferLevel::PRIMARY, command_buffer_count: 1, ..default()})?[0];
			device.begin_command_buffer(command_buffer, &default())?;
			device.cmd_bind_pipeline(command_buffer, PipelineBindPoint::COMPUTE, *pipeline);
			device.cmd_bind_descriptor_sets(command_buffer, PipelineBindPoint::COMPUTE, *layout, 0, &[*descriptor_set], &[]);
			device.cmd_push_constants(command_buffer, *layout, ShaderStageFlags::COMPUTE, 0, std::slice::from_raw_parts(constants.as_ptr() as *const u8, constants.len() * std::mem::size_of::</*T*/f64>()));
			device.cmd_reset_query_pool(command_buffer, *query_pool, 0, 2);
			device.cmd_write_timestamp(command_buffer, PipelineStageFlags::COMPUTE_SHADER, *query_pool, 0);
			device.cmd_dispatch(command_buffer, (len/stride) as u32, 1, 1);
			device.cmd_write_timestamp(command_buffer, PipelineStageFlags::COMPUTE_SHADER, *query_pool, 1);
			device.cmd_pipeline_barrier(command_buffer, PipelineStageFlags::ALL_COMMANDS, PipelineStageFlags::HOST, default(),
				&[MemoryBarrier{src_access_mask: AccessFlags::MEMORY_WRITE, dst_access_mask: AccessFlags::HOST_READ, ..default()}], &[], &[]);
			device.end_command_buffer(command_buffer)?;
			command_buffer
		}
	}
	#[fehler::throws(Result)] pub fn submit_and_wait(&self, command_buffer: CommandBuffer) -> f32 {
		let Self{device, queue, fence, query_pool, ..} = self;
		unsafe {
			device.queue_submit(*queue, &[SubmitInfo::builder().command_buffers(&[command_buffer]).build()], *fence)?;
			device.wait_for_fences(&[*fence], true, !0)?;
			device.reset_fences(&[*fence])?;
			let mut results = vec![0; 2];
			device.get_query_pool_results::<u64>(*query_pool, 0, 2, &mut results, vk::QueryResultFlags::TYPE_64)?;
			(results[1]-results[0]) as f32 * self.timestamp_period * 1e-9
		}
	}
}
