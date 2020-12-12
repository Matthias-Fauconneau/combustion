#include </opt/cuda/include/cuda.h>
#include <fstream>
#include <ios>
#include <vector>
int main() {
	cuInit(0);
	CUdevice device;
	cuDeviceGet(&device, 0);
	CUcontext context;
	cuCtxCreate_v2(&context, CU_CTX_SCHED_BLOCKING_SYNC, device);
	std::ifstream file("main.ptx", std::ios::binary | std::ios::ate);
	std::streamsize fsize = file.tellg();
	file.seekg(0, std::ios::beg);
	std::vector<char> buffer(fsize);
	file.read(buffer.data(), fsize);
	CUmodule module;
	cuModuleLoadData(&module, &buffer[0]);
	CUstream stream;
	cuStreamCreateWithPriority(&stream, CU_STREAM_NON_BLOCKING, 0);
	auto stride = 1;
	auto len = 1*stride;
	CUdeviceptr temperature;
	cuMemAlloc_v2(&temperature, len*sizeof(double));
	auto host_temperature = new double[len];
	for(size_t i=0; i<len; i++) host_temperature[i] = 1000.;
	cuMemcpyHtoD_v2(temperature, host_temperature, len*sizeof(double));
	auto S = 9;
	CUdeviceptr amounts;
	cuMemAlloc_v2(&amounts, (S-1)*len*sizeof(double));
	auto host_amounts = new double[(S-1)*len];
	for(size_t k=0; k<S-1; k++)	for(size_t i=0; i<len; i++) host_amounts[k*len+i] = (double[]){0., 0., 4.874638549881687, 2.4373192749408434, 0., 0., 0., 0.}[k];
	cuMemcpyHtoD_v2(amounts, host_amounts, (S-1)*len*sizeof(double));
	CUfunction function;
	cuModuleGetFunction(&function, module,"dt");
	auto pressure_r = 101325./8.31446261815324;
	cuLaunchKernel(function, len, 1, 1, stride, 1, 1, 0, stream, (void*[]){&pressure_r, &temperature, &amounts}, nullptr);
	cuStreamSynchronize(stream);
	cuMemcpyDtoH_v2(host_temperature, temperature, len*sizeof(double));
	cuMemcpyDtoH_v2(host_amounts, amounts, (S-1)*len*sizeof(double));
	printf("%f ", host_temperature[0]);
	for(size_t k=0; k<S-1; k++) printf("%f ", host_amounts[k*len]);
	printf("\n");
}
