#include <cuda.h>
#include <fstream>
#include <ios>
#include <vector>
#include <cmath>
#include <chrono>
#include <iostream>
using namespace std;
using namespace chrono;

int main(int argc, const char* argv[]) {
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
	auto workgroup_size = std::stoi(argv[1]);
	size_t len = std::stoi(argv[2]); // workgroup_count * workground_size;

	CUdeviceptr temperature;
	cuMemAlloc_v2(&temperature, len*sizeof(double));
	auto host_temperature = new double[len];
	for(size_t i=0; i<len; i++) host_temperature[i] = 1000.;
	cuMemcpyHtoD_v2(temperature, host_temperature, len*sizeof(double));

	CUdeviceptr d_temperature;
	cuMemAlloc_v2(&d_temperature, len*sizeof(double));
	auto host_d_temperature = new double[len];
	for(size_t i=0; i<len; i++) host_d_temperature[i] = NAN;
	cuMemcpyHtoD_v2(d_temperature, host_d_temperature, len*sizeof(double));

	auto S = 9;

	CUdeviceptr amounts;
	cuMemAlloc_v2(&amounts, (S-1)*len*sizeof(double));
	auto host_amounts = new double[(S-1)*len];
	for(size_t k=0; k<S-1; k++) for(size_t i=0; i<len; i++) host_amounts[k*len+i] = (double[]){0., 0., 4.874638549881687, 2.4373192749408434, 0., 0., 0., 0.}[k];
	cuMemcpyHtoD_v2(amounts, host_amounts, (S-1)*len*sizeof(double));

	CUdeviceptr d_amounts;
	cuMemAlloc_v2(&d_amounts, (S-1)*len*sizeof(double));
	auto host_d_amounts = new double[(S-1)*len];
	for(size_t k=0; k<S-1; k++) for(size_t i=0; i<len; i++) host_d_amounts[k*len+i] = NAN;
	cuMemcpyHtoD_v2(d_amounts, host_d_amounts, (S-1)*len*sizeof(double));

	CUfunction function;
	cuModuleGetFunction(&function, module,"rates");
	for(size_t i=0; i<10; i++) {
		auto pressure_r = 101325./8.31446261815324;
		void* args[]{&len, &pressure_r, &temperature, &amounts, &d_temperature, &d_amounts};
		cuStreamSynchronize(stream);
		auto start = high_resolution_clock::now();
		cuLaunchKernel(function, /*workgroup_count*/len/workgroup_size, 1, 1, workgroup_size, 1, 1, 0, stream, args, nullptr);
		cuStreamSynchronize(stream);
		auto stop = high_resolution_clock::now();
		auto mus = duration_cast<microseconds>(stop - start).count();
		cout << len/1000 <<"K in "<< mus/1000 <<"ms = "<< float(len)/(float(mus)/1e6)/1e6 << "M/s"<<endl;
		cuMemcpyDtoH_v2(host_d_temperature, d_temperature, len*sizeof(double));
		cuMemcpyDtoH_v2(host_d_amounts, d_amounts, (S-1)*len*sizeof(double));
		/*printf("T: %f n:", host_temperature[0]);
		for(size_t k=0; k<S-1; k++) printf("%f ", host_amounts[k*len]);
		printf("\n");*/
		printf("%f ", host_d_temperature[0]);
		for(size_t k=0; k<S-1; k++) printf("%f ", host_d_amounts[k*len]);
		printf("\n");
	}
}
