.PHONY: all

all: main

main.ptx: %.ptx: %.cu model.h
	nvcc -DSPECIES=35 -DELEMENTARY=181 -DTHREE_BODY=6 -DPRESSURE_MODIFICATION=1 -DEFFICIENCIES=4 -DFALLOFF=25 --ptx $< -o $@

main: %: %.cpp %.ptx
	LD_LIBRARY_PATH=/opt/cuda/lib64 c++ $< -o $@ -lcuda -g -I/opt/cuda/include

clean :
	rm main.ptx main
