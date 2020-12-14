.PHONY: all

all: main

main.ptx: %.ptx: %.cu system.h
	/opt/cuda/bin/nvcc -DSPECIES=9 -DELEMENTARY=22 -DTHREE_BODY=5 -DFALLOFF=1 --ptx $< -o $@

main: %: %.cpp %.ptx
	LD_LIBRARY_PATH=/opt/cuda/lib64 c++ $< -o $@ -lcuda -g -I/opt/cuda/include
