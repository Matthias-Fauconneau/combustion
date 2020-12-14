.PHONY: all

all: main.ptx

main.ptx: %.ptx: %.cu
	/opt/cuda/bin/nvcc -DSPECIES=9 -DELEMENTARY=22 -DTHREE_BODY=5 -DFALLOFF=1 --ptx $< -o $@

