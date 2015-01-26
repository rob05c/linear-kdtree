CUDA_CC=nvcc
CUDA_FLAGS=-g -std=c++11 -I /usr/local/cuda/include -I ../../cuda/cub
LINK_CC=g++
LINK_FLAGS=-g -L/usr/local/cuda/lib64 -lcuda -lcudart -ltbb
CC_CPP=g++
CPP_FLAGS=-g -std=c++11 -Wall -Wpedantic -Werror -Wfatal-errors

all: lkt
lkt: lkt.o main.o lktcuda.o nocuda.o
	$(LINK_CC) main.o lkt.o lktcuda.o nocuda.o -o lkt -lm $(LINK_FLAGS)
main.o: 
	$(CC_CPP) $(CPP_FLAGS) -c main.cpp -o main.o
lkt.o:
	$(CC_CPP) $(CPP_FLAGS) -c lkt.cpp -o lkt.o
lktcuda.o:
	$(CUDA_CC) $(CUDA_FLAGS) -c lkt.cu -o lktcuda.o
nocuda.o:
	$(CC_CPP) $(CPP_FLAGS) -c nocuda.cpp -o nocuda.o
clean:
	rm -f *o lkt
