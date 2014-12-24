CC=gcc
FLAGS=-g -std=c99 -Wall -Wpedantic -Werror -Wfatal-errors -g
CUDA_CC=nvcc
CUDA_FLAGS=-g -I /usr/local/cuda/include -I ../../cuda/cub
LINK_CC=g++
LINK_FLAGS=-L/usr/local/cuda/lib64 -lcuda -lcudart -ltbb
CC_CPP=g++
CPP_FLAGS= -std=c++11 -Wall -Wpedantic -Werror -Wfatal-errors -g

all: lkt
lkt: lkt.o main.o lktcuda.o nocuda.o
	$(LINK_CC) main.o lkt.o lktcuda.o nocuda.o -o lkt -lm $(LINK_FLAGS)
main.o: 
	$(CC) $(FLAGS) -c main.c -o main.o
lkt.o:
	$(CC) $(FLAGS) -c lkt.c -o lkt.o
lktcuda.o:
	$(CUDA_CC) $(CUDA_FLAGS) -c lkt.cu -o lktcuda.o
nocuda.o:
	$(CC_CPP) $(CPP_FLAGS) -c nocuda.cpp -o nocuda.o
clean:
	rm -f *o lkt
