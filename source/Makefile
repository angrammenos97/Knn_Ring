# define the shell to bash
SHELL := /bin/bash

# define the C/C++ compiler to use,default here is clang
CC = gcc-7
MPICC = mpicc
MPIRUN = mpirun -np 4

# define compile-time flags
CFLAGS = -Wall
LDFLAGS = -lm -lopenblas

# define the source file for the library
SRC = knnring

all: test_sequential test_mpi

extract:
	tar -xvzf code.tar.gz
	rm -f ./$(SRC)/lib/*

tar :
	tar -czvf code.tar.gz ./$(SRC)

test_sequential: extract
	cd $(SRC); make lib; cd ..
	cd $(SRC); cp lib/*.a inc/$(SRC).h ../; cd ..
	$(CC) $(CFLAS) tester.c $(SRC)_sequential.a -o $@ $(LDFLAGS)
	./test_sequential

test_mpi: extract
	cd $(SRC); make lib; cd ..
	cd $(SRC); cp lib/*.a inc/$(SRC).h ../; cd ..
	$(MPICC) $(CFLAS) tester_mpi.c $(SRC)_mpi.a -o $@ $(LDFLAGS)
	$(MPIRUN) ./test_mpi

clean:
	cd $(SRC); make clean; cd ..
	$(RM) -f $(SRC)*.* test_*
