# define the shell to bash
SHELL := /bin/bash

# define the C/C++ compiler to use,default here is clang
CC = gcc-7

# define the source file for the library
SRC = knnring

all: test_sequential

extract:
	tar -xvzf code.tar.gz
	rm -f ./$(SRC)/lib/*

test_sequential: extract
	cd $(SRC); make lib; cd ..
	cd $(SRC); cp lib/*.a inc/$(SRC).h ../; cd ..
	$(CC) tester.c $(SRC)_sequential.a -o $@ -lm
	./test_sequential

clean:
	cd $(SRC); make clean; cd ..
	$(RM) -f $(SRC)*.* test_*
