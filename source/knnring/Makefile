# ####################################################################
#
#			   C/C++ Makefile
#
# Author: Anastasios Grammenos <avgramme@ece.auth.gr>
#
# Adapted from
#  http://www.cs.swarthmore.edu/~newhall/unixhelp/howto_makefiles.html
#
# ####################################################################
#
# 'make'        build executable file 'main'
# 'make lib'		build the libraries .a
# 'make clean'  removes all .o and executable files
#

# define the shell to bash
SHELL := /bin/bash

# define the C/C++ compiler to use,default here is clang
CC = gcc-7
MPICC = mpicc

# setup paths
SRCDIR = src
INCDIR = inc
LIBDIR = lib
MATDIR = matlab

# define compile-time flags
CFLAGS = -Wall
LDFLAGS = -lm -lopenblas

# define any directories containing header files
INCLUDES = -I $(INCDIR)

# define the source file for the library
SRC = knnring

# define the different possible executables
TYPES = sequential mpi

# define the executable file  name
MAIN = main

#
# The following part of the makefile is generic; it can be used to
# build any executable just by changing the definitions above
#

# call everytime
.PRECIOUS: $(LIBDIR)/%.a

all: $(addprefix $(MAIN)_, $(TYPES))

lib: $(addprefix $(LIBDIR)/, $(addsuffix .a, $(addprefix $(SRC)_, $(TYPES))))

$(MAIN)_sequential: $(SRCDIR)/$(MAIN).c $(LIBDIR)/$(SRC)_sequential.a
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ $^ $(LDFLAGS)

$(MAIN)_mpi: $(SRCDIR)/$(MAIN)_mpi.c $(LIBDIR)/$(SRC)_mpi.a
		$(MPICC) $(CFLAGS) $(INCLUDES) -o $@ $^ $(LDFLAGS)

# this is a suffix replacement rule for building .o's from .c's
# it uses automatic variables $<: the name of the prerequisite of
# the rule(a .cpp file) and $@: the name of the target of the rule (a .o file)

$(LIBDIR)/$(SRC)_%.a: $(SRCDIR)/$(SRC)_%.o
	mkdir -p $(LIBDIR)
	ar rcs $@ $<

# (see the gnu make manual section about automatic variables)
$(SRCDIR)/$(SRC)_sequential.o: $(SRCDIR)/$(SRC)_sequential.c
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $< $(LDFLAGS)

$(SRCDIR)/$(SRC)_mpi.o: $(SRCDIR)/$(SRC)_mpi.c
		$(MPICC) $(CFLAGS) $(INCLUDES) -o $@ -c $< $(LDFLAGS)

clean:
	$(RM) $(LIBDIR)/* $(SRCDIR)/*~ $(SRCDIR)/*.o $(INCDIR)/*~ $(MATDIR)/data.m $(addprefix $(MAIN)_, $(TYPES)) *~ data.m
