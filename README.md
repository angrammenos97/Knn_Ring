# Knn_Ring Exercise
Parallel &amp; Distributed Computer Systems Exercise 2

### File explanation:
- Knn_Ring.cpp, knnring.h, knnring.cpp, tester_helper.h, pch.h, pch.cpp are VS temporary files
- source/ contains elearning files, slightly modifided, to evaluate libraries and the knnring folder with the upload files
- source/knnring/inc header file knnring.h and the helping header tester_helper.h
- source/knnring/src source files knnring_sequential|mpi.c, all versions of mpi (s=synchronous, a=asynchronous, r=allreduce) and main_sequential|mpi.c to timer the applications end evaluate them
- source/knnring/Makefile to compile libraries into a lib folder and the main_sequential|mpi_*.c
