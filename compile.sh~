#!/bin/bash

# You need to have the GSL library installed

# to be compiled with a mac without any change in the interpreter
include_path="-I/opt/local/include/ -I/usr/local/include/"
lib_path="-L/opt/local/lib/ -L/usr/local/lib/"

# we compile the included files
gcc -c -O2 -Wall $include_path ./code/include/graph_library.c
gcc -c -O2 -Wall $include_path ./code/include/rewiring.c
gcc -c -O2 -Wall $include_path ./code/include/annealingCk.c
gcc -c -O2 -Wall $include_path ./code/include/annealingCbar.c
gcc -c -O2 -Wall $include_path ./code/include/annealingTRI.c
gcc -c -O2 -Wall $include_path ./code/include/annealingPkkCk.c
gcc -c -O2 -Wall $include_path ./code/include/annealingPkkCbar.c
gcc -c -O2 -Wall $include_path ./code/include/annealingPkkTRI.c
gcc -c -O2 -Wall $include_path ./code/include/annealingKnn.c

# we compile the program
gcc ./code/main.c -o RandGenNet $include_path $lib_path graph_library.o rewiring.o annealingCk.o annealingCbar.o annealingTRI.o annealingPkkCk.o annealingPkkCbar.o annealingPkkTRI.o annealingKnn.o -lm -O2 -Wall -lgsl -lgslcblas

# we remove the objects
rm *.o
