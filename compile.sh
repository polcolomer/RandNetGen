#!/bin/bash

# You need to have the GSL library installed

# to be compiled with a mac without any change in the interpreter
include_path="-I/opt/local/include/ -I/usr/local/include/ -Iinclude/ "
lib_path="-L/opt/local/lib/ -L/usr/local/lib/"

# we compile the included files
gcc -c -O2 -Wall $include_path ./src/graph_library.c
gcc -c -O2 -Wall $include_path ./src/rewiring.c
gcc -c -O2 -Wall $include_path ./src/annealingCk.c
gcc -c -O2 -Wall $include_path ./src/annealingCbar.c
gcc -c -O2 -Wall $include_path ./src/annealingTRI.c
gcc -c -O2 -Wall $include_path ./src/annealingPkkCk.c
gcc -c -O2 -Wall $include_path ./src/annealingPkkCbar.c
gcc -c -O2 -Wall $include_path ./src/annealingPkkTRI.c
gcc -c -O2 -Wall $include_path ./src/annealingKnn.c

# we compile the program
gcc ./src/main.c -o RandGenNet $include_path $lib_path graph_library.o rewiring.o annealingCk.o annealingCbar.o annealingTRI.o annealingPkkCk.o annealingPkkCbar.o annealingPkkTRI.o annealingKnn.o -lm -O2 -Wall -lgsl -lgslcblas

# we remove the objects
rm *.o
