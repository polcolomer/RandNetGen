#Copyright 2014 Pol Colomer de Simon. All rights reserved. Code under License GPLv3.

CC   = gcc                  # compiler
VPATH = src                 # source code directory
LIBS = -lm -lgsl -lgslcblas # libraries
LFLAGS  = -L/opt/local/lib/ -L/usr/local/lib/   # library directories
INCLUDES= -I/opt/local/include/ -I/usr/local/include/ -I./include # include paths
CFLAGS  =-O2 -Wall $(INCLUDES) $(LIBS) $(LFLAGS)  #flags of the compiler

OBJECTS= src/main.o src/graph_library.o src/rewiring.o src/annealingCk.o src/annealingCbar.o src/annealingTRI.o src/annealingPkkTRI.o src/annealingPkkCk.o src/annealingPkkCbar.o src/annealingKnn.o

all: RandNetGen

RandNetGen: $(OBJECTS)
	@echo Linking RandNetGen
	@$(CC) $(OBJECTS) $(CFLAGS) -o RandNetGen

.c.o:
	@echo Compiling $<
	@$(CC) $(CFLAGS) -c $< -o $@

clean:
	@rm src/*.o *~ src/*~ RandNetGen



