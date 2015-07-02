#Copyright 2014 Pol Colomer de Simon. All rights reserved. Code under License GPLv3.

# compiler
CC       = gcc

# Compiler options:
INCLUDES = -I/opt/local/include/ -I/usr/local/include/ -I./include # include paths
CFLAGS   = -O2 -Wall $(INCLUDES)                                   #flags of the compiler

# Linker options:
LIBS     = -lm -lgsl -lgslcblas                         # libraries
LDFLAGS  = -L/opt/local/lib/ -L/usr/local/lib/ $(LIBS)  # library directories

# Files:
CFILES   = main.c graph_library.c rewiring.c annealingCk.c annealingCbar.c \
           annealingTRI.c annealingPkkTRI.c annealingPkkCk.c annealingPkkCbar.c \
           annealingKnn.c

# Actual files are in src directory:
SRCDIR   = src
CSRCS    = $(addprefix $(SRCDIR)/, $(CFILES))

# Objects will be created in $(OBJDIR)/*.o
OBJDIR   = src/obj
OBJECTS  = $(addprefix $(OBJDIR)/, $(CFILES:%.c=%.o))

all: $(OBJDIR) RandNetGen

depend: .depend

.depend: $(CSRCS)
	@rm -f ./.depend
	@$(CC) $(CFLAGS)  -MM $^ | sed 's|[a-zA-Z0-9_-]*\.o|$(OBJDIR)/&|' >>  ./.depend;

include .depend

RandNetGen: $(OBJDIR) $(OBJECTS)
	@echo Linking RandNetGen
	@$(CC) $(OBJECTS) -o RandNetGen $(LDFLAGS)

$(OBJDIR):
	@echo Create object directory: $@
	@mkdir -p $(OBJDIR)


$(OBJECTS): $(OBJDIR)/%.o: $(SRCDIR)/%.c
	@echo Compiling $<
	@$(CC) $(CFLAGS) -c $< -o $@


clean:
	@rm -rf $(OBJDIR) *~ src/*~ RandNetGen 

# By default make understands that each rule refers to a file/directory.
# By saying that "all" and "clean" rules are PHONY we are telling make
# that it should not expect a "./all" nor "./clean" file to be created
.PHONY: all clean depend

