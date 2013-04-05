CC = gcc
TARGET = spm
OMP_TARGET = omp-spm

SRC = khwang_symmetry.c spm_util.c sym_heuristic.c sym_greedy.c sym_exact.c main.c new.c sym_naive.c sym_util.c symm_util.c
OBJ	= $(SRC:.c=.o)
LIBRARY = ./lib/esp.a ./lib/util.a

INCLUDE = ./include

# link standard library
# -lNAME => libNAME.a
# -lm => libm.a
# 
# header files
# -IDICTORY => add dictory for header files
#
# define 
# -DNAME => #define NAME
# -DNAME=VALUE => #define NAME VALUE


# automatic variables
# $^ => list of all the pre-requisites of the rule
# $@ => target
# $< => first pre-requisite


# CFLAGS = -ansi -pedantic -Wall -W -Wconversion -Wshadow -Wcast-qual -Wwrite-strings -I$(INCLUDE)
DEBUG = -g
CFLAGS = -Wall -O3 -std=c99 -I$(INCLUDE)
OPENMP = 

#-------------------------------------------------
all: $(TARGET)

$(TARGET): $(OBJ) $(LIBRARY)
	$(CC) $(CFLAGS) -o $(TARGET) $(OBJ) $(LIBRARY)

openmp: $(OBJ) $(LIBRARY) 
	$(CC) $(CFLAGS) $(OPENMP) -o $(OMP_TARGET) $(OBJ) $(LIBRARY)

#-------------------------------------------------

clean:
	@echo "cleaning..."
	-rm -f $(TARGET) $(OMP_TARGET) *.o *.a

install:

check:
