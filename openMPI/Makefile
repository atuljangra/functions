CC := mpic++
RM := rm -rf
MKDIR := mkdir -p
MV := mv

CFLAGS := -O3 -c
LINK := -o

SRC=.
BIN=bin

sources := $(SRC)/mpi-average.cpp
objects := $(BIN)/mpi-average.o
target := $(BIN)/mpi-average

.PHONY: clean

default:
	$(MKDIR) $(BIN)
	$(CC) $(CFLAGS) $(sources)
	$(MV) *.o $(BIN)
	$(CC) $(LINK) $(target) $(objects)
	
clean: 
	$(RM) $(BIN)
