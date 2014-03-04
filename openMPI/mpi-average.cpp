/*
 MPI program to perform the template 
 (averaging) based computations on a matrix.
*/

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#define N_ARGS 14 // Number of arguments required by the program.
#define FAILURE 1
#define SUCCESS 0

void readInput(int argc, char *argv[]){
  
  if(argc < N_ARGS){
    printf("Please specify proper arguments - ");
    printf("usage:- %s infile outfile m n p q s l a b c d e\n",argv[0]);
    exit(FAILURE);
  }
  return;
}


int main( int argc, char *argv[]){
  readInput(argc,argv);
  return FAILURE;
}
