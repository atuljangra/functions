/*
 MPI program to perform the template 
 (averaging) based computations on a matrix.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>

#define N_ARGS 14 // Number of arguments required by the program.
#define FAILURE 1
#define SUCCESS 0
#define DELM " "

void readMatrix(FILE* infile, float *A, int m, int n){
  char *line = NULL;
  size_t len = 0;
  ssize_t read = -1;
  int i=0, j=0;
  char *token = NULL;
  
  while ((read = getline(&line, &len, infile)) != -1) {
    for( j=0; j<n ; j++, line=NULL){
      token = strtok(line, DELM);
      if (token == NULL)
        break;
      A[i*n+j] = atof(token);
    }
    i++;
  }

  printf("Read Matrix A:- \n");
  for (i=0; i<m; i++) {
    printf("\t| ");
    for (j=0; j<n; j++)
      printf("%4.4f ", A[i*n+j]);
    printf("|\n");
  }
}

void outputMatrix(FILE* outfile, float *A, int m, int n){
  
  int i,j;
  
  for(i=0; i<m; i++){
    for(j=0; j<n; j++){
      fprintf(outfile,"%f",A[i*n+j]);
      if(j<n-1)
        fprintf(outfile," ");
    }
    fprintf(outfile,"\n");
  }
  
  printf("Output Matrix A:- \n");
  for (i=0; i<m; i++) {
    printf("\t| ");
    for (j=0; j<n; j++)
      printf("%4.4f ", A[i*n+j]);
    printf("|\n");
  }
}

int main( int argc, char *argv[]){

  FILE *infile = NULL;
  FILE* outfile = NULL;
  int m, n, p, q, s, l;
  float a, b, c, d, e;
  int myrank, nprocs;
  float *A = NULL;
  
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  
  /*We don't want all the processors to execute these.*/
  if(myrank == 0){
    if(argc < N_ARGS){
      printf("Please specify proper arguments - ");
      printf("usage:- %s infile outfile m n p q s l a b c d e\n",argv[0]);
      exit(FAILURE);
    }
    /*Open files for input output.*/
    infile = fopen(argv[1],"r");
    outfile = fopen(argv[2],"w");
    if( infile == NULL || outfile == NULL){
      printf("Problem opening input output files\n");
      exit(FAILURE);
    }
  }
  
  m = atoi(argv[3]);
  n = atoi(argv[4]);
  p = atoi(argv[5]);
  q = atoi(argv[6]);
  s = atoi(argv[7]);
  l = atoi(argv[8]);
  a = atof(argv[9]);
  b = atof(argv[10]);
  c = atof(argv[11]);
  d = atof(argv[12]);
  e = atof(argv[13]);
  
  if(myrank == 0){
    A = (float *)malloc(m*n*(sizeof(float)));
    readMatrix(infile,A,m,n);
  }
    
  if(myrank == 0){
    outputMatrix(outfile,A,m,n);
    if(fcloseall() != 0)
      printf("Error closing the files\n");
  }
  MPI_Finalize();
  return FAILURE;
}
