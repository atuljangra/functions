/*
 MPI program to perform the template 
 (averaging) based computations on a matrix.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <mpi.h>

#define N_ARGS 14 // Number of arguments required by the program.
#define FAILURE 1
#define SUCCESS 0
#define DELM " "

/* TODO: If needed add some error handling while reading input.*/
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

  printf("\nRead Matrix A:- \n");
  for (i=0; i<m; i++) {
    printf("\t| ");
    for (j=0; j<n; j++)
      printf("%4.4f ", A[i*n+j]);
    printf("|\n");
  }
  printf("\n");
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
  
  printf("\nOutput Matrix A:- \n");
  for (i=0; i<m; i++) {
    printf("\t| ");
    for (j=0; j<n; j++)
      printf("%4.4f ", A[i*n+j]);
    printf("|\n");
  }
  printf("\n");
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
      MPI_Finalize();
      exit(FAILURE);
    }
    /*Open files for input output.*/
    infile = fopen(argv[1],"r");
    outfile = fopen(argv[2],"w");
    if( infile == NULL || outfile == NULL){
      printf("Problem opening input output files\n");
      MPI_Finalize();
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
  
  /*We don't want our program to continue with wrong input.*/
  assert( m>0 && n>0 && p>0 && q>0 && s>0 && l>0 );
  
  /*Leader gets complete matrix A and will transfer blocks to others.*/
  if(myrank == 0){
    A = (float *)malloc(m*n*(sizeof(float)));
    readMatrix(infile,A,m,n);
  }
  
  /*Perform l iterations of averaging over the matrix.*/
  while( (l--)>0 ){
    if(myrank == 0)
      printf("Iteration number %d\n",l+1);
  }
  
  if(myrank == 0){
    outputMatrix(outfile,A,m,n);
    printf("Output written to the file %s\n",argv[2]);
    /*Close all the files before exiting.*/
    if(fcloseall() != 0){
      printf("Error closing the files\n");
      MPI_Finalize();
      exit(FAILURE);
    }
  }
  
  MPI_Finalize();
  return SUCCESS;
}
