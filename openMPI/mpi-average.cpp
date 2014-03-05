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

void readMatrix(FILE* infile, float *A, int m, int n);
void outputMatrix(FILE* outfile, float *A, int m, int n);

int main( int argc, char *argv[]){

  FILE *infile = NULL;
  FILE* outfile = NULL;
  int m, n, p, q, s, l;
  float a, b, c, d, e;
  int myrank, nprocs;
  float *A = NULL;
  MPI_Status status;
  
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  
  /*Executed by leader only.*/
  if(myrank == 0){
    if(argc < N_ARGS){
      printf("ERROR:- Please specify proper arguments - ");
      printf("usage:- %s infile outfile m n p q s l a b c d e\n",argv[0]);
      MPI_Finalize();
      exit(FAILURE);
    }
    infile = fopen(argv[1],"r");
    outfile = fopen(argv[2],"w");
    if( infile == NULL || outfile == NULL){
      printf("ERROR:- Problem opening input output files\n");
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
  A = (float *)calloc(m*n, sizeof(float));
    
  if( !(m>0 && n>0 && p>0 && q>0 && s>0 && l>0) ){
    printf("ERROR:- Invalid Input to the program\n");
    free(A);
    MPI_Finalize();
    exit(FAILURE);
  }
    
  /* TODO:- Well more of a note.
   * One way i'm trying to do this part is for every thread to have
   * space for the A matrix and then leader thread can just send data
   * to every one and each thread can then continue to do its part of
   * the computations, at the end of which they can send the data back
   * to A.
   * 
   * Other way would have been to not reserve space for A on every node.
   * But since this is block cyclic so every thread could be allocated
   * the blocks in some weird manner. Hence for this part we need some
   * kind of tracking mechanism that which blocks do i have and where are
   * they placed in memory, maybe a linked list of blocks.
   * I was not able to find any good MPI-API to help with this and hence
   * for the sake of us having very small time left went with this 
   * approach.
   */

  /************************ Leader thread**************************/
  if(myrank == 0){
    readMatrix(infile,A,m,n);
    
    /*Send the data to all the worker threads in block cyclic fashion.*/
    int i, j, pi, pj, k, offset, dest, tag;
    
    /*Iterate over all the matrix block by block and send the data.*/
    for(i=0; i<m; i+=s){
      for(j=0; j<n; j+=s){
        
        pi = (i/s)%p ; // Row number inside the processors grid.
        pj = (j/s)%q ; // Coulmn number inside the processors grid.
        dest = (pi*q + pj) ; // Rank of the processor where we want to send.
        tag = myrank; // Tag the data with rank of the leader i.e. 0.
        
        //printf("Dest is %d\n",dest);
        if(dest != myrank){
          /*Send the data if this part does not belongs to me.*/
          for(k=i; k<i+s; k++){
            offset = k*n+j;
            MPI_Send( A+offset, s, MPI_FLOAT, dest, tag, MPI_COMM_WORLD);
          }
        }
      }
    }
  }
  /************************ Worker thread**************************/
  if(myrank > 0){
    
    /*Receive the data from the leader thread in block cyclic fashion.*/
    int i, j, pi, pj, k, offset, tag, dest;
    
    /*Iterate over all the matrix block by block and receive the data.*/
    for(i=0; i<m; i+=s){
      for(j=0; j<n; j+=s){
        
        pi = (i/s)%p ; // Row number inside the processors grid.
        pj = (j/s)%q ; // Coulmn number inside the processors grid.
        dest = (pi*q + pj) ; // Rank of the processor where the data belongs.
        tag = 0; // Receiving the data from the leader processor.
        
        if(dest == myrank){
          /*Recv the data if this part belongs to me.*/
          for(k=i; k<i+s; k++){
            offset = k*n+j;
            MPI_Recv( A+offset, s, MPI_FLOAT, 0 , tag, MPI_COMM_WORLD, &status);
          }
        }
      }
    }
  }

/* TODO: Remove this before submitting.*/
/* Uncomment to see the block cyclic distribution of the matrix.*/
/*
  for (int i=0; i<nprocs; i++)
    {
        if (i == myrank) {
            std::cout << "Task " << myrank << std::endl ;
            printf("\nMy Matrix A:- \n");
            for (int j=0; j<m; j++) {
              printf("\t| ");
              for (int k=0; k<n; k++)
                printf("%.2f ", A[j*n+k]);
              printf("|\n");
            }
            printf("\n");
        }
        MPI_Barrier(MPI_COMM_WORLD);
  }
* /

  /*Perform l iterations of averaging*/
  while( (l--)>0 ){
    if(myrank == 0)
      printf("Iteration number %d\n",l+1);
  }
  
  if(myrank == 0){
    outputMatrix(outfile,A,m,n);
    printf("Output written to the file %s\n",argv[2]);
    /*Close all the files before exiting.*/
    if(fcloseall() != 0){
      printf("ERROR:- Problem while closing the files\n");
      free(A);
      MPI_Finalize();
      exit(FAILURE);
    }
  }  
  MPI_Finalize();
  free(A);
  return SUCCESS;
}

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
      printf("%.2f ", A[i*n+j]);
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
      printf("%.2f ", A[i*n+j]);
    printf("|\n");
  }
  printf("\n");
}
