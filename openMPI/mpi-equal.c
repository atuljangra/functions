/* This is more of a hello Mpi type of program.
 * I allocate and initialize array A, 
 * and then I store row wise sum in array B using normal data distribution.
 * Caution: Untested code.
 */
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

#define TAG 1024

int main(int argc, char *argv[]) {
  float **A, **B, **C, *tmp;
  float startTime, endTime;
  int numElements, offset, stripSize, myrank, numnodes, N, i, j, k;
  
  MPI_Init(&argc, &argv);
  
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Comm_size(MPI_COMM_WORLD, &numnodes);
  printf("%d,  %d \n", numnodes, myrank);
  N = atoi(argv[1]);
  
  // allocate A 
  // Master gets whole
  if (myrank == 0) {
    tmp = (float *) malloc (sizeof(float ) * N * N);
    A = (float **) malloc (sizeof(float *) * N);
    for (i = 0; i < N; i++)
      A[i] = &tmp[i * N];
  }
  else {
    tmp = (float *) malloc (sizeof(float ) * N * N / numnodes);
    A = (float **) malloc (sizeof(float *) * N / numnodes);
    for (i = 0; i < N / numnodes; i++)
      A[i] = &tmp[i * N];
  }
  
  // Allocate C
  if (myrank == 0) {
    tmp = (float *) malloc (sizeof(float ) * N * N);
    C = (float **) malloc (sizeof(float *) * N);
    for (i = 0; i < N; i++)
      C[i] = &tmp[i * N];
  }
  else {
    tmp = (float *) malloc (sizeof(float ) * N * N / numnodes);
    C = (float **) malloc (sizeof(float *) * N / numnodes);
    for (i = 0; i < N / numnodes; i++)
      C[i] = &tmp[i * N];
  }

  if (myrank == 0) {
    // initialize A
    for (i=0; i<N; i++) {
      for (j=0; j<N; j++) {
        A[i][j] = N;
        }
    }
  }
  
  // start timer
  if (myrank == 0) {
    startTime = MPI_Wtime();
  }
  
  stripSize = N/numnodes;

// Data distribution of A.
  if (myrank == 0) {
    offset = stripSize;
    numElements = stripSize * N;
    for (i=1; i<numnodes; i++) {
      MPI_Send(A[offset], numElements, MPI_FLOAT, i, TAG, MPI_COMM_WORLD);
      offset += stripSize;
    }
  }
  else {  // receive my part of A
    MPI_Recv(A[0], stripSize * N, MPI_FLOAT, 0, TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }
  
  // everyone initializes C to 0
  for (i=0; i<stripSize; i++) {
    for (j=0; j<N; j++) {
      C[i][j] = 0.0;
    }
  }

  // Actual Computation.
  for (i=0; i<stripSize; i++) {
    for (j=0; j<N; j++) {
	    C[i][j] = A[i][j];
      }
    }

  // Gathering results from workers.
  if (myrank == 0) {
    offset = stripSize; 
    numElements = stripSize * N;
    for (i=1; i<numnodes; i++) {
      MPI_Recv(C[offset], numElements, MPI_FLOAT, i, TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      offset += stripSize;
    }
  }
  else {
    MPI_Send(C[0], stripSize * N, MPI_FLOAT, 0, TAG, MPI_COMM_WORLD);
  }


  // stop timer
  if (myrank == 0) {
    endTime = MPI_Wtime();
    printf("Time is %f\n", endTime-startTime);
  }
  
  // printing results
  if (myrank == 0) {
    for (i=0; i<N; i++) {
      for (j=0; j<N; j++) {
        printf("%.2f ", C[i][j]);
      }
      printf("\n");
    }
  }
  
  MPI_Finalize();
  return 0;
}


