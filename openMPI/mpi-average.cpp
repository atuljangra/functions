/*
 MPI program to perform the template 
 (averaging) based computations on a matrix.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <vector>
#include <mpi.h>
#include <math.h>
#include <unordered_map>

#define N_ARGS 14 // Number of arguments required by the program.
#define FAILURE 1
#define SUCCESS 0
#define DELM " "

/* Different tags for different phases of the program
   used while sending the data from one to another.*/
#define DISTRIBUTE 1
#define COLLECT 2
using namespace std;
void readMatrix(FILE* infile, double *A, int m, int n);
void outputMatrix(FILE* outfile, double *A, int m, int n);

int returnHashKey(int i, int j) {
  return (i*3037 + j*4973);
}

struct data {
  int i;
  int j;
  double value;
};

int main( int argc, char *argv[]){

  FILE *infile = NULL;
  FILE* outfile = NULL;
  int m, n, p, q, s, l, ns, ms;
  double a, b, c, d, e;
  int myrank, nprocs;
  double *A = NULL;
  int i, j, pi, pj, k, offset, tag, dest, *owner;
  int temp_pi, temp_pj, temp_dest;
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
    
  if( !(m>0 && n>0 && p>0 && q>0 && s>0 && l>0 && (m%s)==0 && (n%s)==0 ) ){
    printf("ERROR:- Invalid Input to the program\n");
    MPI_Finalize();
    exit(FAILURE);
  }

  A = (double *)calloc(m*n, sizeof(double));
  owner = (int *)calloc((m/s)*(n/s), sizeof(int));
  if( A==NULL || owner == NULL ){
    printf("ERROR:- Malloc failed on the processor %d\n",myrank);
    free(A);
    free(owner);
    MPI_Finalize();
    exit(FAILURE);
  }
  
  ms = m/s;
  ns = n/s;
  
  /* Initialize the lookup table to see who owns which block. */
  for(i=0; i<ms; i++){
    for(j=0; j<ns; j++){
      pi = i % p ; // Row number inside the processors grid.
      pj = j % q ; // Coulmn number inside the processors grid.
      owner[i*ns+j] = pi*q + pj;
    }
  }
  
  // TODO: Remove this print.
  if(myrank==0) {
    printf("\nOwners matrix is:- \n");
    for (i=0; i<ms; i++) {
      printf("\t| ");
      for (j=0; j<ns; j++)
        printf("%d ", owner[i*ns+j]);
      printf("|\n");
    }
    printf("\n");
  }

  /***************************** COMPUTATION **************************/
  /*Perform l iterations of averaging*/
  while( (l--)>0 ){
    if(myrank == 0)
      printf("Iteration number %d\n",l+1);

  /**************************** DATA DISTRIBUTION *********************/
  if(myrank == 0)
    readMatrix(infile,A,m,n);

  /*Iterate over all the matrix block by block and distribute the data.*/
  for(i=0; i<ms; i++){
    for(j=0; j<ns; j++){
            
      /************************ Leader thread**************************/
      if(owner[i*ns+j] != myrank && myrank == 0){
        /*Send the data if this part does not belongs to me.*/
        for(k=i*s; k<(i+1)*s; k++){
          offset = k*n+j*s;
          tag = DISTRIBUTE;
          MPI_Send( A+offset, s, MPI_DOUBLE, owner[i*ns+j], tag, MPI_COMM_WORLD);
        }
      }
      /************************ Worker thread**************************/ 
      if(owner[i*ns+j] == myrank && myrank != 0){
        /*Recv the data if this part belongs to me.*/
        for(k=i*s; k<(i+1)*s; k++){
          offset = k*n+j*s;
          tag = DISTRIBUTE;
          MPI_Recv( A+offset, s, MPI_DOUBLE, 0 , tag, MPI_COMM_WORLD, &status);
        }
      }
    }
  }
  
  MPI_Bcast (A, m*n, MPI_DOUBLE, 0, MPI_COMM_WORLD); 
      
    //~ for(i=0; i<ms; i++){
      //~ for(j=0; j<ns; j++){
        //~ 
        //~ pi = i % p ; // Row number inside the processors grid.
        //~ pj = j % q ; // Coulmn number inside the processors grid.
        //~ 
        //~ /* Recv data from left processor, Send to the right.*/
        //~ if(owner[i*ns+j] == myrank){
          //~ 
          //~ MPI_Request request;
          //~ 
          //~ temp_pi = pi;
          //~ if( j == 0 ){ // leftmost processor
            //~ temp_pj = (j+1) % q ;
            //~ temp_dest = owner[i*ns+j+1]; //(temp_pi*q + temp_pj);
            //~ 
            //~ printf("\t1 Task %d, i %d,j %d temp-dest %d SENDING\n",myrank,i,j,temp_dest);
            //~ for(k=i*s; k<i*s+s; k++){
              //~ offset = k*n+(j+1)*s-1;
              //~ MPI_Send( A+offset, 1, MPI_DOUBLE, temp_dest, tag, MPI_COMM_WORLD);
            //~ }
            //~ printf("\tTask %d 1 sent succesfully\n",myrank);
          //~ }
          //~ else if( j == ns-1 ){ //rightmost processor
            //~ temp_pj = (j-1)%q ;
            //~ temp_dest = owner[i*ns+j-1]; //(temp_pi*q + temp_pj);
            //~ printf("\t4 Task %d, i %d,j %d temp-dest %d RECVING\n",myrank,i,j,temp_dest);
            //~ for(k=i; k<i*s+s; k++){
              //~ offset = k*n+j*s-1;
              //~ MPI_Recv( A+offset, 1, MPI_DOUBLE, temp_dest , tag, MPI_COMM_WORLD, &status);
            //~ }
            //~ printf("\tTask %d 4 recv succesfully\n",myrank);
          //~ }
          //~ else  {
           //~ temp_pj = (j-1)%q ;
            //~ temp_dest = owner[i*ns+j-1]; //(temp_pi*q + temp_pj);
            //~ printf("\t2 Task %d, i %d,j %d temp-dest %d RECVING\n",myrank,i,j,temp_dest);
            //~ for(k=i*s; k<i*s+s; k++){
              //~ offset = k*n+j*s-1;
              //~ MPI_Recv( A+offset, 1, MPI_DOUBLE, temp_dest , tag, MPI_COMM_WORLD, &status );
            //~ }
            //~ printf("\tTask %d 2 recv succesfully\n",myrank);
        //~ 
            //~ temp_pj = (j+1)%q ;
            //~ temp_dest = owner[i*ns+j+1]; //(temp_pi*q + temp_pj);
            //~ printf("\t3 Task %d, i %d,j %d temp-dest %d Sending\n",myrank,i,j,temp_dest);
            //~ for(k=i*s; k<i*s+s; k++){
              //~ offset = k*n+j*s+s-1;
              //~ MPI_Send( A+offset, 1, MPI_DOUBLE, temp_dest, tag, MPI_COMM_WORLD);
            //~ }
            //~ printf("\tTask %d 3 sent succesfully\n",myrank);
         //~ }
         //~ 
        //~ 
          //~ /* Recv data from the right processor, Send to the left one.*/
          //~ temp_pi = pi;
        //~ 
          //~ if( j+s > n-1 ){ //rightmost processor
            //~ temp_pj = ((j-s)/s)%q ;
            //~ temp_dest = (temp_pi*q + temp_pj);
            //~ 
            //~ for(k=i; k<i+s; k++){
            //~ offset = k*n+j;
              //~ MPI_Send( A+offset, 1, MPI_DOUBLE, temp_dest , tag, MPI_COMM_WORLD);
            //~ }
          //~ }
          //~ else if( j-s < 0 ){ // leftmost processor
            //~ temp_pj = ((j+s)/s) % q ;
            //~ temp_dest = (temp_pi*q + temp_pj);
            //~ 
            //~ for(k=i; k<i+s; k++){
              //~ offset = k*n+j+s;
              //~ MPI_Recv( A+offset, 1, MPI_DOUBLE, temp_dest, tag, MPI_COMM_WORLD, &status);
            //~ }
          //~ }
          //~ else  {
//~ 
            //~ temp_pj = ((j-s)/s)%q ;
            //~ temp_dest = (temp_pi*q + temp_pj);
            //~ 
            //~ for(k=i; k<i+s; k++){
              //~ offset = k*n+j;
              //~ MPI_Send( A+offset, 1, MPI_DOUBLE, temp_dest, tag, MPI_COMM_WORLD);
            //~ }
            //~ 
            //~ temp_pj = ((j+s)/s)%q ;
            //~ temp_dest = (temp_pi*q + temp_pj);
            //~ 
            //~ for(k=i; k<i+s; k++){
              //~ offset = k*n+j+s;
              //~ MPI_Recv( A+offset, 1, MPI_DOUBLE, temp_dest , tag, MPI_COMM_WORLD, &status );
            //~ }
          //~ }
        //~ }
      //~ }
    //~ }
    //~ 
    //~ for(j=0; j<n; j+=s){
      //~ for(i=0; i<m; i+=s){
//~ 
        //~ pi = (i/s) % p ; // Row number inside the processors grid.
        //~ pj = (j/s) % q ; // Coulmn number inside the processors grid.
        //~ owner = (pi*q + pj) ; // Rank of the processor where the data belongs.
        //~ 
        //~ if(myrank == owner){
          //~ 
          //~ /* Recv data from the upper processor, Send to the lower one.*/
          //~ temp_pj = pj;
          //~ 
          //~ if( i-s < 0 ){ // uppermost processor
            //~ temp_pi = ((i+s)/s) % q ;
            //~ temp_dest = (temp_pi*q + temp_pj);
            //~ printf("\t1 Task %d, i %d,j %d Dest %d, temp-dest %d SENDING\n",myrank,i,j,dest,temp_dest);
            //~ offset = (i+s-1)*n+j;
            //~ MPI_Ssend( A+offset, s, MPI_DOUBLE, temp_dest, tag, MPI_COMM_WORLD);
            //~ printf(" 1 sent succesfully\n");
          //~ }
          //~ else if( i+s > m-1 ){ //lowermost processor
            //~ temp_pi = ((i-s)/s)%q ;
            //~ temp_dest = (temp_pi*q + temp_pj);
//~ 
            //~ printf("\t4 Task %d, i %d,j %d Dest %d, temp-dest %d RECVING\n",myrank,i,j,dest,temp_dest);
            //~ offset = (i-1)*n+j;
            //~ tag = COMPUTE;
            //~ MPI_Recv( A+offset, s, MPI_DOUBLE, temp_dest , tag, MPI_COMM_WORLD, &status);
            //~ printf(" 4 recv succesfully\n");
          //~ }
          //~ else  {
            //~ temp_pi = ((i-s)/s)%q ;
            //~ temp_dest = (temp_pi*q + temp_pj);
            //~ 
            //~ printf("\t2 Task %d, i %d,j %d Dest %d, temp-dest %d RECVING\n",myrank,i,j,dest,temp_dest);
            //~ offset = (i-1)*n+j;
            //~ MPI_Recv( A+offset, s, MPI_DOUBLE, temp_dest , tag, MPI_COMM_WORLD, &status );
            //~ printf("2 recv succesfully\n");
            //~ 
            //~ temp_pi = ((i+s)/s)%q ;
            //~ temp_dest = (temp_pi*q + temp_pj);
            //~ 
            //~ printf("\t3 Task %d, i %d,j %d Dest %d, temp-dest %d RECVING\n",myrank,i,j,dest,temp_dest);
            //~ offset = (i+s-1)*n+j;
            //~ MPI_Ssend( A+offset, s, MPI_DOUBLE, temp_dest, tag, MPI_COMM_WORLD);
            //~ printf("3 sent succesfully\n");
          //~ }
          //~ 
        //~ }
      //~ }
    //~ }
    /*Perform the  computation here.*/


    unordered_map <int, data> values;
    /*
     * Insert all changes values into the hash table and after carrying out
     * all the computation, copy values from hash table to the array.
     */ 
    double numer;
    int key;
    double denom;
    for (int iter = 0; iter < m; iter++) {
      for (int jter = 0; jter < n; jter++) {
        pi = (iter/s) % p ; // Row number inside the processors grid.
        pj = (jter/s) % q ; // Coulmn number inside the processors grid.
        int my_owner = (pi*q + pj) ; // Rank of the processor where the data belongs.
        if (my_owner == myrank) {
            numer= c*A[iter*n + jter];
            denom= c;
            numer+= (iter > 0) ? a*A[(iter-1)*n + jter] : 0;
            denom+= (iter > 0)? a : 0;
            numer+= (iter + 1 < m) ? e*A[(iter+1)*n + jter] : 0;
            denom+= (iter + 1 < m) ? e : 0;
            numer+= (jter > 0) ? b*A[iter*n + jter-1] : 0;
            denom+= (jter > 0) ? b : 0;
            numer+= (jter + 1 < n) ? d*A[(iter)*n + jter + 1 ] : 0;
            denom+= (jter + 1 < n) ?d : 0;
            key = returnHashKey(iter, jter);
            //~ cout << "HASH: inserting " << key << " " <<
            //~ " " << iter << " " << jter << " " << numer/denom << endl;
            if (values.find(key) != values.end()) {
              cout << "HASH: PROBLEM\n";
              exit(-1);
            }
            data * newData = new data;
            newData -> i = iter;
            newData -> j = jter;
            newData -> value = numer/denom;
            values[key] = *newData;
//~             A[iter*n + jter] = numer/denom;
        }
        else {
          // Debug
        }
      }
    }

    double valueToPut;
    int iter;
    int jter;
    for (auto it = values.begin(); it != values.end(); it++) {
        iter = it -> second.i;
        jter = it -> second.j;
        valueToPut = it -> second.value;
        //~ cout << "HASH: Adding " << " " <<
            //~ " " << iter << " " << jter << " " <<valueToPut << endl;
        A[iter*n + jter] = valueToPut;
    }

    // Erase the hashmap.
    values.clear();
    assert(values.empty());
    
    
/* TODO: Remove this before submitting.*/
/* Uncomment to see the block cyclic distribution of the matrix.*/
  for (i=0; i<nprocs; i++)
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

  /**************************** DATA COLLECTION ***********************/

  /*Iterate over all the matrix block by block and collect the data.*/
  for(i=0; i<ms; i++){
    for(j=0; j<ns; j++){
      
      /************************ Leader thread**************************/
      if(owner[i*ns+j] != myrank && myrank == 0){
        /*Recv the data if this part does not belongs to me.*/
        for(k=i*s; k<i*s+s; k++){
          offset = k*n+j*s;
          tag = COLLECT;
          MPI_Recv( A+offset, s, MPI_DOUBLE, owner[i*ns+j] , tag, MPI_COMM_WORLD, &status);
        }
      }
      /************************ Worker thread**************************/ 
      if(owner[i*ns+j] == myrank && myrank != 0){
        /*Send the data if this part belongs to me.*/
        for(k=i*s; k<i*s+s; k++){
          offset = k*n+j*s;
          tag = COLLECT;
          MPI_Send( A+offset, s, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD);
        }
      }
    }
  }
 }
 
  /* Print the results and close all the files. */
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

void readMatrix(FILE* infile, double *A, int m, int n){
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

void outputMatrix(FILE* outfile, double *A, int m, int n){
  int i,j;
  
  for(i=0; i<m; i++){
    for(j=0; j<n; j++){
      fprintf(outfile,"%.2f",A[i*n+j]);
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
