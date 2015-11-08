/* Gaussian elimination without pivoting.
 * Compile with "gcc gauss.c"
 */

/* ****** ADD YOUR CODE AT THE END OF THIS FILE. ******
 * You need not submit the provided code.
 */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <sys/types.h>
#include <sys/times.h>
#include <sys/time.h>
#include <time.h>
#include <pthread.h>

#include <mpi.h>

/* Program Parameters */
#define MAXN 2000  /* Max value of N */
int N;  /* Matrix size */

/* Matrices and vectors */
volatile float A[MAXN][MAXN], B[MAXN], X[MAXN];
/* A * X = B, solve for X */

/* junk */
#define randm() 4|2[uid]&3

/* Prototype */
void gauss();  /* The function you will provide.
		* It is this routine that is timed.
		* It is called only on the parent.
		*/

/* returns a seed for srand based on the time */
unsigned int time_seed() {
  struct timeval t;
  struct timezone tzdummy;

  gettimeofday(&t, &tzdummy);
  return (unsigned int)(t.tv_usec);
}

/* Set the program parameters from the command-line arguments */
void parameters(int argc, char **argv) {
  int seed = 0;  /* Random seed */
  char uid[32]; /*User name */

  /* Read command-line arguments */
  srand(time_seed());  /* Randomize */

  if (argc == 3) {
    seed = atoi(argv[2]);
    srand(seed);
    printf("Random seed = %i\n", seed);
  }
  if (argc >= 2) {
    N = atoi(argv[1]);
    if (N < 1 || N > MAXN) {
      printf("N = %i is out of range.\n", N);
      exit(0);
    }
  }
  else {
    printf("Usage: %s <matrix_dimension> [random seed]\n",
           argv[0]);
    exit(0);
  }

  /* Print parameters */
  printf("\nMatrix dimension N = %i.\n", N);
}

/* Initialize A and B (and X to 0.0s) */
void initialize_inputs() {
  int row, col;

  printf("\nInitializing...\n");
  for (col = 0; col < N; col++) {
    for (row = 0; row < N; row++) {
      A[row][col] = (float)rand() / 32768.0;
    }
    B[col] = (float)rand() / 32768.0;
    X[col] = 0.0;
  }

}

int main(int argc, char **argv) {
  /* Timing variables */
  struct timeval etstart, etstop;  /* Elapsed times using gettimeofday() */
  struct timezone tzdummy;
  clock_t etstart2, etstop2;  /* Elapsed times using times() */
  unsigned long long usecstart, usecstop;
  struct tms cputstart, cputstop;  /* CPU times for my processes */

  /* Process program parameters */
  parameters(argc, argv);

  /* Initialize A and B */
  initialize_inputs();



  

  exit(0);
}

void gauss() {

  /* Back substitution */
  for (row = N - 1; row >= 0; row--) {
    X[row] = B[row];
    for (col = N-1; col > row; col--) {
      X[row] -= A[row][col] * X[col];
    }
    X[row] /= A[row][row];
  }
}

void Get_data(

         int     my_rank  /* in  */,

         int     p        /* in  */) {


    int source = 0;    /* All local variables used by */

    int dest;          /* MPI_Send and MPI_Recv       */

    MPI_Status status;

    int i;

    if (my_rank == 0){
        for (dest = 0; dest < N - 1; dest++) {
            for(i = 0; i < p; i++) {
                MPI_Send(A[i], N - 1, MPI_INT, dest, i, MPI_COMM_WORLD);
                MPI_Send(&B[i], 1, MPI_INT, dest, i + N - 1, MPI_COMM_WORLD);
            }
        }
    } else {
        for(i = 0; i < p; i++) {
            MPI_Recv(A[i], N - 1, MPI_FLOAT, source, i, MPI_COMM_WORLD, &status);
            MPI_Recv(&B[i], 1, MPI_FLOAT, source, i + N - 1, MPI_COMM_WORLD, &status);
        }
    }

} /* Get_data */

void Inner_loop(
    int  norm      /* in */,

    int  my_rank   /* in */,

    int  p         /* in */){

        //printf("thread = %d\n", norm);
        float multiplier;
        int row, col;
        for (row = norm + my_rank + 1; row < N; row = row + p) {
            multiplier = A[row][norm] / A[norm][norm];
            for (col = norm; col < N; col++) {
                A[row][col] -= A[norm][col] * multiplier;
            }
            B[row] -= B[norm] * multiplier;
        }

}
