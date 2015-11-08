
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

    int         my_rank;   /* My process rank           */

    int         p;         /* The number of processes   */

    int         norm;      /* The number of rows        */

    int         row;       /* Row number                */

    int         col;       /* Column number             */

    int         source;    /* Process sending matrices  */

    int         dest = 0;  /* All messages go to 0      */

    MPI_Status  status;

    void Get_data(int my_rank, int p);

    void Compute(int norm, int my_rank, int p);

    /* Let the system do what it needs to start up MPI */

    MPI_Init(&argc, &argv);

    /* Get my process rank */

    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    /* Find out how many processes are being used */

    MPI_Comm_size(MPI_COMM_WORLD, &p);

    /* Start Clock */
    printf("\nStarting clock.\n");
    gettimeofday(&etstart, &tzdummy);
    etstart2 = times(&cputstart);

    if(my_rank = 0){
        /* Process program parameters */
        parameters(argc, argv);

        /* Initialize A and B */
        initialize_inputs();

        /* Broadcast N */
        MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);
    } else {

        /* Getting N */
        MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);

        Get_data(my_rank, p);

        /* Gaussian elimination */
        for (norm = 0; norm < N - 1; norm++) {

            Inner_loop(norm, my_rank, p);

            MPI_Barrier(MPI_COMM_WORLD);
        }

    }

    if (my_rank == 0) {

        int row, col;

        /* Back substitution */

        for (row = N - 1; row >= 0; row--) {
            X[row] = B[row];
            for (col = N-1; col > row; col--) {
                X[row] -= A[row][col] * X[col];
            }
            X[row] /= A[row][row];
        }
    }

    MPI_Finalize();

    /* Stop Clock */
    gettimeofday(&etstop, &tzdummy);
    etstop2 = times(&cputstop);
    printf("Stopped clock.\n");
    usecstart = (unsigned long long)etstart.tv_sec * 1000000 + etstart.tv_usec;
    usecstop = (unsigned long long)etstop.tv_sec * 1000000 + etstop.tv_usec;

    printf("\nElapsed time = %g ms.\n",
    (float)(usecstop - usecstart)/(float)1000);

    printf("(CPU times are accurate to the nearest %g ms)\n",
    1.0/(float)CLOCKS_PER_SEC * 1000.0);
    printf("My total CPU time for parent = %g ms.\n",
    (float)( (cputstop.tms_utime + cputstop.tms_stime) -
    (cputstart.tms_utime + cputstart.tms_stime) ) /
    (float)CLOCKS_PER_SEC * 1000);
    printf("My system CPU time for parent = %g ms.\n",
    (float)(cputstop.tms_stime - cputstart.tms_stime) /
    (float)CLOCKS_PER_SEC * 1000);
    printf("My total CPU time for child processes = %g ms.\n",
    (float)( (cputstop.tms_cutime + cputstop.tms_cstime) -
    (cputstart.tms_cutime + cputstart.tms_cstime) ) /
    (float)CLOCKS_PER_SEC * 1000);
    /* Contrary to the man pages, this appears not to include the parent */
    printf("--------------------------------------------\n");



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
