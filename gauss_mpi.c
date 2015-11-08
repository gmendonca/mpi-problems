
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
int N = 2000;  /* Matrix size */

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

    int         my_rank;   /* My process rank           */

    int         p;         /* The number of processes   */

    int         norm;      /* The number of rows        */

    int         row;       /* Row number                */

    int         col;       /* Column number             */

    int         source;    /* Process sending matrices  */

    int         dest = 0;  /* All messages go to 0      */

    MPI_Status  status;

    void Get_data(int my_rank, int p);

    void Inner_loop(int norm, int my_rank, int p);

    /* Let the system do what it needs to start up MPI */

    MPI_Init(&argc, &argv);

    /* Get my process rank */

    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    /* Find out how many processes are being used */

    MPI_Comm_size(MPI_COMM_WORLD, &p);

    double startTime, stopTime;

    if(my_rank = 0) {
        /* Start Clock */
        printf("\nStarting clock.\n");
        startTime = MPI_Wtime();

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

            int i;
            for (i = norm; i < N; i++) {
                MPI_Bcast(A[i], N, MPI_FLOAT, i%p, MPI_COMM_WORLD);
                MPI_Bcast(&B[i], 1, MPI_FLOAT, i%p, MPI_COMM_WORLD);
            }

            Inner_loop(norm, my_rank, p);

            MPI_Barrier(MPI_COMM_WORLD);
        }

    }

    if (my_rank == 0) {

        MPI_Bcast(A[N-1], N, MPI_FLOAT, (N-1)%p, MPI_COMM_WORLD);
        MPI_Bcast(&B[N-1], 1, MPI_FLOAT, (N-1)%p, MPI_COMM_WORLD);

        /* Back substitution */

        for (row = N - 1; row >= 0; row--) {
            X[row] = B[row];
            for (col = N-1; col > row; col--) {
                X[row] -= A[row][col] * X[col];
            }
            X[row] /= A[row][row];
        }

        /* Stop Clock */
        stopTime = MPI_Wtime();

        printf("\nElapsed time = %lf ms.\n",(stopTime - startTime));
        printf("--------------------------------------------\n");
    }

    MPI_Finalize();

    exit(0);

}

void Get_data(

    int     my_rank  /* in  */,

    int     p        /* in  */) {


        int source = 0;    /* All local variables used by */

        int dest;          /* MPI_Send and MPI_Recv       */

        MPI_Status status;

        int i;

        if (my_rank == 0){
            for(i = 1; i < N; i++) {
                dest = i%p;
                if (dest != 0) {
                    MPI_Send(A[i], N - 1, MPI_FLOAT, dest, i, MPI_COMM_WORLD);
                    MPI_Send(&B[i], 1, MPI_FLOAT, dest, i + N - 1, MPI_COMM_WORLD);
                }
            }
        } else {
            for(i = 0; i < N; i++) {
                dest = i%p;
                if (dest == my_rank) {
                    MPI_Recv(A[i], N - 1, MPI_FLOAT, source, i, MPI_COMM_WORLD, &status);
                    MPI_Recv(&B[i], 1, MPI_FLOAT, source, i + N - 1, MPI_COMM_WORLD, &status);
                }
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
            for (row = norm + 1; row < N; row++) {
                if (row % p == my_rank) {
                    multiplier = A[row][norm] / A[norm][norm];
                    for (col = norm; col < N; col++) {
                        A[row][col] -= A[norm][col] * multiplier;
                    }
                    B[row] -= B[norm] * multiplier;
                }
            }
        }
