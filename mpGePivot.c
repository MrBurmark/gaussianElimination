/* 
AUTHOR: 
COMPILE: mpicc mpGePivot.c geutils.c -o mpGe -lm -O3
USAGE: mpirun -n [# number of nodes] ./mpGe [# number of nodes] [file of equation coefficients]
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include "ge.h"

int main(int argc, char** argv) {

    FILE *fp;
    double *eqn, *my_eqn, *checkEqn;
    double *x;
    int size;
    int ok;
    int step;
    int source, dest;

	int i, j;
	int my_rank;
    int num_procs;
    int tag = 0;
    MPI_Status status;

    MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

    if (argc != 3) {
        if (0 == my_rank)
            printf("Usage: mpirun -n [# number of nodes] ./mpGe [# number of nodes]  [file of equation coefficients]\n");
        exit(1);
    }

    /* read file of equations */
    if (0 == my_rank) {
        fp = fopen(argv[2], "r");
        fscanf(fp, "%d", &size);
    }

    MPI_Bcast(&size, 1, MPI_INT, 0, MPI_COMM_WORLD);

    my_eqn = (double *) calloc((size/num_procs+1) * (size+1), sizeof(double));

    if (0 == my_rank) {
        eqn = (double *) calloc(size * (size+1), sizeof(double));
        //checkEqn = (double *) calloc(size * (size+1), sizeof(double));
        x = (double *) calloc(size, sizeof(double));

        readFile(eqn, size, fp);

        /* save copy of matrix for error checking */
        //for (i=0; i<size * (size+1); i++) {
        //    checkEqn[i] = eqn[i];
        //}
    }

    /* Send Data to all processes */
    for (i = 0; i < size; i++) {

        dest = i%num_procs;

        if (dest == my_rank) {

            step = (i/num_procs)*(size+1);

            if (0 == my_rank) {

                memcpy(my_eqn+step,eqn+i*(size+1),(size+1)*sizeof(double));

            } else {

                MPI_Recv(my_eqn+step, size+1, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD, &status);
            }

        } else if (0 == my_rank){

            MPI_Send(eqn+i*(size+1), size+1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
        }
    }



    // int pivotIndex = 0;
    // double f;
    // int maxRow = -1;

    // for (step = 0; step < size-1; step++) {
    //     if (0 == my_rank && step % 50 == 0)
    //         printf("Step %d\n", step);
    //     /* find pivot row */
    //     maxRow = findMaxRow(eqn, step, size);

    //     /* pre-process pivot row */
    //     f = eqn[maxRow * (size+1) + step];

    //     for (j=step; j<size+1; j++) {
    //         eqn[maxRow * (size+1) + j] /= f;
    //     }

    //     /* swap current row with pivot row */
    //     swapRows(eqn, maxRow, step, size);    

    //     pivotIndex = step;

    //     /*
    //     printf("Step %d New pivot: ", step);
    //     printRow(eqn, pivotIndex, size);
    //     */

    //     /* row reduction using pivot row */
    //     for (i=pivotIndex+1; i<size; i++) {
    //         reduce(eqn+i*(size+1), eqn+pivotIndex*(size+1), pivotIndex, size);
    //     }
    // }

    /* print first 23 rows of row-reduced matrix */
    /*
    int endRow = size;
    if (endRow > 23)
    endRow = 23;

    for (i=0; i<endRow; i++) {
    printRow(eqn, i, size);
    }
    */

    /* Retrieve Data from all processes */
    for (i = 0; i < size; i++) {

        source = i%num_procs;

        if (source == my_rank) {

            step = (i/num_procs)*(size+1);

            if (0 == my_rank) {

                memcpy(eqn+i*(size+1), my_eqn+step, (size+1)*sizeof(double));

            } else {

                MPI_Send(my_eqn+step, size+1, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD);
            }

        } else if (0 == my_rank){

            MPI_Recv(eqn+i*(size+1), size+1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
        }
    }

    /* perform back substitution */

    if (0 == my_rank) {
        backSub(x, eqn, size);
        dumpData(x, eqn, size);
    }

    /* check solutions */
    /*
    ok = checkSoln(x, checkEqn, size);
    if (ok == 1)
    printf("All solutions are within error threshold\n");
    else
    printf("Some solutions are not within error threshold\n");
    */

    MPI_Finalize(); // program will die if not last
}
