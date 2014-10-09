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
    double *eqn, *my_eqn, *checkEqn, *pivot;
    double *x;
    int size;
    int ok;
    int step;
    int source, dest;
    int row;
    int winnar;
    int num_rows = 0;
    struct Double_Int my_di, glob_di;

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
    pivot = (double *) calloc(size+1, sizeof(double));

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

            num_rows++;
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



    int pivotIndex = 0;
    double f;
    int my_maxRow = -1;

    for (i = 0; i < size-1; i++) {
        if (0 == my_rank && i % 50 == 0)
            printf("Step %d\n", i);

        /* find pivot row */
        row = i/num_procs;
        if (my_rank < i%num_procs)
            row++;
        my_maxRow = par_findMaxRow(my_eqn, row, i, num_rows, size+1);

        my_di.value = fabs(my_eqn[my_maxRow*(size+1)+i]);
        my_di.row = my_maxRow;

        MPI_Allreduce(&my_di, &glob_di, 1, MPI_DOUBLE_INT, MPI_MAXLOC, MPI_COMM_WORLD);




        // need to copy to pivot and normalize row and swap with proper pivot position


        row = glob_di.row/num_procs;
        winnar = glob_di.row%num_procs;

        if (winnar == my_rank)
            MPI_Bcast(my_eqn+row*(size+1), size+1, MPI_DOUBLE, winnar, MPI_COMM_WORLD);
        else
            MPI_Bcast(pivot, size+1, MPI_DOUBLE, winnar, MPI_COMM_WORLD);





        /* pre-process pivot row */
        f = eqn[my_maxRow * (size+1) + i];

        for (j=i; j<size+1; j++) {
            eqn[my_maxRow * (size+1) + j] /= f;
        }

        /* swap current row with pivot row */
        swapRows(eqn, my_maxRow, i, size);

        pivotIndex = i;

        /*
        printf("Step %d New pivot: ", i);
        printRow(eqn, pivotIndex, size);
        */

        /* row reduction using pivot row */
        for (j=pivotIndex+1; j<size; j++) {
            reduce(eqn+j*(size+1), eqn+pivotIndex*(size+1), pivotIndex, size);
        }
    }

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
