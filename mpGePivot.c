/* 
AUTHOR: 
COMPILE: mpicc mpGePivot.c geutils.c -o mpGe -lm -O3
USAGE: mpirun -n [# number of nodes] ./mpGe [# number of nodes] [file of equation coefficients]
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include "ge.h"

int main(int argc, char** argv) {

    FILE *fp;
    double *eqn, *checkEqn;
    double *x;
    int numRows;
    int ok;
    int step;

	int i, j;
	int my_rank;
    int num_procs;
    int tag = 0;
    MPI_Status status;

    if (argc != 3) {
    	if (0 == my_rank)
    		printf("Usage: mpirun -n [# number of nodes] ./mpGe [# number of nodes]  [file of equation coefficients]\n");
    	exit(1);
    }

    MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);



    /* read file of equations */
    fp = fopen(argv[1], "r");
    fscanf(fp, "%d", &numRows);
    eqn = (double *) calloc(numRows * (numRows+1), sizeof(double));
    checkEqn = (double *) calloc(numRows * (numRows+1), sizeof(double));
    x = (double *) calloc(numRows, sizeof(double));

    if (0 == my_rank)
        readFile(eqn, numRows, fp);



    /* save copy of matrix for error checking */
    for (i=0; i<numRows * (numRows+1); i++) {
        checkEqn[i] = eqn[i];
    }

    int pivotIndex = 0;
    double f;
    int maxRow = -1;

    for (step = 0; step < numRows-1; step++) {
        if (step % 50 == 0)
            printf("Step %d\n", step);
        /* find pivot row */
        maxRow = findMaxRow(eqn, step, numRows);

        /* pre-process pivot row */
        f = eqn[maxRow * (numRows+1) + step];

        for (j=step; j<numRows+1; j++) {
            eqn[maxRow * (numRows+1) + j] /= f;
        }

        /* swap current row with pivot row */
        swapRows(eqn, maxRow, step, numRows);    

        pivotIndex = step;

        /*
        printf("Step %d New pivot: ", step);
        printRow(eqn, pivotIndex, numRows);
        */

        /* row reduction using pivot row */
        for (i=pivotIndex+1; i<numRows; i++) {
            reduce(eqn+i*(numRows+1), eqn+pivotIndex*(numRows+1), pivotIndex, numRows);
        }
    }

    /* print first 23 rows of row-reduced matrix */
    /*
    int endRow = numRows;
    if (endRow > 23)
    endRow = 23;

    for (i=0; i<endRow; i++) {
    printRow(eqn, i, numRows);
    }
    */

    /* perform back substitution */
    backSub(x, eqn, numRows);


    MPI_Finalize();

    dumpData(x, eqn, numRows);

    /* check solutions */
    /*
    ok = checkSoln(x, checkEqn, numRows);
    if (ok == 1)
    printf("All solutions are within error threshold\n");
    else
    printf("Some solutions are not within error threshold\n");
    */

}