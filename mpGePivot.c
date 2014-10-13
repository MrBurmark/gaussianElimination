/* 
AUTHOR: Rachel Beasley, Jason Burmark, Moses Lee
COMPILE: mpicc mpGePivot.c geutils.c -o mpGe -lm -O3
USAGE: mpirun -n [# number of nodes] ./mpGe [# number of nodes] [file of equation coefficients]

This file performs parallelized Gaussian Elimination with partial pivoting using MPI

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
    int stop;
    int num_rows = 0;
    double t1 = 0.0, t2 = 0.0, t3 = 0.0, t4 = 0.0, t_t1 = 0.0, t_t2 = 0.0;
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
        checkEqn = (double *) calloc(size * (size+1), sizeof(double));
        x = (double *) calloc(size, sizeof(double));

        readFile(eqn, size, fp);

        /* save copy of matrix for error checking */
        for (i=0; i<size * (size+1); i++) {
            checkEqn[i] = eqn[i];
        }
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
    stop = size/num_procs;
    if (my_rank < size%num_procs)
        stop++;

    for (i = 0; i < size; i++) {
        // if (0 == my_rank && i % 50 == 0)
        //     printf("Step %d\n", i);


        if (0 == my_rank)
            t_t1 = MPI_Wtime();



        /* find pivot row */
        row = i/num_procs;
        if (my_rank < i%num_procs)
            row++;

        if (row < stop) {
            my_maxRow = par_findMaxRow(my_eqn, row, i, num_rows, size+1);

            my_di.value = fabs(my_eqn[my_maxRow*(size+1)+i]);

        } else {
            my_maxRow = -1;
            my_di.value = -1.0;
        }

        my_di.row = my_maxRow*num_procs+my_rank;

        MPI_Allreduce(&my_di, &glob_di, 1, MPI_DOUBLE_INT, MPI_MAXLOC, MPI_COMM_WORLD);

        // need to copy to pivot and normalize row and swap with proper pivot position

        row = glob_di.row/num_procs;
        winnar = glob_di.row%num_procs;

        if (winnar == my_rank) {
            MPI_Bcast(my_eqn+row*(size+1) + i, size+1 - i, MPI_DOUBLE, winnar, MPI_COMM_WORLD);
            memcpy(pivot + i, my_eqn+row*(size+1) + i, (size+1 - i)*sizeof(double));
        } else {
            MPI_Bcast(pivot + i, size+1 - i, MPI_DOUBLE, winnar, MPI_COMM_WORLD);
        }


        if (glob_di.row != i) {
            if (my_rank == winnar && winnar == i%num_procs) {

                memcpy(my_eqn+row*(size+1) + i, my_eqn+i/num_procs*(size+1) + i, (size+1 - i)*sizeof(double));
                memcpy(my_eqn+i/num_procs*(size+1) + i, pivot + i, (size+1 - i)*sizeof(double));

            } else if (my_rank == winnar) {

                MPI_Recv(my_eqn+row*(size+1) + i, size+1 - i, MPI_DOUBLE, i%num_procs, tag, MPI_COMM_WORLD, &status);

            } else if (my_rank == i%num_procs) {

                MPI_Send(my_eqn+i/num_procs*(size+1) + i, size+1 - i, MPI_DOUBLE, winnar, tag, MPI_COMM_WORLD);
                memcpy(my_eqn+i/num_procs*(size+1) + i, pivot + i, (size+1 - i)*sizeof(double));

            }
        }

        if (0 == my_rank) {
            t_t2 = MPI_Wtime();
            t1 += t_t2 - t_t1;
        }


        /* process pivot row */
        if (my_rank == i%num_procs) {

            f = 1.0 / my_eqn[i/num_procs * (size+1) + i];

            for (j=i; j<size+1; j++) {
                my_eqn[i/num_procs * (size+1) + j] *= f;
            }

        }

        row = i/num_procs;
        if (my_rank <= i%num_procs)
            row++;

        /* row reduction using pivot row */
        for (j = row; j < stop; j++) {
            par_reduce(my_eqn+j*(size+1), pivot, i, size+1);
        }

        if (0 == my_rank) {
            t_t1 = MPI_Wtime();
            t2 += t_t1 - t_t2;
        }


    }

    if (0 == my_rank)
        t3 = MPI_Wtime();

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

    if (0 == my_rank)
        t3 = MPI_Wtime() - t3;

    /* print rows of row-reduced matrix */
    // if (0 == my_rank) 
    //     for (i=0; i<size; i++) 
    //         printRow(eqn, i, size);
    

    /* perform back substitution */

    if (0 == my_rank) {

        t4 = MPI_Wtime();

        backSub(x, eqn, size);

        t4 = MPI_Wtime() - t4;

        printf("nodes: %i\ntotal time, pivot, eliminate, gather, back-sub\n%.9lf\n%.9lf\n%.9lf\n%.9lf\n%.9lf\n", num_procs, t1+t2+t3+t4, t1, t2, t3, t4);

        dumpData(x, eqn, size);

        /* check solutions */
        ok = checkSoln(x, checkEqn, size);
        if (ok == 1)
            printf("All solutions are within error threshold\n");
        else
            printf("Some solutions are not within error threshold\n");
    }

    MPI_Finalize(); // program will die if not last
}
