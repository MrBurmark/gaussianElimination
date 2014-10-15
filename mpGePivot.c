/* 
AUTHOR: Rachel Beasley, Jason Burmark, Moses Lee
COMPILE: mpicc mpGePivot.c geutils.c -o mpGe -lm -O3
USAGE: mpirun -n [# number of nodes] ./mpGe [# number of nodes] [file of equation coefficients]

This file performs parallelized Gaussian Elimination with partial pivoting using MPI
This code normalizes the final row of the matrix

 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include "ge.h"

int main(int argc, char** argv) {

    FILE *fp;
    double *eqn, *my_eqn, *checkEqn, *pivot, *x;
    int size;
    int ok;
    int source;
    int my_row, win_row;
    int pivot_proc;
    int winnar;
    int num_rows = 0;
    double t0 = 0.0, t1 = 0.0, t2 = 0.0, t3 = 0.0, t4 = 0.0, t_t1 = 0.0, t_t2 = 0.0;

	int i, j, k;
    int *all_sizes;
    int *all_offsets;
	int my_rank;
    int my_num_rows;
    int num_procs;
    int tag = 0;
    MPI_Status status;
    struct Double_Int my_di, glob_di;

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
        ok = fscanf(fp, "%d", &size);

        printf("got size\n");
    }

    MPI_Bcast(&size, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if (num_procs > size){
        if (0 == my_rank)
            printf("Too many processes\n");
        exit(1);
    }

    my_eqn = (double *) calloc((size/num_procs+1) * (size+1), sizeof(double));
    pivot = (double *) calloc(size+1, sizeof(double));

    all_sizes = (int *) calloc(num_procs, sizeof(int));
    all_offsets = (int *) calloc(num_procs, sizeof(int));

    for (i=0, j=0, k=0; i < num_procs; i++){

        k = size/num_procs;

        if (i < size%num_procs) k++;

        all_sizes[i] = k*(size+1);
        all_offsets[i] = j*(size+1);
        
        j += k;
    }

    if (0 == my_rank) {
        eqn = (double *) calloc(size * (size+1), sizeof(double));
        checkEqn = (double *) calloc(size * (size+1), sizeof(double));
        x = (double *) calloc(size, sizeof(double));

        readFile(eqn, size, fp);

        printf("read file\n");

        /* save copy of matrix for error checking */
        // for (i=0; i<size * (size+1); i++) {
        //     checkEqn[i] = eqn[i];
        // }
    }

    /* Send Data to all processes */

    MPI_Scatterv(eqn, all_sizes, all_offsets, MPI_DOUBLE, my_eqn, all_sizes[my_rank], MPI_DOUBLE, 0, MPI_COMM_WORLD);

    t0 = MPI_Wtime();

    int pivotIndex = 0;
    double f;
    int my_maxRow = -1;
    my_num_rows = size/num_procs;
    if (my_rank < size%num_procs)
        my_num_rows++;

    for (i = 0; i < size; i++) {
         if (0 == my_rank && i % 50 == 0)
             printf("Step %d\n", i);


        if (0 == my_rank)
            t_t1 = MPI_Wtime();

        my_row = i/num_procs;
        pivot_proc = i%num_procs;
        if (my_rank < pivot_proc)
            my_row++;

        /* find pivot row */
        
        if (my_row < my_num_rows) {
            my_maxRow = par_findMaxRow(my_eqn, my_row, i, my_num_rows, size+1);

            my_di.value = fabs(my_eqn[my_maxRow*(size+1)+i]);

        } else {
            my_maxRow = -1;
            my_di.value = -1.0;
        }

        my_di.row = my_maxRow*num_procs+my_rank;

        MPI_Allreduce(&my_di, &glob_di, 1, MPI_DOUBLE_INT, MPI_MAXLOC, MPI_COMM_WORLD);

        // need to copy to pivot and normalize row and swap with proper pivot position

        win_row = glob_di.row/num_procs;
        winnar = glob_di.row%num_procs;

        if (winnar == my_rank) {
            MPI_Bcast(my_eqn+win_row*(size+1) + i, size+1 - i, MPI_DOUBLE, winnar, MPI_COMM_WORLD);
            memcpy(pivot + i, my_eqn+win_row*(size+1) + i, (size+1 - i)*sizeof(double));
        } else {
            MPI_Bcast(pivot + i, size+1 - i, MPI_DOUBLE, winnar, MPI_COMM_WORLD);
        }


        if (glob_di.row != i) {
            if (my_rank == winnar && winnar == pivot_proc) {

                memcpy(my_eqn+win_row*(size+1) + i, my_eqn+my_row*(size+1) + i, (size+1 - i)*sizeof(double));
                memcpy(my_eqn+my_row*(size+1) + i, pivot + i, (size+1 - i)*sizeof(double));

            } else if (my_rank == winnar) {

                MPI_Recv(my_eqn+win_row*(size+1) + i, size+1 - i, MPI_DOUBLE, pivot_proc, tag, MPI_COMM_WORLD, &status);

            } else if (my_rank == pivot_proc) {

                MPI_Send(my_eqn+my_row*(size+1) + i, size+1 - i, MPI_DOUBLE, winnar, tag, MPI_COMM_WORLD);
                memcpy(my_eqn+my_row*(size+1) + i, pivot + i, (size+1 - i)*sizeof(double));

            }
        }

        if (0 == my_rank) {
            t_t2 = MPI_Wtime();
            t1 += t_t2 - t_t1;
        }


        /* process pivot row */
        if (my_rank == pivot_proc) {

            f = 1.0 / my_eqn[my_row * (size+1) + i];

            for (j=i; j<size+1; j++) {
                my_eqn[my_row * (size+1) + j] *= f;
            }

            my_row++;
        }

        /* row reduction using pivot row */
        for (j = my_row; j < my_num_rows; j++) {
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

    MPI_Gatherv(my_eqn, all_sizes[my_rank], MPI_DOUBLE, checkEqn, all_sizes, all_offsets, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // reorder matrix rows
    if (0 == my_rank){
        for (i = 0; i < size; i++) {

            source = i%num_procs;

            memcpy(eqn+i*(size+1), checkEqn+all_offsets[source], (size+1)*sizeof(double));

            all_offsets[source] += size+1;
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
        t0 = MPI_Wtime() - t0;

        printf("nodes: %i\tsize %d\ttotal time %.9lf\ntotal time, pivot, eliminate, gather, back-sub\n%.9lf\n%.9lf\n%.9lf\n%.9lf\n%.9lf\n", num_procs, size, t0, t1+t2+t3+t4, t1, t2, t3, t4);

        dumpData(x, checkEqn, size);

        /* read in original equation */
        fclose(fp);
        fp = fopen(argv[2], "r");
        ok = fscanf(fp, "%d", &size);
        readFile(checkEqn, size, fp);

        /* read in original equation */
        fclose(fp);
        fp = fopen(argv[2], "r");
        ok = fscanf(fp, "%d", &size);
        readFile(checkEqn, size, fp);

        /* check solutions */
        ok = checkSoln(x, checkEqn, size);
        if (ok == 1)
            printf("All solutions are within error threshold\n");
        else
            printf("Some solutions are not within error threshold\n");
    }

    MPI_Finalize(); // program will die if not last
}
