#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/***
    readFile()

    read file of matrix coefficients 
    data[][]: 2D array to store coefficients (output)
    n: number of rows (input)
    fp: pointer to FILE struct for input file (input)

****/

void readFile(double *data, int n, FILE *fp) {
  int i;
  for (i=0; i< n * (n+1); i++) {
    fscanf(fp, "%lf", &data[i]);
  }
}

/***
    printRow()

    print row r of mat[][]

    mat[][]: nR x nR+1 matrix (input)
    r: row to print
    nR: number of rows

***/

void printRow(double *mat, int r, int nR) {
  int j;

  for (j=0; j<nR+1; j++) {
    printf("%f ", mat[r * (nR+1) + j]);
  }
  printf("\n");
}

void fprintRow(FILE *fp, double *mat, int r, int nR) {
  int j;

  for (j=0; j<nR+1; j++) {
    fprintf(fp, "%f ", mat[r * (nR+1) + j]);
  }
  fprintf(fp, "\n");
}


/***
    reduce()

    row-reduce a row, based on pivot row, from startIndex to endIndex inclusive

    myRow: row to be reduced (input/output)
    pivotRow: pivot row (input)
    startIndex: index of 1st element in row to be reduced (input)
    endIndex: index of last element in row to be reduced  (input)

****/

void reduce(double *myRow, double *pivotRow, int startIndex, int endIndex) {
  int i;
  float factor;
  float oldMr;

  factor = -myRow[startIndex];
  myRow[startIndex] = 0;
  
  for (i=startIndex+1; i<=endIndex; i++) {
    oldMr = myRow[i];
    myRow[i] = myRow[i] + pivotRow[i] * factor;
  }
}

void par_reduce(double *myRow, double *pivot, int start, int end) {
  int i;
  float factor;

  factor = -myRow[start]/pivot[start];
  myRow[start] = 0.0;
  
  for (i = start+1; i < end; i++) {
    myRow[i] = myRow[i] + pivot[i] * factor;
  }
}

/***
    backSub()

    perform back substitution on nR x nR+1 system of equations

    s: array of nR solutions (output)
    mat: upper triangular matrix (nR x nR+1) of coefficients (input)
    nR: number of rows in mat

***/

void backSub(double *s, double *mat, int nR) {
  int i, j;
  double sum = 0;

  for (i=nR-1; i>=0; i--) {
    sum = 0;
    for (j=nR-1; j>i; j--) {
      sum += mat[i*(nR+1) + j] * s[j];
    }
    s[i] = (mat[i*(nR+1) + nR] - sum) / mat[i*(nR+1) + i];
  }
}

/***
    checkSoln()

    check solutions by substituting into original system of equations

    s: array of nR solutions (input)
    mat: matrix of nR x nR+1 coefficients
    nR: number of rows in mat

***/
int checkSoln(double *s, double *mat, int nR) {
  int i, j;
  int ok = 1;
  double sum;
  double diff;
  
  for (i=0; i<nR; i++) {
    sum = 0;
    for (j=0; j<nR; j++) {
      sum += s[j] * mat[i * (nR+1) + j];
    }
    diff = sum - mat[i * (nR+1) + nR];
    if (diff >= .005 || diff < -.005) {
      ok = 0;
      printf("Mismatch at row %d %f %f %f\n", i, diff, sum,  mat[i * (nR+1) + nR]);
    }
  }
  return ok;
}

int checkEqn(double *mat0, double *mat1, int nR) {
  int i, j;
  double diff;
  int ok=1;

  for (i=0; i<nR; i++) {
    for (j=0; j<nR+1; j++) {
      diff = fabs(mat0[i * (nR+1) + j] - mat1[i * (nR+1) + j]);
      if (diff >= .005) {
	ok = 0;
	printf("Mismatch at %d %d %f %f %f\n", i, j, diff,
	       mat0[i * (nR+1) + j], mat1[i * (nR+1) + j]);
      }
    }
  }
  
}

/***
    findMaxRow()

    find row i with max(eqn[][pivotIndex]), return i

    eqn: nR x nR+1 matrix of coefficients (input)
    pivotIndex: index of column for max search (input)
    nR: number of rows in eqn[][] (input)

***/

int findMaxRow(double *eqn, int pivotIndex, int nR) {
  int i, maxI;
  double max;
  max = fabs(eqn[pivotIndex * (nR+1) + pivotIndex]);
  maxI = pivotIndex;
  for (i=pivotIndex+1; i<nR; i++) {
    if (fabs(eqn[i*(nR+1) + pivotIndex]) > max) {
      max = fabs(eqn[i*(nR+1) + pivotIndex]);
      maxI = i;
    }
  }
  return maxI;
  
}

int par_findMaxRow(double *my_eqn, int pivotRow, int pivotCol, int nR, int nC) {
  int i, maxI;
  double max;
  max = fabs(my_eqn[pivotRow * nC + pivotCol]);
  maxI = pivotRow;
  for (i=pivotRow+1; i<nR; i++) {
    if (fabs(my_eqn[i*nC + pivotCol]) > max) {
      max = fabs(my_eqn[i*nC + pivotCol]);
      maxI = i;
    }
  }
  return maxI;
}

/***
    swapRows()

    swap two rows in mat[][]

    mat: nR x nR+1 matrix of coefficients (input/output)
    maxRow: index of one row (input)
    pivotIndex: index of other row (input)
    nR: number of rows in mat[][] (input)

***/

void swapRows(double *mat, int maxRow, int pivotIndex, int nR) {
  double temp;
  int i;

  for (i=0; i<nR+1; i++) {
    temp = mat[maxRow * (nR+1) + i];
    mat[maxRow * (nR+1) + i] = mat[pivotIndex * (nR+1) + i];
    mat[pivotIndex * (nR+1) + i] = temp;
  }
}

/***
    dump file of row-reduced system and solutions 
***/
void dumpData(double *s, double *mat, int nR) {
  FILE *fp;
  int i, j;
  fp = fopen("dump.out", "w");
  fprintf(fp, "%d\n", nR);
  for (i=0; i<nR; i++) {
    fprintf(fp, "%d ", i);
    for (j=0; j<nR+1; j++) {
      fprintf(fp, "%f ", mat[i * (nR+1) + j]);
    }
    fprintf(fp, "\n");
  }
  for (i=0; i<nR; i++) {
    fprintf(fp, "%d %f\n", i, s[i]);
  }
  fclose(fp);
}
