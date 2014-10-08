/****

     File: geMatch.c

     Compile: gcc geMatch.c -o geMatch -lm -O3
     Usage: ./geMatch [file of equation coefficients & solutions]
                      [file of equation coefficients & solutions]

****/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(int argc, char **argv) {

  FILE *fp, *fp1;
  double *eqn, *checkEqn;
  double *x;
  int i, j;
  int numRows, numRows1;
  int ok;
  int inNum, inNum1;
  double inDouble, inDouble1;

  if (argc != 3) {
    printf("Usage: ./geMatchk [coefficients+solutions file] [coefficients/solutions file]\n");
    return(0);
  }

  /* read file of equations */
  fp = fopen(argv[1], "r");
  fscanf(fp, "%d", &numRows);

  fp1 = fopen(argv[2], "r");
  fscanf(fp1, "%d", &numRows1);

  if (numRows != numRows1) {
    printf("number of rows mismatched: %d %d\n", numRows, numRows1);
    return(0);
  }

  printf("Matching upper triangular matrices\n");
  for (i=0; i<numRows; i++) {
    fscanf(fp, "%d", &inNum);
    fscanf(fp1, "%d", &inNum1);
    if (inNum != inNum1) {
      printf("row index mismatched: %d %d\n", inNum, inNum1);
      return(0);
    }
    for (j=0; j<numRows+1; j++) {
      fscanf(fp, "%lf", &inDouble);
      fscanf(fp1, "%lf", &inDouble1);
      if (fabs(inDouble - inDouble1) > .0001) {
	printf("coefficient %d %d mismatched: %f %f\n", i, j, inDouble, inDouble1);
      }
    }
  }

  printf("Matching solutions\n");
  for (i=0; i<numRows; i++) {
    fscanf(fp, "%d", &inNum);
    fscanf(fp1, "%d", &inNum1);
    if (i != inNum || inNum != inNum1) {
      printf("Solution index mismatched %d %d %d\n", i, inNum, inNum1);
    }
    fscanf(fp, "%lf", &inDouble);
    fscanf(fp1, "%lf", &inDouble1);
    if (fabs(inDouble - inDouble1) > .005) {
      printf("solution %d mismatched: %f %f\n", i, inDouble, inDouble1);
    }
  }
    
}
