#ifndef _GE_H
#define _GE_H

void readFile(double *data, int n, FILE *fp);
void printRow(double *mat, int r, int nR);
void fprintRow(FILE *fp, double *mat, int r, int nR);
void reduce(double *myRow, double *pivotRow, int startIndex, int endIndex);
void backSub(double *s, double *mat, int nR);
int checkSoln(double *s, double *mat, int nR);
int checkEqn(double *mat0, double *mat1, int nR);
int findMaxRow(double *eqn, int pivotIndex, int nR);
int par_findMaxRow(double *my_eqn, int pivotRow, int pivotCol, int nR, int nC);
void swapRows(double *mat, int maxRow, int pivotIndex, int nR);
void dumpData(double *s, double *mat, int nR);

struct Double_Int{
double value;
int row;
};

#endif

