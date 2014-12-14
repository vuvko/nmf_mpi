#ifndef __MATRIX_H__
#define __MATRIX_H__
/*
struct Matrix;
typedef struct Matrix Matrix;
*/
typedef struct Matrix
{
    double *data; // stored by rows
    unsigned rows;
    unsigned cols;
    char transposed;
} Matrix;

enum
{
    MAX_ELEMENTS = 1000000000
};

Matrix *make_matrix(double *data, int rows, int cols);
Matrix *make_matrix_zero(int rows, int cols);
Matrix *free_matrix(Matrix *matrix);
Matrix *copy_matrix(Matrix *origin);
Matrix *submatrix(Matrix *matrix, int offset_r, int offset_c, int rows, int cols);
Matrix *dot(Matrix *A, Matrix *B);
void elementwise(Matrix *A, Matrix *B, char op);
void add(Matrix *A, Matrix *B);
void subtract(Matrix *A, Matrix *B);
void prod(Matrix *A, Matrix *B);
void delim(Matrix *A, Matrix *B);
Matrix *transpose(Matrix *matrix);
Matrix *gramm_inner(Matrix *matrix);
Matrix *gramm_outer(Matrix *matrix);
double get(Matrix *matrix, unsigned row, unsigned col);
unsigned midx(Matrix *matrix, unsigned row, unsigned col);

void print_matrix(Matrix *matrix);

double loss(Matrix *V, Matrix *w, Matrix *H);

#endif
