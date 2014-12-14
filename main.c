#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <mpi.h>

#include "common/common.h"
#include "common/matrix.h"
#include "common/parsecfg.h"
#include "common/random.h"

int
main(int argc, char *argv[])
{
    //int rank, nproc;
    
    //MPI_Init(&argc, &argv);
    //MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    //MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    //printf("Hello world from %d of %d processes\n", rank, nproc);
    //MPI_Finalize();
    
    printf("Reading confguration file...\n");
    
    ConfigFile *cfg = config_read("config.txt");
    printf("Getting parameters.\n");
    unsigned rows, cols, topics, max_iter;
    if (!config_get_int(cfg, "rows", &rows)) {
        rows = 10;
    }
    if (!config_get_int(cfg, "cols", &cols)) {
        cols = 10;
    }
    if (!config_get_int(cfg, "topics", &topics)) {
        topics = 10;
    }
    if (!config_get_int(cfg, "max_iter", &max_iter)) {
        max_iter = 10;
    }
    printf("Creating matrix V (%d x %d).\n", rows, cols);
    Matrix *V = make_matrix_zero(rows, cols);
    unsigned i, j;
    for (i = 0; i < V->cols; ++i) {
        for (j = 0; j < V->rows; ++j) {
            if (i == j) {
                V->data[i * V->rows + j] = 10;
            } else {
                V->data[i * V->rows + j] = 1;
            }
        }
    }
    //print_matrix(V);
    printf("Creating matrix W (%d x %d).\n", rows, topics);
    Matrix *W = make_matrix_zero(rows, topics);
    for (i = 0; i < W->cols; ++i) {
        for (j = 0; j < W->rows; ++j) {
            W->data[i * W->rows + j] = (i + j + 1.0) / (rows * cols);
        }
    }
    //print_matrix(W);
    printf("Creating matrix H (%d x %d).\n", topics, cols);
    Matrix *H = make_matrix_zero(topics, cols);
    for (i = 0; i < H->cols; ++i) {
        for (j = 0; j < H->rows; ++j) {
            H->data[i * H->rows + j] = (H->rows * H->cols - (double)(i + j))
                                       / (rows * cols);
        }
    }
    //print_matrix(H);
    
    printf("initial loss = %.5f\n", loss(V, W, H));
    
    unsigned iter;
    for (iter = 0; iter < max_iter; ++iter) {
        printf("Iteration %d:\n", iter + 1);
        printf("  calculating new W.\n");
        Matrix *Ht = transpose(H);
        Matrix *WH = dot(W, H);
        Matrix *VHt = dot(V, Ht);
        Matrix *WHHt = dot(WH, Ht);
        Ht = free_matrix(Ht);
        WH = free_matrix(WH);
        delim(VHt, WHHt);
        WHHt = free_matrix(WHHt);
        prod(W, VHt); // new W here
        VHt = free_matrix(VHt);
        Ht = free_matrix(Ht);
        printf("  calculating new H.\n");
        Matrix *Wt = transpose(W);
        WH = dot(W, H);
        Matrix *WtWH = dot(Wt, WH);
        WH = free_matrix(WH);
        Matrix *WtV = dot(Wt, V);
        Wt = free_matrix(Wt);
        delim(WtV, WtWH);
        WtWH = free_matrix(WtWH);
        prod(H, WtV); // new H here
        WtV = free_matrix(WtV);
        printf("  new loss = %.5f\n", loss(V, W, H));
    }
    V = free_matrix(V);
    W = free_matrix(W);
    H = free_matrix(H);
    cfg = config_free(cfg);
    return 0;
}
