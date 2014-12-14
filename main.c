#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include <mpi.h>
#include <omp.h>

#include "common/common.h"
#include "common/matrix.h"
#include "common/parsecfg.h"
#include "common/random.h"

int
main(int argc, char *argv[])
{
    int rank, nproc;
    //printf("Hello world from %d of %d processes\n", rank, nproc);
    
    printf("Reading confguration file...\n");
    
    ConfigFile *cfg = config_read("config.txt");
    printf("Getting parameters.\n");
    int topics, max_iter, num_threads;
    unsigned div_row, div_col;
    const char *dataset = config_get(cfg, "dataset");
    if (!config_get_int(cfg, "topics", &topics)) {
        topics = 10;
    }
    if (!config_get_int(cfg, "max_iter", &max_iter)) {
        max_iter = 10;
    }
    if (!config_get_int(cfg, "num_threads", &num_threads)) {
        num_threads = 8;
    }
    printf("Initializing MPI.\n");
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Status status;
    Matrix *V = NULL, *W = NULL, *H = NULL;
    Random *rnd = NULL;
    if (rank == 0) {
        printf("Initializing OpenMP with %d threads.\n", num_threads);
        omp_set_num_threads(num_threads);
        printf("Initializing random number generator\n");
        rnd = random_create(cfg);
        if (rnd) {
            srand(rnd->seed);
        }
        printf("Reading matrix V from %s\n", dataset);
        FILE *Vin = NULL;
        if (!(Vin = fopen(dataset, "r"))) {
            fprintf(stderr, "Error opening file %s\n", dataset);
            goto final;
        }
        //V = read_matrix_uci(Vin);
        long pwr = sqrtl(nproc);
        div_row = pwr;
        div_col = pwr;
        if (pwr * pwr > nproc) {
            div_row -= 1;
        } else if (pwr * pwr < nproc) {
            div_col += 1;
        }
        // Sending parameters
        printf("Broadcasting parameters\n");
        MPI_Bcast(&div_row, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
        MPI_Bcast(&div_col, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
        unsigned rows, cols, nnz;
        fscanf(Vin, "%u%u%u", &cols, &rows, &nnz);
        unsigned rows_per_rank = rows;
        int add_row = 0, add_col = 0;
        if (div_row > 0) {
            rows_per_rank /=  div_row;
            add_row = rows % div_row;
        }
        unsigned cols_per_rank = cols;
        if (div_col > 0) {
            cols_per_rank /=  div_col;
            add_col = cols % div_col;
        }
        if (add_col) {
            cols_per_rank += 1;
        }
        int col_mat, row_mat, col_shift = 0, row_shift = 0;
        int col, row;
        double val;
        int i, j, nnz_read = 0;
        char read_val = 1;
        int rec_min = 0, rec_max = div_row, rec_rank = 0;
        for (col_mat = 0; col_mat < div_col; ++col_mat) {
            Matrix **Vs = (Matrix **)calloc(div_row, sizeof(Vs[0]));
            for (col = 0; col < cols_per_rank; ++col) {
                row_shift = 0;
                if (div_row > 0) {
                    add_row = rows % div_row;
                }
                if (add_row) {
                    rows_per_rank += 1;
                }
                for (row_mat = 0; row_mat < div_row; ++row_mat) {
                    if (col == 0) {
                        Vs[row_mat] = make_matrix_zero(rows_per_rank, cols_per_rank);
                    }
                    for (row = 0; row < rows_per_rank; ++row) {
                        if (read_val) {
                            if (nnz_read >= nnz) {
                                break;
                            }
                            fscanf(Vin, "%d%d%lf", &j, &i, &val);
                            printf("Reading [%d %d] -> %lf\n", i, j, val);
                            nnz_read++;
                        }
                        j -= col_shift + 1;
                        i -= row_shift + 1;
                        printf("new %d %d from %d %d (%d %d)\n", i, j, row_shift, col_shift, rows_per_rank, cols_per_rank);
                        if (i >= rows_per_rank || j >= cols_per_rank) {
                            printf("Nope, waiting...\n");
                            j += col_shift + 1;
                            i += row_shift + 1;
                            read_val = 0;
                            break;
                        }
                        read_val = 1;
                        printf("Storing to %d mat.\n", row_mat);
                        Vs[row_mat]->data[midx(Vs[row_mat], i, j)] = val;
                    }
                    row_shift += rows_per_rank;
                    if (add_row) {
                        add_row -= 1;
                        if (!add_row) {
                            rows_per_rank -= 1;
                        }
                    }
                }
            }
            if (col_mat == 0) {
                V = copy_matrix(Vs[0]);
                // Send matrices from here
                printf("Sending matricies V from %d to %d.\n", 1, rec_max);
                for (rec_rank = 1; rec_rank < rec_max; ++rec_rank) {
                    MPI_Send(&Vs[rec_rank]->rows, 1, MPI_UNSIGNED, rec_rank, 
                        0, MPI_COMM_WORLD);
                    MPI_Send(&Vs[rec_rank]->cols, 1, MPI_UNSIGNED, rec_rank, 
                        0, MPI_COMM_WORLD);
                    unsigned size = Vs[rec_rank]->rows * Vs[rec_rank]->cols;
                    MPI_Send(&Vs[rec_rank]->data[0], size, MPI_DOUBLE, rec_rank, 
                        0, MPI_COMM_WORLD);
                }
            } else {
                // Send matrices from here
                printf("Sending matricies V from %d to %d.\n", rec_min, rec_max);
                unsigned idx;
                for (idx = 0; idx < rec_max - rec_min; ++idx) {
                    rec_rank = rec_min + idx;
                    MPI_Send(&Vs[idx]->rows, 1, MPI_UNSIGNED, rec_rank, 
                        0, MPI_COMM_WORLD);
                    MPI_Send(&Vs[idx]->cols, 1, MPI_UNSIGNED, rec_rank, 
                        0, MPI_COMM_WORLD);
                    unsigned size = Vs[idx]->rows * Vs[idx]->cols;
                    MPI_Send(&Vs[idx]->data[0], size, MPI_DOUBLE, rec_rank, 
                        0, MPI_COMM_WORLD);
                }
            }
            for (row_mat = 0; row_mat < div_row; ++row_mat) {
                Vs[row_mat] = free_matrix(Vs[row_mat]);
            }
            free(Vs);
            Vs = NULL;
            col_shift += cols_per_rank;
            if (add_col) {
                add_col -= 1;
                if (!add_col) {
                    cols_per_rank -= 1;
                }
            }
            rec_min += div_row;
            rec_max += div_row;
        }
        fclose(Vin);
        Vin = NULL;
        
        if (div_row > 0) {
            add_row = rows % div_row;
        }
        if (add_row) {
            rows_per_rank += 1;
        }
        if (div_col > 0) {
            add_col = cols % div_col;
        }
        if (add_col) {
            cols_per_rank += 1;
        }
        //int send_rank = 1;
        // Generating and sending Ws
        W = make_matrix_random(rows_per_rank, topics, rnd);
        // Sending W from here
        for (row_mat = 1; row_mat < div_row; ++row_mat) {
            Matrix *Wo = make_matrix_random(rows_per_rank, topics, rnd);
            // Sending W from here
            Wo = free_matrix(Wo);
            if (add_row) {
                add_row -= 1;
                if (!add_row) {
                    rows_per_rank -= 1;
                }
            }
        }
        H = make_matrix_random(topics, cols_per_rank, rnd);
        // Sending H from here
        for (col_mat = 1; col_mat < div_col; ++col_mat) {
            Matrix *Ho = make_matrix_random(topics, cols_per_rank, rnd);
            // Sending H from here
            Ho = free_matrix(Ho);
            if (add_col) {
                add_col -= 1;
                if (!add_col) {
                    cols_per_rank -= 1;
                }
            }
        }
        if (rnd) {
            rnd = rnd->ops->free(rnd);
        }
        printf("initial loss = %.5f\n", loss(V, W, H));
    } else {
        printf("[%d]: Hello world!\n", rank);
        // Receving parameters
        MPI_Bcast(&div_row, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
        printf("[%d]: Received div_row = %d\n", rank, div_row);
        MPI_Bcast(&div_col, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
        printf("[%d]: Received div_col = %d\n", rank, div_col);
        // Receiving V
        unsigned rows = 0, cols = 0;
        MPI_Recv(&rows, 1, MPI_UNSIGNED, 0, 0, MPI_COMM_WORLD, &status);
        MPI_Recv(&cols, 1, MPI_UNSIGNED, 0, 0, MPI_COMM_WORLD, &status);
        double *data = (double *)calloc(rows * cols, sizeof(data[0]));
        MPI_Recv(&data[0], rows * cols, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);
        V = make_matrix(data, rows, cols);
        free(data);
        data = NULL;
        printf("[%d]: Received V of size (%d x %d).\n", rank, rows, cols);
    }
    if (rank == 0) {
        printf("Initializing MPI gorups.\n");
    }
    //int w_group[div_row], h_group[div_col];
    
    unsigned iter;
    for (iter = 0; iter < max_iter; ++iter) {
        if (rank == 0) {
            printf("Iteration %u:\n", iter + 1);
            printf("  calculating new W.\n");
        }
        Matrix *WH = dot(W, H);
        Matrix *WHHt = dott(WH, H);
        WH = free_matrix(WH);
        Matrix *VHt = dott(V, H);
        // Sending and receving matricies VHt here
        delim(VHt, WHHt);
        WHHt = free_matrix(WHHt);
        prod(W, VHt); // new W here
        VHt = free_matrix(VHt);
        if (rank == 0) {
            printf("  calculating new H.\n");
        }
        WH = dot(W, H);
        Matrix *WtWH = tdot(W, WH);
        WH = free_matrix(WH);
        Matrix *WtV = tdot(W, V);
        // Sending and receving matricies WtV here
        delim(WtV, WtWH);
        WtWH = free_matrix(WtWH);
        prod(H, WtV); // new H here
        WtV = free_matrix(WtV);
        if (rank == 0) {
            printf("  new loss = %.5f\n", loss(V, W, H));
        }
    }
final:
    V = free_matrix(V);
    W = free_matrix(W);
    H = free_matrix(H);
    MPI_Finalize();
    cfg = config_free(cfg);
    return 0;
}
