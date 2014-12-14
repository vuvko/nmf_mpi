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
    
    printf("Reading confguration file...\n");
    
    ConfigFile *cfg = config_read("config.txt");
    printf("Getting parameters.\n");
    int topics, max_iter, num_threads, add_col = 0, add_row = 0;
    unsigned rows, cols, div_row, div_col, rows_per_rank = 0, cols_per_rank = 0;
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
    }
    omp_set_num_threads(num_threads);
    if (rank == 0) {
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
        MPI_Bcast(&topics, 1, MPI_INT, 0, MPI_COMM_WORLD);
        unsigned nnz;
        fscanf(Vin, "%u%u%u", &cols, &rows, &nnz);
        rows_per_rank = rows;
        add_row = 0, add_col = 0;
        if (div_row > 0) {
            rows_per_rank /=  div_row;
            add_row = rows % div_row;
        }
        cols_per_rank = cols;
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
                rows_per_rank = rows;
                if (div_row > 0) {
                    rows_per_rank /=  div_row;
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
                            nnz_read++;
                        }
                        j -= col_shift + 1;
                        i -= row_shift + 1;
                        if (i >= rows_per_rank || j >= cols_per_rank) {
                            j += col_shift + 1;
                            i += row_shift + 1;
                            read_val = 0;
                            break;
                        }
                        read_val = 1;
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
    } else {
        // Receving parameters
        MPI_Bcast(&div_row, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
        printf("[%d]: Received div_row = %d\n", rank, div_row);
        MPI_Bcast(&div_col, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
        printf("[%d]: Received div_col = %d\n", rank, div_col);
        MPI_Bcast(&topics, 1, MPI_INT, 0, MPI_COMM_WORLD);
        // Receiving V
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
    MPI_Comm w_comm = 0;
    MPI_Comm h_comm = 0;
    int w_rank = rank % div_row;
    int h_rank = rank / div_row;
    MPI_Comm_split(MPI_COMM_WORLD, w_rank, h_rank, &w_comm);
    MPI_Comm_split(MPI_COMM_WORLD, h_rank, w_rank, &h_comm);
    printf("[%d] my new ranks (%d %d)\n", rank, w_rank, h_rank);
    
    if (rank == 0) {
        rows_per_rank = rows;
        add_row = 0, add_col = 0;
        if (div_row > 0) {
            rows_per_rank /=  div_row;
            add_row = rows % div_row;
        }
        if (add_row) {
            rows_per_rank += 1;
        }
        cols_per_rank = cols;
        if (div_col > 0) {
            cols_per_rank /=  div_col;
            add_col = cols % div_col;
        }
        if (add_col) {
            cols_per_rank += 1;
        }
        unsigned row_mat, col_mat;
        char first = 1;
        // Generating and sending Ws
        for (row_mat = 0; row_mat < div_row; ++row_mat) {
            Matrix *Wo = make_matrix_random(rows_per_rank, topics, rnd);
            // Sending W from here
            unsigned size = Wo->rows * Wo->cols;
            printf("Sending W:\n");
            for (col_mat = 0; col_mat < div_col; ++col_mat) {
                if (first) {
                    W = copy_matrix(Wo);
                    first = 0;
                } else {
                    int send_idx = col_mat * div_row + row_mat;
                    MPI_Send(&Wo->data[0], size, MPI_DOUBLE, 
                        send_idx, 0, MPI_COMM_WORLD);
                }
            }
            //MPI_Bcast(&Wo->data[0], size, MPI_DOUBLE, 0, w_comm);
            Wo = free_matrix(Wo);
            if (add_row) {
                add_row -= 1;
                if (!add_row) {
                    rows_per_rank -= 1;
                }
            }
        }
    } else {
        double *data = (double *)calloc(rows * topics, sizeof(data[0]));
        //MPI_Bcast(&data[0], rows * topics, MPI_DOUBLE, 0, w_comm);
        MPI_Recv(&data[0], rows * topics, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);
        W = make_matrix(data, rows, topics);
        free(data);
        printf("[%d]: Received W of size (%d x %d).\n", rank, rows, topics);
    }
    if (rank == 0) {
        unsigned row_mat, col_mat;
        char first = 1;
        // Generating and sending Hs
        for (col_mat = 0; col_mat < div_col; ++col_mat) {
            Matrix *Ho = make_matrix_random(topics, cols_per_rank, rnd);
            // Sending H from here
            unsigned size = Ho->rows * Ho->cols;
            printf("Sending H:\n");
            for (row_mat = 0; row_mat < div_row; ++row_mat) {
                if (first) {
                    H = copy_matrix(Ho);
                    first = 0;
                } else {
                    int send_idx = col_mat * div_row + row_mat;
                    MPI_Send(&Ho->data[0], size, MPI_DOUBLE, 
                        send_idx, 0, MPI_COMM_WORLD);
                }
            }
            //MPI_Bcast(&Ho->data[0], size, MPI_DOUBLE, 0, h_comm);
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
        rows = V->rows;
        cols = V->cols;
    } else {
        double *data = (double *)calloc(topics * cols, sizeof(data[0]));
        //MPI_Bcast(&data[0], topics * cols, MPI_DOUBLE, 0, h_comm);
        MPI_Recv(&data[0], topics * cols, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);
        H = make_matrix(data, topics, cols);
        free(data);
        printf("[%d]: Received H of size (%d x %d).\n", rank, topics, cols);
    }
    // Sending loss
    double cur_loss = loss(V, W, H);
    unsigned idx, row, col;
    if (rank == 0) {
        double other_loss = 0;
        for (idx = 1; idx < nproc; ++idx) {
            MPI_Recv(&other_loss, 1, MPI_DOUBLE, idx, 0, MPI_COMM_WORLD, &status);
            cur_loss += other_loss;
        }
        printf("initial loss = %.5f\n", cur_loss);
    } else {
        MPI_Send(&cur_loss, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    }
    unsigned iter;
    double init_criteria = 0, tol = 1e-5;
    char stop = 0;
    for (iter = 0; iter < max_iter && !stop; ++iter) {
        if (rank == 0) {
            printf("Iteration %u:\n", iter + 1);
            printf("  calculating new W.\n");
        }
        Matrix *gradW = NULL, *gradH = NULL;
        Matrix *WH = dot(W, H);
        Matrix *WHHt = dott(WH, H);
        // Sending and receving matricies WHHt here
        for (col = 0; col < div_col; ++col) {
            //printf("[%d]: casting %d 1\n", h_rank, col);
            Matrix *WHHt_other = make_matrix_zero(rows, topics);
            MPI_Bcast(&WHHt_other->data[0], rows * topics, MPI_DOUBLE, col, w_comm);
            add(WHHt, WHHt_other);
            WHHt_other = free_matrix(WHHt_other);
        }
        WH = free_matrix(WH);
        Matrix *VHt = dott(V, H);
        // Sending and receving matricies VHt here
        for (col = 0; col < div_col; ++col) {
            //printf("[%d]: casting %d 2\n", rank, col);
            Matrix *VHt_other = make_matrix_zero(rows, topics);
            MPI_Bcast(&VHt_other->data[0], rows * topics, MPI_DOUBLE, col, w_comm);
            add(VHt, VHt_other);
            VHt_other = free_matrix(VHt_other);
        }
        gradW = copy_matrix(WHHt);
        subtract(gradW, VHt);
        delim(VHt, WHHt);
        WHHt = free_matrix(WHHt);
        prod(W, VHt); // new W here
        VHt = free_matrix(VHt);
        if (rank == 0) {
            printf("  calculating new H.\n");
        }
        WH = dot(W, H);
        Matrix *WtWH = tdot(W, WH);
        // Sending and receving matricies WtWH here
        for (row = 0; row < div_row; ++row) {
            //printf("[%d]: casting %d 3\n", rank, row);
            Matrix *WtWH_other = make_matrix_zero(topics, cols);
            MPI_Bcast(&WtWH->data[0], topics * cols, MPI_DOUBLE, row, h_comm);
            add(WtWH, WtWH_other);
            WtWH_other = free_matrix(WtWH_other);
        }
        WH = free_matrix(WH);
        Matrix *WtV = tdot(W, V);
        // Sending and receving matricies WtV here
        for (row = 0; row < div_row; ++row) {
            //printf("[%d]: casting %d 4\n", rank, row);
            Matrix *WtV_other = make_matrix_zero(topics, cols);
            MPI_Bcast(&WtV->data[0], topics * cols, MPI_DOUBLE, row, h_comm);
            add(WtV, WtV_other);
            WtV_other = free_matrix(WtV_other);
        }
        gradH = copy_matrix(WtWH);
        subtract(gradH, WtV);
        delim(WtV, WtWH);
        WtWH = free_matrix(WtWH);
        prod(H, WtV); // new H here
        WtV = free_matrix(WtV);
        // Loss function
        cur_loss = loss(V, W, H);
        if (rank == 0) {
            double other_loss = 0;
            for (idx = 1; idx < nproc; ++idx) {
                MPI_Recv(&other_loss, 1, MPI_DOUBLE, idx, 0, MPI_COMM_WORLD, &status);
                cur_loss += other_loss;
            }
            printf("  new loss = %.5f\n", cur_loss);
        } else {
            MPI_Send(&cur_loss, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
        }
        // Stopping criteria
        double criteria = 0;
        if (iter == 0) {
            init_criteria = stopping(gradW, W) + stopping(gradH, H);
        } else {
            criteria = stopping(gradW, W) + stopping(gradH, H);
            if (rank == 0) {
                double other_criteria = 0;
                for (idx = 1; idx < nproc; ++idx) {
                    MPI_Recv(&other_criteria, 1, MPI_DOUBLE, idx, 0, MPI_COMM_WORLD, &status);
                    criteria += other_criteria;
                }
            } else {
                MPI_Send(&criteria, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
            }
        }
        gradW = free_matrix(gradW);
        gradH = free_matrix(gradH);
        if (iter > 0) {
            if (rank == 0) {
                printf("  Criteria = %.5lf\n", criteria);
                if (criteria < tol * init_criteria) {
                    printf("Stopping by meeting criteria.\n");
                    stop = 1;
                }
                for (idx = 1; idx < nproc; ++idx) {
                    MPI_Send(&stop, 1, MPI_BYTE, idx, 0, MPI_COMM_WORLD);
                }
            } else {
                MPI_Recv(&stop, 1, MPI_BYTE, 0, 0, MPI_COMM_WORLD, &status);
            }
        }
    }
    MPI_Comm_free(&w_comm);
    MPI_Comm_free(&h_comm);
final:
    V = free_matrix(V);
    W = free_matrix(W);
    H = free_matrix(H);
    MPI_Finalize();
    cfg = config_free(cfg);
    return 0;
}
