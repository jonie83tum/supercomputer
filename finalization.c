/**
 * Finalization step - write results and other computational vectors to files
 *
 * @date 22-Oct-2012
 * @author V. Petkov
 */

#include <stdio.h>
#include "util_write_files.h"
#include "mpi.h"
#include "finalization.h"
#include <stdlib.h>

void finalization(char* file_in, char* out_prefix, int total_iters, double residual_ratio,
        int nintci, int nintcf, int points_count, int** points, int* elems, double* var,
        double* cgup, double* su, int num_elem_global, int *local_global_index) {

    int my_rank, num_procs;
    int i;
    double *cgup_global, *cgup_global_all, *var_global, *var_global_all;
    cgup_global = (double*) calloc(sizeof(double), num_elem_global);
    cgup_global_all = (double*) calloc(sizeof(double), num_elem_global);
    var_global = (double*) calloc(sizeof(double), num_elem_global);
    var_global_all = (double*) calloc(sizeof(double), num_elem_global);

    for (i = 0; i <= nintcf; i++) {
        cgup_global[local_global_index[i]] = cgup[i];
        var_global[local_global_index[i]] = var[i];
    }

    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);  /// Get current process id
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);  /// get number of processes

    MPI_Reduce(cgup_global, cgup_global_all, num_elem_global, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(var_global, var_global_all, num_elem_global, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    if (my_rank == 0) {

        char file_out[100];
        sprintf(file_out, "%s_summary.out", out_prefix);

        int status = store_simulation_stats(file_in, file_out, nintci, nintcf, var, total_iters,
                residual_ratio);

        sprintf(file_out, "%s_data.vtk", out_prefix);
        // vtk_write_unstr_grid_header(file_in, file_out, nintci, nintcf, points_count, points, elems);
        vtk_write_unstr_grid_header(file_in, file_out, 0, num_elem_global - 1, points_count, points,
                elems);
        vtk_append_double(file_out, "VAR", 0, num_elem_global, var_global_all);

        if (status != 0)
            fprintf(stderr, "Error when trying to write to file %s\n", file_out);
    }
    free(cgup_global);
    free(cgup_global_all);
}

