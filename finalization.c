/**
 * Finalization step - write results and other computational vectors to files
 *
 * @date 22-Oct-2012
 * @author V. Petkov
 */

#include <stdio.h>
#include "util_write_files.h"
#include "mpi.h"

void finalization(char* file_in, char* out_prefix, int total_iters, double residual_ratio,
        int nintci, int nintcf, int points_count, int** points, int* elems, double* var,
        double* cgup, double* su) {
    /*
     int my_rank, num_procs, i;
     MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);  /// Get current process id
     MPI_Comm_size(MPI_COMM_WORLD, &num_procs);  /// get number of processes
     double *var_g, *cgup_g, *su_g;
     int *num_elem_each;
     num_elem_each = (int*) malloc(num_procs*sizeof(int));

     // prepare for finalization
     printf("nintcf+1 on %d=%d\n", my_rank, nintcf + 1);

     // MPI_Gather(&nintcf, 1, MPI_INT, num_elem_each, 1, MPI_INT, 0, MPI_COMM_WORLD);
     MPI_Barrier(MPI_COMM_WORLD);
     if (my_rank==0){
     for (i=0;i<num_procs;i++){
     // printf("num_elem_each[%d]=%d\n",i,num_elem_each[i]);
     }
     }

     */
    // MPI_Gather(var, nintcf + 1, MPI_DOUBLE, var_g, num_global_elem, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    // MPI_Gather(cgup, nintcf + 1, MPI_DOUBLE, cgup_g, num_global_elem, MPI_DOUBLE, 0,
    //        MPI_COMM_WORLD);
    // MPI_Gather(su, nintcf + 1, MPI_DOUBLE, su_g, num_global_elem, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    /*
     char file_out[100];
     sprintf(file_out, "%s_summary.out", out_prefix);

     int status = store_simulation_stats(file_in, file_out, nintci, nintcf, var, total_iters,
     residual_ratio);

     sprintf(file_out, "%s_data.vtk", out_prefix);
     vtk_write_unstr_grid_header(file_in, file_out, nintci, nintcf, points_count, points, elems);
     vtk_append_double(file_out, "CGUP", nintci, nintcf, cgup);

     if (status != 0)
     fprintf(stderr, "Error when trying to write to file %s\n", file_out);
     */
}

