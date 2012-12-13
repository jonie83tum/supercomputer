/**
 * Functions to test the data distribution and communication lists creation algorithms
 *
 * @date 22-Oct-2012
 * @author V. Petkov
 */
#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include "util_read_files.h"
#include "util_write_files.c"
int test_distribution(char *file_in, char *file_vtk_out, int *local_global_index, int num_elems,
        double *cgup) {
    int i;
    // Read the geometry from the input file

    /** Simulation parameters parsed from the input datasets */
    int nintci, nintcf;  /// internal cells start and end index
    /// external cells start and end index. The external cells are only ghost cells.
    /// They are accessed only through internal cells
    int nextci, nextcf;
    int **lcc;  /// link cell-to-cell array - stores neighboring information
    /// Boundary coefficients for each volume cell (South, East, North, West, High, Low)
    double *bs, *be, *bn, *bw, *bl, *bh;
    double *bp;  /// Pole coefficient
    double *su;  /// Source values
    /** Geometry data */
    int points_count;  /// total number of points that define the geometry
    int** points;  /// coordinates of the points that define the cells - size [points_cnt][3]
    int* elems;  /// definition of the cells using their nodes (points) - each cell has 8 points
    double* distr;

    int read_input = read_binary_geo(file_in, &nintci, &nintcf, &nextci, &nextcf, &lcc, &bs, &be,
            &bn, &bw, &bl, &bh, &bp, &su, &points_count, &points, &elems);

    if (read_input != 0)
        printf("Could not read input file in test distribution");

    // allocate the distr vector and initialize it with zeros
    if ((distr = (double*) calloc(sizeof(double), nintcf - nintci + 1)) == NULL ) {
        fprintf(stderr, "calloc failed to allocate distr");
        return -1;
    }

    // copy the locally initialized CGUP vector to the right place in the distr array
    int ind;
    for (i = 0; i < num_elems; i++) {
        ind = local_global_index[i];
        distr[ind] = cgup[i];
    }

    // write the vtk header
    const char experiment_name[] = "test for distribution";
    vtk_write_unstr_grid_header(experiment_name, file_vtk_out, nintci, nintcf, points_count, points,
            elems);
    // write the values to the vtk file
    vtk_append_double(file_vtk_out, "CGUP", nintci, nintcf, distr);
    return 0;
}

int test_communication(char *file_in, char *file_vtk_out, int *local_global_index, int num_elems,
        int neighbors_count, int* send_count, int** send_list, int* recv_count, int** recv_list) {
    int i, j, id;
    int my_rank, num_procs;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);  // Get current process id
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);  // get number of processes

    // Read the geometry from the input file

    /** Simulation parameters parsed from the input datasets */
    int nintci, nintcf;  /// internal cells start and end index
    /// external cells start and end index. The external cells are only ghost cells.
    /// They are accessed only through internal cells
    int nextci, nextcf;
    int **lcc;  /// link cell-to-cell array - stores neighboring information
    /// Boundary coefficients for each volume cell (South, East, North, West, High, Low)
    double *bs, *be, *bn, *bw, *bl, *bh;
    double *bp;  /// Pole coefficient
    double *su;  /// Source values
    /** Geometry data */
    int points_count;  /// total number of points that define the geometry
    int** points;  /// coordinates of the points that define the cells - size [points_cnt][3]
    int* elems;  /// definition of the cells using their nodes (points) - each cell has 8 points
    double* commlist;
    int ne_g;

    int read_input = read_binary_geo(file_in, &nintci, &nintcf, &nextci, &nextcf, &lcc, &bs, &be,
            &bn, &bw, &bl, &bh, &bp, &su, &points_count, &points, &elems);

    if (read_input != 0)
        printf("Could not read input file in test distribution");

    ne_g = nintcf - nintci + 1;
    // allocate the commlist vector and initialize it with zeros
    if ((commlist = (double*) calloc(sizeof(double), ne_g)) == NULL ) {
        fprintf(stderr, "calloc failed to allocate distr");
        return -1;
    }

    // set all internal cells
    for (i = 0; i < num_elems; i++) {
        id = local_global_index[i];  // get global id of the element
        commlist[id] = 15.0;
    }

    // set the elements which are to be sent
    for (i = 0; i < num_procs; i++) {  // for all neighbors in of this process
        for (j = 0; j < send_count[i]; j++) {  // for all elements in the list
            id = send_list[i][j];  // get global id of the element to be sent
            commlist[id] = 10.0;
        }
    }

    // set the elements which are to be received
    for (i = 0; i < num_procs; i++) {  // for all neighbors in of this process
        for (j = 0; j < recv_count[i]; j++) {  // for all elements in the list
            id = recv_list[i][j];  // get global id of the element to be sent
            commlist[id] = 5.0;
        }
    }

    for (i = 0; i < num_procs; i++) {
        printf("recv_count[%d]=%d\n", i, recv_count[i]);
    }

    // write the vtk header
    const char experiment_name[] = "test for communication";
    vtk_write_unstr_grid_header(experiment_name, file_vtk_out, nintci, nintcf, points_count, points,
            elems);
    // write the values to the vtk file
    vtk_append_double(file_vtk_out, "commlist", nintci, nintcf, commlist);
    /*
     *
     */
    return 0;
}

