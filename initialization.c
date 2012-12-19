/**
 * Initialization step - parse the input file, compute data distribution, initialize LOCAL computational arrays
 *
 * @date 22-Oct-2012
 * @author V. Petkov
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <mpi.h>
#include <metis.h>
#include <math.h>
#include "util_read_files.h"
#include "initialization.h"

int initialization(char* file_in, char* part_type, int* nintci, int* nintcf, int* nextci,
        int* nextcf, int*** lcc, double** bs, double** be, double** bn, double** bw, double** bl,
        double** bh, double** bp, double** su, int* points_count, int*** points, int** elems,
        double** var, double** cgup, double** oc, double** cnorm, int** local_global_index,
        int** global_local_index, int* neighbors_count, int** send_count, int*** send_list,
        int** recv_count, int*** recv_list, int** epart, int** npart, int* objval,
        int* num_global_elem) {
    /********** START INITIALIZATION **********/

    int i = 0;
    int j, k;
    // read-in the input file
    int f_status = read_binary_geo(file_in, &*nintci, &*nintcf, &*nextci, &*nextcf, &*lcc, &*bs,
            &*be, &*bn, &*bw, &*bl, &*bh, &*bp, &*su, &*points_count, &*points, &*elems);

    if (f_status != 0)
        return f_status;
    int my_rank, num_procs;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);  // Get current process id
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);  // get number of processes
    int num_elem_g = *nintcf - *nintci + 1;
    *num_global_elem = num_elem_g;

    // classical partitioning
    if (strcmp(part_type, "classical") == 0) {
        int num_elem = floor(num_elem_g / num_procs);
        if (my_rank == num_procs - 1) {
            num_elem = num_elem_g - num_elem * (num_procs - 1);
        }

        // make the index mapping arrays
        // allocate memory for the mapping arrays
        if ((*global_local_index = (int *) malloc(num_elem_g * sizeof(int))) == NULL ) {
            fprintf(stderr, "malloc failed for global_local_index\n");
            return -1;
        }
        if ((*local_global_index = (int *) malloc(num_elem * sizeof(int))) == NULL ) {
            fprintf(stderr, "malloc failed for local_global_index\n");
            return -1;
        }
        int id;

        if ((*epart = (int*) calloc(sizeof(int), num_elem_g)) == NULL ) {
            fprintf(stderr, "calloc failed for epart\n");
            return -1;
        }
        for (i = 0; i < num_elem; i++) {
            id = my_rank * floor(num_elem_g / num_procs) + i;
            (*local_global_index)[i] = id;
        }

        // fill epart
        for (i = 0; i < num_procs - 1; i++) {  // for all -1 processors
            for (j = 0; j < num_elem; j++) {
                id = i * floor(num_elem_g / num_procs) + j;
                (*epart)[id] = i;
            }
        }
        // for the last processor
        for (j = 0; j < num_elem_g - (num_procs - 1) * floor(num_elem_g / num_procs); j++) {
            id = (num_procs - 1) * floor(num_elem_g / num_procs) + j;
            (*epart)[id] = num_procs - 1;
        }

        // define local arrays
        double *BS_l, *BE_l, *BN_l, *BW_l, *BL_l, *BH_l, *BP_l, *SU_l;
        // allocate memory for the local arrays
        allocate_local_arrays(&BS_l, &BE_l, &BN_l, &BW_l, &BL_l, &BH_l, &BP_l, &SU_l, num_elem);

        // copy the values from the global to the local arrays
        j = 0;
        for (i = 0; i < num_elem_g; i++) {
            if (my_rank == (*epart)[i]) {
                (*local_global_index)[j] = i;
                BS_l[j] = (*bs)[i];
                BE_l[j] = (*be)[i];
                BN_l[j] = (*bn)[i];
                BW_l[j] = (*bw)[i];
                BL_l[j] = (*bl)[i];
                BH_l[j] = (*bh)[i];
                BP_l[j] = (*bp)[i];
                SU_l[j] = (*su)[i];
                j++;
            }
        }
        /*
         int index;
         for (i = 0; i < num_elem; i++) {
         index = my_rank * floor(num_elem_g / num_procs) + i;
         BS_l[i] = (*bs)[index];
         BE_l[i] = (*be)[index];
         BN_l[i] = (*bn)[index];
         BW_l[i] = (*bw)[index];
         BL_l[i] = (*bl)[index];
         BH_l[i] = (*bh)[index];
         BP_l[i] = (*bp)[index];
         SU_l[i] = (*su)[index];
         }
         */

        // free the resources of the global arrays
        free(*bp);
        free(*bh);
        free(*bl);
        free(*bw);
        free(*bn);
        free(*be);
        free(*bs);
        free(*su);

        // set the original pointers to the local pointers
        *bs = BS_l;
        *bh = BH_l;
        *bl = BL_l;
        *bw = BW_l;
        *bn = BN_l;
        *be = BE_l;
        *bp = BP_l;
        *su = SU_l;

        // allocate cgup and initialize it with zeros
        *cgup = (double*) calloc(sizeof(double), num_elem);

        for (i = 0; i < num_elem; i++) {
            (*cgup)[i] = 1.0 / BP_l[i];
        }
        *oc = (double*) calloc(sizeof(double), num_elem);
        *cnorm = (double*) calloc(sizeof(double), num_elem);

        *nintcf = num_elem - 1;

    } else if (strcmp(part_type, "dual") == 0) {
        // distribution with METIS
        // create all arrays in the METIS format
        idx_t ne = *nintcf - *nintci + 1;
        idx_t nn = *points_count;
        idx_t *eptr;
        idx_t *eind;
        idx_t ncommon = 4;
        idx_t nparts = num_procs;
        idx_t objval_m = *objval;
        idx_t *epart_m;
        idx_t *npart_m;

        if ((eptr = (idx_t*) calloc(sizeof(idx_t), ne + 1)) == NULL ) {
            fprintf(stderr, "calloc failed for eptr\n");
            return -1;
        }

        if ((eind = (idx_t*) calloc(sizeof(idx_t), ne * 8)) == NULL ) {
            fprintf(stderr, "calloc failed for eptr\n");
            return -1;
        }
        if ((epart_m = (idx_t*) calloc(sizeof(idx_t), ne)) == NULL ) {
            fprintf(stderr, "calloc failed for epart\n");
            return -1;
        }
        if ((npart_m = (idx_t*) calloc(sizeof(idx_t), nn)) == NULL ) {
            fprintf(stderr, "calloc failed for eptr\n");
            return -1;
        }

        int j = 0;
        for (i = 0; i <= ne; i++) {
            eptr[i] = j;
            j = j + 8;
        }

        for (i = 0; i <= ne * 8; i++) {
            eind[i] = (*elems)[i];
        }

        METIS_PartMeshDual(&ne, &nn, eptr, eind, NULL, NULL, &ncommon, &nparts, NULL, NULL,
                &objval_m, epart_m, npart_m);

        // convert the values back to int
        // allocate memory
        if ((*epart = (int*) calloc(sizeof(int), *nintcf - *nintci + 1)) == NULL ) {
            fprintf(stderr, "calloc failed for epart\n");
            return -1;
        }
        if ((*npart = (int*) calloc(sizeof(int), *points_count)) == NULL ) {
            fprintf(stderr, "calloc failed for eptr\n");
            return -1;
        }
        // copy the values back
        *objval = objval_m;
        for (i = 0; i < ne; i++) {
            (*epart)[i] = epart_m[i];
        }
        for (i = 0; i < nn; i++) {
            (*npart)[i] = npart_m[i];
        }

        // distribute the arrays
        // Calculate the numbers of elements for each process
        int ne_l = 0;
        for (i = 0; i < ne; i++) {
            if (my_rank == (*epart)[i]) {
                ne_l++;
            }
        }

        // make the index mapping arrays
        // allocate memory for the mapping arrays
        if ((*global_local_index = (int *) malloc((*nintcf - *nintci + 1) * sizeof(int))) == NULL ) {
            fprintf(stderr, "malloc failed for global_local_index\n");
            return -1;
        }
        if ((*local_global_index = (int *) malloc(ne_l * sizeof(int))) == NULL ) {
            fprintf(stderr, "malloc failed for local_global_index\n");
            return -1;
        }

        // define local arrays
        double *BS_l, *BE_l, *BN_l, *BW_l, *BL_l, *BH_l, *BP_l, *SU_l;
        // allocate memory for the local arrays
        allocate_local_arrays(&BS_l, &BE_l, &BN_l, &BW_l, &BL_l, &BH_l, &BP_l, &SU_l, ne_l);
        // copy the values from the global to the local arrays
        j = 0;
        for (i = 0; i < ne; i++) {
            if (my_rank == (*epart)[i]) {
                (*local_global_index)[j] = i;
                BS_l[j] = (*bs)[i];
                BE_l[j] = (*be)[i];
                BN_l[j] = (*bn)[i];
                BW_l[j] = (*bw)[i];
                BL_l[j] = (*bl)[i];
                BH_l[j] = (*bh)[i];
                BP_l[j] = (*bp)[i];
                SU_l[j] = (*su)[i];
                j++;
            }
        }

        // free the resources of the global arrays
        free(*bp);
        free(*bh);
        free(*bl);
        free(*bw);
        free(*bn);
        free(*be);
        free(*bs);
        free(*su);

        // set the original pointers to the local pointers
        *bs = BS_l;
        *bh = BH_l;
        *bl = BL_l;
        *bw = BW_l;
        *bn = BN_l;
        *be = BE_l;
        *bp = BP_l;
        *su = SU_l;

        // allocate cgup and initialize it with zeros
        *cgup = (double*) calloc(sizeof(double), ne_l);

        for (i = 0; i < ne_l; i++) {
            (*cgup)[i] = 1.0 / BP_l[i];
        }
        *oc = (double*) calloc(sizeof(double), ne_l);
        *cnorm = (double*) calloc(sizeof(double), ne_l);

        *nintcf = ne_l - 1;
    }
    *nintci = 0;

    // make the global_local_index
    int *gl_lo;
    if ((gl_lo = (int*) calloc(sizeof(int), num_procs)) == NULL ) {
        fprintf(stderr, "calloc failed for gl_lo\n");
        return -1;
    }

    for (i = 0; i < num_elem_g; i++) {
        (*global_local_index)[i] = gl_lo[(*epart)[i]];
        gl_lo[(*epart)[i]]++;
    }

    int **lcc_l;
    // allocating the local LCC array
    if ((lcc_l = (int**) malloc(((*nintcf) + 1) * sizeof(int*))) == NULL ) {
        fprintf(stderr, "malloc failed to allocate first dimension of lcc_l");
        return -1;
    }

    for (i = 0; i <= *nintcf; i++) {
        if ((lcc_l[i] = (int*) malloc(6 * sizeof(int))) == NULL ) {
            fprintf(stderr, "malloc failed to allocate second dimension of lcc_l\n");
            return -1;
        }
    }

    // distribution LCC
    int *dis_lcc;
    if ((dis_lcc = (int*) calloc(sizeof(int), num_procs)) == NULL ) {
        fprintf(stderr, "calloc failed for gl_lo\n");
        return -1;
    }

    // distribute LCC
    for (i = 0; i < num_elem_g; i++) {
        if (my_rank == (*epart)[i]) {
            k = dis_lcc[my_rank];
            for (j = 0; j < 6; j++) {
                lcc_l[k][j] = (*lcc)[i][j];
            }
            dis_lcc[my_rank]++;
        }
    }
    free(dis_lcc);

    // free the memory of the global LCC vector
    for (int i = 0; i < num_elem_g; i++) {
        free((*lcc)[i]);
    }
    free(*lcc);
    // set the original pointer to the local lcc pointer
    *lcc = lcc_l;

    // End partitioning

    // make the communication model

    comm_model(*nintcf + 1, num_elem_g, lcc, local_global_index, global_local_index,
            neighbors_count, send_count, send_list, recv_count, recv_list, epart, nextci, nextcf);

    *var = (double*) calloc(sizeof(double), (*nintcf + 1));
    // cgup oc and cnorm is a local array only
    // *cgup = (double*) calloc(sizeof(double), (*nextcf + 1));

    // *oc = (double*) calloc(sizeof(double), (*nintcf + 1));
    // *cnorm = (double*) calloc(sizeof(double), (*nintcf + 1));

    // initialize the arrays

    for (i = 0; i <= 10; i++) {
        // (*oc)[i] = 0.0;
        (*cnorm)[i] = 1.0;  // only the firt 3 elements have to be one on each proseccor
    }
    /*
     for (i = 0; i < (*nintcf); i++) {

     (*cgup)[i] = 0.0;
     (*var)[i] = 0.0;
     }

     for (i = (*nextci); i <= (*nextcf); i++) {
     (*var)[i] = 0.0;

     (*cgup)[i] = 0.0;
     (*bs)[i] = 0.0;
     (*be)[i] = 0.0;
     (*bn)[i] = 0.0;
     (*bw)[i] = 0.0;
     (*bl)[i] = 0.0;
     (*bh)[i] = 0.0;
     }
     */
    for (i = 0; i <= (*nintcf); i++)
        (*cgup)[i] = 1.0 / ((*bp)[i]);

    return 0;
}
int allocate_local_arrays(double **BS_l, double **BE_l, double **BN_l, double **BW_l, double **BL_l,
        double **BH_l, double **BP_l, double **SU_l, int num_elem) {
// allocate other arrays
    if ((*BS_l = (double *) malloc(num_elem * sizeof(double))) == NULL ) {
        fprintf(stderr, "malloc(BS) failed\n");
        return -1;
    }

    if ((*BE_l = (double *) malloc(num_elem * sizeof(double))) == NULL ) {
        fprintf(stderr, "malloc(BE) failed\n");
        return -1;
    }

    if ((*BN_l = (double *) malloc(num_elem * sizeof(double))) == NULL ) {
        fprintf(stderr, "malloc(BN) failed\n");
        return -1;
    }

    if ((*BW_l = (double *) malloc(num_elem * sizeof(double))) == NULL ) {
        fprintf(stderr, "malloc(BW) failed\n");
        return -1;
    }

    if ((*BL_l = (double *) malloc(num_elem * sizeof(double))) == NULL ) {
        fprintf(stderr, "malloc(BL) failed\n");
        return -1;
    }

    if ((*BH_l = (double *) malloc(num_elem * sizeof(double))) == NULL ) {
        fprintf(stderr, "malloc(BH) failed\n");
        return -1;
    }

    if ((*BP_l = (double *) malloc(num_elem * sizeof(double))) == NULL ) {
        fprintf(stderr, "malloc(BP) failed\n");
        return -1;
    }
    if ((*SU_l = (double *) malloc(num_elem * sizeof(double))) == NULL ) {
        fprintf(stderr, "malloc(SU) failed\n");
        return -1;
    }
    return 0;
}

int comm_model(int ne_l, int ne_g, int*** lcc, int** local_global_index, int** global_local_index,
        int* neighbors_count, int** send_count, int*** send_list, int** recv_count,
        int*** recv_list, int** epart, int* nextci, int* nextcf) {
    int i, j, id, k, total_send_recv, ext_pos;
    int my_rank, num_procs;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);  // Get current process id
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);  // get number of processes

    // allocate memory for send_count and recv_count
    if ((*send_count = (int*) calloc(sizeof(int), num_procs)) == NULL ) {
        fprintf(stderr, "calloc failed for send_count_w\n");
        return -1;
    }
    if ((*recv_count = (int*) calloc(sizeof(int), num_procs)) == NULL ) {
        fprintf(stderr, "calloc failed for recv_count_w\n");
        return -1;
    }
    // allocate memory for send_count_cum (cumulative)
    int *send_count_cum;
    if ((send_count_cum = (int*) calloc(sizeof(int), num_procs)) == NULL ) {
        fprintf(stderr, "calloc failed for send_count_cum\n");
        return -1;
    }

    // calculate the number of neighbors
    int p_n;
    for (i = 0; i < ne_l; i++) {  // for all elements in the process
        for (j = 0; j < 6; j++) {  // for all neighbors of this element
            id = (*lcc)[i][j];  // get the global id for this neighbor
            if (id < ne_g) {  // get whether the neighbor is an external cell
                p_n = (*epart)[id];  // get the process of this neighbor
                if (p_n != my_rank) {
                    (*send_count)[p_n]++;
                }
            }
        }
    }

    // allocating send_list
    if ((*send_list = (int**) malloc(num_procs * sizeof(int*))) == NULL ) {
        fprintf(stderr, "malloc failed to allocate first dimension of LCC");
        return -1;
    }

    for (i = 0; i < num_procs; i++) {
        if (((*send_list)[i] = (int *) malloc((*send_count)[i] * sizeof(int))) == NULL ) {
            fprintf(stderr, "malloc failed to allocate second dimension of LCC\n");
            return -1;
        }
    }

    // allocating recv_list

    if ((*recv_list = (int**) malloc(num_procs * sizeof(int*))) == NULL ) {
        fprintf(stderr, "malloc failed to allocate first dimension of LCC");
        return -1;
    }

    for (i = 0; i < num_procs; i++) {
        if (((*recv_list)[i] = (int *) malloc((*send_count)[i] * sizeof(int))) == NULL ) {
            fprintf(stderr, "malloc failed to allocate second dimension of LCC\n");
            return -1;
        }
    }

    // create a new lcc array
    // this array maps the neighbors in the following order:
    // internal cells[ne_l]; ghost cells[total_send_recv]; external cells[1]
    int **lcc_n;
    // allocating the new LCC array
    if ((lcc_n = (int**) malloc(ne_l * sizeof(int*))) == NULL ) {
        fprintf(stderr, "malloc failed to allocate first dimension of lcc_l");
        return -1;
    }
    for (i = 0; i < ne_l; i++) {
        if ((lcc_n[i] = (int*) malloc(6 * sizeof(int))) == NULL ) {
            fprintf(stderr, "malloc failed to allocate second dimension of lcc_l\n");
            return -1;
        }
    }

    // allocate memory for the send_list_position
    int *send_list_pos;
    if ((send_list_pos = (int*) calloc(sizeof(int), num_procs)) == NULL ) {
        fprintf(stderr, "calloc failed for send_count_w\n");
        return -1;
    }
    // copy send count to receive count because it is the same
    total_send_recv = 0;
    for (i = 0; i < num_procs; i++) {
        (*recv_count)[i] = (*send_count)[i];
        total_send_recv = total_send_recv + (*send_count)[i];
        if (i != 0) {
            send_count_cum[i] = send_count_cum[i - 1] + (*send_count)[i - 1];
        }
    }

    *nextci = ne_l;
    *nextcf = ne_l + total_send_recv - 1;

    // fill the send list
    int p;
    k = 0;
    ext_pos = *nextcf+1;
    for (i = 0; i < ne_l; i++) {  // for all elements in the process
        for (j = 0; j < 6; j++) {  // for all neighbors of this element
            id = (*lcc)[i][j];  // get the global id for this neighbor
            if (id < ne_g) {  // get whether the neighbor is an external cell
                p_n = (*epart)[id];  // get the process of this neighbor
                if (p_n != my_rank) {  // check whether it is an element of another process
                    p = send_list_pos[p_n];
                    (*send_list)[p_n][p] = (*local_global_index)[i];
                    (*recv_list)[p_n][p] = id;
                    send_list_pos[p_n]++;  // increase send_list_count
                    lcc_n[i][j] = *nextci + send_count_cum[p_n] + p;
                } else {
                    lcc_n[i][j] = (*global_local_index)[id];
                }
            } else {
                lcc_n[i][j] = ext_pos;
            }
        }
    }

    // free the memory of the global LCC vector
    for (int i = 0; i < ne_l; i++) {
        free((*lcc)[i]);
    }
    free(*lcc);
    // set the original pointer to the new lcc pointer
    *lcc = lcc_n;

    /*
     // allocating recv_list
     if ((*recv_list = (int**) malloc(num_procs * sizeof(int*))) == NULL ) {
     fprintf(stderr, "malloc failed to allocate first dimension of LCC");
     return -1;
     }

     for (i = 0; i < num_procs; i++) {
     if (((*recv_list)[i] = (int *) malloc((*recv_count)[i] * sizeof(int))) == NULL ) {
     fprintf(stderr, "malloc failed to allocate second dimension of LCC\n");
     return -1;
     }
     }
     */
    /*
     MPI_Request req_s[num_procs];
     MPI_Request req_r[num_procs];
     MPI_Status status_s[num_procs];
     MPI_Status status_r[num_procs];

     for (i = 0; i < num_procs; i++) {

     MPI_Isend(&(*send_list)[i][0], (*send_count)[i], MPI_INT, i, 1, MPI_COMM_WORLD, &req_s[i]);

     }

     for (i = 0; i < num_procs; i++) {

     MPI_Irecv(&(*recv_list)[i][0], (*recv_count)[i], MPI_INT, i, 1, MPI_COMM_WORLD, &req_r[i]);

     }
     MPI_Waitall(num_procs, req_s, status_s);
     MPI_Waitall(num_procs, req_r, status_r);
     */

    return 1;
}
