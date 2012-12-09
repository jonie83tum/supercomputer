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
        int** recv_count, int*** recv_list, int** epart, int** npart, int* objval) {
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
        for (j = 0; j < num_elem_g - floor(num_elem_g / num_procs); j++) {
            id = (num_procs - 1) * floor(num_elem_g / num_procs) + j;
            (*epart)[id] = num_procs - 1;
        }

        // define local arrays
        double *BS_l, *BE_l, *BN_l, *BW_l, *BL_l, *BH_l, *BP_l, *SU_l;
        // allocate memory for the local arrays
        allocate_local_arrays(&BS_l, &BE_l, &BN_l, &BW_l, &BL_l, &BH_l, &BP_l, &SU_l, num_elem);

        // copy the values from the global to the local arrays
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
        // free the resources of the global arrays
        free(*bp);
        free(*bh);
        free(*bl);
        free(*bw);
        free(*bn);
        free(*be);
        free(*bs);
        free(*su);

        // allocate cgup and initialize it with zeros
        *cgup = (double*) calloc(sizeof(double), num_elem);

        for (i = 0; i < num_elem; i++) {
            (*cgup)[i] = 1.0 / BP_l[i];
        }
        *oc = (double*) calloc(sizeof(double), num_elem);
        *cnorm = (double*) calloc(sizeof(double), num_elem);

        *nintcf = num_elem;

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
        if ((*epart = (int*) calloc(sizeof(int), ne)) == NULL ) {
            fprintf(stderr, "calloc failed for epart\n");
            return -1;
        }
        if ((*npart = (int*) calloc(sizeof(int), nn)) == NULL ) {
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
        if ((*global_local_index = (int *) malloc(ne * sizeof(int))) == NULL ) {
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

        // allocate cgup and initialize it with zeros
        *cgup = (double*) calloc(sizeof(double), ne_l);

        for (i = 0; i < ne_l; i++) {
            (*cgup)[i] = 1.0 / BP_l[i];
        }
        *oc = (double*) calloc(sizeof(double), ne_l);
        *cnorm = (double*) calloc(sizeof(double), ne_l);

        *nintcf = ne_l;
    }
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
    if ((lcc_l = (int*) malloc(*nintcf * sizeof(int*))) == NULL ) {
        fprintf(stderr, "malloc failed to allocate first dimension of lcc_l");
        return -1;
    }

    for (i = 0; i < *nintcf; i++) {
        if ((lcc_l[i] = (int) malloc(6 * sizeof(int))) == NULL ) {
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
    comm_model(*nintcf, num_elem_g, lcc, local_global_index, global_local_index, neighbors_count,
            send_count, send_list, recv_count, recv_list, epart);

    *var = (double*) calloc(sizeof(double), (*nextcf + 1));
    // cgup oc and cnorm is a local array only
    // *cgup = (double*) calloc(sizeof(double), (*nextcf + 1));

    // *oc = (double*) calloc(sizeof(double), (*nintcf + 1));
    // *cnorm = (double*) calloc(sizeof(double), (*nintcf + 1));

    // initialize the arrays
    for (i = 0; i <= 10; i++) {
        (*oc)[i] = 0.0;
        (*cnorm)[i] = 1.0;
    }

    for (i = (*nintci); i <= (*nintcf); i++) {
        // cgup is only a local array an it is anyway initialized with calloc
        // (*cgup)[i] = 0.0;
        (*var)[i] = 0.0;
    }

    for (i = (*nextci); i <= (*nextcf); i++) {
        (*var)[i] = 0.0;
        // cgup is only a local array an it is anyway initialized with calloc
        // (*cgup)[i] = 0.0;
        /*do not use the global arrays any more because they are already destroyed
         (*bs)[i] = 0.0;
         (*be)[i] = 0.0;
         (*bn)[i] = 0.0;
         (*bw)[i] = 0.0;
         (*bl)[i] = 0.0;
         (*bh)[i] = 0.0;
         */
    }
    /*do not use bp anymore because it is already destroyed
     for (i = (*nintci); i <= (*nintcf); i++)
     (*cgup)[i] = 1.0 / ((*bp)[i]);
     */
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
        int*** recv_list, int** epart) {

    int i, j, k, l, id;
    int my_rank, num_procs;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);  // Get current process id
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);  // get number of processes

    // helping arrays
    int *send_count_w, *recv_count_w;
    // allocate memory for the helping arrays
    if ((send_count_w = (int*) calloc(sizeof(int), num_procs)) == NULL ) {
        fprintf(stderr, "calloc failed for send_count_w\n");
        return -1;
    }
    if ((recv_count_w = (int*) calloc(sizeof(int), num_procs)) == NULL ) {
        fprintf(stderr, "calloc failed for recv_count_w\n");
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
                    send_count_w[p_n]++;
                }
            }
        }
    }

    // calculate number of neighbors and
    *neighbors_count = 0;
    for (i = 0; i < num_procs; i++) {
        if (send_count_w[i] != 0) {
            *neighbors_count = *neighbors_count + 1;
        }
    }
    // allocate memory and fill send_count
    int process_map_s[*neighbors_count];
    if ((*send_count = (int*) calloc(sizeof(int), *neighbors_count)) == NULL ) {
        fprintf(stderr, "calloc failed for send_count_w\n");
        return -1;
    }
    j = 0;
    for (i = 0; i < num_procs; i++) {
        if (send_count_w[i] != 0) {
            (*send_count)[j] = send_count_w[i];
            process_map_s[j] = i;
            j++;
        }
    }

    // allocate memory for send_list
    if ((*send_list = (int**) calloc(sizeof(int*), num_procs)) == NULL ) {
        fprintf(stderr, "calloc failed for first dimension of send_list\n");
        return -1;
    }
    for (i = 0; i < num_procs; i++) {
        if (((*send_list)[i] = (int*) calloc(sizeof(int), send_count_w[i])) == NULL ) {
            fprintf(stderr, "calloc failed for second dimension of send_list\n");
            return -1;
        }

    }

    // fill the send list
    for (k = 0; k < num_procs; k++) {  // for all neighbors
        if (send_count_w[k] != 0) {
            l = 0;
            for (i = 0; i < ne_l; i++) {  // for all elements in the process
                for (j = 0; j < 6; j++) {  // for all neighbors of this element
                    id = (*lcc)[i][j];  // get the global id for this neighbor
                    if (id < ne_g) {  // get whether the neighbor is an external cell
                        p_n = (*epart)[id];  // get the process of this neighbor
                        if (p_n == process_map_s[k]) {
                            (*send_list)[k][l] = (*local_global_index)[i];
                            l++;
                            if (my_rank == 0 && k == 0) {
                            }
                        }
                    }
                }
            }
        }
    }

    // allocate memory for the recv_count (same size as send_count)
    if ((*recv_count = (int*) calloc(sizeof(int), *neighbors_count)) == NULL ) {
        fprintf(stderr, "calloc failed for send_count_w\n");
        return -1;
    }
/*
    if (my_rank == 0) {
        for (j = 0; j < num_procs; j++) {
            printf("send_count_w[%d]=%d\n", j, send_count_w[j]);
        }
        for (j = 0; j < num_procs; j++) {
            printf("recv_count_w[%d]=%d\n", j, recv_count_w[j]);
        }
    }
    /*
     // scatter all send_count_w (becomes the according recv_count_w)
     for (i = 0; i < num_procs; i++) {
     MPI_Scatter(send_count_w, 1, MPI_INT, &recv_count_w[i], 1, MPI_INT, i, MPI_COMM_WORLD);
     }

     // allocate memory and fill recv_count
     int process_map_r[*neighbors_count];
     if ((*recv_count = (int*) calloc(sizeof(int), *neighbors_count)) == NULL ) {
     fprintf(stderr, "calloc failed for send_count_w\n");
     return -1;
     }
     j = 0;
     for (i = 0; i < num_procs; i++) {
     if (recv_count_w[i] != 0) {
     (*recv_count)[j] = recv_count_w[i];
     process_map_r[j] = i;
     j++;
     }
     }

     // allocate memory for recv_list
     if ((*recv_list = (int**) calloc(sizeof(int*), num_procs)) == NULL ) {
     fprintf(stderr, "calloc failed for first dimension of recv_list\n");
     return -1;
     }
     for (i = 0; i < *neighbors_count; i++) {
     if (((*recv_list)[i] = (int*) calloc(sizeof(int), recv_count_w[i])) == NULL ) {
     fprintf(stderr, "calloc failed for second dimension of recv_list\n");
     return -1;
     }

     }

     // fill recv_list
     MPI_Request req;
     for (i = 0; i < *neighbors_count; i++) {
     MPI_Isend((*send_list)[i], (*send_count)[i], MPI_INT, i, 1, MPI_COMM_WORLD,
     &req);
     }
     for (i = 0; i < *neighbors_count; i++) {
     MPI_Irecv((*recv_list)[i], (*recv_count)[i], MPI_INT, i, 1, MPI_COMM_WORLD,
     &req);
     }
     MPI_Wait(&req, MPI_STATUS_IGNORE);
     /*
     if (my_rank==0){
     for (i=0;i<(*send_count)[0];i++){
     printf("recv_list[0][%d]=%d\n", i, (*recv_list)[0][i]);
     }
     }
     */
    return 1;
}
