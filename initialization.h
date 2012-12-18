/**
 * Initialization step - parse the input file, compute data distribution, initialize LOCAL computational arrays
 *
 * @date 22-Oct-2012
 * @author V. Petkov
 */

#ifndef INITIALIZATION_H_
#define INITIALIZATION_H_

int initialization(char* file_in, char* part_type, int* nintci, int* nintcf, int* nextci,
        int* nextcf, int*** lcc, double** bs, double** be, double** bn, double** bw, double** bl,
        double** bh, double** bp, double** su, int* points_count, int*** points, int** elems,
        double** var, double** cgup, double** oc, double** cnorm, int** local_global_index,
        int** global_local_index, int* neighbors_count, int** send_count, int*** send_list,
        int** recv_count, int*** recv_list, int** epart, int** npart, int* objval,
        int* num_global_elem);

#endif /* INITIALIZATION_H_ */

#include <stdio.h>

int allocate_local_arrays(double **BS, double **BE, double **BN, double **BW, double **BL,
        double **BH, double **BP, double **SU, int num_elem);

int comm_model(int ne_l, int ne_g, int*** lcc, int** local_global_index, int** global_local_index,
        int* neighbors_count, int** send_count, int*** send_list, int** recv_count,
        int*** recv_list, int** epart);
