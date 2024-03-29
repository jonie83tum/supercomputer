/*
 * compute_solution.h
 *
 *  Created on: Oct 21, 2012
 *      Author: petkovve
 */

#ifndef COMPUTE_SOLUTION_H_
#define COMPUTE_SOLUTION_H_

int compute_solution(const int max_iters, int nintci, int nintcf, int nextci, int nextcf, int** lcc,
        double* bp, double* bs, double* bw, double* bl, double* bn, double* be, double* bh,
        double* cnorm, double* var, double *su, double* cgup, double* residual_ratio,
        int* local_global_index, int* global_local_index, int neighbors_count, int* send_count,
        int** send_list, int* recv_count, int** recv_list, int num_global_elem, int num_all_elem);

#endif /* COMPUTE_SOLUTION_H_ */

