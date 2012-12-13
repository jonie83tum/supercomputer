#include "mpi.h"

int main(int argc, char *argv[]) {

    int *send_count, recv_count;
    int **send_list, **recv_list;
    int i, j;
    int my_rank, num_procs;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

    send_count = (int*) calloc(sizeof(int), num_procs);
    recv_count = (int*) calloc(sizeof(int), num_procs);

// some arbitrary numbers for send_count (num_procs=6)

// allocate memeory for the send_list
    send_list = (int**) malloc(num_procs * sizeof(int*));
    for (i = 0; i < num_procs; i++) {
        send_list[i] = (int *) malloc(send_count[i] * sizeof(int));
    }

// fill the send_list somehow
    for (i = 0; i < num_procs; i++) {
        for (j = 0; j < send_count[i]; j++) {
            send_list[i][j] = 1;
        }
    }

// send_count and recv_count are the same in this case
    for (i = 0; i < num_procs; i++) {
        recv_count[i] = send_count[i];
    }

// allocate memeory for the recv_list
    recv_list = (int**) malloc(num_procs * sizeof(int*));
    for (i = 0; i < num_procs; i++) {
        recv_list[i] = (int *) malloc(send_count[i] * sizeof(int));
    }

    MPI_Request req_send, req_recv;
    for (i = 0; i < num_procs; i++) {
        if (my_rank != i) {
            MPI_Isend(send_list[i], send_count[i], MPI_INT, i, 1, MPI_COMM_WORLD, &req_send);
        }
    }
    for (i = 0; i < num_procs; i++) {
        if (my_rank != i) {
            MPI_Irecv(recv_list[i], recv_count[i], MPI_INT, i, 1, MPI_COMM_WORLD, &req_recv);
        }
    }
    MPI_Wait(&req_send, MPI_STATUS_IGNORE);
    MPI_Wait(&req_recv, MPI_STATUS_IGNORE);

    // check the content of the recv_list
    for (i = 0; i < num_procs; i++) {
        for (j = 0; j < recv_count[i]; j++) {
            printf("recv_list[%d][%d]=%d", i, j, recv_list[i][j]);
        }
    }

}
