#include "mpi.h"
#include <stdio.h>

int main(int argc, char *argv[]) {

    int rank, size, i;
    MPI_Datatype type, type2, type_recv;
    int blocklen[3] = {1,1,1};
    int displacement[3] = { 0, 3, 8};
    int buffer[27];
    MPI_Status status;
    int recv[3];

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    if (size < 2) {
        printf("Please run with 2 processes.\n");
        MPI_Finalize();
        return 1;
    }
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    MPI_Type_contiguous(3, MPI_INT, &type2);

    MPI_Type_commit(&type2);
    // MPI_Type_indexed(3, blocklen, displacement, MPI_INT, &type);

    MPI_Type_indexed(3, blocklen, displacement, MPI_INT, &type);
    MPI_Type_commit(&type);


    if (rank == 0) {
        for (i = 0; i < 27; i++)
            buffer[i] = i;
        MPI_Send(buffer, 1, type, 1, 123, MPI_COMM_WORLD);
    }

    if (rank == 1) {
        for (i = 0; i < 27; i++)
            buffer[i] = -1;
        // MPI_Recv(buffer, 1, type, 0, 123, MPI_COMM_WORLD, &status);
        MPI_Recv(recv, 3, MPI_INT, 0, 123, MPI_COMM_WORLD, &status);
        for (i = 0; i < 27; i++)
            printf("buffer[%d] = %d\n", i, buffer[i]);
        for (i=0;i<3;i++){
            printf("recv[%d]=%d\n",i, recv[i]);
        }
        fflush(stdout);
    }

    MPI_Finalize();
    return 0;
}
/*
int rank, size, i;
MPI_Datatype type, type2;
int blocklen[3] = { 2, 3, 1};
int displacement[3] = { 0, 3, 8};
int buffer[27];
MPI_Status status;

MPI_Init(&argc, &argv);
MPI_Comm_size(MPI_COMM_WORLD, &size);
if (size < 2) {
    printf("Please run with 2 processes.\n");
    MPI_Finalize();
    return 1;
}
MPI_Comm_rank(MPI_COMM_WORLD, &rank);

MPI_Type_contiguous(3, MPI_INT, &type2);
MPI_Type_commit(&type2);
MPI_Type_indexed(3, blocklen, displacement, type2, &type);
MPI_Type_commit(&type);

if (rank == 0) {
    for (i = 0; i < 27; i++)
    buffer[i] = i;
    MPI_Send(buffer, 1, type, 1, 123, MPI_COMM_WORLD);
}

if (rank == 1) {
    for (i = 0; i < 27; i++)
    buffer[i] = -1;
    MPI_Recv(buffer, 1, type, 0, 123, MPI_COMM_WORLD, &status);
    for (i = 0; i < 27; i++)
    printf("buffer[%d] = %d\n", i, buffer[i]);
    fflush(stdout);
}

MPI_Finalize();
return 0;
*/
