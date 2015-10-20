#include <mpi.h>
#include <stdio.h>
#include <sys/time.h>
#include <stdlib.h>
#include <math.h>

void main(int argc, char *argv[])
{
	MPI_Status status;
	MPI_Init(&argc, &argv);
	int iproc, nproc;
	MPI_Comm_size(MPI_COMM_WORLD, &nproc);
	MPI_Comm_rank(MPI_COMM_WORLD, &iproc);

	//MPI_COMM communicators[log2(nproc)];
	MPI_Comm currentComm = MPI_COMM_WORLD;
	int numDim = log2(nproc);

	int count;

	//for(count=0;count<log2(nproc);count++)
	//{
	int commSplitMask;
	commSplitMask = (int)pow(2,numDim);

	int color;
	if(iproc & commSplitMask == 0)
		color = 0;
	else color = 1;
	MPI_Comm newComm;
	MPI_Comm_split(currentComm, color, iproc, &newComm);
	//communicators[count] = newComm;
	currentComm = newComm;
	int currentRank;
	MPI_Comm_rank(currentComm, &currentRank);
	printf("Original Rank: %d, New Rank: %d\n", iproc, currentRank);
	//}

}
