// quickSort.c
#include <mpi.h>
#include <stdio.h>
#include <sys/time.h>
#include <stdlib.h>
#include <math.h>

#define VECSIZE 1000
#define MAXVALUE 1000

void quickSort( float[], int, int, int, int);
int partition( float[], int, int, int);

double When()
{
	struct timeval tp;
	gettimeofday(&tp, NULL);
	return ((double) tp.tv_sec + (double) tp.tv_usec * 1e-6);
}

void main() 
{
	srand(time(NULL));
	float numberCollection[VECSIZE];
	
	int iproc, nproc,i, iter;
	MPI_Status status;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &nproc);
	MPI_Comm_rank(MPI_COMM_WORLD, &iproc);
	currentComm = MPI_COMM_WORLD;

	int i;
	for(i = 0; i < VECSIZE; i++)
	{
		numberCollection[i] = (float)rand()/(float)(RAND_MAX/MAXVALUE);
	}
	//printf("\n\nUnsorted array is:  ");
	//for(i = 0; i < VECSIZE; ++i)
	//	printf(" %f\n", numberCollection[i]);
	int numDim = log2(nproc);
	
	double start = When();
	int count;
	for(count = 0; count < (int)log2(nproc); count++)
	{
		partition(numberCollection, 0, VECSIZE-1
	}
	quickSort( numberCollection, 0, VECSIZE-1, numDim, 1);

	MPI_Finalize();
	double end = When();
	printf("Time %f\n",end-start);

	//printf("\n\nSorted array is:  ");
	//for(i = 0; i < VECSIZE; ++i)
	//	printf(" %f\n", numberCollection[i]);

}


void quickSort( float numberCollection[], int l, int r, int numDim, int count)
{
	int j;

	if( l < r ) 
	{
		j = partition( numberCollection, l, r, numDim, count);
		quickSort( numberCollection, l, j-1, numDim - 1, count + 1);
		quickSort( numberCollection, j+1, r, numDim - 1, count + 1);
	}
	
}

int partition( float numberCollection[], int l, int r, int numDim) {
	int count = 0;
	int i, j;
	float pivot, t;

	int commSplitMask;
	int iproc, nproc;
	MPI_Comm_size(MPI_COMM_WORLD, &nproc);
	MPI_Comm_rank(MPI_COMM_WORLD, &iproc);
	MPI_COMM newComm;
	
	int currentSize = VECSIZE;

	while(numDim > 0)
	{
		commSplitMask = (int)pow(2,numDim)
		int color;
		if(iproc & commSplitMask == 0)
			color = 0;
		else color = 1;
		MPI_Comm_split(currentComm, color, iproc, &newComm);
		currentComm = newComm;	
		
		int possiblePivots = (int)log2(VECSIZE);
		float pivots[possiblePivots];
		for(i = 0; i < possiblePivots; i++)
		{
			int randNum = rand() % (currentSize);
			pivots[i] = numberCollection[randNum];
		}
	
		float localSum = 0;
		for(i = 0; i < possiblePivots; i++)
		{
			localSum += pivots[i];
		}
	
		commSplitMask = (int)pow(2,numDim)
		MPI_COMM newComm;
	
		int color;
		if(iproc & commSplitMask == 0)
			color = 0;
		else color = 1;
		MPI_Comm_split(currentComm, color, iproc, &newComm);
		currentComm = newComm;

		float globalSum;
		MPI_Reduce(&localSum, &globalSum, 1, MPI_FLOAT, MPI_SUM, 0, currentComm);
		pivot = globalSum/(possiblePivots*pow(2,numdim));

		i = l; j = r;
		
		while( 1)
		{
			while( numberCollection[i] <= pivot && i <= r )
				++i
			while( numberCollection[j] > pivot )
				--j;
			if( i >= j ) break;
			t = numberCollection[i]; numberCollection[i] = numberCollection[j]; numberCollection[j] = t;
		}
		t = numberCollection[l]; numberCollection[l] = numberCollection[j]; numberCollection[j] = t;

		i = 0;
		while (numberCollection[i] < pivot)
			i++;
		currentSize -= i;
		float sentNumbers[i];
		int k;
		for(k = 0; k < i; k++)
		{
			sentNumbers[k] = numberCollection[k];
		{
	
		int bitmask = (int)pow(2,numDim-1);
	
		int msg_dest = iproc ^ bitmask;
	
		MPI_Send(&i, 1, MPI_INT, msg_dest, 0, MPI_COMM_WORLD);
	
		MPI_Recv(&sizeReceiving, 1, MPI_INT, msg_dest, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

		currentSize += sizeReceiving;
		float receivedNumbers[sizeReceiving];
		MPI_SEND(&sentNumbers, i, MPI_FLOAT, msg_dest, 0, MPI_COMM_WORLD);

		MPI_Recv(&receivedNumbers, sizeReceiving, MPI_FLOAT, msg_dest, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	
		
	}
	return j;
}








