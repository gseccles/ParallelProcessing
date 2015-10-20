// quickSort.c
#include <mpi.h>
#include <stdio.h>
#include <sys/time.h>
#include <stdlib.h>
#include <math.h>

#define VECSIZE 1000
#define MAXVALUE 1000

void quickSort( int[], int, int);
int partition( int[], int, int);
double When()
{
	struct timeval tp;
	gettimeofday(&tp, NULL);
	return ((double) tp.tv_sec + (double) tp.tv_usec * 1e-6);
}

void main(int argc, char *argv[]) 
{
	srand(time(NULL));
	int numberCollection[VECSIZE];
	
	int iproc, nproc,i, iter;
	MPI_Status status;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &nproc);
	MPI_Comm_rank(MPI_COMM_WORLD, &iproc);
	MPI_Comm currentComm = MPI_COMM_WORLD;

	
	MPI_Comm communicators[(int)log2(nproc)];

	for(i = 0; i < VECSIZE; i++)
	{
		numberCollection[i] = rand() % 1000;
	}

	int numDim = log2(nproc);

	int count;

	for(count=0;count<log2(nproc);count++)
	{
		int commSplitMask;
		commSplitMask = (int)pow(2,numDim);

		int commRank;
		MPI_Comm_rank(currentComm, &commRank);

		int color;
		int bitwiseResult = commRank & commSplitMask;
		if(bitwiseResult == 0)
			color = 0;
		else color = 1;
		printf("My rank is %d.  For communicator %d, the Mask is %d, my rank is %d, my result is %d, and my color is %d.\n",iproc,count,commSplitMask,commRank,bitwiseResult, color);
		MPI_Comm newComm;
		MPI_Comm_split(currentComm, color, iproc, &newComm);
		communicators[count] = newComm;
		currentComm = newComm;
		numDim -= 1;
	}
	for(i = 0; i < log2(nproc);i++)
	{
		int commRank;
		MPI_Comm_rank(communicators[i], &commRank);
		//printf("My rank is %d.  For communicator %d, my rank is %d.\n",iproc,i,commRank);
	}

	//quickSort( numberCollection, 0, 8);
	MPI_Finalize();
	
}



void quickSort( int numberCollection[], int l, int r)
{
   int j;

   if( l < r ) 
   {
       j = partition( numberCollection, l, r);
       quickSort( numberCollection, l, j-1);
       quickSort( numberCollection, j+1, r);
   }
	
}



int partition( int numberCollection[], int l, int r) {
   	int pivot, i, j, t;
   	pivot = numberCollection[l];
   	i = l; j = r;
		
   	while( 1)
   	{
		while( numberCollection[i] <= pivot && i <= r )
			++i;
		while( numberCollection[j] > pivot )
		{
			printf("j = %d\n",j);
			printf("numberCollection[j] = %d\n", numberCollection[j]);
			--j;
		}
		if( i >= j ) break;
	   	t = numberCollection[i]; numberCollection[i] = numberCollection[j]; numberCollection[j] = t;
		printf("Current array is:  ");
		for(i = 0; i < 9; ++i)
			printf(" %d ", numberCollection[i]);
		printf("\n");
		printf("j = %d\n",j);
   	}
   	t = numberCollection[l]; numberCollection[l] = numberCollection[j]; numberCollection[j] = t;
   	return j;
}








