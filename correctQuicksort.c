// quickSort.c
#include <mpi.h>
#include <stdio.h>
#include <sys/time.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>


#define VECSIZE 20
#define MAXVALUE 1000

typedef enum { false, true } bool;

void quickSort( int[], int, int);
int partition( int[], int, int);
int parallelPartition( int[], int, int, int);
int findPivot(int[], MPI_Comm, int);
int* sendCollection(int[], int, MPI_Comm,int*, bool);

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
	int *collectionPointer = numberCollection;
	
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
	
	fprintf(stderr,"Rank %d generated list\n", iproc);
	//for(count = 0; count < VECSIZE; count++)
	//{
	//	fprintf(stderr,"My rank is %d. Number[%d] == %d\n",iproc, count, numberCollection[count]);
	//}	

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
		//printf("My rank is %d.  For communicator %d, the Mask is %d, my rank is %d, my result is %d, and my color is %d.\n",iproc,count,commSplitMask,commRank,bitwiseResult, color);
		MPI_Comm newComm;
		MPI_Comm_split(currentComm, color, iproc, &newComm);
		communicators[count] = newComm;
		currentComm = newComm;
		numDim -= 1;
	}
	
	fprintf(stderr,"Rank %d created communicators\n",iproc);

	numDim = log2(nproc);
	int collectionSize = VECSIZE;
	int *sizePointer = &collectionSize;
	
	for(count = 0; count<numDim; count++)
	{
		int pivot = findPivot(numberCollection, communicators[count], *sizePointer);
		int pivotLocation = parallelPartition(collectionPointer, 0, *sizePointer-1, pivot);
		int currentRank;
		MPI_Comm_rank(communicators[count], &currentRank);
		int currentProcesses;
		MPI_Comm_size(communicators[count], &currentProcesses);
		bool sendingHigh;
		if(currentRank >= (currentProcesses/2))
			sendingHigh = false;
		else sendingHigh = true;
		collectionPointer = sendCollection(collectionPointer, pivotLocation, communicators[count], sizePointer, sendingHigh);
	}
	fprintf(stderr,"Rank %d partitioned and traded numbers\n", iproc);
	quickSort( collectionPointer, 0, *sizePointer);
	fprintf(stderr,"Rank %d sorted its list\n",iproc);
	printf("My rank is %d and my size is %d\n",iproc, *sizePointer);	

	//for(count = 0; count < *sizePointer; count++)
	//{
	//	printf("My rank is %d. Number[%d] == %d\n",iproc, count, collectionPointer[count]);
	//}	

	MPI_Finalize();
	
}



void quickSort( int numberCollection[], int l, int r)
{
   int j;

   if( l < r ) 
   {
       j = originalPartition( numberCollection, l, r);
       quickSort( numberCollection, l, j-1);
       quickSort( numberCollection, j+1, r);
   }
	
}



int originalPartition( int numberCollection[], int l, int r) 
{
   	int pivot, i, j, t;
   	pivot = numberCollection[l];
   	i = l; j = r;
		
   	while( 1)
   	{
		while( numberCollection[i] <= pivot && i <= r )
			++i;
		while( numberCollection[j] > pivot )
		{
			//printf("j = %d\n",j);
			//printf("numberCollection[j] = %d\n", numberCollection[j]);
			--j;
		}
		if( i >= j ) break;
	   	t = numberCollection[i]; numberCollection[i] = numberCollection[j]; numberCollection[j] = t;
		//printf("Current array is:  ");
		//for(i = 0; i < 9; ++i)
			//printf(" %d ", numberCollection[i]);
		//printf("\n");
		//printf("j = %d\n",j);
   	}
   	t = numberCollection[l]; numberCollection[l] = numberCollection[j]; numberCollection[j] = t;
   	return j;
}


int parallelPartition( int *numberCollection, int l, int collectionSize, int pivot) 
{
   	int t;
	int left = 0;
	int right = collectionSize-1;
   	while( 1)
   	{
		while( numberCollection[left] <= pivot && left <= right )
			left++;
		while( numberCollection[right] > pivot )
			right--;
		if( left >= right ) 
			{
				if( left > right )
					right++;
				break;
			}
	   	t = numberCollection[left];
		numberCollection[left] = numberCollection[right]; 
		numberCollection[right] = t;
   	}
   	t = numberCollection[l]; 
	numberCollection[l] = numberCollection[right]; 
	numberCollection[right] = t;
	return right;
}


int findPivot(int numberCollection[], MPI_Comm comm, int collectionSize)
{
	int possiblePivots = (int)log2(collectionSize);
	int pivots[possiblePivots];
	int i;
	for(i = 0; i < possiblePivots; i++)
	{
		int randNum = rand() % (collectionSize);
		pivots[i] = numberCollection[randNum];
	}

	int localSum = 0;
	for(i = 0; i < possiblePivots; i++)
	{
		localSum += pivots[i];
	}

	int globalSum;
	MPI_Allreduce(&localSum, &globalSum, 1, MPI_INT, MPI_SUM, comm);
	int globalPivots;
	MPI_Allreduce(&possiblePivots, &globalPivots, 1, MPI_INT, MPI_SUM, comm);

	int currentProcessors;
	MPI_Comm_size(MPI_COMM_WORLD, &currentProcessors);
	return globalSum/globalPivots;
	
}

int* sendCollection(int collection[], int pivotLocation, MPI_Comm comm, int *collectionSize, bool sendingHigh)
{
	int* sentCollection;
	int sentSize;
	int count;
	if(sendingHigh)
	{
		sentSize = *collectionSize-pivotLocation;
		int sentArray[sentSize];
		for(count = 0; count < sentSize; count++)
		{
			sentArray[count] = collection[count + pivotLocation];
		}
		sentCollection = sentArray;
	}
	else
	{
		int sentArray[pivotLocation];
		sentSize = pivotLocation;
		for(count = 0; count < pivotLocation; count++)
		{
			sentArray[count] = collection[count];
		}
		sentCollection = sentArray;
	}
	
	int currentRank;
	MPI_Comm_rank(comm, &currentRank);
	int msg_dest;
	int currentProcesses;
	MPI_Comm_size(comm, &currentProcesses);
	if(sendingHigh)
		msg_dest = currentRank + currentProcesses/2;
	else msg_dest = currentRank - currentProcesses/2;

	MPI_Send(&sentSize, 1, MPI_INT, msg_dest, 0, comm);
	
	int sizeReceiving;
	MPI_Recv(&sizeReceiving, 1, MPI_INT, msg_dest, 0, comm, MPI_STATUS_IGNORE);
	
	MPI_Send(&sentCollection, sentSize, MPI_FLOAT, msg_dest, 0, comm);

	int receivedNumbers[sizeReceiving];
	MPI_Recv(&receivedNumbers, sizeReceiving, MPI_FLOAT, msg_dest, 0, comm, MPI_STATUS_IGNORE);
	
	*collectionSize -= sentSize;
	*collectionSize += sizeReceiving;

	int *newCollectionPtr;
	int newCollection[*collectionSize];
	if(sendingHigh)
	{
		for(count = 0; count < pivotLocation; count++)
		{
			newCollection[count] = collection[count];
		}
		for(count = 0; count < sizeReceiving; count++)
		{
			newCollection[count + pivotLocation] = receivedNumbers[count];
		}
	}
	else
	{
		for(count = 0; count < sizeReceiving; count++)
		{
			newCollection[count] = receivedNumbers[count];
		}
		for(count = 0; count < pivotLocation; count++)
		{
			newCollection[count + sizeReceiving] = collection[count];
		}
	}
	newCollectionPtr = newCollection;
	return newCollectionPtr;
}






