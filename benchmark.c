#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <mpi.h>
#include <math.h>
#define VECSIZE 65536
#define ITERATIONS 50

double When()
{
  struct timeval tp;
  gettimeofday(&tp, NULL);
  return ((double) tp.tv_sec + (double) tp.tv_usec * 1e-6);
}

double* myReduce(int numdim, int rank, double oldValues[], int nproc);
double* myBroadcast(int numdim, int rank, double oldValues[], int nproc);

main(int argc, char *argv[])
{
	int iproc, nproc,i, iter;
	char host[255], message[55];
	MPI_Status status;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &nproc);
	MPI_Comm_rank(MPI_COMM_WORLD, &iproc);

	gethostname(host,253);
	//printf("I am proc %d of %d running on %s\n", iproc, nproc,host);
	// each process has an array of VECSIZE double: ain[VECSIZE]
	double ain[VECSIZE], aout[VECSIZE];
	int  ind[VECSIZE];
	struct {
		double val;
		int   rank;
	} in[VECSIZE], out[VECSIZE];
	int myrank, root = 0;

	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	// Start time here
	srand(myrank+5);
	double start = When();
	for(iter = 0; iter < ITERATIONS; iter++) 
	{
		for(i = 0; i < VECSIZE; i++) 
		{
			ain[i] = rand();
			//printf("init proc %d [%d]=%f\n",myrank,i,ain[i]);
		}
		//printf("Vector Values: \n");
		for (i=0; i<VECSIZE; ++i) 
		{
			in[i].val = ain[i];
			in[i].rank = myrank;
			//printf("\t%f\n",in[i].val);
		}
		int numdim = log(nproc)/log(2) + 0.99;
		//printf("number of dimensions:%d\n",numdim);
		double* output;
		output  = myReduce(numdim, myrank, ain, nproc);
		//MPI_Reduce( in, out, VECSIZE, MPI_DOUBLE_INT, MPI_MAXLOC, root, MPI_COMM_WORLD);
		//At this point, the answer resides on process root
		//if (myrank == root) 
		//{
		/* read ranks out
		*/
			//for (i=0; i<VECSIZE; ++i) 
			//{
				//printf("root out[%d] = %f\n",i,output[i]);
			//}
		//}
		// Now broadcast this max vector to everyone else.
		//MPI_Bcast(out, VECSIZE, MPI_DOUBLE_INT, root, MPI_COMM_WORLD);
		myBroadcast(numdim, myrank,output, nproc);
	}
	MPI_Finalize();
	double end = When();
	if(myrank == root) 
	{
		printf("Time %f\n",end-start);
	}
}

double* myReduce(int numdim, int rank, double oldValues[], int nproc)
{
	//printf("Start Reduction\n");
	int notParticipating = 0;
	int bitmask = 1;
	double values[VECSIZE];
	int i;
	for(i=0; i < VECSIZE; i++){
		values[i] = oldValues[i];
	}
	//double max[] = values;
	double newValues[VECSIZE];
	int curdim;
	for(curdim = 0; curdim < numdim; curdim++)
	{
		int participating = rank & notParticipating;
		if((rank & notParticipating) == 0)
		{
			int bit = rank & bitmask;
			if((rank & bitmask) != 0)
			{
				int msg_dest = rank ^ bitmask;
				MPI_Send(&values, VECSIZE, MPI_DOUBLE, msg_dest, 0, MPI_COMM_WORLD);
			} 
			else 
			{
				int msg_src = rank ^ bitmask;
				//printf("Message source: %d\n", msg_src);
				if(msg_src < nproc)
				{
					//printf("Trying to receive\n");
					MPI_Recv(&newValues, VECSIZE, MPI_DOUBLE, msg_src, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
					int i;
					for(i = 0; i < VECSIZE; i++)
					{
						if(newValues[i] > values[i])
						{
							values[i] = newValues[i];
						}
					}
				}
			}
		}
		notParticipating = notParticipating ^ bitmask;
		bitmask <<= 1;
	}
	static double returnValues[VECSIZE];
	for(i = 0; i < VECSIZE; i++)
	{
		returnValues[i] = values[i];
	}
	return returnValues;
}

double* myBroadcast(int numdim, int rank, double oldValues[], int nproc)
{
	int notparticipating = pow(2,numdim-1)-1;
    int bitmask = pow(2,numdim-1);
    double values[VECSIZE];
    int i;
	for(i=0; i < VECSIZE; i++){
		values[i] = oldValues[i];
	}
	double newValues[VECSIZE];
	int curdim;
    for(curdim = 0; curdim < numdim; curdim++) 
    {
		if ((rank & notparticipating) == 0) 
		{
			int bit = rank & bitmask;
			if ((rank & bitmask) == 0) {
				int msg_dest = rank ^ bitmask;
				if(msg_dest < nproc)
				{
					MPI_Send(&values, VECSIZE, MPI_DOUBLE, msg_dest, 0, MPI_COMM_WORLD);
				}
			} 
			else 
			{
				int msg_src = rank ^ bitmask;
				MPI_Recv(&newValues, VECSIZE, MPI_DOUBLE, msg_src, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				int j;
				for(j = 0; j < VECSIZE; j++)
				{
					values[j] = newValues[j];
				}
			}
		}
		notparticipating >>= 1;
		bitmask >>=1;
	}
	
}

