// quickSort.c
#include <stdio.h>

void quickSort( int[], int, int);
int partition( int[], int, int);


void main() 
{
	int a[] = { 7, 12, 1, -2, 0, 15, 4, 11, 9};

	int i;
	printf("Unsorted array is:  ");
	for(i = 0; i < 9; ++i)
		printf(" %d ", a[i]);
	printf("\n");

	quickSort( a, 0, 8);

	printf("\n\nSorted array is:  ");
	for(i = 0; i < 9; ++i)
		printf(" %d ", a[i]);
	printf("\n");

}



void quickSort( int a[], int l, int r)
{
   int j;

   if( l < r ) 
   {
       j = partition( a, l, r);
       quickSort( a, l, j-1);
       quickSort( a, j+1, r);
   }
	
}



int partition( int a[], int l, int r) {
   	int pivot, i, j, t;
   	pivot = a[l];
   	i = l; j = r;
		
   	while( 1)
   	{
	   	//do ++i; 
		while( a[i] <= pivot && i <= r )
			++i;
	   	//do --j; 
		while( a[j] > pivot )
		{
			printf("j = %d\n",j);
			printf("a[j] = %d\n", a[j]);
			--j;
		}
		if( i >= j ) break;
	   	t = a[i]; a[i] = a[j]; a[j] = t;
		printf("Current array is:  ");
		for(i = 0; i < 9; ++i)
			printf(" %d ", a[i]);
		printf("\n");
		printf("j = %d\n",j);
   	}
   	t = a[l]; a[l] = a[j]; a[j] = t;
   	return j;
}








