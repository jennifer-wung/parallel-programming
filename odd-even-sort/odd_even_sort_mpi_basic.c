#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include <mpi.h>


int checkIsSorted(int nInts, int* buf)
{
	int i;
	for (i = 1; i < nInts; i++)
    		if (buf[i-1] > buf[i])
    			return 0;
	return 1;
}

void swap(int* buf, int i, int j)
{
	int c;
	c = buf[i];
	buf[i] = buf[j];
	buf[j] = c;
}

seqOddEvenSort(int* buf, int nbufInts)
{
	bool sorted = false;
	int i;
	while (!sorted) {
		sorted = true;
		for (i=1;i < nbufInts-1;i+=2) {
			if (buf[i] > buf[i+1]) {swap(buf,i,i+1); sorted = false;}
		}
		for (i=0;i<nbufInts-1;i+=2){
			if (buf[i] > buf[i+1]) {swap(buf,i,i+1); sorted = false;}
		}
	}

}

int main(int argc, char **argv)
{


	/* Initialize MPI and get system information */
	int nprocs, myrank;
	int ierr;
	int *buf, *bcastBuf, recvNum, nints, nbufInts;
	int nints_final;
	bool remainFlag = false;
	bool lessElemts = false;
	int  isSorted = 0;
	int remain = 0;

	clock_t IO_time, commu_time, compute_time, begin;
	double ioSpent_time = 0.0, computeSpent_time = 0.0, commuSpent_time = 0.0;
	double total_IO = 0.0, total_compute = 0.0, total_commu = 0.0;	

	MPI_File fh;
	MPI_Offset filesize;
    	MPI_Offset start;
    	MPI_Offset end;
	MPI_Status status;
	MPI_Init(&argc, &argv); 
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs); 
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	
	/**********************************read file using MPI-IO*************************************/

	if (argc != 4) {
        	if (myrank == 0) 
			fprintf(stderr, "Usage: size infilename outfinename.");
        	MPI_Finalize();
        	exit(1);
    	}

	IO_time = clock();
	/* Open the file */
	ierr = MPI_File_open(MPI_COMM_WORLD,argv[2], MPI_MODE_RDONLY, MPI_INFO_NULL,&fh);
	if (ierr) {
        	if (myrank == 0) 
			fprintf(stderr, "Couldn't open file %s\n", argv[2]);
        	MPI_Finalize();
        	exit(2);
    	}

	/* Get the size of file */
	MPI_File_get_size(fh, &filesize);
	/* Calculate how many elements in total and how many elements in each buf*/
	nints = filesize/sizeof(int);
	if ( nints % nprocs == 0 ) {
		nbufInts = nints/nprocs;
		nints_final = nints;
	}
	else {  
		if (nints >= nprocs) {
			nbufInts = nints/nprocs+1;//result is always floor
			remain = nints%nbufInts;
			remainFlag = true;
		} 
		else {
			nbufInts = 1;
			lessElemts = true;	
		}
		nints_final = nbufInts*nprocs;
		
	}
	/* Calculate starting and ending point of each buf in the file */
	
	start = myrank * nbufInts*sizeof(int);
	//printf("nints:%d; nbufInts:%d\n",nints,nbufInts);
	/* Allocate the buffer to read to, one extra for terminating null int */
	buf = (int *) malloc((nbufInts)*sizeof(int));
	/* everyone reads in their part */
    	MPI_File_read_at_all(fh, start, buf, nbufInts, MPI_INT, MPI_STATUS_IGNORE);
	/*to fill the blanks that doesn't have values from read in file*/
	if (myrank == nprocs-1 && remainFlag == true) {
		int ii; 
		for (ii = nbufInts-1; ii >= remain; ii--)
			buf[ii] = -1;
	}
	else if (lessElemts == true){
		if (myrank >= nints)
			buf[0] = -1;
	}
	MPI_File_close(&fh);

	ioSpent_time = ioSpent_time + ((double)clock() - IO_time);
	//printf("[%d]first time for IO  %f ms\n\n",myrank, ioSpent_time / (CLOCKS_PER_SEC/1000));
	
	/***********************************Odd-even sorting**********************************/
	int oddrank, evenrank;

	compute_time = clock();

	/* Sort the local elements using sequential odd even sort*/
	seqOddEvenSort(buf, nbufInts);

	computeSpent_time = computeSpent_time + ((double)clock() - compute_time);
	
	
	/*Determine odd and even phases for each rank*/
	if (myrank%2 == 0) { 
		oddrank  = myrank-1; //(odd,even)-indexed phase
		evenrank = myrank+1; //(even,odd)-indexed phase
	} 
	else { 
		oddrank  = myrank+1; 
		evenrank = myrank-1; // even phase with process on the left hand side
	} 
	/* Set the ranks of the processors at the end of the linear */ 
	if (oddrank == -1 || oddrank == nprocs) 
		oddrank = MPI_PROC_NULL;
	if (evenrank == -1 || evenrank == nprocs)
		evenrank = MPI_PROC_NULL;
	int i;
	int sendNum;
	int iter;
	nprocs < nbufInts? (iter = nbufInts): (iter = nprocs);

	while(!isSorted) {
	   for (i=0; i < iter; i++) {
		if(i%2 == 1) { // (odd,even)-indexed phase
			if (oddrank != MPI_PROC_NULL) {	
				if (myrank%2 == 0) sendNum = buf[0]; else sendNum = buf[nbufInts-1];
				commu_time = clock();
				MPI_Sendrecv(&sendNum, 1, MPI_INT, oddrank, 1, &recvNum, 1, MPI_INT, oddrank, 1, MPI_COMM_WORLD, &status);
				//MPI_Send(&sendNum, 1, MPI_INT, oddrank, 0, MPI_COMM_WORLD);
				//MPI_Recv(&recvNum, 1, MPI_INT, oddrank, 0, MPI_COMM_WORLD, &status);

				commuSpent_time = commuSpent_time + ((double)clock() - commu_time);

				if (myrank%2 == 0 && recvNum > buf[0]) 
					buf[0] = recvNum; 
				else if(myrank%2 !=0 && recvNum < buf[nbufInts-1]) 
					buf[nbufInts-1] = recvNum;

				compute_time = clock();
				seqOddEvenSort(buf, nbufInts);
			        computeSpent_time = computeSpent_time + ((double)clock() - compute_time);
			}
		}
		else { // (even,odd)-indexed phase 
			if (evenrank != MPI_PROC_NULL) {
				if (myrank%2 == 0) sendNum = buf[nbufInts-1]; else sendNum = buf[0];

				commu_time = clock();
				//MPI_Send(&sendNum, 1, MPI_INT, evenrank, 0, MPI_COMM_WORLD); // rank is odd(even) send to left(right)
                        	//MPI_Recv(&recvNum, 1, MPI_INT, evenrank, 0, MPI_COMM_WORLD, &status); // recieve from the rank that just sent
                        	MPI_Sendrecv(&sendNum, 1, MPI_INT, evenrank, 1, &recvNum, 1, MPI_INT, evenrank, 1, MPI_COMM_WORLD, &status);
				commuSpent_time = commuSpent_time + ((double)clock() - commu_time);

				if (myrank%2 == 0 && recvNum < buf[nbufInts-1]) 
					buf[nbufInts-1] = recvNum; 
				else if (myrank%2 !=0 && recvNum > buf[0])
					buf[0] = recvNum;

				compute_time = clock();
				seqOddEvenSort(buf, nbufInts);
			        computeSpent_time = computeSpent_time + ((double)clock() - compute_time);
			}
		}
	   }
  	   MPI_Barrier(MPI_COMM_WORLD);
	   int *arrayToSort = (int *)malloc((nints_final)*sizeof(int));
	   MPI_Gather(buf,nbufInts,MPI_INT,arrayToSort,nbufInts,MPI_INT,0,MPI_COMM_WORLD);
	   isSorted = checkIsSorted(nints_final,arrayToSort);
	   MPI_Bcast(&isSorted, 1, MPI_INT, 0, MPI_COMM_WORLD);
  	   MPI_Barrier(MPI_COMM_WORLD);
           if(myrank == 0)
	   	printf("[%d] isSorted? %s\n", myrank, isSorted ? "true" : "false");
	}
	
	//printf("[%d]Time for compute  %f ms\n\n",myrank, computeSpent_time / (CLOCKS_PER_SEC/1000));
	//printf("[%d]Time for communication  %f ms\n\n",myrank, commuSpent_time / (CLOCKS_PER_SEC/1000));

	//int *arrayToSort = (int *)malloc((nints_final)*sizeof(int));
	//MPI_Gather(buf,nbufInts,MPI_INT,arrayToSort,nbufInts,MPI_INT,0,MPI_COMM_WORLD);
	
	//if(myrank == 0) {
	//	int i;
	//	int counter = 0;
 	//	//printf("\nSorted array :\n");
 	//	for (i=0;i<nints_final;i++) {
	//		counter++;
	//	//	printf("%d\n",arrayToSort[i]); 
	//	} 
	//	printf("number of sorted numbers: %d\n",counter);
        ////printf("sorted in %f s\n\n",((double)clock() - start) / CLOCKS_PER_SEC);
	//}

	/************************************Write file using MPI-IO************************************************/
	IO_time = clock();

	ierr = MPI_File_open(MPI_COMM_WORLD,argv[3], MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL,&fh);
        if (ierr) {
                if (myrank == 0)
                        fprintf(stderr, "Couldn't open file %s\n", argv[3]);
                MPI_Finalize();
                exit(3);
        }
	if ( nints%nprocs != 0 && nints >= nprocs) {
		if (myrank == 0) {
			int *temp = (int *) malloc((nbufInts)*sizeof(int));
			temp = buf;
			buf = (int *) malloc((remain)*sizeof(int));
			int j;
			for (j = 0; j < remain; j++)
				buf[j] = temp[j+(nbufInts-remain)];
			nbufInts = remain;
		}
		else {
			start = start - (nbufInts-remain)*sizeof(int);
		}
	}
	else if (lessElemts == true) {
		if (myrank < nprocs - nints)
			start = -1;
		else 
			start = start - (nprocs - nints)*sizeof(int);
	}
	MPI_File_write_at(fh, start, buf, nbufInts, MPI_INT, &status);
	MPI_File_close(&fh);

	ioSpent_time = ioSpent_time + ((double)clock()-IO_time);
	
	MPI_Reduce (&ioSpent_time, &total_IO, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce (&computeSpent_time, &total_compute, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce (&commuSpent_time, &total_commu, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	if (myrank == 0) {
		printf("Total time for I/O  %f ms\n", (total_IO/nprocs)/(CLOCKS_PER_SEC/1000));
		printf("Total time for compute  %f ms\n", (total_compute/nprocs)/(CLOCKS_PER_SEC/1000));
		printf("Total time for communication  %f ms\n", (total_commu/nprocs)/(CLOCKS_PER_SEC/1000));
	}

	free(buf);
	MPI_Finalize();
	return 0;
	
}



