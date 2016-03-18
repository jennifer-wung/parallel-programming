#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include <mpi.h>


/* The IncOrder function that is called by qsort is defined as follows */ 
int cmpfunc(const void *e1, const void *e2) 
{ 
	return (*((int*)e1) - *((int*)e2)); 
}

compareSplit(int nbufInts, int* buf, int* recv_buf, int keepsmall)
{
	int i, j;
	int* temp = (int *)malloc(nbufInts*2*sizeof(int));
	for (i=0; i<nbufInts; i++){	
		temp[i] = buf[i];
		temp[i+nbufInts] = recv_buf[i];
	}
	qsort(temp, nbufInts*2, sizeof(int), cmpfunc);
	if (keepsmall) {
		for (i=0; i<nbufInts; i++)
			buf[i] = temp[i];
	} 
	else {
		for (i=nbufInts,j=0; i<nbufInts*2; i++,j++)
                        buf[j] = temp[i];
	}
	free(temp);
}

int main(int argc, char **argv)
{


	/* Initialize MPI and get system information */
	int nprocs, myrank;
	int ierr;
	int *buf, *recv_buf, nints, nbufInts;
	int nints_final;
	bool remainFlag = false;
	bool lessElemts = false;
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
	/* Open the file */
	
	IO_time = clock();

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
    	end   = start + nbufInts*sizeof(int) - 1;
	/* Allocate the buffer to read to, one extra for terminating null int */
	buf = (int *) malloc((nbufInts)*sizeof(int));
	/* everyone reads in their part */
    	MPI_File_read_at_all(fh, start, buf, nbufInts, MPI_INT, MPI_STATUS_IGNORE);
	/*to fill the blanks that doesn't have values from read in file*/
	if (myrank == nprocs-1 && remainFlag == true) {
		int iter; 
		for (iter = nbufInts-1; iter >= remain; iter--)
			buf[iter] = -1;
	}
	else if (lessElemts == true){
		if (myrank+1 > nints)	
			buf[0] = -1;
	}
	MPI_File_close(&fh);
	
	ioSpent_time = ioSpent_time + ((double)clock() - IO_time);

	/***********************************Odd-even sorting**********************************/
	int oddrank, evenrank;
	/* Sort the local elements using the built-in quicksort routine */
	compute_time = clock();

	qsort(buf, nbufInts, sizeof(int), cmpfunc);

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
	for (i=0; i < nprocs; i++) {
		recv_buf = (int *) malloc((nbufInts)*sizeof(int));
		if(i%2 == 1) {// (odd,even)-indexed phase
			if (oddrank != MPI_PROC_NULL) {	
				commu_time = clock();
				//MPI_Send(buf, nbufInts, MPI_INT, oddrank, 0, MPI_COMM_WORLD);
				//MPI_Recv(recv_buf, nbufInts, MPI_INT, oddrank, 0, MPI_COMM_WORLD, &status);
				MPI_Sendrecv(buf, nbufInts, MPI_INT, oddrank, 1, recv_buf, nbufInts, MPI_INT, oddrank, 1, MPI_COMM_WORLD, &status);		
				commuSpent_time = commuSpent_time + ((double)clock() - commu_time);
				

				compute_time = clock();
      				compareSplit(nbufInts, buf, recv_buf, myrank < status.MPI_SOURCE);
			        computeSpent_time = computeSpent_time + ((double)clock() - compute_time);
			}
		}
		else {// (even,odd)-indexed phase 
			if (evenrank != MPI_PROC_NULL) {
				commu_time = clock();
				
				//MPI_Send(buf, nbufInts, MPI_INT, evenrank, 0, MPI_COMM_WORLD); // rank is odd(even) send to left(right)
                        	//MPI_Recv(recv_buf, nbufInts, MPI_INT, evenrank, 0, MPI_COMM_WORLD, &status); // recieve from the rank that just sent 
				MPI_Sendrecv(buf, nbufInts, MPI_INT, evenrank, 1, recv_buf, nbufInts, MPI_INT, evenrank, 1, MPI_COMM_WORLD, &status);		
				commuSpent_time = commuSpent_time + ((double)clock() - commu_time);
				compute_time = clock();
      				compareSplit(nbufInts, buf, recv_buf, myrank < status.MPI_SOURCE);
			        computeSpent_time = computeSpent_time + ((double)clock()- compute_time);
			}
		}
		free(recv_buf);
	}
	//printf("[%d]Time for compute  %f ms\n\n", myrank,computeSpent_time / (CLOCKS_PER_SEC/1000));
	//printf("[%d]Time for communication  %f ms\n\n", myrank,commuSpent_time / (CLOCKS_PER_SEC/1000));
	//int *arrayToSort = (int *)malloc((nints_final)*sizeof(int));
	//MPI_Gather(buf,nbufInts,MPI_INT,arrayToSort,nbufInts,MPI_INT,0,MPI_COMM_WORLD);
	
	//if(myrank == 0) {
	//	//int i;
	//	int counter = 0;
	//	//printf("\nSorted array :\n");
 	//	for (i=0;i<nints_final;i++) {
	//		counter++;
	//		//printf("%d\n",arrayToSort[i]); 
	//	} 
	//	printf("number of sorted numbers: %d\n",counter);
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

	if (nints % nprocs != 0 && nints >= nprocs) {
		if (myrank == 0) {
                	int *temp = (int *) malloc((nbufInts)*sizeof(int));
                	temp = buf;
			free(buf);
                	buf = (int *) malloc(remain*sizeof(int));
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
	//printf("Time for I/O  %f ms\n\n", ioSpent_time/(CLOCKS_PER_SEC/1000));
	free(buf);
	MPI_Finalize();
	//printf("rank: %d --> sorted in %f ms\n\n",myrank, ((double)clock() - begin)/(CLOCKS_PER_SEC/1000));
	return 0;
	
}
