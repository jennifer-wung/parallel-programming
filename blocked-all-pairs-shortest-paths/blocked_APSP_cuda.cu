#include <stdio.h>
#include <stdlib.h>

#include "cuda.h"
#include "cuda_runtime_api.h"

#define INF 	1000000
#define CMCPYHTD cudaMemcpyHostToDevice
#define CMCPYDTH cudaMemcpyDeviceToHost

#define THREAD_WIDTH 2
#define BLOCK_WIDTH 16

#define HANDLE_ERROR(status) \
{ \
	if (status != cudaSuccess) \
	{ \
		fprintf(stderr, "%s failed  at line %d \nError message: %s \n", \
			__FILE__, __LINE__ ,cudaGetErrorString(status)); \
		exit(EXIT_FAILURE); \
	} \
}

const int V = 7000;
bool gDebug = false;
int n, m;       // Number of vertices, edges
int Dist[V*V];


template <int BLOCK_FACTOR> __global__ void cal_phase_one(const unsigned int Round, const unsigned int n, const size_t pitch, int * const d) 
{
	int newPath;
	
	const int tx = threadIdx.x; 
	const int ty = threadIdx.y;

	// To calculate original index of elements in the block, 0 <= ty, tx < BLOCK_FACTOR
	const int v1 = BLOCK_FACTOR*Round + ty; // block_internal_y 
	const int v2 = BLOCK_FACTOR*Round + tx; // block_internal_x

	const int weightID = v1*pitch + v2;

	__shared__ int primary_d[BLOCK_FACTOR][BLOCK_FACTOR];
	//int primary_d[BLOCK_FACTOR][BLOCK_FACTOR];

	if (v1 < n && v2 < n) {
		primary_d[ty][tx] = d[weightID];
	} else {
		primary_d[ty][tx] = INF;
	}

	// Synchronize to make sure the all value are loaded in block
	__syncthreads();
	
	// For each block, it need to compute B times
	#pragma unroll
	for (int i = 0; i < BLOCK_FACTOR; i++) {
		newPath = primary_d[ty][i] + primary_d[i][tx];
		__syncthreads();
		if (newPath < primary_d[ty][tx]) {
			primary_d[ty][tx] = newPath;
		}
		// Synchronize to make sure that all value are current
		__syncthreads();
	}
	
	if (v1 < n && v2 < n) {
		d[weightID] = primary_d[ty][tx];
	}
}

template <int BLOCK_FACTOR> __global__ void cal_phase_two(const unsigned int Round, const unsigned int n, const size_t pitch, int * const d)
{
	// done calculation in phase 1
	if(blockIdx.x == Round) return;		
	int newPath;

        const int tx = threadIdx.x;
        const int ty = threadIdx.y;

        // To calculate original index of elements in the block, 0 <= ty, tx < BLOCK_FACTOR
	// so each cuda block will have its own pivot-block results calculated from phase 1
        int v1 = BLOCK_FACTOR*Round + ty; // block_internal_y
        int v2 = BLOCK_FACTOR*Round + tx; // block_internal_x

	// Shared varialbes are shared within each B*B submatrice
	__shared__ int primary_d[BLOCK_FACTOR][BLOCK_FACTOR];
	__shared__ int current_d[BLOCK_FACTOR][BLOCK_FACTOR];	
	

	// Pivot-block (The result of pivot-row (pivot-column) blocks depends on pivot block in phase 1 and itself)	
	const int cell_primary = v1 * pitch + v2;
	if (v1 < n && v2 < n) {
		primary_d[ty][tx] = d[cell_primary];
	} else {
                primary_d[ty][tx] = INF;
	}
	
	// Pivot-row and pivot-column blocks
	if (blockIdx.y == 0) { // pivot-row blcoks
		v1 = BLOCK_FACTOR*Round + ty;
		v2 = BLOCK_FACTOR*blockIdx.x + tx;
	} else { // pivot-col blocks
		v1 = BLOCK_FACTOR*blockIdx.x + ty;
		v2 = BLOCK_FACTOR*Round + tx;
	}		

	const int cell_current = v1*pitch + v2;
	if (v1 < n && v2 < n) {
		current_d[ty][tx] = d[cell_current];
	} else {
		current_d[ty][tx] = INF;
	} 
	// Synchronize to make sure the all value are loaded in block
	__syncthreads();

	if (blockIdx.y == 0) { // pivot-row blcoks
		#pragma unroll
		for (int i = 0; i < BLOCK_FACTOR; i++) {
			newPath = primary_d[ty][i] + current_d[i][tx];
			__syncthreads();
			if (newPath < current_d[ty][tx]) {
				current_d[ty][tx] = newPath;
			}
			__syncthreads();
		}

	} else { // pivot-col blocks
		#pragma unroll
		for(int i = 0; i < BLOCK_FACTOR; i++) {
			newPath = current_d[ty][i] + primary_d[i][tx];
			__syncthreads();
                        if (newPath < current_d[ty][tx]) {
                                current_d[ty][tx] = newPath;
                        }
			__syncthreads();
		}
	}

	if (v1 < n && v2 < n) {
		d[cell_current] = current_d[ty][tx];
	}
}

template <int BLOCK_SIZE, int THREAD_SIZE> __global__ void cal_phase_three(const unsigned int Round, const unsigned int n, const size_t pitch, int * const d) 
{
	if(blockIdx.x == Round || blockIdx.y == Round) return;
	int newPath;
	int path;

	const int tx = threadIdx.x*THREAD_SIZE;
	const int ty = threadIdx.y*THREAD_SIZE;

	//the ID of the rest of the blocks (no blocks that have been calculated in phase 1 & 2)
  	const int v1 = blockDim.y*blockIdx.y*THREAD_SIZE + ty;
	const int v2 = blockDim.x*blockIdx.x*THREAD_SIZE + tx;
	
	__shared__ int primaryRow_d[BLOCK_SIZE * THREAD_SIZE][BLOCK_SIZE * THREAD_SIZE];
	__shared__ int primaryCol_d[BLOCK_SIZE * THREAD_SIZE][BLOCK_SIZE * THREAD_SIZE];
	
	// pivot-row and pivot-col from phase 1 & 2
	int v1Row = BLOCK_SIZE*Round*THREAD_SIZE + ty; 
	int v2Col = BLOCK_SIZE*Round*THREAD_SIZE + tx;
	
	int idx, idy;
 	int weightID;
	// each cuda block has its own pivot-block results calculated from phase 1 and phase 2 (like a cross shape)
	#pragma unroll
	for (int i = 0; i < THREAD_SIZE; i++) {
		#pragma unroll
		for (int j = 0; j < THREAD_SIZE; j++) {
			idx = tx + j;
			idy = ty + i;
			if(v1Row+i < n && v2+j < n) {
				weightID = (v1Row+i)*pitch + v2 + j;
				primaryRow_d[idy][idx] = d[weightID];
			} else {
				primaryRow_d[idy][idx] = INF;
			}

			if (v1 + i  < n && v2Col + j < n)
			{
				weightID = (v1 + i) * pitch + v2Col + j;
				primaryCol_d[idy][idx] = d[weightID];
			}
			else
			{
				primaryCol_d[idy][idx] = INF;
			}
		}
	}
	__syncthreads();
	#pragma unroll
        for (int i = 0; i < THREAD_SIZE; i++) {
                #pragma unroll
                for (int j = 0; j < THREAD_SIZE; j++) { 
			if(v1+i < n && v2 +j < n) {
				weightID = (v1+i)*pitch + v2 + j;
				path = d[weightID];
				idx = tx + j;
                        	idy = ty + i;
				
				#pragma unroll
				for(int k=0; k < BLOCK_SIZE*THREAD_SIZE; k++) {
					newPath = primaryCol_d[idy][k] + primaryRow_d[k][idx];
					if(path > newPath) {
						path = newPath;
						
					}
					d[weightID] = path;
				}		
				//d[weightID] = path;
			}
		}
	}
}	


int ceil(int a, int b)
{       return (a + b -1)/b;
}

template <int BLOCK_SIZE, int THREAD_SIZE> void block_APSP()
{
	int *dev_d = 0;

	size_t pitch;
	size_t pitch_int;

	// Size of blocking_factor
	const int BLOCK_FACTOR = BLOCK_SIZE * THREAD_SIZE;
	int round = ceil(n, BLOCK_FACTOR);
	printf("FW-Round = %d\n",round);

	cudaError_t cudaStatus;
	cudaStream_t cpyStream;

	// Initialize the grid and block dimensions here
	dim3 dimGridP1(1, 1, 1);
	dim3 dimGridP2(round, 2, 1);
	dim3 dimGridP3(round, round, 1);
	
	dim3 dimBlockP1(BLOCK_FACTOR, BLOCK_FACTOR, 1);
	dim3 dimBlockP2(BLOCK_FACTOR, BLOCK_FACTOR, 1);
	dim3 dimBlockP3(BLOCK_SIZE, BLOCK_SIZE, 1);
	
	if (gDebug)
        {
                printf("|V| %d\n", n);

                printf("Phase 1\n");
                printf("Dim Grid:\nx - %d\ny - %d\nz - %d\n", dimGridP1.x, dimGridP1.y, dimGridP1.z);
                printf("Dim Block::\nx - %d\ny - %d\nz - %d\n", dimBlockP1.x, dimBlockP1.y, dimBlockP1.z);

                printf("\nPhase 2\n");
                printf("Dim Grid:\nx - %d\ny - %d\nz - %d\n", dimGridP2.x, dimGridP2.y, dimGridP2.z);
                printf("Dim Block::\nx - %d\ny - %d\nz - %d\n", dimBlockP2.x, dimBlockP2.y, dimBlockP2.z);

                printf("Phase 3\n");
                printf("Dim Grid:\nx - %d\ny - %d\nz - %d\n", dimGridP3.x, dimGridP3.y, dimGridP3.z);
                printf("Dim Block::\nx - %d\ny - %d\nz - %d\n", dimBlockP3.x, dimBlockP3.y, dimBlockP3.z);
        }

	// Create new stream to copy data	
	cudaStatus = cudaStreamCreate(&cpyStream);
	HANDLE_ERROR(cudaStatus);
	// Allocate GPU buffers for matrix of shortest paths d(G)
	cudaStatus = cudaMallocPitch(&dev_d, &pitch, n*sizeof(int), n);
	HANDLE_ERROR(cudaStatus);	
		
	pitch_int = pitch / sizeof(int);
	//printf("pitch_int = %d\n", pitch_int);
	
	cudaMemcpy2DAsync(dev_d, pitch, Dist, n*sizeof(int), n*sizeof(int), n, CMCPYHTD, cpyStream);

	for(int r = 0; r < round; r++) { 
		cal_phase_one<BLOCK_FACTOR><<<1, dimBlockP1>>>(r, n, pitch_int, dev_d);
		cal_phase_two<BLOCK_FACTOR><<<dimGridP2, dimBlockP2>>>(r, n, pitch_int, dev_d);
		cal_phase_three<BLOCK_SIZE, THREAD_SIZE><<<dimGridP3, dimBlockP3>>>(r, n, pitch_int, dev_d);	
	}

	// Check for any errors launching the kernel
    	cudaStatus = cudaGetLastError();
	HANDLE_ERROR(cudaStatus);

	// cudaDeviceSynchronize waits for the kernel to finish, and returns
    	// any errors encountered during the launch.
    	cudaStatus = cudaDeviceSynchronize();
	HANDLE_ERROR(cudaStatus);

	cudaStatus = cudaMemcpy2D(Dist, n*sizeof(int), dev_d, pitch, n*sizeof(int), n, CMCPYDTH);
	HANDLE_ERROR(cudaStatus);
	
	cudaStatus = cudaFree(dev_d);
	HANDLE_ERROR(cudaStatus);
	
	return;
}

void print_graph()
{
	for(int v1 = 0; v1 < n; v1++)
	{
		for (int v2 = 0; v2 < n; v2++ )
		{	
			if (Dist[v1 * n + v2] < INF)
				printf("%d ", Dist[v1 * n + v2]);
			else
				printf("INF ");
		}
		printf("\n");
	}
	printf("\n");
}

void input(char *inFileName)
{
        FILE *infile = fopen(inFileName, "r");
        fscanf(infile, "%d %d", &n, &m);

        for (int i = 0; i < n*(n+1); ++i)
                Dist[i]=INF;

        for (int i = 0; i < n; ++i)
                Dist[i*n + i ] = 0;
        while (--m >= 0) {
                int v1, v2, w;
                fscanf(infile, "%d %d %d", &v1, &v2, &w);
                Dist[(v1-1) * n + (v2-1)] = w;
        }
}

void output(char *outFileName)
{       FILE *outfile = fopen(outFileName, "w");
        for (int i = 0; i < n; ++i) {
                for (int j = 0; j < n; ++j) {
                        if (Dist[i*n + j] >= INF)  fprintf(outfile, "INF ");
                        else fprintf(outfile, "%d ", Dist[i*n+j]);
                }
                fprintf(outfile, "\n");
        }
}


int init_device()
{	// Task 1: Device Initialization
	int deviceCount;
	//int i;
	cudaGetDeviceCount(&deviceCount);
	if (deviceCount == 0) {
		printf("\n\tNo CUDA device found.\n\n");
		return 0;
	} else {	
		for( int i = 0; i < deviceCount; i++) {
			cudaDeviceProp prop;
			cudaGetDeviceProperties(&prop, i);
			//printf("sharedMemPerBlock: %d\n", prop.sharedMemPerBlock);
			//printf("maxThreadsPerBlock: %d\n", prop.maxThreadsPerBlock);
		}
		//printf("deviceCount: %d\n",deviceCount);
		cudaSetDevice(0);
    		return 0;
	}
}

int main(int argc, char *argv[]) 
{

	init_device();
	input(argv[1]);
	
	// Initialize CUDA Event
	cudaEvent_t start,stop;
	float elapsedTime;
	cudaEventCreate(&start);
	cudaEventCreate(&stop);
	cudaEventRecord(start,0);

	block_APSP<BLOCK_WIDTH, THREAD_WIDTH>();

	// Finish recording
	cudaEventRecord(stop,0);
	cudaEventSynchronize(stop);
	
	// Calculate elasped time
	cudaEventElapsedTime(&elapsedTime,start,stop);
	elapsedTime /= 1000;
	printf ("Time : %f ms\n", elapsedTime*1000);

	output(argv[2]);

	return 0;
}
