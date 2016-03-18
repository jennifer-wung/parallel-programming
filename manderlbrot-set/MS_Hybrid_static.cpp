#include <iostream>
#include <mpi.h>
#include <omp.h>
#include <sys/time.h>
#include <stdlib.h>
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/Xos.h>
#include <string.h>
using namespace std;

// tags for sends/receives
#define DATA_TAG 0
#define RESULT_TAG 1

// complex structure holds real/imag part of num
typedef struct complextype
{
	float real, imag;
} Compl;

typedef struct {
    	int secs;
    	int usecs;
} TIME_DIFF;

void master(int numWorkers, int row_work, int width, int height, int xWin_enable);
void slave(int row_work, int numThreads, int width, int height, double real_min, double real_max, double imag_min, double imag_max,int rank, int nprocs);

double wallclock(void);
TIME_DIFF * my_difftime (struct timeval *, struct timeval *);

struct timeval myTVstart, myTVend;
TIME_DIFF * difference;

int main(int argc, char *argv[]) {

	int nprocs;
	int rank;
	int numThreads;
	int row_work;
	double real_min, real_max, imag_min, imag_max;
	int resol_x, resol_y;
	string xWin_switch;
	int xWin_enable = 0;
	
	/* Initialize for MPI */
    	if (MPI_Init(&argc, &argv) != MPI_SUCCESS) {
        	fprintf(stderr, "MPI initialization error\n");
        	exit(EXIT_FAILURE);
    	}
    	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	
	if (argc == 9) {
		numThreads = abs(atoi(argv[1]));
		real_min = atof(argv[2]);
		real_max = atof(argv[3]);
		imag_min = atof(argv[4]);
		imag_max = atof(argv[5]);
		resol_x  = abs(atoi(argv[6]));
		resol_y  = abs(atoi(argv[7]));
		xWin_switch = argv[8];
	} else {
		cout << "arguments aren't enough!" <<endl;
		exit(1);
	}
	
	if (xWin_switch.compare("enable") == 0) {	
		xWin_enable = 1;
	}	
	
	row_work = resol_y/(nprocs-1);	

	if( rank == 0 ) { // Master
		master(nprocs-1, row_work, resol_x, resol_y, xWin_enable);
	}

	else { // Slave
		if(rank == nprocs -1)
                        row_work = row_work + (resol_y%(nprocs-1));
		slave(row_work, numThreads, resol_x, resol_y, real_min, real_max, imag_min, imag_max,rank, nprocs);
	}

	// finialize mpi stuff
	MPI_Finalize();
	return 0;
}

void master(int numWorkers, int row_work, int width, int height, int xWin_enable) 
{
	Display *display;
	Window window;      //initialization for a window
	GC gc;
	if (xWin_enable) {
		/* Initialize for graphical display */
		int screen;         //which screen 

		/* open connection with the server */ 
		display = XOpenDisplay(NULL);
		if(display == NULL) {
			fprintf(stderr, "cannot open display\n");
			return;
		}

		screen = DefaultScreen(display);

		/* set window position */
		int x = 0;
		int y = 0;
		/* border width in pixels */
		int border_width = 0;
		/* create window */
		window = XCreateSimpleWindow(display, RootWindow(display, screen), x, y, width, height, border_width,
					BlackPixel(display, screen), WhitePixel(display, screen));
		/* create graph */
		XGCValues values;
		long valuemask = 0;

		gc = XCreateGC(display, window, valuemask, &values);
		XSetBackground (display, gc, WhitePixel (display, screen));
		XSetForeground (display, gc, BlackPixel (display, screen));
		XSetBackground(display, gc, 0X0000FF00);
		XSetLineAttributes (display, gc, 1, LineSolid, CapRound, JoinRound);

		/* map(show) the window */
		XMapWindow(display, window);
		XSync(display, 0);
	}

	/* start calculating Mandelbrot set */
	int *colors = (int *) malloc(sizeof(int) * (width+1));
	MPI_Status status;
	int start_row = 0;
	int this_row;
	int id;
	double MPIStartwtime = 0.0, MPIEndwtime;
	MPIStartwtime = wallclock();
	for( int i = 1; i <= numWorkers; i++ ) {
                MPI_Send( &start_row, 1, MPI_INT, i, DATA_TAG, MPI_COMM_WORLD );
                start_row += row_work;
        }

	for (int row = 0; row < height; row++){
		/* Receive a row from a worker */
		MPI_Recv(colors,width+1, MPI_INT, MPI_ANY_SOURCE, RESULT_TAG, MPI_COMM_WORLD, &status);
		
		if(xWin_enable) {
			/* Display received data */
			this_row = colors[0];
			//cout << "this_row: "<<this_row<<endl;
			for(int col = 0; col < width; col++) {
				XSetForeground (display, gc, colors[col+1]);
				XDrawPoint (display, window, gc, col, this_row);
			}
			XFlush(display);
		}	
	}
	MPIEndwtime = wallclock();
        cout << "rank 0->MPI Time elaspsed: "<<MPIEndwtime-MPIStartwtime<<endl;
	//sleep(5);
}

void slave(int row_work, int numThreads, int width, int height, double real_min, double real_max, double imag_min, double imag_max,int rank, int nprocs)
{
	int * colors = (int *)malloc(sizeof(int) * (width+1));
	MPI_Status status;
	int myrow;
	int repeats;
	float temp, lengthsq;
	float scale_real, scale_imag;
	Compl z,c;
	int iam = 0, np = 1;
        char    processor_name[MPI_MAX_PROCESSOR_NAME];
        int     namelen;
        double startwtime = 0.0, endwtime, delta_t = 0.0;
        double MPICommuStartTime = 0.0, MPICommuEndTime, MPI_delta_t = 0.0;
	int chunkSz = width/numThreads;
	/* Compute factors to scale computational region to window */
    	scale_real = (float) (real_max - real_min) / (float) width;
    	scale_imag = (float) (imag_max - imag_min) / (float) height;

	MPICommuStartTime = wallclock();
	MPI_Recv( &myrow, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
	MPICommuEndTime = wallclock();
        MPI_delta_t = MPI_delta_t + (MPICommuEndTime - MPICommuStartTime);

	omp_set_num_threads(numThreads);
	for (int row = myrow; row < myrow+row_work; row++) {
		/* Calculate points and return to master */
		colors[0] = row;
        	startwtime = wallclock();
		#pragma omp parallel for private (repeats,lengthsq,temp, c,z) schedule(static,chunkSz)
		for (int col = 0; col < width; col++) {
			z.real = 0.0;
			z.imag = 0.0;
			repeats = 0;
			lengthsq = 0.0;
			c.imag = imag_min + ((float)row*scale_imag);
			c.real = real_min + ((float)col*scale_real);
			while(repeats < 100000 && lengthsq < 4.0) { 
				temp = z.real*z.real - z.imag*z.imag + c.real;
				z.imag = 2*z.real*z.imag + c.imag;
				z.real = temp;
				lengthsq = z.real*z.real + z.imag*z.imag; 
				repeats++;		
			}
			colors[col+1] = 1024 * 1024 * (repeats % 256);
		}
		endwtime = wallclock();
        	delta_t = delta_t+(endwtime-startwtime);
		MPICommuStartTime = wallclock();
		MPI_Send(colors, width+1, MPI_INT, 0, RESULT_TAG, MPI_COMM_WORLD);
		MPICommuEndTime = wallclock();
                MPI_delta_t = MPI_delta_t + (MPICommuEndTime - MPICommuStartTime);
	}
	/* END OF CHANGED SECTION  */
	//printf("rank %d -> OpenMP wallclock time = %f; MPI communicate time = %f\n", rank, delta_t, MPI_delta_t);
        //fflush(stdout);
    	free(colors);
}

double wallclock(void)
{       struct timeval tv;
        struct timezone tz;
        double t;

        gettimeofday(&tv, &tz);

        t = (double)tv.tv_sec*1000;
        t += ((double)tv.tv_usec)/1000.0;
        t = t/1000;//second
        return t;
}// millisecond

TIME_DIFF *my_difftime (struct timeval * start, struct timeval * end)
{
    	TIME_DIFF * diff = (TIME_DIFF *) malloc ( sizeof (TIME_DIFF) );
 
    	if (start->tv_sec == end->tv_sec) {
        	diff->secs = 0;
        	diff->usecs = end->tv_usec - start->tv_usec;
    	}
    	else {
        	diff->usecs = 1000000 - start->tv_usec;
        	diff->secs = end->tv_sec - (start->tv_sec + 1);
        	diff->usecs += end->tv_usec;
        	if (diff->usecs >= 1000000) {
            		diff->usecs -= 1000000;
            		diff->secs += 1;
        	}
    	}
     
    	return diff;
}
