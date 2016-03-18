#include <iostream>
#include <mpi.h>
#include <sys/time.h>
#include <stdlib.h>
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/Xos.h>
#include <string.h>
using namespace std;

// tags for sends/receives
#define RESULT_TAG 1
#define TERMINATE_TAG 0
#define WORK_TAG 2
#define COUNT_TAG 3
// complex structure holds real/imag part of num
typedef struct complextype
{
	float real, imag;
} Compl;

typedef struct {
    	int secs;
    	int usecs;
} TIME_DIFF;

void master(int numWorkers, int width, int height, int xWin_enable);
void slave(int width, int height, double real_min, double real_max, double imag_min, double imag_max,int rank, int nprocs);

double wallclock(void);
TIME_DIFF * my_difftime (struct timeval *, struct timeval *);

struct timeval myTVstart, myTVend;
TIME_DIFF * difference;

int main(int argc, char *argv[]) {

	int nprocs;
	int rank;
	int numThreads;
	double real_min, real_max, imag_min, imag_max;
	int resol_x, resol_y;
	string xWin_switch;
	int xWin_enable = 0;
	MPI_Status status;
	double totalStartwtime = 0.0, totalEndwtime = 0.0;

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
	if( rank == 0 ) { // Master
		totalStartwtime = wallclock();
		master(nprocs-1, resol_x, resol_y, xWin_enable);
	}

	else { // Slave
		slave(resol_x, resol_y, real_min, real_max, imag_min, imag_max,rank,nprocs);
	}
	if(rank == 0) {
		totalEndwtime = wallclock();
		cout << "Time elaspsed: "<<totalEndwtime-totalStartwtime<<endl;
	}
/*	int pointNum;
	int id;
	if(rank == 0) {
		for(int i = 1; i <= nprocs-1; i++) {
			MPI_Recv(&pointNum,1, MPI_INT, MPI_ANY_SOURCE, COUNT_TAG, MPI_COMM_WORLD, &status);
			id = status.MPI_SOURCE;
			cout << "taskID["<<id<<"] = "<<pointNum<<" rows"<<endl;
		}
	}*/
	// finialize mpi stuff
	MPI_Finalize();
	return 0;
}

void master(int numWorkers, int width, int height, int xWin_enable) 
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
	int this_row, next_row;
	int tasks_not_done;
	int id;
	double totalStartwtime = 0.0, totalEndwtime = 0.0;
	
	/* Set up for dynamic task assignment */

    	next_row = 0;          /* next row to assign to a worker */
    	tasks_not_done = 0;    /* count of workers still working */
	for (int p = 1; p <= numWorkers; p++){
		MPI_Send(&next_row, 1, MPI_INT, p, WORK_TAG, MPI_COMM_WORLD);
		next_row++;
		tasks_not_done++;
	}

	/* Receive results from workers and draw points */
	while(tasks_not_done > 0) {
		/* Receive a row from a worker */
		MPI_Recv(colors,width+1, MPI_INT, MPI_ANY_SOURCE, RESULT_TAG, MPI_COMM_WORLD, &status);
		//cout << "colors[0]: "<<colors[0]<<endl;
		tasks_not_done--;
		id = status.MPI_SOURCE;
		/* More rows? */
        	if (next_row < height) {
			/* If so, give this worker another row to work on */
			//cout << "id: "<<id<<"; next_row: "<<next_row<<endl; 
			MPI_Send(&next_row, 1, MPI_INT, id, WORK_TAG, MPI_COMM_WORLD);
			next_row++;
			tasks_not_done++;
		} else {
			/* Otherwise shut this worker down */
			MPI_Send(&next_row, 0, MPI_INT, id, TERMINATE_TAG, MPI_COMM_WORLD);
		} 
		
		if(xWin_enable) {
			/* Display received data */
			this_row = colors[0];
			for(int col = 0; col < width; col++) {
				XSetForeground (display, gc, colors[col+1]);
				XDrawPoint (display, window, gc, col, this_row);
			}
			XFlush(display);
		}	
	}
	//sleep(5);
}

void slave(int width, int height, double real_min, double real_max, double imag_min, double imag_max,int rank, int nprocs)
{
	int * colors = (int *)malloc(sizeof(int) * (width+1));
	MPI_Status status;
	int row;
	int repeats;
	float temp, lengthsq;
	float scale_real, scale_imag;
	Compl z,c;
	char    processor_name[MPI_MAX_PROCESSOR_NAME];
        int     namelen;
        double startwtime = 0.0, endwtime, delta_t = 0.0;
	int pointCounter; 
	/* Compute factors to scale computational region to window */
    	scale_real = (float) (real_max - real_min) / (float) width;
    	scale_imag = (float) (imag_max - imag_min) / (float) height;
		
	//MPI_Get_processor_name(processor_name,&namelen);
        //fprintf(stdout,"Process %d of %d is on %s\n",rank, nprocs, processor_name);
        //fflush(stdout);
	startwtime = MPI_Wtime();
	MPI_Recv(&row , 1 , MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
	
	/* While master is still sending "work" (not "stop") messages .... */
	while(status.MPI_TAG == WORK_TAG) {
		/* Calculate points and return to master */
		//pointCounter++;
		colors[0] = row;
		c.imag = imag_min + ((float)row*scale_imag);
		for (int col = 0; col < width; col++) {
			z.real = 0.0;
			z.imag = 0.0;
			repeats = 0;
			lengthsq = 0.0;
			c.real = real_min + ((float)col*scale_real);
			while(repeats < 100000 && lengthsq < 4.0) { /* Theorem: If c belongs to M, then |Zn| <= 2. So Zn^2 <= 4 */
				temp = z.real*z.real - z.imag*z.imag + c.real;
				z.imag = 2*z.real*z.imag + c.imag;
				z.real = temp;
				lengthsq = z.real*z.real + z.imag*z.imag; 
				repeats++;		
			}
			colors[col+1] = 1024 * 1024 * (repeats % 256);
		}
		MPI_Send(colors, width+1, MPI_INT, 0, RESULT_TAG, MPI_COMM_WORLD);
		MPI_Recv(&row , 1 , MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
	}
	//endwtime = MPI_Wtime();
	//delta_t = endwtime-startwtime;
	//printf("rank %d -> wall clock time = %f\n", rank, delta_t);
        //fflush(stdout);
	//MPI_Send(&pointCounter, 1, MPI_INT, 0, COUNT_TAG, MPI_COMM_WORLD);
	/* END OF CHANGED SECTION  */
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
}


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
