#include <iostream>
#include <omp.h>
#include <sys/time.h>
#include <stdlib.h>
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/Xos.h>
#include <string.h>
using namespace std;

// tags for sends/receives
#define RESULT_TAG 1

// complex structure holds real/imag part of num
typedef struct complextype
{
	float real, imag;
} Compl;

double wallclock(void);


int main(int argc, char *argv[]) {

	int nprocs;
	int rank;
	int numThreads;
	double real_min, real_max, imag_min, imag_max;
	int width, height;
	string xWin_switch;

	Display *display;
	Window window;      //initialization for a window
	GC gc;
	int repeats;
        float temp, lengthsq;
        float scale_real, scale_imag;
        Compl z,c;
	int iam = 0, np = 1;
	double totalStartwtime = 0.0, totalEndwtime = 0.0;

	if (argc == 9) {
		numThreads = abs(atoi(argv[1]));
		real_min = atof(argv[2]);
		real_max = atof(argv[3]);
		imag_min = atof(argv[4]);
		imag_max = atof(argv[5]);
		width  = abs(atoi(argv[6]));
		height  = abs(atoi(argv[7]));
		xWin_switch = argv[8];
	} else {
		cout << "arguments aren't enough!" <<endl;
		exit(1);
	}

	if (xWin_switch.compare("enable") == 0) {	
		/* Initialize for graphical display */
		int screen;         //which screen 

		/* open connection with the server */ 
		display = XOpenDisplay(NULL);
		if(display == NULL) {
			cout << "cannot open display\n" <<endl;
			return 0;
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

	/* Compute factors to scale computational region to window */
        scale_real = (float) (real_max - real_min) / (float) width;
        scale_imag = (float) (imag_max - imag_min) / (float) height;

	omp_set_num_threads(numThreads);
	double *delta_t = (double *) malloc(sizeof(double)*numThreads);
	double *startwtime = (double *) malloc(sizeof(double)*numThreads);
	double *endwtime = (double *) malloc(sizeof(double)*numThreads);
	int *point = (int *) malloc(sizeof(int)*numThreads);
	for(int k = 0; k < numThreads; k++) {
		point[k] = 0.0;
		delta_t[k] = 0.0;
	}
	#pragma omp parallel private(iam, np)//shared(imag_min, real_min, scale_imag, scale_real)
	{	
		np = omp_get_num_threads();
        	iam = omp_get_thread_num();
		int *colors = (int *) malloc(sizeof(int)*width);
		#pragma omp master
			totalStartwtime = wallclock();	

		#pragma omp for private(repeats,lengthsq,temp, c,z) schedule(dynamic,1)
		for (int row = 0; row < height; row++) {
			startwtime[iam] = wallclock();
			point[iam]++;
			for(int col = 0; col < width; col++) {
				z.real = 0.0;
	                        z.imag = 0.0;
        	                repeats = 0;
                	        lengthsq = 0.0;
				c.imag = imag_min + ((float)row*scale_imag);			
                        	c.real = real_min + ((float)col*scale_real);
                        	while(repeats < 100000 && lengthsq < 4.0) { /* If c belongs to M, then |Zn| <= 2. So Zn^2 <= 4 */
                                	temp = z.real*z.real - z.imag*z.imag + c.real;
                                	z.imag = 2*z.real*z.imag + c.imag;
                                	z.real = temp;
                                	lengthsq = z.real*z.real + z.imag*z.imag;
                                	repeats++;
                        	}
                        	colors[col] = 1024 * 1024 * (repeats % 256);
			}
			#pragma omp critical
			{	
				if(xWin_switch.compare("enable") == 0){
					for(int col = 0; col < width; col++) {
	                                	XSetForeground (display, gc, colors[col]);
        	                        	XDrawPoint (display, window, gc, col, row);
                	        	}
                        		XFlush(display);
				}	
			}

			endwtime[iam] = wallclock();
			delta_t[iam] = delta_t[iam]+(endwtime[iam]-startwtime[iam]);
		}
		#pragma omp master
			totalEndwtime = wallclock();
    		free(colors);
	}
	//totalEndwtime = wallclock();
	//cout << "Time elaspsed: "<<totalEndwtime-totalStartwtime<<"s."<<endl;
	
	//for(int k = 0; k <numThreads; k++) {
	//	cout << "time: threadID["<<k<<"] = "<<delta_t[k]<<endl;
	//	cout << "point: threadID["<<k<<"] = "<<point[k]<<endl;
	//}
	//sleep(5);	

	return 0;
}

double wallclock(void)
{	struct timeval tv;
	struct timezone tz;
	double t;

	gettimeofday(&tv, &tz);

	t = (double)tv.tv_sec*1000;
	t += ((double)tv.tv_usec)/1000.0;
	t = t/1000;//second
	return t;
}

