#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <cstring>
#include <string>
#include <time.h>
#include <ctime>
#include <cmath>
#include <omp.h>    /* for OpenMP */
#include <sys/time.h>     /* gettimeofday() */

#include "xwindow.h"
using namespace std;
/* Gravitational constant */
const double G_CONS=0.0000000000667259;
//const double G_CONS=6.67259;

/* body structure */
struct Body{
	/* position */
        double x,y;
       	/* velocity */
       	double vx,vy;
};

typedef struct {
    	int secs;
    	int usecs;
} TIME_DIFF;
TIME_DIFF * my_difftime (struct timeval *, struct timeval *);

int main(int argc, char* argv[])
{
	/*Timing variables*/
	struct timeval myTVstart, myTVend;
	TIME_DIFF * difference;
	double startTime;
	/***input information***/
	int num_of_threads, num_of_steps;
	double body_mass, duration;
	ifstream infile;
	string xWin_switch; // 1-->enable  0-->disable
	int xWin_size;
	double xWin_length, xWin_xmin, xWin_ymin;
	if (argc > 6) { 
		num_of_threads = abs(atoi(argv[1]));
          	body_mass = abs(atof(argv[2]));
		num_of_steps = abs(atoi(argv[3]));
          	duration = abs(atof(argv[4]));
		infile.open(argv[5]);
		if(!infile.is_open()) {
			cout << "File Could not be open" << endl;
			exit(1);
		} 
		xWin_switch = argv[6];
		if(argc == 11) {
			xWin_xmin = atof(argv[7]);
			xWin_ymin = atof(argv[8]);
			xWin_length = abs(atof(argv[9])); //axis range
			xWin_size = abs(atof(argv[10])); // resolution
			//cout << "xmin: " << xWin_xmin << "; ymin: "<<xWin_ymin<<endl;
			//cout << "xWin_length: " << xWin_length << "; xWin_size" << xWin_size << endl;
		}
	}
	
	/***Read input file (initialize N-Body) and set N-Body variables***/
	string line;
	getline(infile,line);
	cout << "Number of bodies: " << line << endl;
	int num_of_body = abs(atoi(line.c_str()));
	int i,j,k;
    	struct Body Nbody[num_of_body];    /* record the properties of all bodies */
    	double newVx[num_of_body];        /* for calculate each body's velocity */
    	double newVy[num_of_body];
	double force_x[num_of_body];
	double force_y[num_of_body];
	double xWin_scale;
	XWindowDisplay xWin_show;
	if (xWin_switch.compare("enable") == 0) {
		xWin_scale = xWin_size/xWin_length;
		xWin_xmin = xWin_xmin*xWin_scale;
		xWin_ymin = xWin_ymin*xWin_scale;
		xWin_show.initGraph(xWin_size, xWin_size,xWin_scale);
	}

	i = 0;
	while(getline(infile, line))
	{
		istringstream iss(line);
		if (!(iss >> Nbody[i].x >> Nbody[i].y >> Nbody[i].vx >> Nbody[i].vy)) { 
			cout << "Error occurs during reading!" << endl;
			break; 
		} // error
		newVx[i] = Nbody[i].vx;
		newVy[i] = Nbody[i].vy;
		force_x[i] = 0;
		force_y[i] = 0;
		i++;
	}
	infile.close();
	
	/*** start to simulate N-body***/
	startTime = omp_get_wtime();
	gettimeofday (&myTVstart, NULL);
	for (int iter = 0; iter < num_of_steps; iter++) {
	/* set the number of threads */
    	omp_set_num_threads(num_of_threads);
	#pragma omp parallel private (j,k)
	{
		#pragma omp for schedule(static) //collapse(2)
		for (j = 0; j < num_of_body; j++){
			for (k = 0; k < num_of_body; k++){
				if(j == k) 
					continue; //there is no need to calculate the effect from itself
				double delta_x = Nbody[k].x-Nbody[j].x;
				double delta_y = Nbody[k].y-Nbody[j].y;
				double distance = sqrt(pow(delta_x,2)+pow(delta_y,2));
				if(distance==0) 
                			continue; /* if two bodies have the same position, skip the force calculation */
				double accer = G_CONS*body_mass/(distance*distance);
				force_x[j] += accer*delta_x/distance;
				force_y[j] += accer*delta_y/distance;
			}
		}	
		// update the new data
		#pragma omp master
		{	
			if (xWin_switch.compare("enable")==0)
				xWin_show.clear();
		}
		#pragma omp for schedule(static)
		for (int m = 0; m < num_of_body; m++) {
			newVx[m] = Nbody[m].vx + force_x[m]*duration;
			newVy[m] = Nbody[m].vy + force_y[m]*duration;
			Nbody[m].x = Nbody[m].x + newVx[m]*duration;
			Nbody[m].y = Nbody[m].y + newVy[m]*duration;
			Nbody[m].vx = newVx[m];
			Nbody[m].vy = newVy[m];
			force_x[m] = 0;
	                force_y[m] = 0;
		}
		//#pragma omp single
		//cout << "force_x[0] = "<< force_x[0] <<"; Nbody[0].x= "<< (int)xWin_scale*Nbody[0].x-xWin_xmin<<endl; 
		//cout << "Nbody[0].x= "<<Nbody[0].x<<"; Nbody[1].y= "<<Nbody[0].y<<endl;
		//		(int)xWin_scale*Nbody[1].x<<"; Nbody[2].x= "<< xWin_scale*Nbody[2].x<< endl;
		
		#pragma omp master
		{
			if (xWin_switch.compare("enable")==0) {
				for(int n = 0; n< num_of_body;n++) {
					xWin_show.draw((int)xWin_scale*Nbody[n].x-xWin_xmin,(int)xWin_scale*Nbody[n].y+xWin_ymin);
				}
				xWin_show.flush(); 
			}
		}

	}
	//sleep(0.1);
	
	}
	gettimeofday (&myTVend, NULL);
	cout<<"(omp_get_wtime)Time elapsed: "<< omp_get_wtime()-startTime <<"s"<<endl;
	difference = my_difftime (&myTVstart, &myTVend);
	cout << "(gettimeofday)Time elaspsed: "<<difference->secs<<"."<<difference->usecs <<"s."<<endl;
	return 0;
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
