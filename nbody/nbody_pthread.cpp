#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <cstring>
#include <string>
#include <time.h>
#include <unistd.h>
#include <sys/time.h>     /* gettimeofday() */
#include <cmath>
#include <pthread.h>

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
	/* body mass */
	double m;
};

/* calculate info. structure */
struct CalInfo{
	/* an array of NBody old info. */
    	struct Body *NbodyData;
    	/* partition info. */
    	int start,count;
    	/* an array of new velocity info. from bodies in thread per se*/
    	double *vx,*vy;
	/* an array of new force info. from bodies in thread per se*/
    	double *fx,*fy;
	/*duration*/
	double duration;
	/*number of bodies*/
	int numBody;
};

typedef struct {
    	int secs;
    	int usecs;
} TIME_DIFF;
/* function declaration */
/* function for setting the partition for each processor */
void setStartEnd(int size,int rank,int oriStart,int oriEnd,int *istart,int *icount);
/* function for calculating new info of bodies */
void *calculateNewInfo(void *sendInfo);
/* function for update new info of bodies */
void *updateNewInfo(void *sendInfo);
/* calculate time difference*/
TIME_DIFF * my_difftime (struct timeval *, struct timeval *);
int main(int argc, char* argv[])
{
	/*Timing variables*/
	struct timeval myTVstart, myTVend;
	TIME_DIFF * difference;
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
	while(getline(infile, line)) {
		istringstream iss(line);
		if (!(iss >> Nbody[i].x >> Nbody[i].y >> Nbody[i].vx >> Nbody[i].vy)) { 
			cout << "Error occurs during reading!" << endl;
			break; 
		} // error
		newVx[i] = Nbody[i].vx;
		newVy[i] = Nbody[i].vy;
		force_x[i] = 0;
		force_y[i] = 0;
		Nbody[i].m = body_mass;
		i++;
	}
	infile.close();

	/***variables for pthread***/
	if(num_of_threads > num_of_body)
		num_of_threads = num_of_body;
	pthread_t threads[num_of_threads];
    	pthread_attr_t pattr;
    	pthread_attr_init(&pattr);
    	pthread_attr_setdetachstate(&pattr,PTHREAD_CREATE_JOINABLE);
    	struct CalInfo sendInfo[num_of_threads];
	
	for(int t = 0; t < num_of_threads; t++) {
		sendInfo[t].NbodyData = Nbody;
		setStartEnd(num_of_threads, t, 0, num_of_body-1, &sendInfo[t].start,&sendInfo[t].count);
		sendInfo[t].vx=newVx;
        	sendInfo[t].vy=newVy;
		sendInfo[t].fx=force_x;
		sendInfo[t].fy=force_y;
		sendInfo[t].duration = duration;
		sendInfo[t].numBody = num_of_body;
	}

	/*** start to simulate N-body***/
	gettimeofday (&myTVstart, NULL);
	for (int iter = 0; iter < num_of_steps; iter++) {
		for (int t = 0; t < num_of_threads; t++){
			pthread_create(&threads[t],&pattr,calculateNewInfo,&sendInfo[t]);
		}
		for(int t = 0; t < num_of_threads; t++)
			pthread_join(threads[t],NULL);
		/* update the new data */
        	for(int t = 0; t < num_of_threads; t++)
            		pthread_create(&threads[t],NULL,updateNewInfo,&sendInfo[t]);
        	for(int t = 0; t < num_of_threads; t++)
            		pthread_join(threads[t],NULL);

		if (xWin_switch.compare("enable")==0) {
			xWin_show.clear();
			for(int n = 0; n< num_of_body;n++) 
				xWin_show.draw((int)xWin_scale*Nbody[n].x-xWin_xmin,(int)xWin_scale*Nbody[n].y+xWin_ymin);	
			xWin_show.flush();
		}
		//sleep(0.9);
	}
	gettimeofday (&myTVend, NULL);
	difference = my_difftime (&myTVstart, &myTVend);
	cout << "Time elaspsed: "<<difference->secs<<"."<<difference->usecs <<"s."<<endl;
	/* terminate pthread */
    	pthread_attr_destroy(&pattr);
    	pthread_exit(NULL);
	return 0;
}

void setStartEnd(int numThreads,int rank,int oriStart,int oriEnd,int *istart,int *icount){
	int iLength, bpt, iend;
	iLength = oriEnd-oriStart+1;
	bpt = (int)iLength/numThreads; // bodies per thread
	if(iLength > 0) {
		if(rank != numThreads-1) {
			*istart = oriStart + rank*bpt;
			iend = oriStart + (rank+1)*(bpt)-1;
		} else {
			*istart = oriStart + rank*bpt;
			iend = oriEnd;
		}

	}else{
		*istart = 1;
		iend =0;
	}
	*icount = iend-*istart+1;
}

void *calculateNewInfo(void *sendInfo){
	struct CalInfo info = *(struct CalInfo*)(sendInfo);
	struct Body *oldBody = info.NbodyData;
	int start = info.start;
	int count = info.count;
	double *fx = info.fx;
	double *fy = info.fy;
	int num_of_body = info.numBody;

	for (int j = 0; j < count; j++){
                for (int k = 0; k < num_of_body; k++){
                        if(start+j == k)
                        	continue; //there is no need to calculate the effect from itself
                      	double delta_x = oldBody[k].x - oldBody[start+j].x;
                        double delta_y = oldBody[k].y - oldBody[start+j].y;
               		double distance = sqrt(pow(delta_x,2)+pow(delta_y,2));	
                   	if(distance==0)
                    		continue; /* if two bodies have the same position, skip the force calculation */
                  	double accer = G_CONS*(oldBody[start+j].m)/(distance*distance);
                        fx[start+j] += accer*delta_x/distance;
                      	fy[start+j] += accer*delta_y/distance;
		}
	}
	//pthread_exit(NULL);
}

void *updateNewInfo(void *sendInfo){
	struct CalInfo info = *(struct CalInfo*)(sendInfo);
	struct Body *oldBody = info.NbodyData;
	int start = info.start;
	int count = info.count;
	double *fx = info.fx;
	double *fy = info.fy;
    	double *newVx = info.vx;
    	double *newVy = info.vy;
	double duration = info.duration;

    	/* update the position of bodies */
    	for(int i=0;i<count;i++){
		newVx[start+i] = oldBody[start+i].vx + fx[start+i]*duration;
		newVy[start+i] = oldBody[start+i].vy + fy[start+i]*duration;
                oldBody[start+i].x = oldBody[start+i].x + newVx[start+i]*duration;
                oldBody[start+i].y = oldBody[start+i].y + newVy[start+i]*duration;
                oldBody[start+i].vx = newVx[start+i];
                oldBody[start+i].vy = newVy[start+i];
		fx[start+i] = 0;
		fy[start+i] = 0;
    	}
	//pthread_exit(NULL);
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
