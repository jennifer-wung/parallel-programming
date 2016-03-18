#include <iostream>
#include <cstdlib>
#include <cstring>
#include <string>
#include <time.h>
#include <ctime>
#include <cmath>
#include <sys/time.h>     /* gettimeofday() */
#include <omp.h>

#include "particle_system.h"
#include "particle.cpp"
#include "quadtree.h"
#include "barnes_hut.h"
#include "xwindow.h"
using namespace std;

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
	int totalSecs_tree = 0;
	int totalUSecs_tree = 0;
	int totalSecs_compute = 0;
	int totalUSecs_compute = 0;
	/***input information***/
	unsigned int num_of_threads, num_of_steps;
	long double body_mass, duration, theda;
	char* infileName;
	string xWin_switch; // 1-->enable  0-->disable
	int xWin_size;
        double xWin_length, xWin_xmin, xWin_ymin;	
	if (argc > 6) { 
		num_of_threads = abs(atoi(argv[1]));
          	body_mass = abs(atof(argv[2]));
		num_of_steps = abs(atoi(argv[3]));
          	duration = abs(atof(argv[4]));
		infileName = argv[5];
		theda = atof(argv[6]);
		xWin_switch = argv[7];
		if(argc == 12) {
			xWin_xmin = atof(argv[8]);
			xWin_ymin = atof(argv[9]);
			xWin_length = abs(atof(argv[10])); //axis range
			xWin_size = abs(atof(argv[11])); // resolution
			cout << "xmin: " << xWin_xmin << "; ymin: " 
			<<xWin_ymin<<endl;
			cout << "xWin_length: " << xWin_length 
			<< "; xWin_size" << xWin_size << endl;
		}
	}
	
	/***start to do barnes hut algorithm***/
	gettimeofday (&myTVstart, NULL);
	ParticleSystem mPs(infileName, body_mass);
	gettimeofday (&myTVend, NULL);
	difference = my_difftime (&myTVstart, &myTVend);
        cout << "Time elaspsed for IO: "<<difference->secs<<"."<<difference->usecs<<"s."<<endl;

	cout << "particle system built!" << endl;
	if( mPs.getSize() < 1 )
	{
		cout << "No particles in file"<<endl;
		exit(1);
	}
	double xWin_scale;
        XWindowDisplay xWin_show;
        if (xWin_switch.compare("enable") == 0) {
                xWin_scale = xWin_size/xWin_length;
                xWin_xmin = xWin_xmin*xWin_scale;
                xWin_ymin = xWin_ymin*xWin_scale;
                xWin_show.initGraph(xWin_size, xWin_size,xWin_scale);
        }

	gettimeofday (&myTVstart, NULL);
	for (int iter = 0; iter < num_of_steps; iter++) {
		//cout << "iter: "<<iter<<endl<<endl;
		gettimeofday (&myTVstart, NULL);
		Quadtree mQT( &mPs, theda );
		gettimeofday (&myTVend, NULL);
		difference = my_difftime (&myTVstart, &myTVend);
		totalSecs_tree += difference->secs;
		totalUSecs_tree += difference->usecs;
		//cout << "Runnnig Barnes-Hut on all particles" <<endl;
		gettimeofday (&myTVstart, NULL);
		BarnesHut mBH( &mPs, &mQT, num_of_threads, duration);
		mBH.setLast( mPs.getSize() );
		mBH.threadDistri();
		//cout << "BH ends" <<endl;

		Particle *p;
		if (xWin_switch.compare("enable")==0) {
			xWin_show.clear();
			for (int m = 0; m < mPs.getSize(); m++) {
				p = mPs.getParticle(m);
				//cout << "vx="<<p->vx<<", x="<<p->x<<", fx="<<p->fx<<endl;
				p->vx = p->vx + (p->fx)*duration;
				p->vy = p->vy + (p->fy)*duration;
				p->x = p->x + (p->vx)*duration;
				p->y = p->y + (p->vy)*duration;	
				//cout <<"x="<<(int)xWin_scale*(p->x) - xWin_xmin<<", y="<<(int)xWin_scale*(p->y) + xWin_ymin<<endl;
				xWin_show.draw((int)xWin_scale*(p->x) - xWin_xmin,(int)xWin_scale*(p->y) + xWin_ymin);
				xWin_show.flush();
			}
		}
		else {
			//omp_set_num_threads(num_of_threads);
			//#pragma omp parallel
			//{
			//#pragma omp for schedule(static)
			for (int m = 0; m < mPs.getSize(); m++) {
				p = mPs.getParticle(m);
				//cout << "vx="<<p->vx<<", x="<<p->x<<", fx="<<p->fx<<endl;
				p->vx = p->vx + (p->fx)*duration;
				p->vy = p->vy + (p->fy)*duration;
				p->x = p->x + (p->vx)*duration;
				p->y = p->y + (p->vy)*duration;	
			}
			//}
		}
		//sleep(0.5);
		mPs.update();
		gettimeofday (&myTVend, NULL);
		difference = my_difftime (&myTVstart, &myTVend);
                totalSecs_compute += difference->secs;
                totalUSecs_compute += difference->usecs;
		
	}
	//gettimeofday (&myTVend, NULL);
        //difference = my_difftime (&myTVstart, &myTVend);
        cout << "Time elaspsed for building tree: "<<totalSecs_tree<<"."<<totalUSecs_tree<<"s."<<endl;
        cout << "Time elaspsed for computing: "<<totalSecs_compute<<"."<<totalUSecs_compute<<"s."<<endl<<endl<<endl;
	
	
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

