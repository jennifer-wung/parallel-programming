#include <iostream>
#include <cstdlib>
#include <time.h>
#include <ctime>
#include <sys/time.h>     /* gettimeofday() */

#include "roller_coaster.h"
#include "condvar.h"
#include "thread.h"

using namespace std;

// static data variables
static int On_board_passenger_ID[MAX_CAPACITY];  // IDs of on board passengers 
static int On_board_passenger_counter = 0;        // number of passengers on board
static int Release_passenger_counter = 0;
static int CheckIn_counter = 1; // 1 --> allow to check in; 0 -->  not allow to check in
static int Unloading_counter = 0; // 1 --> allow to unload passenger; 0 --> not allow to unload
// static Mutex 
static Mutex mutex_queue;
static Mutex mutex_checkIn;
static Mutex mutex_boarding;
static Mutex mutex_riding;
static Mutex mutex_unloading;

// static CondVar
static CondVar cond_queue(mutex_queue);
static CondVar cond_checkIn(mutex_checkIn);
static CondVar cond_boarding(mutex_boarding);
static CondVar cond_riding(mutex_riding);
static CondVar cond_unloading(mutex_unloading);

typedef struct {
    	int secs;
    	int usecs;
} TIME_DIFF;

TIME_DIFF * my_difftime (struct timeval *, struct timeval *);
/*Timing variables*/
static struct timeval myTVstart, myTVend;
static TIME_DIFF * difference;
static int totalSecs = 0;
static int totalUSecs = 0;
PassengerThread::PassengerThread(int Number, int capacity) : Thread()
{
	passenger_Id = Number;
	Capacity = capacity;
}

//passenger takeks ride in the car
void PassengerThread::TakeRide()
{	
	/***Queue Mutex***/ 
	//passengers come into queue
	mutex_queue.lock();
	while ( Release_passenger_counter == 0) {
		gettimeofday (&myTVstart, NULL);
		cond_queue.wait();
		gettimeofday (&myTVend, NULL);
		difference = my_difftime (&myTVstart, &myTVend);
		totalSecs += difference->secs;
		totalUSecs += difference->usecs;
		//cout << "totalSecs: "<<totalSecs<<"; totalUSecs: "<<totalUSecs<<endl;
		difference = NULL;
	}

	Release_passenger_counter--;
	mutex_queue.unlock();
	
	/***CheckIn Mutex***/ 
 	//time to check in
	mutex_checkIn.lock();
	if (CheckIn_counter == 0) {
		cout << "waiting for checking in!"<<endl;
		cond_checkIn.wait();
	}
	CheckIn_counter--;
	// save my name and increase the counter
	On_board_passenger_ID[On_board_passenger_counter] = passenger_Id;
	On_board_passenger_counter++;
	cout << "Passenger " << On_board_passenger_ID[On_board_passenger_counter-1]<<" on board"<<endl;
	mutex_checkIn.unlock();
	
	/***Boarding Mutex***/
	// if I am the last one to be on board, boarding completes and the car is full
	if (On_board_passenger_counter == Capacity) {
		mutex_boarding.lock();
                cond_boarding.signal();
                mutex_boarding.unlock();
	}

	/***CheckIn Mutex:***/
	// permission for the next passenger to check in 
	mutex_checkIn.lock();
	CheckIn_counter++;
	cond_checkIn.signal();
	mutex_checkIn.unlock();
	
	/***Riding Mutex***/
	// wait till car done running
	mutex_riding.lock();
	cond_riding.wait();
	//cout << "passenger" <<passenger_Id << " is noticed!" << endl;
	mutex_riding.unlock();
	
	/***Unloading Mutex***/
	mutex_unloading.lock();
	Unloading_counter++;
	cond_unloading.signal();
	cout << "passenger " <<passenger_Id <<" got off the car!" << endl;
	mutex_unloading.unlock();
	
}


void* PassengerThread::run()
{
	srand (time(0));
	int wanderTime;
	while(1) {
		TakeRide();
		wanderTime = rand()%10+1;
		cout << "Passenger " << passenger_Id << " wanders for "<< wanderTime<<" secs"<<endl;  
		sleep(wanderTime);
		cout << "[ Passenger " << passenger_Id << " returns for a ride. ]"<<endl;  
	}
	

}


CarThread::CarThread(int capacity, int noOfRides, int ridingTime) : Thread()
{
	Capacity = capacity;
	No_of_rides = noOfRides;
	Riding_time = ridingTime;
}


void* CarThread::run()
{
	for(int i = 1; i <= No_of_rides; i++) {
		for(int j = 1; j <= Capacity; j++){
			/***Queue Mutex***/
			mutex_queue.lock();
			cond_queue.signal();
			Release_passenger_counter++;
			//cout << "Release passenger counter"<<Release_passenger_counter<<endl;
			mutex_queue.unlock();
		}

		/***Boarding Mutex***/
		mutex_boarding.lock();
		cond_boarding.wait();
		mutex_boarding.unlock();
				
		// print on boarding passenger names
		cout << "The car is running for the " << i << " time with passengers:\n ";
		for (int k = 0; k < Capacity; k++) 
               		cout << " " << On_board_passenger_ID[k];
		cout << endl;
		time_t now = time(0);
		cout << "Car departs at: "<< ctime(&now);
		// car runs for a few seconds
		sleep(Riding_time);
		now = time(0);	
		cout << "Car arrives at: "<< ctime(&now);
	
		cout <<"\nFinish riding!" << endl;
		//On_board_passenger_counter = 0;  // start unloading passengers
		for (int jj = 1; jj <= Capacity; jj++) {
			/***Riding Mutex***/
			// notice passenger that the car done running
			cout <<"****Start unloading "<< jj << "th passenger.****"<<endl;
			mutex_riding.lock();
			cond_riding.signal();
			On_board_passenger_counter--;
			mutex_riding.unlock();
			
			/***Unloading Mutex***/
			// unload passenger one after another
			mutex_unloading.lock();
			if (Unloading_counter == 0)
				cond_unloading.wait();
			Unloading_counter--; 
			//cout << jj <<" th passenger unloaded successfully!"<<endl;
			mutex_unloading.unlock();
			
		}
		// definitely, the car is empty and goes for the next run
		cout << "\n\n" <<endl;
	}
	cout << "The car is shot-off for maintenance" << endl;
	cout << "Average time a passenger wait for a ride(C="<<Capacity<<", T="<<No_of_rides<<") is "<< 
		((double)totalSecs/No_of_rides)<<" secs, "<<((double)totalUSecs/No_of_rides)<<" usecs."<<endl;
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
