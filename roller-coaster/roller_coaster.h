#ifndef  _ROLLER_H
#define  _ROLLER_H

#include <thread.h>
#include <condvar.h>

#define MAX_PASSENGERS   10   // Maximum number of passengers 
#define MAX_CAPACITY     10   // Maximum capacity of the car  
#define MAX_NO_RIDES     10   // Maximum number of car rides

//------------------------------------------------------------------------
//// PassengerThread class definition 
////------------------------------------------------------------------------

class PassengerThread: public Thread //derived class for Thread class
{

	public: 
		PassengerThread(int Number, int capacity);
		//~PassengerThread();
		void TakeRide();
		void* run();
	private:
		
		int passenger_Id;
		int Capacity; //capacity of the car
		//void ThreadFunc();
};


//------------------------------------------------------------------------
//// carThread class definition
////------------------------------------------------------------------------

class CarThread: public Thread
{
	public:
          	CarThread(int capacity, int noOfRides, int ridingTime);
		//~CarThread();
		void* run();	
     	private:
          	//void ThreadFunc();
          	int  Capacity;      // capcity of the car  
          	int  No_of_rides;   // number of times the car rides
		int  Riding_time;   // duration of each ride
};

#endif
