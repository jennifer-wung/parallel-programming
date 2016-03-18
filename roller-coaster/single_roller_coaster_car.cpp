#include <iostream>
#include <stdlib.h>

#include "roller_coaster.h"

using namespace std;

int main(int argc, char *argv[])
{
	int Capacity, No_of_rides, Riding_time, No_of_passengers;//number of passengers in the park
	int i;
	PassengerThread *Passenger[MAX_PASSENGERS];

	if (argc != 5) {
		cout << "#_of_passengers Capacity Riding_time #_of_rides" << endl;
	} else {
		No_of_passengers = abs(atoi(argv[1]));
          	Capacity = abs(atoi(argv[2]));
		Riding_time = abs(atoi(argv[3]));
          	No_of_rides = abs(atoi(argv[4]));
	}
	
	if (Capacity > No_of_passengers) { 
          	cout << "Please check your input, car capacity "
               	<< "must be less than the no. of passengers" << endl;
          	exit(1);
     	}

	if (No_of_passengers > MAX_PASSENGERS +1) {
          	cout << "The number of passengers is too large. " 
               	<< "Reset to " << MAX_PASSENGERS << endl;
          	No_of_passengers = MAX_PASSENGERS;
     	}	
     	if (Capacity > MAX_CAPACITY+1) {
          	cout << "Car capacity is too large.  " 
               	<< "Reset to " << MAX_CAPACITY << endl;
          	Capacity = MAX_CAPACITY;
     	}
	
	//cout << "Car thread ready to be generated!"<<endl;	

	CarThread Car(Capacity, No_of_rides, Riding_time); // create and run the car
	Car.start();
	
	for (i = 1; i <= No_of_passengers; i++) { // create and run passengers
		//cout << "passenger" << i<<" thread ready to be generated!"<<endl;
		Passenger[i] = new PassengerThread(i, Capacity);
		Passenger[i]->start();
	}
	
	Car.join(); 
	return 0;
}

