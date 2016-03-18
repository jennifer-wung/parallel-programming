#ifndef PARTICLE_CPP
#define PARTICLE_CPP

/**
 ** Class used to represent a point with an x, y and mass.
 **/

class Particle
{
	public:
		long double x, y;
		long double m;
		long double vx, vy;
		long double fx, fy;

		Particle() : x(0),y(0),m(0),vx(0),vy(0),fx(0),fy(0) {}
		Particle(long double iX, long double iY, long double iM) : 
		x(iM), y(iY), m(iM), vx(0),vy(0),fx(0), fy(0) {}	
};

#endif
