#ifndef PARTICLE_SYSTEM_H
#define PARTICLE_SYSTEM_H

#include "particle.cpp"
//using namespace std;

class ParticleSystem
{
	public:
		/*Construct a ParticleSystem and load particles from a file.*/
		ParticleSystem(char* fileName, long double mass);
		~ParticleSystem();

		/*Load a file into this system*/
		void load(char* fileName, long double mass);
		/*Clear all memory allocation*/
		void clear();

		/*Returns the size in particles of this system.*/
		unsigned int getSize() const;

		/*Returns a long double[7] representing a particle.*/
		Particle* getParticle( unsigned int indice) const;
		void update();
		/*Get the leftmost particle's x*/
		long double getLeft() const;
		/*Get the rightmost partilce's x*/
		long double getRight() const;
		/*Get the topmost particle's y*/
		long double getTop() const;
		/*Get the bottom-most particle's y*/
		long double getBottom() const;

	private:
		unsigned int mSize; // # of particles in the system
		Particle** mParticles;// A Nx7 matrix representing all particles
		
		Particle* minXP;
		Particle* maxXP;
		Particle* minYP;
		Particle* maxYP;

};

#endif
