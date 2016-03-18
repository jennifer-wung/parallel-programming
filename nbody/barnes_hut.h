#ifndef BARNES_HUT_H
#define BARNES_HUT_H

#include "particle_system.h"
#include "quadtree.h"
#include "thread.h"

class BarnesHut : public Thread
{
	public:
      		BarnesHut(ParticleSystem *iPS, Quadtree* iQT, unsigned int threads, long double deltaT);
		/*Run the Barnes-Hut algorithm on this.*/
		void* run();
		/*distribute the particles to threads*/
		void threadDistri();
		/*Get the particle system associated with this*/
		ParticleSystem* getParticleSystem();
		/*Get the quad tree associated wth this.*/
		Quadtree* getQuadtree();
		/*Return the indice of the first particle to act upon*/
		unsigned int getFirst() const;
		/*Return the indice of the last particle to act upon*/
		unsigned int getLast() const;
		/*Set the number of threads to be used*/
		void setNumberOfThreads( unsigned int num );
		/*Set the indice of the first particle to be acted upon.*/
		void setFirst( unsigned int nFirst );
		/*Set the indice of the last particle to be acted upon.*/
		void setLast( unsigned int nLast );
	private:
		ParticleSystem* ps;
		Quadtree* qt;
		unsigned int numThreads;
		unsigned int first;
		unsigned int last;
		long double duration;
};

#endif
