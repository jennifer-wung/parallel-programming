#include <iostream>
#include "barnes_hut.h"
#include "thread.h"

using namespace std;
BarnesHut::BarnesHut( ParticleSystem* iPS, Quadtree* iQT, unsigned int threads, long double deltaT) :
	ps(iPS), qt(iQT), first(0), last(0), Thread(), duration(deltaT) 
{
	numThreads = threads;
}

void* BarnesHut::run()
{
	if(first >= last)
		return NULL;
	if(( this->ps == NULL ) || ( this->qt == NULL ))
	{
		cout << "Tried to run Barnes-Hut with a null ps or qt" << endl;
		return NULL;
	}	
	//cout << "number of threads: "<< this->numThreads<<endl;
	if(this->numThreads == 0)
	{
		//cout << "last: "<<this->last<<endl;
		Particle* p;
		for(unsigned int i = this->first; (i < this->last) && (i < this->ps->getSize()); i++) 
		{
			//cout << "start to update forces"<<endl;
			p = this->ps->getParticle(i);
			this->qt->update(p);
		}
		return NULL;
	}
	
//	//cout << "start allocating treads" <<endl;
//	unsigned int tThreads = (this->numThreads > (this->last - this->first)) ? (this->last - this->first) : this->numThreads;
//	unsigned int ppt = (this->last - this->first) / tThreads;
//	//cout << "tThreads: "<< tThreads<<endl;
//	BarnesHut** workers = new BarnesHut*[ this->numThreads ];
//	
//	for(unsigned int j = 0; j < tThreads; j++)
//	{	
//		//cout << j <<"th workers."<<endl;
//		workers[j] = new BarnesHut(this->ps, this->qt, 0);
//		workers[j]->setFirst(this->first + j*ppt);
//		workers[j]->setLast(this->last + (j+1)*ppt);
//		workers[j]->setNumberOfThreads(0);
//		//cout << "worker's number of thread: "<<workers[j]->numThreads<<endl;
//		workers[j]->start(); // thread starts to run
//		//workers[j]->join();
//		//cout << "thread started!"<<endl;
//		//workers[j]->join();
//	}
//	
//	for( unsigned int i = 0; i < tThreads; i++ )
//	{
//		workers[i]->join();
//		delete workers[i];
//	}
	//delete workers;
}

void BarnesHut::threadDistri()
{
	//cout << "start allocating treads" <<endl;
	unsigned int tThreads = (this->numThreads > (this->last - this->first)) ? (this->last - this->first) : this->numThreads;
	unsigned int ppt = (this->last - this->first) / tThreads;
	//cout << "ppt: "<< ppt<<endl;
	BarnesHut** workers = new BarnesHut*[ this->numThreads ];
	
	for(unsigned int j = 0; j < tThreads; j++)
	{	
		//cout << j <<"th workers."<<endl;
		workers[j] = new BarnesHut(this->ps, this->qt, 0, this->duration);
		workers[j]->setFirst(this->first + j*ppt);
		if( j == (tThreads - 1) )
			workers[j]->setLast( this->last );
		else 
			workers[j]->setLast((j+1)*ppt);
		//cout << "workers[j]->lastFirst" << workers[j]->getLast()<<endl;
			
		workers[j]->setNumberOfThreads(0);
		//cout << "worker's number of thread: "<<workers[j]->numThreads<<endl;
		workers[j]->start(); // thread starts to run
		//workers[j]->join();
		//cout << "thread started!"<<endl;
		//workers[j]->join();
	}
	
	for( unsigned int i = 0; i < tThreads; i++ )
        {
        	workers[i]->join();
       	        delete workers[i];
        }
 
}
ParticleSystem* BarnesHut::getParticleSystem()
{
	return this->ps;
}

Quadtree* BarnesHut::getQuadtree()
{ 
	return this->qt;
} 

unsigned int BarnesHut::getFirst() const
{
	return this->first;
}

unsigned int BarnesHut::getLast() const
{
	return this->last;
} 

void BarnesHut::setNumberOfThreads( unsigned int num )
{
	this->numThreads = num;
} 

void BarnesHut::setFirst( unsigned int nFirst )
{
	this->first = nFirst;
}

void BarnesHut::setLast( unsigned int nLast )
{
	this->last = nLast;
}
