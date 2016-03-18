#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cstring>
#include <string>
#include <cmath>
#include <cstdlib>

#include "particle_system.h"

using namespace std;

ParticleSystem::ParticleSystem(char* fileName, long double mass):
mSize(0), mParticles(NULL), minXP(NULL), maxXP(NULL), minYP(NULL), maxYP(NULL)
{
	this->load(fileName,mass);
}


ParticleSystem::~ParticleSystem()
{
	this->clear();
}

void ParticleSystem::load(char* fileName, long double mass)
{
	ifstream infile(fileName);
	if(!infile.is_open()) {
		cout << "File Could not be open" << endl;
		return;
	} 

	string line;
	getline(infile,line);
	cout << "Number of particles: " << line << endl;
	this->mSize = abs(atoi(line.c_str()));

	this->mParticles = new Particle*[ this->mSize ];

	for(unsigned int i = 0; i < this->mSize; i++)
	{	
		getline(infile, line);
		istringstream iss(line);
		
		this->mParticles[i] = new Particle();

		if(!( iss >> this->mParticles[i]->x >> this->mParticles[i]->y >> this->mParticles[i]->vx >> this->mParticles[i]->vy))
		{
			cout << "Error occurs during reading!" << endl;
			return;
		}

		this->mParticles[i]->m = mass;
		//if (i==0)
		//	cout << "x: "<<this->mParticles[i]->x<<"; m = "<<this->mParticles[i]->m <<endl;
		if((this->minXP == NULL) || (this->mParticles[i]->x < this->minXP->x))
			this->minXP = this->mParticles[i];
		if((this->maxXP == NULL) || ( this->mParticles[i]->x > this->maxXP->x ))
			this->maxXP = this->mParticles[i];
		if((this->minYP == NULL) || ( this->mParticles[i]->y < this->minYP->y ))
			this->minYP = this->mParticles[i];
		if((this->maxYP == NULL) || ( this->mParticles[i]->y > this->maxYP->y ))
			this->maxYP = this->mParticles[i];
	}
	infile.close();
	
}



void ParticleSystem::clear()
{
	if(this->mSize < 1)
		return;
	
	delete[] this->mParticles;
	this->mParticles = NULL;
	this->minXP = NULL;
	this->maxXP = NULL;
	this->minYP = NULL;
	this->maxYP = NULL;
	this->mSize = 0;

}


unsigned int ParticleSystem::getSize() const
{
	return this->mSize;
}

Particle* ParticleSystem::getParticle(unsigned int indice) const
{
	if (indice >= this->mSize )
		return NULL;

	return this->mParticles[indice];
}


void ParticleSystem::update()
{
	this->minXP, this->maxXP, this->minYP, this->maxXP = NULL;
	for( unsigned int i = 0; i < this->mSize; i++)
        {
                this->mParticles[i]->fx = 0;
                this->mParticles[i]->fy = 0;
		//if (i==0)
                //        cout << "x: "<<this->mParticles[i]->x<<"; m = "<<this->mParticles[i]->m <<endl;
                if((this->minXP == NULL) || (this->mParticles[i]->x < this->minXP->x))
                        this->minXP = this->mParticles[i];
                if((this->maxXP == NULL) || ( this->mParticles[i]->x > this->maxXP->x ))
                        this->maxXP = this->mParticles[i];
                if((this->minYP == NULL) || ( this->mParticles[i]->y < this->minYP->y ))
                        this->minYP = this->mParticles[i];
                if((this->maxYP == NULL) || ( this->mParticles[i]->y > this->maxYP->y ))
                        this->maxYP = this->mParticles[i];
        }


} 

long double ParticleSystem::getLeft() const
{ 
	return this->minXP->x;
} 

long double ParticleSystem::getRight() const
{ 
	return this->maxXP->x;
} 

long double ParticleSystem::getBottom() const
{
	return this->minYP->y;
}

long double ParticleSystem::getTop() const
{
	return this->maxYP->y;
} 
