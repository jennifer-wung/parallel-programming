#ifndef QUADTREE_H
#define QUADTREE_H

#include "particle_system.h"
//#include "particle.cpp"

class Quadtree
{
	public:
		/*Create a quadtree based on a particle system.*/
		Quadtree(ParticleSystem* ps, long double thres);
		/**Create a quadtree representing a certain amount of space.
		* @param iL : left hand coordinate
		* @param iR : right hand coordinate
		* @param iB : top coordinate
		* @param iT : bottom coordinate
		* @param iMe : parent node
		* @note : iMe should only be NULL for empty nodes
		*/
		Quadtree(long double iL, long double iR, long double iB, long double iT, Particle* iMe, long double thres);
		
		~Quadtree();
		/*Add a node to this tree*/
		void add(Particle* node);
		/*Add a particle system to this tree*/
		void add(ParticleSystem* ps);
		/*Delete all contents of this tree*/
		void clear();
		/*Updates a particle's forces*/
		void update(Particle* p) const;
		/*Get the parent node this tree represents*/
		Particle* getMe();
		/*recalculate the central mass and positions*/
		void recalculateMe();
		/*to know this node should be pointed to one of the child's tree*/
		unsigned int getQuadrant(Particle* node) const;
		/*return each of the 4 sides*/
		long double getLeft() const;
		long double getRight() const;
		long double getTop() const;
		long double getBottom() const;
		/*Get a pointer to one of my children*/
		Quadtree* getChild(unsigned int indice);
		/*Returns true if I have children*/
		bool isParent() const;
		
	private:
		void makeChildren();//Allocates space for and creates children Quadtrees

		long double left, right;
		long double top, bottom;
		long double theda;
		bool parent;
		Particle* me;
		Quadtree** mChildren;
};

#endif
