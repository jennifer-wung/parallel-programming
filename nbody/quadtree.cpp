#include <iostream>
#include <algorithm>
#include <cmath>
#include <limits>

#include "quadtree.h"
using namespace std;

static const unsigned int NOT_A_QUADRANT = 5;
static const long double QUAD_LEEWAY = 8.0 * numeric_limits<long double>::epsilon();
static const long double G_CONS=0.0000000000667259;

Quadtree::Quadtree(ParticleSystem* ps,long double thres) :
left(0), right(0), top(0), bottom(0), parent(false), me(NULL), mChildren(NULL)
{
	theda = thres;
	this->left = ps->getLeft() - QUAD_LEEWAY;
	this->right = ps->getRight() + QUAD_LEEWAY;
	this->top = ps->getTop() + QUAD_LEEWAY;
	this->bottom = ps->getBottom() - QUAD_LEEWAY;

	long double w = this->right-this->left;
	long double h = this->top - this->bottom;
	
	if(w > h)
	{
		this->bottom -= (w-h)/2.0;
		this->top += (w-h)/2.0;
	}
	else if(h>w)
	{
		this->left -= (h-w)/2.0;
		this->right += (h-w)/2.0;
	}
	//cout << "original box: [" << this->left << ", " <<this->right<<", "<<this->top<<", "<<this->bottom<<"]"<<endl;
	this->add(ps);
}

Quadtree::Quadtree(long double iL, long double iR, long double iB, long double iT, Particle* iMe, long double thres) : 
			left(iL) , right(iR), top(iT), bottom(iB), theda(thres), parent(false), me(iMe), mChildren(NULL)
{
	if( this->left > this->right )
		swap( this->left, this->right );
	if( this->bottom > this->top )
		swap( this->bottom, this->top );
	//cout << "child is made!" <<endl;
}

Quadtree::~Quadtree()
{
	this->clear();
}

void Quadtree::add(Particle* node)
{
	if (node == NULL)
		return;
	if(this->getQuadrant(node) == NOT_A_QUADRANT)
	{	
		cout << "node->x: "<<node->x<<endl;
		cout << "node does not fit in here" << endl;
		return;
	}
	
	if(this->me == NULL)
	{
		//cout << "I'm the only particle!" << endl;
		this->me = node;
		//cout << "this->me->m: "<<this->me->m<<endl;//", this->me->y: "<<this->me->y<<endl;
		return;
	}

	if(!this->parent)
	{	
		this->makeChildren();
		//add current node to children
		//cout << "add->node." << endl;
		//cout << "node quadrant number: "<<this->getQuadrant(node)<<endl;
		this->mChildren[this->getQuadrant(node)]->add(node);
		// and the previos me too
		//cout << "add->me."<<endl;
		//cout << "this->me->x: "<<this->me->x<<", this->me->y: "<<this->me->y<<endl;
		//cout << "box: [" << this->left << ", " <<this->right<<", "<<this->top<<", "<<this->bottom<<"]"<<endl;
		//cout << "quadrant number: "<<this->getQuadrant(this->me)<<endl;
		this->mChildren[this->getQuadrant(this->me)]->add(this->me);
		// now the tree is built
		//this->parent = true;
		// central of mass is changing after adding a particle
		//cout <<"recalculate me."<<endl; 
		this->recalculateMe();
		return;
	}
	//cout << "adding node normally" << endl;
	//cout << "this->me->x: "<<this->me->x<<", this->me->y: "<<this->me->y<<endl;
	//cout << "box: [" << this->left << ", " <<this->right<<", "<<this->top<<", "<<this->bottom<<"]"<<endl;
	//cout << "quadrant number: "<<this->getQuadrant(node)<<endl;
	this->mChildren[this->getQuadrant(node)]->add(node);
	//cout <<"normal calculate me!"<<endl;
	this->recalculateMe();
	return;
}

void Quadtree::add(ParticleSystem* ps)
{
	if(ps==NULL)
		return;
	//cout << "start to build quadtree" <<endl;
	for (unsigned int i=0; i < ps->getSize(); i++)
	{
		//cout << "** " <<i<< "th particle."<<endl;
		//cout << "ps->getParticle->m: "<<ps->getParticle(i)->m<<endl;//", y: "<<ps->getParticle(i)->y<<", m: "<< ps->getParticle(i)->m<<endl;
		this->add( ps->getParticle(i) );
	}
}

void Quadtree::clear()
{
	this->me = NULL;
	if(!this->parent)
		return;
	
	for(unsigned int i=0; i<4; i++)
	{
		delete this->mChildren[i];
		this->mChildren[i] = NULL;
	}	
	this->parent = false;
}

/*recursively check to see if d/r is less than theda, from the topmost layer to the 2nd lower layer*/
void Quadtree::update(Particle* p) const
{
	if(( p == NULL) || (this->me == NULL) || (this->me == p))
		return; // return if empty or only has one particle left
	long double rx = this->me->x - p->x;
	long double ry = this->me->y - p->y;
	long double r2 = rx*rx + ry*ry;
	long double r = sqrt(rx*rx + ry*ry);
	long double r3 = r*r2;
	/*if this->me is not a quadtree, then directly calculate the force exerted on p*/
	if(!this->parent)
	{
		if( fabs( this->me->m ) < 2.0*numeric_limits<long double>::epsilon() )
			return;
		//cout <<"this->me->m: "<<this->me->m<<"; p->m"<<p->m<<endl;
		long double gm = p->m*this->me->m;
		p->fx += G_CONS*gm*rx/r3;
		p->fy += G_CONS*gm*ry/r3;
		//cout << "lowest layer!"<<endl;
		return;
	}

	long double d = this->right - this->left; // size of the cell
	/*if mass is too small or the this->me can't be seen as a node then we jump to the next layer of tree*/
	//cout << "r = "<<r<<"; d = "<<d<<endl;
	//cout << "r/d = "<<r/d<<endl;
	//if(( fabs( this->me->m ) < 2.0*numeric_limits<long double>::epsilon() ) || ((r/d) > this->theda) || (this->getQuadrant(p)!= NOT_A_QUADRANT))
	if((r/d) > this->theda || (this->getQuadrant(p)!= NOT_A_QUADRANT))
	{
		this->mChildren[0]->update(p);
		this->mChildren[1]->update(p);
		this->mChildren[2]->update(p);
		this->mChildren[3]->update(p);
		return;
	}
	else /*this layer of node is far enough to be seen as a node*/
	{
		long double gm = p->m*this->me->m;
		p->fx += G_CONS*gm*rx/r3;
		p->fy += G_CONS*gm*ry/r3;
		//cout << "far enough! p->fx: "<< p->fx <<endl;
		return;
	}	

}

Particle* Quadtree::getMe()
{
	return this->me;
}

void Quadtree::recalculateMe()
{
	//cout << "I'm in the recalculateMe()."<<endl;
	if(!this->parent)
		return;

	Particle* newMe = new Particle();
	//cout << "old me->x = "<< this->me->x<<", old me->y = "<< this->me->y<<", old me->m = "<<this->me->m<<endl; 
	//cout << "this->mChildren->m = "<<this->mChildren[0]->getMe()->m<<", x: "<<this->mChildren[0]->getMe()->x<<endl;
	//this->me->x = 0; this->me->y = 0; 
	//this->me->vx = 0; this->me->vy = 0; this->me->m = 0;
	//cout << "this->mChildren->m = "<<this->mChildren[0]->getMe()->m<<", x: "<<this->mChildren[0]->getMe()->x<<endl;
	for(unsigned int i = 0; i < 4; i++)
	{
		if(this->mChildren[i]->getMe()!=NULL)
		{
			//cout << "add mass, i = " << i << "; this->mChildren->m = "<<this->mChildren[i]->getMe()->m<<", x: "<<this->mChildren[i]->getMe()->x<<endl;
			newMe->m += this->mChildren[i]->getMe()->m;
		}
	}
	//this->me->m = newMe->m;

	//cout << "after adding mass: this->me->m = " << this->me->m <<endl;
	if( fabs( newMe->m ) < 2.0*numeric_limits<long double>::epsilon() )
		return;

	Particle* tChild = NULL;
	for (unsigned int j = 0; j < 4; j++)
	{
		tChild = this->mChildren[j]->getMe();
		if(tChild != NULL) 
		{
			newMe->x += tChild->x * tChild->m;
			newMe->y += tChild->y * tChild->m; 
		}
	}

	newMe->x /= newMe->m;
	newMe->y /= newMe->m;
	this->me = newMe;
	//cout << "new me->x = "<< this->me->x <<", new me->y = "<< this->me->y<<", new me->m"<< this->me->m<<endl;
	return;
}

unsigned int Quadtree::getQuadrant(Particle* node) const
{
	if(node == NULL)
		return NOT_A_QUADRANT;

	//if(( node->x < this->left ) || ( node->x >= this->right ) || ( node->y < this->bottom ) || ( node->y >= this->top ))
	//	return NOT_A_QUADRANT;	

	// 1 0
	// 2 3
	//cout << "get one of the children!" <<endl;
	if(node->x < (this->left+this->right)/2.0)
	{
		if(node->y < (this->top+this->bottom)/2.0)
			return 2;
		return 1; 
	}
	else 
	{
		if(node->y < (this->top+this->bottom)/2.0)
			return 3;
		return 0;
	}
	cout << "does not fit any!"<<endl;
	return NOT_A_QUADRANT;
}

void Quadtree::makeChildren()
{
	if(this->me == NULL)
		return;
		
	if(this->parent) // already made children
	{
		for(unsigned int i = 0; i < 4; i++)
		{
			delete this->mChildren[i];
			this->mChildren[i] = NULL;
		}
	}	
	
	this->mChildren = new Quadtree*[4];
	long double midX = (this->left + this->right)/2.0;
	long double midY = (this->bottom + this->top)/2.0;
	// 1 0
	// 2 3
	this->mChildren[0] = new Quadtree(midX, this->right, midY, this->top, NULL, this->theda);
	this->mChildren[1] = new Quadtree(this->left, midX, midY, this->top, NULL, this->theda);
	this->mChildren[2] = new Quadtree(this->left, midX, this->bottom, midY, NULL, this->theda);
	this->mChildren[3] = new Quadtree(midX, this->right, this->bottom, midY, NULL, this->theda);
	this ->parent = true;

}

long double Quadtree::getLeft() const
{ 
	return this->left;
} 

long double Quadtree::getRight() const
{
	return this->right;
} 

long double Quadtree::getTop() const
{
	return this->top;
}

long double Quadtree::getBottom() const
{
	return this->bottom;
}

Quadtree* Quadtree::getChild(unsigned int indice)
{
	return this->mChildren[indice%4];
}

bool Quadtree::isParent() const
{
	return this->parent;
}

