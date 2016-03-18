#ifndef _thread_h_
#define _thread_h_

#include <pthread.h>

class Thread
{
	public:
		Thread();
		 ~Thread();

		int start();
		int join();
		int detach();
		pthread_t self();

		virtual void* run() = 0 ;

	private:
		pthread_t m_tid;
		int       m_running;
		int       m_detached;


};

#endif
