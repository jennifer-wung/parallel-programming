#include "thread.h"

void* runThread(void* arg)
{
	return ((Thread*)arg)->run();
}

Thread::Thread() : m_tid(0), m_running(0), m_detached(0) {}

Thread::~Thread()
{
	if (m_running == 1 && m_detached == 0)
		pthread_detach(m_tid);

	if (m_running == 1)
		pthread_cancel(m_tid);


}

int Thread::start()
{
	// runThread: The function must accept a void pointer to an object and return a void pointer to an object
	int result = pthread_create(&m_tid, NULL, runThread, this);
	if(result == 0) // thread was successful
		m_running = 1;
	
	return result;
}

int Thread::join()
{
	int result = -1;
	if (m_running == 1) {
		result = pthread_join(m_tid, NULL);
		if (result == 0)
			m_detached = 1; // thread is not detached
	}

	return result;

}

int Thread::detach()
{
	int result = -1;
	if (m_running == 1 && m_detached == 0) {
		result = pthread_detach(m_tid);
		if (result == 0)
			m_detached = 1;
	}

	return result;
}


pthread_t Thread::self() {
	return m_tid;
}
