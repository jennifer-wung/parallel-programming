#ifndef _X_WIN_H_
#define _X_WIN_H_

//#define X_RESN  1000    /* x resolution */
//#define Y_RESN  600    /* y resolution */
#include<X11/Xlib.h>

class XWindowDisplay
{
	public:
		XWindowDisplay();
		void initGraph(int y_resol, int x_resol, double scal);
		void draw(int x, int y);
		void clear();
		void flush();

	private:
		int width;
		int height;
		int x_pos;
		int y_pos;
		double scale;
		GC gc;
		Display *display;
		Window window;      //initialization for a window
		int screen;         //which screen
		Colormap screen_colormap;
		XColor grayColor;
};




#endif
