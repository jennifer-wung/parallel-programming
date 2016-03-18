#include<iostream>
#include<stdio.h>
#include<stdlib.h>
#include<X11/Xlib.h>
#include <unistd.h>
#include<cmath>

#include "xwindow.h"
using namespace std;

XWindowDisplay::XWindowDisplay()
{
	x_pos = 0;
	y_pos = 0;
}

void XWindowDisplay::initGraph(int y_resol, int x_resol, double scal)
{
	height = y_resol;
	width = x_resol;
	scale = scal;
	//char *window_name="N-Body", *display_name=NULL;
	unsigned long valuemask = 0;
    	XGCValues values;
	//XSizeHints size_hints;
    	Pixmap bitmap;
	XSetWindowAttributes attr[1];

	/* open connection with the server */ 
	display = XOpenDisplay(NULL);
	if(display == NULL) {
		fprintf(stderr, "cannot open display\n");
		exit(1);
	}

	/* get display infomation */
	screen = DefaultScreen(display);

	/* set window position */
	//int x = 0;
	//int y = 0;

	/* border width in pixels */
	int border_width = 4;
	/* create window */
	window = XCreateSimpleWindow(display, RootWindow(display, screen), x_pos, y_pos, width, height, 
				border_width, BlackPixel(display, screen), WhitePixel(display, screen));

	
	//size_hints.flags = USPosition|USSize;
    	//size_hints.x = x_pos;
    	//size_hints.y = y_pos;
    	//size_hints.width = width;
    	//size_hints.height = height;
    	//size_hints.min_width = 300;
    	//size_hints.min_height = 300;
    	//XSetNormalHints (display, window, &size_hints);
    	//XStoreName(display, window, window_name);

    	/* create graphics context */
    	gc = XCreateGC (display, window, valuemask, &values);
    	XSetBackground (display, gc, BlackPixel (display, screen));
    	XSetLineAttributes (display, gc, 1, LineSolid, CapRound, JoinRound);
    	attr[0].backing_store = Always;
    	attr[0].backing_planes = 1;
    	attr[0].backing_pixel = BlackPixel(display, screen);
    	XChangeWindowAttributes(display, window, CWBackingStore | CWBackingPlanes | CWBackingPixel, attr);
	
	/* map(show) the window */
	XMapWindow(display, window);
	XSync(display, 0);

    	/* set color */
    	//screen_colormap = DefaultColormap(display, DefaultScreen(display));


	/* draw rectangle */
	XSetForeground(display,gc,BlackPixel(display,screen));
	XFillRectangle(display,window,gc,0,0,width,height);
	XFlush(display);

}

void XWindowDisplay::draw(int x, int y)
{
     	XSetForeground(display,gc,WhitePixel(display,screen));//grayColor.pixel);
      	XDrawPoint (display, window, gc, x, y);
}

void XWindowDisplay::clear()
{
	XSetForeground(display,gc,BlackPixel(display,screen));
	XFillRectangle(display,window,gc,0,0,width,height);
	XFlush(display);
}

void XWindowDisplay::flush()
{
	XFlush(display);	
}
