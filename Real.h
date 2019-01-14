/*======================================================================*
 *   TITLE:   REAL in T.K
 =======================================================================*/
#ifndef __MREAL_H
#define __MREAL_H


#ifdef SUN
#include <stream.h>
inline double sqr(double x) {return(x*x);}
#endif

#ifdef MS_TURBOC_V30
inline double sqr(double x) {return(x*x);}
#endif

#ifdef WIN32
inline double sqr(double x) {return(x*x);}
#endif


//!!typedef double t_t[2][3] ;           // transfo lin 2d complexe 

//typedef double t_precision;
#include <math.h>

#ifdef WIN32
inline double remainder( 
       		double x, 
		double y 
	) 
{
    double fm = fmod(x, y);
    return ( (fm < 0.5 * y ) ? fm : fm - y);
}
#endif


#define arctan(x) atan(x)

extern const double r_zero; 		// = 0.0;
extern const double r_un; 		// = 1.0;
extern const double r_deux; 		// = 2.0;
extern const double r_moins_un; 	// = ( r_zero - r_un );
extern const double aleph; 		// = 1.0E9;
extern const double pi; 		// = 3.14159265;
extern const double r_2pi; 		// = ( r_deux * pi );
extern double eps; 			// = 1.0E-5;


/*
 ! C_reel provides standart functions on real and is automaticaly generated.
 */

#ifdef is_small
#undef is_small
#endif

//#include "c_reel.h"

#endif


