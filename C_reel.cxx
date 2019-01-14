/*======================================================================*
 *   TITLE: Real                                                        *
 =======================================================================*/
#include "stdafx.h"
#include "std_base.h"
#include "real.h"

//#define TRACE_MTH_ARGS
//#include "xtrace.h"

const double v_erreur	= 999999.99E-9;
const double r_zero	= 0.0;
const double r_un	= 1.0;
const double r_deux	= 2.0;
const double r_moins_un=  - 1.0 ;
const double aleph	= 1.0E9;
const double pi	= 3.14159265358979323846;  // see math.h
const double r_2pi	= ( 2.0 * pi );

double eps = 1.0E-5;

static double eps_pile[ 10 ];
static int eps_niv = 0;


int_auto DLL_EXPORT void eps_change( 
	const	double new_eps 
	) 
/*-
 ! Use: 	Change value of eps
 */
{
    eps = new_eps;
}

int_auto DLL_EXPORT void eps_empile( 
	const	double new_eps 
	) 
/*-
 ! Use: 	Change value of eps, and save old value on top of the stack.
 */
{
    precondition(eps_niv>=0);
    precondition(eps_niv<10);

    eps_pile[eps_niv] = eps;
    eps_niv++;
    eps = new_eps;
}

int_auto DLL_EXPORT void eps_empile_rel( 
	const	double coef_eps 
	) 
/*-
 ! Use: 	Multiply current value of eps by coef_eps, 
 !		and save old value on top of the stack.
 */
{
    precondition(eps_niv>=0);
    precondition(eps_niv<10);

    eps_pile[eps_niv] = eps;
    eps_niv++;
    eps = eps * coef_eps;
}

int_auto DLL_EXPORT void eps_depile()
/*-
 ! Use: 	Restore previous value of eps from the stack.
 */
{
    precondition(eps_niv>0);
    precondition(eps_niv<=10);
    eps_niv--;
    eps = eps_pile[eps_niv];
}

int_auto DLL_EXPORT double current_eps( ) 
/*-
 ! Use: 	Return the  current value of eps by.
 */
{
    precondition(eps_niv>=0);
    precondition(eps_niv<10);

    return ( eps );
}

int_auto DLL_EXPORT t_Mbool is_small(const double r)
/*-
 ! Use: 	Return true if fabs(r) is lower than eps.
 */
{
    return ((fabs(r) < eps));
}

int_auto DLL_EXPORT t_Mbool equal(const double r1, const double r2)
/*-
 ! Use: 	Return true if fabs(r1 - r2) is lower than eps.
 */
{
    return (is_small(r2 - r1));
}

int_auto DLL_EXPORT long signe(const double r)
/*-
 ! Use: 	Return 	 1 if r >= 0
 !			-1 if r < 0.
 */
{
    if (r >= 0.0)
	return (1);
    return (-1);
}

int_auto DLL_EXPORT double radian( 
	const	double t 	// angle in degrees.
	) 
/*-
 ! Use: 	Return t in radians.
 */

{
    return (t * ( pi / 180.0 ) );
}

int_auto DLL_EXPORT double degre( 
	const	double t 	// angle in radians
	) 
/*-
 ! Use: 	Return t in degrees.
 */

{
    return (t  * ( 180.0 / pi ));
}

int_auto DLL_EXPORT double angle( 
	const	double u,	// Adjacent side (cosinus)
	const	double v	// Opposite side (sinus)
	) 
/*-
 ! Use: 	Return the angle defined by u and v.
 !		u and v do not need to be normalised.
 ! Result:	in radian, between [0 .. 2pi[
 */
{
    if (u == 0.0)
    {
	if (v < 0.0)
	    return (1.5 * pi);
	return (0.5 * pi);
    }

    if (u <= 0.0)
	return (pi + atan(v / u));
    if (v < 0.0)
	return (r_2pi + atan(v / u));
    return (atan(v / u));
}

int_auto DLL_EXPORT double tang(const double r)
/*-
 ! Use: 	Return tangent of r.
 ! 		Return v_erreur if (is_small(cos(r)).
*/
{
    double c = cos(r);
    if (is_small(c))
	return (v_erreur);
    return (sin(r) / c);
}

int_auto DLL_EXPORT double arcsin(const double r)
/*-
 ! Use: 	Return arc sinus of r.
 ! 		Return v_erreur if (fabs(r) >= 1).
*/
{
    //maj 2935 le test a is_small est trop  brutale donc test a 10e-12
    if ((r > 1.0) || (r < -1.0))
	return (v_erreur);
    if (fabs(r - 1.0) >= 1E-12)
	return ((pi / 2.0) - 2.0 * atan(sqrt((1.0 - r) / (r + 1.0))));
    if ((r > 0.0))
	return (pi / 2.0);
    return (pi * 1.5);
}

int_auto DLL_EXPORT double arccos(const double r)
/*-
 ! Use: 	Return arc cosinus of r.
 ! 		Return v_erreur if (fabs(r) >= 1).
*/
{
    // maj 2935 le test a is_small est trop  brutale donc test a 10e-12
    if ((r > 1.0) || (r < -1.0))
	return (v_erreur);
    if (fabs(r - 1.0) >= 1E-12)
	return (2.0 * atan(sqrt((1.0 - r) / (r + 1.0))));
    if ((r > 0.0))
	return (0.0);
    return (pi);
}

int_auto DLL_EXPORT double substitue( 
	const	double y, 
	const	double a, 
	const	double b, 
	const	double c 
	) 
/*
 !  Use: 	Solve the linear equation ax + by + c = 0  knowing a,b,c,y.
 */

{
    precondition(!is_small(a));
    return (-((y * b + c) / a));
}

int_auto DLL_EXPORT short eq2d( 
	const	double a, 
	const	double b, 
	const	double c, 
		double &x1, 	// First or unique solution
		double &x2 	// Second solution
	) 
/*
 ! Use: 	Compute real solutions (when existing) of equation
 !			   2
 !			a x  + b x + c = 0
 !
 ! Result:	-2 when equation is like  c = 0
 !		 0 when no real solution
 !		 1 when 1 solution
 !		 2 when 2 solutions
 ! PRECONDITION: None
 ! Warning:	Distinction between case 1 and 2 depends on the value of eps.
 */

{
    /*-
     !  treat separately cases where A=B=0, A=0 and other cases
    */
    if (a == 0.0)
    {
	if (!is_small(b))
	{
	    x1 = -c / b;
	    return (1);
	}
	return (-2);
    }

    double bb = b / 2.0;
    double delta = sqr(bb) - a * c;

    if (delta < -eps)
	return (0);

    double ecart = sqrt(fabs(delta)) / a;
    x1 = - bb / a;

    if (is_small(ecart))
	return (1);

    x2 = x1 + ecart;
    x1 -= ecart;
    return (2);
}

int_auto DLL_EXPORT void exchange( 
		double &x, 
		double &y 
	) 
/*-
 ! Use: 	Exchange x and y.
 */
{
    double tp = x;
    x = y;
    y = tp;
}

int_auto DLL_EXPORT double Mmax( 
	const	double x, 
	const	double y 
	) 
/*-
 ! Use: 	Return maximum of x and y.
 */
{
    if (x > y)
	return (x);
    return (y);
}

int_auto DLL_EXPORT double Mmin( 
	const	double x, 
	const	double y 
	) 
/*-
 ! Use: 	Return minimum of x and y.
 */
{
    if (x > y)
	return (y);
    return (x);
}

int_auto DLL_EXPORT void roriente( 
		double &x, 
		double &y 
	) 
/*-
 ! Use: 	Exchange x and y if x > y.
 ! POSTCONDITION: x <= y
 */
{
    if (x > y)
	exchange(x, y);
}

static double r_int32( 
	const	double r 
	) 
{
    const long maxint32 = 2147483647;
    const long minint32 = -2147483647;

    if ((r > maxint32))
	return (maxint32);
    if ((r < minint32))
	return (minint32);
    return (r);
}

static double r_int16(const double r)
{
    const long maxint16 = 32767;
    const long minint16 = -32767;

    if ((r > maxint16))
	return (maxint16);
    if ((r < minint16))
	return (minint16);
    return (r);
}

int_auto DLL_EXPORT double Mftrunc( 
	const	double r, 
 	const	double y 
	) 
/*
 ! Use:		Trunc the real on the given precision.
 !		result = r - i*y with result < r
 */

{
    precondition(y>0.0);
    return (( r > 0 ) ? r - fmod(r,y) : r - fmod(r,y) -y );
}

int_auto DLL_EXPORT double Mfround( 
	const	double r, 
	const	double y 
	) 
/*
 ! Use:		Round the real on the given precision.
 */

{
    precondition(y>0.0);
    return (r - remainder(r,y));
}

int_auto DLL_EXPORT long round( 
	const	double r 
	) 
/*-
 ! Use:		Equivalent of Pascal ROUND function.
 */
{
    return ((long)Mfround( r_int32(r), 1.0));
}

int_auto DLL_EXPORT long r_round32(const double r)
/*
 ! Use:		Portable round function (works on all machines).
 !		Prefer this function to Pascal ROUND function.
 */

{
    return( (long)Mfround( r_int32(r), 1.0));
}

int_auto DLL_EXPORT long r_trunc32(const double r)
/*
 ! Use:		Portable trunc function (works on all machines).
 !		Prefer this function to Pascal TRUNC function.
 */
{
    return( (long)Mftrunc( r_int32(r), 1.0));
}

int_auto DLL_EXPORT short r_round16(double r)
/*
 ! Use:		Portable round function (works on all machines).
 !		Prefer this function to Pascal ROUND function.
 */
{
    return( (short)Mfround( r_int16(r), 1.0));
}

int_auto DLL_EXPORT short r_trunc16(const double r)
/*
 ! Use:		Portable trunc function (works on all machines).
 !		Prefer this function to Pascal TRUNC function.
 */
{
    return( (short)Mftrunc( r_int16(r), 1.0));
}

int_auto DLL_EXPORT double Mfactorial( 
	const	long r
	)

{
    double result = 1;
    long r1 = r;
    while ( r1 > 1 )
    {
	result *= r1;
	r1--;
    }
    return result;
}


int_auto DLL_EXPORT double Mpower( 
	const	double r1,
	const	double r0,
		double& r2
	)

{
    if ( is_small( r1 ) )
      if ( is_small(r0) )
	return( 1);
      else
	r2 = 1;
    else
    {
	if (r0 > 0)
	    if ( r1 > 0)
		r2 = exp(r1 * log( r0));
	    else
		r2 = exp(r1 * log(fabs( r0)));
	else
	{
	    double r2 = r0;
	    if ( !is_small( r1 - 1))
	    {
		int i;
		for (i = 2 ; i <= (long )r1 ; i++)
		{
		    r2 = r2 * r0;
		}
		if ( r1 <= 0)
		{
		    r2 = 1 / r2;
		}
	    }
	}
    }
    return 0;
}
