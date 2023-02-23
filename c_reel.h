int_auto DLL_EXPORT void eps_change( 
	const	double new_eps 
	) 
/*-
 ! Use: 	Change value of eps
 */
;
int_auto DLL_EXPORT void eps_empile( 
	const	double new_eps 
	) 
/*-
 ! Use: 	Change value of eps, and save old value on top of the stack.
 */
;
int_auto DLL_EXPORT void eps_empile_rel( 
	const	double coef_eps 
	) 
/*-
 ! Use: 	Multiply current value of eps by coef_eps, 
 !		and save old value on top of the stack.
 */
;
int_auto DLL_EXPORT void eps_depile()
/*-
 ! Use: 	Restore previous value of eps from the stack.
 */
;
int_auto DLL_EXPORT double current_eps( ) 
/*-
 ! Use: 	Return the  current value of eps by.
 */
;
int_auto DLL_EXPORT bool is_small(const double r)
/*-
 ! Use: 	Return true if fabs(r) is lower than eps.
 */
;
int_auto DLL_EXPORT bool equal(const double r1, const double r2)
/*-
 ! Use: 	Return true if fabs(r1 - r2) is lower than eps.
 */
 ;
int_auto DLL_EXPORT long signe(const double r)
/*-
 ! Use: 	Return 	 1 if r >= 0
 !			-1 if r < 0.
 */
;
int_auto DLL_EXPORT double radian( 
	const	double t 	// angle in degrees.
	) 
/*-
 ! Use: 	Return t in radians.
 */
;
int_auto DLL_EXPORT double degre( 
	const	double t 	// angle in radians
	) 
/*-
 ! Use: 	Return t in degrees.
 */
;
int_auto DLL_EXPORT double angle( 
	const	double u,	// Adjacent side (cosinus)
	const	double v	// Opposite side (sinus)
	) 
/*-
 ! Use: 	Return the angle defined by u and v.
 !		u and v do not need to be normalised.
 ! Result:	in radian, between [0 .. 2pi[
 */
;
int_auto DLL_EXPORT double tang(const double r)
/*-
 ! Use: 	Return tangent of r.
 ! 		Return v_erreur if (is_small(cos(r)).
*/
;
int_auto DLL_EXPORT double arcsin(const double r)
/*-
 ! Use: 	Return arc sinus of r.
 ! 		Return v_erreur if (fabs(r) >= 1).
*/
;
int_auto DLL_EXPORT double arccos(const double r)
/*-
 ! Use: 	Return arc cosinus of r.
 ! 		Return v_erreur if (fabs(r) >= 1).
*/
;
int_auto DLL_EXPORT double substitue( 
	const	double y, 
	const	double a, 
	const	double b, 
	const	double c 
	) 
/*
 !  Use: 	Solve the linear equation ax + by + c = 0  knowing a,b,c,y.
 */
;
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

;
int_auto DLL_EXPORT void exchange( 
		double &x, 
		double &y 
	) 
/*-
 ! Use: 	Exchange x and y.
 */
;
int_auto DLL_EXPORT double Mmax( 
	const	double x, 
	const	double y 
	) 
/*-
 ! Use: 	Return maximum of x and y.
 */
;
int_auto DLL_EXPORT double Mmin( 
	const	double x, 
	const	double y 
	) 
/*-
 ! Use: 	Return minimum of x and y.
 */
;
int_auto DLL_EXPORT void roriente( 
		double &x, 
		double &y 
	) 
;

int_auto DLL_EXPORT double Mftrunc( 
	const	double r, 
 	const	double y 
	) 
/*
 ! Use:		Trunc the real on the given precision.
 !		result = r - i*y with result < r
 */
;
int_auto DLL_EXPORT double Mfround( 
	const	double r, 
	const	double y 
	) 
/*
 ! Use:		Round the real on the given precision.
 */
;
int_auto DLL_EXPORT long l_round( 
	const	double r 
	) 
/*-
 ! Use:		Equivalent of Pascal ROUND function.
 */
;
int_auto DLL_EXPORT long r_round32(const double r)
/*
 ! Use:		Portable round function (works on all machines).
 !		Prefer this function to Pascal ROUND function.
 */
;
int_auto DLL_EXPORT long r_trunc32(const double r)
/*
 ! Use:		Portable trunc function (works on all machines).
 !		Prefer this function to Pascal TRUNC function.
 */
;
int_auto DLL_EXPORT short r_round16(double r)
/*
 ! Use:		Portable round function (works on all machines).
 !		Prefer this function to Pascal ROUND function.
 */
;
int_auto DLL_EXPORT short r_trunc16(const double r)
/*
 ! Use:		Portable trunc function (works on all machines).
 !		Prefer this function to Pascal TRUNC function.
 */
;
int_auto DLL_EXPORT double Mfactorial( 
	const	long r
	)
;

int_auto DLL_EXPORT double Mpower( 
	const	double r1,
	const	double r0,
		double& r2
	)
;