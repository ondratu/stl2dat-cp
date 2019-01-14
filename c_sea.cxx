/*======================================================================*
 *   TITLE: API  for Lines and Plane  Calculations.			*
 =======================================================================*/
//#include "stdafx.h"
#include "std_base.h"
#include "c_reel.h"
#include "ort.h"

//#include "xtrace.h"

long scalplan( 
	const	stl_v& p1, 	// First point defining the plane.
	const	stl_v& p2, 	// Second point defining the plane.
	const	stl_v& p3, 	// third point defining the plane.
   		stl_v &oq, 	// Out: Origin of the plane.
		stl_v &wq, 	// Out: Direction of the plane.
		long &orient 	// Standart orientation of the plane
				// -1: opposite
				//  o: Plane undefined
				//  1: Standard orientation.
	) 
/*-
 ! Use:		Calculate the plane passing through three given points.		
 ! Result:	Status on calculation:	
 ! 			 1: Plane well defined.
 ! 			 0: At least 2 points are identical: no plane
 ! 			-1: three aligned points: no plane.
 */

{
    stl_v uq = p2 - p1;
    stl_v vq = p3 - p1;
    if (uq.null() || vq.null())
	return (0);

    uq.normalize();
    vq.normalize();
    wq = uq * vq;

    if (wq.null())
	return (-1);
    wq.normalize();
    orient = wq.orientation();
    if (orient < 0)
	wq = -wq;
    oq = wq * p1.ps(wq);
    return (1);
}

void Mpl_normalize( 
 	const	stl_v& op, 	//      Non normalized origin of the plane.
	const	stl_v& vp,	//      Non normalized direction of the plane. 
		stl_v& opn, 	// Out: normalized origin of the plane.
		stl_v& vpn 	// Out: normalized direction of the plane.
	) 

/*-
 ! Use:		calculate the normalized form of a plane.
 !		with the origin equal to the projection of < 0 0 0 > on the
 !		plane.
 ! Precond:	vp.normalized()
 */

{
    precondition(vp.normalized());
    vpn = vp;
    opn = vpn * op.ps(vpn);
}

void scal_pl_p_drdr( 
	const	stl_v &od1, 	// Point on the plane.
    	const	stl_v &vd1,	// First direction on the plane. 
	const	stl_v &vd2,	// Second direction on the plane. 
   		stl_v &op, 	// Out: Origin of the plane.
		stl_v &vp 	// Out: Direction of the plane.
	) 
/*-
 ! Use:		Calculate a plane knowing a point and two directions 
 !		lying on it.
 */

{
    precondition(vd1.normalized());
    precondition(vd2.normalized());
    precondition(!vd1.paral(vd2));

    stl_v wq = vd1 * vd2;
    wq.normalize();
    if (wq.orientation() < 0)
	wq = -wq;
    vp = wq;
    op = wq * od1.ps(wq);
}

long scaldr( 
	const	stl_v& p1, 	// First point on the line.
	const	stl_v& p2,	// Second point on the line. 
		stl_v &od,	// Out: origin of the line. 
		stl_v &vd 	// Out: direction of the line.
	) 
/*-
 ! Use:		Calculate a strait line passing through two given points.	
 ! Result: 	 0	Points are identical.
 !            	 1	Line well defined.	
 */

{
    vd = p2 - p1;
    od = p1;
    if (!vd.null())
    {
	vd.normalize();
	if (vd.orientation() < 0)
	    vd = -vd;
	od -= (vd * p1.ps(vd));
	return (1);
    }
    return (0);
}

int si_drdr_copl( 
	const	stl_v& od1,	// Origin of the first line. 
	const	stl_v& vd1,	// Direction of the first line. 
    	const	stl_v& od2,	// Origin of the Second line. 
	const	stl_v& vd2, 	// Direction of the Second line.
		stl_v &p1	// Out: Intersection if found. 
	) 
/*-
 ! Use:		Calculate intersection point between two coplanar lines.
 ! Remarks:	Resolution: intersection point is found for two line with
 !		an angle of 0.001 degrees.
 ! Result:	0  Parallel lines without intersection.
 !             	1  One point of intersection is found.
 !            	2  Parallel lines geometrically identical.
 */

{

    double det_ref = 0.0;
    int i_ref = 0;
    for (int i = 0; i < 3; i++)
    {
	double det = vd1.det(vd2, i);
	if (fabs(det) > fabs(det_ref))
	{
	    i_ref = i;
	    det_ref = det;
	}
    }
    if (!is_small(det_ref))
    {
	p1 = od1 + (vd1 * ((od2 - od1).det(vd2, i_ref) / det_ref));
	return (1);
    }

    stl_v p2 = od2;
    p1 = od1 + (vd1 * vd1.ps(od2 - od1));
    p1 = p1.milieu(p2);
    if (is_small(p1.distance(p2)))
	return (2);
    return (0);
}

int s_dist_drdr( 
	const	stl_v& od1,	// Origin of the first line. 
	const	stl_v& vd1,	// Direction of the first line. 
    	const	stl_v& od2,	// Origin of the Second line. 
	const	stl_v& vd2, 	// Direction of the Second line.
		double &dist,	// Out: Minimal distance between line 1 and 2. 
		stl_v &p1,	// Out: Point on line 1 the nearest of line 2. 
		stl_v &p2 	// Out: Point on line 2 the nearest of line 1.
	) 
/*-
 ! Use:		Calculate de distance between to lines and the nearest points.
 ! Remark:	if lines are parallel p2 is the origine of line 2 and p1 is
 !		the projection of p2 on d1.
 ! Result: 	+1 general case : non coplanar lines
 !		+2 coplanar lines
 !		+3 parallel lines 
 !		+4 lines are identical.
 */

{
    stl_v dir = vd1 * vd2;
    if (dir.null())
    {
	p1 = od1 + (vd1 * vd1.ps(od2 - od1));
	p2 = od2;
	dist = p1.distance(p2);
	if (is_small(dist))
	    return (4);
	return (3);
    }

    /*-
     ! on travail dans le plan passant par
     ! < 0  0 0> de dir DIR
     */
    stl_v p;
    dir.normalize();
    stl_v dec1 = dir * od1.ps(dir);
    stl_v dec2 = dir * od2.ps(dir);
    int result = si_drdr_copl(od1 - dec1, vd1, od2 - dec2, vd2, p);
    if (result > 0)
    {
	p1 = p + dec1;
	p2 = p + dec2;
	dist = p1.distance(p2);
	if (is_small(dist))
	    result = 2;
    }
    return (result);
}

long sidrdr( 
	const	stl_v& od1,	// Origin of the first line. 
	const	stl_v& vd1,	// Direction of the first line. 
    	const	stl_v& od2,	// Origin of the Second line. 
	const	stl_v& vd2, 	// Direction of the Second line.
		stl_v &p1	// Out: Intersection if found. 
	) 
/*-
 ! INTERNAL USE ONLY.
 ! Use:		Calculate the intersection point between 2 lines.		
 ! Warning:	Coplanearity is tested with a tolerance of 0.01 for
 !		compatibility with previous versions
 ! Remarks:	If the lines are geometrically identical the point returned
 !		is the middle of the origins.	
 ! Result:	-2 : non coplanar lines , no intersections
 !		-1 : parallel lines, no intersections
 !		 0 : line geometrically identical. 
 !		 1 : Intersection found.
 */

{

    stl_v p2;
    double dist;
    int result = s_dist_drdr(od1, vd1, od2, vd2, dist, p1, p2);
    p1 = p1.milieu(p2);
    switch (result)
    {
    case 1:
	if ((dist < 0.01))
	    result = 1;
	else
	    result = -2;
	break;
    case 2:
	result = 1;
	break;
    case 3:
	result = -1;
	break;
    case 4:
	result = 0;
	break;
    default:;
	break;
    }
    return (result);
}

double spadr( 
	const	stl_v& pt, 	// Point to be tested.
	const	stl_v& od,	// Origin of the line. 
	const	stl_v& vd 	// Direction of the line.
	) 
/*-
 ! Use:		Calculate the distance between a point and a line.		
 ! Result:	The distance.
 */

{
    return (((pt - od) * vd).norm());
}

double spaplan( 
	const	stl_v& pt,	// Point to be tested. 
	const	stl_v& op,	// Origin of the plane. 
	const	stl_v& vp 	// Direction of the plane.
	) 
/*-
 ! Use:		Calculate the signed distance between a point and a plane.	
 ! Remarks:	positive distance correspond to point lying in the half space
 !		corresponding to the direction of the plane.	
 ! Result:	The signed distance.	
 */

{
    return( vp.ps(pt - op));
}

bool sapptpl( 
	const	stl_v& pt,	// Point to be tested. 
	const	stl_v& op,	// Origin of the plane. 
	const	stl_v& vp 	// Direction of the plane.
	) 
/*-
 ! Use:		Test if a given point is lying on a given plane		
 ! Remarks:	Test based on distance point/plane with a tolerance of +-eps.	
 */

{
    return(  is_small( spaplan( pt, op, vp)));
}

bool sapptdr( 
	const	stl_v& pt, 	// Point to be tested.
	const	stl_v& od,	// Origin of the line. 
	const	stl_v& vd 	// Direction of the line.
	) 
/*-
 ! Use:		Test if a given point is lying on a given line		
 ! Remarks:	Test based on distance point/line with a tolerance of +-eps.	
 */

{
    return( ( vd*( pt - od)).null());
}

long sidrpl( 
	const	stl_v& od,	// Origin of the line. 
	const	stl_v& vd, 	// Direction of the line.
	const	stl_v& op,	// Origin of the plane. 
	const	stl_v& vp, 	// Direction of the plane.
		stl_v &pt 	// .
	) 
/*-
 ! Use:		Calulate the intersection point between a line and a plane.
 ! Result:	-1 line and plane are parallel, no intersection.
 !		 0 line included in the plane.
 !		 1 line and plane non parallel, one intersection point.
 */

{
    /*-
     ! P sur D alors P = OD + l.VD
     !     ! x1         ! x0
     !  VD ! y1     OD  ! y0
     !     ! z1         ! z0
     !  P sur Q alors ux + vy + wz + t = 0
     !  donc          uxo + vyo + wzo + t = (-l)(x1u + y1v + z1w )
     ! L1 est le parallelisme entre le plan et la droite
     */
    long result = -1;
    stl_v odop = op - od;
    double l1 = vd.ps(vp);
    if (is_small(l1))
    {
	if (is_small(odop.ps(vp)))
	    result = 0;
    }
    else
    {
	pt = od + (vd * (odop.ps(vp) / l1));
	result = 1;
    }
    return (result);
}

long siplpl( 
	const	stl_v& op1,	// Origin of the first plane. 
	const	stl_v& vp1, 	// Direction of the first plane.
	const	stl_v& op2,	// Origin of the second plane. 
	const	stl_v& vp2, 	// Direction of the second plane.
		stl_v &od,	// Origin of the intersection line. 
		stl_v &vd 	// Direction of the intersection line.
	) 
/*-
 ! Use:		Calculate the intersection line of two planes.		
 ! Result:	-1 parallel planes, non intersection.
 !		 0 planes are geometrically identical.
 !		 1 non-parallel plane , interserction line is found.
 */

{
     /*------------------------------------------------------------*
      ! On cherche un point sur la droite en resolvant              !
      ! l'intersection des deux plans et d'un plan perpendiculaire  !
      ! aux deux passant par l'origine.                             !
      ! Ce dernier plan est defini par VPS(VDD,X) = 0. On trouve    !
      ! donc directement l'origine de la droite au sens  de MTEL    !
      *-------------------------------------------------------------*/
    precondition(vp1.normalized());
    precondition(vp2.normalized());
    long result = -1;
    vd = vp1*vp2;
    if( ! vd.null())
    {
        vd.normalize();
        double l1 = vp1.ps(op1);
        double l2 = vp2.ps(op2);
        double l3 = vp1.ps(vp2);
        double det =  1.0 - sqr( l3);
        double a1 = ( l1  - l2* l3) / det;
        double a2 = ( l2 - l3* l1) / det;
        od =  ( vp1 * a1 ) + ( vp2 * a2);
        result = 1;
    }
    else
    {
        if( is_small( vp1.ps(op1 - op2))) result = 0;
    }
    return( result);
}

bool sadrpl( 
	const	stl_v& od,	// Origin of the line. 
	const	stl_v& vd, 	// Direction of the line.
	const	stl_v& op,	// Origin of the plane. 
	const	stl_v& vp 	// Direction of the plane.
	) 
/*-
 ! Use:		Test if a line is included in a given plane.	
 */

{
    return(  is_small( vd.ps(vp)) && is_small( op.ps(vp) - od.ps(vp)));
}


double surface3pt(
	const	stl_v& p1,
	const	stl_v& p2,
	const	stl_v& p3
	)
/*-
 ! Use:		Calculates surface of a triangle given by 3 points.	
 */
{
	if (p1 == p2 || p1 == p3 || p2 == p3)
		return(0);
	double b = p1.distance(p2);
	double h = spadr(p3,p1,(p2-p1).normaliz());
	double result = (b * h) / 2;
	return(result);
}

