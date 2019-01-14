/* Automatic API file generated from /disq7/mtelc/std/calc/c_sea.cxx */
/* By user florian on Project /mtelc at 95/09/28 */
#ifndef __C_SEA_H
#define __C_SEA_H
extern DLL_EXPORT
long scalplan(
	const	stl_v& p1,	 // First point defining the plane.
	const	stl_v& p2,	 // Second point defining the plane.
	const	stl_v& p3,	 // third point defining the plane.
		   stl_v &oq,	 // Out: Origin of the plane.
		stl_v &wq,	 // Out: Direction of the plane.
		long &orient	 // Standart orientation of the plane
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
 */;


extern DLL_EXPORT
void Mpl_normalize(
	 const	stl_v& op,	 //      Non normalized origin of the plane.
	const	stl_v& vp,	//      Non normalized direction of the plane. 
		stl_v& opn,	 // Out: normalized origin of the plane.
		stl_v& vpn	 // Out: normalized direction of the plane.
	)
/*-
 ! Use:		calculate the normalized form of a plane.
 !		with the origin equal to the projection of < 0 0 0 > on the
 !		plane.
 ! Precond:	vp.normalized()
 *
 ! precondition(vp.normalized())
 */;

extern DLL_EXPORT
void scal_pl_p_drdr(
	const	stl_v &od1,	 // Point on the plane.
	    const	stl_v &vd1,	// First direction on the plane. 
	const	stl_v &vd2,	// Second direction on the plane. 
		   stl_v &op,	 // Out: Origin of the plane.
		stl_v &vp	 // Out: Direction of the plane.
	)
/*-
 ! Use:		Calculate a plane knowing a point and two directions 
 !		lying on it.
 *
 ! precondition(vd1.normalized())
 ! precondition(vd2.normalized())
 ! precondition(!vd1.paral(vd2))
 */;

extern DLL_EXPORT
long scaldr(
	const	stl_v& p1,	 // First point on the line.
	const	stl_v& p2,	// Second point on the line. 
		stl_v &od,	// Out: origin of the line. 
		stl_v &vd	 // Out: direction of the line.
	)
/*-
 ! Use:		Calculate a strait line passing through two given points.	
 ! Result: 	 0	Points are identical.
 !            	 1	Line well defined.	
 */;


extern DLL_EXPORT
int si_drdr_copl(
	const	stl_v& od1,	// Origin of the first line. 
	const	stl_v& vd1,	// Direction of the first line. 
	    const	stl_v& od2,	// Origin of the Second line. 
	const	stl_v& vd2,	 // Direction of the Second line.
		stl_v &p1	// Out: Intersection if found. 
	)
/*-
 ! Use:		Calculate intersection point between two coplanar lines.
 ! Remarks:	Resolution: intersection point is found for two line with
 !		an angle of 0.001 degrees.
 ! Result:	0  Parallel lines without intersection.
 !             	1  One point of intersection is found.
 !            	2  Parallel lines geometrically identical.
 */;


extern DLL_EXPORT
int s_dist_drdr(
	const	stl_v& od1,	// Origin of the first line. 
	const	stl_v& vd1,	// Direction of the first line. 
	    const	stl_v& od2,	// Origin of the Second line. 
	const	stl_v& vd2,	 // Direction of the Second line.
		double &dist,	// Out: Minimal distance between line 1 and 2. 
		stl_v &p1,	// Out: Point on line 1 the nearest of line 2. 
		stl_v &p2	 // Out: Point on line 2 the nearest of line 1.
	)
/*-
 ! Use:		Calculate de distance between to lines and the nearest points.
 ! Remark:	if lines are parallel p2 is the origine of line 2 and p1 is
 !		the projection of p2 on d1.
 ! Result: 	+1 general case : non coplanar lines
 !		+2 coplanar lines
 !		+3 parallel lines 
 !		+4 lines are identical.
 */;


extern DLL_EXPORT long sidrdr(
	const	stl_v& od1,	// Origin of the first line. 
	const	stl_v& vd1,	// Direction of the first line. 
	    const	stl_v& od2,	// Origin of the Second line. 
	const	stl_v& vd2,	 // Direction of the Second line.
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
 */;


extern DLL_EXPORT double spadr(
	const	stl_v& pt,	 // Point to be tested.
	const	stl_v& od,	// Origin of the line. 
	const	stl_v& vd	 // Direction of the line.
	)
/*-
 ! Use:		Calculate the distance between a point and a line.		
 ! Result:	The distance.
 */;


extern DLL_EXPORT double spaplan(
	const	stl_v& pt,	// Point to be tested. 
	const	stl_v& op,	// Origin of the plane. 
	const	stl_v& vp	 // Direction of the plane.
	)
/*-
 ! Use:		Calculate the signed distance between a point and a plane.	
 ! Remarks:	positive distance correspond to point lying in the half space
 !		corresponding to the direction of the plane.	
 ! Result:	The signed distance.	
 */;


extern DLL_EXPORT bool sapptpl(
	const	stl_v& pt,	// Point to be tested. 
	const	stl_v& op,	// Origin of the plane. 
	const	stl_v& vp	 // Direction of the plane.
	)
/*-
 ! Use:		Test if a given point is lying on a given plane		
 ! Remarks:	Test based on distance point/plane with a tolerance of +-eps.	
 */;


extern DLL_EXPORT
bool sapptdr(
	const	stl_v& pt,	 // Point to be tested.
	const	stl_v& od,	// Origin of the line. 
	const	stl_v& vd	 // Direction of the line.
	)
/*-
 ! Use:		Test if a given point is lying on a given line		
 ! Remarks:	Test based on distance point/line with a tolerance of +-eps.	
 */;


extern DLL_EXPORT
long sidrpl(
	const	stl_v& od,	// Origin of the line. 
	const	stl_v& vd,	 // Direction of the line.
	const	stl_v& op,	// Origin of the plane. 
	const	stl_v& vp,	 // Direction of the plane.
		stl_v &pt	 // .
	)
/*-
 ! Use:		Calulate the intersection point between a line and a plane.
 ! Result:	-1 line and plane are parallel, no intersection.
 !		 0 line included in the plane.
 !		 1 line and plane non parallel, one intersection point.
 */;


extern DLL_EXPORT long siplpl(
	const	stl_v& op1,	// Origin of the first plane. 
	const	stl_v& vp1,	 // Direction of the first plane.
	const	stl_v& op2,	// Origin of the second plane. 
	const	stl_v& vp2,	 // Direction of the second plane.
		stl_v &od,	// Origin of the intersection line. 
		stl_v &vd	 // Direction of the intersection line.
	)
/*-
 ! Use:		Calculate the intersection line of two planes.		
 ! Result:	-1 parallel planes, non intersection.
 !		 0 planes are geometrically identical.
 !		 1 non-parallel plane , interserction line is found.
 *
 ! precondition(vp1.normalized())
 ! precondition(vp2.normalized())
 */;

extern DLL_EXPORT
bool sadrpl(
	const	stl_v& od,	// Origin of the line. 
	const	stl_v& vd,	 // Direction of the line.
	const	stl_v& op,	// Origin of the plane. 
	const	stl_v& vp	 // Direction of the plane.
	)
/*-
 ! Use:		Test if a line is included in a given plane.	
 */;

extern DLL_EXPORT
double surface3pt(
	const	stl_v& p1,
	const	stl_v& p2,
	const	stl_v& p3
	)
/*-
 ! Use:		Calculates surface of a triangle given by 3 points.	
 */;


#endif
