/*======================================================================*
 *   TITLE: Class of point/vectors 2 and 3d.				*
 =======================================================================*/
#ifndef __VECT_H
#define __VECT_H

#include "iostream"
#include "side.h"
#include "real.h"


/*-
 !   forward declaration of stl_v 
 */
class stl_v;


class stl_c
/*-
 ! What:	2D point/vector.
 */

{
friend class stl_square;		// Class of 2d boxes.
friend class t_rep_local;	// Class of 2d coordinates system.
protected:
    double x[2];		// the two coordinates of the 2d point/vector.

public:
   /*-
    ! Constructors.
    !---------------
    */
    
    stl_c ();
    /*-
     ! Use:	 Create a null 2D point/vector coordinates ( 0.0, 0.0).
     */

    
    stl_c ( 
	const	double a, 	// X coordinate of the point/vector.
	const	double b 	// Y coordinate of the point/vector.
	); 
    /*-
     ! Use:	 Create a 2D point/vector knowing its 2 coordinates.
     */

    
    stl_c ( 
	const	double angle 	// polar coordinate.
	); 
    /*-
     ! Use:	Create a 2D normalized Direction knowing its polar
     !		coordinate.
     ! Precond:	( angle >= 0.0 ) && ( angle < r_2pi)
     */

    
    stl_c ( 
		double a[2]	// Array of the 2 coordinates
	); 
    /*-
     ! Use:	 Create a 2D point/vector knowing its 2 coordinates.
     */


   /*-
    !  Methods.
    !----------
    */

    /*-
     !  follow  methods for datas  interrogation.
     */

    
    void x_coord_set( 
	const	double r	// New X coordinate of the point/vector. 
	); 
    /*-
     ! Use:	Set a new value to the X coordinate of the point/vector.
     */

    
    void y_coord_set( 
	const	double r	// New Y coordinate of the point/vector. 
	); 
    /*-
     ! Use:	Set a new value to the Y coordinate of the point/vector.
     */

    
    double x_coord() const;
    /*-
     ! Use:	Return the X coordinate of the point/vector.
     */

    
    double y_coord() const;
    /*-
     ! Use:	Return the Y coordinate of the point/vector.
     */

    
    stl_v to_v( 
		int i_null = 2 	// Index of the null coordinate of
					// 3D point/vector.
	) const ; 
    /*-
     ! Use:	Make a 3D point/vector from a 2D one.
     ! Result:	i_null = 2 : stl_v(   x,   y, 0.0)
     !		i_null = 1 : stl_v(   x, 0.0,   y)
     !		i_null = 0 : stl_v( 0.0,   x,   y)
     ! Precond: (i_null >= 0 ) &&  (i_null  < 3 )
     */

    
    void align( 
	const	stl_c& ref, 	// Point/vector for reference
	const	int icoord	// Index of the coordinate for alignement. 
	); 
    /*-
     ! Use:	Align "this" on ref by identification of one of the 
     ! 		coordinates of "this".
     ! Remarks: coord = 0	X coordinate 
     !		coord = 1 	Y coordinate
     ! Precond: (coord >= 0 ) &&  (coord  < 2 )
     */

    
    stl_locator locate( 
	const	stl_c& pos, 		// Point to be located. 
	const	double acc = 0.0 	// Optionnal tolerance.
	) const; 
    /*-
     ! Use:	Return the location of a given point relative to "this" and
     !		a decomposition of the space in four regions delimited by 
     !		vertical and horizontal axes passing through "this"
     !		with a given accuracy.		
     !Remark:	accuracy is used to detect point ON one of the axes.		
     ! See also:	Side.h	
     */

    /*-
     !  Following overloaded operators are defined accordingly to 
     !  usual math. conventions:
     !
     !	bool	=   stl_c == stl_c	: Equality of 2 point/vector ( eps).
     !	bool	=   stl_c != stl_c	: Non equality of 2 point/vector ( eps).
     !	stl_c&	+=  stl_c		: Addition of two vectors.
     !	stl_c	=   stl_c + stl_c		: Addition of two vectors.
     !	stl_c	=  -stl_c		: Opposite  vector.
     !	stl_c&	-=  stl_c		: Sustraction of two vectors.
     !	stl_c	=   stl_c - stl_c		: Sustraction of two vectors.
     !	stl_c	=   stl_c * stl_c		: Vectorial product of two vectors.
     !	stl_c&	/=  double		: Scalar division of a vector.
     !	stl_c	=   stl_c / double	: Scalar division of a vector.
     !	stl_c&	*=  double		: Scalar product of a vector. 
     !	stl_c	=   stl_c * double	: Scalar product of a vector.
     ! 
     */

    
    bool operator ==( const stl_c& c2) const;
    
    bool operator !=( const stl_c& c2) const;
    
    stl_c& operator +=( const stl_c& c2);
    
    stl_c operator + ( const stl_c& c2)  const;
    
    stl_c  operator -() const;
    
    stl_c& operator -=( const stl_c& c2);
    
    stl_c operator - ( const stl_c& c2)  const;
    
    stl_c& operator *=( const double r);
    
    stl_c operator * ( const double r) const;
    
    stl_c& operator /=( const double r);
    
    stl_c operator / ( const double r) const;
    
    stl_c operator * ( const stl_c& c2)  const;

    
    bool null() const;
    /*-
     ! Use:	Test if a point/vector is a null one.		
     ! Remarks:	A tolerance of eps is used.	
     */

    
    bool normalized() const;
    /*-
     ! Use:	Test if a vector is normalized ( Euclidien Norm ).		
     ! Remarks:	A tolerance of eps is used.	
     */

    
    double norm() const;
    /*-
     ! Use:	Calculate the euclidien norm of a given vector.		
     */

    
    void normalize();
    /*-
     ! Use:	Normalize the vector.		
     ! Precond:	!null();
     ! Postcond: equal( 1.0, norm())
     */

    
    stl_c project( 
	const	stl_c& axis 	// Direction given by a  normalized vector. 
	) const; 
    /*-
     ! Use:	Project a vector on a given direction.		
     ! Warning:	"this" is viewed as a vector and not a point.	
     ! Result:	the projected vector.
     ! Precond: axis.normalized()
     */

    
    stl_c perp( 
	const	stl_c& axis 	// Direction given by a  normalized vector. 
	) const; 
    /*-
     ! Use:	Calculate the componante of the vector pependicular to the
     !		given direction.		
     ! Warning:	"this" is viewed as a vector and not a point.	
     ! Result:	the perpendicular vector.	
     ! Precond: axis.normalized()
     */

    
    void decompose( 
	const	stl_c& axis, 	// In:  Direction given by a  normalized vector.
		stl_c& p_para, 	// Out: parallel componante of "this".
		stl_c& p_ort 	// Out: perpendicular componante of "this".
	) const; 
    /*-
     ! Use:	Calculate parallel and perpendicular componantes on the given
     !		direction.
     ! Warning:	"this" is viewed as a vector and not a point.	
     ! Precond:	axis.normalized()
     ! Postcond: is_small(p_para.ps(p_ort))
     */

    
    bool paral( 
	const	stl_c& c2 	// Vector to compare with "this".
	) const; 
    /*-
     ! Use:	Test if the two vectors are parallel or not.		
     ! Remarks:	A tolerance of eps is used.	
     */

    
    double angle( 
	const	stl_c& c2 	// Second vector.
	) const; 
    /*-
     ! Use:	compute angle from  "this" to v2 .
     ! Result:	The angle.
     ! Postcond: (result >= 0.0) && ( result < r_2pi )
     */

    
    stl_locate entre( 
	const	stl_c& c1, 	// First point.
	const	stl_c& c2 	// Second point.
	) const; 
    /*-
     ! Use:	Test if "this" is lying on the segment defined by the
     !		two given point.		
     ! Result:	Mleft, Mon, Minside, Mright		
     ! Precond: v1 != v2
     */

    
    stl_c orth() const;
    /*-
     ! Use:	Calculate a direction perpendicular to "this"		
     ! Precond:	!null()
     ! Postcond: is_small(ps(result)) && result.normalized()
     */

    
    double ps( 
	const	stl_c& c2 	// Second vector.
	) const; 
    /*-
     ! Use:	Calculate the scalar product of two vectors.		
     */

    
    double distance( 
	const	stl_c& c2 	// Second point.
	) const ; 
    /*-
     ! Use:	Calculate the euclidien distance between  two points.		
     */

    
    stl_c comblin( 
	const	stl_c& c2, 	// Second vector.
	const	double r1, 	// Coeficient for "this".
	const	double r2 	// Coeficient for v2.
	) const; 
    /*-
     ! Use:		 a linear combination of two vectors.		
     */

    
    stl_c milieu( 
	const	stl_c& c2 	// second point.
	) const; 
    /*-
     ! Use:	Calculate the middle point of two points.		
     */

    
    double det( 
	const	stl_c& c2 	// Second vector.
	) const; 
    /*-
     ! Use:	Calculate  the determinants  "this" and v2
     */

    
    double angle() const;
    /*-
     ! Use:	Calculate  of vector "this" from X axis.
     ! Remark:	angle is in radian.
     */

    
    stl_c tourne( 
	const	double a 	// angle of rotation in radian.
	) const; 
    /*-
     ! Use:	Apply a rotation of the given angle and centered on origin
     !		to "this"		
     */

    
    stl_c cartesian_to_polar( ) const;
    /*-
     ! Use:	convert cartesian coordinates to polar one.		
     ! Result:	The equivalent point/vector in polar coordinates.
     */

    
    stl_c polar_to_cartesian() const;
    /*-
     ! Use:	convert polar coordinates to cartesian one.		
     ! Result:	The equivalent point/vector in cartesian coordinates.
     */
 
    
    friend std::ostream & operator << ( 
		std::ostream& s, 	// Ostream where v will be written.
	const	stl_c& c1 	// Point/vector to write in s.
	); 
    /*-
     ! Use:	Write a given point/vector in an ostream.
     ! remark:	Format is < x1 , x2 >		
     */

    
    friend bool lexicaly_ordered( 
	const	stl_c &c1, 	// First point/vector.
	const	stl_c &c2 	// Second point/vector.
	); 
    /*-
     ! Use:	Test if c1 and c2 are lexicaly_ordered ( coordinate x first).
     */

};

extern     
void exchange( 
		stl_c &c1, 	// First point/vector.
		stl_c &c2 	// Second point/vector.
	); 
/*-
 ! Use:	Exchange two point/vectors.		
 */


class stl_lin;

class stl_v
/*-
 ! What:	3D point/vector.
 */

{
friend class stl_lin;	// Class of 3d linear transformations.
friend class stl_box;	// Class of 3 boxes.
friend std::ostream& operator<<( 
		std::ostream& s, 
	const	stl_lin& l1 
	) ;
protected:
    double x[3];	// The three coordinates of the point/vector.

public:
   /*-
    ! Constructors.
    !---------------
    */
    
    stl_v ();
    /*-
     ! Use:	 Create a null point/vector coordinates ( 0.0, 0.0, 0.0)
     */

    
    stl_v ( 
	const	double a, 	// X coordinate of the point/vector.
	const	double b, 	// Y coordinate of the point/vector.
	const	double c 	// Z coordinate of the point/vector.
	); 
    /*-
     ! Use:	 Create a point/vector knowing its 3 coordinates.
     */

    
    stl_v ( 
		double a[3]	// Array of the 3 coordinates.
	); 
    /*-
     ! Use:	 Create a point/vector knowing its 3 coordinates.
     */

    
    stl_v ( 
	const	stl_c& a, 	// 2D point/vector ( XY coordinates).
	const	double b 	// Elevation ( Z coordinate).
	); 
    /*-
     ! Use:	Create a point/vector from a 2D point/vector and an 
     ! 		Elevation.
     */


   /*-
    ! Methods.
    !---------
    */

    /*-
     !  follow  methods for datas  interrogation.
     */
    
    void x_coord_set( 
	const	double r 	// New X coordinate of the point/vector.
	); 
    /*-
     ! Use:	Set a new value to the X coordinate of the point/vector.
     */

    
    void y_coord_set( 
	const	double r	// New Y coordinate of the point/vector. 
	); 
    /*-
     ! Use:	Set a new value to the Y coordinate of the point/vector.
     */

    
    void z_coord_set( 
	const	double r	// New Z coordinate of the point/vector. 
	); 
    /*-
     ! Use:	Set a new value to the Z coordinate of the point/vector.
     */

    
    double x_coord() const;
    /*-
     ! Use:	Return the X coordinate of the point/vector.
     */

    
    double y_coord() const;
    /*-
     ! Use:	Return the Y coordinate of the point/vector.
     */

    
    double z_coord() const;
    /*-
     ! Use:	Return the Z coordinate of the point/vector.
     */

	double coord(int n) const;
    /*-
     ! Use:	Return the nth coordinate of the point/vector.
     */

	void set_coord(int n, double r);
    /*-
     ! Use:	Set the nth coordinate of the point/vector.
     */

    stl_c to_c( 
		int i_null = 2 	// Index of the vanishing coordinate.
	) const ; 
    /*-
     ! Use:	Make a 2D point/vector from a 3D one with two of its
     !		three coordinates.		
     ! Result:	i_null = 2 : stl_c( x, y)
     !		i_null = 1 : stl_c( x, z)
     !		i_null = 0 : stl_c( y, z)
     ! Precond: (i_null >= 0 ) &&  (i_null  < 3 )
     */

    
    void align( 
	const	stl_v& ref, 	// Point/vector for reference.
	const	int icoord	// Index of the coordinate for alignement. 
	); 
    /*-
     ! Use:	Align "this" on ref by identification of one of the 
     ! 		coordinates of "this".
     ! Remarks: coord = 0	X coordinate 
     !		coord = 1 	Y coordinate
     !		coord = 2 	Z coordinate
     ! Precond: (coord >= 0 ) &&  (coord  < 3 )
     */

    
    stl_locator locate( 
	const	stl_v& pos, 		// Point to be tested.
	const	double acc = 0.0 	// Optionnal tolerance.
	) const; 
    /*-
     ! Use:	Return the location of a given point relative to "this" and
     !		a decomposition of the space in height regions delimited by 
     !		X, Y and Z plane  passing through "this"
     !Remark:	accuracy is used to detect point ON one of the plane.		
     ! See also:	Side.h	
     */


    /*-
     !  Following overloaded operators are defined accordingly to 
     !  usual math. conventions:
     !
     !	bool	=   stl_v == stl_v	: Equality of 2 point/vector ( eps).
     !	bool	=   stl_v != stl_v	: Non equality of 2 point/vector ( eps).
     !	stl_v&	+=  stl_v		: Addition of two vectors.
     !	stl_v	=   stl_v + stl_v		: Addition of two vectors.
     !	stl_v	=  -stl_v		: Opposite  vector.
     !	stl_v&	-=  stl_v		: Sustraction of two vectors.
     !	stl_v	=   stl_v - stl_v		: Sustraction of two vectors.
     !	stl_v	=   stl_v * stl_v		: Vectorial product of two vectors.
     !	stl_v&	/=  double		: Scalar division of a vector.
     !	stl_v	=   stl_v / double	: Scalar division of a vector.
     !	stl_v&	*=  double		: Scalar product of a vector. 
     !	stl_v	=   stl_v * double	: Scalar product of a vector. 
     !
     */

    
    bool operator ==( const stl_v& v2) const;
    
    bool operator !=( const stl_v& v2) const;
    
    stl_v& operator +=( const stl_v& v2);
    
    stl_v operator + ( const stl_v& v2) const;
    
    stl_v operator -() const;
    
    stl_v& operator -=( const stl_v& v2);
    
    stl_v operator - ( const stl_v& v2) const;
    
    stl_v operator * ( const stl_v& v2) const;
    
    stl_v& operator /=( const double r);
    
    stl_v operator / ( const double r) const;
    
    stl_v& operator *=( const double r);
    
    stl_v operator * ( const double r) const;

    
    bool null() const;
    /*-
     ! Use:	Test if a point/vector is a null one.		
     ! Remarks:	A tolerance of eps is used.	
     */

    
    bool normalized() const;
    /*-
     ! Use:	Test if a vector is normalized ( Euclidien Norm ).		
     ! Remarks:	A tolerance of eps is used.	
     */

    
    double norm() const;
    /*-
     ! Use:	Calculate the euclidien norm of a given vector.		
     */

    
    void normalize();
    /*-
     ! Use:	Normalize the vector.		
     ! Precond:	!null();
     ! Postcond: equal( 1.0, norm())
     */

	stl_v normaliz();

	stl_v ftrunc(int precision);
    /*-
     ! Use:	truncate to 10 power precision.		
     */
    
    stl_v project( 
	const	stl_v& axis 	// Direction given by a  normalized vector.
	) const; 
    /*-
     ! Use:	Project a vector on a given direction.		
     ! Warning:	"this" is viewed as a vector and not a point.	
     ! Result:	the projected vector.
     ! Precond: axis.normalized()
     */

    
    stl_v perp( 
	const	stl_v& axis 	// Direction given by a  normalized vector.
	) const; 
    /*-
     ! Use:	Calculate the componante of the vector pependicular to the
     !		given direction.		
     ! Warning:	"this" is viewed as a vector and not a point.	
     ! Result:	the perpendicular vector.	
     ! Precond: axis.normalized()
     */

    
    void decompose( 
	const	stl_v& axis ,	// In:  Direction given by a  normalized vector
		stl_v& p_para, 	// Out: parallel componante of "this".
		stl_v& p_ort 	// Out: perpendicular componante of "this".
	) const; 
    /*-
     ! Use:	Calculate parallel and perpendicular componantes on the given
     !		direction.
     ! Warning:	"this" is viewed as a vector and not a point.	
     ! Precond:	axis.normalized()
     ! Postcond: is_small(p_para.ps(p_ort))
     */

    
    bool paral( 
	const	stl_v& v2 	// Vector to compare with "this".
	) const; 
    /*-
     ! Use:	Test if the two vectors are parallel or not.		
     ! Remarks:	A tolerance of eps is used.	
     */

    
    double angle( 
	const	stl_v& v2,	// Second vector. 
	const	stl_v& v_ref 	// Angle is measured around this direction.
	) const; 
    /*-
     ! Use:	compute angle from  "this" to v2 around v_ref in radian.
     ! Result:	The angle.
     ! Precond:	v_ref.normalized()
     ! Postcond: (result >= 0.0) && ( result < r_2pi )
     */


    
    stl_locate entre( 
	const	stl_v& v1, 	// First point.
	const	stl_v& v2 	// Second point.
	) const; 
    /*-
     ! Use:	Test if "this" is lying on the segment defined by the
     !		two given point.
     ! Result:	Mleft, Mon, Minside, Mright		
     ! Precond: v1 != v2
     */

    
    stl_v orth() const;
    /*-
     ! Use:	Calculate a direction perpendicular to "this"		
     ! Precond:	!null()
     ! Postcond: is_small(ps(result)) && result.normalized()
     */

    
    double ps( 
	const	stl_v& v2 	// Second vector.
	) const; 
    /*-
     ! Use:	Calculate the scalar product of two vectors.		
     */

    
    double distance( 
	const	stl_v& v2	// Second point. 
	) const; 
    /*-
     ! Use:	Calculate the euclidien distance between  two points.		
     */

    
    stl_v comblin( 
	const	stl_v& v2 , 	// Second vector.
	const	double r1, 	// Coeficient for "this". 
	const	double r2	// Coeficient for v2. 
	) const; 
    /*-
     ! Use:		 a linear combination of two vectors.		
     */

    
    stl_v milieu( 
	const	stl_v& v2 	// second point..
	) const; 
    /*-
     ! Use:	Calculate the middle point of two points.		
     */

    
    double det( 
	const	stl_v& v2, 		// Second vector.
		int i_det = 0 	// Index of the determinant.
	) const; 
    /*-
     ! Use:	Calculate  one of the 3 determinants of the matrix built on
     !		"this" and v2
     ! Result:	i_det = 0   determinant obtained with first and second row.
     !		i_det = 1   determinant obtained with second and third row.
     !		i_det = 2  determinant obtained with first and third row.
     ! Precond:	(i_det >= 0) &&(i_det < 3);
     */

    
    double mixte( 
	const	stl_v& v2,	// Second Vector. 
	const	stl_v& v3 	// Third Vector.
	) const; 
    /*-
     ! Use:	Calculate  cross product of three vectors.
     */

	bool colinear(const stl_v & v1, const stl_v & v2);

    
    long orientation( ) const;
    /*-
     ! INTERNAL USE ONLY.
     */

    
    long orientation_for_compatibility_only( ) const; 
    /* 
     ! Use:	Ensure compatibility with previous version of Mtel
     !		for direction.
     ! Result:	-1 vector shall be reversed to obtain a standart orientation
     !		 0 null vector
     !		 1 vector is already in standart orientation.
     */
 
    
    bool read( 
	const	std::string & l 	// a given string.
	); 
    /* 
    ! Use: 	Fill itself with a CString containing: < x1 x2 x3 >
    ! Result:	true if the lig was correctly decoded, false otherwise.
    */
    
    
    std::string write() const;
    /*-
     ! Use:	Write the Point/vector in a string and return it		
     ! Result:	"< x1 x2 x3 >"
     */
    
    
    friend std::ostream & operator <<( 
		std::ostream& s, 	// Ostream where v will be written.
	const	stl_v& v 	// Point/vector to write in s.
	); 
    /*-
     ! Use:	Write a given point/vector in an ostream.		
     ! remark:	Format is < x1 , x2 , x3 >		
     */
    friend std::istream & operator >> (std::istream & is, stl_v & v);

    friend double determinant(	stl_v v[3] );
    friend double determinant4(	stl_v v[4] );
};

extern 
void exchange( 
		stl_v &v1, 	// First point/vector.
		stl_v &v2 	// Second point/vector.
	); 
    /*-
     ! Use:	Exchange two point/vectors.		
     */


/*-
 !  Some usefull constants.
 */

extern const stl_c stl_c_zero;	// = stl_c( 0.0, 0.0 );
extern const stl_c stl_c_i;		// = stl_c( 1.0, 0.0 );
extern const stl_c stl_c_j;		// = stl_c( 0.0, 1.0 );
extern const stl_c stl_c_un;	// = stl_c( 1.0, 1.0 );


extern const stl_v stl_v_zero;	// = stl_v( 0.0, 0.0, 0.0 );
extern const stl_v stl_v_i;		// = stl_v( 1.0, 0.0, 0.0 );
extern const stl_v stl_v_j;		// = stl_v( 0.0, 1.0, 0.0 );
extern const stl_v stl_v_k;		// = stl_v( 0.0, 0.0, 1.0 );
extern const stl_v stl_v_one;	// = stl_v( 1.0, 1.0, 1.0 );


#endif

