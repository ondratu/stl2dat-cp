/*======================================================================*
 *   TITLE: 3D Transformations  And Quaternions. 
 =======================================================================*/
#ifndef __ORT_H
#define __ORT_H

#include "vect.h"

/*-
 !   forward class declaration.   
 */
class stl_transf;

class stl_lin
/*-
 ! What:	Class of general 3D affine transformations.
 */

{
protected:
  stl_v l[4];	// The four column defining the transformation.
		// l[ 0]	Transformation of stl_v_i.	
		// l[ 1]	Transformation of stl_v_j.	
		// l[ 2]	Transformation of stl_v_k.	
		// l[ 3]	Transformation of stl_v_zero.	
public:
   /*-
    ! Constructors.
    !---------------
    */
    stl_lin() ;
    /*-
     ! Use:	Create the null transformation.
     */

    stl_lin( 
	const	stl_v& v1 , 	// Transformation of stl_v_i.
	const	stl_v& v2 , 	// Transformation of stl_v_j.
	const	stl_v& v3,	// Transformation of stl_v_k. 
	const	stl_v& v4 	// Transformation of stl_v_zero.
	); 
    /*-
     ! Use:	Create the transformation knowing 4 columns ( v4 is the 
     !		translative part).
     */

    stl_lin( 
	const	stl_v& v1 	// Translative part.
	); 
    /*-
     ! Use:	Create a translation.
     */

    stl_lin( 
	const	double r	// Scale factor. 
	); 
    /*-
     ! Use:	Create an uniform scaling around the origin.
     */

    stl_lin( 
	const	double alpha, 	// Angle of rotation.
	const	stl_v &axis, 	// Direction of axis of rotation.
	const	stl_v &origin 	// Origin of axis of rotation.
	); 
    /*-
     ! Use:	Create a rotation around a given axis.
     ! Precond:	axis.normalized()
     */

    stl_lin( 
	const	double r_x,	// Scale factor along X axis. 
	const	double r_y,	// Scale factor along Y axis. 
	const	double r_z 	// Scale factor along Z axis.
	); 
    /*-
     ! Use:	Create a non-uniform scaling around the origin knowing
     !		the scales factors along the tree major axes.
     */

    stl_lin( 
	const	stl_transf& r	// Transformation as a stl_transf. 
	); 
    /*-
     ! Use:	Create the transformation knowed as a stl_transf one.
     */

    stl_lin( 
	const	stl_lin& l_from, 	// Original coordinates system.
	const	stl_lin& l_to 		// Mapped coordinates system.
	); 
    /*-
     ! Use:	Create the transformation which map a given coordinates
     !		system  on another one. coordinates system are given
     !		as a matrix.
     ! Precond : !is_small(l_from.det())
     ! Remark:	l_to.det() may  be equal to 0.0
     ! Postcond : this * l_from[i] == l_to[i]
     */
    // DECL_SELF_CONSTRUCTOR(stl_lin)

   /*-
    ! Methods.
    !----------
    */

    /*-
     !  follow  methods for datas  interrogation.
     */

    void column_set( 
	const	int i,const 	// Index of the column.
		stl_v& v1 		// New value of the column.
	); 
    /*-
     ! Use:	Set a new value for a given vector column of the transformation.
     */

    void translation_set( 
	const	stl_v& v1 	// New value for the translative part.
	); 
    /*-
     ! Use:	Set a new value for the translative part of the transformation.
     */

    stl_v column( 
	const	int i	// Index of the column. 
	) const; 
    /*-
     ! Use:	Return the value of a column of the transformation.		
     */

    stl_v translation( ) const;
    /*-
     ! Use:	Return the translative part of the transformation.		
     */


    /*-
     ! acessing to datas with  "this"  viewed as a coordinates system.
     */

     void x_axis_set( 
	const	stl_v& v1	//  New value of X axis. 
	); 
    /*-
     ! Use:	Set a new value to the X acis.
     */

    void y_axis_set( 
	const	stl_v& v1 	// New value of Y axis.
	); 
    /*-
     ! Use:	Set a new value to the Y acis.		
     */

     void z_axis_set( 
	const	stl_v& v1	// New value of Z axis. 
	); 
    /*-
     ! Use:	Set a new value to the Z acis.		
     */

     void origin_set( 
	const	stl_v& v1 	// New value of origin.
	); 
    /*-
     ! Use:	Set a new value to the origin.		
     */

     stl_v x_axis() const;
    /*-
     ! Use:	Return the X axis.		
     */

     stl_v y_axis() const;
    /*-
     ! Use:	Return the Y axis.		

     */

     stl_v z_axis() const;
    /*-
     ! Use:	Return the Z axis.		
     */

     stl_v origin() const;
    /*-
     ! Use:	Return the origin.		
     */


    /*-
     !  Following overloaded operators are defined accordingly to 
     !  usual math. conventions:
     !
     !	boolean =  stl_lin == stl_lin 	: test equality of two stl_lin.
     !	stl_lin	=  stl_lin  + stl_lin 	: vectorial addition of stl_lin.
     !	stl_lin	=  stl_lin  - stl_lin	: vectorial susbtraction of stl_lin.
     !	stl_lin	=  stl_lin  * real	: vectorial scaling of a stl_lin
     !	t_v 	=  stl_lin  * t_v	: affine transformation of a t_v
     !	t_v 	=  stl_lin  * stl_lin	: compose two stl_lin
     !	stl_lin	+= stl_lin		: vectorial addition of stl_lin.
     !	stl_lin	-= stl_lin		: vectorial susbtraction of stl_lin.
     !	stl_lin	*= real			: vectorial scaling of a stl_lin
     !	stl_lin	*= stl_lin		: compose two stl_lin
     */

    bool operator == ( const stl_lin& l2) const;
    stl_lin operator + ( const stl_lin& l2)  const;
    stl_lin operator - ( const stl_lin& l2)  const;
    stl_lin operator * ( const double r)   const; 
    stl_v operator * ( const stl_v& v1)  const; 
    stl_lin operator * ( const stl_lin& l2)  const; 
    stl_lin& operator += ( const stl_lin& l2);
    stl_lin& operator -= ( const stl_lin& l2);
    stl_lin& operator *= ( const double r); 
    stl_lin& operator *= ( const stl_lin& l2); 
    
    double    det()        const;  
    /*-
     ! Use:	Calculate the determinant of the vectorial part of the
     !		transformation.		
     */

    double    trace()      const;
    /*-
     ! Use:	Calculate the trace of the vectorial part of the
     !		transformation.
     ! Remark:	Trcae is the sum of the diagonal terms of the vectorial
     !		part of the transformation.	
     */

    bool is_orthogonal_frame() const;
    /*-
     ! Use:	Test if the transformation viewed as a coordinates system is
     !		an orthogonal one.		
     */

    bool is_direct_orthogonal_frame() const;
    /*-
     ! Use:	Test if the transformation viewed as a coordinates system is
     !		an _direct orthogonal one.		
     */

    bool invertible() const;
    /*-
     ! Use:	Test if the transformation is invertible.		
     */

    stl_lin inverse() const;
    /*-
     ! Use:	Calculate the affine inverse of "this". 
     ! Precond:	invertible().                 
     */  

    stl_lin transpose()  const; 
    /*-
     ! Use:	Calculate the transpose of "this". 
     ! Warning:	Work on linear part only.                 
     */ 

    stl_lin change_basis( 
	const	stl_lin& change 	// The coordinates system in which we want
				// the corresponding stl_lin.
	) const; 
    /*-
     ! Use:	Return the stl_lin corresponding to "this". 
     ! 		in the new basis given by change.
     ! Warning:	Work on linear part only.  
     ! Precond:	change.invertible().            
     */  
                                        
    bool pivot( 
		stl_v& sol 	// Out: the solution of the equation.
	) ; 
    /*-
     ! Use:	Solve the equation A sol = B where A is the linear part
     !		of "this and -B the translative part.
     ! Remarks: this is modified : diagonalized by gauss method.
     ! Result:	true if a solution is found.
     */

    double quadratique( 
	const	stl_v& w1, 	// Left vector.
	const	stl_v& w2 	// Right vector.
	) const; 
    /*-
     ! Use:	Calculate the quadratic form R = W1.("this"* W2)  
     ! Warning:	Work with linear part only
     */
    
    stl_v linear ( 
	const	stl_v& v 	// vector to transform.
	) const; 
    /*-
     ! Use:	Applies the linear part of this on v.
     */

    stl_lin linear() const;
    /*-
     ! Use:	Return the linear part of this.
     */

    long to_q( 
		stl_transf &tr 	// The resulting transformation.
	) const; 
    /*-
     ! Use:	Convert in a stl_transf if "this" = A * B
     ! 		where A.is_orthogonal_coordinates system() and B is a scaling
     !		with a scale factor > 0
     ! Warning:	The scaling part B is lost during the conversion.
     ! Result:	1 ok; 0 result not found
     ! Precond:	invertible()
     */
    
    friend
    std::ostream& operator<<( 
		std::ostream& s, 	// Ostream where l1 will be written.
	const	stl_lin& l1 	// transformation to write in s.
	); 
    /*-
     ! Use:	Write a given transformation in an ostream.		
     ! remark:	Format is [ l[0] , l[1] , l[2] | l[3] ]
     */
};


class stl_quaternion
/*-
 ! What:	Class of quaternions for the private use of stl_transf only
 !		so everything are private here and no documentation is 
 !		provided.
 */

{
friend class stl_transf;	// Class of transformation based on quaternions.
private:
    double a;
    stl_v  v;

    stl_quaternion();
    /*-
     ! No public documentation.
     */

public:
    stl_quaternion( 
	const	double r1,
	const	double r2,
	const	double r3, 
	const	double r4
	); 
    /*-
     ! INTERNAL USE ONLY.
     */

    stl_quaternion( 
	const	double r1, 
	const	stl_v& v1
	); 
    /*-
     ! INTERNAL USE ONLY.
     */

    stl_quaternion( 
	const	stl_quaternion& v1
	); 
    /*-
     ! INTERNAL USE ONLY.
     */

    stl_quaternion& operator =( 
	const	stl_quaternion& v1
	); 
    /*-
     ! INTERNAL USE ONLY.
     */

    double real_part() const;	//!! a privatiser
    /*-
     ! INTERNAL USE ONLY.
     */

    stl_v  vector_part() const;	//!! a privatiser
    /*-
     ! INTERNAL USE ONLY.
     */

private:
    /*-
     ! No public documentation for following private methods.
     */

    void real_part_set( 
	const	double r 
	); 
    void vector_part_set( 
	const	stl_v& v1 
	); 
    bool product( 
	const	stl_quaternion & q, 
		stl_quaternion& r, 
		bool a_normaliser 
	) const; 

    long orstd();
    void normalize( );
    friend long stl_lin::to_q( 
		stl_transf & tr
	) const; 
    friend std::ostream& operator<<( 
		std::ostream& s, 
	const	stl_quaternion& q1
	); 
};

class stl_transf
/*-
 ! What:	Class of orthogonal transformations.
 */

{
private:
  stl_quaternion val[3];	// INTERNAL USE ONLY.

public:
   /*-
    ! Constructors.
    !---------------
    */
    stl_transf( 
	const	stl_quaternion& q1,
	const	stl_quaternion& q2, 
	const	stl_quaternion& q3
	); 
    /*-
     ! INTERNAL USE ONLY.
     */

    stl_transf();
    /*-
     ! Use:	Create  null transformation.
     */

    stl_transf( 
	const	stl_transf& v1
	); 

    stl_transf& operator =( 
	const	stl_transf& v1
	); 

    stl_transf( 
	const	double alpha,	// Angle of the rotation. 
	const	stl_v& axis, 	// Direction of the axis of  rotation.
	const	stl_v& origine 	// Origin of the axis of  rotation.
	); 
    /*-
     ! Use:	Create a transformation as a rotation around an axis.
     ! Precond:	axis.normalized()
     */

    stl_transf( 
	const	stl_v& vdp,	// Normal to the plane ( or null vector). 
	const	stl_v& vop 	// Origin of the plane or the point of symetry.
	); 
    /*-
     ! Use:	Create a transformation as a symetry on a point 
     !		( vdp.null() ) or a plane
     ! Precond:   vdp.null() || vdp.normalized()
     */

    stl_transf( 
	const	double r 	//  Scaling factor.
	); 
    /*-
     ! Use:	Create a transformation as a uniform scaling.
     */

    stl_transf( 
	const	stl_v& v 	// vector of translation.
	); 
    /*-
     ! Use:	Create the transform. as a translation
     */     
    
    stl_transf( 
	const	stl_v& op1,	// Origin of first  plane. 
	const	stl_v& vp1,	// Direction of first  plane. 
	const	stl_v& op2,	// Origin of second plane. 
	const	stl_v& vp2 	// Direction of second plane.
	); 
    /*-
     ! Use:	Create a transformation which map a given plane onto another
     !		one. More precisely : map op1 onto op2 ,map vp1 onto vp2
     !		and the intersection  line of the planes is invariant.
     ! Remark:	if vp1 and vp2 are colinear ( parallel planes )
     !  	if vp1 and vp2 are oposite translation of op1 to op2 and symetry
     !		else translation of op1 to op2
     ! Precond:  op[12].normalized().
     ! Postcond: (this * op10 == op2 ) && (  this * vp1 == vp2 )
     ! 		  && ( this.linear( vp1*vp2) = vp1*vp2 )
     */
            
    stl_transf( 
	const	stl_v& x1, 	// First  point.
	const	stl_v& x2, 	// Normalized vectors  orthogonal to x3.
	const	stl_v& x3, 	// Normalized vectors  orthogonal to x2.
	const	stl_v& x4, 	// Second  point.
	const	stl_v& x5,	// Normalized vectors  orthogonal to x6. 
	const	stl_v& x6 	// Normalized vectors  orthogonal to x5.
	); 
    /*-
     ! Use:	Create the transformation which map 3 points on 3 others points
     !		with a rotational linear part.
     ! Precond:	x[2356].nomalized()
     ! Postcond: ( this *x1 == x3 ) && (this.linear(x2) == x5)
     !		 && ( this.linear(x3) == x6).
     */

     stl_transf( 
	const	double alpha1,	// Angle of rotation around X axis. 
	const	double alpha2,	// Angle of rotation around Y axis. 
	const	double alpha3 	// Angle of rotation around Z axis.
	); 
    /*-
     ! Use:	Create the transformation defined by 3 angles of rotation
     !		around canonical axes (radian).
     */

    /*-
     ! Methods.
     !----------
     */

     stl_quaternion quaternion( 
	const	int i
	) const; 					//!! a supprimer 
    /*-
     ! INTERNAL USE ONLY.
     */



private:
     /*-
      ! follow some private method usefull for operator * implementation
      */
     stl_transf prim_compose( 
	const	stl_transf& q_sim2, 
		bool with_normalize 
	) const; 
    /*-
     ! No public documentation.
     */

public:
    /*-
     !  Following overloaded operators are defined accordingly to 
     !  usual math. conventions:
     !
     !	boolean =  stl_transf == stl_transf : test equality of two stl_transf.
     !	t_v 	=  stl_transf  * t_v	  : affine transformation of a t_v
     !	t_v 	=  stl_transf  * t_l	  : compose two stl_transf
     !	t_l	*= stl_transf		  : compose two stl_transf
     */

    bool operator == ( const stl_transf& q2) const;
    stl_transf operator * ( const stl_transf &q) const;
    stl_transf& operator *= ( const stl_transf &q);
    stl_v operator * ( const stl_v &v1) const;
    
    stl_v linear( 
	const	stl_v& v1	// Vector to transform. 
	) const; 
    /*-
     ! Use:	Transform a vector by the linear component of this.
     */

    stl_v scale( 
	const	stl_v& v1	// Vector to transform. 
	) const; 
    /*-
     ! Use:	Transform a vector by a scaling. or the scaling part
     !		of a transformation. 
     */

    stl_transf inverse( ) const;
    /*-
     ! Use:	Calculate the inverse transformation.
     ! Postcond:  a * a.inverse() == stl_q_id;
     */

    bool is_translation() const;
    /*-
     ! Use:	Return yes only if it's a pure translation.
     */

    bool is_rotation() const;
    /*-
     ! Use:	Return yes only if the linear part is a rotation
     !		( with eventual scaling and translation ).
     */

    bool is_plane_reflection() const;
    /*-
     ! Use:	Return yes only if the linear part is a plane reflexion
     !		( with eventual scaling and translation ).
     */

    bool is_point_symetry() const;
    /*-
     ! Use:	Return yes only if the linear part is a point symetry
     !		( with eventual scaling and translation ).
     */

    bool is_scaling() const;
    /*-
     ! Use:	Return yes if a scaling part is present.
     */

    int16_t sign() const;
    /*-
     ! Use:	Return the sign of the determinant.
     */
    
    double scale_factor() const;
    /*-
     ! Use:	Return the scale factor of the transformation.
     */

    stl_v translate_vector() const;
    /*-
     ! Use:	Return the translative part of the transformation. 		
     */

    double angle() const;
    /*-
     ! Use:	Return the angle of rotation part of the transformation.
     ! Precond:	is_rotation().
     */

    void angles_get( 
		double angles[ 3]	// Out: the three angles. 
	) const; 
    /*-
     ! Use:	Calculates the threes angles of rtation of the transformation 
     !		around Z,Y and Z axes respectively.
     ! Postcond:    -r_pi < angles[ 012] <= r_pi
     */

    stl_v axis() const;
    /*-
     ! Use:	Return the  direction of the axis for a rotation.
     !		the normale to  the plane for a plane reflexion.
     ! Precond: is_rotation() || is_plane_reflection()
     */

    stl_v center() const;
    /*-
     ! Use:	Return the center of a scaling or point symetry
     !		Return the origin of the axis for a rotation
     !		Return the origin of plane for a plane reflexion
     ! Precond: is_scaling() || is_point_symetry() || 
     !		is_rotation() || is_plane_reflection()
     */

    friend
    long stl_lin::to_q( 
		stl_transf & tr
	) const; 
    /*-
     ! See also: comment in declaration of this method in class stl_lin above.	
     */

    friend
    std::ostream& operator<<( 
		std::ostream& s,	// Ostream where r1 will be written. 
	const	stl_transf& r1 	// transformation to write in s.
	); 
    /*-
     ! Use:	Write a given transformation in an ostream.
     ! remark:	Format is [  q1 , q2 , q3 ] with qi = < q1.a , q1.v > 
     */		

};


/*-
 ! Some usefull constant.
 */
extern const stl_quaternion stl_q_zero;	// Null quaternion.
extern const stl_transf stl_r_id;		// Identity as a stl_transf.
extern const stl_lin stl_l_zero;		// Null transformation.
extern const stl_lin stl_l_id;		// Identity as a stl_lin.

#endif 

