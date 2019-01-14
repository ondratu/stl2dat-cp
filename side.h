/*======================================================================*
 *   TITLE: Relative Location in 1D, 2D and 3D space 
 =======================================================================*/

#ifndef __SIDE_H
#define __SIDE_H

#include "stdint.h"

enum t_side {
    left_side,
    right_side
};

/*
 ! What	: this enum define location in space in term of
 !	   main axes and/or half plane.
 !
 !  case stl_interval dim = 1 
 !
 !  Mleft   0-----------1	Mright
 !
 !  case t_carre dim = 2 
 !
 !		Mback
 !          1-----------3			Y
 !          |		|			|
 !  Mleft   |		|     Mright		0--X
 !          |		|
 !          |		|
 !          |		|
 !          0-----------2
 !		Mfront
 !
 !  case t_cube dim = 3 
 !						Z
 !						|  Y
 !						| /
 !		    Mtop			0_____X
 !
 !              3---------------7
 !             /.              /|
 !            / .  Mback      / |
 !           /  .            /  |
 !          1---------------5   |
 !          |   .           |   |
 !  Mleft   |   2...........|...6      Mright
 !          |  /            |  /
 !          | /   Mfront    | /
 !          |/              |/
 !          0---------------4
 !
 !		 Mbottom
 */



enum stl_locate
/*-
 ! What	: name regions around a box/point.
 */
{
    Minside	=0,
    Mleft	=1,
    Mright	=2,
    Mfront	=4,
    Mback	=8,
    Mbottom	=16,
    Mtop	=32,
    Mon		=64
};

const int32_t NotMon = ~Mon;

#define IS_INSIDE(x)(((x) & NotMon) == 0)
#define IS_LEFT(x)(((x) & Mleft) == Mleft)
#define IS_RIGHT(x)(((x) & Mright) == Mright)
#define IS_FRONT(x)(((x) & Mfront) == Mfront)
#define IS_BACK(x)(((x) & Mback) == Mback)
#define IS_BOTTOM(x)(((x) & Mbottom) == Mbottom)
#define IS_TOP(x)(((x) & Mtop) == Mtop)
#define IS_ON(x)(((x) & Mon) == Mon)
#define SAME_SIDE(x,y)(((x) & (y) & NotMon) !=0)
typedef int32_t stl_locator;




enum  stl_reg_side
/*-
 ! What:	Used to locate and object relatively to a planar region.
 */

{
    Mreg_unknown 	=0,	// side unknown: problem during the ray fire.
    Mreg_inside		=1,	// Inside the region.
    Mreg_outside	=2,	// outside the region.
    Mreg_on		=3,	// on the border of the region.
    Mreg_cross		=4,	// Crossing the region.
    Mreg_inside_on	=5,	// Inside the region + tangency.
    Mreg_outside_on	=6	// outside the region + tangency.
};


#endif




