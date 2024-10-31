/*======================================================================*
 *   TITLE: Class of vector/points 2 and 3d                             *
 =======================================================================*/
#include "stdafx.h"
#include "std_base.h"
#include "vect.h"
#include "ctype.h"
#include "c_reel.h"

#ifdef MSDOS
#define DLL_EXPORT_VAR2 extern DLL_EXPORT_VAR
#else
#define DLL_EXPORT_VAR2
#endif

const stl_c stl_c_zero	= stl_c( 0.0, 0.0 );
const stl_c stl_c_i		= stl_c( 1.0, 0.0 );
const stl_c stl_c_j		= stl_c( 0.0, 1.0 );
const stl_c stl_c_un	= stl_c( 1.0, 1.0 );

const stl_v stl_v_zero	= stl_v( 0.0, 0.0, 0.0 );
const stl_v stl_v_i		= stl_v( 1.0, 0.0, 0.0 );
const stl_v stl_v_j		= stl_v( 0.0, 1.0, 0.0 );
const stl_v stl_v_k		= stl_v( 0.0, 0.0, 1.0 );
const stl_v stl_v_one	= stl_v( 1.0, 1.0, 1.0 );

stl_c::stl_c ()
{
    x[0] = 0.0;
    x[1] = 0.0;
}

stl_c::stl_c ( 
	const	double a, 
	const	double b 
	) 
{
  x[0] = a;
  x[1] = b;
}

stl_c::stl_c( 
	const	double angle 
	) 
{
  x[0] = cos(angle);
  x[1] = sin(angle);
}

stl_c::stl_c ( 
		double a[2] 
	) 
{
  x[0] = a[0];
  x[1] = a[1];
}

void stl_c::x_coord_set( 
	const	double r 
	) 
{
    x[0] = r;
}

void stl_c::y_coord_set( 
	const	double r 
	) 
{
    x[1] = r;
}

double stl_c::x_coord() const
{
    return( x[0]);
}

double stl_c::y_coord() const
{
    return( x[1]);
}

stl_v stl_c::to_v( 
		int i_null 
	) const 
{
    precondition(i_null > -1);
    precondition(i_null < 3);

    if (i_null == 2)
	return (stl_v (x[0], x[1], 0.0));
    if (i_null == 1)
	return (stl_v (x[0], 0.0, x[1]));
    if (i_null == 0)
	return (stl_v (0.0, x[0], x[1]));
    return (stl_v_zero);
}

void stl_c::align( 
	const	stl_c& ref, 
	const	int coord 
	) 
{
    precondition(coord>=0);
    precondition(coord<2);
    x[coord] = ref.x[coord];
}

stl_locator stl_c::locate( 
	const	stl_c &pos, 
	const	double acc 
	) const 
{
    if (distance(pos)<acc)
    {
	return (Mon);
    }

    stl_locator result = 0;
    if (pos.x[0] < x[0]-acc)
	result |= Mleft;
    else if (pos.x[0] > x[0]+acc)
	result |= Mright;

    if (pos.x[1] < x[1]-acc)
	result |= Mfront;
    else if (pos.x[1] > x[1]+acc)
	result |= Mback;

    return (result);
}

bool stl_c::operator ==( 
	const	stl_c& c2 
	) const 
{
    return( (fabs( x[0]- c2.x[0]) + fabs( x[1] -c2.x[1]) < eps) );
}
bool stl_c::operator !=( 
	const	stl_c& c2 
	) const 
{
    return( !(fabs( x[0]- c2.x[0]) + fabs( x[1] -c2.x[1]) < eps) );
}

stl_c stl_c::operator +( 
	const	stl_c& c2 
	)  const 
{
    stl_c result = *this;
    result += c2;
    return(result);
}

stl_c stl_c::operator -( 
	const	stl_c& c2 
	)  const 
{
    stl_c result = *this;
    result -= c2;
    return(result);
}

stl_c stl_c::operator *( 
	const	double r 
	) const 
{
    stl_c result = *this;
    result *= r;
    return(result);
}

stl_c stl_c::operator /( 
	const	double r 
	) const 
{
    stl_c result = *this;
    result /= r;
    return(result);
}

stl_c stl_c::operator *( 
	const	stl_c& c2 
	)  const 
{
    stl_c result;
    result.x[0] = x[ 0] * c2.x[ 0] - x[ 1] * c2.x[ 1];
    result.x[1] = x[ 0] * c2.x[ 1] + x[ 1] * c2.x[ 0];
    return( result);
}

stl_c& stl_c::operator +=( 
	const	stl_c& c2 
	) 
{
    x[0] += c2.x[0];
    x[1] += c2.x[1];
    return(*this);
}

stl_c& stl_c::operator -=( 
	const	stl_c& c2 
	) 
{
    x[0] -= c2.x[0];
    x[1] -= c2.x[1];
    return(*this);
}

stl_c& stl_c::operator *=( 
	const	double r 
	) 
{
    x[0] *= r;
    x[1] *= r;
    return(*this);
}

stl_c& stl_c::operator /=( 
	const	double r 
	) 
{
    x[0] /= r;
    x[1] /= r;
    return(*this);
}

stl_c stl_c::operator -() const
{
    stl_c result;
    result.x[0] = -x[0];
    result.x[1] = -x[1];
    return(result);
}

bool stl_c::null() const
{
    return( (fabs( x[0]) + fabs( x[1]) < eps) );
}

bool stl_c::normalized() const
{
    return( is_small( norm() - 1.0 ));
}

double stl_c::norm() const
{
    return( sqrt( sqr(x[0]) + sqr(x[1]) ));
}

void stl_c::normalize()
{
    precondition(!null());

     *this /= norm();
}

bool stl_c::paral( 
	const	stl_c &c2 
	) const 
{
    return (is_small(det(c2)));
}

stl_c stl_c::project( 
	const	stl_c& axis 
	) const 
{
    precondition(axis.normalized());
    return( axis * ps(axis));
}

stl_c stl_c::perp(const stl_c& axis) const
{
    precondition(axis.normalized());
    return( *this - axis * ps(axis));
}

void stl_c::decompose( 
	const	stl_c &axis, 
		stl_c &p_para, 
		stl_c &p_ort 
	) const 
{
    stl_c c_loc = *this;
    if (axis == stl_c_zero)
    {
	p_para = stl_c_zero;
	p_ort = c_loc;
	return;
    }
    p_para = axis * (c_loc.ps(axis) / axis.ps( axis));
    p_ort = c_loc - p_para;
}

double stl_c::angle( 
	const	stl_c &c2 
	) const 
{
    if (null())
	return (::angle(c2.x[0], c2.x[1]));
    if (c2.null())
	return (::angle(x[0], 0.0 - x[1]));

    double norma = det(c2);
    double vv = ps(c2);
    if (!is_small(norma))
	return (::angle(vv, norma));
    if (vv >= 0.0)
	return (0.0);
    return (pi);
}

stl_locate stl_c::entre( 
	const	stl_c& c1, 
	const	stl_c& c2 
	) const 
{
    precondition(c1 != c2);

    stl_c v_l = c2 - c1;
    double longueur = v_l.norm();
    stl_c v_mes = *this - c1;
    v_l.normalize();
    double mes = v_mes.ps(v_l);
    if (is_small(mes) || is_small(mes - longueur))
	return (Mon);
    if (mes < 0.0)
	return (Mleft);
    if (mes > longueur)
	return (Mright);
    return (Minside);
}

stl_c stl_c::orth() const
{
    return( stl_c( -x[1],x[0]));
}

double stl_c::ps( 
 	const	stl_c& c2 
	) const 
{
    return( x[0]*c2.x[0] + x[1]*c2.x[1]);
}

double stl_c::distance( const stl_c& c2) const
{
    return( (c2 - *this).norm() );
}

stl_c stl_c::comblin( 
	const	stl_c& c2, 
	const	double r1, 
	const	double r2 
	) const 
{
    return( (*this)*r1 + c2*r2);
}

stl_c stl_c::milieu( 
	const	stl_c& c2 
	) const 
{
    return( ( *this + c2 ) / 2.0 );
}

double stl_c::det( 
 	const	stl_c& c2 
	) const 
{
    return( x[0]*c2.x[1] - x[1]*c2.x[0]);
}

double stl_c::angle() const
{
    return( ::angle( x[0], x[1]));
}

stl_c stl_c::tourne( const double a) const
{
    return ( (*this) * stl_c( a));
}

stl_c stl_c::cartesian_to_polar( ) const
{
     return(stl_c(norm(), angle()));
}

stl_c stl_c::polar_to_cartesian() const
{
    return( stl_c( x[0] * cos( x[1]), x[0] * sin( x[1]) ));
}


std::ostream& operator <<( 
		std::ostream& s, 
	const	stl_c& c1 
	) 
{
    return(s << '<' << c1.x[0] << ',' << c1.x[1] << '>');
}

bool lexicaly_ordered( 
	const	stl_c &c1, 
	const	stl_c &c2 
	) 
{
    /*
     ! this function is neede by dimensioning where
     ! lexical order is  significant.
     */
    stl_locator locate = c1.locate(c2);

    if (IS_LEFT(locate))
	return (false);
    if (!IS_RIGHT(locate) && IS_FRONT(locate))
	return (false);
    return (true);
}

void exchange( 
		stl_c &c1, 
		stl_c &c2 
	) 
{
    stl_c tmp = c2;
	c2  = c1;
	c1  = tmp;
}


stl_v::stl_v()
{
    x[0] = 0.0;
    x[1] = 0.0;
    x[2] = 0.0;
}

stl_v::stl_v( 
	const	double a, 
	const	double b, 
	const	double c 
	) 
{
    x[0] = a;
    x[1] = b;
    x[2] = c;
}

stl_v::stl_v( 
		double a[3] 
	) 
{
    x[0] = a[0];
    x[1] = a[1];
    x[2] = a[2];
}

stl_v::stl_v( 
	const	stl_c &a, 
	const	double b 
	) 
{
    *this = a.to_v();
    x[2] = b;
}

void stl_v::swapyz()
{
    std::swap(x[1], x[2]);
}

void stl_v::x_coord_set( 
	const	double r 
	) 
{
    x[0] = r;
}

void stl_v::y_coord_set( 
	const	double r 
	) 
{
    x[1] = r;
}

void stl_v::z_coord_set( 
	const	double r 
	) 
{
    x[2] = r;
}

double stl_v::x_coord() const
{
    return(x[ 0]);
}

double stl_v::y_coord() const
{
    return(x[ 1]);
}

double stl_v::z_coord() const
{
    return(x[ 2]);
}

double stl_v::coord(int n) const
{
    return(x[ n]);
}

void stl_v::set_coord(int n, double r)
{
	if (n >= 0 && n <= 2)
		x[ n] = r;
}

stl_c stl_v::to_c( 
		int i_null 
	) const 
{
    precondition(i_null >= 0);
    precondition(i_null < 3);

    if (i_null == 2)
	return (stl_c (x[0], x[1]));
    if (i_null == 1)
	return (stl_c (x[0], x[2]));
    if (i_null == 0)
	return (stl_c (x[1], x[2]));
    return( stl_c());
}

void stl_v::align( 
	const	stl_v& ref, 
	const	int coord 
	) 
{
    precondition(coord>=0);
    precondition(coord<3);
    x[coord] = ref.x[coord];
}

stl_locator stl_v::locate( 
	const	stl_v &pos, 
	const	double acc 
	) const 
{
    if (distance(pos)<acc)
    {
	return (Mon);
    }

    stl_locator result = 0;
    if (pos.x[0] < x[0]-acc)
	result |= Mleft;
    else if (pos.x[0] > x[0]+acc)
	result |= Mright;

    if (pos.x[1] < x[1]-acc)
	result |= Mfront;
    else if (pos.x[1] > x[1]+acc)
	result |= Mback;

    if (pos.x[2] < x[2]-acc)
	result |= Mbottom;
    else if (pos.x[2] > x[2]+acc)
	result |= Mtop;

    return (result);
}

bool stl_v::operator == ( 
	const	stl_v &v2 
	) const 
{
    return(fabs(x[0]-v2.x[0]) < eps && fabs(x[1]-v2.x[1]) < eps && fabs(x[2]-v2.x[2]) < eps);
}
bool stl_v::operator != ( 
	const	stl_v &v2 
	) const 
{
    return(!(fabs(x[0]-v2.x[0]) < eps && fabs(x[1]-v2.x[1]) < eps && fabs(x[2]-v2.x[2]) < eps));
}

stl_v stl_v::operator + ( 
	const	stl_v& v2 
	) const 
{
    stl_v result = *this;
    result += v2;
    return(result);
}

stl_v stl_v::operator - ( 
	const	stl_v& v2 
	) const 
{
    stl_v result = *this;
    result -= v2;
    return(result);
}

stl_v stl_v::operator *( 
	const	stl_v &v2 
	) const 
{
    stl_v result;
    result.x[0] = x[1] * v2.x[2] - x[2] * v2.x[1];
    result.x[1] = x[2] * v2.x[0] - x[0] * v2.x[2];
    result.x[2] = x[0] * v2.x[1] - x[1] * v2.x[0];
    return (result);
}

stl_v stl_v::operator * ( 
	const	double r 
	) const 
{
    stl_v result = *this;
    result *= r;
    return(result);
}

stl_v stl_v::operator / ( 
	const	double r 
	) const 
{
    stl_v result = *this;
    result /= r;
    return(result);
}

stl_v& stl_v::operator +=( 
	const	stl_v& v2 
	) 
{
    x[0] += v2.x[0];
    x[1] += v2.x[1];
    x[2] += v2.x[2];
    return(*this);
}

stl_v& stl_v::operator -=( 
	const	stl_v& v2 
	) 
{
    x[0] -= v2.x[0];
    x[1] -= v2.x[1];
    x[2] -= v2.x[2];
    return(*this);
}

stl_v& stl_v::operator *=( 
	const	double r 
	) 
{
    x[0] *= r;
    x[1] *= r;
    x[2] *= r;
    return(*this);
}

stl_v& stl_v::operator /=( 
	const	double r 
	) 
{
    x[0] /= r;
    x[1] /= r;
    x[2] /= r;
    return(*this);
}

stl_v stl_v::operator -() const
{
    stl_v result;
    result.x[0] = -x[0];
    result.x[1] = -x[1];
    result.x[2] = -x[2];
    return(result);
}

bool stl_v::null() const
{
    return ((fabs(x[0]) + fabs(x[1]) + fabs(x[2]) < eps));
}

bool stl_v::normalized() const
{
    return (is_small(norm() - 1.0));
}

double stl_v::norm() const
{
    return (sqrt(sqr(x[0]) + sqr(x[1]) + sqr(x[2])));
}

void stl_v::normalize()
{
	double n = norm();
	if (n != 0)
		*this /= n;
	else
		for( int i =0;i<3;i++)
			x[i] = 0.0;
}

stl_v stl_v::normaliz()
{
	stl_v result = *this;
	result.normalize();
	return(result);
}

double ftrunc( double v, int precision)
{
	double p = pow(10.0,precision);
	if (v < 0) p = -p;
	double a = v + p / 2;
	double b = fmod(a, p);
	double c = b + p / 2;
	double d = v - c;
	return(v - fmod(v + p / 2, p) + p / 2);
}

stl_v stl_v::ftrunc(int precision)
{
	stl_v result = *this;
	for( int i=0;i<3;i++)
		result.x[i] = ::ftrunc(x[i],precision);
	return(result);
}


stl_v stl_v::project( 
	const	stl_v& axis 
	) const 
{
    precondition(axis.normalized());
    return( axis * ps(axis));
}

stl_v stl_v::perp( 
	const	stl_v& axis 
	) const 
{
    precondition(axis.normalized());
    return( *this - axis * ps(axis));
}

void stl_v::decompose( 
	const	stl_v& axis, 
		stl_v &p_para, 
		stl_v &p_ort 
	) const 
{
    stl_v v_loc = *this;
    if (axis == stl_v_zero)
    {
	p_para = stl_v_zero;
	p_ort = v_loc;
	return;
    }
    p_para = axis * (v_loc.ps( axis) / axis.ps( axis));
    p_ort = v_loc - p_para;
}

bool stl_v::paral( 
	const	stl_v& v2 
	) const 
{    
  stl_v v3 = (*this) * v2;
  return( v3.null()); 
}

bool stl_v::colinear(const stl_v & v1, const stl_v & v2)
{
	stl_v d1 = (*this - v1).normaliz();
	stl_v d2 = (*this - v2).normaliz();
	return(d1 == d2 || d1 == -d2);
}


double stl_v::angle( 
	const	stl_v &v2, 
	const	stl_v &v_ref 
	) const 
{
    precondition(v_ref.normalized());
    stl_v norma = (*this) * v2;
    double vv = ps(v2);
    if (!norma.null())
	return (::angle(vv, v_ref.ps( norma)));
    if (vv >= 0.0)
	return (0.0);
    return (pi);
}
  
stl_locate stl_v::entre( 
	const	stl_v& v1, 
	const	stl_v& v2 
	) const 
{                           
    precondition(v1!=v2);

    stl_v  v_l       =   v2  - v1;
    double longueur  =  v_l.norm();
    stl_v  v_mes     = *this - v1;
    v_l.normalize();
    double mes = v_mes.ps(v_l);
    if (is_small(mes) || is_small(mes - longueur))
	return (Mon);
    if (mes < 0.0)
	return (Mleft);
    if (mes > longueur)
	return (Mright);
    return (Minside);
}

stl_v stl_v::orth() const
{
    precondition( !null());

    stl_v result;
    if ( is_small( x[2]) )
        result = stl_v( -x[1], x[0], 0.0);
    else
        result = stl_v( x[2], 0.0,-x[0]);
    result.normalize();

    postcondition( is_small( ps( result)));

    return(result);
}

double stl_v::ps( 
	const	stl_v& v2 
	) const 
{
    return( x[0]*v2.x[0] + x[1]*v2.x[1] + x[2]*v2.x[2]);
}

double stl_v::distance( 
 	const	stl_v& v2 
	) const 
{
    return( (v2 - *this).norm() );
}

stl_v stl_v::comblin( 
	const	stl_v& v2 ,const 
		double r1, 
	const	double r2 
	) const 
{
    return( (*this)*r1 + v2*r2);
}

stl_v stl_v::milieu( const stl_v& v2) const
{
    return( ( *this + v2 ) / 2.0);
}

double stl_v::det( 
	const	stl_v &v2, 
		int i_det 
	) const 
{
    precondition(i_det >= 0);
    precondition(i_det < 3);

    int i = i_det;
    int j = (i_det + 1) % 3;
    return (x[i] * v2.x[j] - x[j] * v2.x[i]);
}

double stl_v::mixte( 
	const	stl_v& v2, 
	const	stl_v& v3 
	) const 
{
   double result  = 0.0;
   for ( int i=0 ;i<3; i++)
   {
       result +=    (x[i] * v2.x[(i+1)%3] * v3.x[(i+2)%3]) 
		-   (v3.x[i]*v2.x[(i+1)%3]*x[(i+2)%3]);
   }
   return(result);
}

long stl_v::orientation() const
{
    if (null())
	return (0);

    bool std = true;

    if (!is_small(x[2]))
	std = (x[2] > 0.0);
    else if (!is_small(x[1]))
	std = (x[1] > 0.0);
    else if (!is_small(x[0]))
	std = (x[0] > 0.0);

    if (std)
	return (1);
    return (-1);
}

long stl_v::orientation_for_compatibility_only() const
{
    return (orientation());
}

bool stl_v::read( 
	const	std::string & l 
	) 
{
/*
    stl_v tmp_v;
    int i = l.skip(' ');
    if (!(l.get(i++) == '<'))
	return (false);
    i = l.skip(' ', i);
    i16 cpt_ind = -1;
    while ((l.get(i) != '>') && (i <= l.len()))
    {
	int i2 = l.find(CString(" >"), i);

	if (i2 <= 0)
	    return (false);
	if (!isdigit(l.get(i)) && l.get(i) != '-')
	    return (false);
	CString l2 = l.mid(i, i2 - 1);

	cpt_ind++;
	if (cpt_ind > 2)
	    return (false);

	tmp_v.x[cpt_ind] = l2.rval();
	i = i2;

	if (l.get(i) == ' ')
	{
	    if (i + 1 > l.len())
		return (false);
	    i = l.skip(' ', ++i);
	}
	else if (!(l.get(i) == '>'))
	    return (false);

    }
    *this = tmp_v;
	*/
    return (true);
}

std::string stl_v::write() const
{
	return " ";
/*
  return(CString("< ") + atos(x[0]) + " " +
                       atos(x[1]) + " " +
                       atos(x[2]) +
         CString(" >") );
*/
}

double small0(double v)
{
	return(is_small(v) ? 0 : v);
}

std::ostream & operator << ( 
		std::ostream & s, 
	const	stl_v &v 
	) 
{
//	char buf[80];
//	sprintf(buf,"%.4g %.4g %.4g", small0(v.x[0]), small0(v.x[1]), small0(v.x[2]));
//    return (s << buf);
	return(s <<  small0(v.x[0]) << " " << small0(v.x[1]) << " " << small0(v.x[2]));
}

void exchange( 
		stl_v &v1, 
		stl_v &v2 
	) 
{
    stl_v tmp = v2;
	v2  = v1;
	v1  = tmp;
}

double determinant(	stl_v v[3] )
{
//	cout << "determinant <" << v[0] << "> <" << v[1] << "> <" << v[2] << ">" << endl;
	return(v[0].x[0] * (v[1].x[1] * v[2].x[2] - v[2].x[1] * v[1].x[2]) +
           v[1].x[0] * (v[2].x[1] * v[0].x[2] - v[0].x[1] * v[2].x[2]) +
           v[2].x[0] * (v[0].x[1] * v[1].x[2] - v[1].x[1] * v[0].x[2]));
}

std::istream & operator >> (std::istream & is, stl_v & v)
{
	is >> v.x[0];
	is >> v.x[1];
	is >> v.x[2];
	return is;
}

#define DIM 4

typedef double ligne[DIM];
typedef ligne matrice[DIM];

//Cette fonction permet de supprimer la ligne et la colonne qui ne sert pas pour le calcul du d?terminant.
void copiesauflc (matrice source, matrice dest, int dim, int ligavirer)
{
	int l,c,ld = 0;
	for (l = 0; l<dim; l++)
		if (l != ligavirer)
		{
			for(c=1;c<dim;c++)
				dest[ld][c-1] = source[l][c];
			ld++;
		}
}


//Dans cette fonction, il appara?t une r?currence permettant de calculer les diff?rents d?terminants form?s par la matrice principale.
double determinant(matrice m, int dim)
{
	matrice sous_m;
	int l,signe = 1;
	double det = 0;
	if (dim == 1)
		return (m[0][0]);
	for (l=0; l<dim; l++)
	{
		copiesauflc(m, sous_m, dim, l);
		det += signe * m[l][0] * determinant(sous_m, dim-1);
		signe =- signe;
	}
	return(det);
}

double determinant4(stl_v v[4] )
{
//	cout << "determinant <" << v[0] << "> <" << v[1] << "> <" << v[2] << ">" << endl;
	matrice mat;
	int l,c;
	for (l = 0; l<4; l++)
	{
		for (c = 0; c<3; c++)
			mat[l][c] = v[l].x[c];
		mat[l][3] = 1;
	}
	return(determinant(mat,4));
}

