/*======================================================================*
 *   TITLE: 3D Transformations   & quaternions                          *
 *   Recap classification                                               *
 *   ---------------------                                              *
 *   -3 : Symetrie et Rotation Perp                                     *
 *   -2 : Symetrie par rapport a un point                               *
 *   -1 : Symetrie par rapport a un plan                                *
 *    1 : Rotation autour d'un axe                                      *
 *    2 : Identite                                                      *
 *    5 : Dilatation                                                    *
 *                                                                      *
 *   Explication                                                        *
 *   -----------                                                        *
 *     Une transformation de l'espace peut etre representee par         *
 *   un quaternion plus un entier                                       *
 *                                                                      *
 *     On peut trier les transformations othogonales par les dimensions *
 *   des leurs espaces propres .                                        *
 *     Les espaces propres peuvent etre de dimension 1 2 ou 3           *
 *     Les valeurs propres peuvent etre 1 -1 ou complexes               *
 *     Un espace propre de dim 1 -> valeur propre = -1 ou 1             *
 *     -  ------ ------ -- --- 2 -> valeur propre = -1    1 ou complexe *
 *     -- ------ ------ -- dim 3 -> valeur propre = -1 ou 1             *
 *     Ceci provient par exemple du fait que les polynomes ( comme le   *
 *     polynome caracteristique de la transfo ) a coefficients reels se *
 *     decompose en produit de polynomes lineaires et de polynomes de   *
 *     degre 2 a determinant negatif ( donc racine complexe ) .         *
 *                                                                      *
 *     Classification des transformations                               *
 *     ----------------------------------                               *
 *     3    : 1/ vp(3) = 1                                              *
 *               c'est l'ID                                             *
 *               correspond a Q[1].A = 2                                *
 *          : 2/ vp(3) = -1                                             *
 *               c'est -ID                                              *
 *               correspond a Q[1].A = -2                               *
 *     2,1  : 1/ vp(2) = 1  alors vp(1) = -1                            *
 *               il s'agit d'une reflexion                              *
 *               correspond a Q[1].A = -1 avec Q[0].A = 0 et            *
 *               Q[0].V axis de la reflexion.                            *
 *          : 2/ vp(2) = -1 alors vp(1) = 1                             *
 *               retournement autour d'un axis                           *
 *               correspond a Q[1].A = 1 et Q[0].A = 0 et Q[0].V = axis du
 *               retournement.                                          *
 *            3/ vp(2) complexe vp(1) = 1                               *
 *               rotation autour de Ep(1) de l'angle indique par le     *
 *               complexe                                               *
 *               Q[1].A = 1 Q[0].A = tgt(angle/2) Q[0].V axis de rotation*
 *            4/ vp(2) complexe vp(1) = -1                              *
 *               rotation uatour de Ep(1) de l'angle indique par le     *
 *               complexe puis reflexion par rapport a l'axis de rotation*
 *               Q[1].A =-3 Q[0].A = tgt(angle/2) Q[0].V axis de rotation*
 *     c'est tout .                                                     *
 *     Remarque: les espaces propres sont toujours perpendiculaires car *
 *     Une base des vecteurs propres doit etre othonormee               *
 *                                                                      *
 * Remarque finale et essentielle                                       *
 * ------------------------------                                       *
 * La classification et tous les calculs sont bases sur le fait que     *
 * l'application OPP de O+(3) vers O-(3) definie par                    *
 *                                                                      *
 *              F(f) = -f                                               *
 * possede les proprietes suivantes                                     *
 *                                                                      *
 * - involutive donc bijective                                          *
 * - Pour un vecteur donne V   F(f)[V] = vopp( f[V] )                   *
 *   ( propriete utilise pour QLSIM )                                   *
 * - Pour la composition : OPP est un homomorphisme sur O+(3)           *
 *   il suffit alors de voir les elements s de O-(3) sous la forme      *
 *   F(x) ou x = F(s) pour en deduire les formules de composition       *
 * - Remarquer que pour une rotation donne de vecteur V et scalaire A   *
 *   r2 = ( A , -V ) est la rotation oppose a R           ( cos , -sin )*
 *   r2 = ( -A,  V ) est la rotation oppose a R           ( cos , -sin )*
 *   r2 = ( -1/A, V ) est la rotation d'angle pi + angle  (-cos , -sin )*
 *   r2 = ( 1/A, V ) est la rotation d'angle pi - angle   (-cos ,  sin )*
 *   n'oubliez pas que 2*pi = 0 pour comprendre la commutativite de ces *
 *   operations                                                         *
 * Calcul : Image par similitude d'un vecteur ( QRIMAGE )               *
 * -------------------------------------------                          *
     si Q = A + V et Q' = A' + V'                                       *
     alors Q.Q' = [ A.A' - V.PS(V') ] + [ A.V' + A'.V + V ^ V' ]        *
     Q =  V                                                             *
     W vecteur                                                          *
                   -1                                                   *
     Evaluons Q*W*Q                                                     *
     puis l'oppose du resultat pour avoir l'image.                      *

     Q*W :                                                              *
     partie reelle:      - Q.V.PS(W)                                    *
     partie vectorielle: Q.V * W                                        *

     On a:                                                              *
          -1          _         2                                       *
     Q*W*Q    =  (Q*W*Q) / !!Q!!                                        *
                                                                        *
          2                                                             *
     !!Q!!   = 1                                                        *
                        _                                               

     On evalue donc Q*W*Q:                                              *

     Partie reelle:       VECTORIEL(Q.V,W) , Q.V)                       *

     Partie vectorielle: + VSCALAIRE(PSV(Q.V,W),Q.V)                    *
                         + VECTORIEL(Q.V,W),Q.V)                        *

     La partie reelle doit etre nulle ( THEOREME )                      *

     On developpe :                                                     *

     - PSV(Q.V,W)*Q.A                                                   *
     + PSV(VSCALAIRE(Q.A,W),Q.V)                                        *
     + PSV(VECTORIEL(Q.V,W),Q.V);                                       *

     avec la linearite de PSV :                                         *

     - PSV(Q.V,W)*Q.A                                                   *
     + PSV(W,Q.V)*Q.A                                                   *
     + PSV(VECTORIEL(Q.V,W),Q.V);                                       *

     VECTORIEL(Q.V,W) est perpendiculaire a Q.V                         *
     donc PSV(VECTORIEL(Q.V,W),Q.V) = 0                                 *

     il reste                                                           *

     - PSV(Q.V,W*//*)*Q.A                                                   *
     + PSV(W,Q.V)*Q.A                                                   *

     qui est bien nul                                                   *

     Developpons la partie vectorielle :                                *

       VSCALAIRE(Q.A, VSCALAIRE(Q.A,W))                                 *
     + VSCALAIRE(Q.A, VECTORIEL(Q.V,W))
     + VSCALAIRE(PSV(Q.V,W),Q.V) 
     - VECTORIEL(VSCALAIRE(Q.A,W),Q.V)
     - VECTORIEL(VECTORIEL(Q.V,W),Q.V)

     puis :
       VSCALAIRE(SQR(Q.A),W))
     + VSCALAIRE(Q.A, VECTORIEL(Q.V,W))
     + VSCALAIRE(PSV(Q.V,W),Q.V)
     - VSCALAIRE(Q.A,VECTORIEL(W,Q.V)
     - VECTORIEL(VECTORIEL(Q.V,W),Q.V)

     puis :

       VSCALAIRE(SQR(Q.A),W)) 
     + VSCALAIRE(Q.A, VECTORIEL(Q.V,W)) 
     + VSCALAIRE(PSV(Q.V,W),Q.V) 
     - VSCALAIRE(Q.A,VECTORIEL(W,Q.V)
     - VECTORIEL(VECTORIEL(Q.V,W),Q.V)

     puis:

       VSCALAIRE(SQR(Q.A),W))
     + VSCALAIRE(Q.A, VECTORIEL(Q.V,W))  \  2*VSCALAIRE(Q.A , VECTORIEL(Q.V,W))
     - VSCALAIRE(Q.A,VECTORIEL(W,Q.V)    /
     + VSCALAIRE(PSV(Q.V,W),Q.V)
     - VECTORIEL(VECTORIEL(Q.V,W),Q.V)

     On developpe le real produit vectoriel:

     VECTORIEL(VECTORIEL(Q.V,W),Q.V) =
     VSCALAIRE(PSV(Q.V,Q.V),W) - VSCALAIRE(PSV(Q.V,W),Q.V)

     On re-integre dans la formule

     + VSCALAIRE(PSV(Q.V,W),Q.V)
     - VSCALAIRE(PSV(Q.V,Q.V),W)
     + VSCALAIRE(PSV(Q.V,W),Q.V)

     Deux lignes se composent

     + 2*VSCALAIRE(PSV(Q.V,W),Q.V)
     - VSCALAIRE(PSV(Q.V,Q.V),W)

     reste (en opposant tout )
     W
     - 2 * VSCALAIRE(PSV(Q.V,W) , Q.V)

     Remarques:

     1/   Donc trois composantes suivant trois vecteurs perp.
          La premiere represente la partie Cosinus
          La deuxieme est la partie constante suivant l'axis de rotation*//*
          La derniere est la partie sinus.

     2/
     Q.A = NORME(Q.V) / TG(ANGLE/2)
     On a toujours NORME(Q.V) = 1
     Supposons Q.V perp a W :

     -------------------------------------------------------------

     Calcul : Image par rotation d'un vecteur
     ----------------------------------------
     si Q = A + V et Q' = A' + V'
     alors Q.Q' = [ A.A' - PSV(V,V') ] + [ A.V' + A'.V + V ^ V' ]

     Q = A + V
     W vecteur
                   -1
     Evaluons Q*W*Q
     Q*W :
     partie reelle:      - PSV(Q.V,W)
     partie vectorielle: VSCALAIRE(Q.A,W) +
                         VECTORIEL(Q.V,W)
     On a:
          -1          _         2
     Q*W*Q    =  (Q*W*Q) / !!Q!!
     On evalue donc Q*W*Q:

     Partie reelle:      - PSV(Q.V,W)*Q.A
                         + PSV( VSCALAIRE(Q.A,W) + VECTORIEL(Q.V,W) , Q.V)

     Partie vectorielle:   VSCALAIRE(Q.A, VSCALAIRE(Q.A,W) + VECTORIEL(Q.V,W))
                         + VSCALAIRE(PSV(Q.V,W),Q.V)
                         - VECTORIEL(VSCALAIRE(Q.A,W) + VECTORIEL(Q.V,W),Q.V)

     Rem: On divisera par la norme a la fin

     La partie reelle doit etre nulle ( THEOREME )

     On developpe :
     - PSV(Q.V,W)*Q.A
     + PSV(VSCALAIRE(Q.A,W),Q.V)
     + PSV(VECTORIEL(Q.V,W),Q.V);

     avec la linearite de PSV :
     - PSV(Q.V,W)*Q.A
     + PSV(W,Q.V)*Q.A
     + PSV(VECTORIEL(Q.V,W),Q.V);

     VECTORIEL(Q.V,W) est perpendiculaire a Q.V
     donc PSV(VECTORIEL(Q.V,W),Q.V) = 0

     il reste
     - PSV(Q.V,W)*Q.A
     + PSV(W,Q.V)*Q.A

     qui est bien nul
     Developpons la partie vectorielle :
       VSCALAIRE(Q.A, VSCALAIRE(Q.A,W))
     + VSCALAIRE(Q.A, VECTORIEL(Q.V,W))
     + VSCALAIRE(PSV(Q.V,W),Q.V)
     - VECTORIEL(VSCALAIRE(Q.A,W),Q.V)
     - VECTORIEL(VECTORIEL(Q.V,W),Q.V)

     puis :
       VSCALAIRE(SQR(Q.A),W))
     + VSCALAIRE(Q.A, VECTORIEL(Q.V,W))
     + VSCALAIRE(PSV(Q.V,W),Q.V)
     - VSCALAIRE(Q.A,VECTORIEL(W,Q.V)
     - VECTORIEL(VECTORIEL(Q.V,W),Q.V)

     puis :
       VSCALAIRE(SQR(Q.A),W))
     + VSCALAIRE(Q.A, VECTORIEL(Q.V,W))
     + VSCALAIRE(PSV(Q.V,W),Q.V)
     - VSCALAIRE(Q.A,VECTORIEL(W,Q.V)
     - VECTORIEL(VECTORIEL(Q.V,W),Q.V)

     puis:
       VSCALAIRE(SQR(Q.A),W))
     + VSCALAIRE(Q.A, VECTORIEL(Q.V,W))  \  2*VSCALAIRE(Q.A , VECTORIEL(Q.V,W))
     - VSCALAIRE(Q.A,VECTORIEL(W,Q.V)    /
     + VSCALAIRE(PSV(Q.V,W),Q.V)
     - VECTORIEL(VECTORIEL(Q.V,W),Q.V)

     On developpe le double produit vectoriel:
     VECTORIEL(VECTORIEL(Q.V,W),Q.V) =
     VSCALAIRE(PSV(Q.V,Q.V),W) - VSCALAIRE(PSV(Q.V,W),Q.V)

     On re-integre dans la formule
       VSCALAIRE(SQR(Q.A),W))
     + 2*VSCALAIRE(Q.A , VECTORIEL(Q.V,W))
     + VSCALAIRE(PSV(Q.V,W),Q.V)
     - VSCALAIRE(PSV(Q.V,Q.V),W)
     + VSCALAIRE(PSV(Q.V,W),Q.V)

     Deux lignes se composent
       VSCALAIRE(SQR(Q.A),W))
     + 2*VSCALAIRE(Q.A , VECTORIEL(Q.V,W))
     + 2*VSCALAIRE(PSV(Q.V,W),Q.V)
     - VSCALAIRE(PSV(Q.V,Q.V),W)

     reste
       VSCALAIRE( [SQR(Q.A) - NORME(Q.V)] , W)
     + 2 * VSCALAIRE(PSV(Q.V,W) , Q.V)
     + 2 * VSCALAIRE( Q.A , VECTORIEL(Q.V,W)

     Remarques:
     1/   Donc trois composantes suivant trois vecteurs perp.
          La premiere represente la partie Cosinus
          La deuxieme est la partie constante suivant l'axis de rotation
          La derniere est la partie sinus.

     2/
     Q.A = NORME(Q.V) / TG(ANGLE/2)
     On a toujours NORME(Q.V) = 1
     Supposons Q.V perp a W :

     ----------------------------------------------------------------

     SI Q.A = 1 = NORME(Q.V) Alors angle = 2*45 = 90 c'est
     le quart de tour or le qurt de tour c'est aussi le produit vectoriel
     et on remarque que la formule devient:

       VSCALAIRE( [SQR(Q.A) - NORME(Q.V)] , W)   = 0
     + 2 * VSCALAIRE(PSV(Q.V,W) , Q.V)           = 0 ( perp)
     + 2 * VSCALAIRE( Q.A , VECTORIEL(Q.V,W)
     le tout divise par NORME(Q)**2 ( qui vaut SQR(A) + 1 = 2)
     donne precisement le produit vectoriel.

     --------------------------------------------

     Si Q.A  = 0 la formule devient :
     -  W  ( divise par SQR(A) + 1 = 1 )
     c'est le demi tour de W autour de v
     ----------------------------------------------
     si W et Q.V sont colineaires
     le produit vectoriel est nul et la formule devient
      VSCALAIRE( [SQR(Q.A) - NORME(Q.V)] , W)
     + 2 * VSCALAIRE(PSV(Q.V,W) , Q.V)
     soit W = b.Q.V

      VSCALAIRE( [SQR(Q.A) - NORME(Q.V)] ,b.Q.V)
     + 2 * VSCALAIRE(b , Q.V)

     ( SQR(A) - 1 + 2 ) * Q.V = ( SQR(A) + 1 ) * Q.V
                              = ( NORME(Q)**2 ) * Q.V

     comme il fallait encore divise par ca on a bien que :
     Pour un vecteur colineaire a l'axis une rotation ne fait rien

 ----------------------------------------------------------------------
                 >        >                    3
  Dilatations: f(V) = k * V avec k element de R

               On ne peut pas composer une dilatation avec une rotation
  ou une symetrie, le resultat n'etant pas representable sous la forme
  quaternionique. On ne peut que les creer, et les appliquer au points
  QDEFDILAT, QDILAT


 *======================================================================*/

#include "stdafx.h"
#include "std_base.h"
#include "c_reel.h"
#include "ort.h"

//#define TRACE_MTH_ARGS

stl_quaternion::stl_quaternion()
{
    a = 0.0;
    v = stl_v_zero;
}

stl_quaternion::stl_quaternion( 
 	const	double r1, 
	const	double r2, 
	const	double r3, 
	const	double r4 
	) 
{
    a = r1;
    v = stl_v( r2,r3,r4);
}

stl_quaternion::stl_quaternion( 
	const	double r1, 
	const	stl_v& v1 
	) 
{
    a = r1;
    v = v1;
}

stl_quaternion::stl_quaternion(
	const	stl_quaternion& q
	) : a( q.a) , v( q.v)
{
}

stl_quaternion& stl_quaternion::operator = ( 
	const	stl_quaternion &q 
	) 
{
    a = q.a;
    v = q.v;
    return (*this);
}

double stl_quaternion::real_part() const
{
    return (a);
}

stl_v stl_quaternion::vector_part() const
{
    return (v);
}

void stl_quaternion::real_part_set( 
	const	double r 
	) 
{
    a = r;
}

void stl_quaternion::vector_part_set( 
	const	stl_v& v1 
	) 
{
    v = v1;
}

bool stl_quaternion::product( 
 	const	stl_quaternion & q, 
		stl_quaternion & r, 
		bool a_normaliser 
	) const 
{
    r.a = (a * q.a) - v.ps(q.v);
    r.v = (q.v * a) + (v * q.a) + (v * q.v);

    if (r.v.null())
    {
	/*-
	 ! Les axes sont / / et angles egaux
	 ! Les axes egaux et les A sont tous deux nuls
	 */
	    r.a = 0.0;
	return (false);
    }

    if (a_normaliser)
	r.normalize();
    return (true);
}

long stl_quaternion::orstd()
{
    long result = v.orientation();
    if (result < 0)
    {
	v = -v;
	a = -a;
    }
    return (result);
}

void stl_quaternion::normalize( )
{
    if (!v.null())
    {
	double norm = v.norm();
	v /= norm;
	a /= norm;
    }
}

ostream& operator<<( 
		ostream& s, 
	const	stl_quaternion& q1 
	) 
{
    return (s << '<' << q1.a << ',' << q1.v << '>');
}


stl_transf::stl_transf()
{
    val[0] = stl_q_zero;
    val[1] = stl_q_zero;
    val[2] = stl_q_zero;
}

stl_transf::stl_transf( 
	const	stl_quaternion& q1, 
  	const	stl_quaternion& q2, 
	const	stl_quaternion& q3 
	) 
{
    val[0] = q1;
    val[1] = q2;
    val[2] = q3;
}
                      
static double normalise_angle( 
	const	double val_a 
	) 
{
    if (is_small(val_a))
	return (0.0);		// 180 degre tg ~ x autour de 0
    if (is_small(val_a - 1.0))
	return (1.0);		// 90 degre ~ 0.0003 degre
    if (is_small(val_a + 1.0))
	return (-1.0);		// -90 degre ~ 0.0003 degre
    return (val_a);
}

stl_transf::stl_transf( 
	const	double alpha, 
	const	stl_v &axis, 
	const	stl_v &origine 
	) 
{
    /*-
     ! On peut avoir alpha = 0	 	-> implique IDENTITE
     ! On peut avoir alpha = 2*PI	 -> implique IDENTITE
     */
    precondition(axis.normalized());

    double alpha2 = alpha / 2.0;
    if (is_small(sin(alpha2)))
    {
	*this = stl_transf(stl_v_zero);
	return;
    }
    val[0].a = normalise_angle(cos(alpha2) / sin(alpha2));
    val[0].v = axis;
    val[1].a = 1.0;
    val[2].a = 1.0;
    val[1].v = origine - linear(origine);
}

stl_transf::stl_transf( 
	const	stl_v& vdp, 
	const	stl_v& vop 
	) 
{
    precondition(vdp.null() || vdp.normalized());
    if (vdp.null())
    {
	val[0] = stl_q_zero;
	val[1].a = -2.0;
	val[1].v = vop * 2.0;
    }
    else
    {
	val[0].v = vdp;
	val[0].a = 0.0;
	val[1].a = -1.0;
	val[1].v = vop - linear(vop);
    }
    val[2].a = 1.0;
}

stl_transf::stl_transf( 
	const	double r 
	) 
{
    *this = stl_r_id;
    val[2].a = r;
}

stl_transf::stl_transf( 
	const	stl_v& v 
	) 
{
    val[0] = stl_q_zero;
    val[1].a = 2.0;
    val[1].v = v;
    val[2].a = 1.0;
}

stl_transf::stl_transf( 
	const	stl_transf& v1
	)
{
    val [0] = v1.val[ 0];
    val [1] = v1.val[ 1];
    val [2] = v1.val[ 2];
}

stl_transf& stl_transf::operator =( 
	const	stl_transf& v1
	)
{
    val [0] = v1.val[ 0];
    val [1] = v1.val[ 1];
    val [2] = v1.val[ 2];
    return( *this);
}

stl_transf::stl_transf( 
	const	stl_v &op1, 
	const	stl_v &vp1, 
	const	stl_v &op2, 
	const	stl_v &vp2 
	) 
{
    precondition(vp1.normalized());
    precondition(vp2.normalized());

    stl_v vd = vp1 * vp2;

    if (vd.null())
    {
	if (vp1.ps(vp2) > 0.0)
	{
	    /*- Cas parallele de meme sens ( identite ) */
	    val[0] = stl_q_zero;
	    val[1].a = 2.0;
	    val[1].v = op2 - op1;
	    val[2].a = 1.0;
	}
	else
	{
	    // Cas parallele oppose ( demi_tour ) : si les origines sont
            // confondues. On prend un vecteur quelconque  non colineaire a VP1

	    vd = vp1 * ( op2 - op1 );

	    if (vd.null())
		vd = vp1.orth();

	    vd.normalize();
	    if (vd.orientation() < 0)
		vd = -vd;
	    val[0].a = 0.0;
	    val[0].v = vd;
	    val[1].a = 1.0;
	    val[2].a = 1.0;
	    val[1].v = op2 - linear(op1);
	}
    }
    else
    {
        // rotation autour de la droite D(O,VD) intersection des plans de 
        // l'angle des normales (VP1, VP2) la'angle = angle des droite
        // mais A := cotgt (alpha/2 ) d'ou le calcul

	double r = vp1.ps(vp2) / (vp1.norm() * vp2.norm());	// = cos alpha

	stl_v tmp_v = vd;
	tmp_v.normalize();
	double rep = tmp_v.orientation();
	if (rep < 0)
	    tmp_v = -tmp_v;
	double rac2 = fabs((r + 1.0) / (r - 1.0));
	double val_a = 0.0;
	if (rac2 >= sqr(eps))
	    val_a = normalise_angle(rep * sqrt(rac2));
	val[0].v = tmp_v;
	val[0].a = val_a;
	val[1].a = 1.0;
	val[2].a = 1.0;
	val[1].v = op2 - linear(op1);
    }
}
        
stl_transf::stl_transf( 
	const	stl_v& x1, 
	const	stl_v& x2, 
	const	stl_v& x3, 
  	const	stl_v& x4, 
	const	stl_v& x5, 
	const	stl_v& x6 
	) 
{
    /*-
     ! Definir la rotation qui amene un plan sur l'autre
     ! Puis Composer avec une rotation qui ramene X3 en X6
     ! Puis Faire correspondre X1 et X4
     */
    precondition(x2.normalized());
    precondition(x3.normalized());
    precondition(x5.normalized());
    precondition(x6.normalized());

    stl_transf rsim = stl_transf(x1, x2, x4, x5);

    stl_v x9 = rsim.linear(x3);

    double tpr =  (!x9.null())  ? x9.angle(x6, x5) : 0.0;
    stl_transf rsim2 = stl_transf(tpr, x5, stl_v_zero);

    stl_transf rsim3 = rsim2 * rsim;

    *this = stl_transf(x4 - (rsim3 * x1)) * rsim3;
}

stl_transf::stl_transf( 
	const	double alpha1, 
 	const	double alpha2, 
	const	double alpha3 
	) 
{
    stl_transf r1 = stl_r_id;
    stl_transf r2 = stl_r_id;
    stl_transf r3 = stl_r_id;

    if (!is_small(alpha1))
	r1 = stl_transf(alpha1, stl_v_k, stl_v_zero);
    if (!is_small(alpha2))
	r2 = stl_transf(alpha2, stl_v_j, stl_v_zero);
    if (!is_small(alpha3))
	r3 = stl_transf(alpha3, stl_v_i, stl_v_zero);
    *this = r1 * r2 * r3;
}

stl_quaternion stl_transf::quaternion( 
	const	int i 
	) const 
{
    precondition(i >= 0);
    precondition(i < 3);
    return (val[i]);
}

stl_transf stl_transf::prim_compose( 
 	const	stl_transf& q_sim2, 
		bool with_normalize 
	) const 
{
    /*-
     ! Compose deux transfo quelconques
     !
     ! Regle de Composition des parties Lineaires:
     ! ---------------------------------
     !
     !  !-----+-----+-----+-----+-----+
     !  !     !     !     !     !     !
     !  !     ! -2  ! -1  !  1  !  2  !
     !  !     !     !     !     !     !
     !  !-----+-----+-----+-----+-----+
     !  !     !     !     !     !     !
     !  ! -2  !  4  !  2  ! -2  !  x  !
     !  !     !     !     !     !     !
     !  !-----+-----+-----+-----+-----+
     !  !     !     !     !     !     !
     !  ! -1  !  2  !  1  ! -1  !  x  !
     !  !     !     !     !     !     !
     !  !-----+-----+-----+-----+-----+
     !  !     !     !     !     !     !
     !  !  1  ! -2  ! -1  !  1  !  x  !
     !  !     !     !     !     !     !
     !  !-----+-----+-----+-----+-----+
     !  !     !     !     !     !     !
     !  !  2  !  x  !  x  !  x  !  x  !
     !  !     !     !     !     !     !
     !  !-----+-----+-----+-----+-----+
     */
    stl_transf result;
    int det1 = (int) val[1].a;
    int det2 = (int) q_sim2.val[1].a;
    if (det1 == 2)
	result = q_sim2;
    else if (det2 == 2)
	result = *this;
    else
    {
	switch (det1 * det2)
	{
	case -2:		/* OK */
	    {
		/*-
                 !   1         *          -2
                 ! Rotation compose par sym/centrale
                 ! Donne une symetrie/rotation de meme
                 ! axis et d'angle ( PI + ancien ) qu'on
                 ! naturellement
                 */
		if (det1 == 1)
		    result.val[0] = val[0];
		else
		    result.val[0] = q_sim2.val[0];
		result.val[1].a = -1.0;
	    }
	    break;
	case -1:		/* NOK */
	    {
		/*-
                 !   1         *          -1
                 ! Rotation compose symetrie/plan
                 ! Donne une symetrie/rotation
                 ! Cas particulier : rota = demi tour
                 ! ou sym/rot et rot d'angle oppose
                 !   alors -identite
                 */
		if (val[0].product(q_sim2.val[0],
				   result.val[0], with_normalize))
		    result.val[1].a = -1.0;
		else
		    result.val[1].a = -2.0;
	    }
	    break;
	case 1:		/* OK */
	    {
		/*-
                 !   1         *          1
                 ! Rotation compose Rotation
                 ! Donne une Rotation
                 ! Cas particulier : Axe egaux et angle
                 ! opposes alors identite
                 !
                 !             OU
                 !   -1         *          -1
                 ! symetrie/plan compose symetrie/plan
                 ! Donne une rotation par rapport a
                 ! l'intersection des plans
                 ! Cas particulier : Axe egaux
                 !   alors -> identite
                 */
		if (val[0].product(q_sim2.val[0],
				   result.val[0], with_normalize))
		    result.val[1].a = 1.0;
		else
		    result.val[1].a = 2.0;
	    }
	    break;
	case 2:		/* OK */
	    {
		/*-
                 !   -1        *           -2
                 ! Sym/Rotation compose par sym/cent
                 ! Donne une rotation de meme
                 ! axis et d'angle ( PI + ancien ) qu'on
                 ! obtient naturellement
                 */
		if (det1 == -1)
		    result.val[0] = val[0];
		else
		    result.val[0] = q_sim2.val[0];
		result.val[1].a = 1.0;
	    }
	    break;
	case 4:		/* OK */
	    {
		/*-
                 !   -2        *          -2
                 ! symetrie/cent compose symetrie/cent
                 ! Donne identite
                 */
		result.val[0] = stl_q_zero;
		result.val[1].a = 2.0;
	    }
	    break;
	default:;
	    break;
	}
    }
    /*-
     ! Orientation standard du vecteur produit
     ! Puis Calcul de la partie AFFINE
     ! Puis Calcul de la partie HOMOTHETIE
     */
    result.val[0].orstd();
    result.val[1].v = (*this) * q_sim2.val[1].v;
    result.val[2].a = val[2].a * q_sim2.val[2].a;
    return (result);
}

bool stl_transf::operator == ( 
	const	stl_transf & q2 
	) const 
{
    if (val[2].a != q2.val[2].a)
	return (false);
    if (val[1].a != q2.val[1].a)
	return (false);
    if (!(val[1].v == q2.val[1].v))
	return (false);
	/*- Id et -id traite a par */
    if (fabs(val[1].a) != 2.0)
	return ((val[0].v == q2.val[0].v)
		&& is_small(val[0].a - q2.val[0].a));
    return (true);
}

stl_transf stl_transf::operator * ( 
	const	stl_transf &q 
	) const 
{
    stl_transf result = *this;
    result *= q;
    return (result);
}

stl_transf& stl_transf::operator *= ( 
	const	stl_transf &q 
	) 
{
    *this = prim_compose(q, true);
    return (*this);
}

stl_v stl_transf::operator * ( 
	const	stl_v &v1 
	) const 
{
    return (scale(linear(v1)) + val[1].v);
}

stl_v stl_transf::linear( 
	const	stl_v &v1 
	) const 
{
    stl_v result = v1;
    int det = (int) val[1].a;
    switch (det)
    {
    case -1:
    case 1:
	{
	    stl_v v_axis = val[0].v * (2.0 * v1.ps(val[0].v));
	    if (is_small(val[0].a))
		result = v_axis - v1 /* ROTATION/DEMI-TOUR */ ;
	    else
	    {
		stl_v vpara = v1 * (sqr(val[0].a) - 1.0);
		stl_v vperp = (val[0].v * v1) * (2.0 * val[0].a);
		result = (vpara + vperp + v_axis) *
		    (1.0 / (sqr(val[0].a) + 1.0));
	    }
	}
	break;
    case -2:
    case 2:
    default:;
	break;
    }
    if (det < 0)
	result = -result;
    return (result);
}
 
stl_v stl_transf::scale( const stl_v& v1) const
{
    if ((val[2].a != 1.0) && !is_small(val[2].a))
	return (v1 * val[2].a);
    return (v1);
}

stl_transf stl_transf::inverse( ) const
{
    /*-
     !  la transformation se calcule par :
     !           V2 = k Q[ 1] * V1 + Q[ 2].V
     ! d'ou l'inverse par
     !         -1      -1         -1      -1
     !   V1 = k   Q[ 1]   * V1 - k   Q[ 1]  * Q[ 2].V
     ! On impose que l'image de Q[ 2].V devienne l'origine
     */
    stl_transf r2 = *this;
    r2.val[0].a = -val[0].a;
    if (!is_small(val[2].a))
	r2.val[2].a = 1.0 / val[2].a;
    r2.val[1].v = -r2.scale(r2.linear(val[1].v));
    return (r2);
}
                                                    
bool stl_transf::is_translation() const
{
    return (val[0].v.null() &&
	    is_small(val[0].a) &&
	    is_small(val[1].a - 2.0) &&
	    is_small(val[2].a - 1.0));
}

bool stl_transf::is_rotation() const
{
    return (is_small(val[1].a - 1.0));
}

bool stl_transf::is_plane_reflection() const
{
    return (is_small(val[1].a + 1.0));
}

bool stl_transf::is_point_symetry() const
{
    return (is_small(val[1].a + 2.0));
}

bool stl_transf::is_scaling() const
{
    return (!is_small(val[2].a - 1.0));
}

i16 stl_transf::sign() const
{
    if (val[1].a < 0)
	return (-1);
    return (1);
}

double stl_transf::scale_factor() const
{
    return (val[2].a);
}

stl_v stl_transf::translate_vector() const
{
    return (val[1].v);
}

double stl_transf::angle() const
{
    precondition(is_rotation());
    /*-
     ! Dans un T_R, l'angle est stoque sous la forme: COTG( Alpha / 2 )
     */
    if (is_small(val[0].a))
	return (pi);
    return (2.0 * arctan(1.0 / val[0].a));
}

static double v_vangle_pi( 
	const	stl_v &v1, 
	const	stl_v &v2, 
	const	stl_v &v3 
	) 
{
    double alpha = v1.angle(v2, v3);
    if (pi < alpha)
	alpha -= r_2pi;
    return (alpha);
}

void stl_transf::angles_get( 
		double angles[ 3] 
	) const 
{
    /*-
     ! Method:
     !  0/ 	The plane of the second rotation is perpendicular to the plane
     !		of the first one ( plane OXY if we work with canonical
     !		frame ( v_i,v_jv_k )).
     !  1/	The normalized projection of  TV_I in the plane  XY is Iaz(V_I)
     !		The angle with  V_I is The angle of the first rotation Az
     !		rem: The plane (Iaz(V_I),TV_I) is the plane of the
     !		second rotation.
     !
     !  2/	The angle of  TV_I with Iaz(V_I) is The angle of the second
     !		rotation Ay
     !
     !  3/	The vector Iaz(V_J) transform of V_J by the first rotation is
     !		invariant by the  second rotation (  axis ).
     !		his angle with  TV_J is the angle of the third rotation Ax
     !
     !  4/	control with  V_K and TV_K
     !  note:	we work with function  V_VANGLE which  normalize vectors
     */
    stl_v tv_i = linear(stl_v_i);
    stl_v tv_j = linear(stl_v_j);
    stl_v v2 = tv_i.perp(stl_v_k);
    if (!v2.null())
    {
	v2.normalize();
	angles[0] = v_vangle_pi(stl_v_i, v2, stl_v_k);
	stl_v v3 = stl_v_i * cos(angles[0] + (pi / 2.0)) +
	stl_v_j * sin(angles[0] + (pi / 2.0));
	angles[1] = v_vangle_pi(v2, tv_i, v3);
	angles[2] = v_vangle_pi(v3, tv_j, tv_i);
    }
    else
    {
	/*-
         !  TV_I est / / a V_K
	 ! angles[0] peut prendre n'importe quelle
	 !   valeur(= 0)
	 */
	angles[0] = 0.0;
	angles[1] = (0.0 - signe(stl_v_k.ps(tv_i))) * pi / 2.0;
	angles[2] = v_vangle_pi(stl_v_j, tv_j, tv_i);
    }
    /*-
     ! A TITRE DE DEBUG :
     ! on regarde si l'image de V_K par
     ! les rot. Ay et Ax est bien TV_K
     */

}

stl_v stl_transf::axis() const
{
    precondition(is_rotation() || is_plane_reflection());
    return (val[0].v);
}

stl_v stl_transf::center() const
{
    precondition(is_scaling() || is_point_symetry()
		 || is_rotation() || is_plane_reflection());
    stl_v result = val[1].v;
    double sc = 0.5;
    if (is_scaling())
	sc = -1.0 / (val[2].a - 1.0);
    if (is_rotation() || is_plane_reflection())
    {
    }
    else
    {
	stl_lin l_r = stl_l_id - stl_lin(*this).linear();
	precondition(l_r.invertible());
	result = l_r.inverse() * result;
    }

    return (result * sc);
}

ostream& operator<<( 
		ostream& s, 
	const	stl_transf& r1 
	) 
{
    return (s << '[' << r1.val[0] << ','
	    << r1.val[1] << ',' << r1.val[2] << ']');
}


stl_lin::stl_lin()
{
    for (int i = 0; i < 4; i++)
	l[i] = stl_v_zero;
}

stl_lin::stl_lin( 
	const	stl_v &v1, 
	const	stl_v &v2, 
	const	stl_v &v3, 
	const	stl_v &v4 
	) 
{
    l[0] = v1;
    l[1] = v2;
    l[2] = v3;
    l[3] = v4;
}

stl_lin::stl_lin( const stl_v& v1)
{
    *this = stl_l_id;
    l[3] = v1;
}

stl_lin::stl_lin( 
	const	double r 
	) 
{
    *this = stl_l_id * r;
}

stl_lin::stl_lin( 
	const	double alpha, 
	const	stl_v &axis, 
	const	stl_v &origin 
	) 
{
    precondition(axis.normalized());
    *this = stl_transf(alpha, axis, origin);
}

stl_lin::stl_lin( 
	const	double r_x, 
	const	double r_y, 
	const	double r_z 
	) 
{
    *this = stl_l_id;
    l[0] *= r_x;
    l[1] *= r_y;
    l[2] *= r_z;
}

stl_lin::stl_lin( 
	const	stl_transf& r 
	) 
{
    l[0] = r.linear(stl_v_i);
    l[1] = r.linear(stl_v_j);
    l[2] = r.linear(stl_v_k);
    l[3] = r * stl_v_zero;
}

stl_lin::stl_lin( 
	const	stl_lin& l_from, 
	const	stl_lin& l_to 
	) 
{
    precondition(l_from.invertible());
    *this = l_to * l_from.inverse();
}

SELF_CONSTRUCTOR(stl_lin)
void stl_lin::column_set( 
	const	int i,const 
		stl_v& v1 
	) 
{
    precondition(i >= 0);
    precondition(i < 4);
    l[i] = v1;
}
void stl_lin::translation_set( 
	const	stl_v& v1 
	) 
{
    l[3] = v1;
}
stl_v stl_lin::column( 
	const	int i 
	) const 
{
    precondition(i >= 0);
    precondition(i < 4);
    return (l[i]);
}
stl_v stl_lin::translation( ) const
{
    return (l[3]);
}
void stl_lin::x_axis_set( 
	const	stl_v &v1 
	) 
{
    l[0] = v1;
}
void stl_lin::y_axis_set( 
	const	stl_v &v1 
	) 
{
    l[1] = v1;
}
void stl_lin::z_axis_set( 
	const	stl_v &v1 
	) 
{
    l[2] = v1;
}
void stl_lin::origin_set( 
	const	stl_v &v1 
	) 
{
    l[3] = v1;
}
stl_v stl_lin::x_axis() const
{
    return (l[0]);
}
stl_v stl_lin::y_axis() const
{
    return (l[1]);
}
stl_v stl_lin::z_axis() const
{
    return (l[2]);
}
stl_v stl_lin::origin() const
{
    return (l[3]);
}

bool stl_lin::operator == ( 
	const	stl_lin & l2 
	) const 
{
    for (int i = 0; i < 4; i++)
	if (!(l[i] == l2.l[i]))
	    return (false);
    return (true);
}

stl_lin stl_lin::operator + ( 
	const	stl_lin & l2 
	) const 
{
    stl_lin result = *this;
    result += l2;
    return (result);
}

stl_lin stl_lin::operator - ( 
	const	stl_lin & l2 
	) const 
{
    stl_lin result = *this;
    result -= l2;
    return (result);
}

stl_lin stl_lin::operator *( 
	const	double r 
	) const 
{
    stl_lin result = *this;
    result *= r;
    return (result);
}

stl_v stl_lin::operator *( 
	const	stl_v &v1 
	) const 
{
    stl_v result = l[3];
    for (int i = 0; i < 3; i++)
	for (int j = 0; j < 3; j++)
	    result.x[i] += l[j].x[i] * v1.x[j];
	/*
	for (int i = 0; i < 3; i++)
	{
		cout << char('x' + i) << "' = l3[" << i << "]";
		for (int j = 0; j < 3; j++)
		{
			result.x[i] += l[j].x[i] * v1.x[j];
			cout << " + l" << j << "[" << i << "] * " << char('x' + i);
		}
		cout << endl;
	}
	*/
    return (result);
}
/*
x1 = l0.x1 
*/

stl_lin stl_lin::operator *( 
	const	stl_lin & l2 
	) const 
{
    stl_lin result;
    result.l[3] = l[3];
    for (int i = 0; i < 3; i++)	// row. 
	for (int j = 0; j < 4; j++)	// collumn.
	    for (int k = 0; k < 3; k++)
		result.l[j].x[i] += l[k].x[i] * l2.l[j].x[k];
    return (result);
}

stl_lin & stl_lin::operator += ( 
	const	stl_lin & l2 
	) 
{
    for (int i = 0; i < 4; i++)
	l[i] += l2.l[i];
    return (*this);
}

stl_lin& stl_lin::operator -= ( 
	const	stl_lin& l2 
	) 
{
    for (int i = 0; i < 4; i++)
	l[i] -= l2.l[i];
    return (*this);
}

stl_lin& stl_lin::operator *= ( 
	const	double r 
	) 
{
    for (int i = 0; i < 4; i++)
	l[i] *= r;
    return (*this);
}
 
stl_lin& stl_lin::operator *= ( 
	const	stl_lin& l2 
	)  
{
    *this = *this * l2;
    return (*this);
}

double stl_lin::det() const
{
    return (l[0].mixte(l[1], l[2]));
}
  
double stl_lin::trace() const
{
    return (l[0].x[0] + l[1].x[1] + l[2].x[2]);
}

bool stl_lin::is_orthogonal_frame() const
{
    /*-
     ! the 3 products are done for symetry and stability
     ! does not work with det() for stability.
     */
    for (int i = 0; i < 3; i++)
    {
	if (!l[i].normalized())
	    return (false);
	if (!l[(i + 2) % 3].paral(l[i] * l[(i + 1) % 3]))
	    return (false);
    }
    return (true);
}

bool stl_lin::is_direct_orthogonal_frame() const
{
    /*-
     ! the 3 products are done for symetry and stability
     ! does not work with det() for stability.
     */
    for (int i = 0; i < 3; i++)
    {
	if (!l[i].normalized())
	    return (false);
	if (!(l[(i + 2) % 3] == l[i] * l[(i + 1) % 3]))
	    return (false);
    }
    return (true);
}
bool stl_lin::invertible() const
{
    return (!is_small(det()));
}

stl_lin stl_lin::inverse() const
{
    precondition(invertible());

    stl_lin result;
    double det1 = det();

    for (int i = 0; i < 3; i++)
    {
	i16 i1 = (i + 1) % 3;
	i16 i2 = (i + 2) % 3;
	for (int j = 0; j < 3; j++)
	{
	    int j1 = (j + 1) % 3;
	    int j2 = (j + 2) % 3;
	    result.l[j].x[i] = (l[i1].x[j1] * l[i2].x[j2]
				- l[i1].x[j2] * l[i2].x[j1]) / det1;
	}
    }
    result.l[3] = result.linear(-l[3]);

    return (result);
}

stl_lin stl_lin::transpose() const
{
    stl_lin result =  *this;
    for (int i = 0; i < 3; i++)
	for (int j = 0; j < i; j++)
	    exchange(result.l[i].x[j] , result.l[j].x[i]);
    return (result);
}
     
stl_lin stl_lin::change_basis( 
	const	stl_lin& change 
	) const 
{
    if (change.is_direct_orthogonal_frame())
	return (change.transpose() * (*this) * change);
    return (change.inverse() * (*this) * change);
}
                               
bool stl_lin::pivot( 
		stl_v &sol 
	)  
{
    //	on permutte la colonne du plus grand nombre, qui  est le pivot
    int i;
    int j;
    int k;
    int n[4];
    double piv;

    for (j = 0; j < 4; j++)
    {
	n[j] = j;
    }

    for ( i = 0; i < 3; i++)
    {
	piv = 0.0;
	for (j = 0; j < 4; j++)
	{
	    if (fabs(l[n[j]].x[i]) > fabs(piv))
	    {
		piv = l[n[j]].x[i];
		int k_tmp = n[j];
		n[j] = n[i];
		n[i] = k_tmp;
	    }
	}

	if (fabs(piv) < eps)
	{
	    return (false);
	}

	for (k = 0; k < 4; k++)
	{
	    l[n[k]].x[i] = (k != i) ? (l[n[k]].x[i] / piv) : 1.0;
	}

	for (j = 0; j < 3; j++)
	{
	    if (j != i)
	    {
		double coef = l[n[i]].x[j];
		for (k = 0; k < 4; k++)
		{
		    l[n[k]].x[j] = (k != i) 
			? (l[n[k]].x[j] - coef * l[n[k]].x[i])
			: 0.0;
		}
	    }
	}
    }

    for (k = 0; k < 3; k++)
    {
	sol.x[n[k]] = l[3].x[k];
    }

    return (true);
}

double stl_lin::quadratique( 
	const	stl_v& w1, 
	const	stl_v& w2 
	) const 
{
    return (w1.ps(linear(w2)));
}

stl_v stl_lin::linear( 
	const	stl_v &v 
	) const 
{
    stl_v result = stl_v_zero;
    for (int i = 0; i < 3; i++)
	for (int j = 0; j < 3; j++)
	    result.x[i] += l[j].x[i] * v.x[j];
    return (result);
}

stl_lin stl_lin::linear() const
{
    return (stl_lin(l[0], l[1], l[2], stl_v_zero));
}

long stl_lin::to_q( 
		stl_transf & tr 
	) const 
{
    precondition(invertible());

    /*-
     ! Warning :
     ! 1- We don't work directly with det() but only with
     !	sign of det() for stability.
     ! 2- we don`t buid tr as a symetry or rotation by call of t_r constructors
     !	  because we work with the explicite translative part and not
     !	  an origin ( of rotation or symetry).
     ! 3- a mixed transformation scaling t_l could be theorically translated
     !	 in a correct t_r we don`t do that today because we forget
     !	 the scaling part so "!!"!.
     !
     ! On evalue le determinant et le cosinus  qui se deduit
     ! de la trace. Ce sont deux invariants de la matrice
     ! Le determinant doit etre nettement # de  0 bien que ce
     ! ne soit pas une condition  suffisante pour etre un
     ! deplacement.
     ! Pour le calcul du cosinus : Voir Les Quaternions de
     ! Casteljau  page 57 formule 4.18
     ! Trace = 1 + 2 * cos_phi
     ! Il peut s'agir de :
     ! ou d'une rotation autour d'un axis
     !   ou d'une rotation/symetrie perp
     ! ou d'un  retournement autour d'un axis
     !   ou d'un reflexion a travers un pl
     ! ou de l'identite ou de  -identite
     */
    stl_lin rot = *this;

    if (!is_orthogonal_frame())
    {
	double scale = (l[0].norm() + l[1].norm() + l[2].norm()) / 3.0;

	if (is_small(scale))
	    return (0);

	rot *= 1.0 / scale;
	if (!rot.is_orthogonal_frame())
	    return (0);
    }

    double trace1 = rot.trace();
    double det1 = signe(rot.det());
    double cos_phi = (trace1 - det1) * det1 / 2.0;

    if (cos_phi < (eps - 1.0))
	cos_phi = -1.0;
    else if (cos_phi > (1.0 - eps))
	cos_phi = 1.0;

    double sin_phi = sqrt(fabs(1.0 - sqr(cos_phi)));

    tr = stl_transf(rot.l[3]);

    if (!is_small(sin_phi))
    {
	/*-
	! C'est une rotation d'angle # 180
	! Voir Les Quaternions de Casteljau page 57 formule 4.19
	! Consequence de T*transposee(T) = I
	! justification fausse
	! Si symetrie AXE=-AXE \ au cas ou
	! et COS_PHI=-COS_PHI  |  on aurait
	! /tout multiplie  par -1
	! Mais comme TR[1].A = -1/TR[1].A o.k
        !
        ! axis cannot be null because we are in a stable
        ! case : rot.is_orthogonal_frame() then rot.det() is close
        ! to 1.0 and 0< sin_phi <=1
	*/
	stl_v axis = stl_v_zero;
	axis.x[0] = (rot.l[2].x[1] - rot.l[1].x[2]) / sin_phi;
	axis.x[1] = (rot.l[0].x[2] - rot.l[2].x[0]) / sin_phi;
	axis.x[2] = (rot.l[1].x[0] - rot.l[0].x[1]) / sin_phi;

	invariant(!axis.null());
	axis.normalize();
	if (axis.orientation() < 0)
	{
	    sin_phi = -sin_phi;
	    axis = -axis;
	}
	double val_a = sqrt((1.0 + cos_phi) / (1.0 - cos_phi));
	/*-
         ! Axis and sinus are known but there sign are not
         ! then we look for sign and we know what is TR[1].A
         */
	double signe = 0.0;
	for (int i = 0; i < 3; i++)
	{
	    stl_v test = rot.l[i] * axis;
	    signe += test.x[i];
	}
	if (signe * det1 < 0.0)
	    tr.val[0].real_part_set(-val_a);
	else
	    tr.val[0].real_part_set(val_a);
	tr.val[0].vector_part_set(axis);
	tr.val[1].real_part_set(det1);
    }
    else if (cos_phi < 0.0)
    {
	/*-
	 !  f(v) + v is colinear to the axis:  rotation of 180 degrees (det1>0)
	 !  f(v) - v  is colinear to the axis: plane-reflexion (det1<0)
	 !
	 !  Then the axis "v" is solution of the linear equation:
	 !   	( f - ( det1 * Id ) ) v = 0
	 */

	for (int i = 0; i < 3; i++)
	    rot.l[i].x[i] += det1;

	for (int i = 0; i < 3; i++)
	    if (rot.l[i].orientation() < 0)
		rot.l[i] = - rot.l[i];

	stl_v axis = rot.l[0] + ( rot.l[1] * 2 ) + rot.l[2];
	invariant(!axis.null());
	axis.normalize();

	tr.val[0].real_part_set(0.0);
	tr.val[0].vector_part_set(axis);
	tr.val[1].real_part_set(det1);
    }
    else
    {
	/*-
	 ! Case of identity/translation
	 ! or point symetry of  center O
	*/
	if (det1 < 0)
	    tr.val[1].real_part_set(-2.0);
    }

    return (1);
}

double small0(double v);

ostream& operator<<( 
		ostream& s, 
	const	stl_lin& l1 
	) 
{
	char buf[80];
    s << l1.l[3];
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
		{
			sprintf(buf," %.3g", small0(l1.l[j].x[i]));
			s << buf;
		}
	return (s);
}


const stl_quaternion stl_q_zero = stl_quaternion(0.0, 0.0, 0.0, 0.0 );
const stl_transf stl_r_id  = stl_transf(	stl_quaternion( 0.0, 0.0, 0.0, 0.0 ),
                        		stl_quaternion( 2.0, 0.0, 0.0, 0.0 ),
                        		stl_quaternion( 1.0, 1.0, 1.0, 1.0 ));

const stl_lin stl_l_zero 	= stl_lin(	stl_v(0.0, 0.0, 0.0),
					stl_v(0.0, 0.0, 0.0),
					stl_v(0.0, 0.0, 0.0),
					stl_v(0.0, 0.0, 0.0));
const stl_lin stl_l_id	= stl_lin(	stl_v(1.0, 0.0, 0.0),
					stl_v(0.0, 1.0, 0.0),
					stl_v(0.0, 0.0, 1.0),
					stl_v(0.0, 0.0, 0.0));

