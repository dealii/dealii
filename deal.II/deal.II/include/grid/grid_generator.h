/*----------------------------   grid_generator.h     ---------------------------*/
/*      $Id$                 */
#ifndef __grid_generator_H
#define __grid_generator_H
/*----------------------------   grid_generator.h     ---------------------------*/


#include <base/forward-declarations.h>
#include <basic/forward-declarations.h>

/**
 * This class offers triangulations of some standard domains such as hypercubes,
 * hyperball and the like. Following is a list of domains that can be generated
 * by the functions of this class:
 * \begin{itemize}
 *    \item Hypercube triangulations: a hypercube triangulation is a
 *       domain which is the tensor product of an interval $[a,b]$ in
 *       the given number of spatial dimensions. If you want to create such
 *       a domain, which is a common test case for model problems, call
 *       #GridGenerator::hyper_cube (tria, a,b)#, which produces a
 *       hypercube domain triangulated with exactly one element. You can
 *       get tensor product meshes by successive refinement of this cell.
 *
 *    \item Generalized L-shape domain:
 *      using the #GridGenerator::hyper_L (tria, a,b)# function produces
 *      the hypercube with the interval $[a,b]$ without the hypercube
 *      made out of the interval $[(a+b)/2,b]$. Let, for example, be $a=-1$
 *      and $b=1$, then the hpyer-L in two dimensions is the region
 *      $[-1,1]^2 - [0,1]^2$. To create a hyper-L in one dimension results in
 *      an error. The function is also implemented for three space dimensions.
 *
 *    \item Hyper balls:
 *      You get the circle or ball (or generalized: hyperball) around origin
 *      #p# and with radius #r# by calling
 *      #Triangulation<dim>::hyper_ball (p, r)#. The circle is triangulated
 *      by five cells, the ball by seven cells. The diameter of the center cell is
 *      chosen so that the aspect ratio of the boundary cells after one refinement
 *      is minimized in some way. To create a hyperball in one dimension results in
 *      an error.
 *
 *      Do not forget to attach a suitable boundary approximation object
 *      to the triangulation object you passed to this function if you later want
 *      the triangulation to be refined at the outer boundaries.
 * \end{itemize}
 *
 * @author Wolfgang Bangerth, 1998, 1999
 */
class GridGenerator 
{
  public:

    				     /**
				      * Initialize the given triangulation with a
				      * hypercube (line in 1D, square in 2D, etc)
				      * consisting of exactly one cell. The
				      * hypercube volume is the tensor product
				      * of the intervall $[left,right]$ in the
				      * present number of dimensions, where
				      * the limits are given as arguments. They
				      * default to zero and unity, then producing
				      * the unit hypercube.
				      *
				      * The triangulation needs to be void
				      * upon calling this function.
				      */
    template <int dim>
    static void hyper_cube (Triangulation<dim> &tria,
			    const double        left = 0.,
			    const double        right= 1.);

				     /**
				      * Initialize the given triangulation with a
				      * hyperball, i.e. a circle or a ball.
				      * See the general documentation for a
				      * more concise description. The center of
				      * the hyperball default to the origin,
				      * the radius defaults to unity.
				      *
				      * The triangulation needs to be void
				      * upon calling this function.
				      */    
    template <int dim>
    static void hyper_ball (Triangulation<dim> &tria,
			    const Point<dim>   &center = Point<dim>(),
			    const double        radius = 1.);

				     /**
				      * Initialize the given triangulation with a
				      * hyper-L consisting of exactly #2^dim-1#
				      * cells. See the general documentation for a
				      * description of the L-region. The limits
				      * default to minus unity and unity.
				      *
				      * The triangulation needs to be void
				      * upon calling this function.
				      */
    template <int dim>
    static void hyper_L (Triangulation<dim> &tria,
			 const double        left = -1.,
			 const double        right= 1.);
};




/*----------------------------   grid_generator.h     ---------------------------*/
/* end of #ifndef __grid_generator_H */
#endif
/*----------------------------   grid_generator.h     ---------------------------*/
