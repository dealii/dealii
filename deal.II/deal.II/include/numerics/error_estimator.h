/*----------------------------   error_estimator.h     ---------------------------*/
/*      $Id$                 */
#ifndef __error_estimator_H
#define __error_estimator_H
/*----------------------------   error_estimator.h     ---------------------------*/


#include <base/exceptions.h>
#include <basic/function.h>
#include <map>


// forward declarations
template <int dim> class DoFHandler;
template <int dim> class Quadrature;
template <int dim> class FiniteElement;
template <int dim> class Boundary;
template <int dim> class Function;
class dVector;




/**
 *  Implementation of the error estimator by Kelly, Gago, Zienkiewicz and
 *  Babuska.
 *  This error estimator tries to approximate the error per cell by integration
 *  of the jump of the gradient of the solution along the faces of each cell.
 *  It can be understood as a gradient recovery estimator; see the survey
 *  of Ainsworth for a complete discussion.
 *
 *  It seem as if this error estimator should only be valid for linear ansatz
 *  spaces, and there are indications that for higher order ansatz spaces the
 *  integrals computed here show superconvergence properties, i.e. they tend
 *  to zero faster than the error itself, thus ruling out the values as error
 *  indicators.
 *  
 *  The error estimator returns a vector of estimated errors per cell which
 *  can be used to feed the #Triangulation<dim>::refine_*# functions.
 *
 *  
 *  \subsection{Implementation}
 *
 *  In principle, the implementation of the error estimation is simple: let
 *  $$ \eta_K^2 =
 *  \frac h{24} \int_{\partial K} \left[\frac{\partial u_h}{\partial n}\right]^2 do
 *  $$
 *  be the error estimator for cell $K$. $[\cdot]$ denotes the jump of the
 *  argument at the face. In the paper of Ainsworth, $h$ is divided by $24$,
 *  but this factor is a bit esoteric, stemming from interpolation estimates
 *  and stability constants which may hold for the Poisson problem, but may
 *  not hold for more general situations. In the implementation, this factor
 *  is considered, but may lead to wrong results. You may scale the vector
 *  appropriately afterwards.
 *
 *  To perform the integration, use is made of the #FEFaceValues# and
 *  #FESubfaceValues# classes. The integration is performed by looping
 *  over all cells and integrating over faces that are not yet treated.
 *  This way we avoid integration on faces twice, once for each time we
 *  visit one of the adjacent cells. In a second loop over all cells, we
 *  sum up the contributions of the faces (which are the integrated
 *  square of the jumps) of each cell and take the square root.
 *
 *  The integration is done using a quadrature formula on the face.
 *  For linear ansatz functions (#FELinear#), the #Gauss2# or even the
 *  #Midpoint# rule will suffice. For higher order elements, it is
 *  necessary to utilize higher order quadrature formulae as well.
 *
 *  We store the contribution of each face in a #map#, as provided by the
 *  C++ standard library, with the iterator pointing to that face being the
 *  key into the map. In fact, we do not store the indicator per face, but
 *  only the integral listed above. When looping the second time over all
 *  cells, we have to sum up the contributions of the faces, multiply them
 *  with $\frac h{24}$ and take the square root. By doing the multiplication
 *  with $h$ in the second loop, we avoid problems to decide with which $h$
 *  to multiply, that of the cell on the one or that of the cell on the other
 *  side of the face.
 *
 *  $h$ is taken to be the greatest length of the diagonals of the cell. For
 *  more or less uniform cells without deformed angles, this coincides with
 *  the diameter of the cell.
 *  
 *
 *  \subsection{Boundary values}
 *  
 *  If the face is at the boundary, i.e. there is no neighboring cell to which
 *  the jump in the gradiend could be computed, there are two possibilities:
 *  \begin{itemize}
 *  \item The face belongs to a Dirichlet boundary. Then the face is not
 *    considered, which can be justified looking at a dual problem technique and
 *    should hold exactly if the boundary can be approximated exactly by the
 *    finite element used (i.e. it is a linear boundary for linear finite elements,
 *    quadratic for isoparametric quadratic elements, etc). For boundaries which
 *    can not be exactly approximated, one should consider the difference
 *    $z-z_h$ on the face, $z$ being a dual problem's solution which is zero at
 *    the true boundary and $z_h$ being an approximation, which in most cases
 *    will be zero on the numerical boundary. Since on the numerical boundary
 *    $z$ will not be zero in general, we would get another term here, but this
 *    one is neglected for practical reasons, in the hope that the error made
 *    here will tend to zero faster than the energy error we wish to estimate.
 *
 *    Though no integration is necessary, in the list of face contributions we
 *    store a zero for this face, which makes summing up the contributions of
 *    the different faces to the cells easier.
 *
 *  \item The face belongs to a Neumann boundary.  In this case, the
 *    contribution of the face $F\in\partial K$ looks like
 *    $$ \int_F \left|g-\frac{\partial u_h}{\partial n}\right| ds $$
 *    where $g$ is the Neumann boundary function.
 *
 *  \item No other boundary conditions are considered.
 *  \end{itemize}
 *
 *  Thanks go to Franz-Theo Suttmeier for clarifications about boundary
 *  conditions.
 *
 *  
 *  \subsection{Handling of hanging nodes}
 *  
 *  The integration along faces with hanging nodes is quite tricky, since one
 *  of the elements has to be shifted one level up or down. See the
 *  documentation for the #FESubfaceValues# class for more information about
 *  technical issues regarding this topic.
 *
 *  In praxi, since we integrate over each face only once, we do this when we
 *  are on the coarser one of the two cells adjacent to a subface (a subface
 *  is defined to be the child of a face; seen from the coarse cell, it is a
 *  subface, while seen from the refined cell it is one of its faces). The
 *  reason is that finding neighborship information is a bit easier then, but
 *  that's all practical reasoning, nothing fundamental.
 *
 *  Since we integrate from the coarse side of the face, we have the mother
 *  face readily at hand and store the result of the integration over that
 *  mother face (being the sum of the integrals along the subfaces) in the
 *  abovementioned map of integrals as well. This consumes some memory more
 *  than needed, but makes the summing up of the face contributions to the
 *  cells easier, since then we have the information from all faces of all
 *  cells at hand and need not think about explicitely determining whether
 *  a face was refined or not. The same applies for boundary faces, see
 *  above.
 *  
 *  @author Wolfgang Bangerth, 1998
 */
template <int dim>
class KellyErrorEstimator {
  public:
    				     /**
				      *	Declare a data type which denotes a
				      *	mapping between a boundary indicator
				      *	and the function denoting the boundary
				      *	values on this part of the boundary.
				      *	Only one boundary function may be given
				      *	for each boundary indicator, which is
				      *	guaranteed by the #map# data type.
				      */
    typedef map<unsigned char,const Function<dim>*> FunctionMap;
    
    static void estimate_error (const DoFHandler<dim>    &dof,
				const Quadrature<dim-1>  &quadrature,
				const FiniteElement<dim> &fe,
				const Boundary<dim>      &boundary,
				const FunctionMap        &neumann_bc,
				const dVector            &solution,
				dVector                  &error);

				     /**
				      * Exception
				      */
    DeclException0 (ExcNotImplemented);
				     /**
				      * Exception
				      */
    DeclException0 (ExcInternalError);
				     /**
				      * Exception
				      */
    DeclException0 (ExcInvalidBoundaryIndicator);

  private:
				     /**
				      * Declare a data type to represent the
				      * mapping between faces and integrated
				      * jumps of gradients. See the general
				      * doc of this class for more information.
				      */
    typedef map<DoFHandler<dim>::face_iterator,double> FaceIntegrals;

				     /**
				      * Redeclare an active cell iterator.
				      * This is simply for convenience.
				      */
    typedef DoFHandler<dim>::active_cell_iterator active_cell_iterator;
    
				     /**
				      * Actually do the computation on a face
				      * which has no hanging nodes (it is
				      * regular),  i.e.
				      * either on the other side there is
				      * nirvana (face is at boundary), or
				      * the other side's refinement level
				      * is the same as that of this side,
				      * then handle the integration of
				      * these both cases together.
				      *
				      * The meaning of the parameters becomes
				      * clear when looking at the source
				      * code. This function is only
				      * externalized from #estimate_error#
				      * to avoid ending up with a function
				      * of 500 lines of code.
				      */
    static void integrate_over_regular_face (const active_cell_iterator &cell,
					     const unsigned int   face_no,
					     const FiniteElement<dim> &fe,
					     const Boundary<dim>      &boundary,
					     const FunctionMap   &neumann_bc,
					     const unsigned int   n_q_points,
					     FEFaceValues<dim>   &fe_face_values_cell,
					     FEFaceValues<dim>   &fe_face_values_neighbor,
					     FaceIntegrals       &face_integrals,
					     const dVector       &solution);

				     /**
				      * The same applies as for the function
				      * above, except that integration is
				      * over face #face_no# of #cell#, where
				      * the respective neighbor is refined,
				      * so that the integration is a bit more
				      * complex.
				      */
    static void integrate_over_irregular_face (const active_cell_iterator &cell,
					       const unsigned int   face_no,
					       const FiniteElement<dim> &fe,
					       const Boundary<dim>      &boundary,
					       const unsigned int    n_q_points,
					       FEFaceValues<dim>    &fe_face_values,
					       FESubfaceValues<dim> &fe_subface_values,
					       FaceIntegrals        &face_integrals,
					       const dVector        &solution);
};



/*----------------------------   error_estimator.h     ---------------------------*/
/* end of #ifndef __error_estimator_H */
#endif
/*----------------------------   error_estimator.h     ---------------------------*/
