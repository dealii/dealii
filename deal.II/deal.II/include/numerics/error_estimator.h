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
   Implementation of the error estimator by Kelly, Gago, Zienkiewicz and
   Babuska.
   This error estimator tries to approximate the error per cell by integration
   of the jump of the gradient of the solution along the faces of each cell.
   It can be understood as a gradient recovery estimator; see the survey
   of Ainsworth for a complete discussion.

   It seem as if this error estimator should only be valid for linear ansatz
   spaces, but no definite answer is given to this question at present.
   
   
   {\bf Implementation}

   In principle, the implementation of the error estimation is simple: let
   $$ \eta_K^2 =
   h \int_{\partial K} \left[\frac{\partial u_h}{\partial n}\right]^2 do
   $$
   be the error estimator for cell $K$. $[\cdot]$ denotes the jump of the
   argument at the face. In the paper of Ainsworth, $h$ is divided by $24$,
   but this factor is a bit esoteric, stemming from interpolation estimates
   and stability constants which may hold for the Poisson problem, but may
   not hold for more general situations. In the implementation, this factor
   is dropped for these reasons.

   To perform the integration, use is made of the #FEFaceValues# class and the
   integration is performed for each cell, i.e. no use is made of the fact, that
   the integration along a face need in principle be done only once for both
   adjacent cells. Clearly there is room for optimization here.

   If the face is at the boundary, i.e. there is no neighboring cell to which
   the jump in the gradiend could be computed, there are two possibilities:
   \begin{itemize}
   \item The face belongs to a Dirichlet boundary. Then the face is not
     considered, which can be justified looking at a dual problem technique and
     should hold exactly if the boundary can be approximated exactly by the
     finite element used (i.e. it is a linear boundary for linear finite elements,
     quadratic for isoparametric quadratic elements, etc). For boundaries which
     can not be exactly approximated, one should consider the difference
     $z-z_h$ on the face, $z$ being a dual problem's solution which is zero at
     the true boundary and $z_h$ being an approximation, which in most cases
     will be zero on the numerical boundary. Since on the numerical boundary
     $z$ will not be zero in general, we would get another term here, but this
     one is neglected for practical reasons, in the hope that the error made
     here will tend to zero faster than the energy error we wish to estimate.

   \item The face belongs to a Neumann boundary.  In this case, the
     contribution of the face $F\in\partial K$ looks like
     $$ \int_F \left|g-\frac{\partial u_h}{\partial n}\right| ds $$
     where $g$ is the Neumann boundary function.

   \item No other boundary conditions are considered.
   \end{itemize}
   
   @author Wolfgang Bangerth, 1998; thanks to Franz-Theo Suttmeier for
     clarifications about boundary conditions.
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
    typedef map<unsigned char,Function<dim>*> FunctionMap;
    
    void estimate_error (const DoFHandler<dim>    &dof,
			 const Quadrature<dim-1>  &quadrature,
			 const FiniteElement<dim> &fe,
			 const Boundary<dim>      &boundary,
			 const FunctionMap        &neumann_bc,
			 const dVector            &solution,
			 dVector                  &error) const;

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
};



/*----------------------------   error_estimator.h     ---------------------------*/
/* end of #ifndef __error_estimator_H */
#endif
/*----------------------------   error_estimator.h     ---------------------------*/
