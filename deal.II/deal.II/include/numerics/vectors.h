/*----------------------------   vectors.h     ---------------------------*/
/*      $Id$                 */
#ifndef __vectors_H
#define __vectors_H
/*----------------------------   vectors.h     ---------------------------*/


template <int dim> class DoFHandler;
template <int dim> class Function;
template <int dim> class Quadrature;
template <int dim> class Boundary;
template <int dim> class FiniteElement;
class ConstraintMatrix;
class dVector;



/**
 * Provide a class which assembles some standard vectors. Among these are
 * interpolations and projections of continuous functions to the finite
 * element space and other operations.
 *
 *
 * \subsection{Description of operations}
 *
 * This collection of methods offers the following operations:
 * \begin{itemize}
 * \item Interpolation: assign each degree of freedom in the vector to be
 *   created the value of the function given as argument. This is identical
 *   to saying that the resulting finite element function (which is isomorphic
 *   to the output vector) has exact function values in all off-points of
 *   ansatz functions. The off-point of an ansatz function is the point where
 *   it assumes its nominal value, e.g. for linear ansatz functions the
 *   off-points are th corners of an element. This function therefore relies
 *   on the assumption that a finite element is used for which the degrees
 *   of freedom are function values (Lagrange elements) rather than gradients,
 *   normal derivatives, second derivatives, etc (Hermite elements, quintic
 *   Argyris element, etc.).
 *
 *   It seems inevitable that some values of the vector to be created are set
 *   twice or even more than that. The reason is that we have to loop over
 *   all cells and get the function values for each of the ansatz functions
 *   located thereon. This applies also to the functions located on faces and
 *   corners which we thus visit more than once. While setting the value
 *   in the vector is not an expensive operation, the evaluation of the
 *   given function may be, taking into account that a virtual function has
 *   to be called.
 *
 * \item Projection: compute the $L_2$-projection of the given function onto
 *   the finite element space. This is done through the solution of the
 *   linear system of equations $M v = f$ where $M$ is the mass matrix
 *   $m_{ij} = \int_\Omega \phi_i(x) \phi_j(x) dx$ and
 *   $f_i = \int_\Omega f(x) \phi_i(x) dx$. The solution vector $v$ then is
 *   the projection.
 *
 *   The solution of the linear system is presently done using a simple CG
 *   method without preconditioning and without multigrid. This is clearly not
 *   too efficient, but sufficient in many cases and simple to implement. This
 *   detail may change in the future.
 * \end{itemize}
 */
template <int dim>
class VectorCreator {
  public:
				     /**
				      * Compute the interpolation of
				      * #function# at the ansatz points to
				      * the finite element space.
				      *
				      * See the general documentation of this
				      * class for further information.
				      */
    static void interpolate (const DoFHandler<dim>    &dof,
			     const FiniteElement<dim> &fe,
			     const Boundary<dim>      &boundary,
			     const Function<dim>      &function,
			     dVector                  &vec);

				     /**
				      * Compute the projection of
				      * #function# to the finite element space.
				      *
				      * See the general documentation of this
				      * class for further information.
				      */
    static void project (const DoFHandler<dim>    &dof,
			 const ConstraintMatrix   &constraints,
			 const FiniteElement<dim> &fe,
			 const Quadrature<dim>    &q,
			 const Boundary<dim>      &boundary,
			 const Function<dim>      &function,
			 dVector                  &vec);
};



/*----------------------------   vectors.h     ---------------------------*/
/* end of #ifndef __vectors_H */
#endif
/*----------------------------   vectors.h     ---------------------------*/
