/*----------------------------   vectors.h     ---------------------------*/
/*      $Id$                 */
#ifndef __vectors_H
#define __vectors_H
/*----------------------------   vectors.h     ---------------------------*/



#include <base/exceptions.h>

template <int dim> class DoFHandler;
template <int dim> class Function;
template <int dim> class Quadrature;
template <int dim> class FiniteElement;
template <int dim> class Boundary;
template <int dim> class StraightBoundary;
class ConstraintMatrix;
class dVector;


/**
 *  Denote which norm/integral is to be computed. The following possibilities
 *  are implemented:
 *  \begin{itemize}
 *  \item #mean#: the function or difference of functions is integrated
 *    on each cell.
 *  \item #L1_norm#: the absolute value of the function is integrated.
 *  \item #L2_norm#: the square of the function is integrated on each
 *    cell; afterwards the root is taken of this value.
 *  \end{itemize}
 */
enum NormType {
      mean,
      L1_norm,
      L2_norm,
      Linfty_norm,
      H1_seminorm,
      H1_norm
};




/**
 * Provide a class which offers some operations on vectors. Amoung these are
 * assemblage of standard vectors, integration of the difference of a
 * finite element solution and a continuous function,
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
 *
 * \item Computing errors:
 *   The function #integrate_difference# performs the calculation of the error
 *   between the finite element solution and a given (continuous) reference
 *   function in different norms. The integration is performed using a given
 *   quadrature formulae and assumes that the given finite element objects equals
 *   that used for the computation of the solution.
 * 
 *   The result ist stored in a vector (named #difference#), where each entry
 *   equals the given norm of the difference on one cell. The order of entries
 *   is the same as a #cell_iterator# takes when started with #begin_active# and
 *   promoted with the #++# operator.
 * 
 *   You can use the #distribute_cell_to_dof_vector# function of the #DoFHandler#
 *   class to convert cell based data to a data vector with values on the degrees
 *   of freedom, which can then be attached to a #DataOut# object to be printed.
 * 
 *   Presently, there is the possibility to compute the following values from the
 *   difference, on each cell: #mean#, #L1_norm#, #L2_norm#, #Linfty_norm#,
 *   #H1_seminorm#.
 *   For the mean difference value, the reference function minus the numerical
 *   solution is computed, not the other way round.
 *
 *   The infinity norm of the difference on a given cell returns the maximum
 *   absolute value of the difference at the quadrature points given by the
 *   quadrature formula parameter. This will in some cases not be too good
 *   an approximation, since for example the Gauss quadrature formulae do
 *   not evaluate the difference at the end or corner points of the cells.
 *   You may want to chose a quadrature formula with more quadrature points
 *   or one with another distribution of the quadrature points in this case.
 *   You should also take into account the superconvergence properties of finite
 *   elements in some points: for example in 1D, the standard finite element
 *   method is a collocation method and should return the exact value at nodal
 *   points. Therefore, the trapezoidal rule should always return a vanishing
 *   L-infinity error. Conversely, in 2D the maximum L-infinity error should
 *   be located at the vertices or at the center of the cell, which would make
 *   it plausible to use the Simpson quadrature rule. On the other hand, there
 *   may be superconvergence at Gauss integration points. These examples are not
 *   intended as a rule of thumb, rather they are though to illustrate that the
 *   use of the wrong quadrature formula may show a significantly wrong result
 *   a nd care should be taken to chose the right formula.
 *
 *   The $H_1$ seminorm is the $L_2$ norm of the gradient of the difference. The
 *   full $H_1$ norm is the sum of the seminorm and the $L_2$ norm.
 * 
 *   To get the {\it global} L_1 error, you have to sum up the entries in
 *   #difference#, e.g. using #dVector::l1_norm# function.
 *   For the global L_2 difference, you have to sum up the squares of the
 *   entries and take the root of the sum, e.g. using #dVector::l2_norm.
 *   These two operations represent the
 *   l_1 and l_2 norms of the vectors, but you need not take the absolute
 *   value of each entry, since the cellwise norms are already positive.
 *  
 *   To get the global mean difference, simply sum up the elements as above.
 *   To get the L_\infty norm, take the maximum of the vector elements, e.g.
 *   using the #dVector::linfty_norm# function.
 *
 *   For the global $H_1$ norm and seminorm, the same rule applies as for the
 *   $L_2$ norm: compute the $l_2$ norm of the cell error vector.
 * \end{itemize}
 */
template <int dim>
class VectorTools {
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

    				     /**
				      * Integrate the difference between
				      * a finite element function and
				      * the reference function, which
				      * is given as a continuous function
				      * object.
				      *
				      * See the general documentation of this
				      * class for more information.
				      */
    static void integrate_difference (const DoFHandler<dim>    &dof,
				      const dVector            &fe_function,
				      const Function<dim>      &exact_solution,
				      dVector                  &difference,
				      const Quadrature<dim>    &q,
				      const FiniteElement<dim> &fe,
				      const NormType           &norm,
				      const Boundary<dim> &boundary=StraightBoundary<dim>());

				     /**
				      * Exception
				      */
    DeclException0 (ExcNotImplemented);
				     /**
				      * Exception
				      */
    DeclException0 (ExcInvalidFE);
};



/*----------------------------   vectors.h     ---------------------------*/
/* end of #ifndef __vectors_H */
#endif
/*----------------------------   vectors.h     ---------------------------*/
