/*----------------------------   vectors.h     ---------------------------*/
/*      $Id$                 */
#ifndef __vectors_H
#define __vectors_H
/*----------------------------   vectors.h     ---------------------------*/



#include <base/exceptions.h>
#include <map>

template <int dim> class DoFHandler;
template <int dim> class Function;
template <int dim> class Quadrature;
template <int dim> class QGauss2;
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
 *   In order to get proper results, it may necessary to treat boundary
 *   conditions right. Below are listed some cases wher this may be needed.
 *   If needed, this is done by $L_2$-projection of the trace of the
 *   given function onto the finite element space restricted to the boundary
 *   of the domain, then taking this information and using it to eliminate
 *   the boundary nodes from the mass matrix of the whole domain, using the
 *   #MatrixTools::apply_boundary_values# function. The projection of the
 *   trace of the function to the boundary is done with the
 *   #VectorTools::project_boundary_values# (see below) function, which is
 *   called with a map of boundary functions in which all boundary indicators
 *   from zero to 254 (255 is used for other purposes, see the #Triangulation#
 *   class documentation) point to the function to be projected. The projection
 *   to the boundary takes place using a second quadrature formula given to
 *   the #project# function.
 *
 *   The projection of the boundary values first, then eliminating them from
 *   the global system of equations is not needed usually. It may be necessary
 *   if you want to enforce special restrictions on the boundary values of the
 *   projected function, for example in time dependant problems: you may want
 *   to project the initial values but need consistency with the boundary
 *   values for later times. Since the latter are projected onto the boundary
 *   in each time step, it is necessary that we also project the boundary
 *   values of the initial values, before projecting them to the whole domain.
 *
 *   Obviously, the results of the two schemes for projection are different.
 *   Usually, when projecting to the boundary first, the $L_2$-norm of the
 *   difference between original function and projection over the whole domain
 *   will be larger (factors of five have been observed) while the $L_2$-norm
 *   of the error integrated over the boundary should of course be less. The
 *   reverse should also hold if no projection to the boundary is performed.
 *
 *   The selection whether the projection to the boundary first is needed is
 *   done with the #project_to_boundary_first# flag passed to the function.
 *   If #false# is given, the additional quadrature formula for faces is
 *   ignored.
 *
 *   You should be aware of the fact that if no projection to the boundary
 *   is requested, a function with zero boundary values may not have zero
 *   boundary values after projection. There is a flag for this especially
 *   important case, which tells the function to enforce zero boundary values
 *   on the respective boundary parts. Since enforced zero boundary values
 *   could also have been reached through projection, but are more economically
 *   obtain using other methods, the #project_to_boundary_first# flag is
 *   ignored if the #enforce_zero_boundary# flag is set.
 *
 *   The solution of the linear system is presently done using a simple CG
 *   method without preconditioning and without multigrid. This is clearly not
 *   too efficient, but sufficient in many cases and simple to implement. This
 *   detail may change in the future.
 *
 * \item Interpolation of boundary values:
 *   The #MatrixTools::apply_boundary_values# function takes a list
 *   of boundary nodes and their values. You can get such a list by interpolation
 *   of a boundary function using the #interpolate_boundary_values# function.
 *   To use it, you have to
 *   specify a list of pairs of boundary indicators (of type #unsigned char#;
 *   see the section in the documentation of the \Ref{Triangulation} class for more
 *   details) and the according functions denoting the dirichlet boundary values
 *   of the nodes on boundary faces with this boundary indicator.
 *
 *   Usually, all other boundary conditions, such as inhomogeneous Neumann values
 *   or mixed boundary conditions are handled in the weak formulation. No attempt
 *   is made to include these into the process of assemblage therefore.
 *
 *   Within this function, boundary values are interpolated, i.e. a node is given
 *   the point value of the boundary function. In some cases, it may be necessary
 *   to use the L2-projection of the boundary function or any other method. For
 *   this purpose to the #VectorTools::project_boundary_values#
 *   function below.
 *
 *   You should be aware that the boundary function may be evaluated at nodes
 *   on the interior of faces. These, however, need not be on the true
 *   boundary, but rather are on the approximation of the boundary represented
 *   by the mapping of the unit cell to the real cell. Since this mapping will
 *   in most cases not be the exact one at the face, the boundary function is
 *   evaluated at points which are not on the boundary and you should make
 *   sure that the returned values are reasonable in some sense anyway.
 *
 * \item Projection of boundary values:
 *   The #project_boundary_values# function acts similar to the
 *   #interpolate_boundary_values# function, apart from the fact that it does
 *   not get the nodal values of boundary nodes by interpolation but rather
 *   through the $L_2$-projection of the trace of the function to the boundary.
 *
 *   The projection takes place on all boundary parts with boundary indicators
 *   listed in the map of boundary functions. These boundary parts may or may
 *   not be contiguous. For these boundary parts, the mass matrix is assembled
 *   using the #MatrixTools::create_boundary_mass_matrix# function, as well as
 *   the appropriate right hand side. Then the resulting system of equations is
 *   solved using a simple CG method (without preconditioning), which is in most
 *   cases sufficient for the present purpose.
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
 *   entries and take the root of the sum, e.g. using #dVector::l2_norm#.
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
				      *	Declare a data type which denotes a
				      *	mapping between a boundary indicator
				      *	and the function denoting the boundary
				      *	values on this part of the boundary.
				      *	Only one boundary function may be given
				      *	for each boundary indicator, which is
				      *	guaranteed by the #map# data type.
				      *	
				      *	See the general documentation of this
				      *	class for more detail.
				      */
    typedef map<unsigned char,const Function<dim>*> FunctionMap;

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
				      * By default, projection to the boundary
				      * and enforcement of zero boundary values
				      * are disabled. The ordering of arguments
				      * to this function is such that you need
				      * not give a second quadrature formula if
				      * you don't want to project to the
				      * boundary first, but that you must if you
				      * want to do so.
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
			 dVector                  &vec,
			 const bool                enforce_zero_boundary = false,
			 const Quadrature<dim-1>  &q_boundary = QGauss2<dim-1>(),
			 const bool                project_to_boundary_first = false);

				     /**
				      * Make up the list of node subject
				      * to Dirichlet boundary conditions
				      * and the values they are to be
				      * assigned, by interpolation around
				      * the boundary. If the
				      * #boundary_values# contained values
				      * before, the new ones are added, or
				      * the old one overwritten if a node
				      * of the boundary part to be prjected
				      * on already was in the variable.
				      *
				      * See the general doc for more
				      * information.
				      */
    static void interpolate_boundary_values (const DoFHandler<dim> &dof,
					     const FunctionMap     &dirichlet_bc,
					     const FiniteElement<dim> &fe,
					     const Boundary<dim> &boundary,
					     map<int,double>     &boundary_values);
    
				     /**
				      * Project #function# to the boundary
				      * of the domain, using the given quadrature
				      * formula for the faces. If the
				      * #boundary_values# contained values
				      * before, the new ones are added, or
				      * the old one overwritten if a node
				      * of the boundary part to be prjected
				      * on already was in the variable.
				      *
				      * See the general documentation of this
				      * class for further information.
				      */
    static void project_boundary_values (const DoFHandler<dim>    &dof,
					 const FunctionMap        &boundary_functions,
					 const FiniteElement<dim> &fe,
					 const Quadrature<dim-1>  &q,
					 const Boundary<dim>      &boundary,
					 map<int,double>          &boundary_values);
    
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
				     /**
				      * Exception
				      */
    DeclException0 (ExcInvalidBoundaryIndicator);
};



/*----------------------------   vectors.h     ---------------------------*/
/* end of #ifndef __vectors_H */
#endif
/*----------------------------   vectors.h     ---------------------------*/
