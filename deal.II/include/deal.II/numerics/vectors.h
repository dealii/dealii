//---------------------------------------------------------------------------
//    $Id$
//
//    Copyright (C) 1998, 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008, 2009, 2010, 2011 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------
#ifndef __deal2__vectors_h
#define __deal2__vectors_h


#include <deal.II/base/config.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/dofs/function_map.h>
#include <deal.II/fe/mapping_q.h>
#include <deal.II/hp/mapping_collection.h>

#include <map>
#include <vector>
#include <set>

DEAL_II_NAMESPACE_OPEN

template <int dim> class Point;
template <int dim> class Function;
template <int dim> class FunctionMap;
template <int dim> class Quadrature;
template <int dim> class QGauss;

template <typename number> class Vector;
template <typename number> class FullMatrix;
template <int dim, int spacedim> class Mapping;
template <int dim, int spacedim> class DoFHandler;
namespace hp
{
  template <int dim, int spacedim> class DoFHandler;
  template <int dim, int spacedim> class MappingCollection;
  template <int dim> class QCollection;
}
class ConstraintMatrix;


//TODO: Move documentation of functions to the functions!
//TODO: (Re)move the basic course on Sobolev spaces

/**
 * Provide a class which offers some operations on vectors. Amoung
 * these are assembling of standard vectors, integration of the
 * difference of a finite element solution and a continuous function,
 * interpolations and projections of continuous functions to the
 * finite element space and other operations.
 *
 * @note There exist two versions of almost each function. One with a
 * Mapping argument and one without. If a code uses a mapping
 * different from MappingQ1 the functions <b>with</b> mapping argument
 * should be used. Code that uses only MappingQ1 may also use the
 * functions without Mapping argument. Each of these latter functions
 * create a MappingQ1 object and just call the respective functions
 * with that object as mapping argument. The functions without Mapping
 * argument still exist to ensure backward compatibility. Nevertheless
 * it is advised to change the user's codes to store a specific
 * Mapping object and to use the functions that take this Mapping
 * object as argument. This gives the possibility to easily extend the
 * user codes to work also on mappings of higher degree, this just by
 * exchanging MappingQ1 by, for example, a MappingQ or another Mapping
 * object of interest.
 *
 * <h3>Description of operations</h3>
 *
 * This collection of methods offers the following operations:
 * <ul>
 * <li> Interpolation: assign each degree of freedom in the vector to be
 *   the value of the function given as argument. This is identical to
 *   saying that the resulting finite element function (which is
 *   isomorphic to the output vector) has exact function values in all
 *   support points of trial functions. The support point of a trial
 *   function is the point where its value equals one, e.g. for linear
 *   trial functions the support points are four corners of an
 *   element. This function therefore relies on the assumption that a
 *   finite element is used for which the degrees of freedom are
 *   function values (Lagrange elements) rather than gradients, normal
 *   derivatives, second derivatives, etc (Hermite elements, quintic
 *   Argyris element, etc.).
 *
 *   It seems inevitable that some values of the vector to be created are set
 *   twice or even more than that. The reason is that we have to loop over
 *   all cells and get the function values for each of the trial functions
 *   located thereon. This applies also to the functions located on faces and
 *   corners which we thus visit more than once. While setting the value
 *   in the vector is not an expensive operation, the evaluation of the
 *   given function may be, taking into account that a virtual function has
 *   to be called.
 *
 * <li> Projection: compute the <i>L</i><sup>2</sup>-projection of the
 * given function onto the finite element space, i.e. if <i>f</i> is
 * the function to be projected, compute <i>f<sub>h</sub></i> in
 * <i>V<sub>h</sub></i> such that
 * (<i>f<sub>h</sub></i>,<i>v<sub>h</sub></i>)=(<i>f</i>,<i>v<sub>h</sub></i>)
 * for all discrete test functions <i>v<sub>h</sub></i>. This is done
 * through the solution of the linear system of equations <i> M v =
 * f</i> where <i>M</i> is the mass matrix $m_{ij} = \int_\Omega
 * \phi_i(x) \phi_j(x) dx$ and $f_i = \int_\Omega f(x) \phi_i(x)
 * dx$. The solution vector $v$ then is the nodal representation of
 * the projection <i>f<sub>h</sub></i>. The project() functions are
 * used in the step-21 and step-23
 * tutorial programs.
 *
 *   In order to get proper results, it be may necessary to treat
 *   boundary conditions right. Below are listed some cases where this
 *   may be needed.  If needed, this is done by <i>L</i><sup>2</sup>-projection of
 *   the trace of the given function onto the finite element space
 *   restricted to the boundary of the domain, then taking this
 *   information and using it to eliminate the boundary nodes from the
 *   mass matrix of the whole domain, using the
 *   MatrixTools::apply_boundary_values() function. The projection of
 *   the trace of the function to the boundary is done with the
 *   VectorTools::project_boundary_values() (see below) function,
 *   which is called with a map of boundary functions FunctioMap in
 *   which all boundary indicators from zero to 254 (255 is used for
 *   other purposes, see the Triangulation class documentation) point
 *   to the function to be projected. The projection to the boundary
 *   takes place using a second quadrature formula on the boundary
 *   given to the project() function. The first quadrature formula is
 *   used to compute the right hand side and for numerical quadrature
 *   of the mass matrix.
 *
 *   The projection of the boundary values first, then eliminating
 *   them from the global system of equations is not needed
 *   usually. It may be necessary if you want to enforce special
 *   restrictions on the boundary values of the projected function,
 *   for example in time dependent problems: you may want to project
 *   the initial values but need consistency with the boundary values
 *   for later times. Since the latter are projected onto the boundary
 *   in each time step, it is necessary that we also project the
 *   boundary values of the initial values, before projecting them to
 *   the whole domain.
 *
 *   Obviously, the results of the two schemes for projection are
 *   different.  Usually, when projecting to the boundary first, the
 *   <i>L</i><sup>2</sup>-norm of the difference between original
 *   function and projection over the whole domain will be larger
 *   (factors of five have been observed) while the
 *   <i>L</i><sup>2</sup>-norm of the error integrated over the
 *   boundary should of course be less. The reverse should also hold
 *   if no projection to the boundary is performed.
 *
 *   The selection whether the projection to the boundary first is
 *   needed is done with the <tt>project_to_boundary_first</tt> flag
 *   passed to the function.  If @p false is given, the additional
 *   quadrature formula for faces is ignored.
 *
 *   You should be aware of the fact that if no projection to the boundary
 *   is requested, a function with zero boundary values may not have zero
 *   boundary values after projection. There is a flag for this especially
 *   important case, which tells the function to enforce zero boundary values
 *   on the respective boundary parts. Since enforced zero boundary values
 *   could also have been reached through projection, but are more economically
 *   obtain using other methods, the @p project_to_boundary_first flag is
 *   ignored if the @p enforce_zero_boundary flag is set.
 *
 *   The solution of the linear system is presently done using a simple CG
 *   method without preconditioning and without multigrid. This is clearly not
 *   too efficient, but sufficient in many cases and simple to implement. This
 *   detail may change in the future.
 *
 * <li> Creation of right hand side vectors:
 *   The create_right_hand_side() function computes the vector
 *   $f_i = \int_\Omega f(x) \phi_i(x) dx$. This is the same as what the
 *   <tt>MatrixCreator::create_*</tt> functions which take a right hand side do,
 *   but without assembling a matrix.
 *
 * <li> Creation of right hand side vectors for point sources:
 *   The create_point_source_vector() function computes the vector
 *   $f_i = \int_\Omega \delta_0(x-x_0) \phi_i(x) dx$.
 *
 * <li> Creation of boundary right hand side vectors: The
 *   create_boundary_right_hand_side() function computes the vector
 *   $f_i = \int_{\partial\Omega} g(x) \phi_i(x) dx$. This is the
 *   right hand side contribution of boundary forces when having
 *   inhomogeneous Neumann boundary values in Laplace's equation or
 *   other second order operators. This function also takes an
 *   optional argument denoting over which parts of the boundary the
 *   integration shall extend.
 *
 * <li> Interpolation of boundary values:
 *   The MatrixTools::apply_boundary_values() function takes a list
 *   of boundary nodes and their values. You can get such a list by interpolation
 *   of a boundary function using the interpolate_boundary_values() function.
 *   To use it, you have to
 *   specify a list of pairs of boundary indicators (of type <tt>unsigned char</tt>;
 *   see the section in the documentation of the Triangulation class for more
 *   details) and the according functions denoting the dirichlet boundary values
 *   of the nodes on boundary faces with this boundary indicator.
 *
 *   Usually, all other boundary conditions, such as inhomogeneous Neumann values
 *   or mixed boundary conditions are handled in the weak formulation. No attempt
 *   is made to include these into the process of matrix and vector assembly therefore.
 *
 *   Within this function, boundary values are interpolated, i.e. a node is given
 *   the point value of the boundary function. In some cases, it may be necessary
 *   to use the L2-projection of the boundary function or any other method. For
 *   this purpose we refer to the project_boundary_values()
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
 *   In 1d the situation is a bit different since there faces (i.e. vertices) have
 *   no boundary indicator. It is assumed that if the boundary indicator zero
 *   is given in the list of boundary functions, the left boundary point is to be
 *   interpolated while the right boundary point is associated with the boundary
 *   index 1 in the map. The respective boundary functions are then evaluated at
 *   the place of the respective boundary point.
 *
 * <li> Projection of boundary values:
 *   The project_boundary_values() function acts similar to the
 *   interpolate_boundary_values() function, apart from the fact that it does
 *   not get the nodal values of boundary nodes by interpolation but rather
 *   through the <i>L</i><sup>2</sup>-projection of the trace of the function to the boundary.
 *
 *   The projection takes place on all boundary parts with boundary
 *   indicators listed in the map (FunctioMap::FunctionMap)
 *   of boundary functions. These boundary parts may or may not be
 *   continuous. For these boundary parts, the mass matrix is
 *   assembled using the
 *   MatrixTools::create_boundary_mass_matrix() function, as
 *   well as the appropriate right hand side. Then the resulting
 *   system of equations is solved using a simple CG method (without
 *   preconditioning), which is in most cases sufficient for the
 *   present purpose.
 *
 * <li> Computing errors:
 *   The function integrate_difference() performs the calculation of
 *   the error between a given (continuous) reference function and the
 *   finite element solution in different norms. The integration is
 *   performed using a given quadrature formula and assumes that the
 *   given finite element objects equals that used for the computation
 *   of the solution.
 *
 *   The result is stored in a vector (named @p difference), where each entry
 *   equals the given norm of the difference on a cell. The order of entries
 *   is the same as a @p cell_iterator takes when started with @p begin_active and
 *   promoted with the <tt>++</tt> operator.
 *
 *   This data, one number per active cell, can be used to generate
 *   graphical output by directly passing it to the DataOut class
 *   through the DataOut::add_data_vector function. Alternatively, it
 *   can be interpolated to the nodal points of a finite element field
 *   using the DoFTools::distribute_cell_to_dof_vector function.
 *
 *   Presently, there is the possibility to compute the following values from the
 *   difference, on each cell: @p mean, @p L1_norm, @p L2_norm, @p Linfty_norm,
 *   @p H1_seminorm and @p H1_norm, see VectorTools::NormType.
 *   For the mean difference value, the reference function minus the numerical
 *   solution is computed, not the other way round.
 *
 *   The infinity norm of the difference on a given cell returns the maximum
 *   absolute value of the difference at the quadrature points given by the
 *   quadrature formula parameter. This will in some cases not be too good
 *   an approximation, since for example the Gauss quadrature formulae do
 *   not evaluate the difference at the end or corner points of the cells.
 *   You may want to choose a quadrature formula with more quadrature points
 *   or one with another distribution of the quadrature points in this case.
 *   You should also take into account the superconvergence properties of finite
 *   elements in some points: for example in 1D, the standard finite element
 *   method is a collocation method and should return the exact value at nodal
 *   points. Therefore, the trapezoidal rule should always return a vanishing
 *   L-infinity error. Conversely, in 2D the maximum L-infinity error should
 *   be located at the vertices or at the center of the cell, which would make
 *   it plausible to use the Simpson quadrature rule. On the other hand, there
 *   may be superconvergence at Gauss integration points. These examples are not
 *   intended as a rule of thumb, rather they are thought to illustrate that the
 *   use of the wrong quadrature formula may show a significantly wrong result
 *   and care should be taken to chose the right formula.
 *
 *   The <i>H</i><sup>1</sup> seminorm is the <i>L</i><sup>2</sup>
 *   norm of the gradient of the difference. The square of the full
 *   <i>H</i><sup>1</sup> norm is the sum of the square of seminorm
 *   and the square of the <i>L</i><sup>2</sup> norm.
 *
 *   To get the global <i>L<sup>1</sup></i> error, you have to sum up the
 *   entries in @p difference, e.g. using
 *   Vector::l1_norm() function.  For the global <i>L</i><sup>2</sup>
 *   difference, you have to sum up the squares of the entries and
 *   take the root of the sum, e.g. using
 *   Vector::l2_norm().  These two operations
 *   represent the <i>l</i><sub>1</sub> and <i>l</i><sub>2</sub> norms of the vectors, but you need
 *   not take the absolute value of each entry, since the cellwise
 *   norms are already positive.
 *
 *   To get the global mean difference, simply sum up the elements as above.
 *   To get the $L_\infty$ norm, take the maximum of the vector elements, e.g.
 *   using the Vector::linfty_norm() function.
 *
 *   For the global <i>H</i><sup>1</sup> norm and seminorm, the same rule applies as for the
 *   <i>L</i><sup>2</sup> norm: compute the <i>l</i><sub>2</sub> norm
 *   of the cell error vector.
 *
 *   Note that, in the codimension one case, if you ask for a norm
 *   that requires the computation of a gradient, then the provided
 *   function is automatically projected along the curve, and the
 *   difference is only computed on the tangential part of the
 *   gradient, since no information is available, on the finite
 *   dimensional one, on the normal component of the gradient.
 * </ul>
 *
 * All functions use the finite element given to the DoFHandler object the last
 * time that the degrees of freedom were distributed over the triangulation. Also,
 * if access to an object describing the exact form of the boundary is needed, the
 * pointer stored within the triangulation object is accessed.
 *
 * @note Instantiations for this template are provided for some vector types,
 * in particular <code>Vector&lt;float&gt;, Vector&lt;double&gt;,
 * BlockVector&lt;float&gt;, BlockVector&lt;double&gt;</code>; others can be
 * generated in application code (see the section on @ref Instantiations in
 * the manual).
 *
 * @ingroup numerics
 * @author Wolfgang Bangerth, Ralf Hartmann, Guido Kanschat, 1998, 1999, 2000, 2001
 */
class VectorTools
{
  public:

				     /**
				      *  Denote which norm/integral is
				      *  to be computed by the
				      *  integrate_difference()
				      *  function of this class. The
				      *  following possibilities are
				      *  implemented:
				     */
  enum NormType
    {
					   /**
					    * The function or
					    * difference of functions
					    * is integrated on each
					    * cell.
					    */
      mean,
					   /**
					    * The absolute value of
					    * the function is
					    * integrated.
					    */
      L1_norm,
					   /**
					    * The square of the
					    * function is integrated
					    * and the the square root
					    * of the result is
					    * computed on each cell.
					    */
      L2_norm,
					   /**
					    * The absolute value to
					    * the <i>p</i>th power is
					    * integrated and the pth
					    * root is computed on each
					    * cell. The exponent
					    * <i>p</i> is the last
					    * parameter of the
					    * function.
					    */
      Lp_norm,
					   /**
					    * The maximum absolute
					    * value of the function.
					    */
      Linfty_norm,
					   /**
					    * #L2_norm of the gradient.
					    */
      H1_seminorm,
					   /**
					    * The square of this norm
					    * is the square of the
					    * #L2_norm plus the square
					    * of the #H1_seminorm.
					    */
      H1_norm,
					   /**
					    * #Lp_norm of the gradient.
					    */
      W1p_seminorm,
					   /**
					    * same as #H1_norm for
					    * <i>L<sup>p</sup></i>.
					    */
      W1p_norm,
					   /**
					    * #Linfty_norm of the gradient.
					    */
      W1infty_seminorm,
					   /**
					    * same as #H1_norm for
					    * <i>L<sup>infty</sup></i>.
					    */
      W1infty_norm

    };
/**
 * @name Interpolation and projection
 */
				     //@{
				     /**
				      * Compute the interpolation of
				      * @p function at the support
				      * points to the finite element
				      * space. It is assumed that the
				      * number of components of
				      * @p function matches that of
				      * the finite element used by
				      * @p dof.
				      *
				      * Note that you may have to call
				      * <tt>hanging_nodes.distribute(vec)</tt>
				      * with the hanging nodes from
				      * space @p dof afterwards, to
				      * make the result continuous
				      * again.
				      *
				      * The template argument <code>DH</code>
				      * may either be of type DoFHandler or
				      * hp::DoFHandler.
				      *
				      * See the general documentation
				      * of this class for further
				      * information.
				      *
				      * @todo The @p mapping argument should be
				      * replaced by a hp::MappingCollection in
				      * case of a hp::DoFHandler.
				      */
  template <class VECTOR, class DH>
    static void interpolate (const Mapping<DH::dimension,DH::space_dimension>    &mapping,
			     const DH              &dof,
			     const Function<DH::space_dimension>   &function,
			     VECTOR                &vec);

				     /**
				      * Calls the @p interpolate()
				      * function above with
				      * <tt>mapping=MappingQ1@<dim>@()</tt>.
				      */
  template <class VECTOR, class DH>
  static void interpolate (const DH              &dof,
			   const Function<DH::space_dimension>   &function,
			   VECTOR                &vec);

				     /**
				      * Interpolate different finite
				      * element spaces. The
				      * interpolation of vector
				      * @p data_1 is executed from the
				      * FE space represented by
				      * @p dof_1 to the vector @p data_2
				      * on FE space @p dof_2. The
				      * interpolation on each cell is
				      * represented by the matrix
				      * @p transfer. Curved boundaries
				      * are neglected so far.
				      *
				      * Note that you may have to call
				      * <tt>hanging_nodes.distribute(data_2)</tt>
				      * with the hanging nodes from
				      * space @p dof_2 afterwards, to
				      * make the result continuous
				      * again.
				      *
				      * @note Instantiations for this template
				      * are provided for some vector types
				      * (see the general documentation of the
				      * class), but only the same vector for
				      * InVector and OutVector. Other
				      * combinations must be instantiated by
				      * hand.
				      */
  template <int dim, class InVector, class OutVector, int spacedim>
  static void interpolate (const DoFHandler<dim,spacedim>    &dof_1,
			   const DoFHandler<dim,spacedim>    &dof_2,
			   const FullMatrix<double> &transfer,
			   const InVector           &data_1,
			   OutVector                &data_2);

				     /**
				      * Compute the projection of
				      * @p function to the finite element space.
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
				      * This function needs the mass
				      * matrix of the finite element
				      * space on the present grid. To
				      * this end, the mass matrix is
				      * assembled exactly using
				      * MatrixTools::create_mass_matrix. This
				      * function performs numerical
				      * quadrature using the given
				      * quadrature rule; you should
				      * therefore make sure that the
				      * given quadrature formula is
				      * also sufficient for the
				      * integration of the mass
				      * matrix.
				      *
				      * See the general documentation of this
				      * class for further information.
				      *
				      * In 1d, the default value of
				      * the boundary quadrature
				      * formula is an invalid object
				      * since integration on the
				      * boundary doesn't happen in
				      * 1d.
				      */
  template <int dim, class VECTOR, int spacedim>
  static void project (const Mapping<dim, spacedim>       &mapping,
		       const DoFHandler<dim,spacedim>    &dof,
		       const ConstraintMatrix   &constraints,
		       const Quadrature<dim>    &quadrature,
		       const Function<spacedim>      &function,
		       VECTOR                   &vec,
		       const bool                enforce_zero_boundary = false,
		       const Quadrature<dim-1>  &q_boundary = (dim > 1 ?
							       QGauss<dim-1>(2) :
							       Quadrature<dim-1>(0)),
		       const bool                project_to_boundary_first = false);

				     /**
				      * Calls the project()
				      * function, see above, with
				      * <tt>mapping=MappingQ1@<dim@>()</tt>.
				      */
  template <int dim, class VECTOR, int spacedim>
  static void project (const DoFHandler<dim,spacedim>    &dof,
		       const ConstraintMatrix   &constraints,
		       const Quadrature<dim>    &quadrature,
		       const Function<spacedim>      &function,
		       VECTOR                   &vec,
		       const bool                enforce_zero_boundary = false,
		       const Quadrature<dim-1>  &q_boundary = (dim > 1 ?
							       QGauss<dim-1>(2) :
							       Quadrature<dim-1>(0)),
		       const bool                project_to_boundary_first = false);

				     /**
				      * Prepare Dirichlet boundary
				      * conditions.  Make up the list
				      * of degrees of freedom subject
				      * to Dirichlet boundary
				      * conditions and the values to
				      * be assigned to them, by
				      * interpolation around the
				      * boundary. If the
				      * @p boundary_values contained
				      * values before, the new ones
				      * are added, or the old ones
				      * overwritten if a node of the
				      * boundary part to be used
				      * was already in the
				      * map of boundary values.
				      *
				      * The parameter
				      * @p boundary_component
				      * corresponds to the number
				      * @p boundary_indicator of the
				      * face.  255 is an illegal
				      * value, since it is reserved
				      * for interior faces.
				      *
				      * The flags in the last
				      * parameter, @p component_mask
				      * denote which components of the
				      * finite element space shall be
				      * interpolated. If it is left as
				      * specified by the default value
				      * (i.e. an empty array), all
				      * components are
				      * interpolated. If it is
				      * different from the default
				      * value, it is assumed that the
				      * number of entries equals the
				      * number of components in the
				      * boundary functions and the
				      * finite element, and those
				      * components in the given
				      * boundary function will be used
				      * for which the respective flag
				      * was set in the component mask.
				      *
				      * It is assumed that the number
				      * of components of the function
				      * in @p boundary_function matches that
				      * of the finite element used by
				      * @p dof.
				      *
				      * If the finite element used has
				      * shape functions that are
				      * non-zero in more than one
				      * component (in deal.II speak:
				      * they are non-primitive), then
				      * these components can presently
				      * not be used for interpolating
				      * boundary values. Thus, the
				      * elements in the component mask
				      * corresponding to the
				      * components of these
				      * non-primitive shape functions
				      * must be @p false.
				      *
				      * See the general doc for more
				      * information.
				      */
  template <class DH>
  static
  void
  interpolate_boundary_values (const Mapping<DH::dimension,DH::space_dimension>            &mapping,
			       const DH                 &dof,
			       const typename FunctionMap<DH::space_dimension>::type &function_map,
			       std::map<unsigned int,double> &boundary_values,
			       const std::vector<bool>       &component_mask = std::vector<bool>());

				     /**
				      * @deprecated This function exists mainly
				      * for backward compatibility.
				      *
				      * Same function as above, but
				      * taking only one pair of
				      * boundary indicator and
				      * corresponding boundary
				      * function. Calls the other
				      * function with remapped
				      * arguments.
				      *
				      */
  template <class DH>
  static
  void
  interpolate_boundary_values (const Mapping<DH::dimension,DH::space_dimension>            &mapping,
			       const DH                 &dof,
			       const unsigned char            boundary_component,
			       const Function<DH::space_dimension>           &boundary_function,
			       std::map<unsigned int,double> &boundary_values,
			       const std::vector<bool>       &component_mask = std::vector<bool>());

				     /**
				      * Calls the other
				      * interpolate_boundary_values()
				      * function, see above, with
				      * <tt>mapping=MappingQ1@<dim@>()</tt>.
				      */
  template <class DH>
  static
  void
  interpolate_boundary_values (const DH        &dof,
			       const unsigned char            boundary_component,
			       const Function<DH::space_dimension>           &boundary_function,
			       std::map<unsigned int,double> &boundary_values,
			       const std::vector<bool>       &component_mask = std::vector<bool>());


				     /**
				      * Calls the other
				      * interpolate_boundary_values()
				      * function, see above, with
				      * <tt>mapping=MappingQ1@<dim@>()</tt>.
				      */
  template <class DH>
  static
  void
  interpolate_boundary_values (const DH                &dof,
			       const typename FunctionMap<DH::space_dimension>::type &function_map,
			       std::map<unsigned int,double>         &boundary_values,
			       const std::vector<bool>               &component_mask = std::vector<bool>());


				     /**
				      * Insert the (algebraic) constraints due
				      * to Dirichlet boundary conditions into
				      * a ConstraintMatrix @p
				      * constraints. This function identifies
				      * the degrees of freedom subject to
				      * Dirichlet boundary conditions, adds
				      * them to the list of constrained DoFs
				      * in @p constraints and sets the
				      * respective inhomogeneity to the value
				      * interpolated around the boundary. If
				      * this routine encounters a DoF that
				      * already is constrained (for instance
				      * by a hanging node constraint, see
				      * below, or any other type of
				      * constraint, e.g. from periodic
				      * boundary conditions), the old setting
				      * of the constraint (dofs the entry is
				      * constrained to, inhomogeneities) is
				      * kept and nothing happens.
				      *
				      * @note When combining adaptively
				      * refined meshes with hanging node
				      * constraints and boundary conditions
				      * like from the current function within
				      * one ConstraintMatrix object, the
				      * hanging node constraints should always
				      * be set first, and then the boundary
				      * conditions since boundary conditions
				      * are not set in the second operation on
				      * degrees of freedom that are already
				      * constrained. This makes sure that the
				      * discretization remains conforming as
				      * is needed. See the discussion on
				      * conflicting constraints in the module
				      * on @ref constraints .
				      *
				      * The parameter @p boundary_component
				      * corresponds to the number @p
				      * boundary_indicator of the face.  255
				      * is an illegal value, since it is
				      * reserved for interior faces.
				      *
				      * The flags in the last parameter, @p
				      * component_mask denote which
				      * components of the finite element
				      * space shall be interpolated. If it
				      * is left as specified by the default
				      * value (i.e. an empty array), all
				      * components are interpolated. If it
				      * is different from the default value,
				      * it is assumed that the number of
				      * entries equals the number of
				      * components in the boundary functions
				      * and the finite element, and those
				      * components in the given boundary
				      * function will be used for which the
				      * respective flag was set in the
				      * component mask.
				      *
				      * It is assumed that the number of
				      * components of the function in @p
				      * boundary_function matches that of
				      * the finite element used by @p dof.
				      *
				      * If the finite element used has shape
				      * functions that are non-zero in more
				      * than one component (in deal.II
				      * speak: they are non-primitive), then
				      * these components can presently not
				      * be used for interpolating boundary
				      * values. Thus, the elements in the
				      * component mask corresponding to the
				      * components of these non-primitive
				      * shape functions must be @p false.
				      *
				      * See the general doc for more
				      * information.
				      *
				      * @ingroup constraints
				      */
  template <class DH>
  static
  void
  interpolate_boundary_values (const Mapping<DH::dimension,DH::space_dimension>            &mapping,
			       const DH                 &dof,
			       const typename FunctionMap<DH::space_dimension>::type &function_map,
			       ConstraintMatrix              &constraints,
			       const std::vector<bool>       &component_mask = std::vector<bool>());

				     /**
				      * @deprecated This function is there
				      * mainly for backward compatibility.
				      *
				      * Same function as above, but taking
				      * only one pair of boundary indicator
				      * and corresponding boundary
				      * function. Calls the other function
				      * with remapped arguments.
				      *
				      * @ingroup constraints
				      */
  template <class DH>
  static
  void
  interpolate_boundary_values (const Mapping<DH::dimension,DH::space_dimension> &mapping,
			       const DH                            &dof,
			       const unsigned char                  boundary_component,
			       const Function<DH::space_dimension> &boundary_function,
			       ConstraintMatrix                    &constraints,
			       const std::vector<bool>             &component_mask = std::vector<bool>());

				     /**
				      * Calls the other
				      * interpolate_boundary_values()
				      * function, see above, with
				      * <tt>mapping=MappingQ1@<dim@>()</tt>.
				      *
				      * @ingroup constraints
				      */
  template <class DH>
  static
  void
  interpolate_boundary_values (const DH                            &dof,
			       const unsigned char                  boundary_component,
			       const Function<DH::space_dimension> &boundary_function,
			       ConstraintMatrix                    &constraints,
			       const std::vector<bool>             &component_mask = std::vector<bool>());


				     /**
				      * Calls the other
				      * interpolate_boundary_values()
				      * function, see above, with
				      * <tt>mapping=MappingQ1@<dim@>()</tt>.
				      *
				      * @ingroup constraints
				      */
  template <class DH>
  static
  void
  interpolate_boundary_values (const DH                &dof,
			       const typename FunctionMap<DH::space_dimension>::type &function_map,
			       ConstraintMatrix        &constraints,
			       const std::vector<bool> &component_mask = std::vector<bool>());


				     /**
				      * Project a function to the boundary
				      * of the domain, using the given
				      * quadrature formula for the faces. If
				      * the @p boundary_values contained
				      * values before, the new ones are
				      * added, or the old one overwritten if
				      * a node of the boundary part to be
				      * projected on already was in the
				      * variable.
				      *
				      * If @p component_mapping is empty, it
				      * is assumed that the number of
				      * components of @p boundary_function
				      * matches that of the finite element
				      * used by @p dof.
				      *
				      * In 1d, projection equals
				      * interpolation. Therefore,
				      * interpolate_boundary_values is
				      * called.
				      *
				      * @arg @p boundary_values: the result
				      * of this function, a map containing
				      * all indices of degrees of freedom at
				      * the boundary (as covered by the
				      * boundary parts in @p
				      * boundary_functions) and the computed
				      * dof value for this degree of
				      * freedom.
				      *
				      * @arg @p component_mapping: if the
				      * components in @p boundary_functions
				      * and @p dof do not coincide, this
				      * vector allows them to be
				      * remapped. If the vector is not
				      * empty, it has to have one entry for
				      * each component in @p dof. This entry
				      * is the component number in @p
				      * boundary_functions that should be
				      * used for this component in @p
				      * dof. By default, no remapping is
				      * applied.
				      */
  template <int dim, int spacedim>
  static void project_boundary_values (const Mapping<dim, spacedim>       &mapping,
				       const DoFHandler<dim,spacedim>    &dof,
				       const typename FunctionMap<spacedim>::type &boundary_functions,
				       const Quadrature<dim-1>  &q,
				       std::map<unsigned int,double> &boundary_values,
 				       std::vector<unsigned int> component_mapping = std::vector<unsigned int>());

				     /**
				      * Calls the project_boundary_values()
				      * function, see above, with
				      * <tt>mapping=MappingQ1@<dim@>()</tt>.
				      */
  template <int dim, int spacedim>
  static void project_boundary_values (const DoFHandler<dim,spacedim>    &dof,
				       const typename FunctionMap<spacedim>::type &boundary_function,
				       const Quadrature<dim-1>  &q,
				       std::map<unsigned int,double> &boundary_values,
				       std::vector<unsigned int> component_mapping = std::vector<unsigned int>());

				     /**
				      * Project a function to the boundary of
				      * the domain, using the given quadrature
				      * formula for the faces. This function
				      * identifies the degrees of freedom
				      * subject to Dirichlet boundary
				      * conditions, adds them to the list of
				      * constrained DoFs in @p constraints and
				      * sets the respective inhomogeneity to
				      * the value resulting from the
				      * projection operation. If this routine
				      * encounters a DoF that already is
				      * constrained (for instance by a hanging
				      * node constraint, see below, or any
				      * other type of constraint, e.g. from
				      * periodic boundary conditions), the old
				      * setting of the constraint (dofs the
				      * entry is constrained to,
				      * inhomogeneities) is kept and nothing
				      * happens.
				      *
				      * @note When combining adaptively
				      * refined meshes with hanging node
				      * constraints and boundary conditions
				      * like from the current function within
				      * one ConstraintMatrix object, the
				      * hanging node constraints should always
				      * be set first, and then the boundary
				      * conditions since boundary conditions
				      * are not set in the second operation on
				      * degrees of freedom that are already
				      * constrained. This makes sure that the
				      * discretization remains conforming as
				      * is needed. See the discussion on
				      * conflicting constraints in the module
				      * on @ref constraints .
				      *
				      * If @p component_mapping is empty, it
				      * is assumed that the number of
				      * components of @p boundary_function
				      * matches that of the finite element
				      * used by @p dof.
				      *
				      * In 1d, projection equals
				      * interpolation. Therefore,
				      * interpolate_boundary_values is
				      * called.
				      *
				      * @arg @p component_mapping: if the
				      * components in @p boundary_functions
				      * and @p dof do not coincide, this
				      * vector allows them to be
				      * remapped. If the vector is not
				      * empty, it has to have one entry for
				      * each component in @p dof. This entry
				      * is the component number in @p
				      * boundary_functions that should be
				      * used for this component in @p
				      * dof. By default, no remapping is
				      * applied.
				      *
				      * @ingroup constraints
				      */
  template <int dim, int spacedim>
  static void project_boundary_values (const Mapping<dim, spacedim>   &mapping,
				       const DoFHandler<dim,spacedim> &dof,
				       const typename FunctionMap<spacedim>::type &boundary_functions,
				       const Quadrature<dim-1>        &q,
				       ConstraintMatrix               &constraints,
 				       std::vector<unsigned int>       component_mapping = std::vector<unsigned int>());

				     /**
				      * Calls the project_boundary_values()
				      * function, see above, with
				      * <tt>mapping=MappingQ1@<dim@>()</tt>.
				      *
				      * @ingroup constraints
				      */
  template <int dim, int spacedim>
  static void project_boundary_values (const DoFHandler<dim,spacedim> &dof,
				       const typename FunctionMap<spacedim>::type &boundary_function,
				       const Quadrature<dim-1>        &q,
				       ConstraintMatrix               &constraints,
				       std::vector<unsigned int>       component_mapping = std::vector<unsigned int>());


				     /**
				      * Compute constraints that correspond to
				      * boundary conditions of the form
				      * $\vec{n}\times\vec{u}=\vec{n}\times\vec{f}$,
				      * i.e. the tangential components of $u$
				      * and $f$ shall coincide.
				      *
				      * If the ConstraintMatrix @p constraints
				      * contained values or other
				      * constraints before, the new ones are
				      * added or the old ones overwritten,
				      * if a node of the boundary part to be
				      * used was already in the list of
				      * constraints. This is handled by
				      * using inhomogeneous constraints. Please
				      * note that when combining adaptive meshes
				      * and this kind of constraints, the
				      * Dirichlet conditions should be set
				      * first, and then completed by hanging
				      * node constraints, in order to make sure
				      * that the discretization remains
				      * consistent. See the discussion on
				      * conflicting constraints in the
				      * module on @ref constraints .
				      *
				      * This function is explecitly written to
				      * use with the FE_Nedelec elements. Thus
				      * it throws an exception, if it is
				      * called with other finite elements.
				      *
				      * The second argument of this function
				      * denotes the first vector component in
				      * the finite element that corresponds to
				      * the vector function that you want to
				      * constrain. For example, if we want to
				      * solve Maxwell's equations in 3d and the
				      * finite element has components
				      * $(E_x,E_y,E_z,B_x,B_y,B_z)$ and we want
				      * the boundary conditions
				      * $\vec{n}\times\vec{B}=\vec{n}\times\vec{f}$,
				      * then @p first_vector_component would
				      * be 3. Vectors are implicitly assumed to
				      * have exactly <code>dim</code> components
				      * that are ordered in the same way as we
				      * usually order the coordinate directions,
				      * i.e. $x$-, $y$-, and finally
				      * $z$-component.
				      *
				      * The parameter @p boundary_component
				      * corresponds to the number
				      * @p boundary_indicator of the face. 255
				      * is an illegal value, since it is
				      * reserved for interior faces.
				      *
				      * The last argument is denoted to compute
				      * the normal vector $\vec{n}$ at the
				      * boundary points.
				      *
				      * <h4>Computing constraints</h4>
				      *
				      * To compute the constraints we use
				      * projection-based interpolation as proposed
				      * in Solin, Segeth and Dolezel
				      * (Higher order finite elements, Chapman&amp;Hall,
				      * 2004) on every face located at the
				      * boundary.
				      *
				      * First one projects $\vec{f}$ on the
				      * lowest-order edge shape functions. Then the
				      * remaining part $(I-P_0)\vec{f}$ of the
				      * function is projected on the remaining
				      * higher-order edge shape functions. In the
				      * last step we project $(I-P_0-P_e)\vec{f}$
				      * on the bubble shape functions defined on
				      * the face.
				      *
				      * @ingroup constraints
				      */
   template <int dim>
   static void project_boundary_values_curl_conforming (const DoFHandler<dim>& dof_handler,
     const unsigned int first_vector_component,
     const Function<dim>& boundary_function,
     const unsigned char boundary_component,
     ConstraintMatrix& constraints,
     const Mapping<dim>& mapping = StaticMappingQ1<dim>::mapping);

				     /**
				      * Same as above for the hp-namespace.
				      *
				      * @ingroup constraints
				      */
   template <int dim>
   static void project_boundary_values_curl_conforming (const hp::DoFHandler<dim>& dof_handler,
     const unsigned int first_vector_component,
     const Function<dim>& boundary_function,
     const unsigned char boundary_component,
     ConstraintMatrix& constraints,
     const hp::MappingCollection<dim, dim>& mapping_collection = hp::StaticMappingQ1<dim>::mapping_collection);


				     /**
				      * Compute constraints that correspond to
				      * boundary conditions of the form
				      * $\vec{n}^T\vec{u}=\vec{n}^T\vec{f}$,
				      * i.e. the normal components of $u$
				      * and $f$ shall coincide.
				      *
				      * If the ConstraintMatrix @p constraints
				      * contained values or other
				      * constraints before, the new ones are
				      * added or the old ones overwritten,
				      * if a node of the boundary part to be
				      * used was already in the list of
				      * constraints. This is handled by
				      * using inhomogeneous constraints. Please
				      * note that when combining adaptive meshes
				      * and this kind of constraints, the
				      * Dirichlet conditions should be set
				      * first, and then completed by hanging
				      * node constraints, in order to make sure
				      * that the discretization remains
				      * consistent. See the discussion on
				      * conflicting constraints in the
				      * module on @ref constraints .
				      *
				      * This function is explecitly written to
				      * use with the FE_RaviartThomas elements.
				      * Thus it throws an exception, if it is
				      * called with other finite elements.
				      *
				      * The second argument of this function
				      * denotes the first vector component in
				      * the finite element that corresponds to
				      * the vector function that you want to
				      * constrain. Vectors are implicitly
				      * assumed to have exactly
				      * <code>dim</code> components that are
				      * ordered in the same way as we
				      * usually order the coordinate directions,
				      * i.e. $x$-, $y$-, and finally
				      * $z$-component.
				      *
				      * The parameter @p boundary_component
				      * corresponds to the number
				      * @p boundary_indicator of the face. 255
				      * is an illegal value, since it is
				      * reserved for interior faces.
				      *
				      * The last argument is denoted to compute
				      * the normal vector $\vec{n}$ at the
				      * boundary points.
				      *
				      * <h4>Computing constraints</h4>
				      *
				      * To compute the constraints we use
				      * interpolation operator proposed
				      * in Brezzi, Fortin (Mixed and Hybrid
				      * (Finite Element Methods, Springer,
				      * 1991) on every face located at the
				      * boundary.
				      *
				      * @ingroup constraints
				      */
   template<int dim>
   static void project_boundary_values_div_conforming (const DoFHandler<dim>& dof_handler,
     const unsigned int first_vector_component,
     const Function<dim>& boundary_function,
     const unsigned char boundary_component,
     ConstraintMatrix& constraints,
     const Mapping<dim>& mapping = StaticMappingQ1<dim>::mapping);

				     /**
				      * Same as above for the hp-namespace.
				      *
				      * @ingroup constraints
				      */
   template<int dim>
   static void project_boundary_values_div_conforming (const hp::DoFHandler<dim>& dof_handler,
     const unsigned int first_vector_component,
     const Function<dim>& boundary_function,
     const unsigned char boundary_component,
     ConstraintMatrix& constraints,
     const hp::MappingCollection<dim, dim>& mapping_collection = hp::StaticMappingQ1<dim>::mapping_collection);
   
   
				     /**
				      * Compute the constraints that
				      * correspond to boundary conditions of
				      * the form $\vec n \cdot \vec u=0$,
				      * i.e. no normal flux if $\vec u$ is a
				      * vector-valued quantity. These
				      * conditions have exactly the form
				      * handled by the ConstraintMatrix class,
				      * so instead of creating a map between
				      * boundary degrees of freedom and
				      * corresponding value, we here create a
				      * list of constraints that are written
				      * into a ConstraintMatrix. This object
				      * may already have some content, for
				      * example from hanging node constraints,
				      * that remains untouched. These
				      * constraints have to be applied to the
				      * linear system like any other such
				      * constraints, i.e. you have to condense
				      * the linear system with the constraints
				      * before solving, and you have to
				      * distribute the solution vector
				      * afterwards.
				      *
				      * The use of this function is
				      * explained in more detail in
				      * step-31. It
				      * doesn't make much sense in 1d,
				      * so the function throws an
				      * exception in that case.
				      *
				      * The second argument of this
				      * function denotes the first
				      * vector component in the finite
				      * element that corresponds to
				      * the vector function that you
				      * want to constrain. For
				      * example, if we were solving a
				      * Stokes equation in 2d and the
				      * finite element had components
				      * $(u,v,p)$, then @p
				      * first_vector_component would
				      * be zero. On the other hand, if
				      * we solved the Maxwell
				      * equations in 3d and the finite
				      * element has components
				      * $(E_x,E_y,E_z,B_x,B_y,B_z)$
				      * and we want the boundary
				      * condition $\vec n\cdot \vec
				      * B=0$, then @p
				      * first_vector_component would
				      * be 3. Vectors are implicitly
				      * assumed to have exactly
				      * <code>dim</code> components
				      * that are ordered in the same
				      * way as we usually order the
				      * coordinate directions,
				      * i.e. $x$-, $y$-, and finally
				      * $z$-component. The function
				      * assumes, but can't check, that
				      * the vector components in the
				      * range
				      * <code>[first_vector_component,first_vector_component+dim)</code>
				      * come from the same base finite
				      * element. For example, in the
				      * Stokes example above, it would
				      * not make sense to use a
				      * <code>FESystem@<dim@>(FE_Q@<dim@>(2),
				      * 1, FE_Q@<dim@>(1), dim)</code>
				      * (note that the first velocity
				      * vector component is a $Q_2$
				      * element, whereas all the other
				      * ones are $Q_1$ elements) as
				      * there would be points on the
				      * boundary where the
				      * $x$-velocity is defined but no
				      * corresponding $y$- or
				      * $z$-velocities.
				      *
				      * The third argument denotes the set of
				      * boundary indicators on which the
				      * boundary condition is to be
				      * enforced. Note that, as explained
				      * below, this is one of the few
				      * functions where it makes a difference
				      * where we call the function multiple
				      * times with only one boundary
				      * indicator, or whether we call the
				      * function onces with the whole set of
				      * boundary indicators at once.
				      *
				      * The last argument is denoted to
				      * compute the normal vector $\vec n$ at
				      * the boundary points.
				      *
				      * @note When combining adaptively
				      * refined meshes with hanging node
				      * constraints and boundary conditions
				      * like from the current function within
				      * one ConstraintMatrix object, the
				      * hanging node constraints should always
				      * be set first, and then the boundary
				      * conditions since boundary conditions
				      * are not set in the second operation on
				      * degrees of freedom that are already
				      * constrained. This makes sure that the
				      * discretization remains conforming as
				      * is needed. See the discussion on
				      * conflicting constraints in the module
				      * on @ref constraints .
				      *
				      *
				      * <h4>Computing constraints in 2d</h4>
				      *
				      * Computing these constraints requires
				      * some smarts. The main question
				      * revolves around the question what the
				      * normal vector is. Consider the
				      * following situation:
				      * <p ALIGN="center">
				      * @image html no_normal_flux_1.png
				      * </p>
				      *
				      * Here, we have two cells that use a
				      * bilinear mapping
				      * (i.e. MappingQ1). Consequently, for
				      * each of the cells, the normal vector
				      * is perpendicular to the straight
				      * edge. If the two edges at the top and
				      * right are meant to approximate a
				      * curved boundary (as indicated by the
				      * dashed line), then neither of the two
				      * computed normal vectors are equal to
				      * the exact normal vector (though they
				      * approximate it as the mesh is refined
				      * further). What is worse, if we
				      * constrain $\vec n \cdot \vec u=0$ at
				      * the common vertex with the normal
				      * vector from both cells, then we
				      * constrain the vector $\vec u$ with
				      * respect to two linearly independent
				      * vectors; consequently, the constraint
				      * would be $\vec u=0$ at this point
				      * (i.e. <i>all</i> components of the
				      * vector), which is not what we wanted.
				      *
				      * To deal with this situation, the
				      * algorithm works in the following way:
				      * at each point where we want to
				      * constrain $\vec u$, we first collect
				      * all normal vectors that adjacent cells
				      * might compute at this point. We then
				      * do not constrain $\vec n \cdot \vec
				      * u=0$ for <i>each</i> of these normal
				      * vectors but only for the
				      * <i>average</i> of the normal
				      * vectors. In the example above, we
				      * therefore record only a single
				      * constraint $\vec n \cdot \vec {\bar
				      * u}=0$, where $\vec {\bar u}$ is the
				      * average of the two indicated normal
				      * vectors.
				      *
				      * Unfortunately, this is not quite
				      * enough. Consider the situation here:
				      *
				      * <p ALIGN="center">
				      * @image html no_normal_flux_2.png
				      * </p>
				      *
				      * If again the top and right edges
				      * approximate a curved boundary, and the
				      * left boundary a separate boundary (for
				      * example straight) so that the exact
				      * boundary has indeed a corner at the
				      * top left vertex, then the above
				      * construction would not work: here, we
				      * indeed want the constraint that $\vec
				      * u$ at this point (because the normal
				      * velocities with respect to both the
				      * left normal as well as the top normal
				      * vector should be zero), not that the
				      * velocity in the direction of the
				      * average normal vector is zero.
				      *
				      * Consequently, we use the following
				      * heuristic to determine whether all
				      * normal vectors computed at one point
				      * are to be averaged: if two normal
				      * vectors for the same point are
				      * computed on <i>different</i> cells,
				      * then they are to be averaged. This
				      * covers the first example above. If
				      * they are computed from the same cell,
				      * then the fact that they are different
				      * is considered indication that they
				      * come from different parts of the
				      * boundary that might be joined by a
				      * real corner, and must not be averaged.
				      *
				      * There is one problem with this
				      * scheme. If, for example, the same
				      * domain we have considered above, is
				      * discretized with the following mesh,
				      * then we get into trouble:
				      *
				      * <p ALIGN="center">
				      * @image html no_normal_flux_2.png
				      * </p>
				      *
				      * Here, the algorithm assumes that the
				      * boundary does not have a corner at the
				      * point where faces $F1$ and $F2$ join
				      * because at that point there are two
				      * different normal vectors computed from
				      * different cells. If you intend for
				      * there to be a corner of the exact
				      * boundary at this point, the only way
				      * to deal with this is to assign the two
				      * parts of the boundary different
				      * boundary indicators and call this
				      * function twice, once for each boundary
				      * indicators; doing so will yield only
				      * one normal vector at this point per
				      * invocation (because we consider only
				      * one boundary part at a time), with the
				      * result that the normal vectors will
				      * not be averaged.
				      *
				      *
				      * <h4>Computing constraints in 3d</h4>
				      *
				      * The situation is more
				      * complicated in 3d. Consider
				      * the following case where we
				      * want to compute the
				      * constraints at the marked
				      * vertex:
				      *
				      * <p ALIGN="center">
				      * @image html no_normal_flux_4.png
				      * </p>
				      *
				      * Here, we get four different
				      * normal vectors, one from each
				      * of the four faces that meet at
				      * the vertex. Even though they
				      * may form a complete set of
				      * vectors, it is not our intent
				      * to constrain all components of
				      * the vector field at this
				      * point. Rather, we would like
				      * to still allow tangential
				      * flow, where the term
				      * "tangential" has to be
				      * suitably defined.
				      *
				      * In a case like this, the
				      * algorithm proceeds as follows:
				      * for each cell that has
				      * computed two tangential
				      * vectors at this point, we
				      * compute the unconstrained
				      * direction as the outer product
				      * of the two tangential vectors
				      * (if necessary multiplied by
				      * minus one). We then average
				      * these tangential
				      * vectors. Finally, we compute
				      * constraints for the two
				      * directions perpendicular to
				      * this averaged tangential
				      * direction.
				      *
				      * There are cases where one cell
				      * contributes two tangential
				      * directions and another one
				      * only one; for example, this
				      * would happen if both top and
				      * front faces of the left cell
				      * belong to the boundary
				      * selected whereas only the top
				      * face of the right cell belongs
				      * to it. This case is not
				      * currently implemented.
				      *
				      *
				      * <h4>Results</h4>
				      *
				      * Because it makes for good
				      * pictures, here are two images
				      * of vector fields on a circle
				      * and on a sphere to which the
				      * constraints computed by this
				      * function have been applied:
				      *
				      * <p ALIGN="center">
				      * @image html no_normal_flux_5.png
				      * @image html no_normal_flux_6.png
				      * </p>
				      *
				      * The vectors fields are not
				      * physically reasonable but the
				      * tangentiality constraint is
				      * clearly enforced. The fact
				      * that the vector fields are
				      * zero at some points on the
				      * boundary is an artifact of the
				      * way it is created, it is not
				      * constrained to be zero at
				      * these points.
				      *
				      * @ingroup constraints
				      */
  template <int dim, template <int, int> class DH, int spacedim>
  static
  void
  compute_no_normal_flux_constraints (const DH<dim,spacedim>         &dof_handler,
				      const unsigned int     first_vector_component,
				      const std::set<unsigned char> &boundary_ids,
				      ConstraintMatrix      &constraints,
				      const Mapping<dim, spacedim>    &mapping = StaticMappingQ1<dim>::mapping);


				     //@}
				     /**
				      * @name Assembling of right hand sides
				      */
				     //@{

				     /**
				      * Create a right hand side
				      * vector. Prior content of the
				      * given @p rhs_vector vector is
				      * deleted.
				      *
				      * See the general documentation of this
				      * class for further information.
				      */
    template <int dim, int spacedim>
    static void create_right_hand_side (const Mapping<dim, spacedim>    &mapping,
					const DoFHandler<dim,spacedim> &dof,
					const Quadrature<dim> &q,
					const Function<spacedim>   &rhs,
					Vector<double>        &rhs_vector);

				     /**
				      * Calls the create_right_hand_side()
				      * function, see above, with
				      * <tt>mapping=MappingQ1@<dim@>()</tt>.
				      */
    template <int dim, int spacedim>
    static void create_right_hand_side (const DoFHandler<dim,spacedim> &dof,
					const Quadrature<dim> &q,
					const Function<spacedim>   &rhs,
					Vector<double>        &rhs_vector);

				     /**
				      * Like the previous set of functions,
				      * but for hp objects.
				      */
    template <int dim, int spacedim>
    static void create_right_hand_side (const hp::MappingCollection<dim,spacedim>    &mapping,
					const hp::DoFHandler<dim,spacedim> &dof,
					const hp::QCollection<dim> &q,
					const Function<spacedim>   &rhs,
					Vector<double>        &rhs_vector);

				     /**
				      * Like the previous set of functions,
				      * but for hp objects.
				      */
    template <int dim, int spacedim>
    static void create_right_hand_side (const hp::DoFHandler<dim,spacedim> &dof,
					const hp::QCollection<dim> &q,
					const Function<spacedim>   &rhs,
					Vector<double>        &rhs_vector);

				     /**
				      * Create a right hand side
				      * vector for a point source at point @p p.
                                      * Prior content of the
				      * given @p rhs_vector vector is
				      * deleted.
				      *
				      * See the general documentation of this
				      * class for further information.
				      */
    template <int dim, int spacedim>
    static void create_point_source_vector(const Mapping<dim,spacedim>    &mapping,
                                           const DoFHandler<dim,spacedim> &dof,
                                           const Point<spacedim>      &p,
                                           Vector<double>        &rhs_vector);

				     /**
				      * Calls the create_point_source_vector()
				      * function, see above, with
				      * <tt>mapping=MappingQ1@<dim@>()</tt>.
				      */
    template <int dim, int spacedim>
    static void create_point_source_vector(const DoFHandler<dim,spacedim> &dof,
                                           const Point<spacedim>      &p,
                                           Vector<double>        &rhs_vector);

				     /**
				      * Like the previous set of functions,
				      * but for hp objects.
				      */
    template <int dim, int spacedim>
    static void create_point_source_vector(const hp::MappingCollection<dim,spacedim>    &mapping,
                                           const hp::DoFHandler<dim,spacedim> &dof,
                                           const Point<spacedim>      &p,
                                           Vector<double>        &rhs_vector);

				     /**
				      * Like the previous set of functions,
				      * but for hp objects. The function uses
				      * the default Q1 mapping object. Note
				      * that if your hp::DoFHandler uses any
				      * active fe index other than zero, then
				      * you need to call the function above
				      * that provides a mapping object for
				      * each active fe index.
				      */
    template <int dim, int spacedim>
    static void create_point_source_vector(const hp::DoFHandler<dim,spacedim> &dof,
                                           const Point<spacedim>      &p,
                                           Vector<double>        &rhs_vector);

                                     /**
				      * Create a right hand side
				      * vector from boundary
				      * forces. Prior content of the
				      * given @p rhs_vector vector is
				      * deleted.
				      *
				      * See the general documentation of this
				      * class for further information.
				      */
    template <int dim, int spacedim>
    static void create_boundary_right_hand_side (const Mapping<dim,spacedim>      &mapping,
						 const DoFHandler<dim,spacedim>   &dof,
						 const Quadrature<dim-1> &q,
						 const Function<spacedim>     &rhs,
						 Vector<double>          &rhs_vector,
						 const std::set<unsigned char> &boundary_indicators = std::set<unsigned char>());

				     /**
				      * Calls the
				      * create_boundary_right_hand_side()
				      * function, see above, with
				      * <tt>mapping=MappingQ1@<dim@>()</tt>.
				      */
    template <int dim, int spacedim>
    static void create_boundary_right_hand_side (const DoFHandler<dim,spacedim>   &dof,
						 const Quadrature<dim-1> &q,
						 const Function<spacedim>     &rhs,
						 Vector<double>          &rhs_vector,
						 const std::set<unsigned char> &boundary_indicators = std::set<unsigned char>());

                                     /**
				      * Same as the set of functions above,
				      * but for hp objects.
				      */
    template <int dim, int spacedim>
    static void create_boundary_right_hand_side (const hp::MappingCollection<dim,spacedim>      &mapping,
						 const hp::DoFHandler<dim,spacedim>   &dof,
						 const hp::QCollection<dim-1> &q,
						 const Function<spacedim>     &rhs,
						 Vector<double>          &rhs_vector,
						 const std::set<unsigned char> &boundary_indicators = std::set<unsigned char>());

				     /**
				      * Calls the
				      * create_boundary_right_hand_side()
				      * function, see above, with a
				      * single Q1 mapping as
				      * collection. This function
				      * therefore will only work if
				      * the only active fe index in
				      * use is zero.
				      */
    template <int dim, int spacedim>
    static void create_boundary_right_hand_side (const hp::DoFHandler<dim,spacedim>   &dof,
						 const hp::QCollection<dim-1> &q,
						 const Function<spacedim>     &rhs,
						 Vector<double>          &rhs_vector,
						 const std::set<unsigned char> &boundary_indicators = std::set<unsigned char>());

				     //@}
				     /**
				      * @name Evaluation of functions
				      * and errors
				      */
				     //@{

    				     /**
				      * Compute the error of the
				      * finite element solution.
				      * Integrate the difference
				      * between a reference function
				      * which is given as a continuous
				      * function object, and a finite
				      * element function.
				      *
				      * The value of @p exponent is
				      * used for computing $L^p$-norms
				      * and $W^{1,p}$-norms.
				      *
				      * The additional argument @p
				      * weight allows to evaluate
				      * weighted norms.  The weight
				      * function may be scalar,
				      * establishing a weight in the
				      * domain for all components
				      * equally. This may be used, for
				      * instance, to only integrates
				      * over parts of the domain.
				      *
				      * The weight function may also
				      * be vector-valued, with as many
				      * components as the finite
				      * element function: Then,
				      * different components get
				      * different weights. A typical
				      * application is when the error
				      * with respect to only one or a
				      * subset of the solution
				      * variables is to be computed,
				      * in which the other components
				      * would have weight values equal
				      * to zero. The
				      * ComponentSelectFunction class
				      * is particularly useful for
				      * this purpose.
				      *
				      * The weight function is
				      * expected to be positive, but
				      * negative values are not
				      * filtered. By default, no
				      * weighting function is given,
				      * i.e. weight=1 in the whole
				      * domain for all vector
				      * components uniformly.
				      *
				      * It is assumed that the number
				      * of components of the function
				      * @p exact_solution matches that
				      * of the finite element used by
				      * @p dof.
				      *
				      * See the general documentation of this
				      * class for more information.
				      *
				      * @note Instantiations for this template
				      * are provided for some vector types
				      * (see the general documentation of the
				      * class), but only for InVectors as in
				      * the documentation of the class,
				      * OutVector only Vector<double> and
				      * Vector<float>.
				      */
    template <int dim, class InVector, class OutVector, int spacedim>
    static void integrate_difference (const Mapping<dim,spacedim>    &mapping,
				      const DoFHandler<dim,spacedim> &dof,
				      const InVector        &fe_function,
				      const Function<spacedim>   &exact_solution,
				      OutVector             &difference,
				      const Quadrature<dim> &q,
				      const NormType        &norm,
				      const Function<spacedim>   *weight=0,
				      const double exponent = 2.);

				     /**
				      * Calls the integrate_difference()
				      * function, see above, with
				      * <tt>mapping=MappingQ1@<dim@>()</tt>.
				      */
    template <int dim, class InVector, class OutVector, int spacedim>
    static void integrate_difference (const DoFHandler<dim,spacedim> &dof,
				      const InVector        &fe_function,
				      const Function<spacedim>   &exact_solution,
				      OutVector             &difference,
				      const Quadrature<dim> &q,
				      const NormType        &norm,
				      const Function<spacedim>   *weight=0,
				      const double exponent = 2.);

    template <int dim, class InVector, class OutVector, int spacedim>
    static void integrate_difference (const hp::MappingCollection<dim,spacedim>    &mapping,
				      const hp::DoFHandler<dim,spacedim> &dof,
				      const InVector        &fe_function,
				      const Function<spacedim>   &exact_solution,
				      OutVector             &difference,
				      const hp::QCollection<dim> &q,
				      const NormType        &norm,
				      const Function<spacedim>   *weight=0,
				      const double exponent = 2.);

				     /**
				      * Calls the integrate_difference()
				      * function, see above, with
				      * <tt>mapping=MappingQ1@<dim@>()</tt>.
				      */
    template <int dim, class InVector, class OutVector, int spacedim>
    static void integrate_difference (const hp::DoFHandler<dim,spacedim> &dof,
				      const InVector        &fe_function,
				      const Function<spacedim>   &exact_solution,
				      OutVector             &difference,
				      const hp::QCollection<dim> &q,
				      const NormType        &norm,
				      const Function<spacedim>   *weight=0,
				      const double exponent = 2.);

				     /**
				      * Point error evaluation. Find
				      * the first cell containing the
				      * given point and compute the
				      * difference of a (possibly
				      * vector-valued) finite element
				      * function and a continuous
				      * function (with as many vector
				      * components as the finite
				      * element) at this point.
				      *
				      * This is a wrapper function
                                      * using a Q1-mapping for cell
                                      * boundaries to call the other
                                      * point_difference() function.
				      */
    template <int dim, class InVector, int spacedim>
    static void point_difference (const DoFHandler<dim,spacedim>& dof,
				  const InVector&        fe_function,
				  const Function<spacedim>&   exact_solution,
				  Vector<double>&        difference,
				  const Point<spacedim>&      point);

				     /**
				      * Point error evaluation. Find
				      * the first cell containing the
				      * given point and compute the
				      * difference of a (possibly
				      * vector-valued) finite element
				      * function and a continuous
				      * function (with as many vector
				      * components as the finite
				      * element) at this point.
				      *
                                      * Compared with the other
                                      * function of the same name,
                                      * this function uses an
                                      * arbitrary mapping to evaluate
                                      * the difference.
				      */
    template <int dim, class InVector, int spacedim>
    static void point_difference (const Mapping<dim, spacedim>    &mapping,
                                  const DoFHandler<dim,spacedim>& dof,
				  const InVector&        fe_function,
				  const Function<spacedim>&   exact_solution,
				  Vector<double>&        difference,
				  const Point<spacedim>&      point);

                                     /**
				      * Evaluate a possibly
				      * vector-valued finite element
				      * function defined by the given
				      * DoFHandler and nodal vector at
				      * the given point, and return
				      * the (vector) value of this
				      * function through the last
				      * argument.
                                      *
                                      * This is a wrapper function
                                      * using a Q1-mapping for cell
                                      * boundaries to call the other
                                      * point_difference() function.
				      */
    template <int dim, class InVector, int spacedim>
    static
    void
    point_value (const DoFHandler<dim,spacedim> &dof,
		 const InVector        &fe_function,
		 const Point<spacedim>      &point,
		 Vector<double>        &value);

				     /**
				      * Evaluate a scalar finite
				      * element function defined by
				      * the given DoFHandler and nodal
				      * vector at the given point, and
				      * return the value of this
				      * function.
                                      *
                                      * Compared with the other
                                      * function of the same name,
                                      * this is a wrapper function using
                                      * a Q1-mapping for cells.
				      *
				      * This function is used in the
				      * "Possibilities for extensions" part of
				      * the results section of @ref step_3
				      * "step-3".
				      */
    template <int dim, class InVector, int spacedim>
    static
    double
    point_value (const DoFHandler<dim,spacedim> &dof,
		 const InVector        &fe_function,
		 const Point<spacedim>      &point);

				     /**
				      * Evaluate a possibly
				      * vector-valued finite element
				      * function defined by the given
				      * DoFHandler and nodal vector at
				      * the given point, and return
				      * the (vector) value of this
				      * function through the last
				      * argument.
                                      *
                                      * Compared with the other
                                      * function of the same name,
                                      * this function uses an arbitrary
                                      * mapping to evaluate the difference.
				      */
    template <int dim, class InVector, int spacedim>
    static
    void
    point_value (const Mapping<dim, spacedim>    &mapping,
                 const DoFHandler<dim,spacedim> &dof,
		 const InVector        &fe_function,
		 const Point<spacedim>      &point,
		 Vector<double>        &value);

				     /**
				      * Evaluate a scalar finite
				      * element function defined by
				      * the given DoFHandler and nodal
				      * vector at the given point, and
				      * return the value of this
				      * function.
                                      *
                                      * Compared with the other
                                      * function of the same name,
                                      * this function uses an arbitrary
                                      * mapping to evaluate the difference.
				      */
    template <int dim, class InVector, int spacedim>
    static
    double
    point_value (const Mapping<dim,spacedim>    &mapping,
                 const DoFHandler<dim,spacedim> &dof,
		 const InVector        &fe_function,
		 const Point<spacedim>      &point);

				     //@}
				     /**
				      * Mean value operations
				      */
				     //@{

                                     /**
				      * Subtract the (algebraic) mean value
				      * from a vector. This function is most
				      * frequently used as a mean-value filter
				      * for Stokes: The pressure in Stokes'
				      * equations with only Dirichlet
				      * boundaries for the velocities is only
				      * determined up to a constant. This
				      * function allows to subtract the mean
				      * value of the pressure. It is usually
				      * called in a preconditioner and
				      * generates updates with mean value
				      * zero. The mean value is computed as
				      * the mean value of the degrees of
				      * freedom values as given by the input
				      * vector; they are not weighted by the
				      * area of cells, i.e. the mean is
				      * computed as $\sum_i v_i$, rather than
				      * as $\int_\Omega v(x) = \int_\Omega
				      * \sum_i v_i \phi_i(x)$. The latter can
				      * be obtained from the
				      * VectorTools::compute_mean_function,
				      * however.
				      *
				      * Apart from the vector @p v to operate
				      * on, this function takes a boolean mask
				      * that has a true entry for
				      * every component for which the mean
				      * value shall be computed and later
				      * subtracted. The argument is used to
				      * denote which components of the
				      * solution vector correspond to the
				      * pressure, and avoid touching all other
				      * components of the vector, such as the
				      * velocity components.
				      *
				      * @note In the context of using this
				      * function to filter out the kernel of
				      * an operator (such as the null space of
				      * the Stokes operator that consists of
				      * the constant pressures), this function
				      * only makes sense for finite elements
				      * for which the null space indeed
				      * consists of the vector
				      * $(1,1,\ldots,1)^T$. This is the case
				      * for example for the usual Lagrange
				      * elements where the sum of all shape
				      * functions equals the function that is
				      * constant one. However, it is not true
				      * for some other functions: for example,
				      * for the FE_DGP element (another valid
				      * choice for the pressure in Stokes
				      * discretizations), the first shape
				      * function on each cell is constant
				      * while further elements are $L_2$
				      * orthogonal to it (on the reference
				      * cell); consequently, the sum of all
				      * shape functions is not equal to one,
				      * and the vector that is associated with
				      * the constant mode is not equal to
				      * $(1,1,\ldots,1)^T$. For such elements,
				      * a different procedure has to be used
				      * when subtracting the mean value.
				      */
    static void subtract_mean_value(Vector<double>          &v,
				    const std::vector<bool> &p_select);

				     /**
				      * Compute the mean value of one
				      * component of the solution.
				      *
				      * This function integrates the
				      * chosen component over the
				      * whole domain and returns the
				      * result, i.e. it computes
				      * $\int_\Omega [u_h(x)]_c \; dx$
				      * where $c$ is the vector component
				      * and $u_h$ is the function
				      * representation of the nodal
				      * vector given as fourth
				      * argument. The integral is evaluated
				      * numerically using the quadrature
				      * formula given as third argument.
				      *
				      * This function is used in the
				      * "Possibilities for extensions" part of
				      * the results section of @ref step_3
				      * "step-3".
				      *
				      * @note The function is most often used
				      * when solving a problem whose solution
				      * is only defined up to a constant, for
				      * example a pure Neumann problem or the
				      * pressure in a Stokes or Navier-Stokes
				      * problem. In both cases, subtracting
				      * the mean value as computed by the
				      * current function, from the nodal
				      * vector does not generally yield the
				      * desired result of a finite element
				      * function with mean value zero. In
				      * fact, it only works for Lagrangian
				      * elements. For all other elements, you
				      * will need to compute the mean value
				      * and subtract it right inside the
				      * evaluation routine.
				      */
    template <int dim, class InVector, int spacedim>
    static double compute_mean_value (const Mapping<dim, spacedim>    &mapping,
				      const DoFHandler<dim,spacedim> &dof,
				      const Quadrature<dim> &quadrature,
				      const InVector        &v,
				      const unsigned int     component);

				     /**
				      * Calls the other compute_mean_value()
				      * function, see above, with
				      * <tt>mapping=MappingQ1@<dim@>()</tt>.
				      */
    template <int dim, class InVector, int spacedim>
    static double compute_mean_value (const DoFHandler<dim,spacedim> &dof,
				      const Quadrature<dim> &quadrature,
				      const InVector        &v,
				      const unsigned int     component);
				     //@}

				     /**
				      * Exception
				      */
    DeclException0 (ExcInvalidBoundaryIndicator);
				     /**
				      * Exception
				      */
    DeclException0 (ExcNonInterpolatingFE);
                                     /**
                                      * Exception
                                      */
    DeclException0 (ExcNoComponentSelected);
};


DEAL_II_NAMESPACE_CLOSE

#endif
