// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 1998 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_vector_tools_h
#define dealii_vector_tools_h


#include <deal.II/base/config.h>

#include <deal.II/numerics/vector_tools_boundary.h>
#include <deal.II/numerics/vector_tools_common.h>
#include <deal.II/numerics/vector_tools_constraints.h>
#include <deal.II/numerics/vector_tools_evaluate.h>
#include <deal.II/numerics/vector_tools_integrate_difference.h>
#include <deal.II/numerics/vector_tools_interpolate.h>
#include <deal.II/numerics/vector_tools_mean_value.h>
#include <deal.II/numerics/vector_tools_point_gradient.h>
#include <deal.II/numerics/vector_tools_point_value.h>
#include <deal.II/numerics/vector_tools_project.h>
#include <deal.II/numerics/vector_tools_rhs.h>


DEAL_II_NAMESPACE_OPEN

// TODO: Move documentation of functions to the functions!

/**
 * Provide a namespace which offers some operations on vectors. Among these
 * are assembling of standard vectors, integration of the difference of a
 * finite element solution and a continuous function, interpolations and
 * projections of continuous functions to the finite element space and other
 * operations.
 *
 * @note There exist two versions of almost all functions, one that takes an
 * explicit Mapping argument and one that does not. The second one generally
 * calls the first with an implicit $Q_1$ argument (i.e., with an argument of
 * kind MappingQ(1)). If your intend your code to use a different
 * mapping than a (bi-/tri-)linear one, then you need to call the functions
 * <b>with</b> mapping argument should be used.
 *
 *
 * <h3>Description of operations</h3>
 *
 * This collection of methods offers the following operations:
 * <ul>
 * <li> Interpolation: assign each degree of freedom in the vector to be the
 * value of the function given as argument. This is identical to saying that
 * the resulting finite element function (which is isomorphic to the output
 * vector) has exact function values in all support points of trial functions.
 * The support point of a trial function is the point where its value equals
 * one, e.g. for linear trial functions the support points are four corners of
 * an element. This function therefore relies on the assumption that a finite
 * element is used for which the degrees of freedom are function values
 * (Lagrange elements) rather than gradients, normal derivatives, second
 * derivatives, etc (Hermite elements, quintic Argyris element, etc.).
 *
 * It seems inevitable that some values of the vector to be created are set
 * twice or even more than that. The reason is that we have to loop over all
 * cells and get the function values for each of the trial functions located
 * thereon. This applies also to the functions located on faces and corners
 * which we thus visit more than once. While setting the value in the vector
 * is not an expensive operation, the evaluation of the given function may be,
 * taking into account that a virtual function has to be called.
 *
 * <li> Projection: compute the <i>L</i><sup>2</sup>-projection of the given
 * function onto the finite element space, i.e. if <i>f</i> is the function to
 * be projected, compute <i>f<sub>h</sub></i> in <i>V<sub>h</sub></i> such
 * that
 * (<i>f<sub>h</sub></i>,<i>v<sub>h</sub></i>)=(<i>f</i>,<i>v<sub>h</sub></i>)
 * for all discrete test functions <i>v<sub>h</sub></i>. This is done through
 * the solution of the linear system of equations <i> M v = f</i> where
 * <i>M</i> is the @ref GlossMassMatrix "mass matrix" $m_{ij} = \int_\Omega \phi_i(x) \phi_j(x) dx$
 * and $f_i = \int_\Omega f(x) \phi_i(x) dx$. The solution vector $v$ then is
 * the nodal representation of the projection <i>f<sub>h</sub></i>. The
 * project() functions are used in the step-21 and step-23 tutorial programs.
 *
 * In order to get proper results, it be may necessary to treat boundary
 * conditions right. Below are listed some cases where this may be needed.  If
 * needed, this is done by <i>L</i><sup>2</sup>-projection of the trace of the
 * given function onto the finite element space restricted to the boundary of
 * the domain, then taking this information and using it to eliminate the
 * boundary nodes from the mass matrix of the whole domain, using the
 * MatrixTools::apply_boundary_values() function. The projection of the trace
 * of the function to the boundary is done with the
 * VectorTools::project_boundary_values() (see below) function, which is
 * called with a map of boundary functions
 * std::map<types::boundary_id, const Function<spacedim,number>*> in which all
 * boundary indicators from zero to numbers::internal_face_boundary_id-1
 * (numbers::internal_face_boundary_id is used for other purposes, see the
 * Triangulation class documentation) point to the function to be projected.
 * The projection to the boundary takes place using a second quadrature
 * formula on the boundary given to the project() function. The first
 * quadrature formula is used to compute the right hand side and for numerical
 * quadrature of the mass matrix.
 *
 * The projection of the boundary values first, then eliminating them from the
 * global system of equations is not needed usually. It may be necessary if
 * you want to enforce special restrictions on the boundary values of the
 * projected function, for example in time dependent problems: you may want to
 * project the initial values but need consistency with the boundary values
 * for later times. Since the latter are projected onto the boundary in each
 * time step, it is necessary that we also project the boundary values of the
 * initial values, before projecting them to the whole domain.
 *
 * Obviously, the results of the two schemes for projection are different.
 * Usually, when projecting to the boundary first, the
 * <i>L</i><sup>2</sup>-norm of the difference between original function and
 * projection over the whole domain will be larger (factors of five have been
 * observed) while the <i>L</i><sup>2</sup>-norm of the error integrated over
 * the boundary should of course be less. The reverse should also hold if no
 * projection to the boundary is performed.
 *
 * The selection whether the projection to the boundary first is needed is
 * done with the <tt>project_to_boundary_first</tt> flag passed to the
 * function.  If @p false is given, the additional quadrature formula for
 * faces is ignored.
 *
 * You should be aware of the fact that if no projection to the boundary is
 * requested, a function with zero boundary values may not have zero boundary
 * values after projection. There is a flag for this especially important
 * case, which tells the function to enforce zero boundary values on the
 * respective boundary parts. Since enforced zero boundary values could also
 * have been reached through projection, but are more economically obtain
 * using other methods, the @p project_to_boundary_first flag is ignored if
 * the @p enforce_zero_boundary flag is set.
 *
 * The solution of the linear system is presently done using a simple CG
 * method without preconditioning and without multigrid. This is clearly not
 * too efficient, but sufficient in many cases and simple to implement. This
 * detail may change in the future.
 *
 * <li> Creation of right hand side vectors: The create_right_hand_side()
 * function computes the vector $f_i = \int_\Omega f(x) \phi_i(x) dx$. This is
 * the same as what the <tt>MatrixCreator::create_*</tt> functions which take
 * a right hand side do, but without assembling a matrix.
 *
 * <li> Creation of right hand side vectors for point sources: The
 * create_point_source_vector() function computes the vector $F_i =
 * \int_\Omega \delta(x-x_0) \phi_i(x) dx$.
 *
 * <li> Creation of boundary right hand side vectors: The
 * create_boundary_right_hand_side() function computes the vector $f_i =
 * \int_{\partial\Omega} g(x) \phi_i(x) dx$. This is the right hand side
 * contribution of boundary forces when having inhomogeneous Neumann boundary
 * values in Laplace's equation or other second order operators. This function
 * also takes an optional argument denoting over which parts of the boundary
 * the integration shall extend. If the default argument is used, it is
 * applied to all boundaries.
 *
 * <li> Interpolation of boundary values: The
 * MatrixTools::apply_boundary_values() function takes a list of boundary
 * nodes and their values. You can get such a list by interpolation of a
 * boundary function using the interpolate_boundary_values() function. To use
 * it, you have to specify a list of pairs of boundary indicators (of type
 * <tt>types::boundary_id</tt>; see the section in the documentation of the
 * Triangulation class for more details) and the according functions denoting
 * the Dirichlet boundary values of the nodes on boundary faces with this
 * boundary indicator.
 *
 * Usually, all other boundary conditions, such as inhomogeneous Neumann
 * values or mixed boundary conditions are handled in the weak formulation. No
 * attempt is made to include these into the process of matrix and vector
 * assembly therefore.
 *
 * Within this function, boundary values are interpolated, i.e. a node is
 * given the point value of the boundary function. In some cases, it may be
 * necessary to use the L2-projection of the boundary function or any other
 * method. For this purpose we refer to the project_boundary_values() function
 * below.
 *
 * You should be aware that the boundary function may be evaluated at nodes on
 * the interior of faces. These, however, need not be on the true boundary,
 * but rather are on the approximation of the boundary represented by the
 * mapping of the unit cell to the real cell. Since this mapping will in most
 * cases not be the exact one at the face, the boundary function is evaluated
 * at points which are not on the boundary and you should make sure that the
 * returned values are reasonable in some sense anyway.
 *
 * In 1d the situation is a bit different since there faces (i.e. vertices)
 * have no boundary indicator. It is assumed that if the boundary indicator
 * zero is given in the list of boundary functions, the left boundary point is
 * to be interpolated while the right boundary point is associated with the
 * boundary index 1 in the map. The respective boundary functions are then
 * evaluated at the place of the respective boundary point.
 *
 * <li> Projection of boundary values: The project_boundary_values() function
 * acts similar to the interpolate_boundary_values() function, apart from the
 * fact that it does not get the nodal values of boundary nodes by
 * interpolation but rather through the <i>L</i><sup>2</sup>-projection of the
 * trace of the function to the boundary.
 *
 * The projection takes place on all boundary parts with boundary indicators
 * listed in the map (std::map<types::boundary_id, const
 * Function<spacedim,number>*>) of boundary functions. These boundary parts may
 * or may not be continuous. For these boundary parts, the @ref GlossMassMatrix "mass matrix" is
 * assembled using the MatrixTools::create_boundary_mass_matrix() function, as
 * well as the appropriate right hand side. Then the resulting system of
 * equations is solved using a simple CG method (without preconditioning), which
 * is in most cases sufficient for the present purpose.
 *
 * <li> Computing errors: The function integrate_difference() performs the
 * calculation of the error between a given (continuous) reference function
 * and the finite element solution in different norms. The integration is
 * performed using a given quadrature formula and assumes that the given
 * finite element objects equals that used for the computation of the
 * solution.
 *
 * The result is stored in a vector (named @p difference), where each entry
 * equals the given norm of the difference on a cell. The order of entries is
 * the same as a @p cell_iterator takes when started with @p begin_active and
 * promoted with the <tt>++</tt> operator.
 *
 * This data, one number per active cell, can be used to generate graphical
 * output by directly passing it to the DataOut class through the
 * DataOut::add_data_vector function. Alternatively, the global error can be
 * computed using VectorTools::compute_global_error(). Finally, the output per
 * cell from VectorTools::integrate_difference() can be interpolated to the
 * nodal points of a finite element field using the
 * DoFTools::distribute_cell_to_dof_vector function.
 *
 * Presently, there is the possibility to compute the following values from
 * the difference, on each cell: @p mean, @p L1_norm, @p L2_norm, @p
 * Linfty_norm, @p H1_seminorm and @p H1_norm, see VectorTools::NormType. For
 * the mean difference value, the reference function minus the numerical
 * solution is computed, not the other way round.
 *
 * The infinity norm of the difference on a given cell returns the maximum
 * absolute value of the difference at the quadrature points given by the
 * quadrature formula parameter. This will in some cases not be too good an
 * approximation, since for example the Gauss quadrature formulae do not
 * evaluate the difference at the end or corner points of the cells. You may
 * want to choose a quadrature formula with more quadrature points or one with
 * another distribution of the quadrature points in this case. You should also
 * take into account the superconvergence properties of finite elements in
 * some points: for example in 1d, the standard finite element method is a
 * collocation method and should return the exact value at nodal points.
 * Therefore, the trapezoidal rule should always return a vanishing L-infinity
 * error. Conversely, in 2d the maximum L-infinity error should be located at
 * the vertices or at the center of the cell, which would make it plausible to
 * use the Simpson quadrature rule. On the other hand, there may be
 * superconvergence at Gauss integration points. These examples are not
 * intended as a rule of thumb, rather they are thought to illustrate that the
 * use of the wrong quadrature formula may show a significantly wrong result
 * and care should be taken to chose the right formula.
 *
 * The <i>H</i><sup>1</sup> seminorm is the <i>L</i><sup>2</sup> norm of the
 * gradient of the difference. The square of the full <i>H</i><sup>1</sup>
 * norm is the sum of the square of seminorm and the square of the
 * <i>L</i><sup>2</sup> norm.
 *
 * To get the global <i>L<sup>1</sup></i> error, you have to sum up the
 * entries in @p difference, e.g. using Vector::l1_norm() function.  For the
 * global <i>L</i><sup>2</sup> difference, you have to sum up the squares of
 * the entries and take the root of the sum, e.g. using Vector::l2_norm().
 * These two operations represent the <i>l</i><sub>1</sub> and
 * <i>l</i><sub>2</sub> norms of the vectors, but you need not take the
 * absolute value of each entry, since the cellwise norms are already
 * positive.
 *
 * To get the global mean difference, simply sum up the elements as above. To
 * get the $L_\infty$ norm, take the maximum of the vector elements, e.g.
 * using the Vector::linfty_norm() function.
 *
 * For the global <i>H</i><sup>1</sup> norm and seminorm, the same rule
 * applies as for the <i>L</i><sup>2</sup> norm: compute the
 * <i>l</i><sub>2</sub> norm of the cell error vector.
 *
 * Note that, in the codimension one case, if you ask for a norm that requires
 * the computation of a gradient, then the provided function is automatically
 * projected along the curve, and the difference is only computed on the
 * tangential part of the gradient, since no information is available on the
 * normal component of the gradient anyway.
 * </ul>
 *
 * All functions use the finite element given to the DoFHandler object the
 * last time that the degrees of freedom were distributed over the
 * triangulation. Also, if access to an object describing the exact form of
 * the boundary is needed, the pointer stored within the triangulation object
 * is accessed.
 *
 * @note Instantiations for this template are provided for some vector types,
 * in particular <code>Vector&lt;float&gt;, Vector&lt;double&gt;,
 * BlockVector&lt;float&gt;, BlockVector&lt;double&gt;</code>; others can be
 * generated in application code (see the section on
 * @ref Instantiations
 * in the manual).
 *
 * @ingroup numerics
 */
namespace VectorTools
{}

DEAL_II_NAMESPACE_CLOSE

#endif
