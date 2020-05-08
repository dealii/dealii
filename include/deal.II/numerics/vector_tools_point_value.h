// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------

#ifndef dealii_vector_tools_point_value_h
#define dealii_vector_tools_point_value_h


#include <deal.II/base/config.h>

DEAL_II_NAMESPACE_OPEN

template <int dim, int spacedim>
class DoFHandler;
template <int dim, typename Number>
class Function;
template <int dim, int spacedim>
class Mapping;
template <int dim, typename Number>
class Point;
template <typename Number>
class Vector;
namespace hp
{
  template <int dim, int spacedim>
  class DoFHandler;
  template <int dim, int spacedim>
  class MappingCollection;
} // namespace hp

namespace VectorTools
{
  /**
   * @name Assembling of right hand sides
   */
  //@{

  /**
   * Create a right hand side vector for a point source at point @p p. In
   * other words, it creates a vector $F$ so that $F_i = \int_\Omega
   * \delta(x-p) \varphi_i(x) dx$ where $\varphi_i$ are the shape functions
   * described by @p dof_handler and @p p is the point at which the delta
   * function is located. Prior content of the given @p rhs_vector
   * vector is deleted. This function is for the case of a scalar finite
   * element.
   *
   * This function is typically used in one of these two contexts:
   * - Let's say you want to solve the same kind of problems many times
   *   over, with different values for right hand sides or coefficients,
   *   and then evaluate the solution at the same point every time. You
   *   could do this by calling VectorTools::point_value() after each
   *   solve, or you could realize that to evaluate the solution $u_h$
   *   at a point $p$, you could rearrange operations like this:
   *   @f{align*}{
   *     u_h(p) &= \sum_j U_j \varphi_j(p) = \sum_j U_j F_j
   *       \\   &= U \cdot F
   *   @f}
   *   with the vector as defined above. In other words, point evaluation
   *   can be achieved with just a single vector-vector product, and the
   *   vector $F$ can be computed once and for all and reused
   *   for each solve, without having to go through the mesh every time
   *   to find out which cell (and where in the cell) the point $p$ is
   *   located.
   * - This function is also useful if you wanted to compute the Green's
   *   function for the problem you are solving. This is because the
   *   Green's function $G(x,p)$ is defined by
   *   @f{align*}{
   *     L G(x,p) &= \delta(x-p)
   *   @f}
   *   where $L$ is the differential operator of your problem. The discrete
   *   version then requires computing the right hand side vector
   *   $F_i = \int_\Omega \varphi_i(x) \delta(x-p)$, which is exactly
   *   the vector computed by the current function.
   *
   * While maybe not relevant for documenting <i>what</i> this
   * function does, it may be interesting to note that delta functions
   * do not exist in reality, and consequently, using this function
   * does not model any real situation. This is, because no real
   * object is able to focus an infinite force density at an
   * infinitesimally small part of the domain (rather, all real
   * devices will spread out the force over a finite area); nor is it
   * possible to measure values at individual points (but all
   * measurements will somehow be averaged over small areas). Only if
   * this area is so small that it cannot be resolved by any mesh does
   * it make sense to model the situation in a way that uses a delta
   * function with the same overall force or sensitivity. On the other
   * hand, a situation that is probably more fruitfully simulated with
   * a delta function is the electric potential of a point source; in
   * this case, the solution is known to have a logarithmic
   * singularity (in 2d) or a $\frac{1}{r}$ singularity (in 3d),
   * neither of which is bounded.
   *
   * Mathematically, the use of delta functions typically leads to exact
   * solutions to which the numerically obtained, approximate solution does
   * not converge. This is because, taking the Laplace equation as an example,
   * the error between exact and numerical solution can be bounded by the
   * expression
   * @f{align*}{
   *   \| u-u_h \|_{L_2} \le C h \| \nabla u \|_{L_2}
   * @f}
   * but when using a delta function on the right hand side, the term
   * $\| \nabla u \|_{L_2} = |u|_{H^1}$ is not finite. This can be seen
   * by using the a-priori bound for solutions of the Laplace equation
   * $-\Delta u = f$ that states that $|u|_{H^1} \le \|f\|_{H^{-1}}$.
   * When using a delta function as right hand side, $f(x)=\delta(x-p)$,
   * one would need to take the $H^{-1}$ norm of a delta function, which
   * however is not finite because $\delta(\cdot-p) \not\in H^{-1}$.
   *
   * The consequence of all of this is that the exact solution of the
   * Laplace equation with a delta function on the right hand side --
   * i.e., the <i>Green's function</i> -- has a singularity at $p$ that
   * is so strong that it cannot be resolved by a finite element
   * solution, and consequently finite element approximations do not
   * converge towards the exact solution in any of the usual norms.
   *
   * All of this is also the case for all of the other usual second-order
   * partial differential equations in dimensions two or higher. (Because
   * in dimension two and higher, $H^1$ functions are not necessarily
   * continuous, and consequently the delta function is not in the dual
   * space $H^{-1}$.)
   */
  template <int dim, int spacedim>
  void
  create_point_source_vector(const Mapping<dim, spacedim> &   mapping,
                             const DoFHandler<dim, spacedim> &dof_handler,
                             const Point<spacedim, double> &  p,
                             Vector<double> &                 rhs_vector);

  /**
   * Call the create_point_source_vector() function, see above, with
   * an implied default $Q_1$ mapping object.
   */
  template <int dim, int spacedim>
  void
  create_point_source_vector(const DoFHandler<dim, spacedim> &dof_handler,
                             const Point<spacedim, double> &  p,
                             Vector<double> &                 rhs_vector);

  /**
   * Like the previous set of functions, but for hp objects.
   */
  template <int dim, int spacedim>
  void
  create_point_source_vector(
    const hp::MappingCollection<dim, spacedim> &mapping,
    const hp::DoFHandler<dim, spacedim> &       dof_handler,
    const Point<spacedim, double> &             p,
    Vector<double> &                            rhs_vector);

  /**
   * Like the previous set of functions, but for hp objects. The function uses
   * an implied default $Q_1$ mapping object. Note that if your hp::DoFHandler
   * uses any active fe index other than zero, then you need to call the
   * function above that provides a mapping object for each active fe index.
   */
  template <int dim, int spacedim>
  void
  create_point_source_vector(const hp::DoFHandler<dim, spacedim> &dof_handler,
                             const Point<spacedim, double> &      p,
                             Vector<double> &                     rhs_vector);

  /**
   * Create a right hand side vector for a point source at point @p p. This
   * variation of the function is meant for vector-valued problems with
   * exactly dim components (it will also work for problems with more than dim
   * components, and in this case simply consider only the first dim
   * components of the shape functions). It computes a right hand side that
   * corresponds to a forcing function that is equal to a delta function times
   * a given direction. In other words, it creates a vector $F$ so that $F_i =
   * \int_\Omega [\mathbf d \delta(x-p)] \cdot \varphi_i(x) dx$. Note here that
   * $\varphi_i$ is a vector-valued function. $\mathbf d$ is the given direction
   * of the source term $\mathbf d \delta(x-p)$ and corresponds to the @p
   * direction argument to be passed to this function.
   *
   * Prior content of the given @p rhs_vector vector is deleted.
   *
   * See the discussion of the first create_point_source_vector() variant for
   * more on the use of delta functions.
   */
  template <int dim, int spacedim>
  void
  create_point_source_vector(const Mapping<dim, spacedim> &   mapping,
                             const DoFHandler<dim, spacedim> &dof_handler,
                             const Point<spacedim, double> &  p,
                             const Point<dim, double> &       direction,
                             Vector<double> &                 rhs_vector);

  /**
   * Call the create_point_source_vector() function for vector-valued finite
   * elements, see above, with an implied default $Q_1$ mapping object.
   */
  template <int dim, int spacedim>
  void
  create_point_source_vector(const DoFHandler<dim, spacedim> &dof_handler,
                             const Point<spacedim, double> &  p,
                             const Point<dim, double> &       direction,
                             Vector<double> &                 rhs_vector);

  /**
   * Like the previous set of functions, but for hp objects.
   */
  template <int dim, int spacedim>
  void
  create_point_source_vector(
    const hp::MappingCollection<dim, spacedim> &mapping,
    const hp::DoFHandler<dim, spacedim> &       dof_handler,
    const Point<spacedim, double> &             p,
    const Point<dim, double> &                  direction,
    Vector<double> &                            rhs_vector);

  /**
   * Like the previous set of functions, but for hp objects. The function uses
   * an implied default $Q_1$ mapping object. Note that if your hp::DoFHandler
   * uses any active fe index other than zero, then you need to call the
   * function above that provides a mapping object for each active fe index.
   */
  template <int dim, int spacedim>
  void
  create_point_source_vector(const hp::DoFHandler<dim, spacedim> &dof_handler,
                             const Point<spacedim, double> &      p,
                             const Point<dim, double> &           direction,
                             Vector<double> &                     rhs_vector);

  // @}

  /**
   * @name Evaluation of functions and errors
   */
  //@{

  /**
   * Point error evaluation. Find the first cell containing the given point
   * and compute the difference of a (possibly vector-valued) finite element
   * function and a continuous function (with as many vector components as the
   * finite element) at this point.
   *
   * This is a wrapper function using a Q1-mapping for cell boundaries to call
   * the other point_difference() function.
   *
   * @note If the cell in which the point is found is not locally owned, an
   * exception of type VectorTools::ExcPointNotAvailableHere is thrown.
   */
  template <int dim, typename VectorType, int spacedim>
  void
  point_difference(
    const DoFHandler<dim, spacedim> &                          dof,
    const VectorType &                                         fe_function,
    const Function<spacedim, typename VectorType::value_type> &exact_solution,
    Vector<typename VectorType::value_type> &                  difference,
    const Point<spacedim, double> &                            point);

  /**
   * Point error evaluation. Find the first cell containing the given point
   * and compute the difference of a (possibly vector-valued) finite element
   * function and a continuous function (with as many vector components as the
   * finite element) at this point.
   *
   * Compared with the other function of the same name, this function uses an
   * arbitrary mapping to evaluate the difference.
   *
   * @note If the cell in which the point is found is not locally owned, an
   * exception of type VectorTools::ExcPointNotAvailableHere is thrown.
   */
  template <int dim, typename VectorType, int spacedim>
  void
  point_difference(
    const Mapping<dim, spacedim> &                             mapping,
    const DoFHandler<dim, spacedim> &                          dof,
    const VectorType &                                         fe_function,
    const Function<spacedim, typename VectorType::value_type> &exact_solution,
    Vector<typename VectorType::value_type> &                  difference,
    const Point<spacedim, double> &                            point);

  /**
   * Evaluate a possibly vector-valued finite element function defined by the
   * given DoFHandler and nodal vector @p fe_function at the given point @p
   * point, and return the (vector) value of this function through the last
   * argument.
   *
   * This function uses a $Q_1$-mapping for the cell the point is evaluated
   * in. If you need to evaluate using a different mapping (for example when
   * using curved boundaries), use the point_difference() function that takes
   * a mapping.
   *
   * This function is not particularly cheap. This is because it first
   * needs to find which cell a given point is in, then find the point
   * on the reference cell that matches the given evaluation point,
   * and then evaluate the shape functions there. You probably do not
   * want to use this function to evaluate the solution at <i>many</i>
   * points. For this kind of application, the FEFieldFunction class
   * offers at least some optimizations. On the other hand, if you
   * want to evaluate <i>many solutions</i> at the same point, you may
   * want to look at the VectorTools::create_point_source_vector()
   * function.
   *
   * @note If the cell in which the point is found is not locally owned, an
   *   exception of type VectorTools::ExcPointNotAvailableHere is thrown.
   *
   * @note This function needs to find the cell within which a point lies,
   *   and this can only be done up to a certain numerical tolerance of course.
   *   Consequently, for points that are on, or close to, the boundary of
   *   a cell, you may get the value of the finite element field either
   *   here or there, depending on which cell the point is found in. This
   *   does not matter (to within the same tolerance) if the finite element
   *   field is continuous. On the other hand, if the finite element in use
   *   is <i>not</i> continuous, then you will get unpredictable values for
   *   points on or close to the boundary of the cell, as one would expect
   *   when trying to evaluate point values of discontinuous functions.
   */
  template <int dim, typename VectorType, int spacedim>
  void
  point_value(const DoFHandler<dim, spacedim> &        dof,
              const VectorType &                       fe_function,
              const Point<spacedim, double> &          point,
              Vector<typename VectorType::value_type> &value);

  /**
   * Same as above for hp.
   *
   * @note If the cell in which the point is found is not locally owned, an
   * exception of type VectorTools::ExcPointNotAvailableHere is thrown.
   *
   * @note This function needs to find the cell within which a point lies,
   *   and this can only be done up to a certain numerical tolerance of course.
   *   Consequently, for points that are on, or close to, the boundary of
   *   a cell, you may get the value of the finite element field either
   *   here or there, depending on which cell the point is found in. This
   *   does not matter (to within the same tolerance) if the finite element
   *   field is continuous. On the other hand, if the finite element in use
   *   is <i>not</i> continuous, then you will get unpredictable values for
   *   points on or close to the boundary of the cell, as one would expect
   *   when trying to evaluate point values of discontinuous functions.
   */
  template <int dim, typename VectorType, int spacedim>
  void
  point_value(const hp::DoFHandler<dim, spacedim> &    dof,
              const VectorType &                       fe_function,
              const Point<spacedim, double> &          point,
              Vector<typename VectorType::value_type> &value);

  /**
   * Evaluate a scalar finite element function defined by the given DoFHandler
   * and nodal vector @p fe_function at the given point @p point, and return
   * the value of this function.
   *
   * This function uses a Q1-mapping for the cell the point is evaluated
   * in. If you need to evaluate using a different mapping (for example when
   * using curved boundaries), use the point_difference() function that takes
   * a mapping.
   *
   * This function is not particularly cheap. This is because it first
   * needs to find which cell a given point is in, then find the point
   * on the reference cell that matches the given evaluation point,
   * and then evaluate the shape functions there. You probably do not
   * want to use this function to evaluate the solution at <i>many</i>
   * points. For this kind of application, the FEFieldFunction class
   * offers at least some optimizations. On the other hand, if you
   * want to evaluate <i>many solutions</i> at the same point, you may
   * want to look at the VectorTools::create_point_source_vector()
   * function.
   *
   * This function is used in the "Possibilities for extensions" part of the
   * results section of
   * @ref step_3 "step-3".
   *
   * @note If the cell in which the point is found is not locally owned, an
   * exception of type VectorTools::ExcPointNotAvailableHere is thrown.
   *
   * @note This function needs to find the cell within which a point lies,
   *   and this can only be done up to a certain numerical tolerance of course.
   *   Consequently, for points that are on, or close to, the boundary of
   *   a cell, you may get the value of the finite element field either
   *   here or there, depending on which cell the point is found in. This
   *   does not matter (to within the same tolerance) if the finite element
   *   field is continuous. On the other hand, if the finite element in use
   *   is <i>not</i> continuous, then you will get unpredictable values for
   *   points on or close to the boundary of the cell, as one would expect
   *   when trying to evaluate point values of discontinuous functions.
   */
  template <int dim, typename VectorType, int spacedim>
  typename VectorType::value_type
  point_value(const DoFHandler<dim, spacedim> &dof,
              const VectorType &               fe_function,
              const Point<spacedim, double> &  point);

  /**
   * Same as above for hp.
   *
   * @note If the cell in which the point is found is not locally owned, an
   * exception of type VectorTools::ExcPointNotAvailableHere is thrown.
   *
   * @note This function needs to find the cell within which a point lies,
   *   and this can only be done up to a certain numerical tolerance of course.
   *   Consequently, for points that are on, or close to, the boundary of
   *   a cell, you may get the value of the finite element field either
   *   here or there, depending on which cell the point is found in. This
   *   does not matter (to within the same tolerance) if the finite element
   *   field is continuous. On the other hand, if the finite element in use
   *   is <i>not</i> continuous, then you will get unpredictable values for
   *   points on or close to the boundary of the cell, as one would expect
   *   when trying to evaluate point values of discontinuous functions.
   */
  template <int dim, typename VectorType, int spacedim>
  typename VectorType::value_type
  point_value(const hp::DoFHandler<dim, spacedim> &dof,
              const VectorType &                   fe_function,
              const Point<spacedim, double> &      point);

  /**
   * Evaluate a possibly vector-valued finite element function defined by the
   * given DoFHandler and nodal vector @p fe_function at the given point @p
   * point, and return the (vector) value of this function through the last
   * argument.
   *
   * Compared with the other function of the same name, this function uses an
   * arbitrary mapping to evaluate the point value.
   *
   * This function is not particularly cheap. This is because it first
   * needs to find which cell a given point is in, then find the point
   * on the reference cell that matches the given evaluation point,
   * and then evaluate the shape functions there. You probably do not
   * want to use this function to evaluate the solution at <i>many</i>
   * points. For this kind of application, the FEFieldFunction class
   * offers at least some optimizations. On the other hand, if you
   * want to evaluate <i>many solutions</i> at the same point, you may
   * want to look at the VectorTools::create_point_source_vector()
   * function.
   *
   * @note If the cell in which the point is found is not locally owned, an
   * exception of type VectorTools::ExcPointNotAvailableHere is thrown.
   *
   * @note This function needs to find the cell within which a point lies,
   *   and this can only be done up to a certain numerical tolerance of course.
   *   Consequently, for points that are on, or close to, the boundary of
   *   a cell, you may get the value of the finite element field either
   *   here or there, depending on which cell the point is found in. This
   *   does not matter (to within the same tolerance) if the finite element
   *   field is continuous. On the other hand, if the finite element in use
   *   is <i>not</i> continuous, then you will get unpredictable values for
   *   points on or close to the boundary of the cell, as one would expect
   *   when trying to evaluate point values of discontinuous functions.
   */
  template <int dim, typename VectorType, int spacedim>
  void
  point_value(const Mapping<dim, spacedim> &           mapping,
              const DoFHandler<dim, spacedim> &        dof,
              const VectorType &                       fe_function,
              const Point<spacedim, double> &          point,
              Vector<typename VectorType::value_type> &value);

  /**
   * Same as above for hp.
   *
   * @note If the cell in which the point is found is not locally owned, an
   * exception of type VectorTools::ExcPointNotAvailableHere is thrown.
   *
   * @note This function needs to find the cell within which a point lies,
   *   and this can only be done up to a certain numerical tolerance of course.
   *   Consequently, for points that are on, or close to, the boundary of
   *   a cell, you may get the value of the finite element field either
   *   here or there, depending on which cell the point is found in. This
   *   does not matter (to within the same tolerance) if the finite element
   *   field is continuous. On the other hand, if the finite element in use
   *   is <i>not</i> continuous, then you will get unpredictable values for
   *   points on or close to the boundary of the cell, as one would expect
   *   when trying to evaluate point values of discontinuous functions.
   */
  template <int dim, typename VectorType, int spacedim>
  void
  point_value(const hp::MappingCollection<dim, spacedim> &mapping,
              const hp::DoFHandler<dim, spacedim> &       dof,
              const VectorType &                          fe_function,
              const Point<spacedim, double> &             point,
              Vector<typename VectorType::value_type> &   value);

  /**
   * Evaluate a scalar finite element function defined by the given DoFHandler
   * and nodal vector @p fe_function at the given point @p point, and return
   * the value of this function.
   *
   * Compared with the other function of the same name, this function uses an
   * arbitrary mapping to evaluate the difference.
   *
   * This function is not particularly cheap. This is because it first
   * needs to find which cell a given point is in, then find the point
   * on the reference cell that matches the given evaluation point,
   * and then evaluate the shape functions there. You probably do not
   * want to use this function to evaluate the solution at <i>many</i>
   * points. For this kind of application, the FEFieldFunction class
   * offers at least some optimizations. On the other hand, if you
   * want to evaluate <i>many solutions</i> at the same point, you may
   * want to look at the VectorTools::create_point_source_vector()
   * function.
   *
   * @note If the cell in which the point is found is not locally owned, an
   * exception of type VectorTools::ExcPointNotAvailableHere is thrown.
   *
   * @note This function needs to find the cell within which a point lies,
   *   and this can only be done up to a certain numerical tolerance of course.
   *   Consequently, for points that are on, or close to, the boundary of
   *   a cell, you may get the value of the finite element field either
   *   here or there, depending on which cell the point is found in. This
   *   does not matter (to within the same tolerance) if the finite element
   *   field is continuous. On the other hand, if the finite element in use
   *   is <i>not</i> continuous, then you will get unpredictable values for
   *   points on or close to the boundary of the cell, as one would expect
   *   when trying to evaluate point values of discontinuous functions.
   */
  template <int dim, typename VectorType, int spacedim>
  typename VectorType::value_type
  point_value(const Mapping<dim, spacedim> &   mapping,
              const DoFHandler<dim, spacedim> &dof,
              const VectorType &               fe_function,
              const Point<spacedim, double> &  point);

  /**
   * Same as above for hp.
   *
   * @note If the cell in which the point is found is not locally owned, an
   * exception of type VectorTools::ExcPointNotAvailableHere is thrown.
   *
   * @note This function needs to find the cell within which a point lies,
   *   and this can only be done up to a certain numerical tolerance of course.
   *   Consequently, for points that are on, or close to, the boundary of
   *   a cell, you may get the value of the finite element field either
   *   here or there, depending on which cell the point is found in. This
   *   does not matter (to within the same tolerance) if the finite element
   *   field is continuous. On the other hand, if the finite element in use
   *   is <i>not</i> continuous, then you will get unpredictable values for
   *   points on or close to the boundary of the cell, as one would expect
   *   when trying to evaluate point values of discontinuous functions.
   */
  template <int dim, typename VectorType, int spacedim>
  typename VectorType::value_type
  point_value(const hp::MappingCollection<dim, spacedim> &mapping,
              const hp::DoFHandler<dim, spacedim> &       dof,
              const VectorType &                          fe_function,
              const Point<spacedim, double> &             point);
  //@}
} // namespace VectorTools

DEAL_II_NAMESPACE_CLOSE

#endif // dealii_vector_tools_point_value_h
