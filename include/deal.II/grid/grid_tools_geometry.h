// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2023 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_grid_tools_geometry_h
#define dealii_grid_tools_geometry_h

#include <deal.II/base/config.h>

#include <deal.II/base/bounding_box.h>
#include <deal.II/base/point.h>
#include <deal.II/base/template_constraints.h>
#include <deal.II/base/tensor.h>

#include <deal.II/fe/mapping.h>

#include <deal.II/grid/manifold.h>
#include <deal.II/grid/tria.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <numeric>
#include <utility>
#include <vector>

DEAL_II_NAMESPACE_OPEN

namespace GridTools
{
  /**
   * @name Information about meshes and cells
   */
  /** @{ */

  /**
   * Return the diameter of a triangulation. The diameter is computed using
   * only the vertices, i.e. if the diameter should be larger than the maximal
   * distance between boundary vertices due to a higher order mapping, then
   * this function will not catch this.
   */
  template <int dim, int spacedim>
  double
  diameter(const Triangulation<dim, spacedim> &tria);

  /**
   * Compute the volume (i.e. the dim-dimensional measure) of the
   * triangulation. We compute the measure using the integral $\sum_K \int_K 1
   * \; dx$ where $K$ are the cells of the given triangulation. The integral
   * is approximated via quadrature. This version of the function uses a
   * linear mapping to compute the JxW values on each cell.
   *
   * If the triangulation is a dim-dimensional one embedded in a higher
   * dimensional space of dimension spacedim, then the value returned is the
   * dim-dimensional measure. For example, for a two-dimensional triangulation
   * in three-dimensional space, the value returned is the area of the surface
   * so described. (This obviously makes sense since the spacedim-dimensional
   * measure of a dim-dimensional triangulation would always be zero if dim @<
   * spacedim).
   *
   * This function also works for objects of type
   * parallel::distributed::Triangulation, in which case the function is a
   * @ref GlossCollectiveOperation "collective operation".
   *
   * @param tria The triangulation.
   * @return The dim-dimensional measure of the domain described by the
   * triangulation, as discussed above.
   */
  template <int dim, int spacedim>
  double
  volume(const Triangulation<dim, spacedim> &tria);

  /**
   * Compute the volume (i.e. the dim-dimensional measure) of the
   * triangulation. We compute the measure using the integral $\sum_K \int_K 1
   * \; dx$ where $K$ are the cells of the given triangulation. The integral
   * is approximated via quadrature for which we use the mapping argument.
   *
   * If the triangulation is a dim-dimensional one embedded in a higher
   * dimensional space of dimension spacedim, then the value returned is the
   * dim-dimensional measure. For example, for a two-dimensional triangulation
   * in three-dimensional space, the value returned is the area of the surface
   * so described. (This obviously makes sense since the spacedim-dimensional
   * measure of a dim-dimensional triangulation would always be zero if dim @<
   * spacedim.
   *
   * This function also works for objects of type
   * parallel::distributed::Triangulation, in which case the function is a
   * @ref GlossCollectiveOperation "collective operation".
   *
   * @param tria The triangulation.
   * @param mapping The Mapping which computes the Jacobians used to
   * approximate the volume via quadrature. Explicitly using a higher-order
   * Mapping (i.e., instead of using the other version of this function) will
   * result in a more accurate approximation of the volume on Triangulations
   * with curvature described by Manifold objects.
   * @return The dim-dimensional measure of the domain described by the
   * triangulation, as discussed above.
   */
  template <int dim, int spacedim>
  double
  volume(const Triangulation<dim, spacedim> &tria,
         const Mapping<dim, spacedim>       &mapping);

  /**
   * Return an approximation of the diameter of the smallest active cell of a
   * triangulation. See step-24 for an example of use of this function.
   *
   * Notice that, even if you pass a non-trivial mapping, the returned value is
   * computed only using information on the vertices of the triangulation,
   * possibly transformed by the mapping. While this is accurate most of the
   * times, it may fail to give the correct result when the triangulation
   * contains very distorted cells.
   */
  template <int dim, int spacedim>
  double
  minimal_cell_diameter(
    const Triangulation<dim, spacedim> &triangulation,
    const Mapping<dim, spacedim>       &mapping =
      (ReferenceCells::get_hypercube<dim>()
#ifndef _MSC_VER
         .template get_default_linear_mapping<dim, spacedim>()
#else
         .ReferenceCell::get_default_linear_mapping<dim, spacedim>()
#endif
         ));

  /**
   * Return an approximation of the diameter of the largest active cell of a
   * triangulation.
   *
   * Notice that, even if you pass a non-trivial mapping to this function, the
   * returned value is computed only using information on the vertices of the
   * triangulation, possibly transformed by the mapping. While this is accurate
   * most of the times, it may fail to give the correct result when the
   * triangulation contains very distorted cells.
   */
  template <int dim, int spacedim>
  double
  maximal_cell_diameter(
    const Triangulation<dim, spacedim> &triangulation,
    const Mapping<dim, spacedim>       &mapping =
      (ReferenceCells::get_hypercube<dim>()
#ifndef _MSC_VER
         .template get_default_linear_mapping<dim, spacedim>()
#else
         .ReferenceCell::get_default_linear_mapping<dim, spacedim>()
#endif
         ));

  /**
   * Given a list of vertices (typically obtained using
   * Triangulation::get_vertices()) as the first, and a list of vertex indices
   * that characterize a single cell as the second argument, return the
   * measure (area, volume) of this cell. If this is a real cell, then you can
   * get the same result using <code>cell-@>measure()</code>, but this
   * function also works for cells that do not exist except that you make it
   * up by naming its vertices from the list.
   *
   * The size of @p vertex_indices, combined with `dim`, implicitly encodes
   * the ReferenceCell type of the provided cell. For example, if `dim == 2` and
   * `vertex_indices.size() == 3` then the cell is a triangle, but if
   * `dim == 2` and `vertex_indices.size() == 4` then the cell is a
   * quadrilateral. A std::vector is implicitly convertible to an ArrayView, so
   * it can be passed directly to this function. See the ArrayView class for
   * more information.
   *
   * @note This function is only implemented for codimension zero objects.
   */
  template <int dim>
  double
  cell_measure(const std::vector<Point<dim>>       &all_vertices,
               const ArrayView<const unsigned int> &vertex_indices);

  /**
   * Return the highest value among ratios between extents in each of the
   * coordinate directions of a @p cell. Moreover, return the dimension
   * relative to the highest elongation.
   *
   * @param[in] cell an iterator pointing to the cell.
   *
   * @return  A std::pair<unsigned int, double> such that the @p first value
   * is the dimension of the highest elongation and the @p second value is the
   * ratio among the dimensions of the @p cell.
   */
  template <int dim, int spacedim>
  std::pair<unsigned int, double>
  get_longest_direction(
    typename Triangulation<dim, spacedim>::active_cell_iterator cell);

  /**
   * This function computes an affine approximation of the map from the unit
   * coordinates to the real coordinates of the form $p_\text{real} = A
   * p_\text{unit} + b $ by a least squares fit of this affine function to the
   * $2^\text{dim}$ vertices representing a quadrilateral or hexahedral cell
   * in `spacedim` dimensions. The result is returned as a pair with the
   * matrix <i>A</i> as the first argument and the vector <i>b</i> describing
   * distance of the plane to the origin.
   *
   * For any valid mesh cell whose geometry is not degenerate, this operation
   * results in a unique affine mapping, even in cases where the actual
   * transformation by a bi-/trilinear or higher order mapping might be
   * singular. The result is exact in case the transformation from the unit to
   * the real cell is indeed affine, such as in one dimension or for Cartesian
   * and affine (parallelogram) meshes in 2d/3d.
   *
   * This approximation is underlying the function
   * TriaAccessor::real_to_unit_cell_affine_approximation() function.
   *
   * For exact transformations to the unit cell, use
   * Mapping::transform_real_to_unit_cell().
   */
  template <int dim, int spacedim>
  std::pair<DerivativeForm<1, dim, spacedim>, Tensor<1, spacedim>>
  affine_cell_approximation(const ArrayView<const Point<spacedim>> &vertices);

  /**
   * Computes an aspect ratio measure for all locally-owned active cells and
   * fills a vector with one entry per cell, given a @p triangulation and
   * @p mapping. The size of the vector that is returned equals the number of
   * active cells. The vector contains zero for non locally-owned cells. The
   * aspect ratio of a cell is defined as the ratio of the maximum to minimum
   * singular value of the Jacobian, taking the maximum over all quadrature
   * points of a quadrature rule specified via @p quadrature. For example, for
   * the special case of rectangular elements in 2d with dimensions $a$ and $b$
   * ($a \geq b$), this function returns the usual aspect ratio definition
   * $a/b$. The above definition using singular values is a generalization to
   * arbitrarily deformed elements. This function is intended to be used for
   * $d=2,3$ space dimensions, but it can also be used for $d=1$ returning a
   * value of 1.
   *
   * @note Inverted elements do not throw an exception. Instead, a value of inf
   * is written into the vector in case of inverted elements.
   *
   * @note Make sure to use enough quadrature points for a precise calculation
   * of the aspect ratio in case of deformed elements.
   *
   * @note In parallel computations the return value will have the length
   * n_active_cells but the aspect ratio is only computed for the cells that
   * are locally owned and placed at index CellAccessor::active_cell_index(),
   * respectively. All other values are set to 0.
   *
   * @note This function can only be used if deal.II was configured with
   * support for LAPACK.
   */
  template <int dim>
  Vector<double>
  compute_aspect_ratio_of_cells(const Mapping<dim>       &mapping,
                                const Triangulation<dim> &triangulation,
                                const Quadrature<dim>    &quadrature);

  /**
   * Computes the maximum aspect ratio by taking the maximum over all cells.
   *
   * @note When running in parallel with a Triangulation that supports MPI,
   * this is a collective call and the return value is the maximum over all
   * processors.
   */
  template <int dim>
  double
  compute_maximum_aspect_ratio(const Mapping<dim>       &mapping,
                               const Triangulation<dim> &triangulation,
                               const Quadrature<dim>    &quadrature);

  /**
   * Compute the smallest box containing the entire triangulation.
   *
   * If the input triangulation is a `parallel::distributed::Triangulation`,
   * then each processor will compute a bounding box enclosing all locally
   * owned, ghost, and artificial cells. In the case of a domain without curved
   * boundaries, these bounding boxes will all agree between processors because
   * the union of the areas occupied by artificial and ghost cells equals the
   * union of the areas occupied by the cells that other processors own.
   * However, if the domain has curved boundaries, this is no longer the case.
   * The bounding box returned may be appropriate for the current processor,
   * but different from the bounding boxes computed on other processors.
   */
  template <int dim, int spacedim>
  BoundingBox<spacedim>
  compute_bounding_box(const Triangulation<dim, spacedim> &triangulation);

  /**
   * Compute and return a bounding box, defined through a pair of points
   * bottom left and top right, that surrounds a subdomain of the @p mesh.
   * Here, the "subdomain" consists of exactly all of those
   * active cells for which the @p predicate returns @p true.
   *
   * For a description of how @p predicate works,
   * see compute_active_cell_halo_layer().
   *
   * @note This function was written before the BoundingBox class was invented.
   *   Consequently, it returns a pair of points, rather than a BoundingBox
   * object as one may expect. However, BoundingBox has a conversion constructor
   * from pairs of points, so the result of this function can still be assigned
   * to a BoundingBox object.
   *
   * @dealiiConceptRequires{concepts::is_triangulation_or_dof_handler<MeshType>}
   */
  template <typename MeshType>
  DEAL_II_CXX20_REQUIRES(concepts::is_triangulation_or_dof_handler<MeshType>)
  std::pair<
    Point<MeshType::space_dimension>,
    Point<MeshType::
            space_dimension>> compute_bounding_box(const MeshType &mesh,
                                                   const std::function<bool(
                                                     const typename MeshType::
                                                       active_cell_iterator &)>
                                                     &predicate);

  /**
   * Return the point on the geometrical object @p object closest to the given
   * point @p trial_point. For example, if @p object is a one-dimensional line
   * or edge, then the returned point will be a point on the geodesic that
   * connects the vertices as the manifold associated with the object sees it
   * (i.e., the geometric line may be curved if it lives in a higher
   * dimensional space). If the iterator points to a quadrilateral in a higher
   * dimensional space, then the returned point lies within the convex hull of
   * the vertices of the quad as seen by the associated manifold.
   *
   * @note This projection is usually not well-posed since there may be
   * multiple points on the object that minimize the distance. The algorithm
   * used in this function is robust (and the output is guaranteed to be on
   * the given @p object) but may only provide a few correct digits if the
   * object has high curvature. If your manifold supports it then the
   * specialized function Manifold::project_to_manifold() may perform better.
   */
  template <typename Iterator>
  Point<Iterator::AccessorType::space_dimension>
  project_to_object(
    const Iterator                                       &object,
    const Point<Iterator::AccessorType::space_dimension> &trial_point);
  /** @} */
} // namespace GridTools

#ifndef DOXYGEN
namespace GridTools
{
  namespace internal
  {
    namespace ProjectToObject
    {
      /**
       * The method GridTools::project_to_object requires taking derivatives
       * along the surface of a simplex. In general these cannot be
       * approximated with finite differences but special differences of the
       * form
       *
       *     df/dx_i - df/dx_j
       *
       * <em>can</em> be approximated. This <code>struct</code> just stores
       * the two derivatives approximated by the stencil (in the case of the
       * example above <code>i</code> and <code>j</code>).
       */
      struct CrossDerivative
      {
        const unsigned int direction_0;
        const unsigned int direction_1;

        CrossDerivative(const unsigned int d0, const unsigned int d1);
      };

      inline CrossDerivative::CrossDerivative(const unsigned int d0,
                                              const unsigned int d1)
        : direction_0(d0)
        , direction_1(d1)
      {}



      /**
       * Standard second-order approximation to the first derivative with a
       * two-point centered scheme. This is used below in a 1d Newton method.
       */
      template <typename F>
      inline auto
      centered_first_difference(const double center,
                                const double step,
                                const F &f) -> decltype(f(center) - f(center))
      {
        return (f(center + step) - f(center - step)) / (2.0 * step);
      }



      /**
       * Standard second-order approximation to the second derivative with a
       * three-point centered scheme. This is used below in a 1d Newton method.
       */
      template <typename F>
      inline auto
      centered_second_difference(const double center,
                                 const double step,
                                 const F &f) -> decltype(f(center) - f(center))
      {
        return (f(center + step) - 2.0 * f(center) + f(center - step)) /
               (step * step);
      }



      /**
       * Fourth order approximation of the derivative
       *
       *     df/dx_i - df/dx_j
       *
       * where <code>i</code> and <code>j</code> are specified by @p
       * cross_derivative. The derivative approximation is at @p center with a
       * step size of @p step and function @p f.
       */
      template <int structdim, typename F>
      inline auto
      cross_stencil(
        const CrossDerivative cross_derivative,
        const Tensor<1, GeometryInfo<structdim>::vertices_per_cell> &center,
        const double                                                 step,
        const F &f) -> decltype(f(center) - f(center))
      {
        Tensor<1, GeometryInfo<structdim>::vertices_per_cell> simplex_vector;
        simplex_vector[cross_derivative.direction_0] = 0.5 * step;
        simplex_vector[cross_derivative.direction_1] = -0.5 * step;
        return (-4.0 * f(center) - 1.0 * f(center + simplex_vector) -
                1.0 / 3.0 * f(center - simplex_vector) +
                16.0 / 3.0 * f(center + 0.5 * simplex_vector)) /
               step;
      }



      /**
       * The optimization algorithm used in GridTools::project_to_object is
       * essentially a gradient descent method. This function computes entries
       * in the gradient of the objective function; see the description in the
       * comments inside GridTools::project_to_object for more information.
       */
      template <int spacedim, int structdim, typename F>
      inline double
      gradient_entry(
        const unsigned int     row_n,
        const unsigned int     dependent_direction,
        const Point<spacedim> &p0,
        const Tensor<1, GeometryInfo<structdim>::vertices_per_cell> &center,
        const double                                                 step,
        const F                                                     &f)
      {
        Assert(row_n < GeometryInfo<structdim>::vertices_per_cell &&
                 dependent_direction <
                   GeometryInfo<structdim>::vertices_per_cell,
               ExcMessage("This function assumes that the last weight is a "
                          "dependent variable (and hence we cannot take its "
                          "derivative directly)."));
        Assert(row_n != dependent_direction,
               ExcMessage(
                 "We cannot differentiate with respect to the variable "
                 "that is assumed to be dependent."));

        const Point<spacedim>     manifold_point = f(center);
        const Tensor<1, spacedim> stencil_value  = cross_stencil<structdim>(
          {row_n, dependent_direction}, center, step, f);
        double entry = 0.0;
        for (unsigned int dim_n = 0; dim_n < spacedim; ++dim_n)
          entry +=
            -2.0 * (p0[dim_n] - manifold_point[dim_n]) * stencil_value[dim_n];
        return entry;
      }

      /**
       * Project onto a d-linear object. This is more accurate than the
       * general algorithm in project_to_object but only works for geometries
       * described by linear, bilinear, or trilinear mappings.
       */
      template <typename Iterator, int spacedim, int structdim>
      Point<spacedim>
      project_to_d_linear_object(const Iterator        &object,
                                 const Point<spacedim> &trial_point)
      {
        // let's look at this for simplicity for a quadrilateral
        // (structdim==2) in a space with spacedim>2 (notate trial_point by
        // y): all points on the surface are given by
        //   x(\xi) = sum_i v_i phi_x(\xi)
        // where v_i are the vertices of the quadrilateral, and
        // \xi=(\xi_1,\xi_2) are the reference coordinates of the
        // quadrilateral. so what we are trying to do is find a point x on the
        // surface that is closest to the point y. there are different ways to
        // solve this problem, but in the end it's a nonlinear problem and we
        // have to find reference coordinates \xi so that J(\xi) = 1/2 ||
        // x(\xi)-y ||^2 is minimal. x(\xi) is a function that is
        // structdim-linear in \xi, so J(\xi) is a polynomial of degree
        // 2*structdim that we'd like to minimize. unless structdim==1, we'll
        // have to use a Newton method to find the answer. This leads to the
        // following formulation of Newton steps:
        //
        // Given \xi_k, find \delta\xi_k so that
        //   H_k \delta\xi_k = - F_k
        // where H_k is an approximation to the second derivatives of J at
        // \xi_k, and F_k is the first derivative of J.  We'll iterate this a
        // number of times until the right hand side is small enough. As a
        // stopping criterion, we terminate if ||\delta\xi||<eps.
        //
        // As for the Hessian, the best choice would be
        //   H_k = J''(\xi_k)
        // but we'll opt for the simpler Gauss-Newton form
        //   H_k = A^T A
        // i.e.
        //   (H_k)_{nm} = \sum_{i,j} v_i*v_j *
        //                   \partial_n phi_i *
        //                   \partial_m phi_j
        // we start at xi=(0.5, 0.5).
        Point<structdim> xi;
        for (unsigned int d = 0; d < structdim; ++d)
          xi[d] = 0.5;

        Point<spacedim> x_k;
        for (const unsigned int i : GeometryInfo<structdim>::vertex_indices())
          x_k += object->vertex(i) *
                 GeometryInfo<structdim>::d_linear_shape_function(xi, i);

        do
          {
            Tensor<1, structdim> F_k;
            for (const unsigned int i :
                 GeometryInfo<structdim>::vertex_indices())
              F_k +=
                (x_k - trial_point) * object->vertex(i) *
                GeometryInfo<structdim>::d_linear_shape_function_gradient(xi,
                                                                          i);

            Tensor<2, structdim> H_k;
            for (const unsigned int i :
                 GeometryInfo<structdim>::vertex_indices())
              for (const unsigned int j :
                   GeometryInfo<structdim>::vertex_indices())
                {
                  Tensor<2, structdim> tmp = outer_product(
                    GeometryInfo<structdim>::d_linear_shape_function_gradient(
                      xi, i),
                    GeometryInfo<structdim>::d_linear_shape_function_gradient(
                      xi, j));
                  H_k += (object->vertex(i) * object->vertex(j)) * tmp;
                }

            const Tensor<1, structdim> delta_xi = -invert(H_k) * F_k;
            xi += delta_xi;

            x_k = Point<spacedim>();
            for (const unsigned int i :
                 GeometryInfo<structdim>::vertex_indices())
              x_k += object->vertex(i) *
                     GeometryInfo<structdim>::d_linear_shape_function(xi, i);

            if (delta_xi.norm() < 1e-7)
              break;
          }
        while (true);

        return x_k;
      }
    } // namespace ProjectToObject

    // We hit an internal compiler error in ICC 15 if we define this as a lambda
    // inside the project_to_object function below.
    template <int structdim>
    inline bool
    weights_are_ok(
      const Tensor<1, GeometryInfo<structdim>::vertices_per_cell> &v)
    {
      // clang has trouble figuring out structdim here, so define it
      // again:
      static const std::size_t n_vertices_per_cell =
        Tensor<1, GeometryInfo<structdim>::vertices_per_cell>::
          n_independent_components;
      std::array<double, n_vertices_per_cell> copied_weights;
      for (unsigned int i = 0; i < n_vertices_per_cell; ++i)
        {
          copied_weights[i] = v[i];
          if (v[i] < 0.0 || v[i] > 1.0)
            return false;
        }

      // check the sum: try to avoid some roundoff errors by summing in order
      std::sort(copied_weights.begin(), copied_weights.end());
      const double sum =
        std::accumulate(copied_weights.begin(), copied_weights.end(), 0.0);
      return std::abs(sum - 1.0) < 1e-10; // same tolerance used in manifold.cc
    }
  } // namespace internal



  template <typename Iterator>
  Point<Iterator::AccessorType::space_dimension>
  project_to_object(
    const Iterator                                       &object,
    const Point<Iterator::AccessorType::space_dimension> &trial_point)
  {
    const int spacedim  = Iterator::AccessorType::space_dimension;
    const int structdim = Iterator::AccessorType::structure_dimension;

    Point<spacedim> projected_point = trial_point;

    if (structdim >= spacedim)
      return projected_point;
    else if (structdim == 1 || structdim == 2)
      {
        using namespace internal::ProjectToObject;
        // Try to use the special flat algorithm for quads (this is better
        // than the general algorithm in 3d). This does not take into account
        // whether projected_point is outside the quad, but we optimize along
        // lines below anyway:
        const int                      dim = Iterator::AccessorType::dimension;
        const Manifold<dim, spacedim> &manifold = object->get_manifold();
        if (structdim == 2 && dynamic_cast<const FlatManifold<dim, spacedim> *>(
                                &manifold) != nullptr)
          {
            projected_point =
              project_to_d_linear_object<Iterator, spacedim, structdim>(
                object, trial_point);
          }
        else
          {
            // We want to find a point on the convex hull (defined by the
            // vertices of the object and the manifold description) that is
            // relatively close to the trial point. This has a few issues:
            //
            // 1. For a general convex hull we are not guaranteed that a unique
            //    minimum exists.
            // 2. The independent variables in the optimization process are the
            //    weights given to Manifold::get_new_point, which must sum to 1,
            //    so we cannot use standard finite differences to approximate a
            //    gradient.
            //
            // There is not much we can do about 1., but for 2. we can derive
            // finite difference stencils that work on a structdim-dimensional
            // simplex and rewrite the optimization problem to use those
            // instead. Consider the structdim 2 case and let
            //
            // F(c0, c1, c2, c3) = Manifold::get_new_point(vertices, {c0, c1,
            // c2, c3})
            //
            // where {c0, c1, c2, c3} are the weights for the four vertices on
            // the quadrilateral. We seek to minimize the Euclidean distance
            // between F(...) and trial_point. We can solve for c3 in terms of
            // the other weights and get, for one coordinate direction
            //
            // d/dc0 ((x0 - F(c0, c1, c2, 1 - c0 - c1 - c2))^2)
            //      = -2(x0 - F(...)) (d/dc0 F(...) - d/dc3 F(...))
            //
            // where we substitute back in for c3 after taking the
            // derivative. We can compute a stencil for the cross derivative
            // d/dc0 - d/dc3: this is exactly what cross_stencil approximates
            // (and gradient_entry computes the sum over the independent
            // variables). Below, we somewhat arbitrarily pick the last
            // component as the dependent one.
            //
            // Since we can now calculate derivatives of the objective
            // function we can use gradient descent to minimize it.
            //
            // Of course, this is much simpler in the structdim = 1 case (we
            // could rewrite the projection as a 1d optimization problem), but
            // to reduce the potential for bugs we use the same code in both
            // cases.
            const double step_size = object->diameter() / 64.0;

            constexpr unsigned int n_vertices_per_cell =
              GeometryInfo<structdim>::vertices_per_cell;

            std::array<Point<spacedim>, n_vertices_per_cell> vertices;
            for (unsigned int vertex_n = 0; vertex_n < n_vertices_per_cell;
                 ++vertex_n)
              vertices[vertex_n] = object->vertex(vertex_n);

            auto get_point_from_weights =
              [&](const Tensor<1, n_vertices_per_cell> &weights)
              -> Point<spacedim> {
              return object->get_manifold().get_new_point(
                make_array_view(vertices.begin(), vertices.end()),
                make_array_view(weights.begin_raw(), weights.end_raw()));
            };

            // pick the initial weights as (normalized) inverse distances from
            // the trial point:
            Tensor<1, n_vertices_per_cell> guess_weights;
            double                         guess_weights_sum = 0.0;
            for (unsigned int vertex_n = 0; vertex_n < n_vertices_per_cell;
                 ++vertex_n)
              {
                const double distance =
                  vertices[vertex_n].distance(trial_point);
                if (distance == 0.0)
                  {
                    guess_weights           = 0.0;
                    guess_weights[vertex_n] = 1.0;
                    guess_weights_sum       = 1.0;
                    break;
                  }
                else
                  {
                    guess_weights[vertex_n] = 1.0 / distance;
                    guess_weights_sum += guess_weights[vertex_n];
                  }
              }
            guess_weights /= guess_weights_sum;
            Assert(internal::weights_are_ok<structdim>(guess_weights),
                   ExcInternalError());

            // The optimization algorithm consists of two parts:
            //
            // 1. An outer loop where we apply the gradient descent algorithm.
            // 2. An inner loop where we do a line search to find the optimal
            //    length of the step one should take in the gradient direction.
            //
            for (unsigned int outer_n = 0; outer_n < 40; ++outer_n)
              {
                const unsigned int dependent_direction =
                  n_vertices_per_cell - 1;
                Tensor<1, n_vertices_per_cell> current_gradient;
                for (unsigned int row_n = 0; row_n < n_vertices_per_cell;
                     ++row_n)
                  {
                    if (row_n != dependent_direction)
                      {
                        current_gradient[row_n] =
                          gradient_entry<spacedim, structdim>(
                            row_n,
                            dependent_direction,
                            trial_point,
                            guess_weights,
                            step_size,
                            get_point_from_weights);

                        current_gradient[dependent_direction] -=
                          current_gradient[row_n];
                      }
                  }

                // We need to travel in the -gradient direction, as noted
                // above, but we may not want to take a full step in that
                // direction; instead, guess that we will go -0.5*gradient and
                // do quasi-Newton iteration to pick the best multiplier. The
                // goal is to find a scalar alpha such that
                //
                // F(x - alpha g)
                //
                // is minimized, where g is the gradient and F is the
                // objective function. To find the optimal value we find roots
                // of the derivative of the objective function with respect to
                // alpha by Newton iteration, where we approximate the first
                // and second derivatives of F(x - alpha g) with centered
                // finite differences.
                double gradient_weight = -0.5;
                auto   gradient_weight_objective_function =
                  [&](const double gradient_weight_guess) -> double {
                  return (trial_point -
                          get_point_from_weights(guess_weights +
                                                 gradient_weight_guess *
                                                   current_gradient))
                    .norm_square();
                };

                for (unsigned int inner_n = 0; inner_n < 10; ++inner_n)
                  {
                    const double update_numerator = centered_first_difference(
                      gradient_weight,
                      step_size,
                      gradient_weight_objective_function);
                    const double update_denominator =
                      centered_second_difference(
                        gradient_weight,
                        step_size,
                        gradient_weight_objective_function);

                    // avoid division by zero. Note that we limit the gradient
                    // weight below
                    if (std::abs(update_denominator) == 0.0)
                      break;
                    gradient_weight =
                      gradient_weight - update_numerator / update_denominator;

                    // Put a fairly lenient bound on the largest possible
                    // gradient (things tend to be locally flat, so the gradient
                    // itself is usually small)
                    if (std::abs(gradient_weight) > 10)
                      {
                        gradient_weight = -10.0;
                        break;
                      }
                  }

                // It only makes sense to take convex combinations with weights
                // between zero and one. If the update takes us outside of this
                // region then rescale the update to stay within the region and
                // try again
                Tensor<1, n_vertices_per_cell> tentative_weights =
                  guess_weights + gradient_weight * current_gradient;

                double new_gradient_weight = gradient_weight;
                for (unsigned int iteration_count = 0; iteration_count < 40;
                     ++iteration_count)
                  {
                    if (internal::weights_are_ok<structdim>(tentative_weights))
                      break;

                    for (unsigned int i = 0; i < n_vertices_per_cell; ++i)
                      {
                        if (tentative_weights[i] < 0.0)
                          {
                            tentative_weights -=
                              (tentative_weights[i] / current_gradient[i]) *
                              current_gradient;
                          }
                        if (tentative_weights[i] < 0.0 ||
                            1.0 < tentative_weights[i])
                          {
                            new_gradient_weight /= 2.0;
                            tentative_weights =
                              guess_weights +
                              new_gradient_weight * current_gradient;
                          }
                      }
                  }

                // the update might still send us outside the valid region, so
                // check again and quit if the update is still not valid
                if (!internal::weights_are_ok<structdim>(tentative_weights))
                  break;

                // if we cannot get closer by traveling in the gradient
                // direction then quit
                if (get_point_from_weights(tentative_weights)
                      .distance(trial_point) <
                    get_point_from_weights(guess_weights).distance(trial_point))
                  guess_weights = tentative_weights;
                else
                  break;
                Assert(internal::weights_are_ok<structdim>(guess_weights),
                       ExcInternalError());
              }
            Assert(internal::weights_are_ok<structdim>(guess_weights),
                   ExcInternalError());
            projected_point = get_point_from_weights(guess_weights);
          }

        // if structdim == 2 and the optimal point is not on the interior then
        // we may be able to get a more accurate result by projecting onto the
        // lines.
        if (structdim == 2)
          {
            std::array<Point<spacedim>, GeometryInfo<structdim>::lines_per_cell>
              line_projections;
            for (unsigned int line_n = 0;
                 line_n < GeometryInfo<structdim>::lines_per_cell;
                 ++line_n)
              {
                line_projections[line_n] =
                  project_to_object(object->line(line_n), trial_point);
              }
            std::sort(line_projections.begin(),
                      line_projections.end(),
                      [&](const Point<spacedim> &a, const Point<spacedim> &b) {
                        return a.distance(trial_point) <
                               b.distance(trial_point);
                      });
            if (line_projections[0].distance(trial_point) <
                projected_point.distance(trial_point))
              projected_point = line_projections[0];
          }
      }
    else
      {
        DEAL_II_NOT_IMPLEMENTED();
        return projected_point;
      }

    return projected_point;
  }
} // namespace GridTools
#endif // DOXYGEN

DEAL_II_NAMESPACE_CLOSE

#endif
