// ---------------------------------------------------------------------
//
// Copyright (C) 2022 by the deal.II authors
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

#ifndef dealii_cgal_utilities_h
#define dealii_cgal_utilities_h

#include <deal.II/base/config.h>

#ifdef DEAL_II_WITH_CGAL
#  include <CGAL/Cartesian.h>
#  include <CGAL/Complex_2_in_triangulation_3.h>
#  include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#  include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#  include <CGAL/Kernel_traits.h>
#  include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#  include <CGAL/Mesh_criteria_3.h>
#  include <CGAL/Mesh_triangulation_3.h>
#  include <CGAL/Polygon_mesh_processing/corefinement.h>
#  include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#  include <CGAL/Polyhedral_mesh_domain_with_features_3.h>
#  include <CGAL/Simple_cartesian.h>
#  include <CGAL/Surface_mesh.h>
#  include <CGAL/Triangulation_3.h>
#  include <CGAL/make_mesh_3.h>
#  include <CGAL/make_surface_mesh.h>

#  include <type_traits>



DEAL_II_NAMESPACE_OPEN
/**
 * Interface to the Computational Geometry Algorithm Library (CGAL).
 *
 * CGAL is a software project that provides easy access to efficient and
 * reliable geometric algorithms. The library offers data structures and
 * algorithms like triangulations, Voronoi diagrams, Boolean operations on
 * polygons and polyhedra, point set processing, arrangements of curves, surface
 * and volume mesh generation, geometry processing,  alpha shapes, convex hull
 * algorithms, shape reconstruction, AABB and KD trees...
 *
 * You can learn more about the CGAL library at https://www.cgal.org/
 */
namespace CGALWrappers
{
#  ifdef CGAL_CONCURRENT_MESH_3
  using ConcurrencyTag = CGAL::Parallel_tag;
#  else
  using ConcurrencyTag = CGAL::Sequential_tag;
#  endif

  /**
   * An enum type given to functions that compute boolean operations between
   * geometrical objects, defined by triangulated surface grids.
   *
   * As an example, we show all supported boolean operations applied to a hedra
   * (the union of two pyramids) and a cube.
   *
   * @image html hedra_cube.png
   */
  enum class BooleanOperation
  {
    /**
     * Given two triangulated surfaces, refine the first surface until its
     * intersection with the second surface is captured exactly by the
     * refinement
     *
     * @image html corefinement.png
     */
    compute_corefinement = 1 << 0,

    /**
     * Given two triangulated surfaces, compute the boolean difference of the
     * first surface minus the second surface
     *
     * @image html boolean_difference.png
     */
    compute_difference = 1 << 1,

    /**
     * Given two triangulated surfaces, compute their intersection
     *
     * @image html boolean_intersection.png
     */
    compute_intersection = 1 << 2,

    /**
     * Given two triangulated surfaces, compute their union
     *
     * @image html boolean_union.png
     */
    compute_union = 1 << 3,
  };

  /**
   * Convert from a deal.II Point to any compatible CGAL point.
   *
   * @tparam CGALPointType Any of the CGAL point types
   * @tparam dim Dimension of the point
   * @param [in] p An input deal.II Point<dim>
   * @return CGALPointType A CGAL point
   */
  template <typename CGALPointType, int dim>
  inline CGALPointType
  dealii_point_to_cgal_point(const dealii::Point<dim> &p);

  /**
   * Convert from various CGAL point types to deal.II Point.
   *
   * @tparam dim Dimension of the point
   * @tparam CGALPointType Any of the CGAL point types
   * @param p An input CGAL point type
   * @return dealii::Point<dim> The corresponding deal.II point.
   */
  template <int dim, typename CGALPointType>
  inline dealii::Point<dim>
  cgal_point_to_dealii_point(const CGALPointType &p);

  /**
   * Given a closed CGAL::Surface_mesh, this function fills the
   * region bounded by the surface with tets, keeping them as coarse as
   * possible.
   *
   * The number of the generated tetrahedrons depends on the refinement level of
   * surface mesh. This function does not attempt to construct a fine
   * triangulation, nor to smooth the final result. If you need finer control on
   * the resulting triangulation, you should consider using directly
   * CGAL::make_mesh_3().
   *
   * @param [in] surface_mesh The (closed) surface mesh bounding the volume that has to be filled.
   * @param [out] triangulation The output triangulation filled with tetrahedra.
   */
  template <typename C3t3>
  void
  cgal_surface_mesh_to_cgal_coarse_triangulation(
    CGAL::Surface_mesh<typename C3t3::Point::Point> &surface_mesh,
    C3t3 &                                           triangulation);

  /**
   * Given two triangulated surface meshes that bound two volumes, execute a
   * boolean operation on them, and store the result in a third surface mesh.
   *
   * Quoting from CGAL documentation
   * (https://doc.cgal.org/latest/Polygon_mesh_processing/index.html#title14):
   *
   * > Given a closed triangulated surface mesh, each connected component splits
   * > the 3D space into two subspaces. The vertex sequence of each face of a
   * > component is seen either clockwise or counterclockwise from these two
   * > subspaces. The subspace that sees the sequence clockwise (resp.
   * > counterclockwise) is on the negative (resp. positive) side of the
   * > component.
   *
   * > Given a closed triangulated surface mesh `surface` with no
   * > self-intersections, we say that `surface` bounds a volume if each
   * > subspace lies exclusively on the positive (or negative) side of all the
   * > incident connected components of `surface`. The volume bounded by
   * > `surface` is the union of all subspaces that are on negative sides of
   * > their incident connected components of `surface`.
   *
   * > There is no restriction on the topology of the input volumes. However,
   * > there are some requirements on the input to guarantee that the operation
   * > is possible. First, the input meshes must not self-intersect. Second, the
   * > operation is possible only if the output can be bounded by a manifold
   * > triangulated surface mesh. In particular this means that the output
   * > volume has no part with zero thickness. Mathematically speaking, the
   * > intersection with an infinitesimally small ball centered in the output
   * > volume is a topological ball. At the surface level this means that no
   * > non-manifold vertex or edge is allowed in the output. For example, it is
   * > not possible to compute the union of two cubes that are disjoint but
   * > sharing an edge.
   *
   * See BooleanOperation for a list of available operations, with the
   * corresponding examples.
   *
   * @param[in] surface_mesh_1 The first surface mesh.
   * @param[in] surface_mesh_2 The second surface mesh.
   * @param[in] boolean_operation See BooleanOperation for the list of the
   * allowed operations.
   * @param[out] output_surface_mesh The surface mesh with the result of the
   * boolean operation.
   */
  template <typename CGALPointType>
  void
  compute_boolean_operation(
    const CGAL::Surface_mesh<CGALPointType> &surface_mesh_1,
    const CGAL::Surface_mesh<CGALPointType> &surface_mesh_2,
    const BooleanOperation &                 boolean_operation,
    CGAL::Surface_mesh<CGALPointType> &      output_surface_mesh);
} // namespace CGALWrappers

#  ifndef DOXYGEN
// Template implementations
namespace CGALWrappers
{
  template <typename CGALPointType, int dim>
  inline CGALPointType
  dealii_point_to_cgal_point(const dealii::Point<dim> &p)
  {
    constexpr int cdim = CGALPointType::Ambient_dimension::value;
    static_assert(dim <= cdim, "Only dim <= cdim supported");
    if constexpr (cdim == 1)
      return CGALPointType(p[0]);
    else if constexpr (cdim == 2)
      return CGALPointType(p[0], dim > 1 ? p[1] : 0);
    else if constexpr (cdim == 3)
      return CGALPointType(p[0], dim > 1 ? p[1] : 0, dim > 2 ? p[2] : 0);
    else
      Assert(false, dealii::ExcNotImplemented());
    return CGALPointType();
  }



  template <int dim, typename CGALPointType>
  inline dealii::Point<dim>
  cgal_point_to_dealii_point(const CGALPointType &p)
  {
    constexpr int cdim = CGALPointType::Ambient_dimension::value;
    if constexpr (dim == 1)
      return dealii::Point<dim>(CGAL::to_double(p.x()));
    else if constexpr (dim == 2)
      return dealii::Point<dim>(CGAL::to_double(p.x()),
                                cdim > 1 ? CGAL::to_double(p.y()) : 0);
    else if constexpr (dim == 3)
      return dealii::Point<dim>(CGAL::to_double(p.x()),
                                cdim > 1 ? CGAL::to_double(p.y()) : 0,
                                cdim > 2 ? CGAL::to_double(p.z()) : 0);
    else
      Assert(false, dealii::ExcNotImplemented());
  }



  template <typename C3t3>
  void
  cgal_surface_mesh_to_cgal_coarse_triangulation(
    CGAL::Surface_mesh<typename C3t3::Point::Point> &surface_mesh,
    C3t3 &                                           triangulation)
  {
    using CGALPointType = typename C3t3::Point::Point;
    Assert(CGAL::is_closed(surface_mesh),
           ExcMessage("The surface mesh must be closed."));

    using K           = typename CGAL::Kernel_traits<CGALPointType>::Kernel;
    using Mesh_domain = CGAL::Polyhedral_mesh_domain_with_features_3<
      K,
      CGAL::Surface_mesh<CGALPointType>>;
    using Tr = typename CGAL::
      Mesh_triangulation_3<Mesh_domain, CGAL::Default, ConcurrencyTag>::type;
    using Mesh_criteria = CGAL::Mesh_criteria_3<Tr>;

    CGAL::Polygon_mesh_processing::triangulate_faces(surface_mesh);
    Mesh_domain domain(surface_mesh);
    domain.detect_features();
    Mesh_criteria criteria;
    // Mesh generation
    triangulation = CGAL::make_mesh_3<C3t3>(domain,
                                            criteria,
                                            CGAL::parameters::no_perturb(),
                                            CGAL::parameters::no_exude());
  }



  template <typename CGALPointType>
  void
  compute_boolean_operation(
    const CGAL::Surface_mesh<CGALPointType> &surface_mesh_1,
    const CGAL::Surface_mesh<CGALPointType> &surface_mesh_2,
    const BooleanOperation &                 boolean_operation,
    CGAL::Surface_mesh<CGALPointType> &      output_surface_mesh)
  {
    Assert(output_surface_mesh.is_empty(),
           ExcMessage(
             "output_surface_mesh must be empty upon calling this function"));
    Assert(CGAL::is_closed(surface_mesh_1),
           ExcMessage(
             "The input surface_mesh_1 must be a closed surface mesh."));
    Assert(CGAL::is_closed(surface_mesh_2),
           ExcMessage(
             "The input surface_mesh_2 must be a closed surface mesh."));
    Assert(CGAL::is_triangle_mesh(surface_mesh_1),
           ExcMessage("The first CGAL mesh must be triangulated."));
    Assert(CGAL::is_triangle_mesh(surface_mesh_2),
           ExcMessage("The second CGAL mesh must be triangulated."));

    bool res      = false;
    auto surf_1   = surface_mesh_1;
    auto surf_2   = surface_mesh_2;
    namespace PMP = CGAL::Polygon_mesh_processing;
    switch (boolean_operation)
      {
        case BooleanOperation::compute_union:
          res = PMP::corefine_and_compute_union(surf_1,
                                                surf_2,
                                                output_surface_mesh);
          break;
        case BooleanOperation::compute_intersection:
          res = PMP::corefine_and_compute_intersection(surf_1,
                                                       surf_2,
                                                       output_surface_mesh);
          break;
        case BooleanOperation::compute_difference:
          res = PMP::corefine_and_compute_difference(surf_1,
                                                     surf_2,
                                                     output_surface_mesh);
          break;
        case BooleanOperation::compute_corefinement:
          PMP::corefine(surf_1,
                        surf_2); // both surfaces are corefined
          output_surface_mesh = std::move(surf_1);
          res                 = true;
          break;
        default:
          output_surface_mesh.clear();
          break;
      }
    Assert(res,
           ExcMessage("The boolean operation was not successfully computed."));
  }
} // namespace CGALWrappers
#  endif

DEAL_II_NAMESPACE_CLOSE

#endif
#endif
