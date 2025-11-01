// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2022 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_cgal_utilities_h
#define dealii_cgal_utilities_h

#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>

#include <deal.II/cgal/additional_data.h>

#ifdef DEAL_II_WITH_CGAL
#  include <deal.II/base/quadrature_lib.h>

#  include <deal.II/cgal/surface_mesh.h>

#  include <deal.II/grid/tria.h>

#  include <boost/hana.hpp>

#  include <CGAL/version.h>
#  if CGAL_VERSION_MAJOR >= 6
#    include <CGAL/Installation/internal/disable_deprecation_warnings_and_errors.h>
#  endif
#  include <CGAL/Cartesian.h>
#  include <CGAL/Complex_2_in_triangulation_3.h>
#  include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#  include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#  include <CGAL/Kernel_traits.h>
#  include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#  include <CGAL/Mesh_criteria_3.h>
#  include <CGAL/Mesh_triangulation_3.h>
// Disable a warning that we get with gcc-13 about a potential uninitialized
// usage of an <anonymous> lambda function in this external CGAL header.
DEAL_II_DISABLE_EXTRA_DIAGNOSTICS
#  include <CGAL/Polygon_mesh_processing/corefinement.h>
DEAL_II_ENABLE_EXTRA_DIAGNOSTICS
#  include <CGAL/Polygon_mesh_processing/measure.h>
#  include <CGAL/Polygon_mesh_processing/remesh.h>
#  include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#  include <CGAL/Polyhedral_mesh_domain_with_features_3.h>
#  include <CGAL/Simple_cartesian.h>
#  include <CGAL/Surface_mesh.h>
#  include <CGAL/Triangulation_3.h>
#  include <CGAL/boost/graph/copy_face_graph.h>
#  include <CGAL/boost/graph/selection.h>
#  include <CGAL/convex_hull_3.h>
#  include <CGAL/make_mesh_3.h>
#  include <CGAL/make_surface_mesh.h>

#  include <fstream>
#  include <limits>
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
   * Given a closed CGAL::Surface_mesh, this function fills the
   * region bounded by the surface with tets.
   *
   *
   * @param [in] surface_mesh The (closed) surface mesh bounding the volume that has to be filled.
   * @param [out] triangulation The output triangulation filled with tetrahedra.
   * @param [in] data AdditionalData object to pass to the CGAL::make_mesh_3 function. See the documentation
   * of that struct for a description of those parameters.
   */
  template <typename C3t3>
  void
  cgal_surface_mesh_to_cgal_triangulation(
    CGAL::Surface_mesh<typename C3t3::Point::Point> &surface_mesh,
    C3t3                                            &triangulation,
    const AdditionalData<3> &data = AdditionalData<3>{});

  /**
   * Given two triangulated surface meshes that bound two volumes, execute a
   * boolean operation on them, and store the result in a third surface mesh.
   *
   * Quoting from CGAL documentation
   * (https://doc.cgal.org/latest/Polygon_mesh_processing/index.html#title14):
   *
   * > Given a closed triangulated surface mesh, each connected component splits
   * > the 3d space into two subspaces. The vertex sequence of each face of a
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
    const BooleanOperation                  &boolean_operation,
    CGAL::Surface_mesh<CGALPointType>       &output_surface_mesh);

  /**
   * Given a CGAL Triangulation describing a polyhedral region, create
   * a Quadrature rule to integrate over the polygon by looping through all the
   * vertices and exploiting QGaussSimplex.
   *
   * @param[in] tria The CGAL triangulation object describing the polyhedral
   * region.
   * @param[in] degree Desired degree of the Quadrature rule on each element of
   * the polyhedral.
   * @return A global Quadrature rule that allows to integrate over the polyhedron.
   */
  template <typename CGALTriangulationType>
  dealii::Quadrature<CGALTriangulationType::Point::Ambient_dimension::value>
  compute_quadrature(const CGALTriangulationType &tria,
                     const unsigned int           degree);

  /**
   * Compute a Quadrature formula over the polygonal/polyhedral region described
   * by a BooleanOperation between two deal.II cell_iterator objects. The degree
   * of the Quadrature
   * rule is specified by the argument @p degree. This function splits the polygon/polyhedron
   * into simplices and collects QGaussSimplex quadrature rules.
   *
   * @param [in] cell0 A cell_iterator to the first deal.II cell.
   * @param [in] cell1 A cell_iterator to the second deal.II cell.
   * @param [in] degree The degree of accuracy you wish to get for the global quadrature formula.
   * @param [in] bool_op The BooleanOperation to be performed.
   * @param [in] mapping0 Mapping object for the first cell.
   * @param [in] mapping1 Mapping object for the first cell.
   * @return The global quadrature rule on the polygon/polyhedron.
   */
  template <int dim0, int dim1, int spacedim>
  dealii::Quadrature<spacedim>
  compute_quadrature_on_boolean_operation(
    const typename dealii::Triangulation<dim0, spacedim>::cell_iterator &cell0,
    const typename dealii::Triangulation<dim1, spacedim>::cell_iterator &cell1,
    const unsigned int                                                   degree,
    const BooleanOperation        &bool_op,
    const Mapping<dim0, spacedim> &mapping0 =
      (ReferenceCells::get_hypercube<dim0>()
         .template get_default_linear_mapping<dim0, spacedim>()),
    const Mapping<dim1, spacedim> &mapping1 =
      (ReferenceCells::get_hypercube<dim1>()
         .template get_default_linear_mapping<dim1, spacedim>()));

  /**
   * A specialization of the function above when the BooleanOperation is an
   * intersection. The rationale behind this specialization is that deal.II
   * affine cells are convex sets, and as the intersection of convex sets is
   * itself convex, this function internally exploits this to use a cheaper way
   * to mesh the inside.
   *
   *
   * @param [in] cell0 A cell_iterator to the first deal.II cell.
   * @param [in] cell1 A cell_iterator to the second deal.II cell.
   * @param [in] mapping0 Mapping object for the first cell.
   * @param [in] mapping1 Mapping object for the first cell.
   * @param [in] degree The degree of accuracy you wish to get for the global quadrature formula.
   * @return The global quadrature rule on the polygon/polyhedron.
   */
  template <int dim0, int dim1, int spacedim>
  dealii::Quadrature<spacedim>
  compute_quadrature_on_intersection(
    const typename dealii::Triangulation<dim0, spacedim>::cell_iterator &cell0,
    const typename dealii::Triangulation<dim1, spacedim>::cell_iterator &cell1,
    const unsigned int                                                   degree,
    const Mapping<dim0, spacedim> &mapping0,
    const Mapping<dim1, spacedim> &mapping1);

  /**
   * Remesh a CGAL::Surface_mesh.
   *
   * If the domain has 1-dimensional exposed features, the criteria includes a
   * sizing field to guide the sampling of 1-dimensional features with
   * protecting balls centers.
   * - CGAL::parameters::edge_size.
   *
   * This utility should be exploited to improve the quality of a mesh coming
   * from boolean operations. As an example, consider two CGAL::Surface_mesh
   * objects describing two hyper_spheres that intersect non-trivially. After
   * computing their corefinement and union with compute_boolean_operation(),
   * the resulting mesh is the following:
   *
   * @image html boolean_union_hyper_spheres.png
   *
   * Clearly, elements (triangles) like the ones arising along the intersection
   * of the two spheres should be avoided as they're quite degenerate and would
   * decrease the accuracy of the results. A call to the present function will
   * try to remesh the surface, according to the given named parameters, to
   * improve its quality. In the present case, the result is the following:
   *
   * @image html boolean_union_hyper_spheres_remeshed.png
   *
   * @param surface_mesh The input CGAL::Surface_mesh.
   * @param [in] data AdditionalData object to pass to the CGAL::make_mesh_3 function. See the documentation
   * of that struct for a description of those parameters.
   */
  template <typename CGALPointType>
  void
  remesh_surface(CGAL::Surface_mesh<CGALPointType> &surface_mesh,
                 const AdditionalData<3>           &data = AdditionalData<3>{});
} // namespace CGALWrappers

#  ifndef DOXYGEN
// Template implementations
namespace CGALWrappers
{
  template <typename C3t3>
  void
  cgal_surface_mesh_to_cgal_triangulation(
    CGAL::Surface_mesh<typename C3t3::Point::Point> &surface_mesh,
    C3t3                                            &triangulation,
    const AdditionalData<3>                         &data)
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
    Mesh_criteria criteria(CGAL::parameters::facet_size  = data.facet_size,
                           CGAL::parameters::facet_angle = data.facet_angle,
                           CGAL::parameters::facet_distance =
                             data.facet_distance,
                           CGAL::parameters::cell_radius_edge_ratio =
                             data.cell_radius_edge_ratio,
                           CGAL::parameters::cell_size = data.cell_size);
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
    const BooleanOperation                  &boolean_operation,
    CGAL::Surface_mesh<CGALPointType>       &output_surface_mesh)
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
    (void)res;
    Assert(res,
           ExcMessage("The boolean operation was not successfully computed."));
  }



  template <typename CGALTriangulationType>
  dealii::Quadrature<CGALTriangulationType::Point::Ambient_dimension::value>
  compute_quadrature(const CGALTriangulationType &tria,
                     const unsigned int           degree)
  {
    Assert(tria.is_valid(), ExcMessage("The triangulation is not valid."));
    Assert(CGALTriangulationType::Point::Ambient_dimension::value == 3,
           ExcNotImplemented());
    Assert(degree > 0,
           ExcMessage("The degree of the Quadrature formula is not positive."));

    constexpr int spacedim =
      CGALTriangulationType::Point::Ambient_dimension::value;
    std::vector<std::array<dealii::Point<spacedim>, spacedim + 1>>
      vec_of_simplices; // tets

    std::array<dealii::Point<spacedim>, spacedim + 1> simplex;
    for (const auto &f : tria.finite_cell_handles())
      {
        for (unsigned int i = 0; i < (spacedim + 1); ++i)
          {
            simplex[i] =
              cgal_point_to_dealii_point<spacedim>(f->vertex(i)->point());
          }

        vec_of_simplices.push_back(simplex);
      }

    return QGaussSimplex<spacedim>(degree).mapped_quadrature(vec_of_simplices);
  }



  template <int dim0, int dim1, int spacedim>
  dealii::Quadrature<spacedim>
  compute_quadrature_on_boolean_operation(
    const typename dealii::Triangulation<dim0, spacedim>::cell_iterator &cell0,
    const typename dealii::Triangulation<dim1, spacedim>::cell_iterator &cell1,
    const unsigned int                                                   degree,
    const BooleanOperation        &bool_op,
    const Mapping<dim0, spacedim> &mapping0,
    const Mapping<dim1, spacedim> &mapping1)
  {
    Assert(dim0 == 3 && dim1 == 3 && spacedim == 3,
           ExcNotImplemented("2d geometries are not yet supported."));
    if (dim1 > dim0)
      {
        return compute_quadrature_on_boolean_operation(
          cell1,
          cell0,
          degree,
          bool_op,
          mapping1,
          mapping0); // This function works for dim1<=dim0, so swap them if this
                     // is not the case.
      }
    if (bool_op == BooleanOperation::compute_intersection)
      {
        return compute_quadrature_on_intersection(
          cell0, cell1, degree, mapping0, mapping1);
      }
    else
      {
        using K         = CGAL::Exact_predicates_inexact_constructions_kernel;
        using CGALPoint = CGAL::Point_3<K>;
        using CGALTriangulation = CGAL::Triangulation_3<K>;
        CGAL::Surface_mesh<CGALPoint> surface_1, surface_2, out_surface;
        dealii_cell_to_cgal_surface_mesh(cell0, mapping0, surface_1);
        dealii_cell_to_cgal_surface_mesh(cell1, mapping1, surface_2);
        // They have to be triangle meshes
        CGAL::Polygon_mesh_processing::triangulate_faces(surface_1);
        CGAL::Polygon_mesh_processing::triangulate_faces(surface_2);
        Assert(CGAL::is_triangle_mesh(surface_1) &&
                 CGAL::is_triangle_mesh(surface_2),
               ExcMessage("The surface must be a triangle mesh."));
        compute_boolean_operation(surface_1, surface_2, bool_op, out_surface);

        CGALTriangulation tria;
        tria.insert(out_surface.points().begin(), out_surface.points().end());
        return compute_quadrature(tria, degree);
      }
  }



  template <int dim0, int dim1, int spacedim>
  dealii::Quadrature<spacedim>
  compute_quadrature_on_intersection(
    const typename dealii::Triangulation<dim0, spacedim>::cell_iterator &cell0,
    const typename dealii::Triangulation<dim1, spacedim>::cell_iterator &cell1,
    const unsigned int                                                   degree,
    const Mapping<dim0, spacedim> &mapping0,
    const Mapping<dim1, spacedim> &mapping1)
  {
    Assert(dim0 == 3 && dim1 == 3 && spacedim == 3,
           ExcNotImplemented("2d geometries are not yet supported."));
    using K         = CGAL::Exact_predicates_inexact_constructions_kernel;
    using CGALPoint = CGAL::Point_3<K>;
    using CGALTriangulation = CGAL::Triangulation_3<K>;

    CGAL::Surface_mesh<CGALPoint> surface_1, surface_2, out_surface;
    dealii_cell_to_cgal_surface_mesh(cell0, mapping0, surface_1);
    dealii_cell_to_cgal_surface_mesh(cell1, mapping1, surface_2);
    // They have to be triangle meshes
    CGAL::Polygon_mesh_processing::triangulate_faces(surface_1);
    CGAL::Polygon_mesh_processing::triangulate_faces(surface_2);

    compute_boolean_operation(surface_1,
                              surface_2,
                              BooleanOperation::compute_intersection,
                              out_surface);
    CGAL::Surface_mesh<CGALPoint> dummy;
    CGALTriangulation             tr;
    CGAL::convex_hull_3(out_surface.points().begin(),
                        out_surface.points().end(),
                        dummy);
    tr.insert(dummy.points().begin(), dummy.points().end());
    return compute_quadrature(tr, degree);
  }



  template <typename CGALPointType>
  void
  remesh_surface(CGAL::Surface_mesh<CGALPointType> &cgal_mesh,
                 const AdditionalData<3>           &data)
  {
    using K           = CGAL::Exact_predicates_inexact_constructions_kernel;
    using Mesh_domain = CGAL::Polyhedral_mesh_domain_with_features_3<K>;
    // Polyhedron type
    using Polyhedron = CGAL::Mesh_polyhedron_3<K>::type;
    // Triangulation
    using Tr = CGAL::Mesh_triangulation_3<Mesh_domain>::type;
    using C3t3 =
      CGAL::Mesh_complex_3_in_triangulation_3<Tr,
                                              Mesh_domain::Corner_index,
                                              Mesh_domain::Curve_index>;
    // Criteria
    using Mesh_criteria = CGAL::Mesh_criteria_3<Tr>;

    Polyhedron poly;
    CGAL::copy_face_graph(cgal_mesh, poly);
    // Create a vector with only one element: the pointer to the
    // polyhedron.
    std::vector<Polyhedron *> poly_ptrs_vector(1, &poly);
    // Create a polyhedral domain, with only one polyhedron,
    // and no "bounding polyhedron", so the volumetric part of the
    // domain will be empty.
    Mesh_domain domain(poly_ptrs_vector.begin(), poly_ptrs_vector.end());
    // Get sharp features
    domain.detect_features(); // includes detection of borders

    // Mesh criteria
    const double default_value_edge_size = std::numeric_limits<double>::max();
    if (data.edge_size > 0 &&
        std::abs(data.edge_size - default_value_edge_size) > 1e-12)
      {
        Mesh_criteria criteria(CGAL::parameters::edge_size   = data.edge_size,
                               CGAL::parameters::facet_size  = data.facet_size,
                               CGAL::parameters::facet_angle = data.facet_angle,
                               CGAL::parameters::facet_distance =
                                 data.facet_distance,
                               CGAL::parameters::cell_radius_edge_ratio =
                                 data.cell_radius_edge_ratio,
                               CGAL::parameters::cell_size = data.cell_size);
        // Mesh generation
        C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain,
                                            criteria,
                                            CGAL::parameters::no_perturb(),
                                            CGAL::parameters::no_exude());
        cgal_mesh.clear();
        CGAL::facets_in_complex_3_to_triangle_mesh(c3t3, cgal_mesh);
      }
    else if (std::abs(data.edge_size - default_value_edge_size) <= 1e-12)
      {
        Mesh_criteria criteria(CGAL::parameters::facet_size  = data.facet_size,
                               CGAL::parameters::facet_angle = data.facet_angle,
                               CGAL::parameters::facet_distance =
                                 data.facet_distance,
                               CGAL::parameters::cell_radius_edge_ratio =
                                 data.cell_radius_edge_ratio,
                               CGAL::parameters::cell_size = data.cell_size);
        // Mesh generation
        C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain,
                                            criteria,
                                            CGAL::parameters::no_perturb(),
                                            CGAL::parameters::no_exude());
        cgal_mesh.clear();
        CGAL::facets_in_complex_3_to_triangle_mesh(c3t3, cgal_mesh);
      }
    else
      {
        DEAL_II_ASSERT_UNREACHABLE();
      }
  }



  /**
   * Resort vertices in deal.II order to vertices in CGAL order.
   *
   * @param structdim Dimension of the entity which is described by the vertices.
   * @param vertices Vertices in deal.II order which are resorted to CGAL order.
   */
  template <int spacedim>
  void
  resort_dealii_vertices_to_cgal_order(const unsigned int            structdim,
                                       std::vector<Point<spacedim>> &vertices)
  {
    // Mark the two arguments as "used" because some compilers complain about
    // arguments used only within an 'if constexpr' block.
    (void)structdim;
    (void)vertices;

    if constexpr (spacedim == 2)
      if (ReferenceCell::n_vertices_to_type(structdim, vertices.size()) ==
          ReferenceCells::Quadrilateral)
        std::swap(vertices[2], vertices[3]);
  }



  /**
   * Get vertices of cell in CGAL ordering.
   *
   * @param cell A cell_iterator to a deal.II cell.
   * @param mapping Mapping object for the cell.
   * @return  Array of vertices in CGAL order.
   */
  template <int dim, int spacedim>
  std::vector<Point<spacedim>>
  get_vertices_in_cgal_order(
    const typename dealii::Triangulation<dim, spacedim>::cell_iterator &cell,
    const Mapping<dim, spacedim>                                       &mapping)
  {
    // Elements have to be rectangular or simplices
    const unsigned int n_vertices = cell->n_vertices();
    Assert((n_vertices == ReferenceCells::get_hypercube<dim>().n_vertices()) ||
             (n_vertices == ReferenceCells::get_simplex<dim>().n_vertices()),
           ExcNotImplemented());

    std::vector<Point<spacedim>> ordered_vertices(n_vertices);
    std::copy_n(mapping.get_vertices(cell).begin(),
                n_vertices,
                ordered_vertices.begin());

    resort_dealii_vertices_to_cgal_order(dim, ordered_vertices);

    return ordered_vertices;
  }



  /**
   * Get vertices of face in CGAL ordering
   *
   * @param cell A cell_iterator to a deal.II cell.
   * @param face_no The face number within the given cell.
   * @param mapping Mapping object for the cell.
   * @return Array of vertices in CGAL order.
   */
  template <int dim, int spacedim>
  std::vector<Point<spacedim>>
  get_vertices_in_cgal_order(
    const typename dealii::Triangulation<dim, spacedim>::cell_iterator &cell,
    const unsigned int                                                  face_no,
    const Mapping<dim, spacedim>                                       &mapping)
  {
    // Elements have to be rectangular or simplices
    const unsigned int n_vertices = cell->face(face_no)->n_vertices();
    Assert(
      (n_vertices == ReferenceCells::get_hypercube<dim - 1>().n_vertices()) ||
        (n_vertices == ReferenceCells::get_simplex<dim - 1>().n_vertices()),
      ExcNotImplemented());

    std::vector<Point<spacedim>> ordered_vertices(n_vertices);
    std::copy_n(mapping.get_vertices(cell, face_no).begin(),
                n_vertices,
                ordered_vertices.begin());

    resort_dealii_vertices_to_cgal_order(dim - 1, ordered_vertices);

    return ordered_vertices;
  }



} // namespace CGALWrappers
#  endif

DEAL_II_NAMESPACE_CLOSE

#else

// Make sure the scripts that create the C++20 module input files have
// something to latch on if the preprocessor #ifdef above would
// otherwise lead to an empty content of the file.
DEAL_II_NAMESPACE_OPEN
DEAL_II_NAMESPACE_CLOSE

#endif
#endif
