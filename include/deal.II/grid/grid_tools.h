// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2001 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_grid_tools_h
#define dealii_grid_tools_h


#include <deal.II/base/config.h>

#include <deal.II/base/bounding_box.h>
#include <deal.II/base/geometry_info.h>
#include <deal.II/base/point.h>

#include <deal.II/boost_adaptors/bounding_box.h>

#include <deal.II/distributed/shared_tria.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping.h>

#include <deal.II/grid/grid_tools_geometry.h>
#include <deal.II/grid/grid_tools_topology.h>
#include <deal.II/grid/manifold.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/lac/petsc_vector.h>
#include <deal.II/lac/sparsity_tools.h>
#include <deal.II/lac/trilinos_vector.h>

#include <deal.II/numerics/rtree.h>

#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/serialization/array.hpp>
#include <boost/serialization/vector.hpp>

#ifdef DEAL_II_WITH_ZLIB
#  include <boost/iostreams/device/back_inserter.hpp>
#  include <boost/iostreams/filter/gzip.hpp>
#  include <boost/iostreams/filtering_stream.hpp>
#  include <boost/iostreams/stream.hpp>
#endif

#include <optional>
#include <set>

#ifdef DEAL_II_HAVE_CXX20
#  include <concepts>
#endif


DEAL_II_NAMESPACE_OPEN

// Forward declarations
#ifndef DOXYGEN
namespace parallel
{
  namespace distributed
  {
    template <int dim, int spacedim>
    DEAL_II_CXX20_REQUIRES((concepts::is_valid_dim_spacedim<dim, spacedim>))
    class Triangulation;
  }
} // namespace parallel

namespace hp
{
  template <int, int>
  class MappingCollection;
}

class SparsityPattern;

namespace GridTools
{
  template <int dim, int spacedim>
  class Cache;
}
#endif

namespace internal
{
  template <int dim, int spacedim, typename MeshType>
  DEAL_II_CXX20_REQUIRES(concepts::is_triangulation_or_dof_handler<MeshType>)
  class ActiveCellIterator
  {
  public:
#ifndef _MSC_VER
    using type = typename MeshType::active_cell_iterator;
#else
    using type = TriaActiveIterator<dealii::CellAccessor<dim, spacedim>>;
#endif
  };

#ifdef _MSC_VER
  template <int dim, int spacedim>
  class ActiveCellIterator<dim, spacedim, DoFHandler<dim, spacedim>>
  {
  public:
    using type =
      TriaActiveIterator<dealii::DoFCellAccessor<dim, spacedim, false>>;
  };
#endif
} // namespace internal

/**
 * This namespace is a collection of algorithms working on triangulations,
 * such as shifting or rotating triangulations, but also finding a cell that
 * contains a given point. See the descriptions of the individual functions
 * for more information.
 *
 * @ingroup grid
 */
namespace GridTools
{
  /**
   * @name Rotating, stretching and otherwise transforming meshes
   */
  /** @{ */

  /**
   * Transform the vertices of the given triangulation by applying the
   * function object provided as first argument to all its vertices.
   *
   * The transformation given as argument is used to transform each vertex.
   * Its respective type has to offer a function-like syntax, i.e. the
   * predicate is either an object of a type that has an <tt>operator()</tt>,
   * or it is a pointer to a non-member function, or it is a lambda function
   * object. In either case, argument and return
   * value have to be of type Point<spacedim>. An example -- a simple
   * transformation that moves the object two units to the
   * right in the $x_1$ direction -- could look like as follows:
   * @code
   *   Triangulation<dim> triangulation;
   *   ... // fill triangulation with something
   *   GridTools::transform ([](const Point<dim> &p) -> Point<dim>
   *                         {
   *                           Point<dim> q = p;
   *                           q[0] += 2;
   *                           return q;
   *                         },
   *                         triangulation);
   * @endcode
   * Here, the transformation is provided by a lambda function that
   * takes a `Point<dim>` as input and returns a `Point<dim>` as output.
   *
   * @note The transformations that make sense to use with this function
   *   should have a Jacobian with a positive determinant. For example,
   *   rotation, shearing, stretching, or scaling all satisfy this (though
   *   there is no requirement that the transformation used actually is
   *   linear, as all of these examples are). On the other hand, reflections
   *   or inversions have a negative determinant of the Jacobian. The
   *   current function has no way of asserting a positive determinant
   *   of the Jacobian, but if you happen to use such a transformation,
   *   the result will be a triangulation in which cells have a negative
   *   volume.
   *
   * @note If you are using a parallel::distributed::Triangulation you will
   * have hanging nodes in your local Triangulation even if your "global" mesh
   * has no hanging nodes. This will cause issues with wrong positioning of
   * hanging nodes in ghost cells if you call the current functions: The
   * vertices of all locally owned cells will be correct, but the vertices of
   * some ghost cells may not. This means that computations like
   * KellyErrorEstimator may give wrong answers.
   *
   * @note This function is in general not compatible with manifolds attached
   * to the triangulation. For example, in order to refine the grid (using
   * manifolds) after the grid transformation, you have to make sure that
   * the original manifold is still valid for the transformed geometry. This
   * does not hold in general, and it is necessary to clear manifolds from
   * the triangulation (for example, using Triangulation::clear_all_manifolds())
   * before the transformation, and then attach new ones after the
   * transformation that are valid for the transformed geometry. There are cases
   * where this is awkward, most notably if you are using a mesh generated by
   * the functions in GridGenerator which generally attach suitable manifolds
   * upon mesh generation; in those cases, you will have to think about how
   * these manifolds were constructed, and create a manifold that is constructed
   * in a similar way but applies to the transformed geometry. As a consequence,
   * if you only care about manifolds for mesh refinement, it is often simpler
   * to just refine the original mesh before transformation as needed, and then
   * simply forget about the manifolds. Of course, manifolds are also used for
   * other cases (e.g., for normal vectors, curved edges and faces, and higher
   * order mappings), and if these are relevant to what you are doing, then
   * there is no alternative to building appropriate manifolds for the
   * transformed geometry and attaching these to the transformed geometry. In
   * general, detaching manifolds from a triangulation and then doing the
   * transformation would look as follows:
   * @code
   * ...
   * triangulation.reset_all_manifolds();
   * GridTools::transform(MyTransformation<dim>(), triangulation);
   * ...
   * @endcode
   *
   * This function is used in the "Possibilities for extensions" section of
   * step-38. It is also used in step-49 and step-53.
   *
   * @dealiiConceptRequires{(std::invocable<Transformation, Point<spacedim>> &&
   *    std::assignable_from<Point<spacedim> &,
   *    std::invoke_result_t<Transformation, Point<spacedim>>>)}
   */
  template <int dim, typename Transformation, int spacedim>
  DEAL_II_CXX20_REQUIRES(
    (std::invocable<Transformation, Point<spacedim>> &&
     std::assignable_from<
       Point<spacedim> &,
       std::invoke_result_t<Transformation, Point<spacedim>>>))
  void transform(const Transformation         &transformation,
                 Triangulation<dim, spacedim> &triangulation);

  /**
   * Shift each vertex of the triangulation by the given shift vector. This
   * function uses the transform() function above, so the requirements on the
   * triangulation stated there hold for this function as well; in particular,
   * this is true about the discussion about manifolds.
   */
  template <int dim, int spacedim>
  void
  shift(const Tensor<1, spacedim>    &shift_vector,
        Triangulation<dim, spacedim> &triangulation);


  /**
   * Rotate all vertices of the given two-dimensional triangulation in
   * counter-clockwise sense around the origin of the coordinate system by the
   * given angle (given in radians, rather than degrees). This function uses
   * the transform() function above, so the requirements on the triangulation
   * stated there hold for this function as well; in particular,
   * this is true about the discussion about manifolds.
   *
   * @note This function is only supported for spacedim=2.
   */
  template <int dim, int spacedim>
  void
  rotate(const double angle, Triangulation<dim, spacedim> &triangulation);

  /**
   * Rotate all vertices of the given @p triangulation in counter-clockwise
   * direction around the @p axis described by a unit vector. Otherwise like the
   * function above; in particular, this function calls the transform() function
   * and so the discussion about manifolds there also applies here.
   *
   * @param[in] angle Angle in radians to rotate the Triangulation by.
   * @param[in] axis A unit vector that defines the axis of rotation.
   * @param[in,out] triangulation The Triangulation object to rotate.
   *
   * @note Implemented for spacedim=3 and dim=1, 2, and 3.
   */
  template <int dim>
  void
  rotate(const Tensor<1, 3, double> &axis,
         const double                angle,
         Triangulation<dim, 3>      &triangulation);

  /**
   * Transform the given triangulation smoothly to a different domain where,
   * typically, each of the vertices at the boundary of the triangulation is
   * mapped to the corresponding points in the @p new_points map.
   *
   * The unknown displacement field $u_d(\mathbf x)$ in direction $d$ is
   * obtained from the minimization problem \f[ \min\, \int \frac{1}{2}
   *   c(\mathbf x)
   *   \mathbf \nabla u_d(\mathbf x) \cdot
   *   \mathbf \nabla u_d(\mathbf x)
   *   \,\rm d x
   * \f]
   * subject to prescribed constraints. The minimizer is obtained by solving the
   * Laplace equation of the dim components of a displacement field that maps
   * the current
   * domain into one described by @p new_points . Linear finite elements with
   * four Gaussian quadrature points in each direction are used. The difference
   * between the vertex positions specified in @p new_points and their current
   * value in @p tria therefore represents the prescribed values of this
   * displacement field at the boundary of the domain, or more precisely at all
   * of those locations for which @p new_points provides values (which may be
   * at part of the boundary, or even in the interior of the domain). The
   * function then evaluates this displacement field at each unconstrained
   * vertex and uses it to place the mapped vertex where the displacement
   * field locates it. Because the solution of the Laplace equation is smooth,
   * this guarantees a smooth mapping from the old domain to the new one.
   *
   * @param[in] new_points The locations where a subset of the existing
   * vertices are to be placed. Typically, this would be a map from the vertex
   * indices of all nodes on the boundary to their new locations, thus
   * completely specifying the geometry of the mapped domain. However, it may
   * also include interior points if necessary and it does not need to include
   * all boundary vertices (although you then lose control over the exact
   * shape of the mapped domain).
   *
   * @param[in,out] tria The Triangulation object. This object is changed
   * in-place, i.e., the previous locations of vertices are overwritten.
   *
   * @param[in] coefficient An optional coefficient for the Laplace problem.
   * Larger values make cells less prone to deformation (effectively
   * increasing their stiffness). The coefficient is evaluated in the
   * coordinate system of the old, undeformed configuration of the
   * triangulation as input, i.e., before the transformation is applied.
   * Should this function be provided, sensible results can only be expected
   * if all coefficients are positive.
   *
   * @param[in] solve_for_absolute_positions If set to <code>true</code>, the
   * minimization problem is formulated with respect to the final vertex
   * positions as opposed to their displacement. The two formulations are
   * equivalent for
   * the homogeneous problem (default value of @p coefficient), but they
   * result in very different mesh motion otherwise. Since in most cases one
   * will be using a non-constant coefficient in displacement formulation, the
   * default value of this parameter is <code>false</code>.
   *
   * @note This function is not currently implemented for the 1d case.
   */
  template <int dim>
  void
  laplace_transform(const std::map<unsigned int, Point<dim>> &new_points,
                    Triangulation<dim>                       &tria,
                    const Function<dim, double> *coefficient = nullptr,
                    const bool solve_for_absolute_positions  = false);

  /**
   * Scale the entire triangulation by the given factor. To preserve the
   * orientation of the triangulation, the factor must be positive.
   *
   * This function uses the transform() function above, so the requirements on
   * the triangulation stated there hold for this function as well.
   */
  template <int dim, int spacedim>
  void
  scale(const double                  scaling_factor,
        Triangulation<dim, spacedim> &triangulation);

  /**
   * Distort the given triangulation by randomly moving around all the
   * vertices of the grid.  The direction of movement of each vertex is
   * random, while the length of the shift vector has a value of @p factor
   * times the minimal length of the active edges adjacent to this vertex.
   * Note that @p factor should obviously be well below <tt>0.5</tt> in order
   * to avoid getting cells that are @ref GlossDistorted "distorted".
   *
   * The function will make sure that vertices on restricted faces
   * (i.e., faces with hanging nodes) will end up in the correct
   * place, i.e. in the middle of the two other vertices of the parent
   * edge, and the analogue in higher space dimensions (vertices on
   * the boundary are not corrected, so don't distort boundary
   * vertices in more than two space dimensions, i.e. in dimensions
   * where boundary vertices can be hanging nodes).
   *
   * If @p keep_boundary is set to @p true (which is the default), then
   * boundary vertices are not moved.
   *
   * @p seed is used for the initialization of the random engine. Its
   * default value initializes the engine with the same state as in
   * previous versions of deal.II.
   *
   * @note If the Triangulation is of distributed kind (derived from
   * parallel::DistributedTriangulationBase) and computations are done in
   * parallel, the new vertex locations will be consistently updated
   * on all ranks.
   *
   */
  template <int dim, int spacedim>
  void
  distort_random(
    const double                  factor,
    Triangulation<dim, spacedim> &triangulation,
    const bool                    keep_boundary = true,
    const unsigned int            seed = boost::random::mt19937::default_seed);

  /**
   * Analyze the boundary cells of a mesh, and if one cell is found at
   * a corner position (with dim adjacent faces on the boundary), and its
   * dim-dimensional angle fraction exceeds @p limit_angle_fraction,
   * refine globally once, and replace the children of such cell
   * with children where the corner is no longer offending the given angle
   * fraction.
   *
   * If no boundary cells exist with two adjacent faces on the boundary, then
   * the triangulation is left untouched. If instead we do have cells with dim
   * adjacent faces on the boundary, then the fraction between the
   * dim-dimensional
   * solid angle and dim*pi/2 is checked against the parameter @p limit_angle_fraction.
   * If it is higher, the grid is refined once, and the children of the
   * offending cell are replaced with some cells that instead respect the limit.
   * After this process the triangulation is flattened, and all Manifold objects
   * are restored as they were in the original triangulation.
   *
   * An example is given by the following mesh, obtained by attaching a
   * SphericalManifold to a mesh generated using GridGenerator::hyper_cube:
   *
   * @code
   * const SphericalManifold<dim> m0;
   * Triangulation<dim> tria;
   * GridGenerator::hyper_cube(tria,-1,1);
   * tria.set_all_manifold_ids_on_boundary(0);
   * tria.set_manifold(0, m0);
   * tria.refine_global(4);
   * @endcode
   *
   * <p ALIGN="center">
   * @image html regularize_mesh_01.png
   * </p>
   *
   * The four cells that were originally the corners of a square will give you
   * some troubles during computations, as the jacobian of the transformation
   * from the reference cell to those cells will go to zero, affecting the error
   * constants of the finite element estimates.
   *
   * Those cells have a corner with an angle that is very close to 180 degrees,
   * i.e., an angle fraction very close to one.
   *
   * The same code, adding a call to regularize_corner_cells:
   * @code
   * const SphericalManifold<dim> m0;
   * Triangulation<dim> tria;
   * GridGenerator::hyper_cube(tria,-1,1);
   * tria.set_all_manifold_ids_on_boundary(0);
   * tria.set_manifold(0, m0);
   * GridTools::regularize_corner_cells(tria);
   * tria.refine_global(2);
   * @endcode
   * generates a mesh that has a much better behavior w.r.t. the jacobian of
   * the Mapping:
   *
   * <p ALIGN="center">
   * @image html regularize_mesh_02.png
   * </p>
   *
   * This mesh is very similar to the one obtained by GridGenerator::hyper_ball.
   * However, using GridTools::regularize_corner_cells one has the freedom to
   * choose when to apply the regularization, i.e., one could in principle first
   * refine a few times, and then call the regularize_corner_cells function:
   *
   * @code
   * const SphericalManifold<dim> m0;
   * Triangulation<dim> tria;
   * GridGenerator::hyper_cube(tria,-1,1);
   * tria.set_all_manifold_ids_on_boundary(0);
   * tria.set_manifold(0, m0);
   * tria.refine_global(2);
   * GridTools::regularize_corner_cells(tria);
   * tria.refine_global(1);
   * @endcode
   *
   * This generates the following mesh:
   *
   * <p ALIGN="center">
   * @image html regularize_mesh_03.png
   * </p>
   *
   * The function is currently implemented only for dim = 2 and
   * will throw an exception if called with dim = 3.
   *
   * @param[in,out] tria Triangulation to regularize.
   *
   * @param[in] limit_angle_fraction Maximum ratio of angle or solid
   * angle that is allowed for a corner element in the mesh.
   */
  template <int dim, int spacedim>
  void
  regularize_corner_cells(Triangulation<dim, spacedim> &tria,
                          const double limit_angle_fraction = .75);

  /** @} */
  /**
   * @name Finding cells and vertices of a triangulation
   */
  /** @{ */

  /**
   * Given a Triangulation's @p cache and a list of @p points, call
   * find_active_cell_around_point() on each element of @p points, and return
   * @p cells, reference positions @p qpoints, and a mapping @p maps from local
   * to global indices into @p points .
   *
   * @param[in] cache The triangulation's GridTools::Cache .
   * @param[in] points A vector of points.
   * @param[in] cell_hint (optional) A cell iterator for a cell which likely
   * contains the first point of @p points.
   *
   * @return A tuple containing the following information:
   *  - @p cells : A vector of all the cells containing at least one of
   *   the @p points .
   *  - @p qpoints : A vector of vectors of points. @p qpoints[i] contains
   *   the reference positions of all points that fall within the cell @p cells[i] .
   *  - @p indices : A vector of vectors of integers, containing the mapping between
   *   local numbering the @p cells array (the first component of the returned
   *   tuple), and global index in the input array @p points . In other words,
   *   the indices stored in the array `indices[c]` correspond to those points
   *   of the input argument `points` that are located on `cells[c]`.
   *
   * If @p points[a] and @p points[b] are the only two points that fall in @p cells[c],
   * then @p qpoints[c][0] and @p qpoints[c][1] are the reference positions of
   * @p points[a] and @p points[b] in @p cells[c], and @p indices[c][0] = a,
   * @p indices[c][1] = b. The function
   * Mapping::transform_unit_to_real(qpoints[c][0])
   * returns @p points[a].
   *
   * The algorithm builds an rtree of @p points to sort them spatially, before
   * attempting to call find_active_cell_around_point().
   *
   * @note This function is not implemented for the codimension one case (<tt>spacedim != dim</tt>).
   *
   * @note If a point is not found inside the mesh, or is lying inside an
   * artificial cell of a parallel::TriangulationBase, the point is silently
   * ignored. If you want to infer for which points the search failed, use the
   * function compute_point_locations_try_all() that also returns a vector of
   * indices indicating the points for which the search failed.
   *
   * @note The actual return type of this function, i.e., the type referenced
   * above as @p return_type, is
   * @code
   * std::tuple<
   *   std::vector<
   *     typename Triangulation<dim, spacedim>::active_cell_iterator>,
   *   std::vector<std::vector<Point<dim>>>,
   *   std::vector<std::vector<unsigned int>>>
   * @endcode
   * The type is abbreviated in the online documentation to improve readability
   * of this page.
   *
   * @note This function optimizes the search by making use of
   * GridTools::Cache::get_cell_bounding_boxes_rtree(), which either returns
   * a cached rtree or builds and stores one. Building an rtree might hinder
   * the performance if the function is called only once on few points.
   */
  template <int dim, int spacedim>
#ifndef DOXYGEN
  std::tuple<
    std::vector<typename Triangulation<dim, spacedim>::active_cell_iterator>,
    std::vector<std::vector<Point<dim>>>,
    std::vector<std::vector<unsigned int>>>
#else
  return_type
#endif
  compute_point_locations(
    const Cache<dim, spacedim>         &cache,
    const std::vector<Point<spacedim>> &points,
    const typename Triangulation<dim, spacedim>::active_cell_iterator
      &cell_hint =
        typename Triangulation<dim, spacedim>::active_cell_iterator());

  /**
   * This function is similar to GridTools::compute_point_locations(),
   * but while compute_point_locations() silently ignores all points for which
   * find_active_cell_around_point() fails, this function also returns a
   * vector containing the indices of the points for which
   * find_active_cell_around_point() failed.
   *
   * @return A tuple containing four elements; the first three
   * are documented in GridTools::compute_point_locations().
   * The last element of the @p return_type contains the
   * indices of points which are neither found inside the mesh
   * nor lie in artificial cells. The @p return_type equals the
   * following tuple type:
   * @code
   *   std::tuple<
   *     std::vector<
   *        typename Triangulation<dim,spacedim>::active_cell_iterator>,
   *     std::vector<std::vector<Point<dim>>>,
   *     std::vector<std::vector<unsigned int>>,
   *     std::vector<unsigned int>
   *   >
   * @endcode
   *
   * @note This function is not implemented for the codimension one case (<tt>spacedim != dim</tt>).
   *
   * @note This function optimizes the search by making use of
   * GridTools::Cache::get_cell_bounding_boxes_rtree(), which either returns
   * a cached rtree or builds and stores one. Building an rtree might hinder
   * the performance if the function is called only once on few points.
   *
   * For a more detailed documentation see
   * GridTools::compute_point_locations().
   */
  template <int dim, int spacedim>
#ifndef DOXYGEN
  std::tuple<
    std::vector<typename Triangulation<dim, spacedim>::active_cell_iterator>,
    std::vector<std::vector<Point<dim>>>,
    std::vector<std::vector<unsigned int>>,
    std::vector<unsigned int>>
#else
  return_type
#endif
  compute_point_locations_try_all(
    const Cache<dim, spacedim>         &cache,
    const std::vector<Point<spacedim>> &points,
    const typename Triangulation<dim, spacedim>::active_cell_iterator
      &cell_hint =
        typename Triangulation<dim, spacedim>::active_cell_iterator());

  /**
   * Given a @p cache and a list of
   * @p local_points for each process, find the points lying on the locally
   * owned part of the mesh and compute the quadrature rules for them.
   * Distributed compute point locations is a function similar to
   * GridTools::compute_point_locations but working for
   * parallel::TriangulationBase objects and, unlike its serial version, also
   * for a distributed triangulation (see parallel::distributed::Triangulation).
   *
   * @param[in] cache a GridTools::Cache object
   * @param[in] local_points the array of points owned by the current process.
   * Every process can have a different array of points which can be empty and
   * not contained within the locally owned part of the triangulation
   * @param[in] global_bboxes a vector of vectors of bounding boxes; it
   * describes the locally owned part of the mesh for each process. The bounding
   * boxes describing which part of the mesh is locally owned by process with
   * rank rk are contained in global_bboxes[rk]. The local description can be
   * obtained from GridTools::compute_mesh_predicate_bounding_box; then the
   * global one can be obtained using either
   * GridTools::exchange_local_bounding_boxes or Utilities::MPI::all_gather
   * @param[in] tolerance Tolerance in terms of unit cell coordinates. Depending
   *   on the problem, it might be necessary to adjust the tolerance in order
   *   to be able to identify a cell. Floating
   *   point arithmetic implies that a point will, in general, not lie exactly
   *   on a vertex, edge, or face. In either case, it is not predictable which
   *   of the cells adjacent to a vertex or an edge/face this function returns.
   *   Consequently, algorithms that call this function need to take into
   *   account that the returned cell will only contain the point approximately.
   * @param[in] enforce_unique_mapping Enforce a one to one mapping between
   points
   *   in real and reference space.
   * @param[in] marked_vertices An array of bools indicating which
   * vertices of @p mesh will be considered within the search
   * as the potentially closest vertex. On receiving a non-empty
   * @p marked_vertices, the function will
   * only search among @p marked_vertices for the closest vertex,
   * otherwise on all vertices in the mesh.

   * @return A tuple containing the quadrature information
   *
   * The elements of the output tuple are:
   * - cells : a vector of all cells containing at least one point.
   * - qpoints : a vector of vector of points; containing in @p qpoints[i]
   *   the reference positions of all points that fall within the cell @p cells[i] .
   * - maps : a vector of vector of integers, containing the mapping between
   *  the numbering in qpoints (previous element of the tuple), and the vector
   *  of local points of the process owning the points.
   * - points : a vector of vector of points. @p points[i][j] is the point in the
   *  real space corresponding.
   *  to @p qpoints[i][j] . Notice @p points are the points lying on the locally
   *  owned part of the mesh; thus these can be either copies of @p local_points
   *  or points received from other processes i.e. local_points for other
   * processes
   * - owners : a vector of vectors; @p owners[i][j] contains the rank of
   *  the process owning the point[i][j] (previous element of the tuple).
   *
   * The function uses the triangulation's mpi communicator: for this reason it
   * throws an assert error if the Triangulation is not derived from
   * parallel::TriangulationBase .
   *
   * In a serial execution the first three elements of the tuple are the same
   * as in GridTools::compute_point_locations .
   *
   * Note: this function is a @ref GlossCollectiveOperation "collective operation".
   *
   * @note The actual return type of this function, i.e., the type referenced
   * above as @p return_type, is
   * @code
   * std::tuple<
   *   std::vector<
   *     typename Triangulation<dim, spacedim>::active_cell_iterator>,
   *   std::vector<std::vector<Point<dim>>>,
   *   std::vector<std::vector<unsigned int>>,
   *   std::vector<std::vector<Point<spacedim>>>,
   *   std::vector<std::vector<unsigned int>>>
   * @endcode
   * The type is abbreviated in the online documentation to improve readability
   * of this page.
   */
  template <int dim, int spacedim>
#ifndef DOXYGEN
  std::tuple<
    std::vector<typename Triangulation<dim, spacedim>::active_cell_iterator>,
    std::vector<std::vector<Point<dim>>>,
    std::vector<std::vector<unsigned int>>,
    std::vector<std::vector<Point<spacedim>>>,
    std::vector<std::vector<unsigned int>>>
#else
  return_type
#endif
  distributed_compute_point_locations(
    const GridTools::Cache<dim, spacedim>                 &cache,
    const std::vector<Point<spacedim>>                    &local_points,
    const std::vector<std::vector<BoundingBox<spacedim>>> &global_bboxes,
    const double                                           tolerance = 1e-10,
    const std::vector<bool>                               &marked_vertices = {},
    const bool enforce_unique_mapping = true);

  namespace internal
  {
    /**
     * Data structure returned by
     * GridTools::internal::distributed_compute_point_locations(). It provides
     * information to perform GridTools::distributed_compute_point_locations()
     * and to set up the communication pattern within
     * Utilities::MPI::RemotePointEvaluation::reinit().
     *
     * @note The name of the fields are chosen with
     *   Utilities::MPI::RemotePointEvaluation in mind. Here, quantities are
     *   computed at specified arbitrary positioned points (and even on remote
     *   processes in the MPI universe) cell by cell and these values are sent
     *   to requesting processes, which receive the result and resort the
     *   result according to the points.
     */
    template <int dim, int spacedim>
    struct DistributedComputePointLocationsInternal
    {
      DistributedComputePointLocationsInternal();

      /**
       * Function which sets up @p send_ranks, @p send_ptrs, @p recv_ranks,
       * and @p recv_ptrs from @p send_components, @p recv_components,
       * and @p n_searched_points. Internally @p send_components and @p recv_components
       * are sorted and enumerated.
       */
      void
      finalize_setup();

      /**
       * Number of searched point locations.
       */
      unsigned int n_searched_points;

      /**
       * Information of each point on sending/evaluation side. The elements of
       * the tuple are as follows: 0) cell level and index, 1) rank of the
       * owning process, 2) local index of the owning process, 3) reference
       * position, 4) real position, 5) permutation index within a send buffer.
       */
      std::vector<std::tuple<std::pair<int, int>,
                             unsigned int,
                             unsigned int,
                             Point<dim>,
                             Point<spacedim>,
                             unsigned int>>
        send_components;

      /**
       * Ranks to send to.
       */
      std::vector<unsigned int> send_ranks;

      /**
       * Pointers of ranges within a send buffer to be sent to the ranks
       * specified by send_ranks. The size of the send buffer is given
       * by send_ptrs.back().
       */
      std::vector<unsigned int> send_ptrs;

      /**
       * Information of each received data value. The elements of the tuple are
       * as follows: 0) rank of sender, 1) local index, 2) enumeration index.
       *
       * @note The vector is sorted according to 1), 0), 2).
       *
       * @note To each point multiple data values might be associated to. This
       *   might be the case if a point coincides with a geometric entity (e.g.,
       *   vertex) that is shared by multiple cells.
       */
      std::vector<std::tuple<unsigned int, unsigned int, unsigned int>>
        recv_components;

      /**
       * Ranks from where data is received.
       */
      std::vector<unsigned int> recv_ranks;

      /**
       * Pointers of ranges within a receive buffer that are filled by ranks
       * specified by recv_ranks. The size of the receive buffer is given by
       * recv_ptrs.back().
       */
      std::vector<unsigned int> recv_ptrs;
    };

    /**
     * A function that fills DistributedComputePointLocationsInternal.
     * If the input argument @p perform_handshake is set to false only
     * the fields needed by
     * GridTools::internal::distributed_compute_point_locations() are filled.
     * If the input argument is set to true additional data structures are
     * set up to be able to set up the communication pattern within
     * Utilities::MPI::RemotePointEvaluation::reinit().
     */
    template <int dim, int spacedim>
    DistributedComputePointLocationsInternal<dim, spacedim>
    distributed_compute_point_locations(
      const GridTools::Cache<dim, spacedim>                 &cache,
      const std::vector<Point<spacedim>>                    &points,
      const std::vector<std::vector<BoundingBox<spacedim>>> &global_bboxes,
      const std::vector<bool>                               &marked_vertices,
      const double                                           tolerance,
      const bool                                             perform_handshake,
      const bool enforce_unique_mapping = false);


    /**
     * Data structure returned by
     * GridTools::internal::distributed_compute_intersection_locations(). It can
     * be converted to
     * GridTools::internal::DistributedComputePointLocationsInternal, which can
     * be used to reinit Utilities::MPI::RemotePointEvaluation.
     */
    template <int structdim, int spacedim>
    struct DistributedComputeIntersectionLocationsInternal
    {
      /**
       * Intersections are assumed to be simplices (as, e.g., provided by CGAL)
       */
      using IntersectionType =
        std::array<dealii::Point<spacedim>, structdim + 1>;

      /**
       * Information of each intersection on sending/evaluation side. The
       * elements of the tuple are as follows: 0) cell level and index, 1) rank
       * of the owning process, 2) local index of the owning process, 3)
       * found intersection.
       *
       * @note The vector is sorted according to 1), 2).
       */
      std::vector<std::tuple<std::pair<int, int>,
                             unsigned int,
                             unsigned int,
                             IntersectionType>>
        send_components;

      /**
       * Information of each received data value. The elements of the tuple are
       * as follows: 0) rank of sender, 1) local index, 2) found intersections.
       *
       * @note The vector is sorted according to 1), 0), 2).
       *
       * @note Multiple intersections between cells can be found
       */
      std::vector<std::tuple<unsigned int, unsigned int, IntersectionType>>
        recv_components;

      /**
       * Pointers of ranges to found intersections for requested intersection.
       */
      std::vector<unsigned int> recv_ptrs;

      /**
       * Distribute quadrature points according to
       * QGaussSimplex<structdim>(n_points_1D) on found intersections and
       * construct GridTools::internal::DistributedComputePointLocationsInternal
       * from class members. This can be done without searching for points again
       * since all information is locally known.
       *
       * @p mapped_quadratures_recv_comp is a pointer to an empty vector of
       * mapped quadratures. By default it is a `nullptr` and the parameter is
       * ignored. Otherwise, the vector is filled with the mapped quadrature
       * rules (in real coordinates) corresponding to recv_components.
       *
       * The parameter @p consistent_numbering_of_sender_and_receiver can be used to ensure
       * points on sender and receiver side are numbered consistently.
       * This parameter is optional if DistributedComputePointLocationsInternal
       * is used to set up RemotePointEvaluation, but might be helpful for
       * debugging or other usage of DistributedComputePointLocationsInternal.
       * Note that setting this parameter true requires an additional
       * communication step during the setup phase.
       */
      template <int dim>
      GridTools::internal::DistributedComputePointLocationsInternal<dim,
                                                                    spacedim>
      convert_to_distributed_compute_point_locations_internal(
        const unsigned int                  n_points_1D,
        const Triangulation<dim, spacedim> &tria,
        const Mapping<dim, spacedim>       &mapping,
        std::vector<Quadrature<spacedim>>  *mapped_quadratures_recv_comp =
          nullptr,
        const bool consistent_numbering_of_sender_and_receiver = false) const;

    private:
      /**
       * Helper function for
       * convert_to_distributed_compute_point_locations_internal(). It sends the
       * indices associated to quadrature points at the receiver side to the
       * sender side, where the information is needed to build
       * GridTools::internal::DistributedComputePointLocationsInternal::send_components
       */
      std::map<unsigned int, std::vector<unsigned int>>
      communicate_indices(
        const std::vector<std::tuple<unsigned int, unsigned int, unsigned int>>
                      &point_recv_components,
        const MPI_Comm comm) const;
    };

    /**
     * A function that fills DistributedComputeIntersectionLocationsInternal.
     * @p intersection_requests are vertices of cells which describe the
     * entities we want to search intersections for.
     * The template parameter @p structdim provides the dimension of the
     * resulting intersection.
     */
    template <int structdim, int dim, int spacedim>
    DistributedComputeIntersectionLocationsInternal<structdim, spacedim>
    distributed_compute_intersection_locations(
      const Cache<dim, spacedim>                      &cache,
      const std::vector<std::vector<Point<spacedim>>> &intersection_requests,
      const std::vector<std::vector<BoundingBox<spacedim>>> &global_bboxes,
      const std::vector<bool>                               &marked_vertices,
      const double                                           tolerance);

  } // namespace internal

  /**
   * Find and return the index of the closest vertex to a given point in the
   * map of vertices passed as the first argument.
   *
   * @param vertices A map of index->vertex, as returned by
   *        GridTools::extract_used_vertices().
   * @param p The target point.
   * @return The index of the vertex that is closest to the target point `p`.
   */
  template <int spacedim>
  unsigned int
  find_closest_vertex(const std::map<unsigned int, Point<spacedim>> &vertices,
                      const Point<spacedim>                         &p);

  /**
   * Find and return the index of the used vertex (or marked vertex) in a
   * given mesh that is located closest to a given point.
   *
   * This function uses the locations of vertices as stored in the
   * triangulation. This is usually sufficient, unless you are using a Mapping
   * that moves the vertices around (for example, MappingQEulerian). In this
   * case, you should call the function with the same name and with an
   * additional Mapping argument.
   *
   * @param mesh A variable of a type that satisfies the requirements of the
   * @ref ConceptMeshType "MeshType concept".
   * @param p The point for which we want to find the closest vertex.
   * @param marked_vertices An array of bools indicating which
   * vertices of @p mesh will be considered within the search
   * as the potentially closest vertex. On receiving a non-empty
   * @p marked_vertices, the function will
   * only search among @p marked_vertices for the closest vertex.
   * The size of this array should be equal to the value returned by
   * Triangulation::n_vertices() for the triangulation underlying the given mesh
   * (as opposed to the value returned by Triangulation::n_used_vertices()).
   * @return The index of the closest vertex found.
   *
   * @dealiiConceptRequires{
   *   concepts::is_triangulation_or_dof_handler<MeshType<dim, spacedim>>}
   */
  template <int dim, template <int, int> class MeshType, int spacedim>
  DEAL_II_CXX20_REQUIRES(
    (concepts::is_triangulation_or_dof_handler<MeshType<dim, spacedim>>))
  unsigned int find_closest_vertex(
    const MeshType<dim, spacedim> &mesh,
    const Point<spacedim>         &p,
    const std::vector<bool>       &marked_vertices = {});

  /**
   * Find and return the index of the used vertex (or marked vertex) in a
   * given mesh that is located closest to a given point. Use the given
   * mapping to compute the actual location of the vertices.
   *
   * If the Mapping does not modify the position of the mesh vertices (like,
   * for example, MappingQEulerian does), then this function is equivalent to
   * the one with the same name, and without the `mapping` argument.
   *
   * @param mapping A mapping used to compute the vertex locations
   * @param mesh A variable of a type that satisfies the requirements of the
   * @ref ConceptMeshType "MeshType concept".
   * @param p The point for which we want to find the closest vertex.
   * @param marked_vertices An array of bools indicating which
   * vertices of @p mesh will be considered within the search
   * as the potentially closest vertex. On receiving a non-empty
   * @p marked_vertices, the function will
   * only search among @p marked_vertices for the closest vertex.
   * The size of this array should be equal to the value returned by
   * Triangulation::n_vertices() for the triangulation underlying the given mesh
   * (as opposed to the value returned by Triangulation::n_used_vertices()).
   * @return The index of the closest vertex found.
   *
   * @dealiiConceptRequires{
   *   concepts::is_triangulation_or_dof_handler<MeshType<dim, spacedim>>}
   */
  template <int dim, template <int, int> class MeshType, int spacedim>
  DEAL_II_CXX20_REQUIRES(
    (concepts::is_triangulation_or_dof_handler<MeshType<dim, spacedim>>))
  unsigned int find_closest_vertex(
    const Mapping<dim, spacedim>  &mapping,
    const MeshType<dim, spacedim> &mesh,
    const Point<spacedim>         &p,
    const std::vector<bool>       &marked_vertices = {});


  /**
   * Find and return a vector of iterators to active cells that surround a
   * given vertex with index @p vertex_index.
   *
   * For locally refined grids, the vertex itself might not be a vertex of all
   * adjacent cells that are returned. However, it will always be either a
   * vertex of a cell or be a hanging node located on a face or an edge of it.
   *
   * @param container A variable of a type that satisfies the requirements of
   * the
   * @ref ConceptMeshType "MeshType concept".
   * @param vertex_index The index of the vertex for which we try to find
   * adjacent cells.
   * @return A vector of cells that lie adjacent to the given vertex.
   *
   * @note It isn't entirely clear at this time whether the function does the
   * right thing with anisotropically refined meshes. It needs to be checked
   * for this case.
   *
   * @dealiiConceptRequires{
   *   concepts::is_triangulation_or_dof_handler<MeshType>}
   */
  template <typename MeshType>
  DEAL_II_CXX20_REQUIRES((concepts::is_triangulation_or_dof_handler<MeshType>))
#ifndef _MSC_VER
  std::vector<typename MeshType::active_cell_iterator>
#else
  std::vector<
    typename dealii::internal::ActiveCellIterator<MeshType::dimension,
                                                  MeshType::space_dimension,
                                                  MeshType>::type>
#endif
    find_cells_adjacent_to_vertex(const MeshType    &container,
                                  const unsigned int vertex_index);

  /**
   * Find an active non-artificial cell that surrounds a given point @p p. The return type
   * is a pair of an iterator to the active cell along with the unit cell
   * coordinates of the point.
   *
   * The algorithm used in this function proceeds by first looking for the
   * vertex located closest to the given point, see
   * GridTools::find_closest_vertex(). Secondly, all adjacent cells to this
   * vertex are found in the mesh, see
   * GridTools::find_cells_adjacent_to_vertex(). Lastly, for each of these
   * cells, the function tests whether the point is inside. This check is
   * performed using the given @p mapping argument to determine whether cells
   * have straight or curved boundaries.
   *
   * If a point lies on the boundary of two or more cells, then the algorithm
   * tries to identify the cell that is of highest refinement level.
   *
   * If the point requested does not lie in a locally-owned or ghost cell,
   * then this function will return the (invalid) MeshType<dim, spacedim>::end()
   * iterator. This case can be handled similarly to the various `std::find()`
   * and `std::lower_bound()` functions.
   *
   * @param mapping The mapping used to determine whether the given point is
   *   inside a given cell.
   * @param mesh A variable of a type that satisfies the requirements of the
   *   @ref ConceptMeshType "MeshType concept".
   * @param p The point for which we want to find the surrounding cell.
   * @param marked_vertices An array of `bool`s indicating whether an
   *   entry in the vertex array should be considered
   *   (and the others must be ignored) as the potentially
   *   closest vertex to the specified point. On specifying a non-default
   *   @p marked_vertices, find_closest_vertex() would
   *   only search among @p marked_vertices for the closest vertex.
   *   The size of this array should be equal to n_vertices() of the
   *   triangulation (as opposed to n_used_vertices() ). The motivation of using
   *   @p marked_vertices is to cut down the search space of vertices if one has
   *   a priori knowledge of a collection of vertices that the point of interest
   *   may be close to.
   * @param tolerance Tolerance in terms of unit cell coordinates. Depending
   *   on the problem, it might be necessary to adjust the tolerance in order
   *   to be able to identify a cell. Floating
   *   point arithmetic implies that a point will, in general, not lie exactly
   *   on a vertex, edge, or face. In either case, it is not predictable which
   *   of the cells adjacent to a vertex or an edge/face this function returns.
   *   Consequently, algorithms that call this function need to take into
   *   account that the returned cell will only contain the point approximately.
   *
   * @return A pair of an iterators into the mesh that points to the
   * surrounding cell, and of the unit cell coordinates of that point. This
   * local position might be located slightly outside an actual unit cell,
   * due to numerical roundoff. Therefore, the point returned by this function
   * should be projected onto the unit cell, using
   * ReferenceCell::closest_point(). This is not automatically performed
   * by the algorithm. The returned cell can be a locally-owned cell or a
   * ghost cell (but not an artificial cell). The returned cell might be a
   * ghost cell even if the given point is a vertex of a locally owned cell.
   * The reason behind is that this is the only way to guarantee that all
   * processors that participate in a parallel triangulation will agree which
   * cell contains a point. For example, if two processors come together
   * at one vertex and the function is called with this vertex, then one
   * processor will return a locally owned cell and the other one a ghost cell.
   *
   * @dealiiConceptRequires{
   *   concepts::is_triangulation_or_dof_handler<MeshType<dim, spacedim>>}
   */
  template <int dim, template <int, int> class MeshType, int spacedim>
  DEAL_II_CXX20_REQUIRES(
    (concepts::is_triangulation_or_dof_handler<MeshType<dim, spacedim>>))
#ifndef _MSC_VER
  std::pair<typename MeshType<dim, spacedim>::active_cell_iterator, Point<dim>>
#else
  std::pair<typename dealii::internal::
              ActiveCellIterator<dim, spacedim, MeshType<dim, spacedim>>::type,
            Point<dim>>
#endif
    find_active_cell_around_point(const Mapping<dim, spacedim>  &mapping,
                                  const MeshType<dim, spacedim> &mesh,
                                  const Point<spacedim>         &p,
                                  const std::vector<bool> &marked_vertices = {},
                                  const double             tolerance = 1.e-10);

  /**
   * A version of the above function that assumes straight boundaries and
   * as a consequence simply calls the above function using MappingQ1 for
   * the mapping argument.
   *
   * @return An iterator into the mesh that points to the surrounding cell.
   *
   * @dealiiConceptRequires{
   *   concepts::is_triangulation_or_dof_handler<MeshType<dim, spacedim>>}
   */
  template <int dim, template <int, int> class MeshType, int spacedim>
  DEAL_II_CXX20_REQUIRES(
    (concepts::is_triangulation_or_dof_handler<MeshType<dim, spacedim>>))
#ifndef _MSC_VER
  typename MeshType<dim, spacedim>::active_cell_iterator
#else
  typename dealii::internal::
    ActiveCellIterator<dim, spacedim, MeshType<dim, spacedim>>::type
#endif
    find_active_cell_around_point(const MeshType<dim, spacedim> &mesh,
                                  const Point<spacedim>         &p,
                                  const std::vector<bool> &marked_vertices = {},
                                  const double             tolerance = 1.e-10);

  /**
   * Another version where we use that mapping on a given
   * cell that corresponds to the active finite element index of that cell.
   * This is obviously only useful for hp-problems, since the active finite
   * element index for all other DoF handlers is always zero.
   */
  template <int dim, int spacedim>
  std::pair<typename DoFHandler<dim, spacedim>::active_cell_iterator,
            Point<dim>>
  find_active_cell_around_point(
    const hp::MappingCollection<dim, spacedim> &mapping,
    const DoFHandler<dim, spacedim>            &mesh,
    const Point<spacedim>                      &p,
    const double                                tolerance = 1.e-10);

  /**
   * Finding an active non-artificial cell around a point can be very expensive
   * in terms of computational costs. This function aims at providing a fast
   * version of the above functions by using a space-tree to speed up the
   * geometry search.
   *
   * @param cache Object with information about the space-tree of a triangulation,
   * see GridTools::Cache.
   * @param p The point for which we want to find the surrounding cell.
   * @param cell_hint Gives a hint for the geometry search, which is beneficial
   * if a-priori knowledge is available regarding the cell on which the point
   * may likely be located. A typical use case would be that this search has
   * to be done for an array of points that are close to each other and where
   * the adjacent cell of the previous point is a good hint for the next point
   * in the array.
   * @param marked_vertices See above.
   * @param tolerance See above.
   *
   *
   * The following code example shows how to use this function:
   *
   * @code
   * GridTools::Cache<dim, dim> cache(triangulation, mapping);
   * auto cell_hint = typename Triangulation<dim, dim>::active_cell_iterator();
   * std::vector<bool> marked_vertices = {};
   * double tolerance = 1.e-10;
   *
   * std::vector<Point<dim>> points; // a vector of many points
   * ...
   *
   * for(auto p : points)
   * {
   *   auto cell_and_ref_point = GridTools::find_active_cell_around_point(
   *     cache, p, cell_hint, marked_vertices, tolerance);
   *
   *   if (cell_and_ref_point.first != triangulation.end())
   *     {
   *      // use current cell as hint for the next point
   *      cell_hint = cell_and_ref_point.first;
   *      // do something with cell_and_ref_point
   *      ...
   *   }
   *  else
   *    {
   *       // The function did not find a locally owned or ghost cell in which
   *       // the point is located. We ought to handle this somehow here.
   *    }
   *   ...
   * }
   * @endcode
   */
  template <int dim, int spacedim>
  std::pair<typename Triangulation<dim, spacedim>::active_cell_iterator,
            Point<dim>>
  find_active_cell_around_point(
    const Cache<dim, spacedim> &cache,
    const Point<spacedim>      &p,
    const typename Triangulation<dim, spacedim>::active_cell_iterator &
      cell_hint = typename Triangulation<dim, spacedim>::active_cell_iterator(),
    const std::vector<bool> &marked_vertices = {},
    const double             tolerance       = 1.e-10);

  /**
   * A version of the previous function that exploits an already existing
   * map between vertices and cells (constructed using the function
   * GridTools::vertex_to_cell_map()), a map of vertex_to_cell_centers (obtained
   * through GridTools::vertex_to_cell_centers_directions()), and
   * optionally an RTree constructed from the used vertices of the
   * Triangulation.
   *
   * @note All of these structures can be queried from a
   * GridTools::Cache object. Note, however, that in this case MeshType
   * has to be Triangulation, so that it might be more appropriate to directly
   * call the function above with argument `cache` in this case.
   *
   * @dealiiConceptRequires{
   *   concepts::is_triangulation_or_dof_handler<MeshType<dim, spacedim>>}
   */
  template <int dim, template <int, int> class MeshType, int spacedim>
  DEAL_II_CXX20_REQUIRES(
    (concepts::is_triangulation_or_dof_handler<MeshType<dim, spacedim>>))
#ifndef _MSC_VER
  std::pair<typename MeshType<dim, spacedim>::active_cell_iterator, Point<dim>>
#else
  std::pair<typename dealii::internal::
              ActiveCellIterator<dim, spacedim, MeshType<dim, spacedim>>::type,
            Point<dim>>
#endif
    find_active_cell_around_point(
      const Mapping<dim, spacedim>  &mapping,
      const MeshType<dim, spacedim> &mesh,
      const Point<spacedim>         &p,
      const std::vector<
        std::set<typename MeshType<dim, spacedim>::active_cell_iterator>>
        &vertex_to_cell_map,
      const std::vector<std::vector<Tensor<1, spacedim>>>
        &vertex_to_cell_centers,
      const typename MeshType<dim, spacedim>::active_cell_iterator &cell_hint =
        typename MeshType<dim, spacedim>::active_cell_iterator(),
      const std::vector<bool> &marked_vertices = {},
      const RTree<std::pair<Point<spacedim>, unsigned int>> &
        used_vertices_rtree = RTree<std::pair<Point<spacedim>, unsigned int>>{},
      const double tolerance = 1.e-10,
      const RTree<
        std::pair<BoundingBox<spacedim>,
                  typename Triangulation<dim, spacedim>::active_cell_iterator>>
        *relevant_cell_bounding_boxes_rtree = nullptr);

  /**
   * As compared to the functions above, this function identifies all active
   * non-artificial cells around a point for a given tolerance level `tolerance`
   * in terms of unit coordinates. Given a first cell with reference coordinates
   * as parameter
   * @p first_cell, e.g. obtained by one of the functions above, all
   * corresponding neighboring cells with points in unit coordinates are also
   * identified.
   *
   * The parameter @p vertex_to_cells allows to accelerate the process of
   * identifying the neighbors of a cell, by first precomputing a map from the
   * vertex indices to the cells. Such data structure is, e.g., provided by
   * GridTools::Cache::get_vertex_to_cell_map().
   *
   * This function is useful e.g. for discontinuous function spaces where, for
   * the case the given point `p` lies on a vertex, edge or face, several
   * cells might hold independent values of the solution that get combined in
   * some way in a user code.
   *
   * This function is used as follows
   * @code
   *   auto first_pair = GridTools::find_active_cell_around_point(...);
   *   auto all_cells  = GridTools::find_all_active_cells_around_point(
   *   			   mapping, mesh, p, tolerance, first_pair);
   * @endcode
   *
   * @dealiiConceptRequires{
   *   concepts::is_triangulation_or_dof_handler<MeshType<dim, spacedim>>}
   */
  template <int dim, template <int, int> class MeshType, int spacedim>
  DEAL_II_CXX20_REQUIRES(
    (concepts::is_triangulation_or_dof_handler<MeshType<dim, spacedim>>))
#ifndef _MSC_VER
  std::vector<std::pair<typename MeshType<dim, spacedim>::active_cell_iterator,
                        Point<dim>>>
#else
  std::vector<std::pair<
    typename dealii::internal::
      ActiveCellIterator<dim, spacedim, MeshType<dim, spacedim>>::type,
    Point<dim>>>
#endif
    find_all_active_cells_around_point(
      const Mapping<dim, spacedim>  &mapping,
      const MeshType<dim, spacedim> &mesh,
      const Point<spacedim>         &p,
      const double                   tolerance,
      const std::pair<typename MeshType<dim, spacedim>::active_cell_iterator,
                      Point<dim>>   &first_cell,
      const std::vector<
        std::set<typename MeshType<dim, spacedim>::active_cell_iterator>>
        *vertex_to_cells = nullptr);

  /**
   * A variant of the previous function that internally calls one of the
   * functions find_active_cell_around_point() to obtain a first cell, and
   * subsequently adds all other active non-artificial cells by calling the
   * function find_all_active_cells_around_point() above.
   *
   * @dealiiConceptRequires{
   *   concepts::is_triangulation_or_dof_handler<MeshType<dim, spacedim>>}
   */
  template <int dim, template <int, int> class MeshType, int spacedim>
  DEAL_II_CXX20_REQUIRES(
    (concepts::is_triangulation_or_dof_handler<MeshType<dim, spacedim>>))
#ifndef _MSC_VER
  std::vector<std::pair<typename MeshType<dim, spacedim>::active_cell_iterator,
                        Point<dim>>>
#else
  std::vector<std::pair<
    typename dealii::internal::
      ActiveCellIterator<dim, spacedim, MeshType<dim, spacedim>>::type,
    Point<dim>>>
#endif
    find_all_active_cells_around_point(
      const Mapping<dim, spacedim>  &mapping,
      const MeshType<dim, spacedim> &mesh,
      const Point<spacedim>         &p,
      const double                   tolerance       = 1e-10,
      const std::vector<bool>       &marked_vertices = {});

  /**
   * Return a list of all descendants of the given cell that are active. For
   * example, if the current cell is once refined but none of its children are
   * any further refined, then the returned list will contain all its
   * children.
   *
   * If the current cell is already active, then the returned list is empty
   * (because the cell has no children that may be active).
   *
   * @tparam MeshType A type that satisfies the requirements of the
   * @ref ConceptMeshType "MeshType concept".
   * @param cell An iterator pointing to a cell of the mesh.
   * @return A list of active descendants of the given cell
   *
   * @note Since in C++ the MeshType template argument can not be deduced from
   * a function call, you will have to specify it after the function name, as
   * for example in
   * @code
   *   GridTools::get_active_child_cells<DoFHandler<dim> > (cell)
   * @endcode
   *
   * @dealiiConceptRequires{concepts::is_triangulation_or_dof_handler<MeshType>}
   */
  template <typename MeshType>
  DEAL_II_CXX20_REQUIRES(concepts::is_triangulation_or_dof_handler<MeshType>)
  std::vector<typename MeshType::active_cell_iterator> get_active_child_cells(
    const typename MeshType::cell_iterator &cell);

  /**
   * Extract the active cells around a given cell @p cell and return them in
   * the vector @p active_neighbors. These neighbors are specifically the
   * <i>face</i> neighbors of a cell or, if that neighbor is further
   * refined, its active children that border on that face. On the other
   * hand, the neighbors returned do not include cells that lie, for
   * example, diagonally opposite to a vertex but are not face neighbors
   * themselves. (In 3d, it also does not include cells that are
   * adjacent to one of the edges of the current cell, but are not
   * face neighbors.)
   *
   * @tparam MeshType A type that satisfies the requirements of the
   * @ref ConceptMeshType "MeshType concept".
   * @param[in] cell An iterator pointing to a cell of the mesh.
   * @param[out] active_neighbors A list of active descendants of the given
   * cell
   *
   * @note Since in C++ the MeshType template argument can not be deduced from
   * a function call, you will have to specify it after the function name, as
   * for example in
   * @code
   *   GridTools::get_active_neighbors<DoFHandler<dim>>(cell, active_neighbors)
   * @endcode
   *
   * @dealiiConceptRequires{concepts::is_triangulation_or_dof_handler<MeshType>}
   */
  template <typename MeshType>
  DEAL_II_CXX20_REQUIRES(concepts::is_triangulation_or_dof_handler<MeshType>)
  void get_active_neighbors(
    const typename MeshType::active_cell_iterator        &cell,
    std::vector<typename MeshType::active_cell_iterator> &active_neighbors);

  /**
   * Extract and return the active cell layer around a subdomain (set of
   * active cells) in the @p mesh (i.e. those that share a common set of
   * vertices with the subdomain but are not a part of it). Here, the
   * "subdomain" consists of exactly all of those cells for which the @p
   * predicate returns @p true.
   *
   * An example of a custom predicate is one that checks for a given material
   * id
   * @code
   * template <int dim>
   * bool
   * pred_mat_id(const typename Triangulation<dim>::active_cell_iterator & cell)
   * {
   *   return cell->material_id() ==  1;
   * }
   * @endcode
   * and we can then extract the layer of cells around this material with the
   * following call:
   * @code
   * GridTools::compute_active_cell_halo_layer(tria, pred_mat_id<dim>);
   * @endcode
   *
   * Predicates that are frequently useful can be found in namespace
   * IteratorFilters. For example, it is possible to extract a layer
   * of cells around all of those cells with a given material id,
   * @code
   * GridTools::compute_active_cell_halo_layer(
   *   tria, IteratorFilters::MaterialIdEqualTo(1, true));
   * @endcode
   * or around all cells with one of a set of active FE indices for a DoFHandler
   * with hp-capabilities
   * @code
   * GridTools::compute_active_cell_halo_layer(
   *   hp_dof_handler, IteratorFilters::ActiveFEIndexEqualTo({1,2}, true));
   * @endcode
   * Note that in the last two examples we ensure that the predicate returns
   * true only for locally owned cells. This means that the halo layer will
   * not contain any artificial cells.
   *
   * @tparam MeshType A type that satisfies the requirements of the
   * @ref ConceptMeshType "MeshType concept".
   * @param[in] mesh A mesh (i.e. objects of type Triangulation or DoFHandler).
   * @param[in] predicate A function  (or object of a type with an operator())
   * defining the subdomain around which the halo layer is to be extracted. It
   * is a function that takes in an active cell and returns a boolean.
   * @return A list of active cells sharing at least one common vertex with
   * the predicated subdomain.
   *
   * @dealiiConceptRequires{concepts::is_triangulation_or_dof_handler<MeshType>}
   */
  template <typename MeshType>
  DEAL_II_CXX20_REQUIRES(concepts::is_triangulation_or_dof_handler<MeshType>)
  std::
    vector<typename MeshType::active_cell_iterator> compute_active_cell_halo_layer(
      const MeshType &mesh,
      const std::function<bool(const typename MeshType::active_cell_iterator &)>
        &predicate);


  /**
   * Extract and return the cell layer around a subdomain (set of
   * cells) on a specified level of the @p mesh (i.e. those cells on
   * that level that share a common set of vertices with the subdomain
   * but are not a part of it). Here, the "subdomain" consists of exactly
   * all of those cells for which the @p predicate returns @p true.
   *
   * @dealiiConceptRequires{concepts::is_triangulation_or_dof_handler<MeshType>}
   */
  template <typename MeshType>
  DEAL_II_CXX20_REQUIRES(concepts::is_triangulation_or_dof_handler<MeshType>)
  std::
    vector<typename MeshType::cell_iterator> compute_cell_halo_layer_on_level(
      const MeshType &mesh,
      const std::function<bool(const typename MeshType::cell_iterator &)>
                        &predicate,
      const unsigned int level);


  /**
   * Extract and return ghost cells which are the active cell layer around all
   * locally owned cells. This is most relevant for
   * parallel::shared::Triangulation where it will return a subset of all
   * ghost cells on a processor, but for parallel::distributed::Triangulation
   * this will return all the ghost cells.
   *
   * @tparam MeshType A type that satisfies the requirements of the
   * @ref ConceptMeshType "MeshType concept".
   * @param[in] mesh A mesh (i.e. objects of type Triangulation or DoFHandler).
   * @return A list of ghost cells
   *
   * @dealiiConceptRequires{concepts::is_triangulation_or_dof_handler<MeshType>}
   */
  template <typename MeshType>
  DEAL_II_CXX20_REQUIRES(concepts::is_triangulation_or_dof_handler<MeshType>)
  std::vector<
    typename MeshType::
      active_cell_iterator> compute_ghost_cell_halo_layer(const MeshType &mesh);

  /**
   * Extract and return the set of active cells within a geometric distance of
   * @p layer_thickness around a subdomain (set of active cells) in the @p mesh.
   * Here, the "subdomain" consists of exactly all of
   * those cells for which the @p predicate returns @p true.
   *
   * The function first computes the cells that form the 'surface' of the
   * subdomain that consists of all of the active cells for which the predicate
   * is true. Using compute_bounding_box(), a bounding box is
   * computed for this subdomain and extended by @p layer_thickness. These
   * cells are called interior subdomain boundary cells.
   * The active cells with all of their vertices outside the extended
   * bounding box are ignored.
   * The cells that are inside the extended bounding box are then checked for
   * their proximity to the interior subdomain boundary cells. This implies
   * checking the distance between a pair of arbitrarily oriented cells,
   * which is not trivial in general. To simplify this, the algorithm checks
   * the distance between the two enclosing spheres of the cells.
   * This will definitely result in slightly more cells being marked but
   * also greatly simplifies the arithmetic complexity of the algorithm.
   *
   * @image html active_cell_layer_within_distance.png
   * The image shows a mesh generated by subdivided_hyper_rectangle(). The cells
   * are marked using three different colors. If the grey colored cells in the
   * image are the cells for which the predicate is true, then the function
   * compute_active_cell_layer_within_distance() will return a set of cell
   * iterators corresponding to the cells colored in red.
   * The red colored cells are the active cells that are within a given
   * distance to the grey colored cells.
   *
   * @tparam MeshType A type that satisfies the requirements of the
   * @ref ConceptMeshType "MeshType concept".
   * @param mesh A mesh (i.e. objects of type Triangulation or DoFHandler).
   * @param predicate A function  (or object of a type with an operator())
   * defining the subdomain around which the halo layer is to be extracted. It
   * is a function that takes in an active cell and returns a boolean.
   * @param layer_thickness specifies the geometric distance within
   * which the function searches for active cells from the predicate domain.
   * If the minimal distance between the enclosing sphere of the an
   * active cell and the enclosing sphere of any of the cells for which
   * the @p predicate returns @p true is less than @p layer_thickness,
   * then the active cell is an \a active_cell_within_distance.
   * @return A list of active cells within a given geometric distance
   * @p layer_thickness from the set of active cells for which the @p predicate
   * returns @p true.
   *
   * See compute_active_cell_halo_layer().
   *
   * @dealiiConceptRequires{concepts::is_triangulation_or_dof_handler<MeshType>}
   */
  template <typename MeshType>
  DEAL_II_CXX20_REQUIRES(concepts::is_triangulation_or_dof_handler<MeshType>)
  std::
    vector<typename MeshType::active_cell_iterator> compute_active_cell_layer_within_distance(
      const MeshType &mesh,
      const std::function<bool(const typename MeshType::active_cell_iterator &)>
                  &predicate,
      const double layer_thickness);

  /**
   * Extract and return a set of ghost cells which are within a
   * @p layer_thickness around all locally owned cells.
   * This is most relevant for parallel::shared::Triangulation
   * where it will return a subset of all ghost cells on a process, but for
   * parallel::distributed::Triangulation this will return all the ghost cells.
   * All the cells for the parallel::shared::Triangulation class that
   * are not owned by the current processor can be considered as ghost cells;
   * in particular, they do not only form a single layer of cells around the
   * locally owned ones.
   *
   * @tparam MeshType A type that satisfies the requirements of the
   * @ref ConceptMeshType "MeshType concept".
   * @param mesh A mesh (i.e. objects of type Triangulation or DoFHandler).
   * @param layer_thickness specifies the geometric distance within
   * which the function searches for active cells from the locally owned cells.
   * @return A subset of ghost cells within a given geometric distance of @p
   * layer_thickness from the locally owned cells of a current process.
   *
   * Also see compute_ghost_cell_halo_layer() and
   * compute_active_cell_layer_within_distance().
   *
   * @dealiiConceptRequires{concepts::is_triangulation_or_dof_handler<MeshType>}
   */
  template <typename MeshType>
  DEAL_II_CXX20_REQUIRES(concepts::is_triangulation_or_dof_handler<MeshType>)
  std::
    vector<typename MeshType::active_cell_iterator> compute_ghost_cell_layer_within_distance(
      const MeshType &mesh,
      const double    layer_thickness);

  /**
   * Compute a collection of bounding boxes so that all active cells for which
   * the given predicate is true, are completely enclosed in at least one of the
   * bounding boxes. Notice the cover is only guaranteed to contain all these
   * active cells but it's not necessarily exact i.e. it can include a bigger
   * area than their union. (This is of course unavoidable in any case if cells
   * are not rectangular or brick-shaped, but it is also true if cells are
   * since it is inefficient to create as many bounding boxes as there are
   * cells; rather, the algorithm here tries to combine the bounding boxes
   * of multiple cells into a cheaper representation, at the cost of a set
   * of bounding boxes that may be larger than the union of the cells -- see
   * the description of the relevant function arguments below.)
   *
   * For each cell at a given refinement level containing active cells for which @p predicate is true,
   * the function creates a bounding box of its children for which @p predicate is true.
   *
   * This results in a cover of all active cells for which @p predicate is true; the parameters
   * @p allow_merge and @p max_boxes are used to reduce the number of cells at a computational cost and
   * covering a bigger n-dimensional volume.
   *
   * @param[in] mesh The mesh object this function is to work on. This
   *   is generally either a triangulation of some kind, or a DoFHandler
   *   object.
   * @param[in] predicate A function-like object that returns true or
   *   false depending on whether the property of the cells to enclose
   *   is satisfied. An example is IteratorFilters::LocallyOwnedCell,
   *   but it can also be a lambda function or anything else that can be
   *   called with a cell as argument.
   *   This predicate is tested only on active cells.
   * @param[in] refinement_level Defines the level at which the
   *  initial bounding boxes are created. The refinement should be set
   *  to a coarse refinement level. A bounding box is created for each
   *  active cell at a coarser level than @p refinement_level; if @p
   *  refinement_level is higher than the number of levels of the
   *  triangulation an exception is thrown.
   * @param[in] allow_merge This flag allows for box merging and, by
   *  default, is false. The algorithm has a cost of O(N^2) where N is
   *  the number of the bounding boxes created from the refinement
   *  level; for this reason, if the flag is set to true, make sure to
   *  choose wisely a coarse enough @p refinement_level.
   * @param[in] max_boxes The maximum number of bounding boxes to
   *  compute. If more are created the smaller ones are merged with
   *  neighbors. By default after merging the boxes which can be
   *  expressed as a single one no more boxes are merged. See the
   *  BoundingBox::get_neighbor_type () function for details.  Notice
   *  only neighboring cells are merged (see the @p get_neighbor_type
   *  function in bounding box class): if the target number of
   *  bounding boxes max_boxes can't be reached by merging neighbors
   *  an exception is thrown.
   *
   * The following image describes an example of the algorithm with @p
   * refinement_level = 2, @p allow_merge = true and @p max_boxes =
   * 1. The cells with the property predicate are in red, the area of
   * a bounding box is slightly orange.
   * @image html bounding_box_predicate.png
   *
   * - 1. In black we can see the cells of the current level.
   * - 2. For each cell containing the red area a bounding box is created: by
   * default these are returned.
   * - 3. Because @p allow_merge = true the number of bounding boxes is reduced while not changing the cover.
   *  If @p max_boxes was left as default or bigger than 1 these two boxes would be returned.
   * - 4. Because @p max_boxes = 1 the smallest bounding box is merged to the bigger one.
   * Notice it is important to choose the parameters wisely. For instance, @p allow_merge = false and
   * @p refinement_level = 1 returns the very same bounding box but with a
   * fraction of the computational cost.
   *
   * This function does not take into account the curvature of cells and thus it
   * is not suited for handling curved geometry: the mapping is assumed to be
   * linear.
   *
   * @dealiiConceptRequires{concepts::is_triangulation_or_dof_handler<MeshType>}
   */
  template <typename MeshType>
  DEAL_II_CXX20_REQUIRES(concepts::is_triangulation_or_dof_handler<MeshType>)
  std::
    vector<BoundingBox<MeshType::space_dimension>> compute_mesh_predicate_bounding_box(
      const MeshType &mesh,
      const std::function<bool(const typename MeshType::active_cell_iterator &)>
                        &predicate,
      const unsigned int refinement_level = 0,
      const bool         allow_merge      = false,
      const unsigned int max_boxes        = numbers::invalid_unsigned_int);

  /**
   * Given an array of points, use the global bounding box description obtained
   * using GridTools::compute_mesh_predicate_bounding_box to guess, for each of
   * them, which process might own it.
   *
   * @param[in] global_bboxes Vector of bounding boxes describing the portion of
   *  mesh with a property for each process.
   * @param[in] points Array of points to test.
   *
   * @return A tuple containing the following information:
   *  - A vector indicized with ranks of processes. For each rank it contains
   *   a vector of the indices of points it might own.
   *  - A map from the index <code>unsigned int</code> of the point in @p points
   *   to the rank of the owner.
   *  - A map from the index <code>unsigned int</code> of the point in @p points
   *   to the ranks of the guessed owners.
   *
   * @note The actual return type of this function, i.e., the type referenced
   * above as @p return_type, is
   * @code
   * std::tuple<std::vector<std::vector<unsigned int>>,
   *            std::map< unsigned int, unsigned int>,
   *            std::map< unsigned int, std::vector<unsigned int>>>
   * @endcode
   * The type is abbreviated in the online documentation to improve readability
   * of this page.
   */
  template <int spacedim>
#ifndef DOXYGEN
  std::tuple<std::vector<std::vector<unsigned int>>,
             std::map<unsigned int, unsigned int>,
             std::map<unsigned int, std::vector<unsigned int>>>
#else
  return_type
#endif
  guess_point_owner(
    const std::vector<std::vector<BoundingBox<spacedim>>> &global_bboxes,
    const std::vector<Point<spacedim>>                    &points);


  /**
   * Given a covering rtree (see GridTools::Cache::get_covering_rtree()), and an
   * array of points, find a superset of processes which, individually,
   * may own the cell containing the points.
   *
   * For further details see GridTools::guess_point_owner; here only
   * different input/output types are reported:
   *
   * @param[in] covering_rtree RTRee which enables us to identify which
   * process(es) in a parallel computation may own the cell that
   * surrounds a given point.
   *
   * @param[in] points A vector of points to consider.
   *
   * @return A tuple containing the following information:
   *  - A map indexed by processor ranks. For each rank it contains
   *   a vector of the indices of points it might own.
   *  - A map from the index <code>unsigned int</code> of the point in @p points
   *   to the rank of the owner; these are points for which a single possible
   *   owner was found.
   *  - A map from the index <code>unsigned int</code> of the point in @p points
   *   to the ranks of the guessed owners; these are points for which multiple
   *   possible owners were found.
   *
   * @note The actual return type of this function, i.e., the type referenced
   * above as @p return_type, is
   * @code
   * std::tuple<std::map<unsigned int, std::vector<unsigned int>>,
   *            std::map<unsigned int, unsigned int>,
   *            std::map<unsigned int, std::vector<unsigned int>>>
   * @endcode
   * The type is abbreviated in the online documentation to improve readability
   * of this page.
   */
  template <int spacedim>
#ifndef DOXYGEN
  std::tuple<std::map<unsigned int, std::vector<unsigned int>>,
             std::map<unsigned int, unsigned int>,
             std::map<unsigned int, std::vector<unsigned int>>>
#else
  return_type
#endif
  guess_point_owner(
    const RTree<std::pair<BoundingBox<spacedim>, unsigned int>> &covering_rtree,
    const std::vector<Point<spacedim>>                          &points);

  /**
   * Return a vector of normalized tensors for each vertex-cell combination of
   * the output of GridTools::vertex_to_cell_map() (which is expected as input
   * parameter for this function). Each tensor represents a geometric vector
   * from the vertex to the respective cell center.
   *
   * An assertion will be thrown if the size of the input vector is not equal to
   * the number of vertices of the triangulation.
   *
   * result[v][c] is a unit Tensor for vertex index v, indicating the direction
   * of the center of the c-th cell with respect to the vertex v.
   */
  template <int dim, int spacedim>
  std::vector<std::vector<Tensor<1, spacedim>>>
  vertex_to_cell_centers_directions(
    const Triangulation<dim, spacedim> &mesh,
    const std::vector<
      std::set<typename Triangulation<dim, spacedim>::active_cell_iterator>>
      &vertex_to_cells);


  /**
   * Return the local vertex index of cell @p cell that is closest to
   * the given location @p position. The location of the vertices is extracted
   * from the (optional) @p mapping argument, to guarantee that the correct
   * answer is returned when the underlying mapping modifies the position of the
   * vertices.
   */
  template <int dim, int spacedim>
  unsigned int
  find_closest_vertex_of_cell(
    const typename Triangulation<dim, spacedim>::active_cell_iterator &cell,
    const Point<spacedim>                                             &position,
    const Mapping<dim, spacedim>                                      &mapping =
      (ReferenceCells::get_hypercube<dim>()
#ifndef _MSC_VER
         .template get_default_linear_mapping<dim, spacedim>()
#else
         .ReferenceCell::get_default_linear_mapping<dim, spacedim>()
#endif
         ));

  /**
   * Compute a globally unique index for each vertex and hanging node
   * associated with a locally owned active cell. The vertices of a ghost cell
   * that are hanging nodes of a locally owned cells have a global index.
   * However, the other vertices of the cells that do not <i>touch</i> an
   * active cell do not have a global index on this processor.
   *
   * The key of the map is the local index of the vertex and the value is the
   * global index. The indices need to be recomputed after refinement or
   * coarsening and may be different.
   */
  template <int dim, int spacedim>
  std::map<unsigned int, types::global_vertex_index>
  compute_local_to_global_vertex_index_map(
    const Triangulation<dim, spacedim> &triangulation);

  /** @} */
  /**
   * @name Partitions and subdomains of triangulations
   */
  /** @{ */

  /**
   * Use graph partitioner to partition the active cells making up the entire
   * domain. After calling this function, the subdomain ids of all active cells
   * will have values
   * between zero and @p n_partitions-1. You can access the subdomain id of a cell by using
   * <tt>cell-@>subdomain_id()</tt>.
   *
   * Use the third argument to select between partitioning algorithms provided
   * by METIS or ZOLTAN. METIS is the default partitioner.
   *
   * If deal.II was not installed with ZOLTAN or METIS, this function will
   * generate an error
   * when the respective partition method is chosen, unless @p n_partitions is one.
   * I.e., you can write a program so that it runs in the single-processor
   * single-partition case without packages installed, and only requires them
   * installed when multiple partitions are required.
   *
   * @note If the `weight` signal has been attached to the @p triangulation,
   * then this will be used and passed to the partitioner.
   */
  template <int dim, int spacedim>
  void
  partition_triangulation(const unsigned int               n_partitions,
                          Triangulation<dim, spacedim>    &triangulation,
                          const SparsityTools::Partitioner partitioner =
                            SparsityTools::Partitioner::metis);

  /**
   * This function performs the same operation as the one above, except that
   * it takes into consideration a specific set of @p cell_weights, which allow the
   * partitioner to balance the graph while taking into consideration the
   * computational effort expended on each cell.
   *
   * @note If the @p cell_weights vector is empty, then no weighting is taken
   * into consideration. If not then the size of this vector must equal to the
   * number of active cells in the triangulation.
   */
  template <int dim, int spacedim>
  void
  partition_triangulation(const unsigned int               n_partitions,
                          const std::vector<unsigned int> &cell_weights,
                          Triangulation<dim, spacedim>    &triangulation,
                          const SparsityTools::Partitioner partitioner =
                            SparsityTools::Partitioner::metis);

  /**
   * This function does the same as the previous one, i.e. it partitions a
   * triangulation using a partitioning algorithm into a number of subdomains
   * identified by the <code>cell-@>subdomain_id()</code> flag.
   *
   * The difference to the previous function is the second argument, a
   * sparsity pattern that represents the connectivity pattern between cells.
   *
   * While the function above builds it directly from the triangulation by
   * considering which cells neighbor each other, this function can take a
   * more refined connectivity graph. The sparsity pattern needs to be of size
   * $N\times N$, where $N$ is the number of active cells in the
   * triangulation. If the sparsity pattern contains an entry at position
   * $(i,j)$, then this means that cells $i$ and $j$ (in the order in which
   * they are traversed by active cell iterators) are to be considered
   * connected; partitioning algorithm will then try to partition the domain in
   * such a way that (i) the subdomains are of roughly equal size, and (ii) a
   * minimal number of connections are broken.
   *
   * This function is mainly useful in cases where connections between cells
   * exist that are not present in the triangulation alone (otherwise the
   * previous function would be the simpler one to use). Such connections may
   * include that certain parts of the boundary of a domain are coupled
   * through symmetric boundary conditions or integrals (e.g. friction contact
   * between the two sides of a crack in the domain), or if a numerical scheme
   * is used that not only connects immediate neighbors but a larger
   * neighborhood of cells (e.g. when solving integral equations).
   *
   * In addition, this function may be useful in cases where the default
   * sparsity pattern is not entirely sufficient. This can happen because the
   * default is to just consider face neighbors, not neighboring cells that
   * are connected by edges or vertices. While the latter couple when using
   * continuous finite elements, they are typically still closely connected in
   * the neighborship graph, and partitioning algorithm
   * will not usually cut important connections in this case. However, if there
   * are vertices in the mesh where many cells (many more than the common 4 or 6
   * in 2d and 3d, respectively) come together, then there will be a significant
   * number of cells that are connected across a vertex, but several degrees
   * removed in the connectivity graph built only using face neighbors. In a
   * case like this, partitioning algorithm may sometimes make bad decisions and
   * you may want to build your own connectivity graph.
   *
   * @note If the `weight` signal has been attached to the @p triangulation,
   * then this will be used and passed to the partitioner.
   */
  template <int dim, int spacedim>
  void
  partition_triangulation(const unsigned int            n_partitions,
                          const SparsityPattern        &cell_connection_graph,
                          Triangulation<dim, spacedim> &triangulation,
                          const SparsityTools::Partitioner partitioner =
                            SparsityTools::Partitioner::metis);

  /**
   * This function performs the same operation as the one above, except that
   * it takes into consideration a specific set of @p cell_weights, which allow the
   * partitioner to balance the graph while taking into consideration the
   * computational effort expended on each cell.
   *
   * @note If the @p cell_weights vector is empty, then no weighting is taken
   * into consideration. If not then the size of this vector must equal to the
   * number of active cells in the triangulation.
   */
  template <int dim, int spacedim>
  void
  partition_triangulation(const unsigned int               n_partitions,
                          const std::vector<unsigned int> &cell_weights,
                          const SparsityPattern        &cell_connection_graph,
                          Triangulation<dim, spacedim> &triangulation,
                          const SparsityTools::Partitioner partitioner =
                            SparsityTools::Partitioner::metis);

  /**
   * Generates a partitioning of the active cells making up the entire domain
   * using the same partitioning scheme as in the p4est library if the flag
   * @p group_siblings is set to true (default behavior of this function).
   * After calling this function, the subdomain ids of all active cells will
   * have values between zero and @p n_partitions-1. You can access the
   * subdomain id of a cell by using <tt>cell-@>subdomain_id()</tt>.
   *
   * @note If the flag @p group_siblings is set to false, children of a
   *       cell might be placed on different processors even though they are all
   *       active, which is an assumption made by p4est. By relaxing this, we
   *       can create partitions owning a single cell (also for refined
   *       meshes).
   */
  template <int dim, int spacedim>
  void
  partition_triangulation_zorder(const unsigned int            n_partitions,
                                 Triangulation<dim, spacedim> &triangulation,
                                 const bool group_siblings = true);

  /**
   * Partitions the cells of a multigrid hierarchy by assigning level subdomain
   * ids using the "youngest child" rule, that is, each cell in the hierarchy is
   * owned by the processor who owns its left most child in the forest, and
   * active cells have the same subdomain id and level subdomain id. You can
   * access the level subdomain id of a cell by using
   * <tt>cell-@>level_subdomain_id()</tt>.
   *
   * Note: This function assumes that the active cells have already been
   * partitioned.
   */
  template <int dim, int spacedim>
  void
  partition_multigrid_levels(Triangulation<dim, spacedim> &triangulation);

  /**
   * This function allows to ask for the owning subdomain of cells identified by
   * CellId objects that do not have to exist on the current process.
   *
   * @note This function has not been implemented yet for
   *   parallel::fullydistributed::Triangulation.
   */
  template <int dim, int spacedim>
  std::vector<types::subdomain_id>
  get_subdomain_association(const Triangulation<dim, spacedim> &triangulation,
                            const std::vector<CellId>          &cell_ids);

  /**
   * For each active cell, return in the output array to which subdomain (as
   * given by the <tt>cell->subdomain_id()</tt> function) it belongs. The
   * output array is supposed to have the right size already when calling this
   * function.
   *
   * This function returns the association of each cell with one subdomain. If
   * you are looking for the association of each @em DoF with a subdomain, use
   * the <tt>DoFTools::get_subdomain_association</tt> function.
   */
  template <int dim, int spacedim>
  void
  get_subdomain_association(const Triangulation<dim, spacedim> &triangulation,
                            std::vector<types::subdomain_id>   &subdomain);

  /**
   * Count how many cells are uniquely associated with the given @p subdomain
   * index.
   *
   * This function may return zero if there are no cells with the given @p
   * subdomain index. This can happen, for example, if you try to partition a
   * coarse mesh into more partitions (one for each processor) than there are
   * cells in the mesh.
   *
   * This function returns the number of cells associated with one subdomain.
   * If you are looking for the association of @em DoFs with this subdomain,
   * use the <tt>DoFTools::count_dofs_with_subdomain_association</tt>
   * function.
   */
  template <int dim, int spacedim>
  unsigned int
  count_cells_with_subdomain_association(
    const Triangulation<dim, spacedim> &triangulation,
    const types::subdomain_id           subdomain);

  /**
   * For a triangulation, return a mask that represents which of its vertices
   * are "owned" by the current process in the same way as we talk about
   * locally owned cells or degrees of freedom (see
   * @ref GlossLocallyOwnedCell
   * and
   * @ref GlossLocallyOwnedDof).
   * For the purpose of this function, we define a locally owned vertex as
   * follows: a vertex is owned by that processor with the smallest subdomain
   * id (which equals the MPI rank of that processor) among all owners of
   * cells adjacent to this vertex. In other words, vertices that are in the
   * interior of a partition of the triangulation are owned by the owner of
   * this partition; for vertices that lie on the boundary between two or more
   * partitions, the owner is the processor with the least subdomain_id among
   * all adjacent subdomains.
   *
   * For sequential triangulations (as opposed to, for example,
   * parallel::distributed::Triangulation), every user vertex is of course
   * owned by the current processor, i.e., the function returns
   * Triangulation::get_used_vertices(). For parallel triangulations, the
   * returned mask is a subset of what Triangulation::get_used_vertices()
   * returns.
   *
   * @param triangulation The triangulation of which the function evaluates
   * which vertices are locally owned.
   * @return The subset of vertices, as described above. The length of the
   * returned array equals Triangulation.n_vertices() and may, consequently,
   * be larger than Triangulation::n_used_vertices().
   */
  template <int dim, int spacedim>
  std::vector<bool>
  get_locally_owned_vertices(const Triangulation<dim, spacedim> &triangulation);

  /** @} */
  /**
   * @name Dealing with distorted cells
   */
  /** @{ */

  /**
   * Given a triangulation and a list of cells whose children have become
   * distorted as a result of mesh refinement, try to fix these cells up by
   * moving the center node around.
   *
   * The function returns a list of cells with distorted children that
   * couldn't be fixed up for whatever reason. The returned list is therefore
   * a subset of the input argument.
   *
   * For a definition of the concept of distorted cells, see the
   * @ref GlossDistorted "glossary entry".
   * The first argument passed to the current function is typically the
   * exception thrown by the Triangulation::execute_coarsening_and_refinement
   * function.
   *
   * @deprecated This function predates deal.II's use of manifolds and use of
   * cell-local transfinite interpolation to place new points and is no longer
   * necessary. See Manifolds::get_default_points_and_weights() for more
   * information.
   */
  template <int dim, int spacedim>
  DEAL_II_DEPRECATED typename Triangulation<dim, spacedim>::DistortedCellList
  fix_up_distorted_child_cells(
    const typename Triangulation<dim, spacedim>::DistortedCellList
                                 &distorted_cells,
    Triangulation<dim, spacedim> &triangulation);



  /** @} */
  /**
   * @name Extracting and creating patches of cells
   *
   * These functions extract and create patches of cells surrounding a single
   * cell, and creating triangulation out of them.
   */
  /** @{ */


  /**
   * This function returns a list of all the active neighbor cells of the
   * given, active cell.  Here, a neighbor is defined as one having at least
   * part of a face in common with the given cell, but not edge (in 3d) or
   * vertex neighbors (in 2d and 3d).
   *
   * The first element of the returned list is the cell provided as argument.
   * The remaining ones are neighbors: The function loops over all faces of
   * that given cell and checks if that face is not on the boundary of the
   * domain. Then, if the neighbor cell does not have any children (i.e., it
   * is either at the same refinement level as the current cell, or coarser)
   * then this neighbor cell is added to the list of cells. Otherwise, if the
   * neighbor cell is refined and therefore has children, then this function
   * loops over all subfaces of current face adds the neighbors behind these
   * sub-faces to the list to be returned.
   *
   * @tparam MeshType A type that satisfies the requirements of the
   * @ref ConceptMeshType "MeshType concept".
   * In C++, the compiler can not determine <code>MeshType</code> from the
   * function call. You need to specify it as an explicit template argument
   * following the function name.
   * @param[in] cell An iterator pointing to a cell of the mesh.
   * @return A list of active cells that form the patch around the given cell
   *
   * @note Patches are often used in defining error estimators that require
   * the solution of a local problem on the patch surrounding each of the
   * cells of the mesh. This also requires manipulating the degrees of freedom
   * associated with the cells of a patch. To this end, there are further
   * functions working on patches in namespace DoFTools.
   *
   * @note In the context of a parallel distributed computation, it only makes
   * sense to call this function on locally owned cells. This is because the
   * neighbors of locally owned cells are either locally owned themselves, or
   * ghost cells. For both, we know that these are in fact the real cells of
   * the complete, parallel triangulation. We can also query the degrees of
   * freedom on these.
   *
   * @dealiiConceptRequires{concepts::is_triangulation_or_dof_handler<MeshType>}
   */
  template <typename MeshType>
  DEAL_II_CXX20_REQUIRES(concepts::is_triangulation_or_dof_handler<MeshType>)
  std::vector<typename MeshType::active_cell_iterator> get_patch_around_cell(
    const typename MeshType::active_cell_iterator &cell);


  /**
   * This function takes a vector of active cells (hereafter named @p
   * patch_cells) as input argument, and returns a vector of their parent
   * cells with the coarsest common level of refinement. In other words, find
   * that set of cells living at the same refinement level so that all cells
   * in the input vector are children of the cells in the set, or are in the
   * set itself.
   *
   * @tparam Container In C++, the compiler can not determine the type of
   * <code>Container</code> from the function call. You need to specify it as
   * an explicit template argument following the function name. This type has
   * to satisfy the requirements of a mesh container (see
   * @ref ConceptMeshType).
   *
   * @param[in] patch_cells A vector of active cells for which this function
   * finds the parents at the coarsest common level. This vector of cells
   * typically results from calling the function
   * GridTools::get_patch_around_cell().
   * @return A list of cells with the coarsest common level of refinement of
   * the input cells.
   */
  template <class Container>
  std::vector<typename Container::cell_iterator>
  get_cells_at_coarsest_common_level(
    const std::vector<typename Container::active_cell_iterator> &patch_cells);

  /**
   * This function constructs a Triangulation (named @p local_triangulation)
   * from a given vector of active cells. This vector (which we think of the
   * cells corresponding to a "patch") contains active cells that are part of
   * an existing global Triangulation. The goal of this function is to build a
   * local Triangulation that contains only the active cells given in @p patch
   * (and potentially a minimum number of additional cells required to form a
   * valid Triangulation). The function also returns a map that allows to
   * identify the cells in the output Triangulation and corresponding cells in
   * the input list.
   *
   * The function copies the location of vertices of cells from the cells of the
   * source triangulation to the triangulation that is built from the list of
   * patch cells.  This adds support for triangulations which have been
   * perturbed or smoothed in some manner which makes the triangulation
   * deviate from the standard deal.II refinement strategy of placing new
   * vertices at midpoints of faces or edges.
   *
   * The operation implemented by this function is frequently used in the
   * definition of error estimators that need to solve "local" problems on
   * each cell and its neighbors. A similar construction is necessary in the
   * definition of the Clement interpolation operator in which one needs to
   * solve a local problem on all cells within the support of a shape
   * function. This function then builds a complete Triangulation from a list
   * of cells that make up such a patch; one can then later attach a
   * DoFHandler to such a Triangulation.
   *
   * If the list of input cells contains only cells at the same refinement
   * level, then the output Triangulation simply consists of a Triangulation
   * containing only exactly these patch cells. On the other hand, if the
   * input cells live on different refinement levels, i.e., the Triangulation
   * of which they are part is adaptively refined, then the construction of
   * the output Triangulation is not so simple because the coarsest level of a
   * Triangulation can not contain hanging nodes. Rather, we first have to
   * find the common refinement level of all input cells, along with their
   * common parents (see GridTools::get_cells_at_coarsest_common_level()),
   * build a Triangulation from those, and then adaptively refine it so that
   * the input cells all also exist in the output Triangulation.
   *
   * A consequence of this procedure is that the output Triangulation may
   * contain more active cells than the ones that exist in the input vector.
   * On the other hand, one typically wants to solve the local problem not on
   * the entire output Triangulation, but only on those cells of it that
   * correspond to cells in the input list.  In this case, a user typically
   * wants to assign degrees of freedom only on cells that are part of the
   * "patch", and somehow ignore those excessive cells. The current function
   * supports this common requirement by setting the user flag for the cells
   * in the output Triangulation that match with cells in the input list.
   * Cells which are not part of the original patch will not have their @p
   * user_flag set; we can then avoid assigning degrees of freedom using the
   * FE_Nothing<dim> element.
   *
   * @tparam Container In C++, the compiler can not determine the type of
   * <code>Container</code> from the function call. You need to specify it as
   * an explicit template argument following the function name. This type that
   * satisfies the requirements of a mesh container (see
   * @ref ConceptMeshType).
   *
   * @param[in] patch A vector of active cells from a common triangulation.
   * These cells may or may not all be at the same refinement level.
   * @param[out] local_triangulation A triangulation whose active cells
   * correspond to the given vector of active cells in @p patch.
   * @param[out] patch_to_global_tria_map A map between the local
   * triangulation which is built as explained above, and the cell iterators
   * in the input list.
   */
  template <class Container>
  void
  build_triangulation_from_patch(
    const std::vector<typename Container::active_cell_iterator> &patch,
    Triangulation<Container::dimension, Container::space_dimension>
      &local_triangulation,
    std::map<
      typename Triangulation<Container::dimension,
                             Container::space_dimension>::active_cell_iterator,
      typename Container::active_cell_iterator> &patch_to_global_tria_map);

  /**
   * This function runs through the degrees of freedom defined by the
   * DoFHandler and for each dof constructs a vector of
   * active_cell_iterators representing the cells of support of the associated
   * basis element at that degree of freedom. This function was originally
   * designed for the implementation of local projections, for instance the
   * Clement interpolant, in conjunction with other local patch functions like
   * GridTools::build_triangulation_from_patch.
   *
   * DoFHandler's built on top of Triangulation or
   * parallel::distributed::Triangulation are supported and handled
   * appropriately.
   *
   * The result is the patch of cells representing the support of the basis
   * element associated to the degree of freedom.  For instance using an FE_Q
   * finite element, we obtain the standard patch of cells touching the degree
   * of freedom and then add other cells that take care of possible hanging node
   * constraints.  Using a FE_DGQ finite element, the degrees of freedom are
   * logically considered to be "interior" to the cells so the patch would
   * consist exclusively of the single cell on which the degree of freedom is
   * located.
   *
   * @param[in] dof_handler The DoFHandler which could be built on a
   * Triangulation or a parallel::distributed::Triangulation with a finite
   * element that has degrees of freedom that are logically associated to a
   * vertex, line, quad, or hex.
   * @return A map from the global_dof_index of
   * degrees of freedom on locally relevant cells to vectors containing
   * DoFHandler::active_cell_iterators of cells in the support of the basis
   * function at that degree of freedom.
   */
  template <int dim, int spacedim>
  std::map<
    types::global_dof_index,
    std::vector<typename DoFHandler<dim, spacedim>::active_cell_iterator>>
  get_dof_to_support_patch_map(DoFHandler<dim, spacedim> &dof_handler);


  /** @} */

  /**
   * @name Dealing with periodic domains
   */
  /** @{ */

  /**
   * Data type that provides all information necessary to create periodicity
   * constraints and a periodic p4est forest with respect to two 'periodic'
   * cell faces.
   *
   * @ingroup reordering
   */
  template <typename CellIterator>
  struct PeriodicFacePair
  {
    /**
     * The cells associated with the two 'periodic' faces.
     */
    CellIterator cell[2];

    /**
     * The local face indices (with respect to the specified cells) of the two
     * 'periodic' faces.
     */
    unsigned int face_idx[2];

    /**
     * The relative orientation of the first face with respect to the second
     * face as described in orthogonal_equality() and
     * DoFTools::make_periodicity_constraints().
     */
    unsigned char orientation;

    /**
     * A @p dim $\times$ @p dim rotation matrix that describes how vector
     * valued DoFs of the first face should be modified prior to constraining
     * to the DoFs of the second face.
     *
     * The rotation matrix is used in DoFTools::make_periodicity_constraints()
     * by applying the rotation to all vector valued blocks listed in the
     * parameter @p first_vector_components of the finite element space. For
     * more details see DoFTools::make_periodicity_constraints() and the
     * glossary
     * @ref GlossPeriodicConstraints "glossary entry on periodic conditions".
     */
    FullMatrix<double> matrix;

    /**
     * Return an estimate, in bytes, for the memory consumption of the object.
     */
    std::size_t
    memory_consumption() const;
  };


  /**
   * An orthogonal equality test for faces.
   *
   * @p face1 and @p face2 are considered equal, if a one to one matching
   * between its vertices can be achieved via an orthogonal equality relation.
   * If no such relation exists then the returned std::optional object is empty
   * (i.e., has_value() will return `false`).
   *
   * Here, two vertices <tt>v_1</tt> and <tt>v_2</tt> are considered equal, if
   * $M\cdot v_1 + offset - v_2$ is parallel to the unit vector in unit
   * direction @p direction. If the parameter @p matrix is a reference to a
   * spacedim x spacedim matrix, $M$ is set to @p matrix, otherwise $M$ is the
   * identity matrix.
   *
   * If the matching was successful, the _relative_ orientation of @p face1 with
   * respect to @p face2 is returned a std::optional<unsigned char>, in which
   * the stored value is the same orientation bit format used elsewhere in the
   * library. More information on that topic can be found in the
   * @ref GlossFaceOrientation "glossary"
   * article.
   */
  template <typename FaceIterator>
  std::optional<unsigned char>
  orthogonal_equality(
    const FaceIterator                                           &face1,
    const FaceIterator                                           &face2,
    const unsigned int                                            direction,
    const Tensor<1, FaceIterator::AccessorType::space_dimension> &offset =
      Tensor<1, FaceIterator::AccessorType::space_dimension>(),
    const FullMatrix<double> &matrix = FullMatrix<double>());

  /**
   * This function will collect periodic face pairs on the coarsest mesh level
   * of the given @p mesh (a Triangulation or DoFHandler) and add them to the
   * vector @p matched_pairs leaving the original contents intact.
   *
   * Define a 'first' boundary as all boundary faces having boundary_id @p
   * b_id1 and a 'second' boundary consisting of all faces belonging to @p
   * b_id2.
   *
   * This function tries to match all faces belonging to the first boundary
   * with faces belonging to the second boundary with the help of
   * orthogonal_equality().
   *
   * The unsigned char that is returned inside of PeriodicFacePair encodes the
   * _relative_ orientation of the first face with respect to the second face,
   * see the documentation of orthogonal_equality() for further details.
   *
   * The @p direction refers to the space direction in which periodicity is
   * enforced. When matching periodic faces this vector component is ignored.
   *
   * The @p offset is a vector tangential to the faces that is added to the
   * location of vertices of the 'first' boundary when attempting to match
   * them to the corresponding vertices of the 'second' boundary. This can be
   * used to implement conditions such as $u(0,y)=u(1,y+1)$.
   *
   * Optionally, a $dim\times dim$ rotation @p matrix can be specified that
   * describes how vector valued DoFs of the first face should be modified
   * prior to constraining to the DoFs of the second face. The @p matrix is
   * used in two places. First, @p matrix will be supplied to
   * orthogonal_equality() and used for matching faces: Two vertices $v_1$ and
   * $v_2$ match if $\text{matrix}\cdot v_1 + \text{offset} - v_2$ is parallel
   * to the unit vector in unit direction @p direction. (For more details see
   * DoFTools::make_periodicity_constraints(), the glossary
   * @ref GlossPeriodicConstraints "glossary entry on periodic conditions"
   * and step-45). Second, @p matrix will be stored in the PeriodicFacePair
   * collection @p matched_pairs for further use.
   *
   * @tparam MeshType A type that satisfies the requirements of the
   * @ref ConceptMeshType "MeshType concept".
   *
   * @note The created std::vector can be used in
   * DoFTools::make_periodicity_constraints() and in
   * parallel::distributed::Triangulation::add_periodicity() to enforce
   * periodicity algebraically.
   *
   * @note Because elements will be added to @p matched_pairs (and existing
   * entries will be preserved), it is possible to call this function several
   * times with different boundary ids to generate a vector with all periodic
   * pairs.
   *
   * @note Since the periodic face pairs are found on the coarsest mesh level,
   * it is necessary to ensure that the coarsest level faces have the correct
   * boundary indicators set. In general, this means that one must first set
   * all boundary indicators on the coarse grid before performing any global
   * or local grid refinement.
   *
   * @dealiiConceptRequires{concepts::is_triangulation_or_dof_handler<MeshType>}
   */
  template <typename MeshType>
  DEAL_II_CXX20_REQUIRES(concepts::is_triangulation_or_dof_handler<MeshType>)
  void collect_periodic_faces(
    const MeshType          &mesh,
    const types::boundary_id b_id1,
    const types::boundary_id b_id2,
    const unsigned int       direction,
    std::vector<PeriodicFacePair<typename MeshType::cell_iterator>>
                                               &matched_pairs,
    const Tensor<1, MeshType::space_dimension> &offset =
      dealii::Tensor<1, MeshType::space_dimension>(),
    const FullMatrix<double> &matrix = FullMatrix<double>());


  /**
   * This compatibility version of collect_periodic_faces() only works on
   * grids with cells in
   * @ref GlossFaceOrientation "standard orientation".
   *
   * Instead of defining a 'first' and 'second' boundary with the help of two
   * boundary_ids this function defines a 'left' boundary as all faces with
   * local face index <code>2*direction</code> and boundary indicator @p b_id
   * and, similarly, a 'right' boundary consisting of all face with local face
   * index <code>2*direction+1</code> and boundary indicator @p b_id. Faces
   * with coordinates only differing in the @p direction component are
   * identified.
   *
   * This function will collect periodic face pairs on the coarsest mesh level
   * and add them to @p matched_pairs leaving the original contents intact.
   *
   * See above function for further details.
   *
   * @note This version of collect_periodic_faces() will not work on
   * meshes with cells not in
   * @ref GlossFaceOrientation "standard orientation".
   *
   * @dealiiConceptRequires{concepts::is_triangulation_or_dof_handler<MeshType>}
   */
  template <typename MeshType>
  DEAL_II_CXX20_REQUIRES(concepts::is_triangulation_or_dof_handler<MeshType>)
  void collect_periodic_faces(
    const MeshType          &mesh,
    const types::boundary_id b_id,
    const unsigned int       direction,
    std::vector<PeriodicFacePair<typename MeshType::cell_iterator>>
                                                       &matched_pairs,
    const dealii::Tensor<1, MeshType::space_dimension> &offset =
      dealii::Tensor<1, MeshType::space_dimension>(),
    const FullMatrix<double> &matrix = FullMatrix<double>());

  /** @} */
  /**
   * @name Dealing with boundary and manifold ids
   */
  /** @{ */

  /**
   * Copy boundary ids to manifold ids on faces and edges at the boundary. The
   * default manifold_id for new Triangulation objects is
   * numbers::flat_manifold_id. This function copies the boundary_ids of
   * the boundary faces and edges to the manifold_ids of the same faces and
   * edges, allowing the user to change the boundary_ids and use them for
   * boundary conditions regardless of the geometry, which will use
   * manifold_ids to create new points. Only active cells will be iterated
   * over. This is a function you'd typically call when there is only one
   * active level on your Triangulation. Mesh refinement will then inherit
   * these indicators to child cells, faces, and edges.
   *
   * The optional parameter @p reset_boundary_ids, indicates whether this
   * function should reset the boundary_ids of boundary faces and edges to its
   * default value 0 after copying its value to the manifold_id. By default,
   * boundary_ids are left untouched.
   *
   * @ingroup manifold
   * @relatesalso boundary
   */
  template <int dim, int spacedim>
  void
  copy_boundary_to_manifold_id(Triangulation<dim, spacedim> &tria,
                               const bool reset_boundary_ids = false);

  /**
   * Map the given boundary ids to the given manifold ids on faces and
   * edges at the boundary.
   *
   * This function copies the boundary ids of the boundary faces and
   * edges that are present in the parameter @p src_boundary_ids to
   * the corresponding manifold id in @p dst_manifold_ids, of the same
   * faces and edges.
   *
   * If the optional parameter @p reset_boundary_ids is non empty,
   * each boundary id in @p src_boundary_ids, is replaced with the
   * corresponding boundary id in @p reset_boundary_ids.
   *
   * An exception is thrown if the size of the input vectors do not
   * match. If a boundary id indicated in @p src_boundary_ids is not
   * present in the triangulation, it is simply ignored during the
   * process.
   *
   * @ingroup manifold
   * @relatesalso boundary
   */
  template <int dim, int spacedim>
  void
  map_boundary_to_manifold_ids(
    const std::vector<types::boundary_id> &src_boundary_ids,
    const std::vector<types::manifold_id> &dst_manifold_ids,
    Triangulation<dim, spacedim>          &tria,
    const std::vector<types::boundary_id> &reset_boundary_ids = {});

  /**
   * Copy material ids to manifold ids. The default manifold_id for new
   * Triangulation objects is numbers::flat_manifold_id. When refinements
   * occurs, the Triangulation asks where to locate new points to the
   * underlying manifold.
   *
   * When reading a Triangulation from a supported input format, typical
   * information that can be stored in a file are boundary conditions for
   * boundary faces (which we store in the boundary_id of the faces), material
   * types for cells (which we store in the material_id of the cells) and in
   * some cases subdomain ids for cells (which we store in the subdomain_id of
   * the cell).
   *
   * If you read one of these grids into a Triangulation, you might still want
   * to use the material_id specified in the input file as a manifold_id
   * description. In this case you can associate a Manifold object to internal
   * cells, and this object will be used by the Triangulation to query
   * Manifold objects for new points. This function iterates over active cells
   * and copies the material_ids to the manifold_ids.
   *
   * The optional parameter @p compute_face_ids, indicates whether this
   * function should also set the manifold_ids of the faces (both for internal
   * faces and for faces on the boundary). If set to true, then each face will
   * get a manifold_id equal to the minimum of the surrounding manifold_ids,
   * ensuring that a unique manifold id is selected for each face of the
   * Triangulation. By default, face manifold_ids are not computed.
   *
   * @ingroup manifold
   */
  template <int dim, int spacedim>
  void
  copy_material_to_manifold_id(Triangulation<dim, spacedim> &tria,
                               const bool compute_face_ids = false);

  /**
   * Propagate manifold indicators associated with the cells of the
   * Triangulation @p tria to their co-dimension one and two objects.
   *
   * This function sets the @p manifold_id of faces and edges (both on the
   * interior and on the boundary) to the value returned by the
   * @p disambiguation_function method, called with the set of
   * manifold indicators of the cells that share the same face or edge.
   *
   * By default, the @p disambiguation_function returns
   * numbers::flat_manifold_id when the set has size greater than one (i.e.,
   * when it is not possible to decide what manifold indicator a face or edge
   * should have according to the manifold indicators of the adjacent cells)
   * and it returns the manifold indicator contained in the set when it has
   * dimension one (i.e., when all adjacent cells and faces have the same
   * manifold indicator).
   *
   * The parameter @p overwrite_only_flat_manifold_ids allows you to specify
   * what to do when a face or an edge already has a manifold indicator
   * different from numbers::flat_manifold_id. If the flag is @p true, the edge
   * or face will maintain its original manifold indicator.
   * If it is @p false, then also the manifold indicator of these faces and edges
   * is set according to the return value of the @p disambiguation_function.
   */
  template <int dim, int spacedim>
  void
  assign_co_dimensional_manifold_indicators(
    Triangulation<dim, spacedim>             &tria,
    const std::function<types::manifold_id(
      const std::set<types::manifold_id> &)> &disambiguation_function =
      [](const std::set<types::manifold_id> &manifold_ids) {
        if (manifold_ids.size() == 1)
          return *manifold_ids.begin();
        else
          return numbers::flat_manifold_id;
      },
    bool overwrite_only_flat_manifold_ids = true);
  /** @} */

  /**
   * Exchange arbitrary data of type @p DataType provided by the function
   * objects from locally owned cells to ghost cells on other processors.
   *
   * After this call, you typically will have received data from @p unpack on
   * every ghost cell as it was given by @p pack on the owning processor.
   * Whether you do or do not receive information to @p unpack on a given
   * ghost cell depends on whether the @p pack function decided that
   * something needs to be sent. It does so using the std::optional
   * mechanism: if the std::optional return object of the @p pack
   * function is empty, then this implies that no data has to be sent for
   * the locally owned cell it was called on. In that case, @p unpack will
   * also not be called on the ghost cell that corresponds to it on the
   * receiving side. On the other hand, if the std::optional object is
   * not empty, then the data stored within it will be sent to the received
   * and the @p unpack function called with it.
   *
   * @tparam DataType The type of the data to be communicated. It is assumed
   *   to be serializable by boost::serialization. In many cases, this
   *   data type can not be deduced by the compiler, e.g., if you provide
   *   lambda functions for the second and third argument
   *   to this function. In this case, you have to explicitly specify
   *   the @p DataType as a template argument to the function call.
   * @tparam MeshType The type of @p mesh.
   *
   * @param mesh A variable of a type that satisfies the requirements of the
   * @ref ConceptMeshType "MeshType concept".
   * @param pack The function that will be called on each locally owned cell
   *   that is a ghost cell somewhere else. As mentioned above, the function
   *   may return a regular data object of type @p DataType to indicate
   *   that data should be sent, or an empty
   *   <code>std::optional@<DataType@></code> to indicate that nothing has
   *   to be sent for this cell.
   * @param unpack The function that will be called for each ghost cell
   *   for which data was sent, i.e., for which the @p pack function
   *   on the sending side returned a non-empty std::optional object.
   *   The @p unpack function is then called with the data sent by the
   *   processor that owns that cell.
   * @param cell_filter Only cells are communicated where this filter function returns
   *   the value `true`. In the default case, the function returns true on all
   * cells and thus, all relevant cells are communicated.
   *
   * <h4> An example </h4>
   *
   * Here is an example that shows how this function is to be used
   * in a concrete context. It is taken from the code that makes
   * sure that the @p active_fe_index (a single unsigned integer) is
   * transported from locally owned cells where one can set it in
   * DoFHandler objects with hp-capabilities, to the corresponding ghost cells
   * on other processors to ensure that one can query the right value also on
   * those processors:
   * @code
   * using active_cell_iterator =
   *   typename DoFHandler<dim,spacedim>::active_cell_iterator;
   * auto pack = [] (const active_cell_iterator &cell) -> unsigned int
   *             {
   *               return cell->active_fe_index();
   *             };
   *
   * auto unpack = [] (const active_cell_iterator &cell,
   *                   const unsigned int active_fe_index) -> void
   *               {
   *                 cell->set_active_fe_index(active_fe_index);
   *               };
   *
   * GridTools::exchange_cell_data_to_ghosts<
   *   unsigned int, DoFHandler<dim,spacedim>> (dof_handler,
   *                                                    pack,
   *                                                    unpack);
   * @endcode
   *
   * You will notice that the @p pack lambda function returns an `unsigned int`,
   * not a `std::optional<unsigned int>`. The former converts
   * automatically to the latter, implying that data will always be transported
   * to the other processor.
   *
   * (In reality, the @p unpack function needs to be a bit more
   * complicated because it is not allowed to call
   * DoFAccessor::set_active_fe_index() on ghost cells. Rather, the
   * @p unpack function directly accesses internal data structures. But
   * you get the idea -- the code could, just as well, have exchanged
   * material ids, user indices, boundary indicators, or any kind of other
   * data with similar calls as the ones above.)
   *
   * @dealiiConceptRequires{concepts::is_triangulation_or_dof_handler<MeshType>}
   */
  template <typename DataType, typename MeshType>
  DEAL_II_CXX20_REQUIRES(concepts::is_triangulation_or_dof_handler<MeshType>)
  void exchange_cell_data_to_ghosts(
    const MeshType                                      &mesh,
    const std::function<std::optional<DataType>(
      const typename MeshType::active_cell_iterator &)> &pack,
    const std::function<void(const typename MeshType::active_cell_iterator &,
                             const DataType &)>         &unpack,
    const std::function<bool(const typename MeshType::active_cell_iterator &)>
      &cell_filter =
        always_return<typename MeshType::active_cell_iterator, bool>{true});

  /**
   * Exchange arbitrary data of type @p DataType provided by the function
   * objects from locally owned level cells to ghost level cells on other
   * processes.
   *
   * In addition to the parameters of exchange_cell_data_to_ghosts(), this
   * function allows to provide a @p cell_filter function, which can be used to only
   * communicate marked cells. In the default case, all relevant cells are
   * communicated.
   *
   * @dealiiConceptRequires{concepts::is_triangulation_or_dof_handler<MeshType>}
   */
  template <typename DataType, typename MeshType>
  DEAL_II_CXX20_REQUIRES(concepts::is_triangulation_or_dof_handler<MeshType>)
  void exchange_cell_data_to_level_ghosts(
    const MeshType                                     &mesh,
    const std::function<std::optional<DataType>(
      const typename MeshType::level_cell_iterator &)> &pack,
    const std::function<void(const typename MeshType::level_cell_iterator &,
                             const DataType &)>        &unpack,
    const std::function<bool(const typename MeshType::level_cell_iterator &)> &
      cell_filter = always_return<typename MeshType::level_cell_iterator, bool>{
        true});

  /* Exchange with all processors of the MPI communicator @p mpi_communicator the vector of bounding
   * boxes @p local_bboxes.
   *
   * This function is meant to exchange bounding boxes describing the locally
   * owned cells in a distributed triangulation obtained with the function
   * GridTools::compute_mesh_predicate_bounding_box .
   *
   * The output vector's size is the number of processes of the MPI
   * communicator:
   * its i-th entry contains the vector @p local_bboxes of the i-th process.
   */
  template <int spacedim>
  std::vector<std::vector<BoundingBox<spacedim>>>
  exchange_local_bounding_boxes(
    const std::vector<BoundingBox<spacedim>> &local_bboxes,
    const MPI_Comm                            mpi_communicator);

  /**
   * In this @ref GlossCollectiveOperation "collective operation" each process provides a vector
   * of bounding boxes and a communicator.
   * All these vectors are gathered on each of the processes,
   * organized in a search tree, and then returned.
   *
   * The idea is that the vector of bounding boxes describes a
   * relevant property of the computations on each process
   * individually, which could also be of use to other processes. An
   * example would be if the input vector of bounding boxes
   * corresponded to a covering of the locally owned partition of a
   * mesh (see
   * @ref GlossLocallyOwnedCell)
   * of a
   * parallel::distributed::Triangulation object. While these may
   * overlap the bounding boxes of other processes, finding which
   * process owns the cell that encloses a given point is vastly
   * easier if the process trying to figure this out has a (relatively
   * small) list of processes whose bounding boxes contain that point.
   *
   * The returned search tree object is an r-tree with packing
   * algorithm, which is provided by boost library. See
   * https://www.boost.org/doc/libs/1_67_0/libs/geometry/doc/html/geometry/spatial_indexes/introduction.html
   * for more information.
   *
   * In the returned tree, each node contains a pair of elements:
   * the first being a bounding box,
   * the second being the rank of the process for which at least some
   * of the locally owned cells overlap with the bounding box.
   *
   * @note This function is a @ref GlossCollectiveOperation "collective operation".
   */
  template <int spacedim>
  RTree<std::pair<BoundingBox<spacedim>, unsigned int>>
  build_global_description_tree(
    const std::vector<BoundingBox<spacedim>> &local_description,
    const MPI_Comm                            mpi_communicator);

  /**
   * Collect for a given triangulation all locally relevant vertices that
   * coincide due to periodicity.
   *
   * Coinciding vertices are put into a group, e.g.: [1, 25, 51], which is
   * labeled by an arbitrary element from it, e.g.: "1". All coinciding vertices
   * store the label to its group, so that they can quickly access all the
   * coinciding vertices in that group: e.g.: 51 ->  "1" -> [1, 25, 51]
   *
   * @param[in] tria Triangulation.
   * @param[out] coinciding_vertex_groups A map of equivalence classes (of
   *             coinciding vertices) labeled by an arbitrary element from them.
   *             Vertices not coinciding are ignored.
   * @param[out] vertex_to_coinciding_vertex_group Map of a vertex to the label
   *             of a group of coinciding vertices. Vertices not contained in
   *             this vector are not coinciding with any other vertex.
   */
  template <int dim, int spacedim>
  void
  collect_coinciding_vertices(
    const Triangulation<dim, spacedim>                &tria,
    std::map<unsigned int, std::vector<unsigned int>> &coinciding_vertex_groups,
    std::map<unsigned int, unsigned int> &vertex_to_coinciding_vertex_group);

  /**
   * Return a map that, for each vertex of the given triangulation,
   * provides a set of all the process subdomain ids whose subdomains
   * are adjacent to that vertex. The set excludes the subdomain id
   * of the current process. As a consequence, for a given vertex,
   * the returned set consists of exactly those subdomain ids that
   * correspond to the
   * @ref GlossGhostCell "ghost cells"
   * adjacent to that vertex, assuming there are any such ghost cells.
   *
   * For vertices that are not adjacent to a ghost cell, the map contains
   * no entries, and this should be interpreted in the same way as if
   * the map contained an entry for a given vertex index, but that the
   * `std::set` associated with that map entry is simply empty.
   * For non-parallel triangulations, the map is consequently empty
   * since no vertex has adjacent ghost cells.
   *
   * @param[in] tria Triangulation.
   */
  template <int dim, int spacedim>
  std::map<unsigned int, std::set<dealii::types::subdomain_id>>
  compute_vertices_with_ghost_neighbors(
    const Triangulation<dim, spacedim> &tria);

  /**
   * An implementation of the marching-square (2d) and marching-cube algorithm
   * for creating data structures (vectors of Point and CellData) to
   * create a linear/bilinear surface mesh on the iso line/contour of a
   * scalar field.
   *
   * To improve the approximation of the iso line/contour and the resulting
   * linear surface mesh, one increases the number of subdivision so that the
   * algorithm is not run on a cell but on subcells with vertex values having
   * been interpolated from the cell values.
   *
   * @note The resulting mesh will contain lines in 2d and triangles in 3d.
   *
   * @note The resulting mesh will not be of high quality, since it might
   *   contain cells with very small diameters if the mesh is cut close to a
   *   vertex.
   *
   * @note Iso lines/contours as a saddle point within a subcell is not
   *       detected by the implemented algorithm.
   */
  template <int dim, typename VectorType>
  class MarchingCubeAlgorithm
  {
  public:
    /**
     * Value type of vector.
     */
    using value_type = typename VectorType::value_type;

    /**
     * Constructor.
     */
    MarchingCubeAlgorithm(const Mapping<dim, dim>       &mapping,
                          const FiniteElement<dim, dim> &fe,
                          const unsigned int             n_subdivisions = 1,
                          const double                   tolerance = 1e-10);

    /**
     * Process all locally-owned cells and fill @p vertices and @p cells for all
     * cells that are cut.
     *
     * @note This function is only implemented for dim>1. Use
     * process(background_dof_handler, ls_vector, iso_level, vertices) for
     * dim==1.
     *
     * @note Duplicate vertices are not deleted.
     */
    void
    process(const DoFHandler<dim>   &background_dof_handler,
            const VectorType        &ls_vector,
            const double             iso_level,
            std::vector<Point<dim>> &vertices,
            std::vector<CellData<dim == 1 ? 1 : dim - 1>> &cells) const;

    /**
     * Process all locally-owned cells and fill @p vertices for all cells that
     * are cut.
     */
    void
    process(const DoFHandler<dim>   &background_dof_handler,
            const VectorType        &ls_vector,
            const double             iso_level,
            std::vector<Point<dim>> &vertices) const;

    /**
     * Process the provided cell and fill @p vertices and @p cells for all cells
     * that are cut.
     *
     * @note The resulting vectors are empty if the cell is not cut.
     *
     * @note This function is only implemented for dim>1. Use
     * process_cell(cell, ls_vector, iso_level, vertices) for dim==1.
     */
    void
    process_cell(const typename DoFHandler<dim>::active_cell_iterator &cell,
                 const VectorType                              &ls_vector,
                 const double                                   iso_level,
                 std::vector<Point<dim>>                       &vertices,
                 std::vector<CellData<dim == 1 ? 1 : dim - 1>> &cells) const;
    /**
     * Process the provided cell and fill @p vertices for all cells that are cut.
     *
     * @note The resulting vector is empty if the cell is not cut.
     */
    void
    process_cell(const typename DoFHandler<dim>::active_cell_iterator &cell,
                 const VectorType        &ls_vector,
                 const double             iso_level,
                 std::vector<Point<dim>> &vertices) const;

  private:
    /**
     * Internal function to create a quadrature rule with n_subdivisions+1
     * equally-positioned quadrature points.
     */
    static Quadrature<dim>
    create_quadrature_rule(const unsigned int n_subdivisions);

    /**
     * Process a cell.
     */
    void
    process_cell(std::vector<value_type>                       &ls_values,
                 const std::vector<Point<dim>>                 &points,
                 const double                                   iso_level,
                 std::vector<Point<dim>>                       &vertices,
                 std::vector<CellData<dim == 1 ? 1 : dim - 1>> &cells,
                 const bool write_back_cell_data = true) const;

    /**
     * Dummy function for 1d processing a sub-cell.
     */
    void
    process_sub_cell(const std::vector<value_type> &,
                     const std::vector<Point<1>> &,
                     const std::vector<unsigned int> &,
                     const double,
                     std::vector<Point<1>> &,
                     std::vector<CellData<1>> &,
                     const bool) const
    {
      AssertThrow(false, ExcNotImplemented());
    }

    /**
     * Process a sub-cell (2d).
     *
     * @note Subcells with saddle points are ignored. Please increase the number
     *   of subdivisions in this case.
     */
    void
    process_sub_cell(const std::vector<value_type>   &ls_values,
                     const std::vector<Point<2>>     &points,
                     const std::vector<unsigned int> &mask,
                     const double                     iso_level,
                     std::vector<Point<2>>           &vertices,
                     std::vector<CellData<1>>        &cells,
                     const bool write_back_cell_data) const;

    /**
     * Process a sub-cell (3d).
     */
    void
    process_sub_cell(const std::vector<value_type>   &ls_values,
                     const std::vector<Point<3>>     &points,
                     const std::vector<unsigned int> &mask,
                     const double                     iso_level,
                     std::vector<Point<3>>           &vertices,
                     std::vector<CellData<2>>        &cells,
                     const bool write_back_cell_data) const;

    /**
     * Number of subdivisions each cell is subdivided into in each direction to
     * improve the approximation.
     */
    const unsigned int n_subdivisions;

    /**
     * Absolute tolerance specifying the minimum distance between a vertex and
     * the cut point so that a line is considered cut.
     */
    const double tolerance;

    /**
     * FEValues used internally and set up with a quadrature rule with the
     * correct number of subdivisions.
     */
    mutable FEValues<dim> fe_values;
  };



  /**
   * @name Exceptions
   */
  /** @{ */

  /**
   * Exception
   */
  DeclException1(ExcInvalidNumberOfPartitions,
                 int,
                 << "The number of partitions you gave is " << arg1
                 << ", but must be greater than zero.");
  /**
   * Exception
   */
  DeclException1(ExcNonExistentSubdomain,
                 int,
                 << "The subdomain id " << arg1
                 << " has no cells associated with it.");
  /**
   * Exception
   */
  DeclException0(ExcTriangulationHasBeenRefined);

  /**
   * Exception
   */
  DeclException1(ExcScalingFactorNotPositive,
                 double,
                 << "The scaling factor must be positive, but it is " << arg1
                 << '.');

  /**
   * Exception
   */
  DeclException1(ExcVertexNotUsed,
                 unsigned int,
                 << "The given vertex with index " << arg1
                 << " is not used in the given triangulation.");

  /** @} */

} /*namespace GridTools*/

/* ----------------- Template function --------------- */

#ifndef DOXYGEN

namespace GridTools
{
  template <int dim>
  double
  cell_measure(
    const std::vector<Point<dim>> &all_vertices,
    const unsigned int (&indices)[GeometryInfo<dim>::vertices_per_cell])
  {
    // We forward call to the ArrayView version:
    const ArrayView<const unsigned int> view(
      indices, GeometryInfo<dim>::vertices_per_cell);
    return cell_measure(all_vertices, view);
  }



  // This specialization is defined here so that the general template in the
  // source file doesn't need to have further 1d overloads for the internal
  // functions it calls.
  template <>
  inline Triangulation<1, 1>::DistortedCellList
  fix_up_distorted_child_cells(const Triangulation<1, 1>::DistortedCellList &,
                               Triangulation<1, 1> &)
  {
    return {};
  }



  template <int dim, typename Transformation, int spacedim>
  DEAL_II_CXX20_REQUIRES(
    (std::invocable<Transformation, Point<spacedim>> &&
     std::assignable_from<
       Point<spacedim> &,
       std::invoke_result_t<Transformation, Point<spacedim>>>))
  void transform(const Transformation         &transformation,
                 Triangulation<dim, spacedim> &triangulation)
  {
    std::vector<bool> treated_vertices(triangulation.n_vertices(), false);

    // loop over all active cells, and
    // transform those vertices that
    // have not yet been touched. note
    // that we get to all vertices in
    // the triangulation by only
    // visiting the active cells.
    typename Triangulation<dim, spacedim>::active_cell_iterator
      cell = triangulation.begin_active(),
      endc = triangulation.end();
    for (; cell != endc; ++cell)
      for (const unsigned int v : cell->vertex_indices())
        if (treated_vertices[cell->vertex_index(v)] == false)
          {
            // transform this vertex
            cell->vertex(v) = transformation(cell->vertex(v));
            // and mark it as treated
            treated_vertices[cell->vertex_index(v)] = true;
          };


    // now fix any vertices on hanging nodes so that we don't create any holes
    if (dim == 2)
      {
        typename Triangulation<dim, spacedim>::active_cell_iterator
          cell = triangulation.begin_active(),
          endc = triangulation.end();
        for (; cell != endc; ++cell)
          for (const unsigned int face : cell->face_indices())
            if (cell->face(face)->has_children() &&
                !cell->face(face)->at_boundary())
              {
                Assert(cell->reference_cell() ==
                         ReferenceCells::get_hypercube<dim>(),
                       ExcNotImplemented());

                // this line has children
                cell->face(face)->child(0)->vertex(1) =
                  (cell->face(face)->vertex(0) + cell->face(face)->vertex(1)) /
                  2.0;
              }
      }
    else if (dim == 3)
      {
        typename Triangulation<dim, spacedim>::active_cell_iterator
          cell = triangulation.begin_active(),
          endc = triangulation.end();
        for (; cell != endc; ++cell)
          for (const unsigned int face : cell->face_indices())
            if (cell->face(face)->has_children() &&
                !cell->face(face)->at_boundary())
              {
                if (static_cast<uint8_t>(cell->face(face)->refinement_case()) ==
                    RefinementCase<dim - 1>::isotropic_refinement)
                  {
                    Assert(cell->reference_cell() ==
                             ReferenceCells::get_hypercube<dim>(),
                           ExcNotImplemented());

                    // this face has hanging nodes
                    cell->face(face)->child(0)->vertex(1) =
                      (cell->face(face)->vertex(0) +
                       cell->face(face)->vertex(1)) /
                      2.0;
                    cell->face(face)->child(0)->vertex(2) =
                      (cell->face(face)->vertex(0) +
                       cell->face(face)->vertex(2)) /
                      2.0;
                    cell->face(face)->child(1)->vertex(3) =
                      (cell->face(face)->vertex(1) +
                       cell->face(face)->vertex(3)) /
                      2.0;
                    cell->face(face)->child(2)->vertex(3) =
                      (cell->face(face)->vertex(2) +
                       cell->face(face)->vertex(3)) /
                      2.0;

                    // center of the face
                    cell->face(face)->child(0)->vertex(3) =
                      (cell->face(face)->vertex(0) +
                       cell->face(face)->vertex(1) +
                       cell->face(face)->vertex(2) +
                       cell->face(face)->vertex(3)) /
                      4.0;
                  }
                else
                  {
                    // Special case for anisotropic refinement
                    for (unsigned int line = 0;
                         line < GeometryInfo<dim - 1>::faces_per_cell;
                         line++)
                      if (cell->face(face)->line(line)->has_children())
                        cell->face(face)->line(line)->child(0)->vertex(1) =
                          (cell->face(face)->line(line)->vertex(0) +
                           cell->face(face)->line(line)->vertex(1)) /
                          2.0;
                  }
              }
      }

    // Make sure FEValues notices that the mesh has changed
    triangulation.signals.mesh_movement();
  }



  template <typename MeshType>
  DEAL_II_CXX20_REQUIRES(concepts::is_triangulation_or_dof_handler<MeshType>)
  std::vector<typename MeshType::active_cell_iterator> get_active_child_cells(
    const typename MeshType::cell_iterator &cell)
  {
    std::vector<typename MeshType::active_cell_iterator> child_cells;

    if (cell->has_children())
      {
        for (unsigned int child = 0; child < cell->n_children(); ++child)
          if (cell->child(child)->has_children())
            {
              const std::vector<typename MeshType::active_cell_iterator>
                children = get_active_child_cells<MeshType>(cell->child(child));
              child_cells.insert(child_cells.end(),
                                 children.begin(),
                                 children.end());
            }
          else
            child_cells.push_back(cell->child(child));
      }

    return child_cells;
  }



  template <typename MeshType>
  DEAL_II_CXX20_REQUIRES(concepts::is_triangulation_or_dof_handler<MeshType>)
  void get_active_neighbors(
    const typename MeshType::active_cell_iterator        &cell,
    std::vector<typename MeshType::active_cell_iterator> &active_neighbors)
  {
    active_neighbors.clear();
    for (const unsigned int n : cell->face_indices())
      if (!cell->at_boundary(n))
        {
          if (MeshType::dimension == 1)
            {
              // check children of neighbor. note
              // that in 1d children of the neighbor
              // may be further refined. In 1d the
              // case is simple since we know what
              // children bound to the present cell
              typename MeshType::cell_iterator neighbor_child =
                cell->neighbor(n);
              if (!neighbor_child->is_active())
                {
                  while (neighbor_child->has_children())
                    neighbor_child = neighbor_child->child(n == 0 ? 1 : 0);

                  Assert(neighbor_child->neighbor(n == 0 ? 1 : 0) == cell,
                         ExcInternalError());
                }
              active_neighbors.push_back(neighbor_child);
            }
          else
            {
              if (cell->face(n)->has_children())
                // this neighbor has children. find
                // out which border to the present
                // cell
                for (unsigned int c = 0;
                     c < cell->face(n)->n_active_descendants();
                     ++c)
                  active_neighbors.push_back(
                    cell->neighbor_child_on_subface(n, c));
              else
                {
                  // the neighbor must be active
                  // himself
                  Assert(cell->neighbor(n)->is_active(), ExcInternalError());
                  active_neighbors.push_back(cell->neighbor(n));
                }
            }
        }
  }



  template <typename CellIterator>
  std::size_t
  PeriodicFacePair<CellIterator>::memory_consumption() const
  {
    return sizeof(*this) + matrix.memory_consumption();
  }



  namespace internal
  {
    template <typename DataType,
              typename MeshType,
              typename MeshCellIteratorType>
    DEAL_II_CXX20_REQUIRES(concepts::is_triangulation_or_dof_handler<MeshType>)
    inline void exchange_cell_data(
      const MeshType &mesh,
      const std::function<std::optional<DataType>(const MeshCellIteratorType &)>
        &pack,
      const std::function<void(const MeshCellIteratorType &, const DataType &)>
                                                                 &unpack,
      const std::function<bool(const MeshCellIteratorType &)>    &cell_filter,
      const std::function<void(
        const std::function<void(const MeshCellIteratorType &,
                                 const types::subdomain_id)> &)> &process_cells,
      const std::function<std::set<types::subdomain_id>(
        const parallel::TriangulationBase<MeshType::dimension,
                                          MeshType::space_dimension> &)>
        &compute_ghost_owners)
    {
#  ifndef DEAL_II_WITH_MPI
      (void)mesh;
      (void)pack;
      (void)unpack;
      (void)cell_filter;
      (void)process_cells;
      (void)compute_ghost_owners;
      Assert(false, ExcNeedsMPI());
#  else
      constexpr int dim      = MeshType::dimension;
      constexpr int spacedim = MeshType::space_dimension;
      auto          tria =
        dynamic_cast<const parallel::TriangulationBase<dim, spacedim> *>(
          &mesh.get_triangulation());
      Assert(
        tria != nullptr,
        ExcMessage(
          "The function exchange_cell_data_to_ghosts() only works with parallel triangulations."));

      if (const auto tria = dynamic_cast<
            const parallel::shared::Triangulation<dim, spacedim> *>(
            &mesh.get_triangulation()))
        {
          Assert(
            tria->with_artificial_cells(),
            ExcMessage(
              "The functions GridTools::exchange_cell_data_to_ghosts() and "
              "GridTools::exchange_cell_data_to_level_ghosts() can only "
              "operate on a single layer of ghost cells. However, you have "
              "given a Triangulation object of type "
              "parallel::shared::Triangulation without artificial cells "
              "resulting in an arbitrary number of ghost layers. "
              "To use this function for a Triangulation object of type "
              "parallel::shared::Triangulation, make sure to create the "
              "Triangulation object with allow_artificial_cells set to true. "
              "This results in a parallel::shared::Triangulation with only "
              "a single layer of ghost cells."));
        }

      // build list of cells to request for each neighbor
      std::set<dealii::types::subdomain_id> ghost_owners =
        compute_ghost_owners(*tria);
      std::map<dealii::types::subdomain_id,
               std::vector<typename CellId::binary_type>>
        neighbor_cell_list;

      for (const auto ghost_owner : ghost_owners)
        neighbor_cell_list[ghost_owner] = {};

      process_cells([&](const auto &cell, const auto key) -> void {
        if (cell_filter(cell))
          {
            constexpr int spacedim = MeshType::space_dimension;
            neighbor_cell_list[key].emplace_back(
              cell->id().template to_binary<spacedim>());
          }
      });

      Assert(ghost_owners.size() == neighbor_cell_list.size(),
             ExcInternalError());


      // Before sending & receiving, make sure we protect this section with
      // a mutex:
      static Utilities::MPI::CollectiveMutex      mutex;
      Utilities::MPI::CollectiveMutex::ScopedLock lock(
        mutex, tria->get_mpi_communicator());

      const int mpi_tag =
        Utilities::MPI::internal::Tags::exchange_cell_data_request;
      const int mpi_tag_reply =
        Utilities::MPI::internal::Tags::exchange_cell_data_reply;

      // send our requests
      std::vector<MPI_Request> requests(ghost_owners.size());
      {
        unsigned int idx = 0;
        for (const auto &it : neighbor_cell_list)
          {
            // send the data about the relevant cells
            const int ierr = MPI_Isend(it.second.data(),
                                       it.second.size() * sizeof(it.second[0]),
                                       MPI_BYTE,
                                       it.first,
                                       mpi_tag,
                                       tria->get_mpi_communicator(),
                                       &requests[idx]);
            AssertThrowMPI(ierr);
            ++idx;
          }
      }

      // receive requests and reply with the results
      std::vector<MPI_Request>       reply_requests(ghost_owners.size());
      std::vector<std::vector<char>> sendbuffers(ghost_owners.size());

      for (unsigned int idx = 0; idx < ghost_owners.size(); ++idx)
        {
          MPI_Status status;
          int        ierr = MPI_Probe(MPI_ANY_SOURCE,
                               mpi_tag,
                               tria->get_mpi_communicator(),
                               &status);
          AssertThrowMPI(ierr);

          int len;
          ierr = MPI_Get_count(&status, MPI_BYTE, &len);
          AssertThrowMPI(ierr);
          Assert(len % sizeof(typename CellId::binary_type) == 0,
                 ExcInternalError());

          const unsigned int n_cells =
            len / sizeof(typename CellId::binary_type);
          std::vector<typename CellId::binary_type> cells_with_requests(
            n_cells);
          std::vector<DataType> data_to_send;
          data_to_send.reserve(n_cells);
          std::vector<bool> cell_carries_data(n_cells, false);

          ierr = MPI_Recv(cells_with_requests.data(),
                          len,
                          MPI_BYTE,
                          status.MPI_SOURCE,
                          status.MPI_TAG,
                          tria->get_mpi_communicator(),
                          &status);
          AssertThrowMPI(ierr);

          // store data for each cell
          for (unsigned int c = 0; c < static_cast<unsigned int>(n_cells); ++c)
            {
              const auto cell =
                tria->create_cell_iterator(CellId(cells_with_requests[c]));

              MeshCellIteratorType mesh_it(tria,
                                           cell->level(),
                                           cell->index(),
                                           &mesh);

              const std::optional<DataType> data = pack(mesh_it);
              if (data)
                {
                  data_to_send.emplace_back(*data);
                  cell_carries_data[c] = true;
                }
            }

          // collect data for sending the reply in a buffer

          // (a) make room for storing the local offsets in case we receive
          // other data
          sendbuffers[idx].resize(sizeof(std::size_t));

          // (b) append the actual data and store how much memory it
          // corresponds to, which we then insert into the leading position of
          // the sendbuffer
          std::size_t size_of_send =
            Utilities::pack(data_to_send,
                            sendbuffers[idx],
                            /*enable_compression*/ false);
          std::memcpy(sendbuffers[idx].data(),
                      &size_of_send,
                      sizeof(std::size_t));

          // (c) append information of certain cells that got left out in case
          // we need it
          if (data_to_send.size() < n_cells)
            Utilities::pack(cell_carries_data,
                            sendbuffers[idx],
                            /*enable_compression*/ false);

          // send data
          ierr = MPI_Isend(sendbuffers[idx].data(),
                           sendbuffers[idx].size(),
                           MPI_BYTE,
                           status.MPI_SOURCE,
                           mpi_tag_reply,
                           tria->get_mpi_communicator(),
                           &reply_requests[idx]);
          AssertThrowMPI(ierr);
        }

      // finally receive the replies
      std::vector<char> receive;
      for (unsigned int id = 0; id < neighbor_cell_list.size(); ++id)
        {
          MPI_Status status;
          int        ierr = MPI_Probe(MPI_ANY_SOURCE,
                               mpi_tag_reply,
                               tria->get_mpi_communicator(),
                               &status);
          AssertThrowMPI(ierr);

          int len;
          ierr = MPI_Get_count(&status, MPI_BYTE, &len);
          AssertThrowMPI(ierr);

          receive.resize(len);

          ierr = MPI_Recv(receive.data(),
                          len,
                          MPI_BYTE,
                          status.MPI_SOURCE,
                          status.MPI_TAG,
                          tria->get_mpi_communicator(),
                          &status);
          AssertThrowMPI(ierr);

          // (a) first determine the length of the data section in the
          // received buffer
          auto        data_iterator = receive.begin();
          std::size_t size_of_received_data =
            Utilities::unpack<std::size_t>(data_iterator,
                                           data_iterator + sizeof(std::size_t));
          data_iterator += sizeof(std::size_t);

          // (b) unpack the data section in the indicated region
          auto received_data = Utilities::unpack<std::vector<DataType>>(
            data_iterator,
            data_iterator + size_of_received_data,
            /*enable_compression*/ false);
          data_iterator += size_of_received_data;

          // (c) check if the received data contained fewer entries than the
          // number of cells we identified in the beginning, in which case we
          // need to extract the boolean vector with the relevant information
          const std::vector<typename CellId::binary_type> &this_cell_list =
            neighbor_cell_list[status.MPI_SOURCE];
          AssertIndexRange(received_data.size(), this_cell_list.size() + 1);
          std::vector<bool> cells_with_data;
          if (received_data.size() < this_cell_list.size())
            {
              cells_with_data = Utilities::unpack<std::vector<bool>>(
                data_iterator, receive.end(), /*enable_compression*/ false);
              AssertDimension(cells_with_data.size(), this_cell_list.size());
            }

          // (d) go through the received data and call the user-provided
          // unpack function
          auto received_data_iterator = received_data.begin();
          for (unsigned int c = 0; c < this_cell_list.size(); ++c)
            if (cells_with_data.empty() || cells_with_data[c])
              {
                const typename Triangulation<dim, spacedim>::cell_iterator
                  tria_cell = tria->create_cell_iterator(this_cell_list[c]);

                MeshCellIteratorType cell(tria,
                                          tria_cell->level(),
                                          tria_cell->index(),
                                          &mesh);

                unpack(cell, *received_data_iterator);
                ++received_data_iterator;
              }
        }

      // make sure that all communication is finished
      // when we leave this function.
      if (requests.size() > 0)
        {
          const int ierr =
            MPI_Waitall(requests.size(), requests.data(), MPI_STATUSES_IGNORE);
          AssertThrowMPI(ierr);
        }
      if (reply_requests.size() > 0)
        {
          const int ierr = MPI_Waitall(reply_requests.size(),
                                       reply_requests.data(),
                                       MPI_STATUSES_IGNORE);
          AssertThrowMPI(ierr);
        }


#  endif // DEAL_II_WITH_MPI
    }

  } // namespace internal

  template <typename DataType, typename MeshType>
  DEAL_II_CXX20_REQUIRES(concepts::is_triangulation_or_dof_handler<MeshType>)
  inline void exchange_cell_data_to_ghosts(
    const MeshType                                      &mesh,
    const std::function<std::optional<DataType>(
      const typename MeshType::active_cell_iterator &)> &pack,
    const std::function<void(const typename MeshType::active_cell_iterator &,
                             const DataType &)>         &unpack,
    const std::function<bool(const typename MeshType::active_cell_iterator &)>
      &cell_filter)
  {
#  ifndef DEAL_II_WITH_MPI
    (void)mesh;
    (void)pack;
    (void)unpack;
    (void)cell_filter;
    Assert(false, ExcNeedsMPI());
#  else
    internal::exchange_cell_data<DataType,
                                 MeshType,
                                 typename MeshType::active_cell_iterator>(
      mesh,
      pack,
      unpack,
      cell_filter,
      [&](const auto &process) {
        for (const auto &cell : mesh.active_cell_iterators())
          if (cell->is_ghost())
            process(cell, cell->subdomain_id());
      },
      [](const auto &tria) { return tria.ghost_owners(); });
#  endif
  }



  template <typename DataType, typename MeshType>
  DEAL_II_CXX20_REQUIRES(concepts::is_triangulation_or_dof_handler<MeshType>)
  inline void exchange_cell_data_to_level_ghosts(
    const MeshType                                     &mesh,
    const std::function<std::optional<DataType>(
      const typename MeshType::level_cell_iterator &)> &pack,
    const std::function<void(const typename MeshType::level_cell_iterator &,
                             const DataType &)>        &unpack,
    const std::function<bool(const typename MeshType::level_cell_iterator &)>
      &cell_filter)
  {
#  ifndef DEAL_II_WITH_MPI
    (void)mesh;
    (void)pack;
    (void)unpack;
    (void)cell_filter;
    Assert(false, ExcNeedsMPI());
#  else
    internal::exchange_cell_data<DataType,
                                 MeshType,
                                 typename MeshType::level_cell_iterator>(
      mesh,
      pack,
      unpack,
      cell_filter,
      [&](const auto &process) {
        for (const auto &cell : mesh.cell_iterators())
          if (cell->is_ghost_on_level())
            process(cell, cell->level_subdomain_id());
      },
      [](const auto &tria) { return tria.level_ghost_owners(); });
#  endif
  }
} // namespace GridTools

#endif // DOXYGEN

DEAL_II_NAMESPACE_CLOSE

#endif
