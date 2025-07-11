// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2019 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_mapping_q_cache_h
#define dealii_mapping_q_cache_h


#include <deal.II/base/config.h>

#include <deal.II/base/function.h>
#include <deal.II/base/mg_level_object.h>

#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/tria.h>

#include <boost/container/small_vector.hpp>
#include <boost/signals2/connection.hpp>


DEAL_II_NAMESPACE_OPEN

// Forward declarations
#ifndef DOXYGEN
template <int dim, int spacedim>
DEAL_II_CXX20_REQUIRES((concepts::is_valid_dim_spacedim<dim, spacedim>))
class DoFHandler;
#endif


/**
 * @addtogroup mapping
 * @{
 */

/**
 * This class implements a caching strategy for objects of the MappingQ family
 * in terms of the MappingQ::compute_mapping_support_points() function,
 * which is used in all operations of MappingQ. The information of the
 * mapping is pre-computed by the MappingQCache::initialize() function.
 *
 * The use of this class is discussed extensively in step-65.
 */
template <int dim, int spacedim = dim>
class MappingQCache : public MappingQ<dim, spacedim>
{
public:
  /**
   * Constructor. @p polynomial_degree denotes the polynomial degree of the
   * polynomials that are used to map cells from the reference to the real
   * cell.
   */
  explicit MappingQCache(const unsigned int polynomial_degree);

  /**
   * Copy constructor.
   */
  explicit MappingQCache(const MappingQCache<dim, spacedim> &mapping);

  /**
   * Destructor.
   */
  ~MappingQCache();

  /**
   * clone() functionality. For documentation, see Mapping::clone().
   */
  virtual std::unique_ptr<Mapping<dim, spacedim>>
  clone() const override;

  /**
   * Returns @p false because the preservation of vertex locations depends on
   * the mapping handed to the reinit() function.
   */
  virtual bool
  preserves_vertex_locations() const override;

  /**
   * Initialize the data cache by computing the mapping support points for all
   * cells (on all levels) of the given triangulation.
   *
   * @note The cache is invalidated upon the signal
   * Triangulation::Signals::any_change of the underlying triangulation.
   */
  void
  initialize(const Mapping<dim, spacedim>       &mapping,
             const Triangulation<dim, spacedim> &triangulation);

  /**
   * Initialize the data cache by letting the function given as an argument
   * provide the mapping support points for all cells (on all levels) of the
   * given triangulation. The function must return a vector of
   * `Point<spacedim>` whose length is the same as the size of the polynomial
   * space, $(p+1)^\text{dim}$, where $p$ is the polynomial degree of the
   * mapping, and it must be in the order the mapping or FE_Q sort their
   * points, i.e., all $2^\text{dim}$ vertex points first, then the points on
   * the lines, quads, and hexes according to the usual hierarchical
   * numbering. No attempt is made to validate these points internally, except
   * for the number of given points.
   *
   * @note If multiple threads are enabled, this function will run in
   * parallel, invoking the function passed in several times. Thus, in case
   * MultithreadInfo::n_threads()>1, the user code must make sure that the
   * function, typically a lambda, does not write into data shared with other
   * threads.
   *
   * @note The cache is invalidated upon the signal
   * Triangulation::Signals::any_change of the underlying triangulation.
   */
  void
  initialize(const Triangulation<dim, spacedim> &triangulation,
             const std::function<std::vector<Point<spacedim>>(
               const typename Triangulation<dim, spacedim>::cell_iterator &)>
               &compute_points_on_cell);

  /**
   * Initialize the data cache by computing the mapping support points for all
   * cells (on all levels) of the given triangulation and a given @p mapping
   * and transforming these points via the function @p transformation_function.
   *
   * The bool @p function_describes_relative_displacement indicates that
   * the function @p transformation_function maps to absolute coordinates.
   * If the parameter is set to true, the return value of the function is
   * interpreted as relative deformation and the result is eventually added
   * to the original point for the support points eventually used by this class.
   *
   * This function calls the previous function so the comments regarding
   * threading listed above apply also here.
   *
   * @note The cache is invalidated upon the signal
   * Triangulation::Signals::any_change of the underlying triangulation.
   */
  void
  initialize(const Mapping<dim, spacedim>       &mapping,
             const Triangulation<dim, spacedim> &tria,
             const std::function<Point<spacedim>(
               const typename Triangulation<dim, spacedim>::cell_iterator &,
               const Point<spacedim> &)>        &transformation_function,
             const bool function_describes_relative_displacement);

  /**
   * The same as above but taking a dealii::Function object.
   */
  void
  initialize(const Mapping<dim, spacedim>       &mapping,
             const Triangulation<dim, spacedim> &tria,
             const Function<spacedim>           &transformation_function,
             const bool function_describes_relative_displacement);

  /**
   * Initialize the data cache of the active cells by a discrete field
   * (specified
   * by @p dof_handler and @p vector) that describes the absolute or
   * relative position of each support point.
   *
   * @note By using this function for reinitialization, this class behaves like
   *   MappingFEField (vector_describes_relative_displacement == false) or
   *   MappingQEulerian (vector_describes_relative_displacement == true), but
   *   with much more efficient operations internally.
   */
  template <typename VectorType>
  void
  initialize(const Mapping<dim, spacedim>    &mapping,
             const DoFHandler<dim, spacedim> &dof_handler,
             const VectorType                &vector,
             const bool vector_describes_relative_displacement);

  /**
   * Initialize the data cache of all non-artificial cells by a solution
   * (specified by @p dof_handler and a set of @p vectors on all levels of the
   * triangulation) that describes the absolute or relative position of each
   * support point.
   *
   * @note By using this function for reinitialization, this class behaves like
   *   MappingFEField (vector_describes_relative_displacement == false) or
   *   MappingQEulerian (vector_describes_relative_displacement == true), but
   *   with much more efficient operations internally.
   */
  template <typename VectorType>
  void
  initialize(const Mapping<dim, spacedim>    &mapping,
             const DoFHandler<dim, spacedim> &dof_handler,
             const MGLevelObject<VectorType> &vectors,
             const bool vector_describes_relative_displacement);

  /**
   * @copydoc Mapping::get_vertices()
   */
  virtual boost::container::small_vector<Point<spacedim>,
#ifndef _MSC_VER
                                         ReferenceCells::max_n_vertices<dim>()
#else
                                         GeometryInfo<dim>::vertices_per_cell
#endif
                                         >
  get_vertices(const typename Triangulation<dim, spacedim>::cell_iterator &cell)
    const override;

  /**
   * Return the memory consumption (in bytes) of the cache.
   */
  std::size_t
  memory_consumption() const;

protected:
  /**
   * This is the main function overridden from the base class MappingQ.
   */
  virtual std::vector<Point<spacedim>>
  compute_mapping_support_points(
    const typename Triangulation<dim, spacedim>::cell_iterator &cell)
    const override;

private:
  /**
   * The point cache filled upon calling initialize(). It is made a shared
   * pointer to allow several instances (created via clone()) to share this
   * cache.
   */
  std::shared_ptr<std::vector<std::vector<std::vector<Point<spacedim>>>>>
    support_point_cache;

  /**
   * The connection to Triangulation::signals::any that must be reset once
   * this class goes out of scope.
   */
  boost::signals2::connection clear_signal;

  /**
   * Specifies if support_point_cache has been set up for the cells on the
   * levels.
   */
  bool uses_level_info;
};

/** @} */

DEAL_II_NAMESPACE_CLOSE

#endif
