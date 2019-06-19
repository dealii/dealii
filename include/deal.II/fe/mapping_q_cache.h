// ---------------------------------------------------------------------
//
// Copyright (C) 2019 by the deal.II authors
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

#ifndef dealii_mapping_q_cache_h
#define dealii_mapping_q_cache_h


#include <deal.II/base/config.h>

#include <deal.II/fe/mapping_q_generic.h>

#include <deal.II/grid/tria.h>


DEAL_II_NAMESPACE_OPEN



/*!@addtogroup mapping */
/*@{*/


/**
 * This class implements a caching strategy for objects of the MappingQ family
 * in terms of the MappingQGeneric::compute_mapping_support_points() function,
 * which is used in all operations of MappingQGeneric. The information of the
 * mapping is pre-computed by the MappingQCache::initialize() function.
 *
 * @author Martin Kronbichler, 2019
 */
template <int dim, int spacedim = dim>
class MappingQCache : public MappingQGeneric<dim, spacedim>
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
   * cells (on all levels) of the given triangulation. Note that the cache is
   * invalidated upon the signal Triangulation::Signals::any_change of the
   * underlying triangulation.
   */
  void
  initialize(const Triangulation<dim, spacedim> &  triangulation,
             const MappingQGeneric<dim, spacedim> &mapping);

  /**
   * Return the memory consumption (in bytes) of the cache.
   */
  std::size_t
  memory_consumption() const;

protected:
  /**
   * This is the main function overriden from the base class MappingQGeneric.
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
};

DEAL_II_NAMESPACE_CLOSE

#endif
