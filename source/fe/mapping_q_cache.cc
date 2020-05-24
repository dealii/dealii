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

#include <deal.II/base/memory_consumption.h>
#include <deal.II/base/work_stream.h>

#include <deal.II/fe/mapping_q_cache.h>

#include <functional>

DEAL_II_NAMESPACE_OPEN

template <int dim, int spacedim>
MappingQCache<dim, spacedim>::MappingQCache(
  const unsigned int polynomial_degree)
  : MappingQGeneric<dim, spacedim>(polynomial_degree)
{}



template <int dim, int spacedim>
MappingQCache<dim, spacedim>::MappingQCache(
  const MappingQCache<dim, spacedim> &mapping)
  : MappingQGeneric<dim, spacedim>(mapping)
  , support_point_cache(mapping.support_point_cache)
{}



template <int dim, int spacedim>
MappingQCache<dim, spacedim>::~MappingQCache()
{
  // When this object goes out of scope, we want the cache to get cleared and
  // free its memory before the signal is disconnected in order to not work on
  // invalid memory that has been left back by freeing an object of this
  // class.
  support_point_cache.reset();
  clear_signal.disconnect();
}



template <int dim, int spacedim>
std::unique_ptr<Mapping<dim, spacedim>>
MappingQCache<dim, spacedim>::clone() const
{
  return std::make_unique<MappingQCache<dim, spacedim>>(*this);
}



template <int dim, int spacedim>
bool
MappingQCache<dim, spacedim>::preserves_vertex_locations() const
{
  return false;
}



template <int dim, int spacedim>
void
MappingQCache<dim, spacedim>::initialize(
  const Triangulation<dim, spacedim> &  triangulation,
  const MappingQGeneric<dim, spacedim> &mapping)
{
  AssertDimension(this->get_degree(), mapping.get_degree());
  initialize(triangulation,
             [&mapping](const typename Triangulation<dim>::cell_iterator &cell)
               -> std::vector<Point<spacedim>> {
               return mapping.compute_mapping_support_points(cell);
             });
}



template <int dim, int spacedim>
void
MappingQCache<dim, spacedim>::initialize(
  const Triangulation<dim, spacedim> &triangulation,
  const std::function<std::vector<Point<spacedim>>(
    const typename Triangulation<dim, spacedim>::cell_iterator &)>
    &compute_points_on_cell)
{
  clear_signal = triangulation.signals.any_change.connect(
    [&]() -> void { this->support_point_cache.reset(); });

  support_point_cache =
    std::make_shared<std::vector<std::vector<std::vector<Point<spacedim>>>>>(
      triangulation.n_levels());
  for (unsigned int l = 0; l < triangulation.n_levels(); ++l)
    (*support_point_cache)[l].resize(triangulation.n_raw_cells(l));

  WorkStream::run(
    triangulation.begin(),
    triangulation.end(),
    [&](const typename Triangulation<dim, spacedim>::cell_iterator &cell,
        void *,
        void *) {
      (*support_point_cache)[cell->level()][cell->index()] =
        compute_points_on_cell(cell);
      AssertDimension(
        (*support_point_cache)[cell->level()][cell->index()].size(),
        Utilities::pow(this->get_degree() + 1, dim));
    },
    /* copier */ std::function<void(void *)>(),
    /* scratch_data */ nullptr,
    /* copy_data */ nullptr,
    2 * MultithreadInfo::n_threads(),
    /* chunk_size = */ 1);
}



template <int dim, int spacedim>
std::size_t
MappingQCache<dim, spacedim>::memory_consumption() const
{
  if (support_point_cache.get() != nullptr)
    return sizeof(*this) +
           MemoryConsumption::memory_consumption(*support_point_cache);
  else
    return sizeof(*this);
}



template <int dim, int spacedim>
std::vector<Point<spacedim>>
MappingQCache<dim, spacedim>::compute_mapping_support_points(
  const typename Triangulation<dim, spacedim>::cell_iterator &cell) const
{
  Assert(support_point_cache.get() != nullptr,
         ExcMessage("Must call MappingQCache::initialize() before "
                    "using it or after mesh has changed!"));

  AssertIndexRange(cell->level(), support_point_cache->size());
  AssertIndexRange(cell->index(), (*support_point_cache)[cell->level()].size());
  return (*support_point_cache)[cell->level()][cell->index()];
}



//--------------------------- Explicit instantiations -----------------------
#include "mapping_q_cache.inst"


DEAL_II_NAMESPACE_CLOSE
