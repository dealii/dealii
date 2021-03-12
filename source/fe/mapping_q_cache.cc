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

#include <deal.II/fe/fe_nothing.h>
#include <deal.II/fe/fe_tools.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q1.h>
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
  const Mapping<dim, spacedim> &      mapping,
  const Triangulation<dim, spacedim> &triangulation)
{
  // FE and FEValues in the case they are needed
  FE_Nothing<dim, spacedim> fe;
  Threads::ThreadLocalStorage<std::unique_ptr<FEValues<dim, spacedim>>>
    fe_values_all;

  this->initialize(
    triangulation,
    [&](const typename Triangulation<dim, spacedim>::cell_iterator &cell) {
      const auto mapping_q_generic =
        dynamic_cast<const MappingQGeneric<dim, spacedim> *>(&mapping);
      if (mapping_q_generic != nullptr &&
          this->get_degree() == mapping_q_generic->get_degree())
        {
          return mapping_q_generic->compute_mapping_support_points(cell);
        }
      else
        {
          // get FEValues (thread-safe); in the case that this thread has not
          // created a an FEValues object yet, this helper-function also
          // creates one with the right quadrature rule
          auto &fe_values = fe_values_all.get();
          if (fe_values.get() == nullptr)
            {
              QGaussLobatto<dim> quadrature_gl(this->polynomial_degree + 1);

              std::vector<Point<dim>> quadrature_points;
              for (const auto i :
                   FETools::hierarchic_to_lexicographic_numbering<dim>(
                     this->polynomial_degree))
                quadrature_points.push_back(quadrature_gl.point(i));
              Quadrature<dim> quadrature(quadrature_points);

              fe_values = std::make_unique<FEValues<dim, spacedim>>(
                mapping, fe, quadrature, update_quadrature_points);
            }

          fe_values->reinit(cell);
          return fe_values->get_quadrature_points();
        }
    });
}



template <int dim, int spacedim>
void
MappingQCache<dim, spacedim>::initialize(
  const Triangulation<dim, spacedim> &  triangulation,
  const MappingQGeneric<dim, spacedim> &mapping)
{
  this->initialize(mapping, triangulation);
}



template <int dim, int spacedim>
void
MappingQCache<dim, spacedim>::initialize(
  const Triangulation<dim, spacedim> &triangulation,
  const std::function<std::vector<Point<spacedim>>(
    const typename Triangulation<dim, spacedim>::cell_iterator &)>
    &compute_points_on_cell)
{
  clear_signal.disconnect();
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
void
MappingQCache<dim, spacedim>::initialize(
  const Mapping<dim, spacedim> &      mapping,
  const Triangulation<dim, spacedim> &tria,
  const std::function<Point<spacedim>(
    const typename Triangulation<dim, spacedim>::cell_iterator &,
    const Point<spacedim> &)> &       transformation_function,
  const bool                          function_describes_relative_displacement)
{
  // FE and FEValues in the case they are needed
  FE_Nothing<dim, spacedim> fe;
  Threads::ThreadLocalStorage<std::unique_ptr<FEValues<dim, spacedim>>>
    fe_values_all;

  this->initialize(
    tria,
    [&](const typename Triangulation<dim, spacedim>::cell_iterator &cell) {
      std::vector<Point<spacedim>> points;

      const auto mapping_q_generic =
        dynamic_cast<const MappingQGeneric<dim, spacedim> *>(&mapping);

      if (mapping_q_generic != nullptr &&
          this->get_degree() == mapping_q_generic->get_degree())
        {
          points = mapping_q_generic->compute_mapping_support_points(cell);
        }
      else
        {
          // get FEValues (thread-safe); in the case that this thread has not
          // created a an FEValues object yet, this helper-function also
          // creates one with the right quadrature rule
          auto &fe_values = fe_values_all.get();
          if (fe_values.get() == nullptr)
            {
              QGaussLobatto<dim> quadrature_gl(this->polynomial_degree + 1);

              std::vector<Point<dim>> quadrature_points;
              for (const auto i :
                   FETools::hierarchic_to_lexicographic_numbering<dim>(
                     this->polynomial_degree))
                quadrature_points.push_back(quadrature_gl.point(i));
              Quadrature<dim> quadrature(quadrature_points);

              fe_values = std::make_unique<FEValues<dim, spacedim>>(
                mapping, fe, quadrature, update_quadrature_points);
            }

          fe_values->reinit(cell);
          points = fe_values->get_quadrature_points();
        }

      for (auto &p : points)
        if (function_describes_relative_displacement)
          p += transformation_function(cell, p);
        else
          p = transformation_function(cell, p);

      return points;
    });
}



template <int dim, int spacedim>
void
MappingQCache<dim, spacedim>::initialize(
  const Mapping<dim, spacedim> &      mapping,
  const Triangulation<dim, spacedim> &tria,
  const Function<spacedim> &          transformation_function,
  const bool                          function_describes_relative_displacement)
{
  AssertDimension(transformation_function.n_components, spacedim);

  this->initialize(
    mapping,
    tria,
    [&](const typename Triangulation<dim, spacedim>::cell_iterator &,
        const Point<spacedim> &point) -> Point<spacedim> {
      Point<spacedim> new_point;
      for (int c = 0; c < spacedim; ++c)
        new_point[c] = transformation_function.value(point, c);
      return new_point;
    },
    function_describes_relative_displacement);
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
