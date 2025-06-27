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

#include <deal.II/base/memory_consumption.h>
#include <deal.II/base/utilities.h>
#include <deal.II/base/work_stream.h>

#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_nothing.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_tools.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q_cache.h>

#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/lac/petsc_vector.h>
#include <deal.II/lac/trilinos_vector.h>

#include <functional>

DEAL_II_NAMESPACE_OPEN

template <int dim, int spacedim>
MappingQCache<dim, spacedim>::MappingQCache(
  const unsigned int polynomial_degree)
  : MappingQ<dim, spacedim>(polynomial_degree)
  , uses_level_info(false)
{}



template <int dim, int spacedim>
MappingQCache<dim, spacedim>::MappingQCache(
  const MappingQCache<dim, spacedim> &mapping)
  : MappingQ<dim, spacedim>(mapping)
  , support_point_cache(mapping.support_point_cache)
  , uses_level_info(mapping.uses_level_info)
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
  const Mapping<dim, spacedim>       &mapping,
  const Triangulation<dim, spacedim> &triangulation)
{
  // FE and FEValues in the case they are needed
  FE_Nothing<dim, spacedim> fe;
  Threads::ThreadLocalStorage<std::unique_ptr<FEValues<dim, spacedim>>>
    fe_values_all;

  this->initialize(
    triangulation,
    [&](const typename Triangulation<dim, spacedim>::cell_iterator &cell) {
      const auto mapping_q =
        dynamic_cast<const MappingQ<dim, spacedim> *>(&mapping);
      if (mapping_q != nullptr && this->get_degree() == mapping_q->get_degree())
        {
          return mapping_q->compute_mapping_support_points(cell);
        }
      else
        {
          // get FEValues (thread-safe); in the case that this thread has not
          // created a an FEValues object yet, this helper-function also
          // creates one with the right quadrature rule
          auto &fe_values = fe_values_all.get();
          if (fe_values.get() == nullptr)
            {
              const QGaussLobatto<dim> quadrature_gl(this->polynomial_degree +
                                                     1);

              std::vector<Point<dim>> quadrature_points;
              for (const auto i :
                   FETools::hierarchic_to_lexicographic_numbering<dim>(
                     this->polynomial_degree))
                quadrature_points.push_back(quadrature_gl.point(i));
              const Quadrature<dim> quadrature(quadrature_points);

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
      // Do not use `this` in Assert because nvcc when using C++20 assumes that
      // `this` is an integer and we get the following error: invalid type
      // argument of unary '*' (have 'int')
      [[maybe_unused]] const unsigned int d = this->get_degree() + 1;
      AssertDimension(
        (*support_point_cache)[cell->level()][cell->index()].size(),
        Utilities::pow(d, dim));
    },
    /* copier */ std::function<void(void *)>(),
    /* scratch_data */ nullptr,
    /* copy_data */ nullptr,
    2 * MultithreadInfo::n_threads(),
    /* chunk_size = */ 1);

  uses_level_info = true;
}



template <int dim, int spacedim>
void
MappingQCache<dim, spacedim>::initialize(
  const Mapping<dim, spacedim>       &mapping,
  const Triangulation<dim, spacedim> &tria,
  const std::function<Point<spacedim>(
    const typename Triangulation<dim, spacedim>::cell_iterator &,
    const Point<spacedim> &)>        &transformation_function,
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

      const auto mapping_q =
        dynamic_cast<const MappingQ<dim, spacedim> *>(&mapping);

      if (mapping_q != nullptr && this->get_degree() == mapping_q->get_degree())
        {
          points = mapping_q->compute_mapping_support_points(cell);
        }
      else
        {
          // get FEValues (thread-safe); in the case that this thread has not
          // created a an FEValues object yet, this helper-function also
          // creates one with the right quadrature rule
          auto &fe_values = fe_values_all.get();
          if (fe_values.get() == nullptr)
            {
              const QGaussLobatto<dim> quadrature_gl(this->polynomial_degree +
                                                     1);

              std::vector<Point<dim>> quadrature_points;
              for (const auto i :
                   FETools::hierarchic_to_lexicographic_numbering<dim>(
                     this->polynomial_degree))
                quadrature_points.push_back(quadrature_gl.point(i));
              const Quadrature<dim> quadrature(quadrature_points);

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

  uses_level_info = true;
}



template <int dim, int spacedim>
void
MappingQCache<dim, spacedim>::initialize(
  const Mapping<dim, spacedim>       &mapping,
  const Triangulation<dim, spacedim> &tria,
  const Function<spacedim>           &transformation_function,
  const bool                          function_describes_relative_displacement)
{
  AssertDimension(transformation_function.n_components, spacedim);

  this->initialize(
    mapping,
    tria,
    [&](const auto &, const auto &point) {
      Point<spacedim> new_point;
      for (unsigned int c = 0; c < spacedim; ++c)
        new_point[c] = transformation_function.value(point, c);
      return new_point;
    },
    function_describes_relative_displacement);

  uses_level_info = true;
}



namespace
{
  template <typename VectorType>
  void
  copy_locally_owned_data_from(
    const VectorType &vector,
    LinearAlgebra::distributed::Vector<typename VectorType::value_type>
      &vector_ghosted)
  {
    LinearAlgebra::ReadWriteVector<typename VectorType::value_type> temp;
    temp.reinit(vector.locally_owned_elements());
    temp.import_elements(vector, VectorOperation::insert);
    vector_ghosted.import_elements(temp, VectorOperation::insert);
  }
} // namespace



template <int dim, int spacedim>
template <typename VectorType>
void
MappingQCache<dim, spacedim>::initialize(
  const Mapping<dim, spacedim>    &mapping,
  const DoFHandler<dim, spacedim> &dof_handler,
  const VectorType                &vector,
  const bool                       vector_describes_relative_displacement)
{
  AssertDimension(dof_handler.get_fe_collection().size(), 1);
  const FiniteElement<dim, spacedim> &fe = dof_handler.get_fe();
  AssertDimension(fe.n_base_elements(), 1);
  AssertDimension(fe.element_multiplicity(0), spacedim);

  const unsigned int is_fe_q =
    dynamic_cast<const FE_Q<dim, spacedim> *>(&fe.base_element(0)) != nullptr;
  const unsigned int is_fe_dgq =
    dynamic_cast<const FE_DGQ<dim, spacedim> *>(&fe.base_element(0)) != nullptr;

  const auto lexicographic_to_hierarchic_numbering =
    Utilities::invert_permutation(
      FETools::hierarchic_to_lexicographic_numbering<spacedim>(
        this->get_degree()));

  // Step 1: copy global vector so that the ghost values are such that the
  // cache can be set up for all ghost cells
  LinearAlgebra::distributed::Vector<typename VectorType::value_type>
                 vector_ghosted;
  const IndexSet locally_relevant_dofs =
    DoFTools::extract_locally_relevant_dofs(dof_handler);
  vector_ghosted.reinit(dof_handler.locally_owned_dofs(),
                        locally_relevant_dofs,
                        dof_handler.get_mpi_communicator());
  copy_locally_owned_data_from(vector, vector_ghosted);
  vector_ghosted.update_ghost_values();

  // FE and FEValues in the case they are needed
  FE_Nothing<dim, spacedim> fe_nothing;
  Threads::ThreadLocalStorage<std::unique_ptr<FEValues<dim, spacedim>>>
    fe_values_all;

  // Interpolation of values is needed if we cannot just read off locations
  // from the solution vectors (as in the case of FE_Q and FE_DGQ with the
  // same polynomial degree as this class has).
  const bool interpolation_of_values_is_needed =
    ((is_fe_q || is_fe_dgq) && fe.degree == this->get_degree()) == false;

  // Step 2: loop over all cells
  this->initialize(
    dof_handler.get_triangulation(),
    [&](const typename Triangulation<dim, spacedim>::cell_iterator &cell_tria)
      -> std::vector<Point<spacedim>> {
      const bool is_active_non_artificial_cell =
        (cell_tria->is_active() == true) &&
        (cell_tria->is_artificial() == false);

      const typename DoFHandler<dim, spacedim>::cell_iterator cell_dofs(
        &cell_tria->get_triangulation(),
        cell_tria->level(),
        cell_tria->index(),
        &dof_handler);

      const auto mapping_q =
        dynamic_cast<const MappingQ<dim, spacedim> *>(&mapping);

      // Step 2a) set up and reinit FEValues (if needed)
      if (
        ((vector_describes_relative_displacement ||
          (is_active_non_artificial_cell == false)) &&
         ((mapping_q != nullptr &&
           this->get_degree() == mapping_q->get_degree()) ==
          false)) /*condition 1: points need to be computed via FEValues*/
        ||
        (is_active_non_artificial_cell && interpolation_of_values_is_needed) /*condition 2: interpolation of values is needed*/)
        {
          // get FEValues (thread-safe); in the case that this thread has
          // not created a an FEValues object yet, this helper-function also
          // creates one with the right quadrature rule
          auto &fe_values = fe_values_all.get();
          if (fe_values.get() == nullptr)
            {
              const QGaussLobatto<dim> quadrature_gl(this->polynomial_degree +
                                                     1);

              std::vector<Point<dim>> quadrature_points;
              for (const auto i :
                   FETools::hierarchic_to_lexicographic_numbering<dim>(
                     this->polynomial_degree))
                quadrature_points.push_back(quadrature_gl.point(i));
              const Quadrature<dim> quadrature(quadrature_points);

              fe_values = std::make_unique<FEValues<dim, spacedim>>(
                mapping,
                interpolation_of_values_is_needed ?
                  fe :
                  static_cast<const FiniteElement<dim, spacedim> &>(fe_nothing),
                quadrature,
                update_quadrature_points | update_values);
            }

          if (interpolation_of_values_is_needed)
            fe_values->reinit(cell_dofs);
          else
            fe_values->reinit(cell_tria);
        }

      std::vector<Point<spacedim>> result;

      // Step 2b) read of quadrature points in the relative displacement case
      // note: we also take this path for non-active or artificial cells so that
      // these cells are filled with some useful data
      if (vector_describes_relative_displacement ||
          is_active_non_artificial_cell == false)
        {
          if (mapping_q != nullptr &&
              this->get_degree() == mapping_q->get_degree())
            result = mapping_q->compute_mapping_support_points(cell_tria);
          else
            result = fe_values_all.get()->get_quadrature_points();

          // for non-active or artificial cells we are done here and return
          // the absolute positions, since the provided vector cannot contain
          // any useful information for these cells
          if (is_active_non_artificial_cell == false)
            return result;
        }
      else
        {
          result.resize(
            Utilities::pow<unsigned int>(this->get_degree() + 1, dim));
        }

      // Step 2c) read global vector and adjust points accordingly
      if (interpolation_of_values_is_needed == false)
        {
          // case 1: FE_Q or FE_DGQ with same degree as this class has; this
          // is the simple case since no interpolation is needed
          std::vector<types::global_dof_index> dof_indices(
            fe.n_dofs_per_cell());
          cell_dofs->get_dof_indices(dof_indices);

          for (unsigned int i = 0; i < dof_indices.size(); ++i)
            {
              const auto id = fe.system_to_component_index(i);

              if (is_fe_q)
                {
                  // case 1a: FE_Q
                  if (vector_describes_relative_displacement)
                    result[id.second][id.first] +=
                      vector_ghosted(dof_indices[i]);
                  else
                    result[id.second][id.first] =
                      vector_ghosted(dof_indices[i]);
                }
              else
                {
                  // case 1b: FE_DGQ
                  if (vector_describes_relative_displacement)
                    result[lexicographic_to_hierarchic_numbering[id.second]]
                          [id.first] += vector_ghosted(dof_indices[i]);
                  else
                    result[lexicographic_to_hierarchic_numbering[id.second]]
                          [id.first] = vector_ghosted(dof_indices[i]);
                }
            }
        }
      else
        {
          // case 2: general case; interpolation is needed
          // note: the following code could be optimized for tensor-product
          // elements via application of sum factorization as is done on
          // MatrixFree/FEEvaluation
          auto &fe_values = fe_values_all.get();

          std::vector<Vector<typename VectorType::value_type>> values(
            fe_values->n_quadrature_points,
            Vector<typename VectorType::value_type>(spacedim));

          fe_values->get_function_values(vector_ghosted, values);

          for (unsigned int q = 0; q < fe_values->n_quadrature_points; ++q)
            for (unsigned int c = 0; c < spacedim; ++c)
              if (vector_describes_relative_displacement)
                result[q][c] += values[q][c];
              else
                result[q][c] = values[q][c];
        }

      return result;
    });

  uses_level_info = false;
}



template <int dim, int spacedim>
template <typename VectorType>
void
MappingQCache<dim, spacedim>::initialize(
  const Mapping<dim, spacedim>    &mapping,
  const DoFHandler<dim, spacedim> &dof_handler,
  const MGLevelObject<VectorType> &vectors,
  const bool                       vector_describes_relative_displacement)
{
  AssertDimension(dof_handler.get_fe_collection().size(), 1);
  const FiniteElement<dim, spacedim> &fe = dof_handler.get_fe();
  AssertDimension(fe.n_base_elements(), 1);
  AssertDimension(fe.element_multiplicity(0), spacedim);
  AssertDimension(0, vectors.min_level());
  AssertDimension(dof_handler.get_triangulation().n_global_levels() - 1,
                  vectors.max_level());

  const unsigned int is_fe_q =
    dynamic_cast<const FE_Q<dim, spacedim> *>(&fe.base_element(0)) != nullptr;
  const unsigned int is_fe_dgq =
    dynamic_cast<const FE_DGQ<dim, spacedim> *>(&fe.base_element(0)) != nullptr;

  const auto lexicographic_to_hierarchic_numbering =
    Utilities::invert_permutation(
      FETools::hierarchic_to_lexicographic_numbering<spacedim>(
        this->get_degree()));

  // Step 1: copy global vector so that the ghost values are such that the
  // cache can be set up for all ghost cells
  MGLevelObject<
    LinearAlgebra::distributed::Vector<typename VectorType::value_type>>
    vectors_ghosted(vectors.min_level(), vectors.max_level());

  for (unsigned int l = vectors.min_level(); l <= vectors.max_level(); ++l)
    {
      const IndexSet locally_relevant_dofs =
        DoFTools::extract_locally_relevant_level_dofs(dof_handler, l);
      vectors_ghosted[l].reinit(dof_handler.locally_owned_mg_dofs(l),
                                locally_relevant_dofs,
                                dof_handler.get_mpi_communicator());
      copy_locally_owned_data_from(vectors[l], vectors_ghosted[l]);
      vectors_ghosted[l].update_ghost_values();
    }

  // FE and FEValues in the case they are needed
  FE_Nothing<dim, spacedim> fe_nothing;
  Threads::ThreadLocalStorage<std::unique_ptr<FEValues<dim, spacedim>>>
    fe_values_all;

  // Interpolation of values is needed if we cannot just read off locations
  // from the solution vectors (as in the case of FE_Q and FE_DGQ with the
  // same polynomial degree as this class has).
  const bool interpolation_of_values_is_needed =
    ((is_fe_q || is_fe_dgq) && fe.degree == this->get_degree()) == false;

  // Step 2: loop over all cells
  this->initialize(
    dof_handler.get_triangulation(),
    [&](const typename Triangulation<dim, spacedim>::cell_iterator &cell_tria)
      -> std::vector<Point<spacedim>> {
      const bool is_non_artificial_cell =
        cell_tria->level_subdomain_id() != numbers::artificial_subdomain_id;

      const typename DoFHandler<dim, spacedim>::level_cell_iterator cell_dofs(
        &cell_tria->get_triangulation(),
        cell_tria->level(),
        cell_tria->index(),
        &dof_handler);

      const auto mapping_q =
        dynamic_cast<const MappingQ<dim, spacedim> *>(&mapping);

      // Step 2a) set up and reinit FEValues (if needed)
      if (
        ((vector_describes_relative_displacement ||
          (is_non_artificial_cell == false)) &&
         ((mapping_q != nullptr &&
           this->get_degree() == mapping_q->get_degree()) ==
          false)) /*condition 1: points need to be computed via FEValues*/
        ||
        (is_non_artificial_cell == true && interpolation_of_values_is_needed) /*condition 2: interpolation of values is needed*/)
        {
          // get FEValues (thread-safe); in the case that this thread has
          // not created a an FEValues object yet, this helper-function also
          // creates one with the right quadrature rule
          auto &fe_values = fe_values_all.get();
          if (fe_values.get() == nullptr)
            {
              const QGaussLobatto<dim> quadrature_gl(this->polynomial_degree +
                                                     1);

              std::vector<Point<dim>> quadrature_points;
              for (const auto i :
                   FETools::hierarchic_to_lexicographic_numbering<dim>(
                     this->polynomial_degree))
                quadrature_points.push_back(quadrature_gl.point(i));
              const Quadrature<dim> quadrature(quadrature_points);

              fe_values = std::make_unique<FEValues<dim, spacedim>>(
                mapping,
                interpolation_of_values_is_needed ?
                  fe :
                  static_cast<const FiniteElement<dim, spacedim> &>(fe_nothing),
                quadrature,
                update_quadrature_points | update_values);
            }

          if (interpolation_of_values_is_needed)
            fe_values->reinit(cell_dofs);
          else
            fe_values->reinit(cell_tria);
        }

      std::vector<Point<spacedim>> result;

      // Step 2b) read of quadrature points in the relative displacement case
      // note: we also take this path for non-active or artificial cells so that
      // these cells are filled with some useful data
      if (vector_describes_relative_displacement ||
          (is_non_artificial_cell == false))
        {
          if (mapping_q != nullptr &&
              this->get_degree() == mapping_q->get_degree())
            result = mapping_q->compute_mapping_support_points(cell_tria);
          else
            result = fe_values_all.get()->get_quadrature_points();

          // for non-active or artificial cells we are done here and return
          // the absolute positions, since the provided vector cannot contain
          // any useful information for these cells
          if (is_non_artificial_cell == false)
            return result;
        }
      else
        {
          result.resize(
            Utilities::pow<unsigned int>(this->get_degree() + 1, dim));
        }

      // Step 2c) read global vector and adjust points accordingly
      if (interpolation_of_values_is_needed == false)
        {
          // case 1: FE_Q or FE_DGQ with same degree as this class has; this
          // is the simple case since no interpolation is needed
          std::vector<types::global_dof_index> dof_indices(
            fe.n_dofs_per_cell());
          cell_dofs->get_mg_dof_indices(dof_indices);

          for (unsigned int i = 0; i < dof_indices.size(); ++i)
            {
              const auto id = fe.system_to_component_index(i);

              if (is_fe_q)
                {
                  // case 1a: FE_Q
                  if (vector_describes_relative_displacement)
                    result[id.second][id.first] +=
                      vectors_ghosted[cell_tria->level()](dof_indices[i]);
                  else
                    result[id.second][id.first] =
                      vectors_ghosted[cell_tria->level()](dof_indices[i]);
                }
              else
                {
                  // case 1b: FE_DGQ
                  if (vector_describes_relative_displacement)
                    result[lexicographic_to_hierarchic_numbering[id.second]]
                          [id.first] +=
                      vectors_ghosted[cell_tria->level()](dof_indices[i]);
                  else
                    result[lexicographic_to_hierarchic_numbering[id.second]]
                          [id.first] =
                            vectors_ghosted[cell_tria->level()](dof_indices[i]);
                }
            }
        }
      else
        {
          // case 2: general case; interpolation is needed
          // note: the following code could be optimized for tensor-product
          // elements via application of sum factorization as is done on
          // MatrixFree/FEEvaluation
          auto &fe_values = fe_values_all.get();

          std::vector<types::global_dof_index> dof_indices(
            fe.n_dofs_per_cell());
          cell_dofs->get_mg_dof_indices(dof_indices);

          std::vector<typename VectorType::value_type> dof_values(
            fe.n_dofs_per_cell());

          for (unsigned int i = 0; i < fe.n_dofs_per_cell(); ++i)
            dof_values[i] = vectors_ghosted[cell_tria->level()](dof_indices[i]);

          for (unsigned int c = 0; c < spacedim; ++c)
            for (unsigned int i = 0; i < fe.n_dofs_per_cell(); ++i)
              for (unsigned int q = 0; q < fe_values->n_quadrature_points; ++q)
                if (vector_describes_relative_displacement == false && i == 0)
                  result[q][c] =
                    dof_values[i] * fe_values->shape_value_component(i, q, c);
                else
                  result[q][c] +=
                    dof_values[i] * fe_values->shape_value_component(i, q, c);
        }

      return result;
    });

  uses_level_info = true;
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

  Assert(uses_level_info || cell->is_active(), ExcInternalError());

  AssertIndexRange(cell->level(), support_point_cache->size());
  AssertIndexRange(cell->index(), (*support_point_cache)[cell->level()].size());
  return (*support_point_cache)[cell->level()][cell->index()];
}



template <int dim, int spacedim>
boost::container::small_vector<Point<spacedim>,
#ifndef _MSC_VER
                               ReferenceCells::max_n_vertices<dim>()
#else
                               GeometryInfo<dim>::vertices_per_cell
#endif
                               >
MappingQCache<dim, spacedim>::get_vertices(
  const typename Triangulation<dim, spacedim>::cell_iterator &cell) const
{
  Assert(support_point_cache.get() != nullptr,
         ExcMessage("Must call MappingQCache::initialize() before "
                    "using it or after mesh has changed!"));

  Assert(uses_level_info || cell->is_active(), ExcInternalError());

  AssertIndexRange(cell->level(), support_point_cache->size());
  AssertIndexRange(cell->index(), (*support_point_cache)[cell->level()].size());
  const auto ptr = (*support_point_cache)[cell->level()][cell->index()].begin();
  return boost::container::small_vector<Point<spacedim>,
#ifndef _MSC_VER
                                        ReferenceCells::max_n_vertices<dim>()
#else
                                        GeometryInfo<dim>::vertices_per_cell
#endif
                                        >(ptr, ptr + cell->n_vertices());
}



//--------------------------- Explicit instantiations -----------------------
#include "fe/mapping_q_cache.inst"


DEAL_II_NAMESPACE_CLOSE
