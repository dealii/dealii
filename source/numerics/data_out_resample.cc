// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2021 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q_iso_q1.h>
#include <deal.II/fe/fe_tools.h>

#include <deal.II/lac/la_parallel_block_vector.h>
#include <deal.II/lac/la_parallel_vector.h>

#include <deal.II/numerics/data_out_dof_data.templates.h>
#include <deal.II/numerics/data_out_resample.h>
#include <deal.II/numerics/vector_tools.h>

#include <sstream>

DEAL_II_NAMESPACE_OPEN


template <int dim, int patch_dim, int spacedim>
DataOutResample<dim, patch_dim, spacedim>::DataOutResample(
  const Triangulation<patch_dim, spacedim> &patch_tria,
  const Mapping<patch_dim, spacedim>       &patch_mapping)
  : patch_dof_handler(patch_tria)
  , patch_mapping(&patch_mapping)
{}



template <int dim, int patch_dim, int spacedim>
void
DataOutResample<dim, patch_dim, spacedim>::update_mapping(
  const Mapping<dim, spacedim> &mapping,
  const unsigned int            n_subdivisions)
{
  this->mapping = &mapping;
  this->point_to_local_vector_indices.clear();

  FE_Q_iso_Q1<patch_dim, spacedim> fe(
    std::max<unsigned int>(1, n_subdivisions));
  patch_dof_handler.distribute_dofs(fe);

  std::vector<std::pair<types::global_dof_index, Point<spacedim>>> points_all;

  const Quadrature<patch_dim> quadrature(fe.get_unit_support_points());

  FEValues<patch_dim, spacedim> fe_values(*patch_mapping,
                                          fe,
                                          quadrature,
                                          update_quadrature_points);

  std::vector<types::global_dof_index> dof_indices(fe.n_dofs_per_cell());

  const IndexSet active_dofs =
    DoFTools::extract_locally_active_dofs(patch_dof_handler);
  partitioner = std::make_shared<Utilities::MPI::Partitioner>(
    patch_dof_handler.locally_owned_dofs(), active_dofs, MPI_COMM_WORLD);

  for (const auto &cell : patch_dof_handler.active_cell_iterators() |
                            IteratorFilters::LocallyOwnedCell())
    {
      fe_values.reinit(cell);

      cell->get_dof_indices(dof_indices);

      for (unsigned int i = 0; i < dof_indices.size(); ++i)
        points_all.emplace_back(partitioner->global_to_local(dof_indices[i]),
                                fe_values.quadrature_point(i));
    }

  std::sort(points_all.begin(),
            points_all.end(),
            [](const auto &a, const auto &b) { return a.first < b.first; });
  points_all.erase(std::unique(points_all.begin(),
                               points_all.end(),
                               [](const auto &a, const auto &b) {
                                 return a.first == b.first;
                               }),
                   points_all.end());

  std::vector<Point<spacedim>> points;

  for (const auto &i : points_all)
    {
      point_to_local_vector_indices.push_back(i.first);
      points.push_back(i.second);
    }

  rpe.reinit(points, *this->triangulation, *this->mapping);
}



template <int dim, int patch_dim, int spacedim>
void
DataOutResample<dim, patch_dim, spacedim>::build_patches(
  const Mapping<dim, spacedim>                                 &mapping,
  const unsigned int                                            n_subdivisions,
  const typename DataOut<patch_dim, spacedim>::CurvedCellRegion curved_region)
{
  this->update_mapping(mapping, n_subdivisions);
  this->build_patches(curved_region);
}



template <int dim, int patch_dim, int spacedim>
void
DataOutResample<dim, patch_dim, spacedim>::build_patches(
  const typename DataOut<patch_dim, spacedim>::CurvedCellRegion curved_region)
{
  patch_data_out.clear();

  if (rpe.is_ready() == false)
    {
      Assert(
        this->mapping,
        ExcMessage(
          "Mapping is not valid anymore! Please register a new mapping via "
          "update_mapping() or the other build_patches() function."));
      update_mapping(*this->mapping, patch_dof_handler.get_fe().degree);
    }

  std::vector<std::shared_ptr<LinearAlgebra::distributed::Vector<double>>>
    vectors;

  patch_data_out.attach_dof_handler(patch_dof_handler);

  unsigned int counter = 0;

  for (const auto &data : this->dof_data)
    {
      const auto data_ptr = dynamic_cast<
        internal::DataOutImplementation::DataEntry<dim, spacedim, double> *>(
        data.get());

      Assert(data_ptr, ExcNotImplemented());

      const auto &dh = *data_ptr->dof_handler;

      if constexpr (running_in_debug_mode())
        {
          for (const auto &fe : dh.get_fe_collection())
            Assert(
              fe.n_base_elements() == 1,
              ExcMessage(
                "This class currently only supports scalar elements and elements "
                "with a single base element."));
        }

      for (unsigned int comp = 0; comp < dh.get_fe_collection().n_components();
           ++comp)
        {
          const auto values = VectorTools::point_values<1>(
            rpe, dh, data_ptr->vector, VectorTools::EvaluationFlags::avg, comp);

          vectors.emplace_back(
            std::make_shared<LinearAlgebra::distributed::Vector<double>>(
              partitioner));

          for (unsigned int j = 0; j < values.size(); ++j)
            vectors.back()->local_element(point_to_local_vector_indices[j]) =
              values[j];

          vectors.back()->set_ghost_state(true);

          // we can give the vectors arbitrary names ("temp_*") here, since
          // these are only used internally (by patch_data_out) but not later on
          // during the actual output to file
          patch_data_out.add_data_vector(
            *vectors.back(),
            std::string("temp_" + std::to_string(counter)),
            DataOut_DoFData<patch_dim, patch_dim, spacedim, spacedim>::
              DataVectorType::type_dof_data);

          ++counter;
        }
    }

  patch_data_out.build_patches(*patch_mapping,
                               patch_dof_handler.get_fe().degree,
                               curved_region);
}



template <int dim, int patch_dim, int spacedim>
const std::vector<typename DataOutBase::Patch<patch_dim, spacedim>> &
DataOutResample<dim, patch_dim, spacedim>::get_patches() const
{
  return patch_data_out.get_patches();
}


// explicit instantiations
#include "numerics/data_out_resample.inst"

DEAL_II_NAMESPACE_CLOSE
