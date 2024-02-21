// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2009 - 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_integration_info_templates_h
#define dealii_integration_info_templates_h


#include <deal.II/base/config.h>

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/meshworker/dof_info.h>
#include <deal.II/meshworker/integration_info.h>

DEAL_II_NAMESPACE_OPEN


namespace MeshWorker
{
  template <int dim, int sdim>
  void
  IntegrationInfo<dim, sdim>::initialize_data(
    const std::shared_ptr<VectorDataBase<dim, sdim>> &data)
  {
    global_data            = data;
    const unsigned int nqp = fevalv[0]->n_quadrature_points;

    values.resize(global_data->n_values());
    // For all selected finite element functions
    for (auto &function_values : values)
      {
        function_values.resize(n_components);
        // For all components
        for (auto &component_values : function_values)
          {
            component_values.resize(nqp);
          }
      }

    gradients.resize(global_data->n_gradients());
    // For all selected finite element functions
    for (auto &function_gradients : gradients)
      {
        function_gradients.resize(n_components);
        // For all components
        for (auto &component_gradients : function_gradients)
          {
            component_gradients.resize(nqp);
          }
      }

    hessians.resize(global_data->n_hessians());
    // For all selected finite element functions
    for (auto &function_hessians : hessians)
      {
        function_hessians.resize(n_components);
        // For all components
        for (auto &component_hessians : function_hessians)
          {
            component_hessians.resize(nqp);
          }
      }
  }


  template <int dim, int sdim>
  void
  IntegrationInfo<dim, sdim>::clear()
  {
    fevalv.resize(0);
  }



  template <int dim, int sdim>
  template <typename number>
  void
  IntegrationInfo<dim, sdim>::fill_local_data(
    const DoFInfo<dim, sdim, number> &info,
    bool                              split_fevalues)
  {
    if (split_fevalues)
      {
        unsigned int comp = 0;
        // Loop over all blocks
        for (unsigned int b = 0; b < info.block_info->local().size(); ++b)
          {
            const unsigned int fe_no = info.block_info->base_element(b);
            const FEValuesBase<dim, sdim> &fe     = this->fe_values(fe_no);
            const unsigned int             n_comp = fe.get_fe().n_components();
            const unsigned int             block_start =
              info.block_info->local().block_start(b);
            const unsigned int block_size =
              info.block_info->local().block_size(b);

            if (info.level_cell)
              this->global_data->mg_fill(values,
                                         gradients,
                                         hessians,
                                         fe,
                                         info.cell->level(),
                                         info.indices,
                                         comp,
                                         n_comp,
                                         block_start,
                                         block_size);
            else
              this->global_data->fill(values,
                                      gradients,
                                      hessians,
                                      fe,
                                      info.indices,
                                      comp,
                                      n_comp,
                                      block_start,
                                      block_size);
            comp += n_comp;
          }
      }
    else
      {
        const FEValuesBase<dim, sdim> &fe     = this->fe_values(0);
        const unsigned int             n_comp = fe.get_fe().n_components();
        if (info.level_cell)
          this->global_data->mg_fill(values,
                                     gradients,
                                     hessians,
                                     fe,
                                     info.cell->level(),
                                     info.indices,
                                     0,
                                     n_comp,
                                     0,
                                     info.indices.size());
        else
          this->global_data->fill(values,
                                  gradients,
                                  hessians,
                                  fe,
                                  info.indices,
                                  0,
                                  n_comp,
                                  0,
                                  info.indices.size());
      }
  }


  template <int dim, int sdim>
  std::size_t
  IntegrationInfo<dim, sdim>::memory_consumption() const
  {
    std::size_t mem = sizeof(*this) +
                      MemoryConsumption::memory_consumption(fevalv) -
                      sizeof(fevalv);
    for (unsigned int i = 0; i < fevalv.size(); ++i)
      mem += fevalv[i]->memory_consumption();
    return mem;
  }

  //----------------------------------------------------------------------//

  template <int dim, int sdim>
  IntegrationInfoBox<dim, sdim>::IntegrationInfoBox()
  {
    cell_flags     = update_default;
    boundary_flags = update_default;
    face_flags     = update_default;
    neighbor_flags = update_default;
  }


  template <int dim, int sdim>
  void
  IntegrationInfoBox<dim, sdim>::initialize_update_flags(bool neighbor_geometry)
  {
    cell_flags |= update_JxW_values;
    boundary_flags |= UpdateFlags(update_JxW_values | update_normal_vectors);
    face_flags |= boundary_flags;
    neighbor_flags |= neighbor_geometry ? boundary_flags : update_default;

    if (cell_selector.has_values())
      cell_flags |= update_values;
    if (cell_selector.has_gradients())
      cell_flags |= update_gradients;
    if (cell_selector.has_hessians())
      cell_flags |= update_hessians;

    if (boundary_selector.has_values())
      boundary_flags |= update_values;
    if (boundary_selector.has_gradients())
      boundary_flags |= update_gradients;
    if (boundary_selector.has_hessians())
      boundary_flags |= update_hessians;

    if (face_selector.has_values())
      face_flags |= update_values;
    if (face_selector.has_gradients())
      face_flags |= update_gradients;
    if (face_selector.has_hessians())
      face_flags |= update_hessians;

    if (face_selector.has_values())
      neighbor_flags |= update_values;
    if (face_selector.has_gradients())
      neighbor_flags |= update_gradients;
    if (face_selector.has_hessians())
      neighbor_flags |= update_hessians;
  }


  template <int dim, int sdim>
  void
  IntegrationInfoBox<dim, sdim>::add_update_flags(const UpdateFlags flags,
                                                  bool              cell,
                                                  bool              boundary,
                                                  bool              face,
                                                  bool              neighbor)
  {
    if (cell)
      cell_flags |= flags;
    if (boundary)
      boundary_flags |= flags;
    if (face)
      face_flags |= flags;
    if (neighbor)
      neighbor_flags |= flags;
  }


  template <int dim, int sdim>
  std::size_t
  IntegrationInfoBox<dim, sdim>::memory_consumption() const
  {
    std::size_t mem =
      sizeof(*this) + MemoryConsumption::memory_consumption(cell_quadrature) -
      sizeof(cell_quadrature) +
      MemoryConsumption::memory_consumption(boundary_quadrature) -
      sizeof(boundary_quadrature) +
      MemoryConsumption::memory_consumption(face_quadrature) -
      sizeof(face_quadrature) +
      MemoryConsumption::memory_consumption(cell_selector) -
      sizeof(cell_selector) +
      MemoryConsumption::memory_consumption(boundary_selector) -
      sizeof(boundary_selector) +
      MemoryConsumption::memory_consumption(face_selector) -
      sizeof(face_selector) + MemoryConsumption::memory_consumption(cell) -
      sizeof(cell) + MemoryConsumption::memory_consumption(boundary) -
      sizeof(boundary) + MemoryConsumption::memory_consumption(face) -
      sizeof(face) + MemoryConsumption::memory_consumption(subface) -
      sizeof(subface) + MemoryConsumption::memory_consumption(neighbor) -
      sizeof(neighbor);
    //   if (cell_data != 0)
    //     mem += MemoryConsumption::memory_consumption(*cell_data);
    //   if (boundary_data != 0)
    //     mem += MemoryConsumption::memory_consumption(*boundary_data);
    //   if (face_data != 0)
    //     mem += MemoryConsumption::memory_consumption(*face_data);

    return mem;
  }
} // namespace MeshWorker


DEAL_II_NAMESPACE_CLOSE

#endif
