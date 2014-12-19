// ---------------------------------------------------------------------
//
// Copyright (C) 2009 - 2013 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------

#include <deal.II/meshworker/dof_info.h>
#include <deal.II/meshworker/integration_info.h>
#include <deal.II/base/quadrature_lib.h>

DEAL_II_NAMESPACE_OPEN


namespace MeshWorker
{
  template<int dim, int sdim>
  void
  IntegrationInfo<dim,sdim>::initialize_data(
    const std_cxx11::shared_ptr<VectorDataBase<dim,sdim> > &data)
  {
    global_data = data;
    const unsigned int nqp = fevalv[0]->n_quadrature_points;

    values.resize(global_data->n_values());
//    deallog << "values: " << values.size() << " [";
    // For all selected finite
    // element functions
    for (unsigned int i=0; i<values.size(); ++i)
      {
        values[i].resize(n_components);
//      deallog << ' ' << values[i].size() << " {";
        // For all components
        for (unsigned int j=0; j<values[i].size(); ++j)
          {
            values[i][j].resize(nqp);
//          deallog << ' ' << values[i][j].size();
          }
//      deallog << " }";
      }
//    deallog << " ]" << std::endl;

    gradients.resize(global_data->n_gradients());
    // For all selected finite
    // element functions
    for (unsigned int i=0; i<gradients.size(); ++i)
      {
        gradients[i].resize(n_components);
        // For all components
        for (unsigned int j=0; j<gradients[i].size(); ++j)
          {
            gradients[i][j].resize(nqp);
          }
      }

    hessians.resize(global_data->n_hessians());
    // For all selected finite
    // element functions
    for (unsigned int i=0; i<hessians.size(); ++i)
      {
        hessians[i].resize(n_components);
        // For all components
        for (unsigned int j=0; j<hessians[i].size(); ++j)
          {
            hessians[i][j].resize(nqp);
          }
      }
  }


  template<int dim, int sdim>
  void
  IntegrationInfo<dim,sdim>::clear()
  {
    fevalv.resize(0);
  }



  template<int dim, int sdim>
  template <typename number>
  void
  IntegrationInfo<dim,sdim>::fill_local_data(const DoFInfo<dim, sdim, number> &info, bool split_fevalues)
  {
    if (split_fevalues)
      {
        unsigned int comp = 0;
        // Loop over all blocks
        for (unsigned int b=0; b<info.block_info->local().size(); ++b)
          {
            const unsigned int fe_no = info.block_info->base_element(b);
            const FEValuesBase<dim,sdim> &fe = this->fe_values(fe_no);
            const unsigned int n_comp = fe.get_fe().n_components();
            const unsigned int block_start = info.block_info->local().block_start(b);
            const unsigned int block_size = info.block_info->local().block_size(b);

            if (info.level_cell)
              this->global_data->mg_fill(values, gradients, hessians, fe, info.cell->level(), info.indices,
                                         comp, n_comp, block_start, block_size);
            else
              this->global_data->fill(values, gradients, hessians, fe, info.indices,
                                      comp, n_comp, block_start, block_size);
            comp += n_comp;
          }
      }
    else
      {
        const FEValuesBase<dim,sdim> &fe = this->fe_values(0);
        const unsigned int n_comp = fe.get_fe().n_components();
        if (info.level_cell)
          this->global_data->mg_fill(values, gradients, hessians, fe, info.cell->level(), info.indices,
                                     0, n_comp, 0, info.indices.size());
        else
          this->global_data->fill(values, gradients, hessians, fe, info.indices,
                                  0, n_comp, 0, info.indices.size());
      }
  }


  template<int dim, int sdim>
  std::size_t
  IntegrationInfo<dim,sdim>::memory_consumption () const
  {
    std::size_t mem = sizeof(*this)
                      + MemoryConsumption::memory_consumption(fevalv)
                      - sizeof (fevalv);
    for (unsigned int i=0; i<fevalv.size(); ++i)
      mem += fevalv[i]->memory_consumption();
    return mem;
  }

//----------------------------------------------------------------------//

  template<int dim, int sdim>
  IntegrationInfoBox<dim,sdim>::IntegrationInfoBox()
  {
    cell_flags = update_default;
    boundary_flags = update_default;
    face_flags = update_default;
    neighbor_flags = update_default;
  }


  template<int dim, int sdim>
  void
  IntegrationInfoBox<dim,sdim>::initialize_update_flags (bool neighbor_geometry)
  {
    cell_flags |= update_JxW_values;
    boundary_flags |= UpdateFlags(update_JxW_values | update_normal_vectors);
    face_flags |= boundary_flags;
    neighbor_flags |= neighbor_geometry
                      ? boundary_flags
                      : update_default;

    if (cell_selector.has_values() != 0) cell_flags |= update_values;
    if (cell_selector.has_gradients() != 0) cell_flags |= update_gradients;
    if (cell_selector.has_hessians() != 0) cell_flags |= update_hessians;

    if (boundary_selector.has_values() != 0) boundary_flags |= update_values;
    if (boundary_selector.has_gradients() != 0) boundary_flags |= update_gradients;
    if (boundary_selector.has_hessians() != 0) boundary_flags |= update_hessians;

    if (face_selector.has_values() != 0) face_flags |= update_values;
    if (face_selector.has_gradients() != 0) face_flags |= update_gradients;
    if (face_selector.has_hessians() != 0) face_flags |= update_hessians;

    if (face_selector.has_values() != 0) neighbor_flags |= update_values;
    if (face_selector.has_gradients() != 0) neighbor_flags |= update_gradients;
    if (face_selector.has_hessians() != 0) neighbor_flags |= update_hessians;
  }


  template <int dim, int sdim>
  void
  IntegrationInfoBox<dim,sdim>::add_update_flags(
    const UpdateFlags flags,
    bool cell,
    bool boundary,
    bool face,
    bool neighbor)
  {
    if (cell) cell_flags |= flags;
    if (boundary) boundary_flags |= flags;
    if (face) face_flags |= flags;
    if (neighbor) neighbor_flags |= flags;
  }


  template<int dim, int sdim>
  std::size_t
  IntegrationInfoBox<dim,sdim>::memory_consumption () const
  {
    std::size_t mem = sizeof(*this)
                      + MemoryConsumption::memory_consumption(cell_quadrature)
                      - sizeof (cell_quadrature)
                      + MemoryConsumption::memory_consumption(boundary_quadrature)
                      - sizeof (boundary_quadrature)
                      + MemoryConsumption::memory_consumption(face_quadrature)
                      - sizeof (face_quadrature)
                      + MemoryConsumption::memory_consumption(cell_selector)
                      -sizeof (cell_selector)
                      + MemoryConsumption::memory_consumption(boundary_selector)
                      -sizeof (boundary_selector)
                      + MemoryConsumption::memory_consumption(face_selector)
                      -sizeof (face_selector)
                      + MemoryConsumption::memory_consumption(cell)
                      - sizeof(cell)
                      + MemoryConsumption::memory_consumption(boundary)
                      - sizeof(boundary)
                      + MemoryConsumption::memory_consumption(face)
                      - sizeof(face)
                      + MemoryConsumption::memory_consumption(subface)
                      - sizeof(subface)
                      + MemoryConsumption::memory_consumption(neighbor)
                      - sizeof(neighbor);
//   if (cell_data != 0)
//     mem += MemoryConsumption::memory_consumption(*cell_data);
//   if (boundary_data != 0)
//     mem += MemoryConsumption::memory_consumption(*boundary_data);
//   if (face_data != 0)
//     mem += MemoryConsumption::memory_consumption(*face_data);

    return mem;
  }
}


DEAL_II_NAMESPACE_CLOSE

