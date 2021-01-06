// ---------------------------------------------------------------------
//
// Copyright (C) 2020 by the deal.II authors
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


#ifndef dealii_vector_repartitioner_templates_h
#define dealii_vector_repartitioner_templates_h

#include <deal.II/multigrid/mg_transfer_global_coarsening.templates.h>

#include <deal.II/numerics/vector_repartitioner.h>


DEAL_II_NAMESPACE_OPEN

namespace
{
  template <typename MeshType>
  class PermutationFineDoFHandlerView
    : public internal::FineDoFHandlerView<MeshType>
  {
  public:
    PermutationFineDoFHandlerView(const MeshType &dof_handler_dst,
                                  const MeshType &dof_handler_src)
      : internal::FineDoFHandlerView<MeshType>(dof_handler_dst, dof_handler_src)
    {
      // get reference to triangulations
      const auto &tria_dst = dof_handler_dst.get_triangulation();
      const auto &tria_src = dof_handler_src.get_triangulation();

      // create index sets
      IndexSet is_dst_locally_owned(this->cell_id_translator.size());
      IndexSet is_dst_remote(this->cell_id_translator.size());
      IndexSet is_src_locally_owned(this->cell_id_translator.size());

      for (auto cell : tria_dst.active_cell_iterators())
        if (!cell->is_artificial() && cell->is_locally_owned())
          is_dst_locally_owned.add_index(
            this->cell_id_translator.translate(cell));


      for (auto cell : tria_src.active_cell_iterators())
        if (!cell->is_artificial() && cell->is_locally_owned())
          {
            is_src_locally_owned.add_index(
              this->cell_id_translator.translate(cell));
            is_dst_remote.add_index(this->cell_id_translator.translate(cell));
          }

      this->reinit(is_dst_locally_owned,
                   is_dst_remote,
                   is_src_locally_owned,
                   false);
    }
  };
} // namespace



template <typename Number>
void
VectorRepartitioner<LinearAlgebra::distributed::Vector<Number>>::
  update_forwards(LinearAlgebra::distributed::Vector<Number> &      dst,
                  const LinearAlgebra::distributed::Vector<Number> &src) const
{
  // create new source vector with matching ghost values
  LinearAlgebra::distributed::Vector<Number> src_extended(extended_partitioner);

  // copy locally owned values from original source vector
  src_extended.copy_locally_owned_data_from(src);

  // update ghost values
  src_extended.update_ghost_values();

  // copy locally owned values from temporal array to destination vector
  for (unsigned int i = 0; i < indices.size(); ++i)
    dst.local_element(i) = src_extended.local_element(indices[i]);
}



template <typename Number>
void
VectorRepartitioner<LinearAlgebra::distributed::Vector<Number>>::
  update_backwards(LinearAlgebra::distributed::Vector<Number> &      dst,
                   const LinearAlgebra::distributed::Vector<Number> &src) const
{
  // create new source vector with matching ghost values
  LinearAlgebra::distributed::Vector<Number> dst_extended(extended_partitioner);

  // copy locally owned values from temporal array to destination vector
  for (unsigned int i = 0; i < indices.size(); ++i)
    dst_extended.local_element(indices[i]) = src.local_element(i);

  // update ghost values
  dst_extended.compress(VectorOperation::values::add);

  // copy locally owned values from original source vector
  dst.copy_locally_owned_data_from(dst_extended);
}



template <typename Number>
template <int dim, int spacedim>
void
VectorRepartitioner<LinearAlgebra::distributed::Vector<Number>>::reinit(
  const DoFHandler<dim, spacedim> &dof_handler_dst,
  const DoFHandler<dim, spacedim> &dof_handler_src)
{
  const PermutationFineDoFHandlerView<DoFHandler<dim, spacedim>> view(
    dof_handler_src, dof_handler_dst);

  // 1) setup intermediate vector
  this->extended_partitioner = std::make_shared<Utilities::MPI::Partitioner>(
    view.locally_owned_dofs(),
    view.locally_relevant_dofs(),
    internal::get_mpi_comm(
      dof_handler_src) /*TODO: generalize for different comms*/);

  // 2) setup indices for fast copying between vector
  this->indices.resize(dof_handler_dst.locally_owned_dofs().n_elements());

  std::vector<types::global_dof_index> indices_src;
  std::vector<types::global_dof_index> indices_dst;

  // loop over all active locally owned cells
  for (const auto &cell_dst : dof_handler_dst.active_cell_iterators())
    {
      if (cell_dst->is_artificial() || cell_dst->is_locally_owned() == false)
        continue;

      const unsigned int dofs_per_cell =
        dof_handler_src.get_fe_collection()[cell_dst->active_fe_index()]
          .n_dofs_per_cell();

      indices_src.resize(dofs_per_cell);
      indices_dst.resize(dofs_per_cell);

      // get view on source cell (possibly remote)
      const auto cell_src = view.get_cell(cell_dst);

      // get indices
      cell_dst->get_dof_indices(indices_dst);
      cell_src.get_dof_indices(indices_src);

      // use src/dst indices to set up indices
      for (unsigned int j = 0; j < dofs_per_cell; ++j)
        if (dof_handler_dst.locally_owned_dofs().is_element(indices_dst[j]))
          this->indices[dof_handler_dst.locally_owned_dofs().index_within_set(
            indices_dst[j])] =
            this->extended_partitioner->global_to_local(indices_src[j]);
    }
}


DEAL_II_NAMESPACE_CLOSE

#endif
