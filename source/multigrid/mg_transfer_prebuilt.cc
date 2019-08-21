// ---------------------------------------------------------------------
//
// Copyright (C) 2003 - 2019 by the deal.II authors
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


#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/block_sparse_matrix.h>
#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/la_parallel_block_vector.h>
#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/sparsity_tools.h>
#include <deal.II/lac/trilinos_epetra_vector.h>
#include <deal.II/lac/trilinos_parallel_block_vector.h>
#include <deal.II/lac/trilinos_vector.h>
#include <deal.II/lac/vector.h>

#include <deal.II/multigrid/mg_tools.h>
#include <deal.II/multigrid/mg_transfer.h>

#include <algorithm>

DEAL_II_NAMESPACE_OPEN


template <typename VectorType>
MGTransferPrebuilt<VectorType>::MGTransferPrebuilt(
  const MGConstrainedDoFs &mg_c)
{
  this->mg_constrained_dofs = &mg_c;
}



template <typename VectorType>
MGTransferPrebuilt<VectorType>::MGTransferPrebuilt(
  const AffineConstraints<double> & /*c*/,
  const MGConstrainedDoFs &mg_c)
{
  this->mg_constrained_dofs = &mg_c;
}



template <typename VectorType>
void
MGTransferPrebuilt<VectorType>::initialize_constraints(
  const MGConstrainedDoFs &mg_c)
{
  this->mg_constrained_dofs = &mg_c;
}



template <typename VectorType>
void
MGTransferPrebuilt<VectorType>::initialize_constraints(
  const AffineConstraints<double> & /*c*/,
  const MGConstrainedDoFs &mg_c)
{
  initialize_constraints(mg_c);
}



template <typename VectorType>
void
MGTransferPrebuilt<VectorType>::clear()
{
  MGLevelGlobalTransfer<VectorType>::clear();
  prolongation_matrices.resize(0);
  prolongation_sparsities.resize(0);
  interface_dofs.resize(0);
}



template <typename VectorType>
void
MGTransferPrebuilt<VectorType>::prolongate(const unsigned int to_level,
                                           VectorType &       dst,
                                           const VectorType & src) const
{
  Assert((to_level >= 1) && (to_level <= prolongation_matrices.size()),
         ExcIndexRange(to_level, 1, prolongation_matrices.size() + 1));

#ifdef DEBUG
  if (this->mg_constrained_dofs != nullptr)
    Assert(this->mg_constrained_dofs->get_user_constraint_matrix(to_level - 1)
               .get_local_lines()
               .size() == 0,
           ExcNotImplemented());
#endif

  prolongation_matrices[to_level - 1]->vmult(dst, src);
}



template <typename VectorType>
void
MGTransferPrebuilt<VectorType>::restrict_and_add(const unsigned int from_level,
                                                 VectorType &       dst,
                                                 const VectorType & src) const
{
  Assert((from_level >= 1) && (from_level <= prolongation_matrices.size()),
         ExcIndexRange(from_level, 1, prolongation_matrices.size() + 1));
  (void)from_level;

  prolongation_matrices[from_level - 1]->Tvmult_add(dst, src);
}


namespace
{
  /**
   * Helper function for build_matrices. Checks for identity constrained dofs
   * and replace with the indices of the dofs to which they are constrained
   */
  void
  replace(const MGConstrainedDoFs *             mg_constrained_dofs,
          const unsigned int                    level,
          std::vector<types::global_dof_index> &dof_indices)
  {
    if (mg_constrained_dofs != nullptr &&
        mg_constrained_dofs->get_level_constraints(level).n_constraints() > 0)
      for (auto &ind : dof_indices)
        if (mg_constrained_dofs->get_level_constraints(level)
              .is_identity_constrained(ind))
          {
            Assert(mg_constrained_dofs->get_level_constraints(level)
                       .get_constraint_entries(ind)
                       ->size() == 1,
                   ExcInternalError());
            ind = mg_constrained_dofs->get_level_constraints(level)
                    .get_constraint_entries(ind)
                    ->front()
                    .first;
          }
  }
} // namespace

template <typename VectorType>
template <int dim, int spacedim>
void
MGTransferPrebuilt<VectorType>::build_matrices(
  const DoFHandler<dim, spacedim> &mg_dof)
{
  const unsigned int n_levels = mg_dof.get_triangulation().n_global_levels();
  const unsigned int dofs_per_cell = mg_dof.get_fe().dofs_per_cell;

  this->sizes.resize(n_levels);
  for (unsigned int l = 0; l < n_levels; ++l)
    this->sizes[l] = mg_dof.n_dofs(l);

  // reset the size of the array of
  // matrices. call resize(0) first,
  // in order to delete all elements
  // and clear their memory. then
  // repopulate these arrays
  //
  // note that on resize(0), the
  // shared_ptr class takes care of
  // deleting the object it points to
  // by itself
  prolongation_matrices.resize(0);
  prolongation_sparsities.resize(0);
  prolongation_matrices.reserve(n_levels - 1);
  prolongation_sparsities.reserve(n_levels - 1);

  for (unsigned int i = 0; i < n_levels - 1; ++i)
    {
      prolongation_sparsities.emplace_back(
        new typename internal::MatrixSelector<VectorType>::Sparsity);
      prolongation_matrices.emplace_back(
        new typename internal::MatrixSelector<VectorType>::Matrix);
    }

  // two fields which will store the
  // indices of the multigrid dofs
  // for a cell and one of its children
  std::vector<types::global_dof_index> dof_indices_parent(dofs_per_cell);
  std::vector<types::global_dof_index> dof_indices_child(dofs_per_cell);
  std::vector<types::global_dof_index> entries(dofs_per_cell);

  // for each level: first build the sparsity
  // pattern of the matrices and then build the
  // matrices themselves. note that we only
  // need to take care of cells on the coarser
  // level which have children
  for (unsigned int level = 0; level < n_levels - 1; ++level)
    {
      // reset the dimension of the structure.  note that for the number of
      // entries per row, the number of parent dofs coupling to a child dof is
      // necessary. this, of course, is the number of degrees of freedom per
      // cell
      //
      // increment dofs_per_cell since a useless diagonal element will be
      // stored
      IndexSet level_p1_relevant_dofs;
      DoFTools::extract_locally_relevant_level_dofs(mg_dof,
                                                    level + 1,
                                                    level_p1_relevant_dofs);
      DynamicSparsityPattern                  dsp(this->sizes[level + 1],
                                 this->sizes[level],
                                 level_p1_relevant_dofs);
      typename DoFHandler<dim>::cell_iterator cell, endc = mg_dof.end(level);
      for (cell = mg_dof.begin(level); cell != endc; ++cell)
        if (cell->has_children() &&
            (mg_dof.get_triangulation().locally_owned_subdomain() ==
               numbers::invalid_subdomain_id ||
             cell->level_subdomain_id() ==
               mg_dof.get_triangulation().locally_owned_subdomain()))
          {
            cell->get_mg_dof_indices(dof_indices_parent);

            replace(this->mg_constrained_dofs, level, dof_indices_parent);

            Assert(cell->n_children() ==
                     GeometryInfo<dim>::max_children_per_cell,
                   ExcNotImplemented());
            for (unsigned int child = 0; child < cell->n_children(); ++child)
              {
                // set an alias to the prolongation matrix for this child
                const FullMatrix<double> &prolongation =
                  mg_dof.get_fe().get_prolongation_matrix(
                    child, cell->refinement_case());

                Assert(prolongation.n() != 0, ExcNoProlongation());

                cell->child(child)->get_mg_dof_indices(dof_indices_child);

                replace(this->mg_constrained_dofs,
                        level + 1,
                        dof_indices_child);

                // now tag the entries in the
                // matrix which will be used
                // for this pair of parent/child
                for (unsigned int i = 0; i < dofs_per_cell; ++i)
                  {
                    entries.resize(0);
                    for (unsigned int j = 0; j < dofs_per_cell; ++j)
                      if (prolongation(i, j) != 0)
                        entries.push_back(dof_indices_parent[j]);
                    dsp.add_entries(dof_indices_child[i],
                                    entries.begin(),
                                    entries.end());
                  }
              }
          }

#ifdef DEAL_II_WITH_MPI
      if (internal::MatrixSelector<
            VectorType>::requires_distributed_sparsity_pattern)
        {
          // Since PETSc matrices do not offer the functionality to fill up in-
          // complete sparsity patterns on their own, the sparsity pattern must
          // be manually distributed.

          // Retrieve communicator from triangulation if it is parallel
          const parallel::TriangulationBase<dim, spacedim> *dist_tria =
            dynamic_cast<const parallel::TriangulationBase<dim, spacedim> *>(
              &(mg_dof.get_triangulation()));

          MPI_Comm communicator = dist_tria != nullptr ?
                                    dist_tria->get_communicator() :
                                    MPI_COMM_SELF;

          // Compute # of locally owned MG dofs / processor for distribution
          const std::vector<::dealii::IndexSet>
            locally_owned_mg_dofs_per_processor =
              mg_dof.compute_locally_owned_mg_dofs_per_processor(level + 1);
          std::vector<::dealii::types::global_dof_index>
            n_locally_owned_mg_dofs_per_processor(
              locally_owned_mg_dofs_per_processor.size(), 0);

          for (std::size_t index = 0;
               index < n_locally_owned_mg_dofs_per_processor.size();
               ++index)
            {
              n_locally_owned_mg_dofs_per_processor[index] =
                locally_owned_mg_dofs_per_processor[index].n_elements();
            }

          // Distribute sparsity pattern
          ::dealii::SparsityTools::distribute_sparsity_pattern(
            dsp,
            n_locally_owned_mg_dofs_per_processor,
            communicator,
            dsp.row_index_set());
        }
#endif

      internal::MatrixSelector<VectorType>::reinit(
        *prolongation_matrices[level],
        *prolongation_sparsities[level],
        level,
        dsp,
        mg_dof);
      dsp.reinit(0, 0);

      // In the end, the entries in this object will only be real valued.
      // Nevertheless, we have to take the underlying scalar type of the
      // vector we want to use this class with. The global matrix the entries
      // of this matrix are copied into has to be able to perform a
      // matrix-vector multiplication and this is in general only implemented if
      // the scalar type for matrix and vector is the same. Furthermore,
      // copying entries between this local object and the global matrix is only
      // implemented if the objects have the same scalar type.
      FullMatrix<typename VectorType::value_type> prolongation;

      // now actually build the matrices
      for (cell = mg_dof.begin(level); cell != endc; ++cell)
        if (cell->has_children() &&
            (mg_dof.get_triangulation().locally_owned_subdomain() ==
               numbers::invalid_subdomain_id ||
             cell->level_subdomain_id() ==
               mg_dof.get_triangulation().locally_owned_subdomain()))
          {
            cell->get_mg_dof_indices(dof_indices_parent);

            replace(this->mg_constrained_dofs, level, dof_indices_parent);

            Assert(cell->n_children() ==
                     GeometryInfo<dim>::max_children_per_cell,
                   ExcNotImplemented());
            for (unsigned int child = 0; child < cell->n_children(); ++child)
              {
                // set an alias to the prolongation matrix for this child
                prolongation = mg_dof.get_fe().get_prolongation_matrix(
                  child, cell->refinement_case());

                if (this->mg_constrained_dofs != nullptr &&
                    this->mg_constrained_dofs->have_boundary_indices())
                  for (unsigned int j = 0; j < dofs_per_cell; ++j)
                    if (this->mg_constrained_dofs->is_boundary_index(
                          level, dof_indices_parent[j]))
                      for (unsigned int i = 0; i < dofs_per_cell; ++i)
                        prolongation(i, j) = 0.;

                cell->child(child)->get_mg_dof_indices(dof_indices_child);

                replace(this->mg_constrained_dofs,
                        level + 1,
                        dof_indices_child);

                // now set the entries in the matrix
                for (unsigned int i = 0; i < dofs_per_cell; ++i)
                  prolongation_matrices[level]->set(dof_indices_child[i],
                                                    dofs_per_cell,
                                                    dof_indices_parent.data(),
                                                    &prolongation(i, 0),
                                                    true);
              }
          }
      prolongation_matrices[level]->compress(VectorOperation::insert);
    }

  this->fill_and_communicate_copy_indices(mg_dof);
}



template <typename VectorType>
void
MGTransferPrebuilt<VectorType>::print_matrices(std::ostream &os) const
{
  for (unsigned int level = 0; level < prolongation_matrices.size(); ++level)
    {
      os << "Level " << level << std::endl;
      prolongation_matrices[level]->print(os);
      os << std::endl;
    }
}



template <typename VectorType>
std::size_t
MGTransferPrebuilt<VectorType>::memory_consumption() const
{
  std::size_t result = MGLevelGlobalTransfer<VectorType>::memory_consumption();
  for (unsigned int i = 0; i < prolongation_matrices.size(); ++i)
    result += prolongation_matrices[i]->memory_consumption() +
              prolongation_sparsities[i]->memory_consumption();

  return result;
}


// explicit instantiation
#include "mg_transfer_prebuilt.inst"


DEAL_II_NAMESPACE_CLOSE
