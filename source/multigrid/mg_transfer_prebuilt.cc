// ---------------------------------------------------------------------
//
// Copyright (C) 2003 - 2016 by the deal.II authors
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


#include <deal.II/base/logstream.h>
#include <deal.II/base/function.h>

#include <deal.II/lac/vector.h>
#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/parallel_vector.h>
#include <deal.II/lac/parallel_block_vector.h>
#include <deal.II/lac/petsc_vector.h>
#include <deal.II/lac/petsc_block_vector.h>
#include <deal.II/lac/trilinos_vector.h>
#include <deal.II/lac/trilinos_block_vector.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/block_sparse_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/multigrid/mg_tools.h>
#include <deal.II/multigrid/mg_transfer.h>

#include <algorithm>

DEAL_II_NAMESPACE_OPEN


template<typename VectorType>
MGTransferPrebuilt<VectorType>::MGTransferPrebuilt ()
{}



template<typename VectorType>
MGTransferPrebuilt<VectorType>::MGTransferPrebuilt (const ConstraintMatrix &c, const MGConstrainedDoFs &mg_c)
  :
  constraints(&c)
{
  this->mg_constrained_dofs = &mg_c;
}



template <typename VectorType>
MGTransferPrebuilt<VectorType>::~MGTransferPrebuilt ()
{}



template <typename VectorType>
void MGTransferPrebuilt<VectorType>::initialize_constraints
(const ConstraintMatrix &c, const MGConstrainedDoFs &mg_c)
{
  constraints = &c;
  this->mg_constrained_dofs = &mg_c;
}



template <typename VectorType>
void MGTransferPrebuilt<VectorType>::clear ()
{
  MGLevelGlobalTransfer<VectorType>::clear();
  prolongation_matrices.resize(0);
  prolongation_sparsities.resize(0);
  interface_dofs.resize(0);
  constraints = 0;
}



template <typename VectorType>
void MGTransferPrebuilt<VectorType>::prolongate (const unsigned int to_level,
                                                 VectorType        &dst,
                                                 const VectorType  &src) const
{
  Assert ((to_level >= 1) && (to_level<=prolongation_matrices.size()),
          ExcIndexRange (to_level, 1, prolongation_matrices.size()+1));

  prolongation_matrices[to_level-1]->vmult (dst, src);
}



template <typename VectorType>
void MGTransferPrebuilt<VectorType>::restrict_and_add (const unsigned int from_level,
                                                       VectorType        &dst,
                                                       const VectorType  &src) const
{
  Assert ((from_level >= 1) && (from_level<=prolongation_matrices.size()),
          ExcIndexRange (from_level, 1, prolongation_matrices.size()+1));
  (void)from_level;

  prolongation_matrices[from_level-1]->Tvmult_add (dst, src);
}



template <typename VectorType>
template <int dim, int spacedim>
void MGTransferPrebuilt<VectorType>::build_matrices
(const DoFHandler<dim,spacedim>  &mg_dof)
{
  const unsigned int n_levels      = mg_dof.get_triangulation().n_global_levels();
  const unsigned int dofs_per_cell = mg_dof.get_fe().dofs_per_cell;

  this->sizes.resize(n_levels);
  for (unsigned int l=0; l<n_levels; ++l)
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
  prolongation_matrices.resize (0);
  prolongation_sparsities.resize (0);

  for (unsigned int i=0; i<n_levels-1; ++i)
    {
      prolongation_sparsities.push_back
      (std_cxx11::shared_ptr<typename internal::MatrixSelector<VectorType>::Sparsity> (new typename internal::MatrixSelector<VectorType>::Sparsity));
      prolongation_matrices.push_back
      (std_cxx11::shared_ptr<typename internal::MatrixSelector<VectorType>::Matrix> (new typename internal::MatrixSelector<VectorType>::Matrix));
    }

  // two fields which will store the
  // indices of the multigrid dofs
  // for a cell and one of its children
  std::vector<types::global_dof_index> dof_indices_parent (dofs_per_cell);
  std::vector<types::global_dof_index> dof_indices_child (dofs_per_cell);
  std::vector<types::global_dof_index> entries (dofs_per_cell);

  // for each level: first build the sparsity
  // pattern of the matrices and then build the
  // matrices themselves. note that we only
  // need to take care of cells on the coarser
  // level which have children
  for (unsigned int level=0; level<n_levels-1; ++level)
    {
      // reset the dimension of the structure.  note that for the number of
      // entries per row, the number of parent dofs coupling to a child dof is
      // necessary. this, of course, is the number of degrees of freedom per
      // cell
      //
      // increment dofs_per_cell since a useless diagonal element will be
      // stored
      IndexSet level_p1_relevant_dofs;
      DoFTools::extract_locally_relevant_level_dofs(mg_dof, level+1,
                                                    level_p1_relevant_dofs);
      DynamicSparsityPattern dsp (this->sizes[level+1],
                                  this->sizes[level],
                                  level_p1_relevant_dofs);
      typename DoFHandler<dim>::cell_iterator cell, endc = mg_dof.end(level);
      for (cell=mg_dof.begin(level); cell != endc; ++cell)
        if (cell->has_children() &&
            ( mg_dof.get_triangulation().locally_owned_subdomain()==numbers::invalid_subdomain_id
              || cell->level_subdomain_id()==mg_dof.get_triangulation().locally_owned_subdomain()
            ))
          {
            cell->get_mg_dof_indices (dof_indices_parent);

            Assert(cell->n_children()==GeometryInfo<dim>::max_children_per_cell,
                   ExcNotImplemented());
            for (unsigned int child=0; child<cell->n_children(); ++child)
              {
                // set an alias to the prolongation matrix for this child
                const FullMatrix<double> &prolongation
                  = mg_dof.get_fe().get_prolongation_matrix (child,
                                                             cell->refinement_case());

                Assert (prolongation.n() != 0, ExcNoProlongation());

                cell->child(child)->get_mg_dof_indices (dof_indices_child);

                // now tag the entries in the
                // matrix which will be used
                // for this pair of parent/child
                for (unsigned int i=0; i<dofs_per_cell; ++i)
                  {
                    entries.resize(0);
                    for (unsigned int j=0; j<dofs_per_cell; ++j)
                      if (prolongation(i,j) != 0)
                        entries.push_back (dof_indices_parent[j]);
                    dsp.add_entries (dof_indices_child[i],
                                     entries.begin(), entries.end());
                  }
              }
          }

      internal::MatrixSelector<VectorType>::reinit(*prolongation_matrices[level],
                                                   *prolongation_sparsities[level],
                                                   level,
                                                   dsp,
                                                   mg_dof);
      dsp.reinit(0,0);

      FullMatrix<double> prolongation;

      // now actually build the matrices
      for (cell=mg_dof.begin(level); cell != endc; ++cell)
        if (cell->has_children() &&
            (mg_dof.get_triangulation().locally_owned_subdomain()==numbers::invalid_subdomain_id
             || cell->level_subdomain_id()==mg_dof.get_triangulation().locally_owned_subdomain())
           )
          {
            cell->get_mg_dof_indices (dof_indices_parent);

            Assert(cell->n_children()==GeometryInfo<dim>::max_children_per_cell,
                   ExcNotImplemented());
            for (unsigned int child=0; child<cell->n_children(); ++child)
              {
                // set an alias to the prolongation matrix for this child
                prolongation
                  = mg_dof.get_fe().get_prolongation_matrix (child,
                                                             cell->refinement_case());

                if (this->mg_constrained_dofs != 0 &&
                    this->mg_constrained_dofs->have_boundary_indices())
                  for (unsigned int j=0; j<dofs_per_cell; ++j)
                    if (this->mg_constrained_dofs->is_boundary_index(level, dof_indices_parent[j]))
                      for (unsigned int i=0; i<dofs_per_cell; ++i)
                        prolongation(i,j) = 0.;

                cell->child(child)->get_mg_dof_indices (dof_indices_child);

                // now set the entries in the matrix
                for (unsigned int i=0; i<dofs_per_cell; ++i)
                  prolongation_matrices[level]->set (dof_indices_child[i],
                                                     dofs_per_cell,
                                                     &dof_indices_parent[0],
                                                     &prolongation(i,0),
                                                     true);
              }
          }
      prolongation_matrices[level]->compress(VectorOperation::insert);
    }

  this->fill_and_communicate_copy_indices(mg_dof);
}



template <typename VectorType>
void
MGTransferPrebuilt<VectorType>::print_matrices (std::ostream &os) const
{
  for (unsigned int level = 0; level<prolongation_matrices.size(); ++level)
    {
      os << "Level " << level << std::endl;
      prolongation_matrices[level]->print(os);
      os << std::endl;
    }
}



template <typename VectorType>
std::size_t
MGTransferPrebuilt<VectorType>::memory_consumption () const
{
  std::size_t result = MGLevelGlobalTransfer<VectorType>::memory_consumption();
  for (unsigned int i=0; i<prolongation_matrices.size(); ++i)
    result += prolongation_matrices[i]->memory_consumption()
              + prolongation_sparsities[i]->memory_consumption();

  return result;
}


// explicit instantiation
#include "mg_transfer_prebuilt.inst"


DEAL_II_NAMESPACE_CLOSE
