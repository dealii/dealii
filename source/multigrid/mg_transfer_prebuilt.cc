// ---------------------------------------------------------------------
//
// Copyright (C) 2003 - 2013 by the deal.II authors
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
#include <deal.II/lac/compressed_simple_sparsity_pattern.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe.h>
#include <deal.II/multigrid/mg_dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/multigrid/mg_tools.h>
#include <deal.II/multigrid/mg_transfer.h>
#include <deal.II/multigrid/mg_transfer.templates.h>

#include <algorithm>

DEAL_II_NAMESPACE_OPEN


template<class VECTOR>
MGTransferPrebuilt<VECTOR>::MGTransferPrebuilt ()
{}


template<class VECTOR>
MGTransferPrebuilt<VECTOR>::MGTransferPrebuilt (const ConstraintMatrix &c, const MGConstrainedDoFs &mg_c)
  :
  constraints(&c),
  mg_constrained_dofs(&mg_c)
{}


template <class VECTOR>
MGTransferPrebuilt<VECTOR>::~MGTransferPrebuilt ()
{}


template <class VECTOR>
void MGTransferPrebuilt<VECTOR>::initialize_constraints (
  const ConstraintMatrix &c, const MGConstrainedDoFs &mg_c)
{
  constraints = &c;
  mg_constrained_dofs = &mg_c;
}


template <class VECTOR>
void MGTransferPrebuilt<VECTOR>::clear ()
{
  sizes.resize(0);
  prolongation_matrices.resize(0);
  prolongation_sparsities.resize(0);
  copy_indices.resize(0);
  copy_indices_to_me.resize(0);
  copy_indices_from_me.resize(0);
  component_to_block_map.resize(0);
  interface_dofs.resize(0);
  constraints = 0;
  mg_constrained_dofs = 0;
}


template <class VECTOR>
void MGTransferPrebuilt<VECTOR>::prolongate (
  const unsigned int to_level,
  VECTOR            &dst,
  const VECTOR      &src) const
{
  Assert ((to_level >= 1) && (to_level<=prolongation_matrices.size()),
          ExcIndexRange (to_level, 1, prolongation_matrices.size()+1));

  prolongation_matrices[to_level-1]->vmult (dst, src);
}


template <class VECTOR>
void MGTransferPrebuilt<VECTOR>::restrict_and_add (
  const unsigned int   from_level,
  VECTOR       &dst,
  const VECTOR &src) const
{
  Assert ((from_level >= 1) && (from_level<=prolongation_matrices.size()),
          ExcIndexRange (from_level, 1, prolongation_matrices.size()+1));

  prolongation_matrices[from_level-1]->Tvmult_add (dst, src);
}


template <typename VECTOR>
template <int dim, int spacedim>
void MGTransferPrebuilt<VECTOR>::build_matrices (
  const DoFHandler<dim,spacedim>  &mg_dof)
{
  const unsigned int n_levels      = mg_dof.get_tria().n_global_levels();
  const unsigned int dofs_per_cell = mg_dof.get_fe().dofs_per_cell;

  sizes.resize(n_levels);
  for (unsigned int l=0; l<n_levels; ++l)
    sizes[l] = mg_dof.n_dofs(l);

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
      (std_cxx11::shared_ptr<typename internal::MatrixSelector<VECTOR>::Sparsity> (new typename internal::MatrixSelector<VECTOR>::Sparsity));
      prolongation_matrices.push_back
      (std_cxx11::shared_ptr<typename internal::MatrixSelector<VECTOR>::Matrix> (new typename internal::MatrixSelector<VECTOR>::Matrix));
    }

  // two fields which will store the
  // indices of the multigrid dofs
  // for a cell and one of its children
  std::vector<types::global_dof_index> dof_indices_parent (dofs_per_cell);
  std::vector<types::global_dof_index> dof_indices_child (dofs_per_cell);

  // for each level: first build the sparsity
  // pattern of the matrices and then build the
  // matrices themselves. note that we only
  // need to take care of cells on the coarser
  // level which have children
  for (unsigned int level=0; level<n_levels-1; ++level)
    {

      // reset the dimension of the structure.
      // note that for the number of entries
      // per row, the number of parent dofs
      // coupling to a child dof is
      // necessary. this, of course, is the
      // number of degrees of freedom per
      // cell
      // increment dofs_per_cell
      // since a useless diagonal
      // element will be stored
      CompressedSimpleSparsityPattern csp (sizes[level+1],
                                           sizes[level]);
      std::vector<types::global_dof_index> entries (dofs_per_cell);
      for (typename DoFHandler<dim,spacedim>::cell_iterator cell=mg_dof.begin(level);
           cell != mg_dof.end(level); ++cell)
        if (cell->has_children() &&
            ( mg_dof.get_tria().locally_owned_subdomain()==numbers::invalid_subdomain_id
              || cell->level_subdomain_id()==mg_dof.get_tria().locally_owned_subdomain()
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
                    csp.add_entries (dof_indices_child[i],
                                     entries.begin(), entries.end());
                  }
              }
          }

      internal::MatrixSelector<VECTOR>::reinit(*prolongation_matrices[level],
                                               *prolongation_sparsities[level],
                                               level,
                                               csp,
                                               mg_dof);
      csp.reinit(0,0);

      FullMatrix<double> prolongation;

      // now actually build the matrices
      for (typename DoFHandler<dim,spacedim>::cell_iterator cell=mg_dof.begin(level);
           cell != mg_dof.end(level); ++cell)
        if (cell->has_children() &&
            (mg_dof.get_tria().locally_owned_subdomain()==numbers::invalid_subdomain_id
             || cell->level_subdomain_id()==mg_dof.get_tria().locally_owned_subdomain())
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

                if (mg_constrained_dofs != 0 && mg_constrained_dofs->set_boundary_values())
                  for (unsigned int j=0; j<dofs_per_cell; ++j)
                    if (mg_constrained_dofs->is_boundary_index(level, dof_indices_parent[j]))
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

  // Now we are filling the variables copy_indices*, which are essentially
  // maps from global to mgdof for each level stored as a std::vector of
  // pairs. We need to split this map on each level depending on the ownership
  // of the global and mgdof, so that we later not access non-local elements
  // in copy_to/from_mg.
  // We keep track in the bitfield dof_touched which global dof has
  // been processed already (on the current level). This is the same as
  // the multigrid running in serial.
  // Only entering on the finest level gives wrong results (why?)

  copy_indices.resize(n_levels);
  copy_indices_from_me.resize(n_levels);
  copy_indices_to_me.resize(n_levels);
  IndexSet globally_relevant;
  DoFTools::extract_locally_relevant_dofs(mg_dof, globally_relevant);

  std::vector<types::global_dof_index> global_dof_indices (dofs_per_cell);
  std::vector<types::global_dof_index> level_dof_indices  (dofs_per_cell);
  //  for (int level=mg_dof.get_tria().n_levels()-1; level>=0; --level)
  for (unsigned int level=0; level<mg_dof.get_tria().n_levels(); ++level)
    {
      std::vector<bool> dof_touched(globally_relevant.n_elements(), false);
      copy_indices[level].clear();
      copy_indices_from_me[level].clear();
      copy_indices_to_me[level].clear();

      typename DoFHandler<dim,spacedim>::active_cell_iterator
      level_cell = mg_dof.begin_active(level);
      const typename DoFHandler<dim,spacedim>::active_cell_iterator
      level_end  = mg_dof.end_active(level);

      for (; level_cell!=level_end; ++level_cell)
        {
          if (mg_dof.get_tria().locally_owned_subdomain()!=numbers::invalid_subdomain_id
              &&  (level_cell->level_subdomain_id()==numbers::artificial_subdomain_id
                   ||  level_cell->subdomain_id()==numbers::artificial_subdomain_id)
             )
            continue;

          // get the dof numbers of this cell for the global and the level-wise
          // numbering
          level_cell->get_dof_indices (global_dof_indices);
          level_cell->get_mg_dof_indices (level_dof_indices);

          for (unsigned int i=0; i<dofs_per_cell; ++i)
            {
              // we need to ignore if the DoF is on a refinement edge (hanging node)
              if (mg_constrained_dofs != 0
                  && mg_constrained_dofs->at_refinement_edge(level, level_dof_indices[i]))
                continue;
              unsigned int global_idx = globally_relevant.index_within_set(global_dof_indices[i]);
              //skip if we did this global dof already (on this or a coarser level)
              if (dof_touched[global_idx])
                continue;
              bool global_mine = mg_dof.locally_owned_dofs().is_element(global_dof_indices[i]);
              bool level_mine = mg_dof.locally_owned_mg_dofs(level).is_element(level_dof_indices[i]);

              if (global_mine && level_mine)
                copy_indices[level].push_back(
                  std::pair<unsigned int, unsigned int> (global_dof_indices[i], level_dof_indices[i]));
              else if (level_mine)
                copy_indices_from_me[level].push_back(
                  std::pair<unsigned int, unsigned int> (global_dof_indices[i], level_dof_indices[i]));
              else if (global_mine)
                copy_indices_to_me[level].push_back(
                  std::pair<unsigned int, unsigned int> (global_dof_indices[i], level_dof_indices[i]));
              else
                continue;

              dof_touched[global_idx] = true;
            }
        }
    }

  // If we are in debugging mode, we order the copy indices, so we get
  // more reliable output for regression texts
#ifdef DEBUG
  std::less<std::pair<types::global_dof_index, unsigned int> > compare;
  for (unsigned int level=0; level<copy_indices.size(); ++level)
    std::sort(copy_indices[level].begin(), copy_indices[level].end(), compare);
  for (unsigned int level=0; level<copy_indices_from_me.size(); ++level)
    std::sort(copy_indices_from_me[level].begin(), copy_indices_from_me[level].end(), compare);
  for (unsigned int level=0; level<copy_indices_to_me.size(); ++level)
    std::sort(copy_indices_to_me[level].begin(), copy_indices_to_me[level].end(), compare);
#endif
}


template <class VECTOR>
void
MGTransferPrebuilt<VECTOR>::print_matrices (std::ostream &os) const
{
  for (unsigned int level = 0; level<prolongation_matrices.size(); ++level)
    {
      os << "Level " << level << std::endl;
      prolongation_matrices[level]->print(os);
      os << std::endl;
    }
}

template <class VECTOR>
void
MGTransferPrebuilt<VECTOR>::print_indices (std::ostream &os) const
{
  for (unsigned int level = 0; level<copy_indices.size(); ++level)
    {
      for (unsigned int i=0; i<copy_indices[level].size(); ++i)
        os << "copy_indices[" << level
           << "]\t" << copy_indices[level][i].first << '\t' << copy_indices[level][i].second << std::endl;
    }

  for (unsigned int level = 0; level<copy_indices_from_me.size(); ++level)
    {
      for (unsigned int i=0; i<copy_indices_from_me[level].size(); ++i)
        os << "copy_ifrom  [" << level
           << "]\t" << copy_indices_from_me[level][i].first << '\t' << copy_indices_from_me[level][i].second << std::endl;
    }
  for (unsigned int level = 0; level<copy_indices_to_me.size(); ++level)
    {
      for (unsigned int i=0; i<copy_indices_to_me[level].size(); ++i)
        os << "copy_ito    [" << level
           << "]\t" << copy_indices_to_me[level][i].first << '\t' << copy_indices_to_me[level][i].second << std::endl;
    }
}


// explicit instantiation
#include "mg_transfer_prebuilt.inst"


DEAL_II_NAMESPACE_CLOSE
