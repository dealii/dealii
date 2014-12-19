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
#include <deal.II/lac/block_indices.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/block_sparse_matrix.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe.h>
#include <deal.II/multigrid/mg_dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/multigrid/mg_transfer_component.h>
#include <deal.II/multigrid/mg_transfer_component.templates.h>
#include <deal.II/multigrid/mg_tools.h>

#include <algorithm>
#include <numeric>
#include <iostream>

DEAL_II_NAMESPACE_OPEN


namespace
{
  /**
   * Adjust block-vectors on all
   * levels to correct size.  Count
   * the numbers of degrees of
   * freedom on each level
   * component-wise. Then, assign
   * each block of @p vector the
   * corresponding size.
   *
   * The boolean field @p selected
   * allows restricting this
   * operation to certain
   * components. In this case, @p
   * vector will only have as many
   * blocks as there are true
   * values in @p selected (no
   * blocks of length zero are
   * padded in). If this argument
   * is omitted, all blocks will be
   * considered.
   *
   * Degrees of freedom must be
   * sorted by component in order
   * to obtain reasonable results
   * from this function.
   *
   * The argument
   * @p target_component allows to
   * re-sort and group components
   * as in
   * DoFRenumbering::component_wise.
   *
   *
   */
  template <int dim, typename number, int spacedim>
  void
  reinit_vector_by_components (
    const dealii::DoFHandler<dim,spacedim> &mg_dof,
    MGLevelObject<BlockVector<number> > &v,
    const std::vector<bool> &sel,
    const std::vector<unsigned int> &target_comp,
    std::vector<std::vector<types::global_dof_index> > &ndofs)
  {
    std::vector<bool> selected=sel;
    std::vector<unsigned int> target_component=target_comp;
    const unsigned int ncomp = mg_dof.get_fe().n_components();

    // If the selected and
    // target_component have size 0,
    // they must be replaced by default
    // values.
    //
    // Since we already made copies
    // directly after this function was
    // called, we use the arguments
    // directly.
    if (target_component.size() == 0)
      {
        target_component.resize(ncomp);
        for (unsigned int i=0; i<ncomp; ++i)
          target_component[i] = i;
      }

    // If selected is an empty vector,
    // all components are selected.
    if (selected.size() == 0)
      {
        selected.resize(target_component.size());
        std::fill_n (selected.begin(), ncomp, false);
        for (unsigned int i=0; i<target_component.size(); ++i)
          selected[target_component[i]] = true;
      }

    Assert (selected.size() == target_component.size(),
            ExcDimensionMismatch(selected.size(), target_component.size()));

    // Compute the number of blocks needed
    const unsigned int n_selected
      = std::accumulate(selected.begin(),
                        selected.end(),
                        0U);

    if (ndofs.size() == 0)
      {
        std::vector<std::vector<types::global_dof_index> >
        new_dofs(mg_dof.get_tria().n_levels(),
                 std::vector<types::global_dof_index>(target_component.size()));
        std::swap(ndofs, new_dofs);
        MGTools::count_dofs_per_block (mg_dof, ndofs, target_component);
      }

    for (unsigned int level=v.min_level();
         level<=v.max_level(); ++level)
      {
        v[level].reinit(n_selected, 0);
        unsigned int k=0;
        for (unsigned int i=0; i<selected.size() && (k<v[level].n_blocks()); ++i)
          {
            if (selected[i])
              {
                v[level].block(k++).reinit(ndofs[level][i]);
              }
            v[level].collect_sizes();
          }
      }
  }


  /**
   * Adjust vectors on all levels
   * to correct size.  Count the
   * numbers of degrees of freedom
   * on each level component-wise
   * in a single component. Then,
   * assign @p vector the
   * corresponding size.
   *
   * The boolean field @p selected
   * may be nonzero in a single
   * component, indicating the
   * block of a block vector the
   * argument @p v corresponds to.
   *
   * Degrees of freedom must be
   * sorted by component in order
   * to obtain reasonable results
   * from this function.
   *
   * The argument
   * @p target_component allows to
   * re-sort and group components
   * as in
   * DoFRenumbering::component_wise.
   */
  template <int dim, typename number, int spacedim>
  void
  reinit_vector_by_components (
    const dealii::DoFHandler<dim,spacedim> &mg_dof,
    MGLevelObject<dealii::Vector<number> > &v,
    const ComponentMask &component_mask,
    const std::vector<unsigned int> &target_component,
    std::vector<std::vector<types::global_dof_index> > &ndofs)
  {
    Assert (component_mask.represents_n_components(target_component.size()),
            ExcMessage ("The component mask does not have the correct size."));

    unsigned int selected_block = 0;
    for (unsigned int i=0; i<target_component.size(); ++i)
      if (component_mask[i])
        selected_block = target_component[i];

    if (ndofs.size() == 0)
      {
        std::vector<std::vector<types::global_dof_index> >
        new_dofs(mg_dof.get_tria().n_levels(),
                 std::vector<types::global_dof_index>(target_component.size()));
        std::swap(ndofs, new_dofs);
        MGTools::count_dofs_per_block (mg_dof, ndofs,
                                       target_component);
      }

    for (unsigned int level=v.min_level();
         level<=v.max_level(); ++level)
      {
        v[level].reinit(ndofs[level][selected_block]);
      }
  }
}


template <typename number>
template <int dim, class InVector, int spacedim>
void
MGTransferSelect<number>::do_copy_to_mg (
  const DoFHandler<dim,spacedim>        &mg_dof_handler,
  MGLevelObject<Vector<number> > &dst,
  const InVector                 &src) const
{
  dst=0;

  Assert(sizes.size()==mg_dof_handler.get_tria().n_levels(),
         ExcMatricesNotBuilt());

  reinit_vector_by_components(mg_dof_handler, dst,
                              mg_component_mask,
                              mg_target_component, sizes);

  // traverse the grid top-down
  // (i.e. starting with the most
  // refined grid). this way, we can
  // always get that part of one
  // level of the output vector which
  // corresponds to a region which is
  // more refined, by restriction of
  // the respective vector on the
  // next finer level, which we then
  // already have built.

  bool first = true;
  for (unsigned int level=mg_dof_handler.get_tria().n_levels(); level!=0;)
    {
      --level;

      typedef std::vector<std::pair<types::global_dof_index, unsigned int> >::const_iterator IT;
      for (IT i=copy_to_and_from_indices[level].begin();
           i != copy_to_and_from_indices[level].end(); ++i)
        dst[level](i->second) = src(i->first);
      // for that part of the level
      // which is further refined:
      // get the defect by
      // restriction of the defect on
      // one level higher
      if (!first)
        restrict_and_add (level+1, dst[level], dst[level+1]);
      first = false;
    }
}


template <int dim, int spacedim>
void MGTransferComponentBase::build_matrices (
  const DoFHandler<dim,spacedim> &,
  const DoFHandler<dim,spacedim> &mg_dof)
{
  // Fill target component with
  // standard values (identity) if it
  // is empty
  if (target_component.size() == 0)
    {
      target_component.resize(mg_dof.get_fe().n_components());
      for (unsigned int i=0; i<target_component.size(); ++i)
        target_component[i] = i;
    }
  else
    {
      // otherwise, check it for consistency
      Assert (target_component.size() == mg_dof.get_fe().n_components(),
              ExcDimensionMismatch(target_component.size(),
                                   mg_dof.get_fe().n_components()));

      for (unsigned int i=0; i<target_component.size(); ++i)
        {
          Assert(i<target_component.size(),
                 ExcIndexRange(i,0,target_component.size()));
        }
    }
  // Do the same for the multilevel
  // components. These may be
  // different.
  if (mg_target_component.size() == 0)
    {
      mg_target_component.resize(mg_dof.get_fe().n_components());
      for (unsigned int i=0; i<mg_target_component.size(); ++i)
        mg_target_component[i] = target_component[i];
    }
  else
    {
      Assert (mg_target_component.size() == mg_dof.get_fe().n_components(),
              ExcDimensionMismatch(mg_target_component.size(),
                                   mg_dof.get_fe().n_components()));

      for (unsigned int i=0; i<mg_target_component.size(); ++i)
        {
          Assert(i<mg_target_component.size(),
                 ExcIndexRange(i,0,mg_target_component.size()));
        }
    }

  const FiniteElement<dim> &fe = mg_dof.get_fe();

  // Effective number of components
  // is the maximum entry in
  // mg_target_component. This
  // assumes that the values in that
  // vector don't have holes.
  const unsigned int n_components  =
    *std::max_element(mg_target_component.begin(), mg_target_component.end()) + 1;
  const unsigned int dofs_per_cell = fe.dofs_per_cell;
  const unsigned int n_levels      = mg_dof.get_tria().n_levels();

  Assert (mg_component_mask.represents_n_components(fe.n_components()),
          ExcMessage ("Component mask has wrong size."));

  // Compute the lengths of all blocks
  sizes.resize(n_levels);
  for (unsigned int l=0; l<n_levels; ++l)
    sizes[l].resize(n_components);

  MGTools::count_dofs_per_block(mg_dof, sizes, mg_target_component);

  // Fill some index vectors
  // for later use.
  mg_component_start = sizes;
  // Compute start indices from sizes
  for (unsigned int l=0; l<mg_component_start.size(); ++l)
    {
      types::global_dof_index k=0;
      for (unsigned int i=0; i<mg_component_start[l].size(); ++i)
        {
          const types::global_dof_index t=mg_component_start[l][i];
          mg_component_start[l][i] = k;
          k += t;
        }
    }

  component_start.resize(*std::max_element (target_component.begin(),
                                            target_component.end()) + 1);
  DoFTools::
  count_dofs_per_block (mg_dof, component_start, target_component);

  types::global_dof_index k=0;
  for (unsigned int i=0; i<component_start.size(); ++i)
    {
      const types::global_dof_index t=component_start[i];
      component_start[i] = k;
      k += t;
    }

  // Build index vectors for
  // copy_to_mg and
  // copy_from_mg. These vectors must
  // be prebuilt, since the
  // get_dof_indices functions are
  // too slow

  copy_to_and_from_indices.resize(n_levels);

// Building the prolongation matrices starts here!

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
      prolongation_sparsities
      .push_back (std_cxx11::shared_ptr<BlockSparsityPattern> (new BlockSparsityPattern));
      prolongation_matrices
      .push_back (std_cxx11::shared_ptr<BlockSparseMatrix<double> > (new BlockSparseMatrix<double>));
    }

  // two fields which will store the
  // indices of the multigrid dofs
  // for a cell and one of its children
  std::vector<types::global_dof_index> dof_indices_parent (dofs_per_cell);
  std::vector<types::global_dof_index> dof_indices_child (dofs_per_cell);

  // for each level: first build the
  // sparsity pattern of the matrices
  // and then build the matrices
  // themselves. note that we only
  // need to take care of cells on
  // the coarser level which have
  // children
  for (unsigned int level=0; level<n_levels-1; ++level)
    {
      // reset the dimension of the
      // structure.  note that for
      // the number of entries per
      // row, the number of parent
      // dofs coupling to a child dof
      // is necessary. this, is the
      // number of degrees of freedom
      // per cell
      prolongation_sparsities[level]->reinit (n_components, n_components);
      for (unsigned int i=0; i<n_components; ++i)
        for (unsigned int j=0; j<n_components; ++j)
          if (i==j)
            prolongation_sparsities[level]->block(i,j)
            .reinit(sizes[level+1][i],
                    sizes[level][j],
                    dofs_per_cell+1);
          else
            prolongation_sparsities[level]->block(i,j)
            .reinit(sizes[level+1][i],
                    sizes[level][j],
                    0);

      prolongation_sparsities[level]->collect_sizes();

      for (typename DoFHandler<dim,spacedim>::cell_iterator cell=mg_dof.begin(level);
           cell != mg_dof.end(level); ++cell)
        if (cell->has_children())
          {
            cell->get_mg_dof_indices (dof_indices_parent);

            for (unsigned int child=0; child<cell->n_children(); ++child)
              {
                // set an alias to the
                // prolongation matrix for
                // this child
                const FullMatrix<double> &prolongation
                  = mg_dof.get_fe().get_prolongation_matrix (child, cell->refinement_case());

                cell->child(child)->get_mg_dof_indices (dof_indices_child);

                // now tag the entries in the
                // matrix which will be used
                // for this pair of parent/child
                for (unsigned int i=0; i<dofs_per_cell; ++i)
                  for (unsigned int j=0; j<dofs_per_cell; ++j)
                    if (prolongation(i,j) != 0)
                      {
                        const unsigned int icomp
                          = fe.system_to_component_index(i).first;
                        const unsigned int jcomp
                          = fe.system_to_component_index(j).first;
                        if ((icomp==jcomp) && mg_component_mask[icomp])
                          prolongation_sparsities[level]->add(dof_indices_child[i],
                                                              dof_indices_parent[j]);
                      };
              };
          };
      prolongation_sparsities[level]->compress ();

      prolongation_matrices[level]->reinit (*prolongation_sparsities[level]);
      // now actually build the matrices
      for (typename DoFHandler<dim,spacedim>::cell_iterator cell=mg_dof.begin(level);
           cell != mg_dof.end(level); ++cell)
        if (cell->has_children())
          {
            cell->get_mg_dof_indices (dof_indices_parent);

            for (unsigned int child=0; child<cell->n_children(); ++child)
              {
                // set an alias to the
                // prolongation matrix for
                // this child
                const FullMatrix<double> &prolongation
                  = mg_dof.get_fe().get_prolongation_matrix (child, cell->refinement_case());

                cell->child(child)->get_mg_dof_indices (dof_indices_child);

                // now set the entries in the
                // matrix
                for (unsigned int i=0; i<dofs_per_cell; ++i)
                  for (unsigned int j=0; j<dofs_per_cell; ++j)
                    if (prolongation(i,j) != 0)
                      {
                        const unsigned int icomp = fe.system_to_component_index(i).first;
                        const unsigned int jcomp = fe.system_to_component_index(j).first;
                        if ((icomp==jcomp) && mg_component_mask[icomp])
                          prolongation_matrices[level]->set(dof_indices_child[i],
                                                            dof_indices_parent[j],
                                                            prolongation(i,j));
                      }
              }
          }
    }
  // impose boundary conditions
  // but only in the column of
  // the prolongation matrix
  //TODO: this way is not very efficient

  if (boundary_indices.size() != 0)
    {
      std::vector<std::vector<types::global_dof_index> >
      dofs_per_component(mg_dof.get_tria().n_levels(),
                         std::vector<types::global_dof_index>(n_components));

      MGTools::count_dofs_per_block (mg_dof, dofs_per_component, mg_target_component);
      for (unsigned int level=0; level<n_levels-1; ++level)
        {
          if (boundary_indices[level].size() == 0)
            continue;

          for (unsigned int iblock=0; iblock<n_components; ++iblock)
            for (unsigned int jblock=0; jblock<n_components; ++jblock)
              if (iblock==jblock)
                {
                  const types::global_dof_index n_dofs = prolongation_matrices[level]->block(iblock,jblock).m();
                  for (types::global_dof_index i=0; i<n_dofs; ++i)
                    {
                      SparseMatrix<double>::iterator anfang = prolongation_matrices[level]->block(iblock,jblock).begin(i),
                                                     ende = prolongation_matrices[level]->block(iblock,jblock).end(i);
                      for (; anfang != ende; ++anfang)
                        {
                          const types::global_dof_index column_number = anfang->column();

                          //convert global indices into local ones
                          const BlockIndices block_indices_coarse (dofs_per_component[level]);
                          const types::global_dof_index global_j = block_indices_coarse.local_to_global(iblock, column_number);

                          std::set<types::global_dof_index>::const_iterator found_dof =
                            boundary_indices[level].find(global_j);

                          const bool is_boundary_index =
                            (found_dof != boundary_indices[level].end());

                          if (is_boundary_index)
                            {
                              prolongation_matrices[level]->block(iblock,jblock)
                              .set(i,column_number,0);
                            }
                        }
                    }
                }
        }
    }
}


template <typename number>
template <int dim, int spacedim>
void MGTransferSelect<number>::build_matrices (
  const DoFHandler<dim,spacedim> &dof,
  const DoFHandler<dim,spacedim> &mg_dof,
  unsigned int select,
  unsigned int mg_select,
  const std::vector<unsigned int> &t_component,
  const std::vector<unsigned int> &mg_t_component,
  const std::vector<std::set<types::global_dof_index> > &bdry_indices)
{
  const FiniteElement<dim> &fe = mg_dof.get_fe();
  unsigned int ncomp = mg_dof.get_fe().n_components();

  target_component = t_component;
  mg_target_component = mg_t_component;
  boundary_indices = bdry_indices;

  selected_component = select;
  mg_selected_component = mg_select;

  {
    std::vector<bool> tmp(ncomp, false);
    for (unsigned int c=0; c<ncomp; ++c)
      if (t_component[c] == selected_component)
        tmp[c] = true;
    component_mask = ComponentMask(tmp);
  }

  {
    std::vector<bool> tmp(ncomp, false);
    for (unsigned int c=0; c<ncomp; ++c)
      if (mg_t_component[c] == mg_selected_component)
        tmp[c] = true;
    mg_component_mask = ComponentMask(tmp);
  }

  // If components are renumbered,
  // find the first original
  // component corresponding to the
  // target component.
  for (unsigned int i=0; i<target_component.size(); ++i)
    {
      if (target_component[i] == select)
        {
          selected_component = i;
          break;
        }
    }

  for (unsigned int i=0; i<mg_target_component.size(); ++i)
    {
      if (mg_target_component[i] == mg_select)
        {
          mg_selected_component = i;
          break;
        }
    }

  MGTransferComponentBase::build_matrices (dof, mg_dof);

  interface_dofs.resize(mg_dof.get_tria().n_levels());
  for (unsigned int l=0; l<mg_dof.get_tria().n_levels(); ++l)
    interface_dofs[l].resize(mg_dof.n_dofs(l));
  MGTools::extract_inner_interface_dofs(mg_dof, interface_dofs);

  // use a temporary vector to create the
  // relation between global and level dofs
  std::vector<types::global_dof_index> temp_copy_indices;
  std::vector<types::global_dof_index> global_dof_indices (fe.dofs_per_cell);
  std::vector<types::global_dof_index> level_dof_indices  (fe.dofs_per_cell);
  for (int level=dof.get_tria().n_levels()-1; level>=0; --level)
    {
      copy_to_and_from_indices[level].clear();
      typename DoFHandler<dim,spacedim>::active_cell_iterator
      level_cell = mg_dof.begin_active(level);
      const typename DoFHandler<dim,spacedim>::active_cell_iterator
      level_end  = mg_dof.end_active(level);

      temp_copy_indices.resize (0);
      temp_copy_indices.resize (mg_dof.n_dofs(level), numbers::invalid_dof_index);

      // Compute coarse level right hand side
      // by restricting from fine level.
      for (; level_cell!=level_end; ++level_cell)
        {
          // get the dof numbers of
          // this cell for the global
          // and the level-wise
          // numbering
          level_cell->get_dof_indices(global_dof_indices);
          level_cell->get_mg_dof_indices (level_dof_indices);

          for (unsigned int i=0; i<fe.dofs_per_cell; ++i)
            {
              const unsigned int component
                = fe.system_to_component_index(i).first;
              if (component_mask[component] &&
                  !interface_dofs[level][level_dof_indices[i]])
                {
                  const types::global_dof_index level_start
                    = mg_component_start[level][mg_target_component[component]];
                  const types::global_dof_index global_start
                    = component_start[target_component[component]];
                  temp_copy_indices[level_dof_indices[i]-level_start] =
                    global_dof_indices[i] - global_start;
                }
            }
        }

      // write indices from vector into the map from
      // global to level dofs
      const types::global_dof_index n_active_dofs =
        std::count_if (temp_copy_indices.begin(), temp_copy_indices.end(),
                       std::bind2nd(std::not_equal_to<types::global_dof_index>(),
                                    numbers::invalid_dof_index));
      copy_to_and_from_indices[level].resize (n_active_dofs);
      types::global_dof_index counter = 0;
      for (types::global_dof_index i=0; i<temp_copy_indices.size(); ++i)
        if (temp_copy_indices[i] != numbers::invalid_dof_index)
          copy_to_and_from_indices[level][counter++] =
            std::pair<types::global_dof_index, unsigned int> (temp_copy_indices[i], i);
      Assert (counter == n_active_dofs, ExcInternalError());
    }
}



// explicit instantiations
#include "mg_transfer_component.inst"


DEAL_II_NAMESPACE_CLOSE
