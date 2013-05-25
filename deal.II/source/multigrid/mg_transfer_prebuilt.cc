//---------------------------------------------------------------------------
//    $Id$
//
//    Copyright (C) 2003, 2004, 2005, 2006, 2007, 2008, 2009, 2010, 2012, 2013 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------

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
  const unsigned int n_levels      = mg_dof.get_tria().n_levels();
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
      (std_cxx1x::shared_ptr<typename internal::MatrixSelector<VECTOR>::Sparsity> (new typename internal::MatrixSelector<VECTOR>::Sparsity));
      prolongation_matrices.push_back
      (std_cxx1x::shared_ptr<typename internal::MatrixSelector<VECTOR>::Matrix> (new typename internal::MatrixSelector<VECTOR>::Matrix));
    }

  // two fields which will store the
  // indices of the multigrid dofs
  // for a cell and one of its children
  std::vector<unsigned int> dof_indices_parent (dofs_per_cell);
  std::vector<unsigned int> dof_indices_child (dofs_per_cell);

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
      std::vector<unsigned int> entries (dofs_per_cell);
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
                // set an alias to the
                // prolongation matrix for
                // this child
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
                // set an alias to the
                // prolongation matrix for
                // this child
                const FullMatrix<double> &prolongation
                  = mg_dof.get_fe().get_prolongation_matrix (child,
                                                             cell->refinement_case());

                cell->child(child)->get_mg_dof_indices (dof_indices_child);

                // now set the entries in the
                // matrix
                for (unsigned int i=0; i<dofs_per_cell; ++i)
                  prolongation_matrices[level]->set (dof_indices_child[i],
                                                     dofs_per_cell,
                                                     &dof_indices_parent[0],
                                                     &prolongation(i,0),
                                                     true);
              }
          }
    }


  // impose boundary conditions
  // but only in the column of
  // the prolongation matrix
  if (mg_constrained_dofs != 0)
    if (mg_constrained_dofs->set_boundary_values())
      {
        std::vector<unsigned int> constrain_indices;
        for (int level=n_levels-2; level>=0; --level)
          {
            if (mg_constrained_dofs->get_boundary_indices()[level].size() == 0)
              continue;

            // need to delete all the columns in the
            // matrix that are on the boundary. to achieve
            // this, create an array as long as there are
            // matrix columns, and find which columns we
            // need to filter away.
            constrain_indices.resize (0);
            constrain_indices.resize (prolongation_matrices[level]->n(), 0);
            std::set<unsigned int>::const_iterator dof
            = mg_constrained_dofs->get_boundary_indices()[level].begin(),
            endd = mg_constrained_dofs->get_boundary_indices()[level].end();
            for (; dof != endd; ++dof)
              constrain_indices[*dof] = 1;

            const unsigned int n_dofs = prolongation_matrices[level]->m();
            for (unsigned int i=0; i<n_dofs; ++i)
              {
                typename internal::MatrixSelector<VECTOR>::Matrix::iterator
                start_row = prolongation_matrices[level]->begin(i),
                end_row   = prolongation_matrices[level]->end(i);
                for (; start_row != end_row; ++start_row)
                  {
                    if (constrain_indices[start_row->column()] == 1)
                      start_row->value() = 0;
                  }
              }
          }
      }

  // to find the indices that describe the
  // relation between global dofs and local
  // numbering on the individual level, first
  // create a temp vector where the ith level
  // entry contains the respective global
  // entry. this gives a neat way to find those
  // indices. in a second step, actually build
  // the std::vector<std::pair<uint,uint> > that
  // only contains the active dofs on the
  // levels.

  copy_indices.resize(n_levels);
  std::vector<unsigned int> temp_copy_indices;
  std::vector<unsigned int> global_dof_indices (dofs_per_cell);
  std::vector<unsigned int> level_dof_indices  (dofs_per_cell);
  for (int level=mg_dof.get_tria().n_levels()-1; level>=0; --level)
    {
      copy_indices[level].clear();
      typename DoFHandler<dim,spacedim>::active_cell_iterator
      level_cell = mg_dof.begin_active(level);
      const typename DoFHandler<dim,spacedim>::active_cell_iterator
      level_end  = mg_dof.end_active(level);

      temp_copy_indices.resize (0);
      temp_copy_indices.resize (mg_dof.n_dofs(level), numbers::invalid_unsigned_int);

      // Compute coarse level right hand side
      // by restricting from fine level.
      for (; level_cell!=level_end; ++level_cell)
        {
          if (mg_dof.get_tria().locally_owned_subdomain()!=numbers::invalid_subdomain_id
              &&  level_cell->level_subdomain_id()!=mg_dof.get_tria().locally_owned_subdomain())
            continue;

          // get the dof numbers of
          // this cell for the global
          // and the level-wise
          // numbering
          level_cell->get_dof_indices(global_dof_indices);
          level_cell->get_mg_dof_indices (level_dof_indices);

          for (unsigned int i=0; i<dofs_per_cell; ++i)
            {
              if (mg_constrained_dofs != 0)
                {
                  if (!mg_constrained_dofs->at_refinement_edge(level,level_dof_indices[i]))
                    temp_copy_indices[level_dof_indices[i]] = global_dof_indices[i];
                }
              else
                temp_copy_indices[level_dof_indices[i]] = global_dof_indices[i];
            }
        }

      // now all the active dofs got a valid entry,
      // the other ones have an invalid entry. Count
      // the invalid entries and then resize the
      // copy_indices object. Then, insert the pairs
      // of global index and level index into
      // copy_indices.
      const unsigned int n_active_dofs =
        std::count_if (temp_copy_indices.begin(), temp_copy_indices.end(),
                       std::bind2nd(std::not_equal_to<unsigned int>(),
                                    numbers::invalid_unsigned_int));
      copy_indices[level].resize (n_active_dofs);
      unsigned int counter = 0;
      for (unsigned int i=0; i<temp_copy_indices.size(); ++i)
        if (temp_copy_indices[i] != numbers::invalid_unsigned_int)
          copy_indices[level][counter++] =
            std::pair<unsigned int, unsigned int> (temp_copy_indices[i], i);
      Assert (counter == n_active_dofs, ExcInternalError());
    }
}


template <class VECTOR>
void
MGTransferPrebuilt<VECTOR>::print_matrices (std::ostream& os) const
{
  for (unsigned int level = 0;level<prolongation_matrices.size();++level)
    {
      os << "Level " << level << std::endl;
      prolongation_matrices[level]->print(os);
      os << std::endl;
    }
}


// explicit instantiation
#include "mg_transfer_prebuilt.inst"


DEAL_II_NAMESPACE_CLOSE
