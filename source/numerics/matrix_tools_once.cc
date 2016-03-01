// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2015 by the deal.II authors
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

#include <deal.II/base/function.h>
#include <deal.II/base/quadrature.h>
#include <deal.II/base/work_stream.h>
#include <deal.II/base/geometry_info.h>
#include <deal.II/base/quadrature.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q1.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/hp/fe_values.h>
#include <deal.II/hp/mapping_collection.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/block_sparse_matrix.h>

#ifdef DEAL_II_WITH_PETSC
#  include <deal.II/lac/petsc_parallel_sparse_matrix.h>
#  include <deal.II/lac/petsc_sparse_matrix.h>
#  include <deal.II/lac/petsc_parallel_vector.h>
#  include <deal.II/lac/petsc_vector.h>
#  include <deal.II/lac/petsc_parallel_block_sparse_matrix.h>
#endif

#ifdef DEAL_II_WITH_TRILINOS
#  include <deal.II/lac/trilinos_sparse_matrix.h>
#  include <deal.II/lac/trilinos_vector.h>
#  include <deal.II/lac/trilinos_block_sparse_matrix.h>
#  include <deal.II/lac/trilinos_block_vector.h>
#endif

#include <algorithm>
#include <set>
#include <cmath>


DEAL_II_NAMESPACE_OPEN

namespace MatrixTools
{

#ifdef DEAL_II_WITH_PETSC

  namespace internal
  {
    namespace PETScWrappers
    {
      template <typename PETScMatrix, typename PETScVector>
      void
      apply_boundary_values (const std::map<types::global_dof_index,PetscScalar> &boundary_values,
                             PETScMatrix      &matrix,
                             PETScVector      &solution,
                             PETScVector      &right_hand_side,
                             const bool        eliminate_columns)
      {
        (void)eliminate_columns;
        Assert (eliminate_columns == false, ExcNotImplemented());

        Assert (matrix.n() == right_hand_side.size(),
                ExcDimensionMismatch(matrix.n(), right_hand_side.size()));
        Assert (matrix.n() == solution.size(),
                ExcDimensionMismatch(matrix.n(), solution.size()));

        // if no boundary values are to be applied, then
        // jump straight to the compress() calls that we still have
        // to perform because they are collective operations
        if (boundary_values.size() > 0)
          {
            const std::pair<types::global_dof_index, types::global_dof_index> local_range
              = matrix.local_range();
            Assert (local_range == right_hand_side.local_range(),
                    ExcInternalError());
            Assert (local_range == solution.local_range(),
                    ExcInternalError());

            // determine the first nonzero diagonal
            // entry from within the part of the
            // matrix that we can see. if we can't
            // find such an entry, take one
            PetscScalar average_nonzero_diagonal_entry = 1;
            for (types::global_dof_index i=local_range.first; i<local_range.second; ++i)
              if (matrix.diag_element(i) != PetscScalar ())
                {
                  average_nonzero_diagonal_entry = std::abs(matrix.diag_element(i));
                  break;
                }

            // figure out which rows of the matrix we
            // have to eliminate on this processor
            std::vector<types::global_dof_index> constrained_rows;
            for (std::map<types::global_dof_index,PetscScalar>::const_iterator
                 dof  = boundary_values.begin();
                 dof != boundary_values.end();
                 ++dof)
              if ((dof->first >= local_range.first) &&
                  (dof->first < local_range.second))
                constrained_rows.push_back (dof->first);

            // then eliminate these rows and set
            // their diagonal entry to what we have
            // determined above. note that for petsc
            // matrices interleaving read with write
            // operations is very expensive. thus, we
            // here always replace the diagonal
            // element, rather than first checking
            // whether it is nonzero and in that case
            // preserving it. this is different from
            // the case of deal.II sparse matrices
            // treated in the other functions.
            matrix.clear_rows (constrained_rows, average_nonzero_diagonal_entry);

            std::vector<types::global_dof_index> indices;
            std::vector<PetscScalar>  solution_values;
            for (std::map<types::global_dof_index,PetscScalar>::const_iterator
                 dof  = boundary_values.begin();
                 dof != boundary_values.end();
                 ++dof)
              if ((dof->first >= local_range.first) &&
                  (dof->first < local_range.second))
                {
                  indices.push_back (dof->first);
                  solution_values.push_back (dof->second);
                }
            solution.set (indices, solution_values);

            // now also set appropriate values for
            // the rhs
            for (unsigned int i=0; i<solution_values.size(); ++i)
              solution_values[i] *= average_nonzero_diagonal_entry;

            right_hand_side.set (indices, solution_values);
          }
        else
          {
            // clear_rows() is a collective operation so we still have to call
            // it:
            std::vector<types::global_dof_index> constrained_rows;
            matrix.clear_rows (constrained_rows, 1.);
          }

        // clean up
        solution.compress (VectorOperation::insert);
        right_hand_side.compress (VectorOperation::insert);
      }
    }
  }



  void
  apply_boundary_values (const std::map<types::global_dof_index,PetscScalar> &boundary_values,
                         PETScWrappers::SparseMatrix   &matrix,
                         PETScWrappers::Vector   &solution,
                         PETScWrappers::Vector   &right_hand_side,
                         const bool        eliminate_columns)
  {
    // simply redirect to the generic function
    // used for both petsc matrix types
    internal::PETScWrappers::apply_boundary_values (boundary_values, matrix, solution,
                                                    right_hand_side, eliminate_columns);
  }



  void
  apply_boundary_values (const std::map<types::global_dof_index,PetscScalar> &boundary_values,
                         PETScWrappers::MPI::SparseMatrix   &matrix,
                         PETScWrappers::MPI::Vector   &solution,
                         PETScWrappers::MPI::Vector   &right_hand_side,
                         const bool        eliminate_columns)
  {
    // simply redirect to the generic function
    // used for both petsc matrix types
    internal::PETScWrappers::apply_boundary_values (boundary_values, matrix, solution,
                                                    right_hand_side, eliminate_columns);
  }


  void
  apply_boundary_values (const std::map<types::global_dof_index,PetscScalar>  &boundary_values,
                         PETScWrappers::MPI::BlockSparseMatrix &matrix,
                         PETScWrappers::MPI::BlockVector        &solution,
                         PETScWrappers::MPI::BlockVector        &right_hand_side,
                         const bool                            eliminate_columns)
  {
    Assert (matrix.n() == right_hand_side.size(),
            ExcDimensionMismatch(matrix.n(), right_hand_side.size()));
    Assert (matrix.n() == solution.size(),
            ExcDimensionMismatch(matrix.n(), solution.size()));
    Assert (matrix.n_block_rows() == matrix.n_block_cols(),
            ExcNotQuadratic());

    const unsigned int n_blocks = matrix.n_block_rows();

    // We need to find the subdivision
    // into blocks for the boundary values.
    // To this end, generate a vector of
    // maps with the respective indices.
    std::vector<std::map<dealii::types::global_dof_index,PetscScalar> > block_boundary_values(n_blocks);
    {
      int block = 0;
      dealii::types::global_dof_index offset = 0;
      for (std::map<types::global_dof_index,PetscScalar>::const_iterator
           dof  = boundary_values.begin();
           dof != boundary_values.end();
           ++dof)
        {
          if (dof->first >= matrix.block(block,0).m() + offset)
            {
              offset += matrix.block(block,0).m();
              block++;
            }
          const types::global_dof_index index = dof->first - offset;
          block_boundary_values[block].insert(std::pair<types::global_dof_index, PetscScalar> (index,dof->second));
        }
    }

    // Now call the non-block variants on
    // the diagonal subblocks and the
    // solution/rhs.
    for (unsigned int block=0; block<n_blocks; ++block)
      internal::PETScWrappers::apply_boundary_values(block_boundary_values[block],
                                                     matrix.block(block,block),
                                                     solution.block(block),
                                                     right_hand_side.block(block),
                                                     eliminate_columns);

    // Finally, we need to do something
    // about the off-diagonal matrices. This
    // is luckily not difficult. Just clear
    // the whole row.
    for (unsigned int block_m=0; block_m<n_blocks; ++block_m)
      {
        const std::pair<types::global_dof_index, types::global_dof_index> local_range
          = matrix.block(block_m,0).local_range();

        std::vector<types::global_dof_index> constrained_rows;
        for (std::map<types::global_dof_index,PetscScalar>::const_iterator
             dof  = block_boundary_values[block_m].begin();
             dof != block_boundary_values[block_m].end();
             ++dof)
          if ((dof->first >= local_range.first) &&
              (dof->first < local_range.second))
            constrained_rows.push_back (dof->first);

        for (unsigned int block_n=0; block_n<n_blocks; ++block_n)
          if (block_m != block_n)
            matrix.block(block_m,block_n).clear_rows(constrained_rows);
      }
  }

#endif



#ifdef DEAL_II_WITH_TRILINOS

  namespace internal
  {
    namespace TrilinosWrappers
    {
      template <typename TrilinosMatrix, typename TrilinosVector>
      void
      apply_boundary_values (const std::map<types::global_dof_index,TrilinosScalar> &boundary_values,
                             TrilinosMatrix      &matrix,
                             TrilinosVector      &solution,
                             TrilinosVector      &right_hand_side,
                             const bool           eliminate_columns)
      {
        Assert (eliminate_columns == false, ExcNotImplemented());
        (void)eliminate_columns;

        Assert (matrix.n() == right_hand_side.size(),
                ExcDimensionMismatch(matrix.n(), right_hand_side.size()));
        Assert (matrix.n() == solution.size(),
                ExcDimensionMismatch(matrix.m(), solution.size()));

        // if no boundary values are to be applied, then
        // jump straight to the compress() calls that we still have
        // to perform because they are collective operations
        if (boundary_values.size() > 0)
          {
            const std::pair<types::global_dof_index, types::global_dof_index> local_range
              = matrix.local_range();
            Assert (local_range == right_hand_side.local_range(),
                    ExcInternalError());
            Assert (local_range == solution.local_range(),
                    ExcInternalError());

            // determine the first nonzero diagonal
            // entry from within the part of the
            // matrix that we can see. if we can't
            // find such an entry, take one
            TrilinosScalar average_nonzero_diagonal_entry = 1;
            for (types::global_dof_index i=local_range.first; i<local_range.second; ++i)
              if (matrix.diag_element(i) != 0)
                {
                  average_nonzero_diagonal_entry = std::fabs(matrix.diag_element(i));
                  break;
                }

            // figure out which rows of the matrix we
            // have to eliminate on this processor
            std::vector<types::global_dof_index> constrained_rows;
            for (std::map<types::global_dof_index,TrilinosScalar>::const_iterator
                 dof  = boundary_values.begin();
                 dof != boundary_values.end();
                 ++dof)
              if ((dof->first >= local_range.first) &&
                  (dof->first < local_range.second))
                constrained_rows.push_back (dof->first);

            // then eliminate these rows and
            // set their diagonal entry to
            // what we have determined
            // above. if the value already is
            // nonzero, it will be preserved,
            // in accordance with the basic
            // matrix classes in deal.II.
            matrix.clear_rows (constrained_rows, average_nonzero_diagonal_entry);

            std::vector<types::global_dof_index> indices;
            std::vector<TrilinosScalar>  solution_values;
            for (std::map<types::global_dof_index,TrilinosScalar>::const_iterator
                 dof  = boundary_values.begin();
                 dof != boundary_values.end();
                 ++dof)
              if ((dof->first >= local_range.first) &&
                  (dof->first < local_range.second))
                {
                  indices.push_back (dof->first);
                  solution_values.push_back (dof->second);
                }
            solution.set (indices, solution_values);

            // now also set appropriate
            // values for the rhs
            for (unsigned int i=0; i<solution_values.size(); ++i)
              solution_values[i] *= matrix.diag_element(indices[i]);

            right_hand_side.set (indices, solution_values);
          }
        else
          {
            // clear_rows() is a collective operation so we still have to call
            // it:
            std::vector<types::global_dof_index> constrained_rows;
            matrix.clear_rows (constrained_rows, 1.);
          }

        // clean up
        matrix.compress (VectorOperation::insert);
        solution.compress (VectorOperation::insert);
        right_hand_side.compress (VectorOperation::insert);
      }



      template <typename TrilinosMatrix, typename TrilinosBlockVector>
      void
      apply_block_boundary_values (const std::map<types::global_dof_index,TrilinosScalar> &boundary_values,
                                   TrilinosMatrix      &matrix,
                                   TrilinosBlockVector &solution,
                                   TrilinosBlockVector &right_hand_side,
                                   const bool          eliminate_columns)
      {
        Assert (eliminate_columns == false, ExcNotImplemented());

        Assert (matrix.n() == right_hand_side.size(),
                ExcDimensionMismatch(matrix.n(), right_hand_side.size()));
        Assert (matrix.n() == solution.size(),
                ExcDimensionMismatch(matrix.n(), solution.size()));
        Assert (matrix.n_block_rows() == matrix.n_block_cols(),
                ExcNotQuadratic());

        const unsigned int n_blocks = matrix.n_block_rows();

        // We need to find the subdivision
        // into blocks for the boundary values.
        // To this end, generate a vector of
        // maps with the respective indices.
        std::vector<std::map<types::global_dof_index,TrilinosScalar> > block_boundary_values(n_blocks);
        {
          int block=0;
          types::global_dof_index offset = 0;
          for (std::map<types::global_dof_index,TrilinosScalar>::const_iterator
               dof  = boundary_values.begin();
               dof != boundary_values.end();
               ++dof)
            {
              if (dof->first >= matrix.block(block,0).m() + offset)
                {
                  offset += matrix.block(block,0).m();
                  block++;
                }
              const types::global_dof_index index = dof->first - offset;
              block_boundary_values[block].insert(
                std::pair<types::global_dof_index, TrilinosScalar> (index,dof->second));
            }
        }

        // Now call the non-block variants on
        // the diagonal subblocks and the
        // solution/rhs.
        for (unsigned int block=0; block<n_blocks; ++block)
          TrilinosWrappers::apply_boundary_values(block_boundary_values[block],
                                                  matrix.block(block,block),
                                                  solution.block(block),
                                                  right_hand_side.block(block),
                                                  eliminate_columns);

        // Finally, we need to do something
        // about the off-diagonal matrices. This
        // is luckily not difficult. Just clear
        // the whole row.
        for (unsigned int block_m=0; block_m<n_blocks; ++block_m)
          {
            const std::pair<types::global_dof_index, types::global_dof_index> local_range
              = matrix.block(block_m,0).local_range();

            std::vector<types::global_dof_index> constrained_rows;
            for (std::map<types::global_dof_index,TrilinosScalar>::const_iterator
                 dof  = block_boundary_values[block_m].begin();
                 dof != block_boundary_values[block_m].end();
                 ++dof)
              if ((dof->first >= local_range.first) &&
                  (dof->first < local_range.second))
                constrained_rows.push_back (dof->first);

            for (unsigned int block_n=0; block_n<n_blocks; ++block_n)
              if (block_m != block_n)
                matrix.block(block_m,block_n).clear_rows(constrained_rows);
          }
      }
    }
  }




  void
  apply_boundary_values (const std::map<types::global_dof_index,TrilinosScalar> &boundary_values,
                         TrilinosWrappers::SparseMatrix   &matrix,
                         TrilinosWrappers::Vector         &solution,
                         TrilinosWrappers::Vector         &right_hand_side,
                         const bool        eliminate_columns)
  {
    // simply redirect to the generic function
    // used for both trilinos matrix types
    internal::TrilinosWrappers::apply_boundary_values (boundary_values, matrix, solution,
                                                       right_hand_side, eliminate_columns);
  }



  void
  apply_boundary_values (const std::map<types::global_dof_index,TrilinosScalar> &boundary_values,
                         TrilinosWrappers::SparseMatrix   &matrix,
                         TrilinosWrappers::MPI::Vector    &solution,
                         TrilinosWrappers::MPI::Vector    &right_hand_side,
                         const bool        eliminate_columns)
  {
    // simply redirect to the generic function
    // used for both trilinos matrix types
    internal::TrilinosWrappers::apply_boundary_values (boundary_values, matrix, solution,
                                                       right_hand_side, eliminate_columns);
  }



  void
  apply_boundary_values (const std::map<types::global_dof_index,TrilinosScalar>  &boundary_values,
                         TrilinosWrappers::BlockSparseMatrix &matrix,
                         TrilinosWrappers::BlockVector        &solution,
                         TrilinosWrappers::BlockVector        &right_hand_side,
                         const bool                            eliminate_columns)
  {
    internal::TrilinosWrappers::apply_block_boundary_values (boundary_values, matrix,
                                                             solution, right_hand_side,
                                                             eliminate_columns);
  }



  void
  apply_boundary_values (const std::map<types::global_dof_index,TrilinosScalar>  &boundary_values,
                         TrilinosWrappers::BlockSparseMatrix &matrix,
                         TrilinosWrappers::MPI::BlockVector   &solution,
                         TrilinosWrappers::MPI::BlockVector   &right_hand_side,
                         const bool                            eliminate_columns)
  {
    internal::TrilinosWrappers::apply_block_boundary_values (boundary_values, matrix,
                                                             solution, right_hand_side,
                                                             eliminate_columns);
  }

#endif

} // namespace MatrixTools

DEAL_II_NAMESPACE_CLOSE
