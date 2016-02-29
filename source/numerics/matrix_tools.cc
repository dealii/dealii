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


#include <algorithm>
#include <set>
#include <cmath>


DEAL_II_NAMESPACE_OPEN




namespace MatrixTools
{
  namespace
  {
    template <typename Iterator>
    bool column_less_than(const typename Iterator::value_type p,
                          const unsigned int column)
    {
      return (p.column() < column);
    }
  }

//TODO:[WB] I don't think that the optimized storage of diagonals is needed (GK)
  template <typename number>
  void
  apply_boundary_values (const std::map<types::global_dof_index,number> &boundary_values,
                         SparseMatrix<number>  &matrix,
                         Vector<number>   &solution,
                         Vector<number>   &right_hand_side,
                         const bool        eliminate_columns)
  {
    Assert (matrix.n() == right_hand_side.size(),
            ExcDimensionMismatch(matrix.n(), right_hand_side.size()));
    Assert (matrix.n() == solution.size(),
            ExcDimensionMismatch(matrix.n(), solution.size()));
    Assert (matrix.n() == matrix.m(),
            ExcDimensionMismatch(matrix.n(), matrix.m()));

    // if no boundary values are to be applied
    // simply return
    if (boundary_values.size() == 0)
      return;


    const types::global_dof_index n_dofs = matrix.m();

    // if a diagonal entry is zero
    // later, then we use another
    // number instead. take it to be
    // the first nonzero diagonal
    // element of the matrix, or 1 if
    // there is no such thing
    number first_nonzero_diagonal_entry = 1;
    for (unsigned int i=0; i<n_dofs; ++i)
      if (matrix.diag_element(i) != number())
        {
          first_nonzero_diagonal_entry = matrix.diag_element(i);
          break;
        }


    typename std::map<types::global_dof_index,number>::const_iterator dof  = boundary_values.begin(),
                                                                      endd = boundary_values.end();
    for (; dof != endd; ++dof)
      {
        Assert (dof->first < n_dofs, ExcInternalError());

        const types::global_dof_index dof_number = dof->first;
        // for each boundary dof:

        // set entries of this line to zero except for the diagonal
        // entry
        for (typename SparseMatrix<number>::iterator
             p = matrix.begin(dof_number);
             p != matrix.end(dof_number); ++p)
          if (p->column() != dof_number)
            p->value() = 0.;

        // set right hand side to
        // wanted value: if main diagonal
        // entry nonzero, don't touch it
        // and scale rhs accordingly. If
        // zero, take the first main
        // diagonal entry we can find, or
        // one if no nonzero main diagonal
        // element exists. Normally, however,
        // the main diagonal entry should
        // not be zero.
        //
        // store the new rhs entry to make
        // the gauss step more efficient
        number new_rhs;
        if (matrix.diag_element(dof_number) != number())
          {
            new_rhs = dof->second * matrix.diag_element(dof_number);
            right_hand_side(dof_number) = new_rhs;
          }
        else
          {
            matrix.set (dof_number, dof_number,
                        first_nonzero_diagonal_entry);
            new_rhs = dof->second * first_nonzero_diagonal_entry;
            right_hand_side(dof_number) = new_rhs;
          }


        // if the user wants to have
        // the symmetry of the matrix
        // preserved, and if the
        // sparsity pattern is
        // symmetric, then do a Gauss
        // elimination step with the
        // present row
        if (eliminate_columns)
          {
            // store the only nonzero entry
            // of this line for the Gauss
            // elimination step
            const number diagonal_entry = matrix.diag_element(dof_number);

            // we have to loop over all rows of the matrix which have
            // a nonzero entry in the column which we work in
            // presently. if the sparsity pattern is symmetric, then
            // we can get the positions of these rows cheaply by
            // looking at the nonzero column numbers of the present
            // row. we need not look at the first entry of each row,
            // since that is the diagonal element and thus the present
            // row
            for (typename SparseMatrix<number>::iterator
                 q = matrix.begin(dof_number)+1;
                 q != matrix.end(dof_number); ++q)
              {
                const types::global_dof_index row = q->column();

                // find the position of
                // element
                // (row,dof_number)
                bool (*comp)(const typename SparseMatrix<number>::iterator::value_type p,
                             const unsigned int column)
                  = &column_less_than<typename SparseMatrix<number>::iterator>;
                const typename SparseMatrix<number>::iterator
                p = Utilities::lower_bound(matrix.begin(row)+1,
                                           matrix.end(row),
                                           dof_number,
                                           comp);

                // check whether this line has an entry in the
                // regarding column (check for ==dof_number and !=
                // next_row, since if row==dof_number-1, *p is a
                // past-the-end pointer but points to dof_number
                // anyway...)
                //
                // there should be such an entry! we know this because
                // we have assumed that the sparsity pattern is
                // symmetric and we only walk over those rows for
                // which the current row has a column entry
                Assert ((p != matrix.end(row))
                        &&
                        (p->column() == dof_number),
                        ExcMessage("This function is trying to access an element of the "
                                   "matrix that doesn't seem to exist. Are you using a "
                                   "nonsymmetric sparsity pattern? If so, you are not "
                                   "allowed to set the eliminate_column argument of this "
                                   "function, see the documentation."));

                // correct right hand side
                right_hand_side(row) -= static_cast<number>(p->value()) /
                                        diagonal_entry * new_rhs;

                // set matrix entry to zero
                p->value() = 0.;
              }
          }

        // preset solution vector
        solution(dof_number) = dof->second;
      }
  }



  template <typename number>
  void
  apply_boundary_values (const std::map<types::global_dof_index,number> &boundary_values,
                         BlockSparseMatrix<number>  &matrix,
                         BlockVector<number>   &solution,
                         BlockVector<number>   &right_hand_side,
                         const bool             eliminate_columns)
  {
    const unsigned int blocks = matrix.n_block_rows();

    Assert (matrix.n() == right_hand_side.size(),
            ExcDimensionMismatch(matrix.n(), right_hand_side.size()));
    Assert (matrix.n() == solution.size(),
            ExcDimensionMismatch(matrix.n(), solution.size()));
    Assert (matrix.n_block_rows() == matrix.n_block_cols(),
            ExcNotQuadratic());
    Assert (matrix.get_sparsity_pattern().get_row_indices() ==
            matrix.get_sparsity_pattern().get_column_indices(),
            ExcNotQuadratic());
    Assert (matrix.get_sparsity_pattern().get_column_indices() ==
            solution.get_block_indices (),
            ExcBlocksDontMatch ());
    Assert (matrix.get_sparsity_pattern().get_row_indices() ==
            right_hand_side.get_block_indices (),
            ExcBlocksDontMatch ());

    // if no boundary values are to be applied
    // simply return
    if (boundary_values.size() == 0)
      return;


    const types::global_dof_index n_dofs = matrix.m();

    // if a diagonal entry is zero
    // later, then we use another
    // number instead. take it to be
    // the first nonzero diagonal
    // element of the matrix, or 1 if
    // there is no such thing
    number first_nonzero_diagonal_entry = 0;
    for (unsigned int diag_block=0; diag_block<blocks; ++diag_block)
      {
        for (unsigned int i=0; i<matrix.block(diag_block,diag_block).n(); ++i)
          if (matrix.block(diag_block,diag_block).diag_element(i) != 0)
            {
              first_nonzero_diagonal_entry
                = matrix.block(diag_block,diag_block).diag_element(i);
              break;
            }
        // check whether we have found
        // something in the present
        // block
        if (first_nonzero_diagonal_entry != 0)
          break;
      }
    // nothing found on all diagonal
    // blocks? if so, use 1.0 instead
    if (first_nonzero_diagonal_entry == 0)
      first_nonzero_diagonal_entry = 1;


    typename std::map<types::global_dof_index,number>::const_iterator dof  = boundary_values.begin(),
                                                                      endd = boundary_values.end();
    const BlockSparsityPattern &
    sparsity_pattern = matrix.get_sparsity_pattern();

    // pointer to the mapping between
    // global and block indices. since
    // the row and column mappings are
    // equal, store a pointer on only
    // one of them
    const BlockIndices &
    index_mapping = sparsity_pattern.get_column_indices();

    // now loop over all boundary dofs
    for (; dof != endd; ++dof)
      {
        Assert (dof->first < n_dofs, ExcInternalError());
        (void)n_dofs;

        // get global index and index
        // in the block in which this
        // dof is located
        const types::global_dof_index dof_number = dof->first;
        const std::pair<unsigned int,types::global_dof_index>
        block_index = index_mapping.global_to_local (dof_number);

        // for each boundary dof:

        // set entries of this line
        // to zero except for the diagonal
        // entry. Note that the diagonal
        // entry is always the first one
        // in a row for square matrices
        for (unsigned int block_col=0; block_col<blocks; ++block_col)
          for (typename SparseMatrix<number>::iterator
               p = (block_col == block_index.first ?
                    matrix.block(block_index.first,block_col).begin(block_index.second) + 1 :
                    matrix.block(block_index.first,block_col).begin(block_index.second));
               p != matrix.block(block_index.first,block_col).end(block_index.second);
               ++p)
            p->value() = 0;

        // set right hand side to
        // wanted value: if main diagonal
        // entry nonzero, don't touch it
        // and scale rhs accordingly. If
        // zero, take the first main
        // diagonal entry we can find, or
        // one if no nonzero main diagonal
        // element exists. Normally, however,
        // the main diagonal entry should
        // not be zero.
        //
        // store the new rhs entry to make
        // the gauss step more efficient
        number new_rhs;
        if (matrix.block(block_index.first, block_index.first)
            .diag_element(block_index.second) != 0.0)
          new_rhs = dof->second *
                    matrix.block(block_index.first, block_index.first)
                    .diag_element(block_index.second);
        else
          {
            matrix.block(block_index.first, block_index.first)
            .diag_element(block_index.second)
              = first_nonzero_diagonal_entry;
            new_rhs = dof->second * first_nonzero_diagonal_entry;
          }
        right_hand_side.block(block_index.first)(block_index.second)
          = new_rhs;


        // if the user wants to have
        // the symmetry of the matrix
        // preserved, and if the
        // sparsity pattern is
        // symmetric, then do a Gauss
        // elimination step with the
        // present row. this is a
        // little more complicated for
        // block matrices.
        if (eliminate_columns)
          {
            // store the only nonzero entry
            // of this line for the Gauss
            // elimination step
            const number diagonal_entry
              = matrix.block(block_index.first,block_index.first)
                .diag_element(block_index.second);

            // we have to loop over all
            // rows of the matrix which
            // have a nonzero entry in
            // the column which we work
            // in presently. if the
            // sparsity pattern is
            // symmetric, then we can
            // get the positions of
            // these rows cheaply by
            // looking at the nonzero
            // column numbers of the
            // present row.
            //
            // note that if we check
            // whether row @p{row} in
            // block (r,c) is non-zero,
            // then we have to check
            // for the existence of
            // column @p{row} in block
            // (c,r), i.e. of the
            // transpose block
            for (unsigned int block_row=0; block_row<blocks; ++block_row)
              {
                // get pointers to the sparsity patterns of this block and of
                // the transpose one
                const SparsityPattern &this_sparsity
                  = sparsity_pattern.block (block_row, block_index.first);

                SparseMatrix<number> &this_matrix
                  = matrix.block(block_row, block_index.first);
                SparseMatrix<number> &transpose_matrix
                  = matrix.block(block_index.first, block_row);

                // traverse the row of the transpose block to find the
                // interesting rows in the present block.  don't use the
                // diagonal element of the diagonal block
                for (typename SparseMatrix<number>::iterator
                     q = (block_index.first == block_row ?
                          transpose_matrix.begin(block_index.second)+1 :
                          transpose_matrix.begin(block_index.second));
                     q != transpose_matrix.end(block_index.second);
                     ++q)
                  {
                    // get the number of the column in this row in which a
                    // nonzero entry is. this is also the row of the transpose
                    // block which has an entry in the interesting row
                    const types::global_dof_index row = q->column();

                    // find the position of element (row,dof_number) in this
                    // block (not in the transpose one). note that we have to
                    // take care of special cases with square sub-matrices
                    bool (*comp)(typename SparseMatrix<number>::iterator::value_type p,
                                 const unsigned int column)
                      = &column_less_than<typename SparseMatrix<number>::iterator>;

                    typename SparseMatrix<number>::iterator p = this_matrix.end();

                    if (this_sparsity.n_rows() == this_sparsity.n_cols())
                      {
                        if (this_matrix.begin(row)->column()
                            ==
                            block_index.second)
                          p = this_matrix.begin(row);
                        else
                          p = Utilities::lower_bound(this_matrix.begin(row)+1,
                                                     this_matrix.end(row),
                                                     block_index.second,
                                                     comp);
                      }
                    else
                      p = Utilities::lower_bound(this_matrix.begin(row),
                                                 this_matrix.end(row),
                                                 block_index.second,
                                                 comp);

                    // check whether this line has an entry in the
                    // regarding column (check for ==dof_number and !=
                    // next_row, since if row==dof_number-1, *p is a
                    // past-the-end pointer but points to dof_number
                    // anyway...)
                    //
                    // there should be such an entry! we know this because
                    // we have assumed that the sparsity pattern is
                    // symmetric and we only walk over those rows for
                    // which the current row has a column entry
                    Assert ((p->column() == block_index.second) &&
                            (p != this_matrix.end(row)),
                            ExcInternalError());

                    // correct right hand side
                    right_hand_side.block(block_row)(row)
                    -= p->value() /
                       diagonal_entry * new_rhs;

                    // set matrix entry to zero
                    p->value() = 0.;
                  }
              }
          }

        // preset solution vector
        solution.block(block_index.first)(block_index.second) = dof->second;
      }
  }





  template <typename number>
  void
  local_apply_boundary_values (const std::map<types::global_dof_index,number> &boundary_values,
                               const std::vector<types::global_dof_index> &local_dof_indices,
                               FullMatrix<number> &local_matrix,
                               Vector<number>     &local_rhs,
                               const bool          eliminate_columns)
  {
    Assert (local_dof_indices.size() == local_matrix.m(),
            ExcDimensionMismatch(local_dof_indices.size(),
                                 local_matrix.m()));
    Assert (local_dof_indices.size() == local_matrix.n(),
            ExcDimensionMismatch(local_dof_indices.size(),
                                 local_matrix.n()));
    Assert (local_dof_indices.size() == local_rhs.size(),
            ExcDimensionMismatch(local_dof_indices.size(),
                                 local_rhs.size()));

    // if there is nothing to do, then exit
    // right away
    if (boundary_values.size() == 0)
      return;

    // otherwise traverse all the dofs used in
    // the local matrices and vectors and see
    // what's there to do

    // if we need to treat an entry, then we
    // set the diagonal entry to its absolute
    // value. if it is zero, we used to set it
    // to one, which is a really terrible
    // choice that can lead to hours of
    // searching for bugs in programs (I
    // experienced this :-( ) if the matrix
    // entries are otherwise very large. this
    // is so since iterative solvers would
    // simply not correct boundary nodes for
    // their correct values since the residual
    // contributions of their rows of the
    // linear system is almost zero if the
    // diagonal entry is one. thus, set it to
    // the average absolute value of the
    // nonzero diagonal elements.
    //
    // we only compute this value lazily the
    // first time we need it.
    number average_diagonal = 0;
    const unsigned int n_local_dofs = local_dof_indices.size();
    for (unsigned int i=0; i<n_local_dofs; ++i)
      {
        const typename std::map<types::global_dof_index, number>::const_iterator
        boundary_value = boundary_values.find (local_dof_indices[i]);
        if (boundary_value != boundary_values.end())
          {
            // remove this row, except for the
            // diagonal element
            for (unsigned int j=0; j<n_local_dofs; ++j)
              if (i != j)
                local_matrix(i,j) = 0;

            // replace diagonal entry by its
            // absolute value to make sure that
            // everything remains positive, or
            // by the average diagonal value if
            // zero
            if (local_matrix(i,i) == 0.)
              {
                // if average diagonal hasn't
                // yet been computed, do so now
                if (average_diagonal == 0.)
                  {
                    unsigned int nonzero_diagonals = 0;
                    for (unsigned int k=0; k<n_local_dofs; ++k)
                      if (local_matrix(k,k) != 0.)
                        {
                          average_diagonal += std::fabs(local_matrix(k,k));
                          ++nonzero_diagonals;
                        }
                    if (nonzero_diagonals != 0)
                      average_diagonal /= nonzero_diagonals;
                    else
                      average_diagonal = 0;
                  }

                // only if all diagonal entries
                // are zero, then resort to the
                // last measure: choose one
                if (average_diagonal == 0.)
                  average_diagonal = 1.;

                local_matrix(i,i) = average_diagonal;
              }
            else
              local_matrix(i,i) = std::fabs(local_matrix(i,i));

            // and replace rhs entry by correct
            // value
            local_rhs(i) = local_matrix(i,i) * boundary_value->second;

            // finally do the elimination step
            // if requested
            if (eliminate_columns == true)
              {
                for (unsigned int row=0; row<n_local_dofs; ++row)
                  if (row != i)
                    {
                      local_rhs(row) -= local_matrix(row,i) * boundary_value->second;
                      local_matrix(row,i) = 0;
                    }
              }
          }
      }
  }
}



// explicit instantiations
#include "matrix_tools.inst"


DEAL_II_NAMESPACE_CLOSE
