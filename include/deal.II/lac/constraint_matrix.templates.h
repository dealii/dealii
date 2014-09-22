// ---------------------------------------------------------------------
//
// Copyright (C) 1999 - 2013 by the deal.II authors
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


#ifndef __deal2__constraint_matrix_templates_h
#define __deal2__constraint_matrix_templates_h


#include <deal.II/lac/constraint_matrix.h>

#include <deal.II/base/table.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparsity_pattern.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/block_sparsity_pattern.h>
#include <deal.II/lac/block_sparse_matrix.h>
#include <deal.II/lac/parallel_vector.h>
#include <deal.II/lac/parallel_block_vector.h>
#include <deal.II/lac/petsc_parallel_vector.h>
#include <deal.II/lac/petsc_vector.h>
#include <deal.II/lac/trilinos_vector.h>

#include <iomanip>

DEAL_II_NAMESPACE_OPEN


template<typename number>
void
ConstraintMatrix::condense (const SparseMatrix<number> &uncondensed,
                            SparseMatrix<number>       &condensed) const
{
  // create two dummy vectors and enter the
  // other function
  Vector<number> dummy (0);
  condense (uncondensed, dummy, condensed, dummy);
}



template<typename number>
void
ConstraintMatrix::condense (SparseMatrix<number> &uncondensed) const
{
  Vector<number> dummy (0);
  condense (uncondensed, dummy);
}



template <typename number>
void
ConstraintMatrix::condense (BlockSparseMatrix<number> &uncondensed) const
{
  BlockVector<number> dummy (0);
  condense (uncondensed, dummy);
}



template<class VectorType>
void
ConstraintMatrix::condense (const VectorType &vec_ghosted,
                            VectorType       &vec) const
{
  Assert (sorted == true, ExcMatrixNotClosed());

  // if this is called with different arguments, we need to copy the data over:
  if (&vec != &vec_ghosted)
    vec = vec_ghosted;

  // distribute all entries, and set them to zero. do so in
  // two loops because in the first one we need to add to elements
  // and in the second one we need to set elements to zero. for
  // parallel vectors, this can only work if we can put a compress()
  // in between, but we don't want to call compress() twice per entry
  for (std::vector<ConstraintLine>::const_iterator
       constraint_line = lines.begin();
       constraint_line!=lines.end(); ++constraint_line)
    {
      // in case the constraint is
      // inhomogeneous, this function is not
      // appropriate. Throw an exception.
      Assert (constraint_line->inhomogeneity == 0.,
              ExcMessage ("Inhomogeneous constraint cannot be condensed "
                          "without any matrix specified."));

      const typename VectorType::value_type old_value = vec_ghosted(constraint_line->line);
      for (size_type q=0; q!=constraint_line->entries.size(); ++q)
        if (vec.in_local_range(constraint_line->entries[q].first) == true)
          vec(constraint_line->entries[q].first)
          += (static_cast<typename VectorType::value_type>
              (old_value) *
              constraint_line->entries[q].second);
    }

  vec.compress(VectorOperation::add);

  for (std::vector<ConstraintLine>::const_iterator
       constraint_line = lines.begin();
       constraint_line!=lines.end(); ++constraint_line)
    if (vec.in_local_range(constraint_line->line) == true)
      vec(constraint_line->line) = 0.;

  vec.compress(VectorOperation::insert);
}



template <class VectorType>
void
ConstraintMatrix::condense (VectorType &vec) const
{
  condense(vec, vec);
}



template<typename number, class VectorType>
void
ConstraintMatrix::condense (const SparseMatrix<number> &uncondensed,
                            const VectorType           &uncondensed_vector,
                            SparseMatrix<number>       &condensed,
                            VectorType                 &condensed_vector) const
{
  // check whether we work on real vectors
  // or we just used a dummy when calling
  // the other function above.
  const bool use_vectors = (uncondensed_vector.size() == 0 &&
                            condensed_vector.size() == 0) ? false : true;

  const SparsityPattern &uncondensed_struct = uncondensed.get_sparsity_pattern ();

  Assert (sorted == true, ExcMatrixNotClosed());
  Assert (uncondensed_struct.is_compressed() == true, ExcMatrixNotClosed());
  Assert (condensed.get_sparsity_pattern().is_compressed() == true, ExcMatrixNotClosed());
  Assert (uncondensed_struct.n_rows() == uncondensed_struct.n_cols(),
          ExcNotQuadratic());
  Assert (condensed.n() == condensed.m(),
          ExcNotQuadratic());
  AssertDimension (condensed.n()+n_constraints(), uncondensed.n());
  if (use_vectors == true)
    {
      AssertDimension (condensed_vector.size()+n_constraints(),
                       uncondensed_vector.size());
      AssertDimension (condensed_vector.size(), condensed.m());
    }

  // store for each line of the matrix
  // its new line number
  // after compression. If the shift is
  // -1, this line will be condensed away
  std::vector<int> new_line;

  new_line.reserve (uncondensed_struct.n_rows());

  std::vector<ConstraintLine>::const_iterator next_constraint = lines.begin();
  size_type                                   shift           = 0;
  const size_type n_rows = uncondensed_struct.n_rows();

  if (next_constraint == lines.end())
    // if no constraint is to be handled
    for (size_type row=0; row!=n_rows; ++row)
      new_line.push_back (row);
  else
    for (size_type row=0; row!=n_rows; ++row)
      if (row == next_constraint->line)
        {
          // this line is constrained
          new_line.push_back (-1);
          // note that @p lines is ordered
          ++shift;
          ++next_constraint;
          if (next_constraint == lines.end())
            // nothing more to do; finish rest
            // of loop
            {
              for (size_type i=row+1; i<n_rows; ++i)
                new_line.push_back (i-shift);
              break;
            };
        }
      else
        new_line.push_back (row-shift);


  next_constraint = lines.begin();

  // note: in this loop we need not check
  // whether @p next_constraint is a valid
  // iterator, since @p next_constraint is
  // only evaluated so often as there are
  // entries in new_line[*] which tells us
  // which constraints exist
  for (size_type row=0; row<uncondensed_struct.n_rows(); ++row)
    if (new_line[row] != -1)
      {
        // line not constrained
        // copy entries if column will not
        // be condensed away, distribute
        // otherwise
        for (typename SparseMatrix<number>::const_iterator
             p = uncondensed.begin(row);
             p != uncondensed.end(row); ++p)
          if (new_line[p->column()] != -1)
            condensed.add (new_line[row],
                           new_line[p->column()],
                           p->value());
          else
            {
              // let c point to the
              // constraint of this column
              std::vector<ConstraintLine>::const_iterator c = lines.begin();
              while (c->line != p->column())
                ++c;

              for (size_type q=0; q!=c->entries.size(); ++q)
                // distribute to rows with
                // appropriate weight
                condensed.add (new_line[row], new_line[c->entries[q].first],
                               p->value() * c->entries[q].second);

              // take care of inhomogeneity:
              // need to subtract this element from the
              // vector. this corresponds to an
              // explicit elimination in the respective
              // row of the inhomogeneous constraint in
              // the matrix with Gauss elimination
              if (use_vectors == true)
                condensed_vector(new_line[row]) -= p->value() *
                                                   c->inhomogeneity;
            }

        if (use_vectors == true)
          condensed_vector(new_line[row]) += uncondensed_vector(row);
      }
    else
      // line must be distributed
      {
        for (typename SparseMatrix<number>::const_iterator
             p = uncondensed.begin(row);
             p != uncondensed.end(row); ++p)
          // for each column: distribute
          if (new_line[p->column()] != -1)
            // column is not constrained
            for (size_type q=0; q!=next_constraint->entries.size(); ++q)
              condensed.add (new_line[next_constraint->entries[q].first],
                             new_line[p->column()],
                             p->value() *
                             next_constraint->entries[q].second);

          else
            // not only this line but
            // also this col is constrained
            {
              // let c point to the constraint
              // of this column
              std::vector<ConstraintLine>::const_iterator c = lines.begin();
              while (c->line != p->column())
                ++c;

              for (size_type r=0; r!=c->entries.size(); ++r)
                for (size_type q=0; q!=next_constraint->entries.size(); ++q)
                  condensed.add (new_line[next_constraint->entries[q].first],
                                 new_line[c->entries[r].first],
                                 p->value() *
                                 next_constraint->entries[r].second *
                                 c->entries[r].second);

              if (use_vectors == true)
                for (size_type q=0; q!=next_constraint->entries.size(); ++q)
                  condensed_vector (new_line[next_constraint->entries[q].first])
                  -= p->value() *
                     next_constraint->entries[q].second *
                     c->inhomogeneity;
            }

        // condense the vector
        if (use_vectors == true)
          for (size_type q=0; q!=next_constraint->entries.size(); ++q)
            condensed_vector(new_line[next_constraint->entries[q].first])
            +=
              uncondensed_vector(row) * next_constraint->entries[q].second;

        ++next_constraint;
      };
}



template<typename number, class VectorType>
void
ConstraintMatrix::condense (SparseMatrix<number> &uncondensed,
                            VectorType           &vec) const
{
  // check whether we work on real vectors
  // or we just used a dummy when calling
  // the other function above.
  const bool use_vectors = vec.size() == 0 ? false : true;

  const SparsityPattern &sparsity = uncondensed.get_sparsity_pattern ();

  Assert (sorted == true, ExcMatrixNotClosed());
  Assert (sparsity.is_compressed() == true, ExcMatrixNotClosed());
  Assert (sparsity.n_rows() == sparsity.n_cols(),
          ExcNotQuadratic());
  if (use_vectors == true)
    AssertDimension (vec.size(), sparsity.n_rows());

  double average_diagonal = 0;
  for (size_type i=0; i<uncondensed.m(); ++i)
    average_diagonal += std::fabs (uncondensed.diag_element(i));
  average_diagonal /= uncondensed.m();

  // store for each index whether it must be
  // distributed or not. If entry is
  // invalid_size_type, no distribution is
  // necessary.  otherwise, the number states
  // which line in the constraint matrix
  // handles this index
  std::vector<size_type> distribute (sparsity.n_rows(),
                                     numbers::invalid_size_type);

  for (size_type c=0; c<lines.size(); ++c)
    distribute[lines[c].line] = c;

  const size_type n_rows = sparsity.n_rows();
  for (size_type row=0; row<n_rows; ++row)
    {
      if (distribute[row] == numbers::invalid_size_type)
        // regular line. loop over cols
        {
          for (typename SparseMatrix<number>::iterator
               entry = uncondensed.begin(row);
               entry != uncondensed.end(row); ++entry)
            {
              const size_type column = entry->column();

              // end of row reached?
              // this should not
              // happen, since we only
              // operate on compressed
              // matrices!
              Assert (column != SparsityPattern::invalid_entry,
                      ExcMatrixNotClosed());

              if (distribute[column] != numbers::invalid_size_type)
                // distribute entry at
                // regular row @p row
                // and irregular column
                // sparsity.get_column_numbers()[j];
                // set old entry to
                // zero
                {
                  for (size_type q=0;
                       q!=lines[distribute[column]].entries.size(); ++q)
                    uncondensed.add (row,
                                     lines[distribute[column]].entries[q].first,
                                     entry->value() *
                                     lines[distribute[column]].entries[q].second);

                  // need to subtract this element from the
                  // vector. this corresponds to an
                  // explicit elimination in the respective
                  // row of the inhomogeneous constraint in
                  // the matrix with Gauss elimination
                  if (use_vectors == true)
                    vec(row) -=
                      entry->value() * lines[distribute[column]].inhomogeneity;

                  // set old value to zero
                  entry->value() = 0.;
                }
            }
        }
      else
        // row must be distributed
        {
          for (typename SparseMatrix<number>::iterator
               entry = uncondensed.begin(row);
               entry != uncondensed.end(row); ++entry)
            {
              const size_type column = entry->column();

              // end of row reached?
              // this should not
              // happen, since we only
              // operate on compressed
              // matrices!
              Assert (column != SparsityPattern::invalid_entry,
                      ExcMatrixNotClosed());

              if (distribute[column] == numbers::invalid_size_type)
                // distribute entry at
                // irregular row
                // @p row and regular
                // column
                // column. set
                // old entry to zero
                {
                  for (size_type q=0;
                       q!=lines[distribute[row]].entries.size(); ++q)
                    uncondensed.add (lines[distribute[row]].entries[q].first,
                                     column,
                                     entry->value() *
                                     lines[distribute[row]].entries[q].second);

                  // set old entry to zero
                  entry->value() = 0.;
                }
              else
                // distribute entry at
                // irregular row @p row and
                // irregular column
                // @p column set old entry
                // to one on main
                // diagonal, zero otherwise
                {
                  for (size_type p=0; p!=lines[distribute[row]].entries.size(); ++p)
                    {
                      for (size_type q=0;
                           q!=lines[distribute[column]].entries.size(); ++q)
                        uncondensed.add (lines[distribute[row]].entries[p].first,
                                         lines[distribute[column]].entries[q].first,
                                         entry->value() *
                                         lines[distribute[row]].entries[p].second *
                                         lines[distribute[column]].entries[q].second);

                      if (use_vectors == true)
                        vec(lines[distribute[row]].entries[p].first) -=
                          entry->value() * lines[distribute[row]].entries[p].second *
                          lines[distribute[column]].inhomogeneity;
                    }

                  // set old entry to correct
                  // value
                  entry->value() = (row == column ? average_diagonal : 0. );
                }
            }

          // take care of vector
          if (use_vectors == true)
            {
              for (size_type q=0; q!=lines[distribute[row]].entries.size(); ++q)
                vec(lines[distribute[row]].entries[q].first)
                += (vec(row) * lines[distribute[row]].entries[q].second);

              vec(lines[distribute[row]].line) = 0.;
            }
        }
    }
}



template <typename number, class BlockVectorType>
void
ConstraintMatrix::condense (BlockSparseMatrix<number> &uncondensed,
                            BlockVectorType           &vec) const
{
  // check whether we work on real vectors
  // or we just used a dummy when calling
  // the other function above.
  const bool use_vectors = vec.n_blocks() == 0 ? false : true;

  const size_type blocks = uncondensed.n_block_rows();

  const BlockSparsityPattern &
  sparsity = uncondensed.get_sparsity_pattern ();

  Assert (sorted == true, ExcMatrixNotClosed());
  Assert (sparsity.is_compressed() == true, ExcMatrixNotClosed());
  Assert (sparsity.n_rows() == sparsity.n_cols(),
          ExcNotQuadratic());
  Assert (sparsity.n_block_rows() == sparsity.n_block_cols(),
          ExcNotQuadratic());
  Assert (sparsity.n_block_rows() == sparsity.n_block_cols(),
          ExcNotQuadratic());
  Assert (sparsity.get_column_indices() == sparsity.get_row_indices(),
          ExcNotQuadratic());

  if (use_vectors == true)
    {
      AssertDimension (vec.size(), sparsity.n_rows());
      AssertDimension (vec.n_blocks(), sparsity.n_block_rows());
    }

  double average_diagonal = 0;
  for (size_type b=0; b<uncondensed.n_block_rows(); ++b)
    for (size_type i=0; i<uncondensed.block(b,b).m(); ++i)
      average_diagonal += std::fabs (uncondensed.block(b,b).diag_element(i));
  average_diagonal /= uncondensed.m();

  const BlockIndices &
  index_mapping = sparsity.get_column_indices();

  // store for each index whether it must be
  // distributed or not. If entry is
  // numbers::invalid_size_type,
  // no distribution is necessary.
  // otherwise, the number states which line
  // in the constraint matrix handles this
  // index
  std::vector<size_type> distribute (sparsity.n_rows(),
                                     numbers::invalid_size_type);

  for (size_type c=0; c<lines.size(); ++c)
    distribute[lines[c].line] = c;

  const size_type n_rows = sparsity.n_rows();
  for (size_type row=0; row<n_rows; ++row)
    {
      // get index of this row
      // within the blocks
      const std::pair<size_type,size_type>
      block_index = index_mapping.global_to_local(row);
      const size_type block_row = block_index.first;

      if (distribute[row] == numbers::invalid_size_type)
        // regular line. loop over
        // all columns and see
        // whether this column must
        // be distributed
        {

          // to loop over all entries
          // in this row, we have to
          // loop over all blocks in
          // this blockrow and the
          // corresponding row
          // therein
          for (size_type block_col=0; block_col<blocks; ++block_col)
            {
              for (typename SparseMatrix<number>::iterator
                   entry = uncondensed.block(block_row, block_col).begin(block_index.second);
                   entry != uncondensed.block(block_row, block_col).end(block_index.second);
                   ++entry)
                {
                  const size_type global_col
                    = index_mapping.local_to_global(block_col,entry->column());

                  if (distribute[global_col] != numbers::invalid_size_type)
                    // distribute entry at
                    // regular row @p row
                    // and irregular column
                    // global_col; set old
                    // entry to zero
                    {
                      const double old_value = entry->value ();

                      for (size_type q=0;
                           q!=lines[distribute[global_col]].entries.size(); ++q)
                        uncondensed.add (row,
                                         lines[distribute[global_col]].entries[q].first,
                                         old_value *
                                         lines[distribute[global_col]].entries[q].second);

                      // need to subtract this element from the
                      // vector. this corresponds to an
                      // explicit elimination in the respective
                      // row of the inhomogeneous constraint in
                      // the matrix with Gauss elimination
                      if (use_vectors == true)
                        vec(row) -= entry->value() *
                                    lines[distribute[global_col]].inhomogeneity;

                      entry->value() = 0.;
                    }
                }
            }
        }
      else
        {
          // row must be
          // distributed. split the
          // whole row into the
          // chunks defined by the
          // blocks
          for (size_type block_col=0; block_col<blocks; ++block_col)
            {
              for (typename SparseMatrix<number>::iterator
                   entry = uncondensed.block(block_row, block_col).begin(block_index.second);
                   entry != uncondensed.block(block_row, block_col).end(block_index.second);
                   ++entry)
                {
                  const size_type global_col
                    = index_mapping.local_to_global (block_col, entry->column());

                  if (distribute[global_col] ==
                      numbers::invalid_size_type)
                    // distribute
                    // entry at
                    // irregular
                    // row @p row
                    // and regular
                    // column
                    // global_col. set
                    // old entry to
                    // zero
                    {
                      const double old_value = entry->value();

                      for (size_type q=0;
                           q!=lines[distribute[row]].entries.size(); ++q)
                        uncondensed.add (lines[distribute[row]].entries[q].first,
                                         global_col,
                                         old_value *
                                         lines[distribute[row]].entries[q].second);

                      entry->value() = 0.;
                    }
                  else
                    // distribute entry at
                    // irregular row @p row
                    // and irregular column
                    // @p global_col set old
                    // entry to one if on
                    // main diagonal, zero
                    // otherwise
                    {
                      const double old_value = entry->value ();

                      for (size_type p=0; p!=lines[distribute[row]].entries.size(); ++p)
                        {
                          for (size_type q=0; q!=lines[distribute[global_col]].entries.size(); ++q)
                            uncondensed.add (lines[distribute[row]].entries[p].first,
                                             lines[distribute[global_col]].entries[q].first,
                                             old_value *
                                             lines[distribute[row]].entries[p].second *
                                             lines[distribute[global_col]].entries[q].second);

                          if (use_vectors == true)
                            vec(lines[distribute[row]].entries[p].first) -=
                              old_value * lines[distribute[row]].entries[p].second *
                              lines[distribute[global_col]].inhomogeneity;
                        }

                      entry->value() = (row == global_col ? average_diagonal : 0. );
                    }
                }
            }

          // take care of vector
          if (use_vectors == true)
            {
              for (size_type q=0; q!=lines[distribute[row]].entries.size(); ++q)
                vec(lines[distribute[row]].entries[q].first)
                += (vec(row) * lines[distribute[row]].entries[q].second);

              vec(lines[distribute[row]].line) = 0.;
            }
        }
    }
}


//TODO: I'm sure the followng could be made more elegant by using a bit of
//introspection using static member variables of the various vector
//classes to dispatch between the different functions, rather than using
//knowledge of the individual types

// number of functions to select the right implementation for set_zero().
namespace internal
{
  namespace ConstraintMatrix
  {
    namespace
    {
      typedef types::global_dof_index size_type;

      template<class VEC>
      void set_zero_parallel(const std::vector<size_type> &cm, VEC &vec, size_type shift = 0)
      {
        Assert(!vec.has_ghost_elements(), ExcInternalError());
        IndexSet locally_owned = vec.locally_owned_elements();
        for (typename std::vector<size_type>::const_iterator it = cm.begin();
             it != cm.end(); ++it)
          {
            // If shift>0 then we are working on a part of a BlockVector
            // so vec(i) is actually the global entry i+shift.
            // We first make sure the line falls into the range of vec,
            // then check if is part of the local part of the vector, before
            // finally setting it to 0.
            if ((*it)<shift)
              continue;
            size_type idx = *it - shift;
            if (idx<vec.size() && locally_owned.is_element(idx))
              vec(idx) = 0.;
          }
      }

      template<typename Number>
      void set_zero_parallel(const std::vector<size_type> &cm, parallel::distributed::Vector<Number> &vec, size_type shift = 0)
      {
        for (typename std::vector<size_type>::const_iterator it = cm.begin();
             it != cm.end(); ++it)
          {
            // If shift>0 then we are working on a part of a BlockVector
            // so vec(i) is actually the global entry i+shift.
            // We first make sure the line falls into the range of vec,
            // then check if is part of the local part of the vector, before
            // finally setting it to 0.
            if ((*it)<shift)
              continue;
            size_type idx = *it - shift;
            if (vec.in_local_range(idx))
              vec(idx) = 0.;
          }
        vec.zero_out_ghosts();
      }

      template<class VEC>
      void set_zero_in_parallel(const std::vector<size_type> &cm, VEC &vec, internal::bool2type<false>)
      {
        set_zero_parallel(cm, vec, 0);
      }

      // in parallel for BlockVectors
      template<class VEC>
      void set_zero_in_parallel(const std::vector<size_type> &cm, VEC &vec, internal::bool2type<true>)
      {
        size_type start_shift = 0;
        for (size_type j=0; j<vec.n_blocks(); ++j)
          {
            set_zero_parallel(cm, vec.block(j), start_shift);
            start_shift += vec.block(j).size();
          }
      }

      template<class VEC>
      void set_zero_serial(const std::vector<size_type> &cm, VEC &vec)
      {
        for (typename std::vector<size_type>::const_iterator it = cm.begin();
             it != cm.end(); ++it)
          vec(*it) = 0.;
      }

      template<class VEC>
      void set_zero_all(const std::vector<size_type> &cm, VEC &vec)
      {
        set_zero_in_parallel<VEC>(cm, vec, internal::bool2type<IsBlockVector<VEC>::value>());
        vec.compress(VectorOperation::insert);
      }


      template<class T>
      void set_zero_all(const std::vector<size_type> &cm, dealii::Vector<T> &vec)
      {
        set_zero_serial(cm, vec);
      }

      template<class T>
      void set_zero_all(const std::vector<size_type> &cm, dealii::BlockVector<T> &vec)
      {
        set_zero_serial(cm, vec);
      }
    }
  }
}


template <class VectorType>
void
ConstraintMatrix::set_zero (VectorType &vec) const
{
  // since we lines is a private member, we cannot pass it to the functions
  // above. therefore, copy the content which is cheap
  std::vector<size_type> constrained_lines(lines.size());
  for (unsigned int i=0; i<lines.size(); ++i)
    constrained_lines[i] = lines[i].line;
  internal::ConstraintMatrix::set_zero_all(constrained_lines, vec);
}




template <typename VectorType>
void
ConstraintMatrix::
distribute_local_to_global (const Vector<double>            &local_vector,
                            const std::vector<size_type> &local_dof_indices,
                            VectorType                      &global_vector,
                            const FullMatrix<double>        &local_matrix) const
{
  Assert (sorted == true, ExcMatrixNotClosed());
  AssertDimension (local_vector.size(), local_dof_indices.size());
  AssertDimension (local_matrix.m(), local_dof_indices.size());
  AssertDimension (local_matrix.n(), local_dof_indices.size());

  const size_type n_local_dofs = local_vector.size();
  if (lines.empty())
    global_vector.add(local_dof_indices, local_vector);
  else
    for (size_type i=0; i<n_local_dofs; ++i)
      {
        // check whether the current index is
        // constrained. if not, just write the entry
        // into the vector. otherwise, need to resolve
        // the constraint
        if (is_constrained(local_dof_indices[i]) == false)
          {
            global_vector(local_dof_indices[i]) += local_vector(i);
            continue;
          }

        // find the constraint line to the given
        // global dof index
        const size_type line_index = calculate_line_index (local_dof_indices[i]);
        const ConstraintLine *position =
          lines_cache.size() <= line_index ? 0 : &lines[lines_cache[line_index]];

        // Gauss elimination of the matrix columns
        // with the inhomogeneity. Go through them one
        // by one and again check whether they are
        // constrained. If so, distribute the constraint
        const double val = position->inhomogeneity;
        if (val != 0)
          for (size_type j=0; j<n_local_dofs; ++j)
            if (is_constrained(local_dof_indices[j]) == false)
              global_vector(local_dof_indices[j]) -= val * local_matrix(j,i);
            else
              {
                const double matrix_entry = local_matrix(j,i);
                if (matrix_entry == 0)
                  continue;

                const ConstraintLine &position_j =
                  lines[lines_cache[calculate_line_index(local_dof_indices[j])]];
                for (size_type q=0; q<position_j.entries.size(); ++q)
                  {
                    Assert (!(!local_lines.size()
                              || local_lines.is_element(position_j.entries[q].first))
                            || is_constrained(position_j.entries[q].first) == false,
                            ExcMessage ("Tried to distribute to a fixed dof."));
                    global_vector(position_j.entries[q].first)
                    -= val * position_j.entries[q].second * matrix_entry;
                  }
              }

        // now distribute the constraint,
        // but make sure we don't touch
        // the entries of fixed dofs
        for (size_type j=0; j<position->entries.size(); ++j)
          {
            Assert (!(!local_lines.size()
                      || local_lines.is_element(position->entries[j].first))
                    || is_constrained(position->entries[j].first) == false,
                    ExcMessage ("Tried to distribute to a fixed dof."));
            global_vector(position->entries[j].first)
            += local_vector(i) * position->entries[j].second;
          }
      }
}



template<class VectorType>
void
ConstraintMatrix::distribute (const VectorType &condensed,
                              VectorType       &uncondensed) const
{
  Assert (sorted == true, ExcMatrixNotClosed());
  AssertDimension (condensed.size()+n_constraints(), uncondensed.size());

  // store for each line of the new vector
  // its old line number before
  // distribution. If the shift is
  // -1, this line was condensed away
  std::vector<int> old_line;

  old_line.reserve (uncondensed.size());

  std::vector<ConstraintLine>::const_iterator next_constraint = lines.begin();
  size_type                                   shift           = 0;
  size_type n_rows = uncondensed.size();

  if (next_constraint == lines.end())
    // if no constraint is to be handled
    for (size_type row=0; row!=n_rows; ++row)
      old_line.push_back (row);
  else
    for (size_type row=0; row!=n_rows; ++row)
      if (row == next_constraint->line)
        {
          // this line is constrained
          old_line.push_back (-1);
          // note that @p lines is ordered
          ++shift;
          ++next_constraint;
          if (next_constraint == lines.end())
            // nothing more to do; finish rest
            // of loop
            {
              for (size_type i=row+1; i<n_rows; ++i)
                old_line.push_back (i-shift);
              break;
            }
        }
      else
        old_line.push_back (row-shift);


  next_constraint = lines.begin();
  // note: in this loop we need not check
  // whether @p next_constraint is a valid
  // iterator, since @p next_constraint is
  // only evaluated so often as there are
  // entries in new_line[*] which tells us
  // which constraints exist
  for (size_type line=0; line<uncondensed.size(); ++line)
    if (old_line[line] != -1)
      // line was not condensed away
      uncondensed(line) = condensed(old_line[line]);
    else
      {
        // line was condensed away,
        // create it newly. first set
        // it to zero
        uncondensed(line) = next_constraint->inhomogeneity;
        // then add the different
        // contributions
        for (size_type i=0; i<next_constraint->entries.size(); ++i)
          uncondensed(line) += (condensed(old_line[next_constraint->entries[i].first]) *
                                next_constraint->entries[i].second);
        ++next_constraint;
      };
}


namespace internal
{
  namespace
  {
    // create an output vector that consists of the input vector's locally owned
    // elements plus some ghost elements that need to be imported from elsewhere
    //
    // this is an operation that is different for all vector types and so we
    // need a few overloads
#ifdef DEAL_II_WITH_TRILINOS
    void
    import_vector_with_ghost_elements (const TrilinosWrappers::MPI::Vector &vec,
                                       const IndexSet                      &/*locally_owned_elements*/,
                                       const IndexSet                      &needed_elements,
                                       TrilinosWrappers::MPI::Vector       &output,
                                       const internal::bool2type<false>     /*is_block_vector*/)
    {
      Assert(!vec.has_ghost_elements(),
             TrilinosWrappers::VectorBase::ExcGhostsPresent());
#ifdef DEAL_II_WITH_MPI
      const Epetra_MpiComm *mpi_comm
        = dynamic_cast<const Epetra_MpiComm *>(&vec.trilinos_vector().Comm());

      Assert (mpi_comm != 0, ExcInternalError());
      output.reinit (needed_elements, mpi_comm->GetMpiComm());
#else
      output.reinit (needed_elements, MPI_COMM_WORLD);
#endif
      output = vec;
    }
#endif

#ifdef DEAL_II_WITH_PETSC
    void
    import_vector_with_ghost_elements (const PETScWrappers::MPI::Vector &vec,
                                       const IndexSet                   &locally_owned_elements,
                                       const IndexSet                   &needed_elements,
                                       PETScWrappers::MPI::Vector       &output,
                                       const internal::bool2type<false>  /*is_block_vector*/)
    {
      output.reinit (vec.get_mpi_communicator(), locally_owned_elements, needed_elements);
      output = vec;
    }
#endif

    template <typename number>
    void
    import_vector_with_ghost_elements (const parallel::distributed::Vector<number> &vec,
                                       const IndexSet                              &locally_owned_elements,
                                       const IndexSet                              &needed_elements,
                                       parallel::distributed::Vector<number>       &output,
                                       const internal::bool2type<false>             /*is_block_vector*/)
    {
      // TODO: the in vector might already have all elements. need to find a
      // way to efficiently avoid the copy then
      const_cast<parallel::distributed::Vector<number>&>(vec).zero_out_ghosts();
      output.reinit (locally_owned_elements, needed_elements, vec.get_mpi_communicator());
      output = vec;
      output.update_ghost_values();
    }


    // all other vector non-block vector types are sequential and we should
    // not have this function called at all -- so throw an exception
    template <typename Vector>
    void
    import_vector_with_ghost_elements (const Vector                     &/*vec*/,
                                       const IndexSet                   &/*locally_owned_elements*/,
                                       const IndexSet                   &/*needed_elements*/,
                                       Vector                           &/*output*/,
                                       const internal::bool2type<false>  /*is_block_vector*/)
    {
      Assert (false, ExcMessage ("We shouldn't even get here!"));
    }


    // for block vectors, simply dispatch to the individual blocks
    template <class VectorType>
    void
    import_vector_with_ghost_elements (const VectorType                &vec,
                                       const IndexSet                  &locally_owned_elements,
                                       const IndexSet                  &needed_elements,
                                       VectorType                      &output,
                                       const internal::bool2type<true>  /*is_block_vector*/)
    {
      output.reinit (vec.n_blocks());

      types::global_dof_index block_start = 0;
      for (unsigned int b=0; b<vec.n_blocks(); ++b)
        {
          import_vector_with_ghost_elements (vec.block(b),
                                             locally_owned_elements.get_view (block_start, block_start+vec.block(b).size()),
                                             needed_elements.get_view (block_start, block_start+vec.block(b).size()),
                                             output.block(b),
                                             internal::bool2type<false>());
          block_start += vec.block(b).size();
        }

      output.collect_sizes ();
    }
  }
}


template <class VectorType>
void
ConstraintMatrix::distribute (VectorType &vec) const
{
  Assert (sorted==true, ExcMatrixNotClosed());

  // if the vector type supports parallel storage and if the vector actually
  // does store only part of the vector, distributing is slightly more
  // complicated. we might be able to skip the complicated part if one
  // processor owns everything and pretend that this is a sequential vector,
  // but it is difficult for the other processors to know whether they should
  // not do anything or if other processors will create a temporary vector,
  // exchange data (requiring communication, maybe even with the processors
  // that do not own anything because of that particular parallel model), and
  // call compress() finally. the first case here is for the complicated case,
  // the last else is for the simple case (sequential vector)
  const IndexSet vec_owned_elements = vec.locally_owned_elements();
  if (vec.supports_distributed_data == true)
    {
      // This processor owns only part of the vector. one may think that
      // every processor should be able to simply communicate those elements
      // it owns and for which it knows that they act as sources to constrained
      // DoFs to the owner of these DoFs. This would lead to a scheme where all
      // we need to do is to add some local elements to (possibly non-local) ones
      // and then call compress().
      //
      // Alas, this scheme does not work as evidenced by the disaster of bug #51,
      // see http://code.google.com/p/dealii/issues/detail?id=51 and the
      // reversion of one attempt that implements this in r29662. Rather, we
      // need to get a vector that has all the *sources* or constraints we
      // own locally, possibly as ghost vector elements, then read from them,
      // and finally throw away the ghosted vector. Implement this in the following.
      IndexSet needed_elements = vec_owned_elements;

      typedef std::vector<ConstraintLine>::const_iterator constraint_iterator;
      for (constraint_iterator it = lines.begin();
           it != lines.end(); ++it)
        if (vec_owned_elements.is_element(it->line))
          for (unsigned int i=0; i<it->entries.size(); ++i)
            if (!vec_owned_elements.is_element(it->entries[i].first))
              needed_elements.add_index(it->entries[i].first);

      VectorType ghosted_vector;
      internal::import_vector_with_ghost_elements (vec,
                                                   vec_owned_elements, needed_elements,
                                                   ghosted_vector,
                                                   internal::bool2type<IsBlockVector<VectorType>::value>());

      for (constraint_iterator it = lines.begin();
           it != lines.end(); ++it)
        if (vec_owned_elements.is_element(it->line))
          {
            typename VectorType::value_type
            new_value = it->inhomogeneity;
            for (unsigned int i=0; i<it->entries.size(); ++i)
              new_value += (static_cast<typename VectorType::value_type>
                            (ghosted_vector(it->entries[i].first)) *
                            it->entries[i].second);
            Assert(numbers::is_finite(new_value), ExcNumberNotFinite());
            vec(it->line) = new_value;
          }

      // now compress to communicate the entries that we added to
      // and that weren't to local processors to the owner
      //
      // this shouldn't be strictly necessary but it probably doesn't
      // hurt either
      vec.compress (VectorOperation::insert);
    }
  else
    // purely sequential vector (either because the type doesn't
    // support anything else or because it's completely stored
    // locally)
    {
      std::vector<ConstraintLine>::const_iterator next_constraint = lines.begin();
      for (; next_constraint != lines.end(); ++next_constraint)
        {
          // fill entry in line
          // next_constraint.line by adding the
          // different contributions
          typename VectorType::value_type
          new_value = next_constraint->inhomogeneity;
          for (unsigned int i=0; i<next_constraint->entries.size(); ++i)
            new_value += (static_cast<typename VectorType::value_type>
                          (vec(next_constraint->entries[i].first)) *
                          next_constraint->entries[i].second);
          Assert(numbers::is_finite(new_value), ExcNumberNotFinite());
          vec(next_constraint->line) = new_value;
        }
    }
}



// Some helper definitions for the local_to_global functions.
namespace internals
{
  typedef types::global_dof_index size_type;

  // this struct contains all the information we need to store about each of
  // the global entries (global_row): are they obtained directly by some local
  // entry (local_row) or some constraints (constraint_position). This is not
  // directly used in the user code, but accessed via the GlobalRowsFromLocal.
  //
  // The actions performed here correspond to reshaping the constraint
  // information from global degrees of freedom to local ones (i.e.,
  // cell-related DoFs), and also transforming the constraint information from
  // compressed row storage (each local dof that is constrained has a list of
  // constraint entries associated to it) into compressed column storage based
  // on the cell-related DoFs (we have a list of global degrees of freedom,
  // and to each we have a list of local rows where the entries come from). To
  // increase the speed, we additionally store whether an entry is generated
  // directly from the local degrees of freedom or whether it comes from a
  // constraint.
  struct Distributing
  {
    Distributing (const size_type global_row = numbers::invalid_size_type,
                  const size_type local_row = numbers::invalid_size_type);
    Distributing (const Distributing &in);
    Distributing &operator = (const Distributing &in);
    bool operator < (const Distributing &in) const
    {
      return global_row<in.global_row;
    };

    size_type global_row;
    size_type local_row;
    mutable size_type constraint_position;
  };

  inline
  Distributing::Distributing (const size_type global_row,
                              const size_type local_row) :
    global_row (global_row),
    local_row (local_row),
    constraint_position (numbers::invalid_size_type) {}

  inline
  Distributing::Distributing (const Distributing &in)
    :
    constraint_position (numbers::invalid_size_type)
  {
    *this = (in);
  }

  inline
  Distributing &Distributing::operator = (const Distributing &in)
  {
    global_row = in.global_row;
    local_row = in.local_row;
    // the constraints pointer should not contain any data here.
    Assert (constraint_position == numbers::invalid_size_type,
            ExcInternalError());

    if (in.constraint_position != numbers::invalid_size_type)
      {
        constraint_position = in.constraint_position;
        in.constraint_position = numbers::invalid_size_type;
      }
    return *this;
  }



  // this is a cache for constraints that are encountered on a local level.
  // The functionality is similar to
  // std::vector<std::vector<std::pair<uint,double> > >, but tuned so that
  // frequent memory allocation for each entry is avoided. The data is put
  // into a std::vector<std::pair<uint,double> > and the row length is kept
  // fixed at row_length. Both the number of rows and the row length can
  // change is this structure is filled. In that case, the data is
  // rearranged. This is not directly used in the user code, but accessed via
  // the GlobalRowsFromLocal.
  struct DataCache
  {
    DataCache ()
      :
      row_length (8)
    {}

    void reinit ()
    {
      individual_size.resize(0);
      data.resize(0);
    }

    size_type insert_new_index (const std::pair<size_type,double> &pair)
    {
      Assert(row_length > 0, ExcInternalError());
      const unsigned int index = individual_size.size();
      individual_size.push_back(1);
      data.resize(individual_size.size()*row_length);
      data[index*row_length] = pair;
      individual_size[index] = 1;
      return index;
    }

    void append_index (const size_type index,
                       const std::pair<size_type,double> &pair)
    {
      AssertIndexRange (index, individual_size.size());
      const size_type my_length = individual_size[index];
      if (my_length == row_length)
        {
          AssertDimension(data.size(), individual_size.size()*row_length);
          // no space left in this row, need to double row_length and
          // rearrange the data items. Move all items to the right except the
          // first one, starting at the back. Since individual_size contains
          // at least one element when we get here, subtracting 1 works fine.
          data.resize(2*data.size());
          for (size_type i=individual_size.size()-1; i>0; --i)
            std::memmove(&data[i*row_length*2], &data[i*row_length],
                         individual_size[i]*
                         sizeof(std::pair<size_type,double>));
          row_length *= 2;
        }
      data[index*row_length+my_length] = pair;
      individual_size[index] = my_length + 1;
    }

    size_type
    get_size (const size_type index) const
    {
      return individual_size[index];
    }

    const std::pair<size_type,double> *
    get_entry (const size_type index) const
    {
      return &data[index*row_length];
    }

    size_type row_length;

    std::vector<std::pair<size_type,double> > data;

    std::vector<size_type> individual_size;
  };



  // collects all the global rows from a local contribution (cell) and their
  // origin (direct/constraint). this is basically a vector consisting of
  // "Distributing" structs using access via the DataCache. Provides some
  // specialized sort and insert functions.
  //
  // in case there are no constraints, this is basically a list of pairs
  // <uint,unit> with the first index being the global index and the second
  // index the local index. The list is sorted with respect to the global
  // index.
  //
  // in case there are constraints, a global dof might get a contribution also
  // because it gets data from a constrained dof. This means that a global dof
  // might also have indirect contributions from a local dof via a constraint,
  // besides the direct ones.
  //
  // The actions performed here correspond to reshaping the constraint
  // information from global degrees of freedom to local ones (i.e.,
  // cell-related DoFs), and also transforming the constraint information from
  // compressed row storage (each local dof that is constrained has a list of
  // constraint entries associated to it) into compressed column storage based
  // on the cell-related DoFs (we have a list of global degrees of freedom,
  // and to each we have a list of local rows where the entries come from). To
  // increase the speed, we additionally store whether an entry is generated
  // directly from the local degrees of freedom or whether it comes from a
  // constraint.
  class GlobalRowsFromLocal
  {
  public:
    GlobalRowsFromLocal ()
      :
      n_active_rows (0),
      n_inhomogeneous_rows (0)
    {}

    void reinit (const size_type n_local_rows)
    {
      total_row_indices.resize(n_local_rows);
      for (unsigned int i=0; i<n_local_rows; ++i)
        total_row_indices[i].constraint_position = numbers::invalid_size_type;
      n_active_rows = n_local_rows;
      n_inhomogeneous_rows = 0;
      data_cache.reinit();
    }

    // implemented below
    void insert_index (const size_type global_row,
                       const size_type local_row,
                       const double       constraint_value);
    void sort ();

    // Print object for debugging purpose
    void print(std::ostream &os)
    {
      os << "Active rows " << n_active_rows << std::endl
         << "Constr rows " << n_constraints() << std::endl
         << "Inhom  rows " << n_inhomogeneous_rows << std::endl
         << "Local: ";
      for (size_type i=0 ; i<total_row_indices.size() ; ++i)
        os << ' ' << std::setw(4) << total_row_indices[i].local_row;
      os << std::endl
         << "Global:";
      for (size_type i=0 ; i<total_row_indices.size() ; ++i)
        os << ' ' << std::setw(4) << total_row_indices[i].global_row;
      os << std::endl
         << "ConPos:";
      for (size_type i=0 ; i<total_row_indices.size() ; ++i)
        os << ' ' << std::setw(4) << total_row_indices[i].constraint_position;
      os << std::endl;
    }


    // return all kind of information on the constraints

    // returns the number of global indices in the struct
    size_type size () const
    {
      return n_active_rows;
    }

    // returns the number of constraints that are associated to the
    // counter_index-th entry in the list
    size_type size (const size_type counter_index) const
    {
      return (total_row_indices[counter_index].constraint_position ==
              numbers::invalid_size_type ?
              0 :
              data_cache.get_size(total_row_indices[counter_index].
                                  constraint_position));
    }

    // returns the global row of the counter_index-th entry in the list
    size_type global_row (const size_type counter_index) const
    {
      return total_row_indices[counter_index].global_row;
    }

    // returns the global row of the counter_index-th entry in the list
    size_type &global_row (const size_type counter_index)
    {
      return total_row_indices[counter_index].global_row;
    }

    // returns the local row in the cell matrix associated with the
    // counter_index-th entry in the list. Returns invalid_size_type for
    // constrained rows
    size_type local_row (const size_type counter_index) const
    {
      return total_row_indices[counter_index].local_row;
    }

    // writable index
    size_type &local_row (const size_type counter_index)
    {
      return total_row_indices[counter_index].local_row;
    }

    // returns the local row in the cell matrix associated with the
    // counter_index-th entry in the list in the index_in_constraint-th
    // position of constraints
    size_type local_row (const size_type counter_index,
                         const size_type index_in_constraint) const
    {
      return (data_cache.get_entry(total_row_indices[counter_index].constraint_position)
              [index_in_constraint]).first;
    }

    // returns the value of the constraint in the counter_index-th entry in
    // the list in the index_in_constraint-th position of constraints
    double constraint_value (const size_type counter_index,
                             const size_type index_in_constraint) const
    {
      return (data_cache.get_entry(total_row_indices[counter_index].constraint_position)
              [index_in_constraint]).second;
    }

    // returns whether there is one row with indirect contributions (i.e.,
    // there has been at least one constraint with non-trivial ConstraintLine)
    bool have_indirect_rows () const
    {
      return data_cache.individual_size.empty() == false;
    }

    // append an entry that is constrained. This means that there is one less
    // nontrivial row
    void insert_constraint (const size_type constrained_local_dof)
    {
      --n_active_rows;
      total_row_indices[n_active_rows].local_row = constrained_local_dof;
      total_row_indices[n_active_rows].global_row = numbers::invalid_size_type;
    }

    // returns the number of constrained dofs in the structure. Constrained
    // dofs do not contribute directly to the matrix, but are needed in order
    // to set matrix diagonals and resolve inhomogeneities
    size_type n_constraints () const
    {
      return total_row_indices.size()-n_active_rows;
    }

    // returns the number of constrained dofs in the structure that have an
    // inhomogeneity
    size_type n_inhomogeneities () const
    {
      return n_inhomogeneous_rows;
    }

    // tells the structure that the ith constraint is
    // inhomogeneous. inhomogeneous constraints contribute to right hand
    // sides, so to have fast access to them, put them before homogeneous
    // constraints
    void set_ith_constraint_inhomogeneous (const size_type i)
    {
      Assert (i >= n_inhomogeneous_rows, ExcInternalError());
      std::swap (total_row_indices[n_active_rows+i],
                 total_row_indices[n_active_rows+n_inhomogeneous_rows]);
      n_inhomogeneous_rows++;
    }

    // the local row where constraint number i was detected, to find that row
    // easily when the GlobalRowsToLocal has been set up
    size_type constraint_origin (size_type i) const
    {
      return total_row_indices[n_active_rows+i].local_row;
    }

    // a vector that contains all the global ids and the corresponding local
    // ids as well as a pointer to that data where we store how to resolve
    // constraints.
    std::vector<Distributing> total_row_indices;

  private:
    // holds the actual data from the constraints
    DataCache                 data_cache;

    // how many rows there are, constraints disregarded
    size_type                 n_active_rows;

    // the number of rows with inhomogeneous constraints
    size_type                 n_inhomogeneous_rows;
  };

  // a function that appends an additional row to the list of values, or
  // appends a value to an already existing row. Similar functionality as for
  // std::map<size_type,Distributing>, but here done for a
  // std::vector<Distributing>, much faster for short lists as we have them
  // here
  inline
  void
  GlobalRowsFromLocal::insert_index (const size_type global_row,
                                     const size_type local_row,
                                     const double    constraint_value)
  {
    typedef std::vector<Distributing>::iterator index_iterator;
    index_iterator pos, pos1;
    Distributing row_value (global_row);
    std::pair<size_type,double> constraint (local_row, constraint_value);

    // check whether the list was really sorted before entering here
    for (size_type i=1; i<n_active_rows; ++i)
      Assert (total_row_indices[i-1] < total_row_indices[i], ExcInternalError());

    pos = Utilities::lower_bound (total_row_indices.begin(),
                                  total_row_indices.begin()+n_active_rows,
                                  row_value);
    if (pos->global_row == global_row)
      pos1 = pos;
    else
      {
        pos1 = total_row_indices.insert(pos, row_value);
        ++n_active_rows;
      }

    if (pos1->constraint_position == numbers::invalid_size_type)
      pos1->constraint_position = data_cache.insert_new_index (constraint);
    else
      data_cache.append_index (pos1->constraint_position, constraint);
  }

  // this sort algorithm sorts std::vector<Distributing>, but does not take
  // the constraints into account. this means that in case that constraints
  // are already inserted, this function does not work as expected. Use
  // shellsort, which is very fast in case the indices are already sorted
  // (which is the usual case with DG elements), and not too slow in other
  // cases
  inline
  void
  GlobalRowsFromLocal::sort ()
  {
    size_type i, j, j2, temp, templ, istep;
    size_type step;

    // check whether the constraints are really empty.
    const size_type length = size();

    // make sure that we are in the range of the vector
    AssertIndexRange (length, total_row_indices.size()+1);
    for (size_type i=0; i<length; ++i)
      Assert (total_row_indices[i].constraint_position ==
              numbers::invalid_size_type,
              ExcInternalError());

    step = length/2;
    while (step > 0)
      {
        for (i=step; i < length; i++)
          {
            istep = step;
            j = i;
            j2 = j-istep;
            temp = total_row_indices[i].global_row;
            templ = total_row_indices[i].local_row;
            if (total_row_indices[j2].global_row > temp)
              {
                while ((j >= istep) && (total_row_indices[j2].global_row > temp))
                  {
                    total_row_indices[j].global_row = total_row_indices[j2].global_row;
                    total_row_indices[j].local_row = total_row_indices[j2].local_row;
                    j = j2;
                    j2 -= istep;
                  }
                total_row_indices[j].global_row = temp;
                total_row_indices[j].local_row = templ;
              }
          }
        step = step>>1;
      }
  }



  /**
   * Scratch data that is used during calls to distribute_local_to_global and
   * add_entries_local_to_global. In order to avoid frequent memory
   * allocation, we keep the data alive from one call to the next in a static
   * variable. Since we want to allow for different number types in matrices,
   * this is a template.
   *
   * Since each thread gets its private version of scratch data out of the
   * ThreadLocalStorage, no conflicting access can occur. For this to be
   * valid, we need to make sure that no call within
   * distribute_local_to_global is made that by itself can spawn
   * tasks. Otherwise, we might end up in a situation where several threads
   * fight for the data.
   *
   * Access to the scratch data is only through the accessor class which
   * handles the access as well as marking the data as used.
   */
  template <typename Number>
  class ConstraintMatrixData
  {
  public:
    struct ScratchData
    {
      /**
       * Constructor, does nothing.
       */
      ScratchData ()
        :
        in_use (false)
      {}

      /**
       * Copy constructor, does nothing
       */
      ScratchData (const ScratchData &)
        :
        in_use (false)
      {}

      /**
       * Stores whether the data is currently in use.
       */
      bool in_use;

      /**
       * Temporary array for column indices
       */
      std::vector<size_type> columns;

      /**
       * Temporary array for column values
       */
      std::vector<Number>    values;

      /**
       * Temporary array for block start indices
       */
      std::vector<size_type> block_starts;

      /**
       * Temporary array for vector indices
       */
      std::vector<size_type> vector_indices;

      /**
       * Data array for reorder row/column indices. Use a shared ptr to
       * global_rows to avoid defining in the .h file
       */
      GlobalRowsFromLocal global_rows;

      /**
       * Data array for reorder row/column indices. Use a shared ptr to
       * global_rows to avoid defining in the .h file
       */
      GlobalRowsFromLocal global_columns;
    };

    /**
     * Accessor class to guard access to scratch_data
     */
    class ScratchDataAccessor
    {
    public:
      /**
       * Constructor. Grabs a scratch data object on the current thread and
       * mark it as used
       */
      ScratchDataAccessor()
        :
        my_scratch_data(&ConstraintMatrixData::scratch_data.get())
      {
        Assert(my_scratch_data->in_use == false,
               ExcMessage("Access to thread-local scratch data tried, but it is already "
                          "in use"));
        my_scratch_data->in_use = true;
      }

      /**
       * Destructor. Mark scratch data as available again.
       */
      ~ScratchDataAccessor()
      {
        my_scratch_data->in_use = false;
      }

      /**
       * Dereferencing operator.
       */
      ScratchData &operator* ()
      {
        return *my_scratch_data;
      }

      /**
       * Dereferencing operator.
       */
      ScratchData *operator-> ()
      {
        return my_scratch_data;
      }

    private:
      ScratchData *my_scratch_data;
    };

  private:
    /**
     * The actual data object that contains a scratch data for each thread.
     */
    static Threads::ThreadLocalStorage<ScratchData> scratch_data;
  };



  // function for block matrices: Find out where in the list of local dofs
  // (sorted according to global ids) the individual blocks start. Transform
  // the global indices to block-local indices in order to be able to use
  // functions like vector.block(1)(block_local_id), instead of
  // vector(global_id). This avoids transforming indices one-by-one later on.
  template <class BlockType>
  inline
  void
  make_block_starts (const BlockType        &block_object,
                     GlobalRowsFromLocal    &global_rows,
                     std::vector<size_type> &block_starts)
  {
    AssertDimension (block_starts.size(), block_object.n_block_rows()+1);

    typedef std::vector<Distributing>::iterator row_iterator;
    row_iterator block_indices = global_rows.total_row_indices.begin();

    const size_type num_blocks = block_object.n_block_rows();
    const size_type n_active_rows = global_rows.size();

    // find end of rows.
    block_starts[0] = 0;
    for (size_type i=1; i<num_blocks; ++i)
      {
        row_iterator first_block =
          Utilities::lower_bound (block_indices,
                                  global_rows.total_row_indices.begin()+n_active_rows,
                                  Distributing(block_object.get_row_indices().block_start(i)));
        block_starts[i] = first_block - global_rows.total_row_indices.begin();
        block_indices = first_block;
      }
    block_starts[num_blocks] = n_active_rows;

    // transform row indices to block-local index space
    for (size_type i=block_starts[1]; i<n_active_rows; ++i)
      global_rows.global_row(i) = block_object.get_row_indices().
                                  global_to_local(global_rows.global_row(i)).second;
  }



  // same as before, but for std::vector<uint> instead of
  // GlobalRowsFromLocal. Used in functions for sparsity patterns.
  template <class BlockType>
  inline
  void
  make_block_starts (const BlockType        &block_object,
                     std::vector<size_type> &row_indices,
                     std::vector<size_type> &block_starts)
  {
    AssertDimension (block_starts.size(), block_object.n_block_rows()+1);

    typedef std::vector<size_type>::iterator row_iterator;
    row_iterator col_indices = row_indices.begin();

    const size_type num_blocks = block_object.n_block_rows();

    // find end of rows.
    block_starts[0] = 0;
    for (size_type i=1; i<num_blocks; ++i)
      {
        row_iterator first_block =
          Utilities::lower_bound (col_indices,
                                  row_indices.end(),
                                  block_object.get_row_indices().block_start(i));
        block_starts[i] = first_block - row_indices.begin();
        col_indices = first_block;
      }
    block_starts[num_blocks] = row_indices.size();

    // transform row indices to local index space
    for (size_type i=block_starts[1]; i<row_indices.size(); ++i)
      row_indices[i] = block_object.get_row_indices().
                       global_to_local(row_indices[i]).second;
  }



  // resolves constraints of one column at the innermost loop. goes through
  // the origin of each global entry and finds out which data we need to
  // collect.
  static inline
  double resolve_matrix_entry (const GlobalRowsFromLocal &global_rows,
                               const GlobalRowsFromLocal &global_cols,
                               const size_type            i,
                               const size_type            j,
                               const size_type            loc_row,
                               const FullMatrix<double> &local_matrix)
  {
    const size_type loc_col = global_cols.local_row(j);
    double col_val;

    // case 1: row has direct contribution in local matrix. decide whether col
    // has a direct contribution. if not, set the value to zero.
    if (loc_row != numbers::invalid_size_type)
      {
        col_val = ((loc_col != numbers::invalid_size_type) ?
                   local_matrix(loc_row, loc_col) : 0);

        // account for indirect contributions by constraints in column
        for (size_type p=0; p<global_cols.size(j); ++p)
          col_val += (local_matrix(loc_row, global_cols.local_row(j,p)) *
                      global_cols.constraint_value(j,p));
      }

    // case 2: row has no direct contribution in local matrix
    else
      col_val = 0;

    // account for indirect contributions by constraints in row, going trough
    // the direct and indirect references in the given column.
    for (size_type q=0; q<global_rows.size(i); ++q)
      {
        double add_this = (loc_col != numbers::invalid_size_type)
                          ? local_matrix(global_rows.local_row(i,q), loc_col) : 0;

        for (size_type p=0; p<global_cols.size(j); ++p)
          add_this += (local_matrix(global_rows.local_row(i,q),
                                    global_cols.local_row(j,p))
                       *
                       global_cols.constraint_value(j,p));
        col_val += add_this * global_rows.constraint_value(i,q);
      }
    return col_val;
  }



  // computes all entries that need to be written into global_rows[i]. Lists
  // the resulting values in val_ptr, and the corresponding column indices in
  // col_ptr.
  template <typename number>
  inline
  void
  resolve_matrix_row (const GlobalRowsFromLocal &global_rows,
                      const GlobalRowsFromLocal &global_cols,
                      const size_type            i,
                      const size_type            column_start,
                      const size_type            column_end,
                      const FullMatrix<double>  &local_matrix,
                      size_type                *&col_ptr,
                      number                   *&val_ptr)
  {
    if (column_end == column_start)
      return;

    AssertIndexRange (column_end-1, global_cols.size());
    const size_type loc_row = global_rows.local_row(i);

    // fast function if there are no indirect references to any of the local
    // rows at all on this set of dofs (saves a lot of checks). the only check
    // we actually need to perform is whether the matrix element is zero.
    if (global_rows.have_indirect_rows() == false &&
        global_cols.have_indirect_rows() == false)
      {
        AssertIndexRange(loc_row, local_matrix.m());
        const double *matrix_ptr = &local_matrix(loc_row, 0);

        for (size_type j=column_start; j<column_end; ++j)
          {
            const size_type loc_col = global_cols.local_row(j);
            AssertIndexRange(loc_col, local_matrix.n());
            const double col_val = matrix_ptr[loc_col];
            if (col_val != 0.)
              {
                *val_ptr++ = static_cast<number> (col_val);
                *col_ptr++ = global_cols.global_row(j);
              }
          }
      }

    // more difficult part when there are indirect references and when we need
    // to do some more checks.
    else
      {
        for (size_type j=column_start; j<column_end; ++j)
          {
            double col_val = resolve_matrix_entry (global_rows, global_cols, i, j,
                                                   loc_row, local_matrix);

            // if we got some nontrivial value, append it to the array of
            // values.
            if (col_val != 0.)
              {
                *val_ptr++ = static_cast<number> (col_val);
                *col_ptr++ = global_cols.global_row(j);
              }
          }
      }
  }



  // specialized function that can write into the row of a
  // SparseMatrix<number>.
  namespace dealiiSparseMatrix
  {
    template <typename SparseMatrixIterator>
    static inline
    void add_value (const double          value,
                    const size_type       row,
                    const size_type       column,
                    SparseMatrixIterator &matrix_values)
    {
      if (value != 0.)
        {
          while (matrix_values->column() < column)
            ++matrix_values;
          Assert (matrix_values->column() == column,
                  typename SparseMatrix<typename SparseMatrixIterator::MatrixType::value_type>::ExcInvalidIndex(row, column));
          matrix_values->value() += value;
        }
    }
  }


  // similar as before, now with shortcut for deal.II sparse matrices. this
  // lets us avoid using extra arrays, and does all the operations just in
  // place, i.e., in the respective matrix row
  template <typename number>
  inline
  void
  resolve_matrix_row (const GlobalRowsFromLocal &global_rows,
                      const size_type            i,
                      const size_type            column_start,
                      const size_type            column_end,
                      const FullMatrix<double>  &local_matrix,
                      SparseMatrix<number>      *sparse_matrix)
  {
    if (column_end == column_start)
      return;

    AssertIndexRange (column_end-1, global_rows.size());
    const SparsityPattern &sparsity = sparse_matrix->get_sparsity_pattern();

    if (sparsity.n_nonzero_elements() == 0)
      return;

    const size_type row = global_rows.global_row(i);
    const size_type loc_row = global_rows.local_row(i);

    typename SparseMatrix<number>::iterator
    matrix_values = sparse_matrix->begin(row);
    const bool optimize_diagonal = sparsity.n_rows() == sparsity.n_cols();

    // distinguish three cases about what can happen for checking whether the
    // diagonal is the first element of the row. this avoids if statements at
    // the innermost loop positions

    if (!optimize_diagonal) // case 1: no diagonal optimization in matrix
      {
        if (global_rows.have_indirect_rows() == false)
          {
            AssertIndexRange (loc_row, local_matrix.m());
            const double *matrix_ptr = &local_matrix(loc_row, 0);

            for (size_type j=column_start; j<column_end; ++j)
              {
                const size_type loc_col = global_rows.local_row(j);
                const double col_val = matrix_ptr[loc_col];
                dealiiSparseMatrix::add_value (col_val, row,
                                               global_rows.global_row(j),
                                               matrix_values);
              }
          }
        else
          {
            for (size_type j=column_start; j<column_end; ++j)
              {
                double col_val = resolve_matrix_entry (global_rows, global_rows, i, j,
                                                       loc_row, local_matrix);
                dealiiSparseMatrix::add_value (col_val, row,
                                               global_rows.global_row(j),
                                               matrix_values);
              }
          }
      }
    else if (i>=column_start && i<column_end) // case 2: can split loop
      {
        ++matrix_values; // jump over diagonal element
        if (global_rows.have_indirect_rows() == false)
          {
            AssertIndexRange (loc_row, local_matrix.m());
            const double *matrix_ptr = &local_matrix(loc_row, 0);

            sparse_matrix->begin(row)->value() += matrix_ptr[loc_row];
            for (size_type j=column_start; j<i; ++j)
              {
                const size_type loc_col = global_rows.local_row(j);
                const double col_val = matrix_ptr[loc_col];
                dealiiSparseMatrix::add_value(col_val, row,
                                              global_rows.global_row(j),
                                              matrix_values);
              }
            for (size_type j=i+1; j<column_end; ++j)
              {
                const size_type loc_col = global_rows.local_row(j);
                const double col_val = matrix_ptr[loc_col];
                dealiiSparseMatrix::add_value(col_val, row,
                                              global_rows.global_row(j),
                                              matrix_values);
              }
          }
        else
          {
            sparse_matrix->begin(row)->value() +=
              resolve_matrix_entry (global_rows, global_rows, i, i,
                                    loc_row, local_matrix);
            for (size_type j=column_start; j<i; ++j)
              {
                double col_val = resolve_matrix_entry (global_rows, global_rows, i, j,
                                                       loc_row, local_matrix);
                dealiiSparseMatrix::add_value (col_val, row,
                                               global_rows.global_row(j),
                                               matrix_values);
              }
            for (size_type j=i+1; j<column_end; ++j)
              {
                double col_val = resolve_matrix_entry (global_rows, global_rows, i, j,
                                                       loc_row, local_matrix);
                dealiiSparseMatrix::add_value (col_val, row,
                                               global_rows.global_row(j),
                                               matrix_values);
              }
          }
      }
    // case 3: can't say - need to check inside the loop
    else if (global_rows.have_indirect_rows() == false)
      {
        ++matrix_values; // jump over diagonal element
        AssertIndexRange (loc_row, local_matrix.m());
        const double *matrix_ptr = &local_matrix(loc_row, 0);

        for (size_type j=column_start; j<column_end; ++j)
          {
            const size_type loc_col = global_rows.local_row(j);
            const double col_val = matrix_ptr[loc_col];
            if (row==global_rows.global_row(j))
              sparse_matrix->begin(row)->value() += col_val;
            else
              dealiiSparseMatrix::add_value(col_val, row,
                                            global_rows.global_row(j),
                                            matrix_values);
          }
      }
    else
      {
        ++matrix_values; // jump over diagonal element
        for (size_type j=column_start; j<column_end; ++j)
          {
            double col_val = resolve_matrix_entry (global_rows, global_rows, i,
                                                   j, loc_row, local_matrix);
            if (row==global_rows.global_row(j))
              sparse_matrix->begin(row)->value() += col_val;
            else
              dealiiSparseMatrix::add_value (col_val, row,
                                             global_rows.global_row(j),
                                             matrix_values);
          }
      }
  }



  // Same function to resolve all entries that will be added to the given
  // global row global_rows[i] as before, now for sparsity pattern
  inline
  void
  resolve_matrix_row (const GlobalRowsFromLocal     &global_rows,
                      const size_type                i,
                      const size_type                column_start,
                      const size_type                column_end,
                      const Table<2,bool>           &dof_mask,
                      std::vector<size_type>::iterator &col_ptr)
  {
    if (column_end == column_start)
      return;

    const size_type loc_row = global_rows.local_row(i);

    // fast function if there are no indirect references to any of the local
    // rows at all on this set of dofs
    if (global_rows.have_indirect_rows() == false)
      {
        Assert(loc_row < dof_mask.n_rows(),
               ExcInternalError());

        for (size_type j=column_start; j<column_end; ++j)
          {
            const size_type loc_col = global_rows.local_row(j);
            Assert(loc_col < dof_mask.n_cols(), ExcInternalError());

            if (dof_mask(loc_row,loc_col) == true)
              *col_ptr++ = global_rows.global_row(j);
          }
      }

    // slower functions when there are indirect references and when we need to
    // do some more checks.
    else
      {
        for (size_type j=column_start; j<column_end; ++j)
          {
            const size_type loc_col = global_rows.local_row(j);
            if (loc_row != numbers::invalid_size_type)
              {
                Assert (loc_row < dof_mask.n_rows(), ExcInternalError());
                if (loc_col != numbers::invalid_size_type)
                  {
                    Assert (loc_col < dof_mask.n_cols(), ExcInternalError());
                    if (dof_mask(loc_row,loc_col) == true)
                      goto add_this_index;
                  }

                for (size_type p=0; p<global_rows.size(j); ++p)
                  if (dof_mask(loc_row,global_rows.local_row(j,p)) == true)
                    goto add_this_index;
              }

            for (size_type q=0; q<global_rows.size(i); ++q)
              {
                if (loc_col != numbers::invalid_size_type)
                  {
                    Assert (loc_col < dof_mask.n_cols(), ExcInternalError());
                    if (dof_mask(global_rows.local_row(i,q),loc_col) == true)
                      goto add_this_index;
                  }

                for (size_type p=0; p<global_rows.size(j); ++p)
                  if (dof_mask(global_rows.local_row(i,q),
                               global_rows.local_row(j,p)) == true)
                    goto add_this_index;
              }

            continue;
            // if we got some nontrivial value, append it to the array of
            // values.
add_this_index:
            *col_ptr++ = global_rows.global_row(j);
          }
      }
  }



  // to make sure that the global matrix remains invertible, we need to do
  // something with the diagonal elements. add the absolute value of the local
  // matrix, so the resulting entry will always be positive and furthermore be
  // in the same order of magnitude as the other elements of the matrix
  //
  // note that this also captures the special case that a dof is both
  // constrained and fixed (this can happen for hanging nodes in 3d that also
  // happen to be on the boundary). in that case, following the program flow
  // in distribute_local_to_global, it is realized that when distributing the
  // row and column no elements of the matrix are actually touched if all the
  // degrees of freedom to which this dof is constrained are also constrained
  // (the usual case with hanging nodes in 3d). however, in the line below, we
  // do actually do something with this dof
  template <typename MatrixType, typename VectorType>
  inline void
  set_matrix_diagonals (const internals::GlobalRowsFromLocal &global_rows,
                        const std::vector<size_type>         &local_dof_indices,
                        const FullMatrix<double>             &local_matrix,
                        const ConstraintMatrix               &constraints,
                        MatrixType                           &global_matrix,
                        VectorType                           &global_vector,
                        bool                                 use_inhomogeneities_for_rhs)
  {
    if (global_rows.n_constraints() > 0)
      {
        double average_diagonal = 0;
        for (size_type i=0; i<local_matrix.m(); ++i)
          average_diagonal += std::fabs (local_matrix(i,i));
        average_diagonal /= static_cast<double>(local_matrix.m());

        for (size_type i=0; i<global_rows.n_constraints(); i++)
          {
            const size_type local_row = global_rows.constraint_origin(i);
            const size_type global_row = local_dof_indices[local_row];
            const typename MatrixType::value_type new_diagonal
              = (std::fabs(local_matrix(local_row,local_row)) != 0 ?
                 std::fabs(local_matrix(local_row,local_row)) : average_diagonal);
            global_matrix.add(global_row, global_row, new_diagonal);

            // if the use_inhomogeneities_for_rhs flag is set to true, the
            // inhomogeneities are used to create the global vector. instead
            // of fill in a zero in the ith components with an inhomogeneity,
            // we set those to: inhomogeneity(i)*global_matrix (i,i).
            if (use_inhomogeneities_for_rhs == true)
              global_vector(global_row) += constraints.get_inhomogeneity(global_row) * new_diagonal;
          }
      }
  }



  // similar function as the one above for setting matrix diagonals, but now
  // doing that for sparsity patterns when setting them up using
  // add_entries_local_to_global. In case we keep constrained entries, add all
  // the rows and columns related to the constrained dof, otherwise just add
  // the diagonal
  template <typename SparsityType>
  inline void
  set_sparsity_diagonals (const internals::GlobalRowsFromLocal &global_rows,
                          const std::vector<size_type>         &local_dof_indices,
                          const Table<2,bool>                  &dof_mask,
                          const bool                            keep_constrained_entries,
                          SparsityType                         &sparsity_pattern)
  {
    // if we got constraints, need to add the diagonal element and, if the
    // user requested so, also the rest of the entries in rows and columns
    // that have been left out above
    if (global_rows.n_constraints() > 0)
      {
        for (size_type i=0; i<global_rows.n_constraints(); i++)
          {
            const size_type local_row = global_rows.constraint_origin(i);
            const size_type global_row = local_dof_indices[local_row];
            if (keep_constrained_entries == true)
              {
                for (size_type j=0; j<local_dof_indices.size(); ++j)
                  {
                    if (dof_mask(local_row,j) == true)
                      sparsity_pattern.add(global_row,
                                           local_dof_indices[j]);
                    if (dof_mask(j,local_row) == true)
                      sparsity_pattern.add(local_dof_indices[j],
                                           global_row);
                  }
              }
            else
              // don't keep constrained entries - just add the diagonal.
              sparsity_pattern.add(global_row,global_row);
          }
      }
  }

} // end of namespace internals



// Basic idea of setting up a list of
// all global dofs: first find all rows and columns
// that we are going to write touch,
// and then go through the
// lines and collect all the local rows that
// are related to it.
void
ConstraintMatrix::
make_sorted_row_list (const std::vector<size_type>   &local_dof_indices,
                      internals::GlobalRowsFromLocal &global_rows) const
{
  const size_type n_local_dofs = local_dof_indices.size();
  AssertDimension (n_local_dofs, global_rows.size());

  // when distributing the local data to the global matrix, we can quite
  // cheaply sort the indices (obviously, this introduces the need for
  // allocating some memory on the way, but we need to do this only for rows,
  // whereas the distribution process itself goes over rows and columns). This
  // has the advantage that when writing into the global matrix, we can make
  // use of the sortedness.

  // so the first step is to create a sorted list of all row values that are
  // possible. these values are either the rows from unconstrained dofs, or
  // some indices introduced by dofs constrained to a combination of some
  // other dofs. regarding the data type, choose an STL vector of a pair of
  // unsigned ints (for global columns) and internal data (containing local
  // columns + possible jumps from constraints). Choosing an STL map or
  // anything else M.K. knows of would be much more expensive here!

  // cache whether we have to resolve any indirect rows generated from
  // resolving constrained dofs.
  size_type added_rows = 0;

  // first add the indices in an unsorted way and only keep track of the
  // constraints that appear. They are resolved in a second step.
  for (size_type i = 0; i<n_local_dofs; ++i)
    {
      if (is_constrained(local_dof_indices[i]) == false)
        {
          global_rows.global_row(added_rows)  = local_dof_indices[i];
          global_rows.local_row(added_rows++) = i;
        }
      else
        global_rows.insert_constraint(i);
    }
  global_rows.sort();

  const size_type n_constrained_rows = n_local_dofs-added_rows;
  for (size_type i=0; i<n_constrained_rows; ++i)
    {
      const size_type local_row = global_rows.constraint_origin(i);
      AssertIndexRange(local_row, n_local_dofs);
      const size_type global_row = local_dof_indices[local_row];
      Assert (is_constrained(global_row), ExcInternalError());
      const ConstraintLine &position =
        lines[lines_cache[calculate_line_index(global_row)]];
      if (position.inhomogeneity != 0)
        global_rows.set_ith_constraint_inhomogeneous (i);
      for (size_type q=0; q<position.entries.size(); ++q)
        global_rows.insert_index (position.entries[q].first,
                                  local_row,
                                  position.entries[q].second);
    }
}



// Same function as before, but now do only extract the global indices that
// come from the local ones without storing their origin. Used for sparsity
// pattern generation.
inline
void
ConstraintMatrix::
make_sorted_row_list (const std::vector<size_type> &local_dof_indices,
                      std::vector<size_type>       &active_dofs) const
{
  const size_type n_local_dofs = local_dof_indices.size();
  size_type added_rows = 0;
  for (size_type i = 0; i<n_local_dofs; ++i)
    {
      if (is_constrained(local_dof_indices[i]) == false)
        {
          active_dofs[added_rows++] = local_dof_indices[i];
          continue;
        }

      active_dofs[n_local_dofs-i+added_rows-1] = i;
    }
  std::sort (active_dofs.begin(), active_dofs.begin()+added_rows);

  const size_type n_constrained_dofs = n_local_dofs-added_rows;
  for (size_type i=n_constrained_dofs; i>0; --i)
    {
      const size_type local_row = active_dofs.back();

      // remove constrained entry since we are going to resolve it in place
      active_dofs.pop_back();
      const size_type global_row = local_dof_indices[local_row];
      const ConstraintLine &position =
        lines[lines_cache[calculate_line_index(global_row)]];
      for (size_type q=0; q<position.entries.size(); ++q)
        {
          const size_type new_index = position.entries[q].first;
          if (active_dofs[active_dofs.size()-i] < new_index)
            active_dofs.insert(active_dofs.end()-i+1,new_index);

          // make binary search to find where to put the new index in order to
          // keep the list sorted
          else
            {
              std::vector<size_type>::iterator it =
                Utilities::lower_bound(active_dofs.begin(),
                                       active_dofs.end()-i+1,
                                       new_index);
              if (*it != new_index)
                active_dofs.insert(it, new_index);
            }
        }
    }
}



// Resolve the constraints from the vector and apply inhomogeneities.
inline
double
ConstraintMatrix::
resolve_vector_entry (const size_type                       i,
                      const internals::GlobalRowsFromLocal &global_rows,
                      const Vector<double>                 &local_vector,
                      const std::vector<size_type>         &local_dof_indices,
                      const FullMatrix<double>             &local_matrix) const
{
  const size_type loc_row = global_rows.local_row(i);
  const size_type n_inhomogeneous_rows = global_rows.n_inhomogeneities();
  double val = 0;
  // has a direct contribution from some local entry. If we have inhomogeneous
  // constraints, compute the contribution of the inhomogeneity in the current
  // row.
  if (loc_row != numbers::invalid_size_type)
    {
      val = local_vector(loc_row);
      for (size_type i=0; i<n_inhomogeneous_rows; ++i)
        val -= (lines[lines_cache[calculate_line_index(local_dof_indices
                                                       [global_rows.constraint_origin(i)])]].
                inhomogeneity *
                local_matrix(loc_row, global_rows.constraint_origin(i)));
    }

  // go through the indirect contributions
  for (size_type q=0; q<global_rows.size(i); ++q)
    {
      const size_type loc_row_q = global_rows.local_row(i,q);
      double add_this = local_vector (loc_row_q);
      for (size_type k=0; k<n_inhomogeneous_rows; ++k)
        add_this -= (lines[lines_cache[calculate_line_index
                                       (local_dof_indices
                                        [global_rows.constraint_origin(k)])]].
                     inhomogeneity *
                     local_matrix(loc_row_q,global_rows.constraint_origin(k)));
      val += add_this * global_rows.constraint_value(i,q);
    }
  return val;
}


// internal implementation for distribute_local_to_global for standard
// (non-block) matrices
template <typename MatrixType, typename VectorType>
void
ConstraintMatrix::distribute_local_to_global (
  const FullMatrix<double>        &local_matrix,
  const Vector<double>            &local_vector,
  const std::vector<size_type>    &local_dof_indices,
  MatrixType                      &global_matrix,
  VectorType                      &global_vector,
  bool                            use_inhomogeneities_for_rhs,
  internal::bool2type<false>) const
{
  // check whether we work on real vectors or we just used a dummy when
  // calling the other function above.
  const bool use_vectors = (local_vector.size() == 0 &&
                            global_vector.size() == 0) ? false : true;
  typedef typename MatrixType::value_type number;
  const bool use_dealii_matrix =
    types_are_equal<MatrixType,SparseMatrix<number> >::value;

  AssertDimension (local_matrix.n(), local_dof_indices.size());
  AssertDimension (local_matrix.m(), local_dof_indices.size());
  Assert (global_matrix.m() == global_matrix.n(), ExcNotQuadratic());
  if (use_vectors == true)
    {
      AssertDimension (local_matrix.m(), local_vector.size());
      AssertDimension (global_matrix.m(), global_vector.size());
    }
  Assert (lines.empty() || sorted == true, ExcMatrixNotClosed());

  const size_type n_local_dofs = local_dof_indices.size();

  typename internals::ConstraintMatrixData<number>::ScratchDataAccessor
  scratch_data;

  internals::GlobalRowsFromLocal &global_rows = scratch_data->global_rows;
  global_rows.reinit(n_local_dofs);
  make_sorted_row_list (local_dof_indices, global_rows);

  const size_type n_actual_dofs = global_rows.size();

  // create arrays for the column data (indices and values) that will then be
  // written into the matrix. Shortcut for deal.II sparse matrix. We can use
  // the scratch data if we have a double matrix. Otherwise, we need to create
  // an array in any case since we cannot know about the actual data type in
  // the ConstraintMatrix class (unless we do cast). This involves a little
  // bit of logic to determine the type of the matrix value.
  std::vector<size_type> &cols = scratch_data->columns;
  std::vector<number>     &vals = scratch_data->values;
  SparseMatrix<number> *sparse_matrix
    = dynamic_cast<SparseMatrix<number> *>(&global_matrix);
  if (use_dealii_matrix == false)
    {
      cols.resize (n_actual_dofs);
      vals.resize (n_actual_dofs);
    }
  else
    Assert (sparse_matrix != 0, ExcInternalError());

  // now do the actual job. go through all the global rows that we will touch
  // and call resolve_matrix_row for each of those.
  for (size_type i=0; i<n_actual_dofs; ++i)
    {
      const size_type row = global_rows.global_row(i);

      // calculate all the data that will be written into the matrix row.
      if (use_dealii_matrix == false)
        {
          size_type *col_ptr = &cols[0];
          // cast is uncritical here and only used to avoid compiler
          // warnings. We never access a non-double array
          number *val_ptr = &vals[0];
          internals::resolve_matrix_row (global_rows, global_rows, i, 0,
                                         n_actual_dofs,
                                         local_matrix, col_ptr, val_ptr);
          const size_type n_values = col_ptr - &cols[0];
          if (n_values > 0)
            global_matrix.add(row, n_values, &cols[0], &vals[0], false,
                              true);
        }
      else
        internals::resolve_matrix_row (global_rows, i, 0, n_actual_dofs,
                                       local_matrix, sparse_matrix);

      // now to the vectors. besides doing the same job as we did above (i.e.,
      // distribute the content of the local vector into the global one), need
      // to account for inhomogeneities here: thie corresponds to eliminating
      // the respective column in the local matrix with value on the right
      // hand side.
      if (use_vectors == true)
        {
          const double val = resolve_vector_entry (i, global_rows,
                                                   local_vector,
                                                   local_dof_indices,
                                                   local_matrix);

          if (val != 0)
            global_vector(row) += static_cast<typename VectorType::value_type>(val);
        }
    }

  internals::set_matrix_diagonals (global_rows, local_dof_indices,
                                   local_matrix, *this,
                                   global_matrix, global_vector, use_inhomogeneities_for_rhs);
}



template <typename MatrixType>
void
ConstraintMatrix::distribute_local_to_global (
  const FullMatrix<double>     &local_matrix,
  const std::vector<size_type> &row_indices,
  const std::vector<size_type> &col_indices,
  MatrixType                   &global_matrix) const
{
  typedef double number;

  AssertDimension (local_matrix.m(), row_indices.size());
  AssertDimension (local_matrix.n(), col_indices.size());
  //Assert (sorted == true, ExcMatrixNotClosed());

  const size_type n_local_row_dofs = row_indices.size();
  const size_type n_local_col_dofs = col_indices.size();

  typename internals::ConstraintMatrixData<number>::ScratchDataAccessor
  scratch_data;
  internals::GlobalRowsFromLocal &global_rows = scratch_data->global_rows;
  global_rows.reinit(n_local_row_dofs);
  internals::GlobalRowsFromLocal &global_cols = scratch_data->global_columns;
  global_cols.reinit(n_local_col_dofs);
  make_sorted_row_list (row_indices, global_rows);
  make_sorted_row_list (col_indices, global_cols);

  const size_type n_actual_row_dofs = global_rows.size();
  const size_type n_actual_col_dofs = global_cols.size();

  // create arrays for the column data (indices and values) that will then be
  // written into the matrix. Shortcut for deal.II sparse matrix
  std::vector<size_type> &cols = scratch_data->columns;
  std::vector<number>     &vals = scratch_data->values;
  cols.resize(n_actual_col_dofs);
  vals.resize(n_actual_col_dofs);

  // now do the actual job.
  for (size_type i=0; i<n_actual_row_dofs; ++i)
    {
      const size_type row = global_rows.global_row(i);

      // calculate all the data that will be written into the matrix row.
      size_type *col_ptr = &cols[0];
      number    *val_ptr = &vals[0];
      internals::resolve_matrix_row (global_rows, global_cols, i, 0,
                                     n_actual_col_dofs,
                                     local_matrix, col_ptr, val_ptr);
      const size_type n_values = col_ptr - &cols[0];
      if (n_values > 0)
        global_matrix.add(row, n_values, &cols[0], &vals[0], false, true);
    }
}


// similar function as above, but now specialized for block matrices. See the
// other function for additional comments.
template <typename MatrixType, typename VectorType>
void
ConstraintMatrix::
distribute_local_to_global (const FullMatrix<double>     &local_matrix,
                            const Vector<double>         &local_vector,
                            const std::vector<size_type> &local_dof_indices,
                            MatrixType                   &global_matrix,
                            VectorType                   &global_vector,
                            bool                          use_inhomogeneities_for_rhs,
                            internal::bool2type<true>) const
{
  const bool use_vectors = (local_vector.size() == 0 &&
                            global_vector.size() == 0) ? false : true;
  typedef typename MatrixType::value_type number;
  const bool use_dealii_matrix =
    types_are_equal<MatrixType,BlockSparseMatrix<number> >::value;

  AssertDimension (local_matrix.n(), local_dof_indices.size());
  AssertDimension (local_matrix.m(), local_dof_indices.size());
  Assert (global_matrix.m() == global_matrix.n(), ExcNotQuadratic());
  Assert (global_matrix.n_block_rows() == global_matrix.n_block_cols(),
          ExcNotQuadratic());
  if (use_vectors == true)
    {
      AssertDimension (local_matrix.m(), local_vector.size());
      AssertDimension (global_matrix.m(), global_vector.size());
    }
  Assert (sorted == true, ExcMatrixNotClosed());

  typename internals::ConstraintMatrixData<number>::ScratchDataAccessor
  scratch_data;

  const size_type n_local_dofs = local_dof_indices.size();
  internals::GlobalRowsFromLocal &global_rows = scratch_data->global_rows;
  global_rows.reinit(n_local_dofs);

  make_sorted_row_list (local_dof_indices, global_rows);
  const size_type n_actual_dofs = global_rows.size();

  std::vector<size_type> &global_indices = scratch_data->vector_indices;
  if (use_vectors == true)
    {
      global_indices.resize(n_actual_dofs);
      for (size_type i=0; i<n_actual_dofs; ++i)
        global_indices[i] = global_rows.global_row(i);
    }

  // additional construct that also takes care of block indices.
  const size_type num_blocks   = global_matrix.n_block_rows();
  std::vector<size_type> &block_starts = scratch_data->block_starts;
  block_starts.resize(num_blocks+1);
  internals::make_block_starts (global_matrix, global_rows, block_starts);

  std::vector<size_type> &cols = scratch_data->columns;
  std::vector<number>     &vals = scratch_data->values;
  if (use_dealii_matrix == false)
    {
      cols.resize (n_actual_dofs);
      vals.resize (n_actual_dofs);
    }

  // the basic difference to the non-block variant from now onwards is that we
  // go through the blocks of the matrix separately, which allows us to set
  // the block entries individually
  for (size_type block=0; block<num_blocks; ++block)
    {
      const size_type next_block = block_starts[block+1];
      for (size_type i=block_starts[block]; i<next_block; ++i)
        {
          const size_type row = global_rows.global_row(i);

          for (size_type block_col=0; block_col<num_blocks; ++block_col)
            {
              const size_type start_block = block_starts[block_col],
                              end_block = block_starts[block_col+1];
              if (use_dealii_matrix == false)
                {
                  size_type *col_ptr = &cols[0];
                  number *val_ptr = &vals[0];
                  internals::resolve_matrix_row (global_rows, global_rows, i,
                                                 start_block, end_block,
                                                 local_matrix, col_ptr, val_ptr);
                  const size_type n_values = col_ptr - &cols[0];
                  if (n_values > 0)
                    global_matrix.block(block, block_col).add(row, n_values,
                                                              &cols[0], &vals[0],
                                                              false, true);
                }
              else
                {
                  SparseMatrix<number> *sparse_matrix
                    = dynamic_cast<SparseMatrix<number> *>(&global_matrix.block(block,
                                                           block_col));
                  Assert (sparse_matrix != 0, ExcInternalError());
                  internals::resolve_matrix_row (global_rows, i, start_block,
                                                 end_block, local_matrix, sparse_matrix);
                }
            }

          if (use_vectors == true)
            {
              const double val = resolve_vector_entry (i, global_rows,
                                                       local_vector,
                                                       local_dof_indices,
                                                       local_matrix);

              if (val != 0)
                global_vector(global_indices[i]) +=
                  static_cast<typename VectorType::value_type>(val);
            }
        }
    }

  internals::set_matrix_diagonals (global_rows, local_dof_indices,
                                   local_matrix, *this,
                                   global_matrix, global_vector, use_inhomogeneities_for_rhs);
}



template <typename SparsityType>
void
ConstraintMatrix::
add_entries_local_to_global (const std::vector<size_type> &local_dof_indices,
                             SparsityType                 &sparsity_pattern,
                             const bool                    keep_constrained_entries,
                             const Table<2,bool>          &dof_mask,
                             internal::bool2type<false> ) const
{
  Assert (sparsity_pattern.n_rows() == sparsity_pattern.n_cols(), ExcNotQuadratic());

  const size_type n_local_dofs = local_dof_indices.size();
  bool dof_mask_is_active = false;
  if (dof_mask.n_rows() == n_local_dofs)
    {
      dof_mask_is_active = true;
      AssertDimension (dof_mask.n_cols(), n_local_dofs);
    }

  internals::ConstraintMatrixData<double>::ScratchDataAccessor scratch_data;

  // if the dof mask is not active, all we have to do is to add some indices
  // in a matrix format. To do this, we first create an array of all the
  // indices that are to be added. these indices are the local dof indices
  // plus some indices that come from constraints.
  if (dof_mask_is_active == false)
    {
      std::vector<size_type> &actual_dof_indices = scratch_data->columns;
      actual_dof_indices.resize(n_local_dofs);
      make_sorted_row_list (local_dof_indices, actual_dof_indices);
      const size_type n_actual_dofs = actual_dof_indices.size();

      // now add the indices we collected above to the sparsity pattern. Very
      // easy here - just add the same array to all the rows...
      for (size_type i=0; i<n_actual_dofs; ++i)
        sparsity_pattern.add_entries(actual_dof_indices[i],
                                     actual_dof_indices.begin(),
                                     actual_dof_indices.end(),
                                     true);

      // need to add the whole row and column structure in case we keep
      // constrained entries. Unfortunately, we can't use the nice matrix
      // structure we use elsewhere, so manually add those indices one by one.
      for (size_type i=0; i<n_local_dofs; i++)
        if (is_constrained(local_dof_indices[i]))
          {
            if (keep_constrained_entries == true)
              for (size_type j=0; j<n_local_dofs; j++)
                {
                  sparsity_pattern.add (local_dof_indices[i], local_dof_indices[j]);
                  sparsity_pattern.add (local_dof_indices[j], local_dof_indices[i]);
                }
            else
              sparsity_pattern.add (local_dof_indices[i], local_dof_indices[i]);
          }

      return;
    }


  // complicated case: we need to filter out some indices. then the function
  // gets similar to the function for distributing matrix entries, see there
  // for additional comments.
  internals::GlobalRowsFromLocal &global_rows = scratch_data->global_rows;
  global_rows.reinit(n_local_dofs);
  make_sorted_row_list (local_dof_indices, global_rows);
  const size_type n_actual_dofs = global_rows.size();

  // create arrays for the column indices that will then be written into the
  // sparsity pattern.
  std::vector<size_type> &cols = scratch_data->columns;
  cols.resize(n_actual_dofs);

  for (size_type i=0; i<n_actual_dofs; ++i)
    {
      std::vector<size_type>::iterator col_ptr = cols.begin();
      const size_type row = global_rows.global_row(i);
      internals::resolve_matrix_row (global_rows, i, 0, n_actual_dofs,
                                     dof_mask, col_ptr);

      // finally, write all the information that accumulated under the given
      // process into the global matrix row and into the vector
      if (col_ptr != cols.begin())
        sparsity_pattern.add_entries(row, cols.begin(), col_ptr,
                                     true);
    }
  internals::set_sparsity_diagonals (global_rows, local_dof_indices,
                                     dof_mask, keep_constrained_entries,
                                     sparsity_pattern);
}




template <typename SparsityType>
void
ConstraintMatrix::
add_entries_local_to_global (const std::vector<size_type> &row_indices,
                             const std::vector<size_type> &col_indices,
                             SparsityType                 &sparsity_pattern,
                             const bool                    keep_constrained_entries,
                             const Table<2,bool>          &dof_mask) const
{
  const size_type n_local_rows = row_indices.size();
  const size_type n_local_cols = col_indices.size();
  bool dof_mask_is_active = false;
  if (dof_mask.n_rows() == n_local_rows && dof_mask.n_cols() == n_local_cols)
    dof_mask_is_active = true;

  // if constrained entries should be kept, need to add rows and columns of
  // those to the sparsity pattern
  if (keep_constrained_entries == true)
    {
      for (size_type i=0; i<row_indices.size(); i++)
        if (is_constrained(row_indices[i]))
          for (size_type j=0; j<col_indices.size(); j++)
            sparsity_pattern.add (row_indices[i], col_indices[j]);
      for (size_type i=0; i<col_indices.size(); i++)
        if (is_constrained(col_indices[i]))
          for (size_type j=0; j<row_indices.size(); j++)
            sparsity_pattern.add (row_indices[j], col_indices[i]);
    }

  // if the dof mask is not active, all we have to do is to add some indices
  // in a matrix format. To do this, we first create an array of all the
  // indices that are to be added. these indices are the local dof indices
  // plus some indices that come from constraints.
  if (dof_mask_is_active == false)
    {
      std::vector<size_type> actual_row_indices (n_local_rows);
      std::vector<size_type> actual_col_indices (n_local_cols);
      make_sorted_row_list (row_indices, actual_row_indices);
      make_sorted_row_list (col_indices, actual_col_indices);
      const size_type n_actual_rows = actual_row_indices.size();

      // now add the indices we collected above to the sparsity pattern. Very
      // easy here - just add the same array to all the rows...
      for (size_type i=0; i<n_actual_rows; ++i)
        sparsity_pattern.add_entries(actual_row_indices[i],
                                     actual_col_indices.begin(),
                                     actual_col_indices.end(),
                                     true);
      return;
    }


  // TODO: implement this
  Assert (false, ExcNotImplemented());
}




template <typename SparsityType>
void
ConstraintMatrix::
add_entries_local_to_global (const std::vector<size_type> &local_dof_indices,
                             SparsityType                 &sparsity_pattern,
                             const bool                    keep_constrained_entries,
                             const Table<2,bool>          &dof_mask,
                             internal::bool2type<true> ) const
{
  // just as the other add_entries_local_to_global function, but now
  // specialized for block matrices.
  Assert (sparsity_pattern.n_rows() == sparsity_pattern.n_cols(), ExcNotQuadratic());
  Assert (sparsity_pattern.n_block_rows() == sparsity_pattern.n_block_cols(),
          ExcNotQuadratic());

  const size_type n_local_dofs = local_dof_indices.size();
  const size_type num_blocks = sparsity_pattern.n_block_rows();

  internals::ConstraintMatrixData<double>::ScratchDataAccessor scratch_data;

  bool dof_mask_is_active = false;
  if (dof_mask.n_rows() == n_local_dofs)
    {
      dof_mask_is_active = true;
      AssertDimension (dof_mask.n_cols(), n_local_dofs);
    }

  if (dof_mask_is_active == false)
    {
      std::vector<size_type> &actual_dof_indices = scratch_data->columns;
      actual_dof_indices.resize(n_local_dofs);
      make_sorted_row_list (local_dof_indices, actual_dof_indices);
      const size_type n_actual_dofs = actual_dof_indices.size();

      // additional construct that also takes care of block indices.
      std::vector<size_type> &block_starts = scratch_data->block_starts;
      block_starts.resize(num_blocks+1);
      internals::make_block_starts (sparsity_pattern, actual_dof_indices,
                                    block_starts);

      for (size_type block=0; block<num_blocks; ++block)
        {
          const size_type next_block = block_starts[block+1];
          for (size_type i=block_starts[block]; i<next_block; ++i)
            {
              Assert (i<n_actual_dofs, ExcInternalError());
              const size_type row = actual_dof_indices[i];
              Assert (row < sparsity_pattern.block(block,0).n_rows(),
                      ExcInternalError());
              std::vector<size_type>::iterator index_it = actual_dof_indices.begin();
              for (size_type block_col = 0; block_col<num_blocks; ++block_col)
                {
                  const size_type next_block_col = block_starts[block_col+1];
                  sparsity_pattern.block(block,block_col).
                  add_entries(row,
                              index_it,
                              actual_dof_indices.begin() + next_block_col,
                              true);
                  index_it = actual_dof_indices.begin() + next_block_col;
                }
            }
        }

      for (size_type i=0; i<n_local_dofs; i++)
        if (is_constrained(local_dof_indices[i]))
          {
            if (keep_constrained_entries == true)
              for (size_type j=0; j<n_local_dofs; j++)
                {
                  sparsity_pattern.add (local_dof_indices[i], local_dof_indices[j]);
                  sparsity_pattern.add (local_dof_indices[j], local_dof_indices[i]);
                }
            else
              sparsity_pattern.add (local_dof_indices[i], local_dof_indices[i]);
          }

      return;
    }

  // difficult case with dof_mask, similar to the distribute_local_to_global
  // function for block matrices
  internals::GlobalRowsFromLocal &global_rows = scratch_data->global_rows;
  global_rows.reinit(n_local_dofs);
  make_sorted_row_list (local_dof_indices, global_rows);
  const size_type n_actual_dofs = global_rows.size();

  // additional construct that also takes care of block indices.
  std::vector<size_type> &block_starts = scratch_data->block_starts;
  block_starts.resize(num_blocks+1);
  internals::make_block_starts(sparsity_pattern, global_rows, block_starts);

  std::vector<size_type> &cols = scratch_data->columns;
  cols.resize(n_actual_dofs);

  // the basic difference to the non-block variant from now onwards is that we
  // go through the blocks of the matrix separately.
  for (size_type block=0; block<num_blocks; ++block)
    {
      const size_type next_block = block_starts[block+1];
      for (size_type i=block_starts[block]; i<next_block; ++i)
        {
          const size_type row = global_rows.global_row(i);
          for (size_type block_col=0; block_col<num_blocks; ++block_col)
            {
              const size_type begin_block = block_starts[block_col],
                              end_block = block_starts[block_col+1];
              std::vector<size_type>::iterator col_ptr = cols.begin();
              internals::resolve_matrix_row (global_rows, i, begin_block,
                                             end_block, dof_mask, col_ptr);

              sparsity_pattern.block(block, block_col).add_entries(row,
                                                                   cols.begin(),
                                                                   col_ptr,
                                                                   true);
            }
        }
    }

  internals::set_sparsity_diagonals (global_rows, local_dof_indices,
                                     dof_mask, keep_constrained_entries,
                                     sparsity_pattern);
}


DEAL_II_NAMESPACE_CLOSE

#endif
