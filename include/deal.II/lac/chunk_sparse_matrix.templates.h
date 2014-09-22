// ---------------------------------------------------------------------
//
// Copyright (C) 2008 - 2014 by the deal.II authors
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

#ifndef __deal2__chunk_sparse_matrix_templates_h
#define __deal2__chunk_sparse_matrix_templates_h


#include <deal.II/base/template_constraints.h>
#include <deal.II/base/parallel.h>
#include <deal.II/lac/chunk_sparse_matrix.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>


#include <ostream>
#include <iomanip>
#include <algorithm>
#include <functional>
#include <cmath>

#include <vector>
#include <numeric>

#include <deal.II/base/thread_management.h>

DEAL_II_NAMESPACE_OPEN


namespace internal
{
//TODO: the goal of the ChunkSparseMatrix class is to stream data and use
// the vectorization features of modern processors. to make this happen,
// we will have to vectorize the functions in the following namespace, either
// by hand or by using, for example, optimized BLAS versions for them.
  namespace ChunkSparseMatrix
  {
    /**
     * Declare type for container size.
     */
    typedef types::global_dof_index size_type;

    /**
     * Add the result of multiplying a chunk
     * of size chunk_size times chunk_size by
     * a source vector fragment of size
     * chunk_size to the destination vector
     * fragment.
     */
    template <typename MatrixIterator,
              typename SrcIterator,
              typename DstIterator>
    inline
    void
    chunk_vmult_add (const size_type      chunk_size,
                     const MatrixIterator matrix,
                     const SrcIterator    src,
                     DstIterator          dst)
    {
      MatrixIterator matrix_row = matrix;

      for (size_type i=0; i<chunk_size;
           ++i,  matrix_row += chunk_size)
        {
          typename std::iterator_traits<DstIterator>::value_type
          sum = 0;

          for (size_type j=0; j<chunk_size; ++j)
            sum += matrix_row[j] * src[j];

          dst[i] += sum;
        }
    }



    /**
     * Like the previous function, but subtract. We need this for computing
     * the residual.
     */
    template <typename MatrixIterator,
              typename SrcIterator,
              typename DstIterator>
    inline
    void
    chunk_vmult_subtract (const size_type chunk_size,
                          const MatrixIterator matrix,
                          const SrcIterator    src,
                          DstIterator          dst)
    {
      MatrixIterator matrix_row = matrix;

      for (size_type i=0; i<chunk_size;
           ++i,  matrix_row += chunk_size)
        {
          typename std::iterator_traits<DstIterator>::value_type
          sum = 0;

          for (size_type j=0; j<chunk_size; ++j)
            sum += matrix_row[j] * src[j];

          dst[i] -= sum;
        }
    }


    /**
     * Add the result of multiplying the transpose of a chunk of size
     * chunk_size times chunk_size by a source vector fragment of size
     * chunk_size to the destination vector fragment.
     */
    template <typename MatrixIterator,
              typename SrcIterator,
              typename DstIterator>
    inline
    void
    chunk_Tvmult_add (const size_type      chunk_size,
                      const MatrixIterator matrix,
                      const SrcIterator    src,
                      DstIterator          dst)
    {
      for (size_type i=0; i<chunk_size; ++i)
        {
          typename std::iterator_traits<DstIterator>::value_type
          sum = 0;

          for (size_type j=0; j<chunk_size; ++j)
            sum += matrix[j*chunk_size+i] * src[j];

          dst[i] += sum;
        }
    }


    /**
     * Produce the result of the matrix scalar product $u^TMv$ for an
     * individual chunk.
     */
    template <typename result_type,
              typename MatrixIterator,
              typename SrcIterator1,
              typename SrcIterator2>
    inline
    result_type
    chunk_matrix_scalar_product (const size_type      chunk_size,
                                 const MatrixIterator matrix,
                                 const SrcIterator1   u,
                                 const SrcIterator2   v)
    {
      result_type result = 0;

      MatrixIterator matrix_row = matrix;

      for (size_type i=0; i<chunk_size;
           ++i,  matrix_row += chunk_size)
        {
          typename std::iterator_traits<SrcIterator2>::value_type
          sum = 0;

          for (size_type j=0; j<chunk_size; ++j)
            sum += matrix_row[j] * v[j];

          result += u[i] * sum;
        }

      return result;
    }



    /**
     * Perform a vmult_add using the ChunkSparseMatrix data structures, but
     * only using a subinterval of the matrix rows.
     *
     * In the sequential case, this function is called on all rows, in the
     * parallel case it may be called on a subrange, at the discretion of the
     * task scheduler.
     */
    template <typename number,
              typename InVector,
              typename OutVector>
    void vmult_add_on_subrange (const ChunkSparsityPattern &cols,
                                const unsigned int  begin_row,
                                const unsigned int  end_row,
                                const number       *values,
                                const std::size_t  *rowstart,
                                const size_type    *colnums,
                                const InVector     &src,
                                OutVector          &dst)
    {
      const size_type m = cols.n_rows();
      const size_type n = cols.n_cols();
      const size_type chunk_size = cols.get_chunk_size();

      // loop over all chunks. note that we need to treat the last chunk row
      // and column differently if they have padding elements
      const size_type n_filled_last_rows = m % chunk_size;
      const size_type n_filled_last_cols = n % chunk_size;

      const size_type last_regular_row = n_filled_last_rows > 0 ?
                                         std::min(m/chunk_size,
                                                  static_cast<size_type>(end_row)) :
                                         end_row;
      const size_type irregular_col = n/chunk_size;

      typename OutVector::iterator dst_ptr = dst.begin()+chunk_size*begin_row;
      const number *val_ptr= &values[rowstart[begin_row]*chunk_size*chunk_size];
      const size_type *colnum_ptr = &colnums[rowstart[begin_row]];
      for (unsigned int chunk_row=begin_row; chunk_row<last_regular_row;
           ++chunk_row)
        {
          const number *const val_end_of_row = &values[rowstart[chunk_row+1] *
                                                       chunk_size * chunk_size];
          while (val_ptr != val_end_of_row)
            {
              if (*colnum_ptr != irregular_col)
                chunk_vmult_add (chunk_size,
                                 val_ptr,
                                 src.begin() + *colnum_ptr * chunk_size,
                                 dst_ptr);
              else
                // we're at a chunk column that has padding
                for (size_type r=0; r<chunk_size; ++r)
                  for (size_type c=0; c<n_filled_last_cols; ++c)
                    dst_ptr[r] += (val_ptr[r*chunk_size + c] *
                                   src(*colnum_ptr * chunk_size + c));

              ++colnum_ptr;
              val_ptr += chunk_size * chunk_size;
            }

          dst_ptr += chunk_size;
        }

      // now deal with last chunk row if necessary
      if (n_filled_last_rows > 0 && end_row == (m/chunk_size+1))
        {
          const size_type chunk_row = last_regular_row;

          const number *const val_end_of_row = &values[rowstart[chunk_row+1] *
                                                       chunk_size * chunk_size];
          while (val_ptr != val_end_of_row)
            {
              if (*colnum_ptr != irregular_col)
                {
                  // we're at a chunk row but not column that has padding
                  for (size_type r=0; r<n_filled_last_rows; ++r)
                    for (size_type c=0; c<chunk_size; ++c)
                      dst_ptr[r]
                      += (val_ptr[r*chunk_size + c] *
                          src(*colnum_ptr * chunk_size + c));
                }
              else
                // we're at a chunk row and column that has padding
                for (size_type r=0; r<n_filled_last_rows; ++r)
                  for (size_type c=0; c<n_filled_last_cols; ++c)
                    dst_ptr[r]
                    += (val_ptr[r*chunk_size + c] *
                        src(*colnum_ptr * chunk_size + c));

              ++colnum_ptr;
              val_ptr += chunk_size * chunk_size;
            }
        }
      Assert(std::size_t(colnum_ptr-&colnums[0]) == rowstart[end_row],
             ExcInternalError());
      Assert(std::size_t(val_ptr-&values[0]) ==
             rowstart[end_row] * chunk_size * chunk_size,
             ExcInternalError());
    }
  }
}



template <typename number>
ChunkSparseMatrix<number>::ChunkSparseMatrix ()
  :
  cols(0, "ChunkSparseMatrix"),
  val(0),
  max_len(0)
{}



template <typename number>
ChunkSparseMatrix<number>::ChunkSparseMatrix (const ChunkSparseMatrix &m)
  :
  Subscriptor (m),
  cols(0, "ChunkSparseMatrix"),
  val(0),
  max_len(0)
{
  Assert (m.cols==0, ExcInvalidConstructorCall());
  Assert (m.val==0, ExcInvalidConstructorCall());
  Assert (m.max_len==0, ExcInvalidConstructorCall());
}



template <typename number>
ChunkSparseMatrix<number> &
ChunkSparseMatrix<number>::operator = (const ChunkSparseMatrix<number> &m)
{
  Assert (m.cols==0, ExcInvalidConstructorCall());
  Assert (m.val==0, ExcInvalidConstructorCall());
  Assert (m.max_len==0, ExcInvalidConstructorCall());

  return *this;
}



template <typename number>
ChunkSparseMatrix<number>::ChunkSparseMatrix (const ChunkSparsityPattern &c)
  :
  cols(0, "ChunkSparseMatrix"),
  val(0),
  max_len(0)
{
  reinit (c);
}



template <typename number>
ChunkSparseMatrix<number>::ChunkSparseMatrix (const ChunkSparsityPattern &c,
                                              const IdentityMatrix  &id)
  :
  cols(0, "ChunkSparseMatrix"),
  val(0),
  max_len(0)
{
  Assert (c.n_rows() == id.m(), ExcDimensionMismatch (c.n_rows(), id.m()));
  Assert (c.n_cols() == id.n(), ExcDimensionMismatch (c.n_cols(), id.n()));

  reinit (c);
  for (size_type i=0; i<n(); ++i)
    this->set(i,i,1.);
}



template <typename number>
ChunkSparseMatrix<number>::~ChunkSparseMatrix ()
{
  cols = 0;

  if (val != 0)
    delete[] val;
}



namespace internal
{
  namespace ChunkSparseMatrix
  {
    template<typename T>
    void zero_subrange (const unsigned int begin,
                        const unsigned int end,
                        T *dst)
    {
      std::memset (dst+begin,0,(end-begin)*sizeof(T));
    }
  }
}



template <typename number>
ChunkSparseMatrix<number> &
ChunkSparseMatrix<number>::operator = (const double d)
{
  Assert (d==0, ExcScalarAssignmentOnlyForZeroValue());

  Assert (cols != 0, ExcNotInitialized());
  Assert (cols->sparsity_pattern.compressed || cols->empty(),
          ChunkSparsityPattern::ExcNotCompressed());

  // do initial zeroing of elements in parallel. Try to achieve a similar
  // layout as when doing matrix-vector products, as on some NUMA systems, a
  // memory block is assigned to memory banks where the first access is
  // generated. For sparse matrices, the first operations is usually the
  // operator=. The grain size is chosen to reflect the number of rows in
  // minimum_parallel_grain_size, weighted by the number of nonzero entries
  // per row on average.
  const unsigned int matrix_size = cols->sparsity_pattern.n_nonzero_elements()
                                   * cols->chunk_size * cols->chunk_size;
  const unsigned int grain_size =
    internal::SparseMatrix::minimum_parallel_grain_size *
    (matrix_size+m()) / m();
  if (matrix_size>grain_size)
    parallel::apply_to_subranges (0U, matrix_size,
                                  std_cxx11::bind(&internal::ChunkSparseMatrix::template
                                                  zero_subrange<number>,
                                                  std_cxx11::_1, std_cxx11::_2,
                                                  val),
                                  grain_size);
  else if (matrix_size > 0)
    std::memset (&val[0], 0, matrix_size*sizeof(number));

  return *this;
}



template <typename number>
ChunkSparseMatrix<number> &
ChunkSparseMatrix<number>::operator= (const IdentityMatrix  &id)
{
  Assert (cols->n_rows() == id.m(),
          ExcDimensionMismatch (cols->n_rows(), id.m()));
  Assert (cols->n_cols() == id.n(),
          ExcDimensionMismatch (cols->n_cols(), id.n()));

  *this = 0;
  for (size_type i=0; i<n(); ++i)
    this->set(i,i,1.);

  return *this;
}



template <typename number>
void
ChunkSparseMatrix<number>::reinit (const ChunkSparsityPattern &sparsity)
{
  cols = &sparsity;

  if (cols->empty())
    {
      if (val != 0)
        delete[] val;
      val = 0;
      max_len = 0;
      return;
    }

  // allocate not just m() * n() elements but enough so that we can store full
  // chunks. this entails some padding elements
  const size_type chunk_size = cols->get_chunk_size();
  const size_type N = cols->sparsity_pattern.n_nonzero_elements() *
                      chunk_size * chunk_size;
  if (N > max_len || max_len == 0)
    {
      if (val != 0)
        delete[] val;
      val = new number[N];
      max_len = N;
    }

  // fill with zeros. do not just fill N elements but all that we allocated to
  // ensure that also the padding elements are zero and not left at previous
  // values
  this->operator=(0.);
}



template <typename number>
void
ChunkSparseMatrix<number>::clear ()
{
  cols = 0;
  if (val) delete[] val;
  val = 0;
  max_len = 0;
}



template <typename number>
bool
ChunkSparseMatrix<number>::empty () const
{
  if (cols == 0)
    return true;
  else
    return cols->empty();
}



template <typename number>
typename ChunkSparseMatrix<number>::size_type
ChunkSparseMatrix<number>::n_nonzero_elements () const
{
  Assert (cols != 0, ExcNotInitialized());
  return cols->n_nonzero_elements ();
}



template <typename number>
typename ChunkSparseMatrix<number>::size_type
ChunkSparseMatrix<number>::n_actually_nonzero_elements () const
{
  Assert (cols != 0, ExcNotInitialized());

  // count those elements that are nonzero, even if they lie in the padding
  // around the matrix. since we have the invariant that padding elements are
  // zero, nothing bad can happen here
  const size_type chunk_size = cols->get_chunk_size();
  return std::count_if(&val[0],
                       &val[cols->sparsity_pattern.n_nonzero_elements () *
                            chunk_size * chunk_size],
                       std::bind2nd(std::not_equal_to<double>(), 0));
}



template <typename number>
void
ChunkSparseMatrix<number>::symmetrize ()
{
  Assert (cols != 0, ExcNotInitialized());
  Assert (cols->rows == cols->cols, ExcNotQuadratic());

  Assert (false, ExcNotImplemented());
}



template <typename number>
template <typename somenumber>
ChunkSparseMatrix<number> &
ChunkSparseMatrix<number>::copy_from (const ChunkSparseMatrix<somenumber> &matrix)
{
  Assert (cols != 0, ExcNotInitialized());
  Assert (val != 0, ExcNotInitialized());
  Assert (cols == matrix.cols, ExcDifferentChunkSparsityPatterns());

  // copy everything, including padding elements
  const size_type chunk_size = cols->get_chunk_size();
  std::copy (&matrix.val[0],
             &matrix.val[cols->sparsity_pattern.n_nonzero_elements()
                         * chunk_size * chunk_size],
             &val[0]);

  return *this;
}



template <typename number>
template <typename somenumber>
void
ChunkSparseMatrix<number>::copy_from (const FullMatrix<somenumber> &matrix)
{
  // first delete previous content
  *this = 0;

  // then copy old matrix
  for (size_type row=0; row<matrix.m(); ++row)
    for (size_type col=0; col<matrix.n(); ++col)
      if (matrix(row,col) != 0)
        set (row, col, matrix(row,col));
}



template <typename number>
template <typename somenumber>
void
ChunkSparseMatrix<number>::add (const number factor,
                                const ChunkSparseMatrix<somenumber> &matrix)
{
  Assert (cols != 0, ExcNotInitialized());
  Assert (val != 0, ExcNotInitialized());
  Assert (cols == matrix.cols, ExcDifferentChunkSparsityPatterns());

  // add everything, including padding elements
  const size_type     chunk_size = cols->get_chunk_size();
  number             *val_ptr    = &val[0];
  const somenumber   *matrix_ptr = &matrix.val[0];
  const number *const end_ptr    = &val[cols->sparsity_pattern.n_nonzero_elements()
                                        * chunk_size * chunk_size];

  while (val_ptr != end_ptr)
    *val_ptr++ += factor **matrix_ptr++;
}



template <typename number>
void
ChunkSparseMatrix<number>::extract_row_copy (const size_type row,
                                             const size_type array_length,
                                             size_type &row_length,
                                             size_type *column_indices,
                                             number *values) const
{
  AssertIndexRange(cols->row_length(row), array_length+1);
  AssertIndexRange(row, m());
  const unsigned int chunk_size = cols->get_chunk_size();
  const size_type reduced_row = row/chunk_size;

  SparsityPattern::iterator it = cols->sparsity_pattern.begin(reduced_row),
                            itend = cols->sparsity_pattern.end(reduced_row);
  const number *val_ptr = &val[(it-cols->sparsity_pattern.begin(0))*chunk_size*chunk_size
                               +(row%chunk_size)*chunk_size];

  // find out if we did padding and if this row is affected by it
  if (cols->n_cols() % chunk_size == 0)
    {
      for ( ; it < itend; ++it)
        {
          for (unsigned int c=0; c<chunk_size; ++c)
            {
              *values++ = val_ptr[c];
              *column_indices++ = it->column()*chunk_size+c;
            }
          val_ptr += chunk_size*chunk_size;
        }
      row_length = chunk_size * (cols->sparsity_pattern.row_length(reduced_row));
    }
  else
    {
      const unsigned int last_chunk_size = cols->n_cols() % chunk_size;
      row_length = 0;
      for ( ; it < itend; ++it)
        {
          const unsigned int next_chunk_size =
            (it->column()==cols->sparsity_pattern.n_cols()-1) ?
            last_chunk_size : chunk_size;
          for (unsigned int c=0; c<next_chunk_size; ++c, ++row_length)
            {
              *values++ = val_ptr[c];
              *column_indices++ = it->column()*chunk_size+c;
            }
          val_ptr += chunk_size*chunk_size;
        }
    }
}



template <typename number>
template <class OutVector, class InVector>
void
ChunkSparseMatrix<number>::vmult (OutVector &dst,
                                  const InVector &src) const
{
  Assert (cols != 0, ExcNotInitialized());
  Assert (val != 0, ExcNotInitialized());
  Assert(m() == dst.size(), ExcDimensionMismatch(m(),dst.size()));
  Assert(n() == src.size(), ExcDimensionMismatch(n(),src.size()));

  Assert (!PointerComparison::equal(&src, &dst), ExcSourceEqualsDestination());

  // set the output vector to zero and then add to it the contributions of
  // vmults from individual chunks. this is what vmult_add does
  dst = 0;
  vmult_add (dst, src);
}



template <typename number>
template <class OutVector, class InVector>
void
ChunkSparseMatrix<number>::Tvmult (OutVector &dst,
                                   const InVector &src) const
{
  Assert (val != 0, ExcNotInitialized());
  Assert (cols != 0, ExcNotInitialized());
  Assert(n() == dst.size(), ExcDimensionMismatch(n(),dst.size()));
  Assert(m() == src.size(), ExcDimensionMismatch(m(),src.size()));

  Assert (!PointerComparison::equal(&src, &dst), ExcSourceEqualsDestination());

  Assert (cols != 0, ExcNotInitialized());
  Assert (val != 0, ExcNotInitialized());
  Assert(m() == dst.size(), ExcDimensionMismatch(m(),dst.size()));
  Assert(n() == src.size(), ExcDimensionMismatch(n(),src.size()));

  Assert (!PointerComparison::equal(&src, &dst), ExcSourceEqualsDestination());

  // set the output vector to zero and then add to it the contributions of
  // vmults from individual chunks. this is what vmult_add does
  dst = 0;
  Tvmult_add (dst, src);
}



template <typename number>
template <class OutVector, class InVector>
void
ChunkSparseMatrix<number>::vmult_add (OutVector &dst,
                                      const InVector &src) const
{
  Assert (cols != 0, ExcNotInitialized());
  Assert (val != 0, ExcNotInitialized());
  Assert(m() == dst.size(), ExcDimensionMismatch(m(),dst.size()));
  Assert(n() == src.size(), ExcDimensionMismatch(n(),src.size()));

  Assert (!PointerComparison::equal(&src, &dst), ExcSourceEqualsDestination());
  parallel::apply_to_subranges (0U, cols->sparsity_pattern.n_rows(),
                                std_cxx11::bind (&internal::ChunkSparseMatrix::vmult_add_on_subrange
                                                 <number,InVector,OutVector>,
                                                 std_cxx11::cref(*cols),
                                                 std_cxx11::_1, std_cxx11::_2,
                                                 val,
                                                 cols->sparsity_pattern.rowstart,
                                                 cols->sparsity_pattern.colnums,
                                                 std_cxx11::cref(src),
                                                 std_cxx11::ref(dst)),
                                internal::SparseMatrix::minimum_parallel_grain_size/cols->chunk_size+1);

}


template <typename number>
template <class OutVector, class InVector>
void
ChunkSparseMatrix<number>::Tvmult_add (OutVector &dst,
                                       const InVector &src) const
{
  Assert (cols != 0, ExcNotInitialized());
  Assert (val != 0, ExcNotInitialized());
  Assert(m() == dst.size(), ExcDimensionMismatch(m(),dst.size()));
  Assert(n() == src.size(), ExcDimensionMismatch(n(),src.size()));

  Assert (!PointerComparison::equal(&src, &dst), ExcSourceEqualsDestination());

  const size_type n_chunk_rows = cols->sparsity_pattern.n_rows();

  // loop over all chunks. note that we need to treat the last chunk row and
  // column differently if they have padding elements
  const bool rows_have_padding = (m() % cols->chunk_size != 0),
             cols_have_padding = (n() % cols->chunk_size != 0);

  const size_type n_regular_chunk_rows
    = (rows_have_padding ?
       n_chunk_rows-1 :
       n_chunk_rows);

  // like in vmult_add, but don't keep an iterator into dst around since we're
  // not traversing it sequentially this time
  const number    *val_ptr    = val;
  const size_type *colnum_ptr = cols->sparsity_pattern.colnums;

  for (size_type chunk_row=0; chunk_row<n_regular_chunk_rows; ++chunk_row)
    {
      const number *const val_end_of_row = &val[cols->sparsity_pattern.rowstart[chunk_row+1]
                                                * cols->chunk_size
                                                * cols->chunk_size];
      while (val_ptr != val_end_of_row)
        {
          if ((cols_have_padding == false)
              ||
              (*colnum_ptr != cols->sparsity_pattern.n_cols()-1))
            internal::ChunkSparseMatrix::chunk_Tvmult_add
            (cols->chunk_size,
             val_ptr,
             src.begin() + chunk_row * cols->chunk_size,
             dst.begin() + *colnum_ptr * cols->chunk_size);
          else
            // we're at a chunk column that has padding
            for (size_type r=0; r<cols->chunk_size; ++r)
              for (size_type c=0; c<n() % cols->chunk_size; ++c)
                dst(*colnum_ptr * cols->chunk_size + c)
                += (val_ptr[r*cols->chunk_size + c] *
                    src(chunk_row * cols->chunk_size + r));

          ++colnum_ptr;
          val_ptr += cols->chunk_size * cols->chunk_size;
        }
    }

  // now deal with last chunk row if necessary
  if (rows_have_padding)
    {
      const size_type chunk_row = n_chunk_rows - 1;

      const number *const val_end_of_row = &val[cols->sparsity_pattern.rowstart[chunk_row+1]
                                                * cols->chunk_size
                                                * cols->chunk_size];
      while (val_ptr != val_end_of_row)
        {
          if ((cols_have_padding == false)
              ||
              (*colnum_ptr != cols->sparsity_pattern.n_cols()-1))
            {
              // we're at a chunk row but not column that has padding
              for (size_type r=0; r<m() % cols->chunk_size; ++r)
                for (size_type c=0; c<cols->chunk_size; ++c)
                  dst(*colnum_ptr * cols->chunk_size + c)
                  += (val_ptr[r*cols->chunk_size + c] *
                      src(chunk_row * cols->chunk_size + r));
            }
          else
            // we're at a chunk row and column that has padding
            for (size_type r=0; r<m() % cols->chunk_size; ++r)
              for (size_type c=0; c<n() % cols->chunk_size; ++c)
                dst(*colnum_ptr * cols->chunk_size + c)
                += (val_ptr[r*cols->chunk_size + c] *
                    src(chunk_row * cols->chunk_size + r));

          ++colnum_ptr;
          val_ptr += cols->chunk_size * cols->chunk_size;
        }
    }
}


template <typename number>
template <typename somenumber>
somenumber
ChunkSparseMatrix<number>::matrix_norm_square (const Vector<somenumber> &v) const
{
  Assert (cols != 0, ExcNotInitialized());
  Assert (val != 0, ExcNotInitialized());
  Assert(m() == v.size(), ExcDimensionMismatch(m(),v.size()));
  Assert(n() == v.size(), ExcDimensionMismatch(n(),v.size()));

  somenumber result = 0;

  ////////////////
  // like matrix_scalar_product, except that the two vectors are now the same

  const size_type n_chunk_rows = cols->sparsity_pattern.n_rows();

  // loop over all chunks. note that we need to treat the last chunk row and
  // column differently if they have padding elements
  const bool rows_have_padding = (m() % cols->chunk_size != 0),
             cols_have_padding = (n() % cols->chunk_size != 0);

  const size_type n_regular_chunk_rows
    = (rows_have_padding ?
       n_chunk_rows-1 :
       n_chunk_rows);

  const number    *val_ptr    = val;
  const size_type *colnum_ptr = cols->sparsity_pattern.colnums;
  typename Vector<somenumber>::const_iterator v_ptr = v.begin();

  for (size_type chunk_row=0; chunk_row<n_regular_chunk_rows; ++chunk_row)
    {
      const number *const val_end_of_row = &val[cols->sparsity_pattern.rowstart[chunk_row+1]
                                                * cols->chunk_size
                                                * cols->chunk_size];
      while (val_ptr != val_end_of_row)
        {
          if ((cols_have_padding == false)
              ||
              (*colnum_ptr != cols->sparsity_pattern.n_cols()-1))
            result +=
              internal::ChunkSparseMatrix::
              chunk_matrix_scalar_product<somenumber>
              (cols->chunk_size,
               val_ptr,
               v_ptr,
               v.begin() + *colnum_ptr * cols->chunk_size);
          else
            // we're at a chunk column that has padding
            for (size_type r=0; r<cols->chunk_size; ++r)
              for (size_type c=0; c<n() % cols->chunk_size; ++c)
                result
                +=
                  v(chunk_row * cols->chunk_size + r)
                  * (val_ptr[r*cols->chunk_size + c] *
                     v(*colnum_ptr * cols->chunk_size + c));

          ++colnum_ptr;
          val_ptr += cols->chunk_size * cols->chunk_size;
        }


      v_ptr += cols->chunk_size;
    }

  // now deal with last chunk row if necessary
  if (rows_have_padding)
    {
      const size_type chunk_row = n_chunk_rows - 1;

      const number *const val_end_of_row = &val[cols->sparsity_pattern.rowstart[chunk_row+1]
                                                * cols->chunk_size
                                                * cols->chunk_size];
      while (val_ptr != val_end_of_row)
        {
          if ((cols_have_padding == false)
              ||
              (*colnum_ptr != cols->sparsity_pattern.n_cols()-1))
            {
              // we're at a chunk row but not column that has padding
              for (size_type r=0; r<m() % cols->chunk_size; ++r)
                for (size_type c=0; c<cols->chunk_size; ++c)
                  result
                  +=
                    v(chunk_row * cols->chunk_size + r)
                    * (val_ptr[r*cols->chunk_size + c] *
                       v(*colnum_ptr * cols->chunk_size + c));
            }
          else
            // we're at a chunk row and column that has padding
            for (size_type r=0; r<m() % cols->chunk_size; ++r)
              for (size_type c=0; c<n() % cols->chunk_size; ++c)
                result
                +=
                  v(chunk_row * cols->chunk_size + r)
                  * (val_ptr[r*cols->chunk_size + c] *
                     v(*colnum_ptr * cols->chunk_size + c));

          ++colnum_ptr;
          val_ptr += cols->chunk_size * cols->chunk_size;
        }
    }

  return result;
}



template <typename number>
template <typename somenumber>
somenumber
ChunkSparseMatrix<number>::matrix_scalar_product (const Vector<somenumber> &u,
                                                  const Vector<somenumber> &v) const
{
  Assert (cols != 0, ExcNotInitialized());
  Assert (val != 0, ExcNotInitialized());
  Assert(m() == u.size(), ExcDimensionMismatch(m(),u.size()));
  Assert(n() == v.size(), ExcDimensionMismatch(n(),v.size()));

  // the following works like the vmult_add function
  somenumber result = 0;

  const size_type n_chunk_rows = cols->sparsity_pattern.n_rows();

  // loop over all chunks. note that we need to treat the last chunk row and
  // column differently if they have padding elements
  const bool rows_have_padding = (m() % cols->chunk_size != 0),
             cols_have_padding = (n() % cols->chunk_size != 0);

  const size_type n_regular_chunk_rows
    = (rows_have_padding ?
       n_chunk_rows-1 :
       n_chunk_rows);

  const number    *val_ptr    = val;
  const size_type *colnum_ptr = cols->sparsity_pattern.colnums;
  typename Vector<somenumber>::const_iterator u_ptr = u.begin();

  for (size_type chunk_row=0; chunk_row<n_regular_chunk_rows; ++chunk_row)
    {
      const number *const val_end_of_row = &val[cols->sparsity_pattern.rowstart[chunk_row+1]
                                                * cols->chunk_size
                                                * cols->chunk_size];
      while (val_ptr != val_end_of_row)
        {
          if ((cols_have_padding == false)
              ||
              (*colnum_ptr != cols->sparsity_pattern.n_cols()-1))
            result +=
              internal::ChunkSparseMatrix::
              chunk_matrix_scalar_product<somenumber>
              (cols->chunk_size,
               val_ptr,
               u_ptr,
               v.begin() + *colnum_ptr * cols->chunk_size);
          else
            // we're at a chunk column that has padding
            for (size_type r=0; r<cols->chunk_size; ++r)
              for (size_type c=0; c<n() % cols->chunk_size; ++c)
                result
                +=
                  u(chunk_row * cols->chunk_size + r)
                  * (val_ptr[r*cols->chunk_size + c] *
                     v(*colnum_ptr * cols->chunk_size + c));

          ++colnum_ptr;
          val_ptr += cols->chunk_size * cols->chunk_size;
        }


      u_ptr += cols->chunk_size;
    }

  // now deal with last chunk row if necessary
  if (rows_have_padding)
    {
      const size_type chunk_row = n_chunk_rows - 1;

      const number *const val_end_of_row = &val[cols->sparsity_pattern.rowstart[chunk_row+1]
                                                * cols->chunk_size
                                                * cols->chunk_size];
      while (val_ptr != val_end_of_row)
        {
          if ((cols_have_padding == false)
              ||
              (*colnum_ptr != cols->sparsity_pattern.n_cols()-1))
            {
              // we're at a chunk row but not column that has padding
              for (size_type r=0; r<m() % cols->chunk_size; ++r)
                for (size_type c=0; c<cols->chunk_size; ++c)
                  result
                  +=
                    u(chunk_row * cols->chunk_size + r)
                    * (val_ptr[r*cols->chunk_size + c] *
                       v(*colnum_ptr * cols->chunk_size + c));
            }
          else
            // we're at a chunk row and column that has padding
            for (size_type r=0; r<m() % cols->chunk_size; ++r)
              for (size_type c=0; c<n() % cols->chunk_size; ++c)
                result
                +=
                  u(chunk_row * cols->chunk_size + r)
                  * (val_ptr[r*cols->chunk_size + c] *
                     v(*colnum_ptr * cols->chunk_size + c));

          ++colnum_ptr;
          val_ptr += cols->chunk_size * cols->chunk_size;
        }
    }

  return result;
}



template <typename number>
typename ChunkSparseMatrix<number>::real_type
ChunkSparseMatrix<number>::l1_norm () const
{
  Assert (cols != 0, ExcNotInitialized());
  Assert (val != 0, ExcNotInitialized());

  const size_type n_chunk_rows = cols->sparsity_pattern.n_rows();

  // loop over all rows and columns; it is safe to also loop over the padding
  // elements (they are zero) if we make sure that the vector into which we
  // sum column sums is large enough
  Vector<real_type> column_sums(cols->sparsity_pattern.n_cols() *
                                cols->chunk_size);

  for (size_type chunk_row=0; chunk_row<n_chunk_rows; ++chunk_row)
    for (size_type j=cols->sparsity_pattern.rowstart[chunk_row];
         j<cols->sparsity_pattern.rowstart[chunk_row+1] ; ++j)
      for (size_type r=0; r<cols->chunk_size; ++r)
        for (size_type s=0; s<cols->chunk_size; ++s)
          column_sums(cols->sparsity_pattern.colnums[j] *
                      cols->chunk_size + s) +=
                        numbers::NumberTraits<number>::abs(val[j * cols->chunk_size *
                                                               cols->chunk_size +
                                                               r * cols->chunk_size +
                                                               s]);

  return column_sums.linfty_norm();
}



template <typename number>
typename ChunkSparseMatrix<number>::real_type
ChunkSparseMatrix<number>::linfty_norm () const
{
  Assert (cols != 0, ExcNotInitialized());
  Assert (val != 0, ExcNotInitialized());

  // this function works like l1_norm(). it can be made more efficient
  // (without allocating a temporary vector) as is done in the SparseMatrix
  // class but since it is rarely called in time critical places it is
  // probably not worth it
  const size_type n_chunk_rows = cols->sparsity_pattern.n_rows();

  // loop over all rows and columns; it is safe to also loop over the padding
  // elements (they are zero) if we make sure that the vector into which we
  // sum column sums is large enough
  Vector<real_type> row_sums(cols->sparsity_pattern.n_rows() *
                             cols->chunk_size);

  for (size_type chunk_row=0; chunk_row<n_chunk_rows; ++chunk_row)
    for (size_type j=cols->sparsity_pattern.rowstart[chunk_row];
         j<cols->sparsity_pattern.rowstart[chunk_row+1] ; ++j)
      for (size_type r=0; r<cols->chunk_size; ++r)
        for (size_type s=0; s<cols->chunk_size; ++s)
          row_sums(chunk_row * cols->chunk_size + r) +=
            numbers::NumberTraits<number>::abs(val[j * cols->chunk_size *
                                                   cols->chunk_size +
                                                   r * cols->chunk_size +
                                                   s]);

  return row_sums.linfty_norm();
}



template <typename number>
typename ChunkSparseMatrix<number>::real_type
ChunkSparseMatrix<number>::frobenius_norm () const
{
  // simply add up all entries in the sparsity pattern, without taking any
  // reference to rows or columns
  //
  // padding elements are zero, so we can add them up as well
  real_type norm_sqr = 0;
  for (const number *ptr = &val[0]; ptr != &val[max_len]; ++ptr)
    norm_sqr +=  numbers::NumberTraits<number>::abs_square(*ptr);

  return std::sqrt (norm_sqr);
}



template <typename number>
template <typename somenumber>
somenumber
ChunkSparseMatrix<number>::residual (Vector<somenumber>       &dst,
                                     const Vector<somenumber> &u,
                                     const Vector<somenumber> &b) const
{
  Assert (cols != 0, ExcNotInitialized());
  Assert (val != 0, ExcNotInitialized());
  Assert(m() == dst.size(), ExcDimensionMismatch(m(),dst.size()));
  Assert(m() == b.size(), ExcDimensionMismatch(m(),b.size()));
  Assert(n() == u.size(), ExcDimensionMismatch(n(),u.size()));

  Assert (&u != &dst, ExcSourceEqualsDestination());

  // set dst=b, then subtract the result of A*u from it. since the purpose of
  // the current class is to promote streaming of data rather than more random
  // access patterns, breaking things up into two loops may be reasonable
  dst = b;

  /////////
  // the rest of this function is like vmult_add, except that we subtract
  // rather than add A*u
  /////////
  const size_type n_chunk_rows = cols->sparsity_pattern.n_rows();

  // loop over all chunks. note that we need to treat the last chunk row and
  // column differently if they have padding elements
  const bool rows_have_padding = (m() % cols->chunk_size != 0),
             cols_have_padding = (n() % cols->chunk_size != 0);

  const size_type n_regular_chunk_rows
    = (rows_have_padding ?
       n_chunk_rows-1 :
       n_chunk_rows);

  const number       *val_ptr    = val;
  const size_type *colnum_ptr = cols->sparsity_pattern.colnums;
  typename Vector<somenumber>::iterator dst_ptr = dst.begin();

  for (size_type chunk_row=0; chunk_row<n_regular_chunk_rows; ++chunk_row)
    {
      const number *const val_end_of_row = &val[cols->sparsity_pattern.rowstart[chunk_row+1]
                                                * cols->chunk_size
                                                * cols->chunk_size];
      while (val_ptr != val_end_of_row)
        {
          if ((cols_have_padding == false)
              ||
              (*colnum_ptr != cols->sparsity_pattern.n_cols()-1))
            internal::ChunkSparseMatrix::chunk_vmult_subtract
            (cols->chunk_size,
             val_ptr,
             u.begin() + *colnum_ptr * cols->chunk_size,
             dst_ptr);
          else
            // we're at a chunk column that has padding
            for (size_type r=0; r<cols->chunk_size; ++r)
              for (size_type c=0; c<n() % cols->chunk_size; ++c)
                dst(chunk_row * cols->chunk_size + r)
                -= (val_ptr[r*cols->chunk_size + c] *
                    u(*colnum_ptr * cols->chunk_size + c));

          ++colnum_ptr;
          val_ptr += cols->chunk_size * cols->chunk_size;
        }


      dst_ptr += cols->chunk_size;
    }

  // now deal with last chunk row if necessary
  if (rows_have_padding)
    {
      const size_type chunk_row = n_chunk_rows - 1;

      const number *const val_end_of_row = &val[cols->sparsity_pattern.rowstart[chunk_row+1]
                                                * cols->chunk_size
                                                * cols->chunk_size];
      while (val_ptr != val_end_of_row)
        {
          if ((cols_have_padding == false)
              ||
              (*colnum_ptr != cols->sparsity_pattern.n_cols()-1))
            {
              // we're at a chunk row but not column that has padding
              for (size_type r=0; r<m() % cols->chunk_size; ++r)
                for (size_type c=0; c<cols->chunk_size; ++c)
                  dst(chunk_row * cols->chunk_size + r)
                  -= (val_ptr[r*cols->chunk_size + c] *
                      u(*colnum_ptr * cols->chunk_size + c));
            }
          else
            // we're at a chunk row and column that has padding
            for (size_type r=0; r<m() % cols->chunk_size; ++r)
              for (size_type c=0; c<n() % cols->chunk_size; ++c)
                dst(chunk_row * cols->chunk_size + r)
                -= (val_ptr[r*cols->chunk_size + c] *
                    u(*colnum_ptr * cols->chunk_size + c));

          ++colnum_ptr;
          val_ptr += cols->chunk_size * cols->chunk_size;
        }


      dst_ptr += cols->chunk_size;
    }

  // finally compute the norm
  return dst.l2_norm();
}



template <typename number>
template <typename somenumber>
void
ChunkSparseMatrix<number>::precondition_Jacobi (Vector<somenumber>       &dst,
                                                const Vector<somenumber> &src,
                                                const number              /*om*/) const
{
  Assert (cols != 0, ExcNotInitialized());
  Assert (val != 0, ExcNotInitialized());
  Assert (m() == n(), ExcMessage("This operation is only valid on square matrices."));

  Assert (dst.size() == n(), ExcDimensionMismatch (dst.size(), n()));
  Assert (src.size() == n(), ExcDimensionMismatch (src.size(), n()));

  Assert (false, ExcNotImplemented());
}



template <typename number>
template <typename somenumber>
void
ChunkSparseMatrix<number>::precondition_SSOR (Vector<somenumber>       &dst,
                                              const Vector<somenumber> &src,
                                              const number              /*om*/) const
{
  // to understand how this function works you may want to take a look at the
  // CVS archives to see the original version which is much clearer...
  Assert (cols != 0, ExcNotInitialized());
  Assert (val != 0, ExcNotInitialized());
  Assert (m() == n(), ExcMessage("This operation is only valid on square matrices."));

  Assert (dst.size() == n(), ExcDimensionMismatch (dst.size(), n()));
  Assert (src.size() == n(), ExcDimensionMismatch (src.size(), n()));

  Assert (false, ExcNotImplemented());
}


template <typename number>
template <typename somenumber>
void
ChunkSparseMatrix<number>::precondition_SOR (Vector<somenumber> &dst,
                                             const Vector<somenumber> &src,
                                             const number om) const
{
  Assert (cols != 0, ExcNotInitialized());
  Assert (val != 0, ExcNotInitialized());
  Assert (m() == n(), ExcMessage("This operation is only valid on square matrices."));


  dst = src;
  SOR(dst,om);
}


template <typename number>
template <typename somenumber>
void
ChunkSparseMatrix<number>::precondition_TSOR (Vector<somenumber> &dst,
                                              const Vector<somenumber> &src,
                                              const number om) const
{
  Assert (cols != 0, ExcNotInitialized());
  Assert (val != 0, ExcNotInitialized());
  Assert (m() == n(), ExcMessage("This operation is only valid on square matrices."));


  dst = src;
  TSOR(dst,om);
}


template <typename number>
template <typename somenumber>
void
ChunkSparseMatrix<number>::SOR (Vector<somenumber> &dst,
                                const number /*om*/) const
{
  Assert (cols != 0, ExcNotInitialized());
  Assert (val != 0, ExcNotInitialized());
  Assert (m() == n(), ExcMessage("This operation is only valid on square matrices."));
  Assert (m() == dst.size(), ExcDimensionMismatch(m(),dst.size()));

  Assert (false, ExcNotImplemented());
}


template <typename number>
template <typename somenumber>
void
ChunkSparseMatrix<number>::TSOR (Vector<somenumber> &dst,
                                 const number /*om*/) const
{
  Assert (cols != 0, ExcNotInitialized());
  Assert (val != 0, ExcNotInitialized());
  Assert (m() == n(), ExcMessage("This operation is only valid on square matrices."));
  Assert (m() == dst.size(), ExcDimensionMismatch(m(),dst.size()));

  Assert (false, ExcNotImplemented());
}


template <typename number>
template <typename somenumber>
void
ChunkSparseMatrix<number>::PSOR (Vector<somenumber> &dst,
                                 const std::vector<size_type> &permutation,
                                 const std::vector<size_type> &inverse_permutation,
                                 const number /*om*/) const
{
  Assert (cols != 0, ExcNotInitialized());
  Assert (val != 0, ExcNotInitialized());
  Assert (m() == n(), ExcMessage("This operation is only valid on square matrices."));

  Assert (m() == dst.size(), ExcDimensionMismatch(m(),dst.size()));
  Assert (m() == permutation.size(),
          ExcDimensionMismatch(m(), permutation.size()));
  Assert (m() == inverse_permutation.size(),
          ExcDimensionMismatch(m(), inverse_permutation.size()));

  Assert (false, ExcNotImplemented());
}


template <typename number>
template <typename somenumber>
void
ChunkSparseMatrix<number>::TPSOR (Vector<somenumber> &dst,
                                  const std::vector<size_type> &permutation,
                                  const std::vector<size_type> &inverse_permutation,
                                  const number /*om*/) const
{
  Assert (cols != 0, ExcNotInitialized());
  Assert (val != 0, ExcNotInitialized());
  Assert (m() == n(), ExcMessage("This operation is only valid on square matrices."));

  Assert (m() == dst.size(), ExcDimensionMismatch(m(),dst.size()));
  Assert (m() == permutation.size(),
          ExcDimensionMismatch(m(), permutation.size()));
  Assert (m() == inverse_permutation.size(),
          ExcDimensionMismatch(m(), inverse_permutation.size()));

  Assert (false, ExcNotImplemented());
}



template <typename number>
template <typename somenumber>
void
ChunkSparseMatrix<number>::SOR_step (Vector<somenumber> &v,
                                     const Vector<somenumber> &b,
                                     const number        /*om*/) const
{
  Assert (cols != 0, ExcNotInitialized());
  Assert (val != 0, ExcNotInitialized());
  Assert (m() == n(), ExcMessage("This operation is only valid on square matrices."));

  Assert (m() == v.size(), ExcDimensionMismatch(m(),v.size()));
  Assert (m() == b.size(), ExcDimensionMismatch(m(),b.size()));

  Assert (false, ExcNotImplemented());
}



template <typename number>
template <typename somenumber>
void
ChunkSparseMatrix<number>::TSOR_step (Vector<somenumber> &v,
                                      const Vector<somenumber> &b,
                                      const number        /*om*/) const
{
  Assert (cols != 0, ExcNotInitialized());
  Assert (val != 0, ExcNotInitialized());
  Assert (m() == n(), ExcMessage("This operation is only valid on square matrices."));

  Assert (m() == v.size(), ExcDimensionMismatch(m(),v.size()));
  Assert (m() == b.size(), ExcDimensionMismatch(m(),b.size()));

  Assert (false, ExcNotImplemented());
}



template <typename number>
template <typename somenumber>
void
ChunkSparseMatrix<number>::SSOR_step (Vector<somenumber> &v,
                                      const Vector<somenumber> &b,
                                      const number        om) const
{
  SOR_step(v,b,om);
  TSOR_step(v,b,om);
}



template <typename number>
template <typename somenumber>
void
ChunkSparseMatrix<number>::SSOR (Vector<somenumber> &dst,
                                 const number /*om*/) const
{
  Assert (cols != 0, ExcNotInitialized());
  Assert (val != 0, ExcNotInitialized());
  Assert (m() == n(), ExcMessage("This operation is only valid on square matrices."));

  Assert (m() == dst.size(), ExcDimensionMismatch(m(),dst.size()));

  Assert (false, ExcNotImplemented());
}



template <typename number>
void ChunkSparseMatrix<number>::print (std::ostream &out) const
{
  AssertThrow (out, ExcIO());

  Assert (cols != 0, ExcNotInitialized());
  Assert (val != 0, ExcNotInitialized());

  Assert (false, ExcNotImplemented());

  AssertThrow (out, ExcIO());
}


template <typename number>
void ChunkSparseMatrix<number>::print_formatted (std::ostream &out,
                                                 const unsigned int precision,
                                                 const bool scientific,
                                                 const unsigned int width_,
                                                 const char *zero_string,
                                                 const double denominator) const
{
  AssertThrow (out, ExcIO());

  Assert (cols != 0, ExcNotInitialized());
  Assert (val != 0, ExcNotInitialized());

  unsigned int width = width_;

  Assert (false, ExcNotImplemented());

  std::ios::fmtflags old_flags = out.flags();
  unsigned int old_precision = out.precision (precision);

  if (scientific)
    {
      out.setf (std::ios::scientific, std::ios::floatfield);
      if (!width)
        width = precision+7;
    }
  else
    {
      out.setf (std::ios::fixed, std::ios::floatfield);
      if (!width)
        width = precision+2;
    }

  for (size_type i=0; i<m(); ++i)
    {
      for (size_type j=0; j<n(); ++j)
        if (cols->sparsity_pattern(i,j) != SparsityPattern::invalid_entry)
          out << std::setw(width)
              << val[cols->sparsity_pattern(i,j)] * denominator << ' ';
        else
          out << std::setw(width) << zero_string << ' ';
      out << std::endl;
    };
  AssertThrow (out, ExcIO());

  // reset output format
  out.precision(old_precision);
  out.flags (old_flags);
}



template <typename number>
void ChunkSparseMatrix<number>::print_pattern (std::ostream &out,
                                               const double threshold) const
{
  AssertThrow (out, ExcIO());

  Assert (cols != 0, ExcNotInitialized());
  Assert (val != 0, ExcNotInitialized());

  const size_type chunk_size = cols->get_chunk_size();

  // loop over all chunk rows and columns, and each time we find something
  // repeat it chunk_size times in both directions
  for (size_type i=0; i<cols->sparsity_pattern.n_rows(); ++i)
    {
      for (size_type d=0; d<chunk_size; ++d)
        for (size_type j=0; j<cols->sparsity_pattern.n_cols(); ++j)
          if (cols->sparsity_pattern(i,j) == SparsityPattern::invalid_entry)
            {
              for (size_type e=0; e<chunk_size; ++e)
                out << '.';
            }
          else if (std::fabs(val[cols->sparsity_pattern(i,j)]) > threshold)
            {
              for (size_type e=0; e<chunk_size; ++e)
                out << '*';
            }
          else
            {
              for (size_type e=0; e<chunk_size; ++e)
                out << ':';
            }
      out << std::endl;
    }

  AssertThrow (out, ExcIO());
}



template <typename number>
void
ChunkSparseMatrix<number>::block_write (std::ostream &out) const
{
  AssertThrow (out, ExcIO());

  // first the simple objects, bracketed in [...]
  out << '[' << max_len << "][";
  // then write out real data
  out.write (reinterpret_cast<const char *>(&val[0]),
             reinterpret_cast<const char *>(&val[max_len])
             - reinterpret_cast<const char *>(&val[0]));
  out << ']';

  AssertThrow (out, ExcIO());
}



template <typename number>
void
ChunkSparseMatrix<number>::block_read (std::istream &in)
{
  AssertThrow (in, ExcIO());

  char c;

  // first read in simple data
  in >> c;
  AssertThrow (c == '[', ExcIO());
  in >> max_len;

  in >> c;
  AssertThrow (c == ']', ExcIO());
  in >> c;
  AssertThrow (c == '[', ExcIO());

  // reallocate space
  delete[] val;
  val  = new number[max_len];

  // then read data
  in.read (reinterpret_cast<char *>(&val[0]),
           reinterpret_cast<char *>(&val[max_len])
           - reinterpret_cast<char *>(&val[0]));
  in >> c;
  AssertThrow (c == ']', ExcIO());
}



template <typename number>
std::size_t
ChunkSparseMatrix<number>::memory_consumption () const
{
  return sizeof(*this) + max_len*sizeof(number);
}


DEAL_II_NAMESPACE_CLOSE

#endif
