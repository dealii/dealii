// ---------------------------------------------------------------------
//
// Copyright (C) 2008 - 2013 by the deal.II authors
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


#include <deal.II/lac/chunk_sparsity_pattern.h>
#include <deal.II/lac/compressed_sparsity_pattern.h>
#include <deal.II/lac/compressed_set_sparsity_pattern.h>
#include <deal.II/lac/compressed_simple_sparsity_pattern.h>
#include <deal.II/lac/full_matrix.h>


DEAL_II_NAMESPACE_OPEN


ChunkSparsityPattern::ChunkSparsityPattern ()
{
  reinit (0,0,0,1);
}



ChunkSparsityPattern::ChunkSparsityPattern (const ChunkSparsityPattern &s)
  :
  Subscriptor(),
  chunk_size (s.chunk_size),
  sparsity_pattern(s.sparsity_pattern)
{
  Assert (s.rows == 0, ExcInvalidConstructorCall());
  Assert (s.cols == 0, ExcInvalidConstructorCall());

  reinit (0,0,0,0);
}



ChunkSparsityPattern::ChunkSparsityPattern (const size_type m,
                                            const size_type n,
                                            const size_type max_per_row,
                                            const size_type chunk_size,
                                            const bool)
{
  Assert (chunk_size > 0, ExcInvalidNumber (chunk_size));

  reinit (m,n,max_per_row, chunk_size);
}


ChunkSparsityPattern::ChunkSparsityPattern (const size_type m,
                                            const size_type n,
                                            const size_type max_per_row,
                                            const size_type chunk_size)
{
  Assert (chunk_size > 0, ExcInvalidNumber (chunk_size));

  reinit (m,n,max_per_row, chunk_size);
}



ChunkSparsityPattern::ChunkSparsityPattern (
  const size_type m,
  const size_type n,
  const std::vector<size_type> &row_lengths,
  const size_type chunk_size,
  const bool)
{
  Assert (chunk_size > 0, ExcInvalidNumber (chunk_size));

  reinit (m, n, row_lengths, chunk_size);
}



ChunkSparsityPattern::ChunkSparsityPattern (
  const size_type m,
  const size_type n,
  const std::vector<size_type> &row_lengths,
  const size_type chunk_size)
{
  Assert (chunk_size > 0, ExcInvalidNumber (chunk_size));

  reinit (m, n, row_lengths, chunk_size);
}



ChunkSparsityPattern::ChunkSparsityPattern (const size_type n,
                                            const size_type max_per_row,
                                            const size_type chunk_size)
{
  reinit (n, n, max_per_row, chunk_size);
}



ChunkSparsityPattern::ChunkSparsityPattern (
  const size_type m,
  const std::vector<size_type> &row_lengths,
  const size_type chunk_size,
  const bool)
{
  Assert (chunk_size > 0, ExcInvalidNumber (chunk_size));

  reinit (m, m, row_lengths, chunk_size);
}



ChunkSparsityPattern::ChunkSparsityPattern (
  const size_type                m,
  const std::vector<size_type > &row_lengths,
  const size_type chunk_size)
{
  Assert (chunk_size > 0, ExcInvalidNumber (chunk_size));

  reinit (m, m, row_lengths, chunk_size);
}



ChunkSparsityPattern::~ChunkSparsityPattern ()
{}



ChunkSparsityPattern &
ChunkSparsityPattern::operator = (const ChunkSparsityPattern &s)
{
  Assert (s.rows == 0, ExcInvalidConstructorCall());
  Assert (s.cols == 0, ExcInvalidConstructorCall());

  // perform the checks in the underlying object as well
  sparsity_pattern = s.sparsity_pattern;

  return *this;
}



void
ChunkSparsityPattern::reinit (const size_type m,
                              const size_type n,
                              const size_type max_per_row,
                              const size_type chunk_size,
                              const bool)
{
  reinit (m, n, max_per_row, chunk_size);
}


void
ChunkSparsityPattern::reinit (const size_type m,
                              const size_type n,
                              const size_type max_per_row,
                              const size_type chunk_size)
{
  Assert (chunk_size > 0, ExcInvalidNumber (chunk_size));

  // simply map this function to the other @p{reinit} function
  const std::vector<size_type> row_lengths (m, max_per_row);
  reinit (m, n, row_lengths, chunk_size);
}



void
ChunkSparsityPattern::reinit (
  const size_type m,
  const size_type n,
  const VectorSlice<const std::vector<size_type> > &row_lengths,
  const size_type chunk_size,
  const bool)
{
  reinit (m, n, row_lengths, chunk_size);
}


void
ChunkSparsityPattern::reinit (
  const size_type m,
  const size_type n,
  const VectorSlice<const std::vector<size_type> > &row_lengths,
  const size_type chunk_size)
{
  Assert (row_lengths.size() == m, ExcInvalidNumber (m));
  Assert (chunk_size > 0, ExcInvalidNumber (chunk_size));

  rows = m;
  cols = n;

  this->chunk_size = chunk_size;

  // pass down to the necessary information to the underlying object. we need
  // to calculate how many chunks we need: we need to round up (m/chunk_size)
  // and (n/chunk_size). rounding up in integer arithmetic equals
  // ((m+chunk_size-1)/chunk_size):
  const size_type m_chunks = (m+chunk_size-1) / chunk_size,
                  n_chunks = (n+chunk_size-1) / chunk_size;

  // compute the maximum number of chunks in each row. the passed array
  // denotes the number of entries in each row of the big matrix -- in the
  // worst case, these are all in independent chunks, so we have to calculate
  // it as follows (as an example: let chunk_size==2, row_lengths={2,2,...},
  // and entries in row zero at columns {0,2} and for row one at {4,6} -->
  // we'll need 4 chunks for the first chunk row!) :
  std::vector<unsigned int> chunk_row_lengths (m_chunks, 0);
  for (size_type i=0; i<m; ++i)
    chunk_row_lengths[i/chunk_size] += row_lengths[i];

  // for the case that the reduced sparsity pattern optimizes the diagonal but
  // the actual sparsity pattern does not, need to take one more entry in the
  // row to fit the user-required entry
  if (m != n && m_chunks == n_chunks)
    for (unsigned int i=0; i<m_chunks; ++i)
      ++chunk_row_lengths[i];

  sparsity_pattern.reinit (m_chunks,
                           n_chunks,
                           chunk_row_lengths);
}



void
ChunkSparsityPattern::compress ()
{
  sparsity_pattern.compress ();
}



template <typename SparsityType>
void
ChunkSparsityPattern::copy_from (const SparsityType &csp,
                                 const size_type     chunk_size,
                                 const bool)
{
  copy_from (csp, chunk_size);
}



namespace internal
{
  namespace
  {
    /**
     * Declare type for container size.
     */
    typedef types::global_dof_index size_type;

    // distinguish between compressed sparsity types that define row_begin()
    // and SparsityPattern that uses begin() as iterator type
    template <typename Sparsity>
    void copy_row (const Sparsity                  &csp,
                   const size_type                  row,
                   const unsigned int               chunk_size,
                   CompressedSimpleSparsityPattern &dst)
    {
      typename Sparsity::row_iterator col_num = csp.row_begin (row);
      const size_type reduced_row = row/chunk_size;
      for (; col_num != csp.row_end (row); ++col_num)
        dst.add (reduced_row, *col_num/chunk_size);
    }

    void copy_row (const SparsityPattern           &csp,
                   const size_type                  row,
                   const unsigned int               chunk_size,
                   CompressedSimpleSparsityPattern &dst)
    {
      SparsityPattern::iterator col_num = csp.begin (row);
      const size_type reduced_row = row/chunk_size;
      for (; col_num != csp.end (row); ++col_num)
        dst.add (reduced_row, col_num->column()/chunk_size);
    }
  }
}


template <typename SparsityType>
void
ChunkSparsityPattern::copy_from (const SparsityType &csp,
                                 const size_type     chunk_size)
{
  Assert (chunk_size > 0, ExcInvalidNumber (chunk_size));
  this->chunk_size = chunk_size;
  rows = csp.n_rows();
  cols = csp.n_cols();

  // simple case: just use the other sparsity pattern
  if (chunk_size == 1)
    {
      sparsity_pattern.copy_from (csp);
      return;
    }

  // create a temporary compressed sparsity pattern that collects all entries
  // from the input sparsity pattern and then initialize the underlying small
  // sparsity pattern
  const size_type m_chunks = (csp.n_rows()+chunk_size-1) / chunk_size,
                  n_chunks = (csp.n_cols()+chunk_size-1) / chunk_size;
  CompressedSimpleSparsityPattern temporary_sp(m_chunks, n_chunks);

  for (size_type row = 0; row<csp.n_rows(); ++row)
    internal::copy_row(csp, row, chunk_size, temporary_sp);

  sparsity_pattern.copy_from (temporary_sp);
}




template <typename number>
void ChunkSparsityPattern::copy_from (const FullMatrix<number> &matrix,
                                      const size_type chunk_size,
                                      const bool)
{
  copy_from (matrix, chunk_size);
}


template <typename number>
void ChunkSparsityPattern::copy_from (const FullMatrix<number> &matrix,
                                      const size_type chunk_size)
{
  Assert (chunk_size > 0, ExcInvalidNumber (chunk_size));

  // count number of entries per row, then initialize the underlying sparsity
  // pattern. remember to also allocate space for the diagonal entry (if that
  // hasn't happened yet) if m==n since we always allocate that for diagonal
  // matrices
  std::vector<size_type> entries_per_row (matrix.m(), 0);
  for (size_type row=0; row<matrix.m(); ++row)
    {
      for (size_type col=0; col<matrix.n(); ++col)
        if (matrix(row,col) != 0)
          ++entries_per_row[row];

      if ((matrix.m() == matrix.n())
          &&
          (matrix(row,row) == 0))
        ++entries_per_row[row];
    }

  reinit (matrix.m(), matrix.n(),
          entries_per_row,
          chunk_size);

  // then actually fill it
  for (size_type row=0; row<matrix.m(); ++row)
    for (size_type col=0; col<matrix.n(); ++col)
      if (matrix(row,col) != 0)
        add (row,col);

  // finally compress
  compress ();
}


void
ChunkSparsityPattern::reinit (
  const size_type m,
  const size_type n,
  const std::vector<size_type> &row_lengths,
  const size_type chunk_size,
  const bool)
{
  reinit (m, n, row_lengths, chunk_size);
}


void
ChunkSparsityPattern::reinit (
  const size_type m,
  const size_type n,
  const std::vector<size_type> &row_lengths,
  const size_type chunk_size)
{
  Assert (chunk_size > 0, ExcInvalidNumber (chunk_size));

  reinit(m, n, make_slice(row_lengths), chunk_size);
}



namespace internal
{
  namespace
  {
    template <typename SP>
    void copy_sparsity (const SP        &src,
                        SparsityPattern &dst)
    {
      dst.copy_from(src);
    }

    void copy_sparsity (const SparsityPattern &src,
                        SparsityPattern       &dst)
    {
      dst = src;
    }
  }
}



template <typename Sparsity>
void
ChunkSparsityPattern::create_from
(const unsigned int m,
 const unsigned int n,
 const Sparsity    &sparsity_pattern_for_chunks,
 const unsigned int chunk_size_in,
 const bool)
{
  Assert (m > (sparsity_pattern_for_chunks.n_rows()-1) * chunk_size_in &&
          m <= sparsity_pattern_for_chunks.n_rows() * chunk_size_in,
          ExcMessage("Number of rows m is not compatible with chunk size "
                     "and number of rows in sparsity pattern for the chunks."));
  Assert (n > (sparsity_pattern_for_chunks.n_cols()-1) * chunk_size_in &&
          n <= sparsity_pattern_for_chunks.n_cols() * chunk_size_in,
          ExcMessage("Number of columns m is not compatible with chunk size "
                     "and number of columns in sparsity pattern for the chunks."));

  internal::copy_sparsity(sparsity_pattern_for_chunks, sparsity_pattern);
  chunk_size = chunk_size_in;
  rows = m;
  cols = n;
}



bool
ChunkSparsityPattern::empty () const
{
  return sparsity_pattern.empty();
}



ChunkSparsityPattern::size_type
ChunkSparsityPattern::max_entries_per_row () const
{
  return sparsity_pattern.max_entries_per_row() * chunk_size;
}



void
ChunkSparsityPattern::add (const size_type i,
                           const size_type j)
{
  Assert (i<rows, ExcInvalidIndex(i,rows));
  Assert (j<cols, ExcInvalidIndex(j,cols));

  sparsity_pattern.add (i/chunk_size, j/chunk_size);
}


bool
ChunkSparsityPattern::exists (const size_type i,
                              const size_type j) const
{
  Assert (i<rows, ExcIndexRange(i,0,rows));
  Assert (j<cols, ExcIndexRange(j,0,cols));

  return sparsity_pattern.exists (i/chunk_size,
                                  j/chunk_size);
}



void
ChunkSparsityPattern::symmetrize ()
{
  // matrix must be square. note that the for some matrix sizes, the current
  // sparsity pattern may not be square even if the underlying sparsity
  // pattern is (e.g. a 10x11 matrix with chunk_size 4)
  Assert (rows==cols, ExcNotQuadratic());

  sparsity_pattern.symmetrize ();
}



ChunkSparsityPattern::size_type
ChunkSparsityPattern::row_length (const size_type i) const
{
  Assert (i<rows, ExcIndexRange(i,0,rows));

  // find out if we did padding and if this row is affected by it
  if (n_cols() % chunk_size == 0)
    return sparsity_pattern.row_length (i/chunk_size) * chunk_size;
  else
    // if columns don't align, then just iterate over all chunks and see
    // what this leads to
    {
      SparsityPattern::const_iterator p = sparsity_pattern.begin(i/chunk_size),
                                      end = sparsity_pattern.end(i/chunk_size);
      unsigned int n = 0;
      for ( ; p != end; ++p)
        if (p->column() != sparsity_pattern.n_cols() - 1)
          n += chunk_size;
        else
          n += (n_cols() % chunk_size);
      return n;
    }
}



ChunkSparsityPattern::size_type
ChunkSparsityPattern::n_nonzero_elements () const
{
  if ((n_rows() % chunk_size == 0)
      &&
      (n_cols() % chunk_size == 0))
    return (sparsity_pattern.n_nonzero_elements() *
            chunk_size *
            chunk_size);
  else
    // some of the chunks reach beyond the extent of this matrix. this
    // requires a somewhat more complicated computations, in particular if the
    // columns don't align
    {
      if ((n_rows() % chunk_size != 0)
          &&
          (n_cols() % chunk_size == 0))
        {
          // columns align with chunks, but not rows
          size_type n = sparsity_pattern.n_nonzero_elements() *
                        chunk_size *
                        chunk_size;
          n -= (sparsity_pattern.n_rows() * chunk_size - n_rows()) *
               sparsity_pattern.row_length(sparsity_pattern.n_rows()-1) *
               chunk_size;
          return n;
        }

      else
        {
          // if columns don't align, then just iterate over all chunks and see
          // what this leads to. follow the advice in the documentation of the
          // sparsity pattern iterators to do the loop over individual rows,
          // rather than all elements
          size_type n = 0;

          for (size_type row = 0; row < sparsity_pattern.n_rows(); ++row)
            {
              SparsityPattern::const_iterator p = sparsity_pattern.begin(row);
              for (; p!=sparsity_pattern.end(row); ++p)
                if ((row != sparsity_pattern.n_rows() - 1)
                    &&
                    (p->column() != sparsity_pattern.n_cols() - 1))
                  n += chunk_size * chunk_size;
                else if ((row == sparsity_pattern.n_rows() - 1)
                         &&
                         (p->column() != sparsity_pattern.n_cols() - 1))
                  // last chunk row, but not last chunk column. only a smaller
                  // number (n_rows % chunk_size) of rows actually exist
                  n += (n_rows() % chunk_size) * chunk_size;
                else if ((row != sparsity_pattern.n_rows() - 1)
                         &&
                         (p->column() == sparsity_pattern.n_cols() - 1))
                  // last chunk column, but not row
                  n += (n_cols() % chunk_size) * chunk_size;
                else
                  // bottom right chunk
                  n += (n_cols() % chunk_size) *
                       (n_rows() % chunk_size);
            }

          return n;
        }
    }
}



void
ChunkSparsityPattern::print (std::ostream &out) const
{
  Assert ((sparsity_pattern.rowstart!=0) && (sparsity_pattern.colnums!=0),
          ExcEmptyObject());

  AssertThrow (out, ExcIO());

  for (size_type i=0; i<sparsity_pattern.rows; ++i)
    for (size_type d=0;
         (d<chunk_size) && (i*chunk_size + d < n_rows());
         ++d)
      {
        out << '[' << i *chunk_size+d;
        for (size_type j=sparsity_pattern.rowstart[i];
             j<sparsity_pattern.rowstart[i+1]; ++j)
          if (sparsity_pattern.colnums[j] != sparsity_pattern.invalid_entry)
            for (size_type e=0;
                 ((e<chunk_size) &&
                  (sparsity_pattern.colnums[j]*chunk_size + e < n_cols()));
                 ++e)
              out << ',' << sparsity_pattern.colnums[j]*chunk_size+e;
        out << ']' << std::endl;
      }

  AssertThrow (out, ExcIO());
}



void
ChunkSparsityPattern::print_gnuplot (std::ostream &out) const
{
  Assert ((sparsity_pattern.rowstart!=0) &&
          (sparsity_pattern.colnums!=0), ExcEmptyObject());

  AssertThrow (out, ExcIO());

  // for each entry in the underlying sparsity pattern, repeat everything
  // chunk_size x chunk_size times
  for (size_type i=0; i<sparsity_pattern.rows; ++i)
    for (size_type j=sparsity_pattern.rowstart[i];
         j<sparsity_pattern.rowstart[i+1]; ++j)
      if (sparsity_pattern.colnums[j] != sparsity_pattern.invalid_entry)
        for (size_type d=0;
             ((d<chunk_size) &&
              (sparsity_pattern.colnums[j]*chunk_size+d < n_cols()));
             ++d)
          for (size_type e=0;
               (e<chunk_size) && (i*chunk_size + e < n_rows());
               ++e)
            // while matrix entries are usually written (i,j), with i vertical
            // and j horizontal, gnuplot output is x-y, that is we have to
            // exchange the order of output
            out << sparsity_pattern.colnums[j]*chunk_size+d << " "
                << -static_cast<signed int>(i*chunk_size+e)
                << std::endl;

  AssertThrow (out, ExcIO());
}



ChunkSparsityPattern::size_type
ChunkSparsityPattern::bandwidth () const
{
  // calculate the bandwidth from that of the underlying sparsity
  // pattern. note that even if the bandwidth of that is zero, then the
  // bandwidth of the chunky pattern is chunk_size-1, if it is 1 then the
  // chunky pattern has chunk_size+(chunk_size-1), etc
  //
  // we'll cut it off at max(n(),m())
  return std::min (sparsity_pattern.bandwidth()*chunk_size
                   + (chunk_size-1),
                   std::max(n_rows(), n_cols()));
}



bool
ChunkSparsityPattern::stores_only_added_elements () const
{
  if (chunk_size == 1)
    return sparsity_pattern.stores_only_added_elements ();
  else
    return false;
}



void
ChunkSparsityPattern::block_write (std::ostream &out) const
{
  AssertThrow (out, ExcIO());

  // first the simple objects, bracketed in [...]
  out << '['
      << rows << ' '
      << cols << ' '
      << chunk_size << ' '
      << "][";
  // then the underlying sparsity pattern
  sparsity_pattern.block_write (out);
  out << ']';

  AssertThrow (out, ExcIO());
}



void
ChunkSparsityPattern::block_read (std::istream &in)
{
  AssertThrow (in, ExcIO());

  char c;

  // first read in simple data
  in >> c;
  AssertThrow (c == '[', ExcIO());
  in >> rows
     >> cols
     >> chunk_size;

  in >> c;
  AssertThrow (c == ']', ExcIO());
  in >> c;
  AssertThrow (c == '[', ExcIO());

  // then read the underlying sparsity pattern
  sparsity_pattern.block_read (in);

  in >> c;
  AssertThrow (c == ']', ExcIO());
}



std::size_t
ChunkSparsityPattern::memory_consumption () const
{
  return (sizeof(*this) +
          sparsity_pattern.memory_consumption());
}



// explicit instantiations
template
void ChunkSparsityPattern::copy_from<SparsityPattern> (const SparsityPattern &,
                                                       const size_type,
                                                       const bool);
template
void ChunkSparsityPattern::copy_from<CompressedSparsityPattern> (const CompressedSparsityPattern &,
    const size_type ,
    const bool);
template
void ChunkSparsityPattern::copy_from<CompressedSetSparsityPattern> (const CompressedSetSparsityPattern &,
    const size_type ,
    const bool);
template
void ChunkSparsityPattern::copy_from<CompressedSimpleSparsityPattern> (const CompressedSimpleSparsityPattern &,
    const size_type ,
    const bool);
template
void ChunkSparsityPattern::copy_from<CompressedSparsityPattern> (const CompressedSparsityPattern &,
    const size_type);
template
void ChunkSparsityPattern::copy_from<CompressedSetSparsityPattern> (const CompressedSetSparsityPattern &,
    const size_type);
template
void ChunkSparsityPattern::copy_from<CompressedSimpleSparsityPattern> (const CompressedSimpleSparsityPattern &,
    const size_type);
template
void ChunkSparsityPattern::create_from<SparsityPattern>
(const unsigned int,
 const unsigned int,
 const SparsityPattern &,
 const unsigned int,
 const bool);
template
void ChunkSparsityPattern::create_from<CompressedSparsityPattern>
(const unsigned int,
 const unsigned int,
 const CompressedSparsityPattern &,
 const unsigned int,
 const bool);
template
void ChunkSparsityPattern::create_from<CompressedSetSparsityPattern>
(const unsigned int,
 const unsigned int,
 const CompressedSetSparsityPattern &,
 const unsigned int,
 const bool);
template
void ChunkSparsityPattern::create_from<CompressedSimpleSparsityPattern>
(const unsigned int,
 const unsigned int,
 const CompressedSimpleSparsityPattern &,
 const unsigned int,
 const bool);
template
void ChunkSparsityPattern::copy_from<float> (const FullMatrix<float> &,
                                             const size_type ,
                                             const bool);
template
void ChunkSparsityPattern::copy_from<double> (const FullMatrix<double> &,
                                              const size_type ,
                                              const bool);
template
void ChunkSparsityPattern::copy_from<float> (const FullMatrix<float> &,
                                             const size_type);
template
void ChunkSparsityPattern::copy_from<double> (const FullMatrix<double> &,
                                              const size_type);

DEAL_II_NAMESPACE_CLOSE
