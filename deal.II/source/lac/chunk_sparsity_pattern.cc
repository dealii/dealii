//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2008, 2011 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------

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

  reinit (0,0,0,0, false);
}



ChunkSparsityPattern::ChunkSparsityPattern (const unsigned int m,
					    const unsigned int n,
					    const unsigned int max_per_row,
					    const unsigned int chunk_size,
					    const bool optimize_diag)
{
  Assert (chunk_size > 0, ExcInvalidNumber (chunk_size));

  reinit (m,n,max_per_row, chunk_size, optimize_diag);
}



ChunkSparsityPattern::ChunkSparsityPattern (
  const unsigned int m,
  const unsigned int n,
  const std::vector<unsigned int>& row_lengths,
  const unsigned int chunk_size,
  const bool optimize_diag)
{
  Assert (chunk_size > 0, ExcInvalidNumber (chunk_size));

  reinit (m, n, row_lengths, chunk_size, optimize_diag);
}



ChunkSparsityPattern::ChunkSparsityPattern (const unsigned int n,
					    const unsigned int max_per_row,
					    const unsigned int chunk_size)
{
  reinit (n, n, max_per_row, chunk_size, true);
}



ChunkSparsityPattern::ChunkSparsityPattern (
  const unsigned int               m,
  const std::vector<unsigned int>& row_lengths,
  const unsigned int chunk_size,
  const bool optimize_diag)
{
  Assert (chunk_size > 0, ExcInvalidNumber (chunk_size));

  reinit (m, m, row_lengths, chunk_size, optimize_diag);
}



ChunkSparsityPattern::~ChunkSparsityPattern ()
{}



ChunkSparsityPattern &
ChunkSparsityPattern::operator = (const ChunkSparsityPattern &s)
{
  Assert (s.rows == 0, ExcInvalidConstructorCall());
  Assert (s.cols == 0, ExcInvalidConstructorCall());

				   // perform the checks in the underlying
				   // object as well
  sparsity_pattern = s.sparsity_pattern;

  return *this;
}



void
ChunkSparsityPattern::reinit (const unsigned int m,
			      const unsigned int n,
			      const unsigned int max_per_row,
			      const unsigned int chunk_size,
			      const bool optimize_diag)
{
  Assert (chunk_size > 0, ExcInvalidNumber (chunk_size));

				   // simply map this function to the
				   // other @p{reinit} function
  const std::vector<unsigned int> row_lengths (m, max_per_row);
  reinit (m, n, row_lengths, chunk_size, optimize_diag);
}



void
ChunkSparsityPattern::reinit (
  const unsigned int m,
  const unsigned int n,
  const VectorSlice<const std::vector<unsigned int> >&row_lengths,
  const unsigned int chunk_size,
  const bool optimize_diag)
{
  Assert (row_lengths.size() == m, ExcInvalidNumber (m));
  Assert (chunk_size > 0, ExcInvalidNumber (chunk_size));

  rows = m;
  cols = n;

  this->chunk_size = chunk_size;

				   // pass down to the necessary information
				   // to the underlying object. we need to
				   // calculate how many chunks we need: we
				   // need to round up (m/chunk_size) and
				   // (n/chunk_size). rounding up in integer
				   // arithmetic equals
				   // ((m+chunk_size-1)/chunk_size):
  const unsigned int m_chunks = (m+chunk_size-1) / chunk_size,
		     n_chunks = (n+chunk_size-1) / chunk_size;

				   // compute the maximum number of chunks in
				   // each row. the passed array denotes the
				   // number of entries in each row of the big
				   // matrix -- in the worst case, these are
				   // all in independent chunks, so we have to
				   // calculate it as follows (as an example:
				   // let chunk_size==2,
				   // row_lengths={2,2,...}, and entries in
				   // row zero at columns {0,2} and for row
				   // one at {4,6} --> we'll need 4 chunks for
				   // the first chunk row!) :
  std::vector<unsigned int> chunk_row_lengths (m_chunks, 0);
  for (unsigned int i=0; i<m; ++i)
    chunk_row_lengths[i/chunk_size] += row_lengths[i];

  sparsity_pattern.reinit (m_chunks,
			   n_chunks,
			   chunk_row_lengths,
			   optimize_diag);
}



void
ChunkSparsityPattern::compress ()
{
  sparsity_pattern.compress ();
}



template <typename SparsityType>
void
ChunkSparsityPattern::copy_from (const SparsityType &csp,
				 const unsigned int  chunk_size,
				 const bool          optimize_diag)
{
  Assert (chunk_size > 0, ExcInvalidNumber (chunk_size));

				   // count number of entries per row, then
				   // initialize the underlying sparsity
				   // pattern
  std::vector<unsigned int> entries_per_row (csp.n_rows(), 0);
  for (unsigned int row = 0; row<csp.n_rows(); ++row)
    entries_per_row[row] = csp.row_length(row);

  reinit (csp.n_rows(), csp.n_cols(),
	  entries_per_row,
	  chunk_size, optimize_diag);

				   // then actually fill it
  for (unsigned int row = 0; row<csp.n_rows(); ++row)
    {
      typename SparsityType::row_iterator col_num = csp.row_begin (row);

      for (; col_num != csp.row_end (row); ++col_num)
	add (row, *col_num);
    }

				   // finally compress
  compress ();
}




template <typename number>
void ChunkSparsityPattern::copy_from (const FullMatrix<number> &matrix,
				      const unsigned int chunk_size,
				      const bool         optimize_diag)
{
  Assert (chunk_size > 0, ExcInvalidNumber (chunk_size));

				   // count number of entries per row, then
				   // initialize the underlying sparsity
				   // pattern
  std::vector<unsigned int> entries_per_row (matrix.m(), 0);
  for (unsigned int row=0; row<matrix.m(); ++row)
    for (unsigned int col=0; col<matrix.n(); ++col)
      if (matrix(row,col) != 0)
	++entries_per_row[row];

  reinit (matrix.m(), matrix.n(),
	  entries_per_row,
	  chunk_size, optimize_diag);

				   // then actually fill it
  for (unsigned int row=0; row<matrix.m(); ++row)
    for (unsigned int col=0; col<matrix.n(); ++col)
      if (matrix(row,col) != 0)
	add (row,col);

				   // finally compress
  compress ();
}


void
ChunkSparsityPattern::reinit (
  const unsigned int m,
  const unsigned int n,
  const std::vector<unsigned int>& row_lengths,
  const unsigned int chunk_size,
  const bool optimize_diag)
{
  Assert (chunk_size > 0, ExcInvalidNumber (chunk_size));

  reinit(m, n, make_slice(row_lengths), chunk_size, optimize_diag);
}



bool
ChunkSparsityPattern::empty () const
{
  return sparsity_pattern.empty();
}



unsigned int
ChunkSparsityPattern::max_entries_per_row () const
{
  return sparsity_pattern.max_entries_per_row() * chunk_size;
}



void
ChunkSparsityPattern::add (const unsigned int i,
			   const unsigned int j)
{
  Assert (i<rows, ExcInvalidIndex(i,rows));
  Assert (j<cols, ExcInvalidIndex(j,cols));

  sparsity_pattern.add (i/chunk_size, j/chunk_size);
}


bool
ChunkSparsityPattern::exists (const unsigned int i,
			      const unsigned int j) const
{
  Assert (i<rows, ExcIndexRange(i,0,rows));
  Assert (j<cols, ExcIndexRange(j,0,cols));

  return sparsity_pattern.exists (i/chunk_size,
				  j/chunk_size);
}



unsigned int
ChunkSparsityPattern::row_length (const unsigned int i) const
{
  Assert (i<rows, ExcIndexRange(i,0,rows));

  return sparsity_pattern.row_length (i/chunk_size) * chunk_size;
}



void
ChunkSparsityPattern::symmetrize ()
{
				   // matrix must be square. note that the for
				   // some matrix sizes, the current sparsity
				   // pattern may not be square even if the
				   // underlying sparsity pattern is (e.g. a
				   // 10x11 matrix with chunk_size 4)
  Assert (rows==cols, ExcNotQuadratic());

  sparsity_pattern.symmetrize ();
}



unsigned int
ChunkSparsityPattern::n_nonzero_elements () const
{
  if ((n_rows() % chunk_size == 0)
      &&
      (n_cols() % chunk_size == 0))
    return (sparsity_pattern.n_nonzero_elements() *
	    chunk_size *
	    chunk_size);
  else
				     // some of the chunks reach beyond the
				     // extent of this matrix. this requires a
				     // somewhat more complicated
				     // computations, in particular if the
				     // columns don't align
    {
      if ((n_rows() % chunk_size != 0)
	  &&
	  (n_cols() % chunk_size == 0))
	{
					   // columns align with chunks, but
					   // not rows
	  unsigned int n = sparsity_pattern.n_nonzero_elements() *
			   chunk_size *
			   chunk_size;
	  n -= (sparsity_pattern.n_rows() * chunk_size - n_rows()) *
	       sparsity_pattern.row_length(sparsity_pattern.n_rows()-1) *
	       chunk_size;
	  return n;
	}

      else
	{
					   // if columns don't align, then
					   // just iterate over all chunks and
					   // see what this leads to
	  SparsityPattern::const_iterator p = sparsity_pattern.begin();
	  unsigned int n = 0;
	  for (; p!=sparsity_pattern.end(); ++p)
	    if ((p->row() != sparsity_pattern.n_rows() - 1)
		&&
		(p->column() != sparsity_pattern.n_cols() - 1))
	      n += chunk_size * chunk_size;
	    else
	      if ((p->row() == sparsity_pattern.n_rows() - 1)
		  &&
		  (p->column() != sparsity_pattern.n_cols() - 1))
						 // last chunk row, but not
						 // last chunk column. only a
						 // smaller number (n_rows %
						 // chunk_size) of rows
						 // actually exist
		n += (n_rows() % chunk_size) * chunk_size;
	      else
		if ((p->row() != sparsity_pattern.n_rows() - 1)
		    &&
		    (p->column() == sparsity_pattern.n_cols() - 1))
						   // last chunk column, but
						   // not row
		  n += (n_cols() % chunk_size) * chunk_size;
		else
						   // bottom right chunk
		  n += (n_cols() % chunk_size) *
		       (n_rows() % chunk_size);

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

  for (unsigned int i=0; i<sparsity_pattern.rows; ++i)
    for (unsigned int d=0;
	 (d<chunk_size) && (i*chunk_size + d < n_rows());
	  ++d)
      {
	out << '[' << i*chunk_size+d;
	for (unsigned int j=sparsity_pattern.rowstart[i];
	     j<sparsity_pattern.rowstart[i+1]; ++j)
	  if (sparsity_pattern.colnums[j] != sparsity_pattern.invalid_entry)
	    for (unsigned int e=0;
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

				   // for each entry in the underlying
				   // sparsity pattern, repeat everything
				   // chunk_size x chunk_size times
  for (unsigned int i=0; i<sparsity_pattern.rows; ++i)
    for (unsigned int j=sparsity_pattern.rowstart[i];
	 j<sparsity_pattern.rowstart[i+1]; ++j)
      if (sparsity_pattern.colnums[j] != sparsity_pattern.invalid_entry)
	for (unsigned int d=0;
	     ((d<chunk_size) &&
	      (sparsity_pattern.colnums[j]*chunk_size+d < n_cols()));
	     ++d)
	  for (unsigned int e=0;
	       (e<chunk_size) && (i*chunk_size + e < n_rows());
	       ++e)
					     // while matrix entries are
					     // usually written (i,j), with i
					     // vertical and j horizontal,
					     // gnuplot output is x-y, that is
					     // we have to exchange the order
					     // of output
	    out << sparsity_pattern.colnums[j]*chunk_size+d << " "
		<< -static_cast<signed int>(i*chunk_size+e)
		<< std::endl;

  AssertThrow (out, ExcIO());
}



unsigned int
ChunkSparsityPattern::bandwidth () const
{
				   // calculate the bandwidth from that of the
				   // underlying sparsity pattern. note that
				   // even if the bandwidth of that is zero,
				   // then the bandwidth of the chunky pattern
				   // is chunk_size-1, if it is 1 then the
				   // chunky pattern has
				   // chunk_size+(chunk_size-1), etc
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

                                   // first the simple objects,
                                   // bracketed in [...]
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

				   // then read the underlying sparsity
				   // pattern
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
void ChunkSparsityPattern::copy_from<CompressedSparsityPattern> (const CompressedSparsityPattern &,
								 const unsigned int,
								 const bool);
template
void ChunkSparsityPattern::copy_from<CompressedSetSparsityPattern> (const CompressedSetSparsityPattern &,
								    const unsigned int,
								    const bool);
template
void ChunkSparsityPattern::copy_from<CompressedSimpleSparsityPattern> (const CompressedSimpleSparsityPattern &,
								       const unsigned int,
								       const bool);
template
void ChunkSparsityPattern::copy_from<float> (const FullMatrix<float> &,
					     const unsigned int,
					     const bool);
template
void ChunkSparsityPattern::copy_from<double> (const FullMatrix<double> &,
					      const unsigned int,
					      const bool);

DEAL_II_NAMESPACE_CLOSE
