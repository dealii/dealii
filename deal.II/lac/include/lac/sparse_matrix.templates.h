//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008, 2009, 2010 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------

#ifndef __deal2__sparse_matrix_templates_h
#define __deal2__sparse_matrix_templates_h


#include <base/config.h>
#include <base/template_constraints.h>
#include <base/parallel.h>
#include <base/thread_management.h>
#include <base/multithread_info.h>
#include <lac/sparse_matrix.h>
#include <lac/vector.h>
#include <lac/full_matrix.h>
#include <lac/compressed_simple_sparsity_pattern.h>


// we only need output streams, but older compilers did not provide
// them in a separate include file
#ifdef HAVE_STD_OSTREAM_HEADER
#  include <ostream>
#else
#  include <iostream>
#endif

#include <iomanip>
#include <algorithm>
#include <functional>
#include <cmath>
#include <vector>
#include <numeric>
#include <base/std_cxx1x/bind.h>



DEAL_II_NAMESPACE_OPEN


template <typename number>
SparseMatrix<number>::SparseMatrix ()
                :
		cols(0, "SparseMatrix"),
		val(0),
		max_len(0)
{}



template <typename number>
SparseMatrix<number>::SparseMatrix (const SparseMatrix &m)
                :
		Subscriptor (m),
		cols(0, "SparseMatrix"),
		val(0),
		max_len(0)
{
  Assert (m.cols==0, ExcInvalidConstructorCall());
  Assert (m.val==0, ExcInvalidConstructorCall());
  Assert (m.max_len==0, ExcInvalidConstructorCall());
}



template <typename number>
SparseMatrix<number>&
SparseMatrix<number>::operator = (const SparseMatrix<number> &m)
{
  Assert (m.cols==0, ExcInvalidConstructorCall());
  Assert (m.val==0, ExcInvalidConstructorCall());
  Assert (m.max_len==0, ExcInvalidConstructorCall());

  return *this;
}



template <typename number>
SparseMatrix<number>::SparseMatrix (const SparsityPattern &c)
                :
		cols(0, "SparseMatrix"),
		val(0),
		max_len(0)
{
  reinit (c);
}



template <typename number>
SparseMatrix<number>::SparseMatrix (const SparsityPattern &c,
				    const IdentityMatrix  &id)
                :
		cols(0, "SparseMatrix"),
		val(0),
		max_len(0)
{
  Assert (c.n_rows() == id.m(), ExcDimensionMismatch (c.n_rows(), id.m()));
  Assert (c.n_cols() == id.n(), ExcDimensionMismatch (c.n_cols(), id.n()));

  reinit (c);
  for (unsigned int i=0; i<n(); ++i)
    this->set(i,i,1.);
}



template <typename number>
SparseMatrix<number>::~SparseMatrix ()
{
  cols = 0;

  if (val != 0)
    delete[] val;
}



template <typename number>
SparseMatrix<number> &
SparseMatrix<number>::operator = (const double d)
{
  Assert (d==0, ExcScalarAssignmentOnlyForZeroValue());

  Assert (cols != 0, ExcNotInitialized());
  Assert (cols->compressed || cols->empty(), SparsityPattern::ExcNotCompressed());

  if (val != 0)
    memset (&val[0], 0, cols->n_nonzero_elements()*sizeof(number));

  return *this;
}



template <typename number>
SparseMatrix<number> &
SparseMatrix<number>::operator= (const IdentityMatrix  &id)
{
  Assert (cols->n_rows() == id.m(),
	  ExcDimensionMismatch (cols->n_rows(), id.m()));
  Assert (cols->n_cols() == id.n(),
	  ExcDimensionMismatch (cols->n_cols(), id.n()));

  *this = 0;
  for (unsigned int i=0; i<n(); ++i)
    this->set(i,i,1.);

  return *this;
}



template <typename number>
void
SparseMatrix<number>::reinit (const SparsityPattern &sparsity)
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

  const std::size_t N = cols->n_nonzero_elements();
  if (N > max_len || max_len == 0)
    {
      if (val != 0)
        delete[] val;
      val = new number[N];
      max_len = N;
    }

  *this = 0.;
}



template <typename number>
void
SparseMatrix<number>::clear ()
{
  cols = 0;
  if (val) delete[] val;
  val = 0;
  max_len = 0;
}



template <typename number>
bool
SparseMatrix<number>::empty () const
{
  if (cols == 0)
    return true;
  else
    return cols->empty();
}



template <typename number>
unsigned int
SparseMatrix<number>::n_nonzero_elements () const
{
  Assert (cols != 0, ExcNotInitialized());
  return cols->n_nonzero_elements ();
}



template <typename number>
unsigned int
SparseMatrix<number>::n_actually_nonzero_elements (const double threshold) const
{
  Assert (cols != 0, ExcNotInitialized());
  Assert (threshold >= 0, ExcMessage ("Negative threshold!"));
  unsigned int nnz = 0;
  for (unsigned int i=0; i<n_nonzero_elements(); ++i)
    if (std::fabs(val[i]) > threshold)
      ++nnz;
  return nnz;
}



template <typename number>
void
SparseMatrix<number>::symmetrize ()
{
  Assert (cols != 0, ExcNotInitialized());
  Assert (cols->rows == cols->cols, ExcNotQuadratic());

  const unsigned int n_rows = m();
  for (unsigned int row=0; row<n_rows; ++row)
    {
				       // first skip diagonal entry
      number             *val_ptr = &val[cols->rowstart[row]];
      if (cols->optimize_diagonal())
	  ++val_ptr;
      const unsigned int *colnum_ptr = &cols->colnums[cols->rowstart[row]+1];
      const number       *const val_end_of_row = &val[cols->rowstart[row+1]];

				       // treat lower left triangle
      while ((val_ptr != val_end_of_row) && (*colnum_ptr<row))
	{
					   // compute the mean of this
					   // and the transpose value
	  const number mean_value = (*val_ptr +
				     val[(*cols)(*colnum_ptr,row)]) / 2.0;
					   // set this value and the
					   // transpose one to the
					   // mean
	  *val_ptr = mean_value;
	  set (*colnum_ptr, row, mean_value);

					   // advance pointers
	  ++val_ptr;
	  ++colnum_ptr;
	};
    };
}



template <typename number>
template <typename somenumber>
SparseMatrix<number> &
SparseMatrix<number>::copy_from (const SparseMatrix<somenumber> &matrix)
{
  Assert (cols != 0, ExcNotInitialized());
  Assert (val != 0, ExcNotInitialized());
  Assert (cols == matrix.cols, ExcDifferentSparsityPatterns());

  std::copy (&matrix.val[0], &matrix.val[cols->n_nonzero_elements()],
	     &val[0]);

  return *this;
}



template <typename number>
template <typename somenumber>
void
SparseMatrix<number>::copy_from (const FullMatrix<somenumber> &matrix)
{
				   // first delete previous content
  *this = 0;

				   // then copy old matrix
  for (unsigned int row=0; row<matrix.m(); ++row)
    for (unsigned int col=0; col<matrix.n(); ++col)
      if (matrix(row,col) != 0)
	set (row, col, matrix(row,col));
}



template <typename number>
template <typename somenumber>
void
SparseMatrix<number>::add (const number factor,
			   const SparseMatrix<somenumber> &matrix)
{
  Assert (cols != 0, ExcNotInitialized());
  Assert (val != 0, ExcNotInitialized());
  Assert (cols == matrix.cols, ExcDifferentSparsityPatterns());

  number             *val_ptr    = &val[0];
  const somenumber   *matrix_ptr = &matrix.val[0];
  const number *const end_ptr    = &val[cols->n_nonzero_elements()];

  while (val_ptr != end_ptr)
    *val_ptr++ += factor * *matrix_ptr++;
}


namespace internal
{
  namespace SparseMatrix
  {
				     /**
				      * Perform a vmult using the SparseMatrix
				      * data structures, but only using a
				      * subinterval for the row indices.
				      *
				      * In the sequential case, this function
				      * is called on all rows, in the parallel
				      * case it may be called on a subrange,
				      * at the discretion of the task
				      * scheduler.
				      */
    template <typename number,
	      typename InVector,
	      typename OutVector>
    void vmult_on_subrange (const unsigned int  begin_row,
			    const unsigned int  end_row,
			    const number       *values,
			    const std::size_t  *rowstart,
			    const unsigned int *colnums,
			    const InVector     &src,
			    OutVector          &dst,
			    const bool          add)
    {
      const number       *val_ptr    = &values[rowstart[begin_row]];
      const unsigned int *colnum_ptr = &colnums[rowstart[begin_row]];
      typename OutVector::iterator dst_ptr = dst.begin() + begin_row;

      if (add == false)
	for (unsigned int row=begin_row; row<end_row; ++row)
	  {
	    typename OutVector::value_type s = 0.;
	    const number *const val_end_of_row = &values[rowstart[row+1]];
	    while (val_ptr != val_end_of_row)
	      s += *val_ptr++ * src(*colnum_ptr++);
	    *dst_ptr++ = s;
	  }
      else
	for (unsigned int row=begin_row; row<end_row; ++row)
	  {
	    typename OutVector::value_type s = *dst_ptr;
	    const number *const val_end_of_row = &values[rowstart[row+1]];
	    while (val_ptr != val_end_of_row)
	      s += *val_ptr++ * src(*colnum_ptr++);
	    *dst_ptr++ = s;
	  }
    }
  }
}


template <typename number>
template <typename number2>
void
SparseMatrix<number>::add (const unsigned int  row,
			   const unsigned int  n_cols,
			   const unsigned int *col_indices,
			   const number2      *values,
			   const bool          elide_zero_values,
			   const bool          col_indices_are_sorted)
{
  Assert (cols != 0, ExcNotInitialized());

				   // if we have sufficiently many columns
				   // and sorted indices it is faster to
				   // just go through the column indices and
				   // look whether we found one, rather than
				   // doing many binary searches
  if (elide_zero_values == false && col_indices_are_sorted == true &&
      n_cols > 3)
    {
				   // check whether the given indices are
				   // really sorted
#ifdef DEBUG
      for (unsigned int i=1; i<n_cols; ++i)
	Assert (col_indices[i] > col_indices[i-1],
		ExcMessage("List of indices not sorted or with duplicates."));
#endif

      const unsigned int * this_cols =
	&cols->get_column_numbers()[cols->get_rowstart_indices()[row]];
      number * val_ptr = &val[cols->get_rowstart_indices()[row]];

      if (cols->optimize_diagonal() == true)
	{

				   // find diagonal and add it if found
	  Assert (this_cols[0] == row, ExcInternalError());
	  const unsigned int * diag_pos =
	    internals::SparsityPatternTools::optimized_lower_bound (col_indices,
								    col_indices+n_cols,
								    row);
	  const unsigned int diag = diag_pos - col_indices;
	  unsigned int post_diag = diag;
	  if (diag != n_cols && *diag_pos == row)
	    {
	      val_ptr[0] += *(values+(diag_pos-col_indices));
	      ++post_diag;
	    }

				   // add indices before diagonal
	  unsigned int counter = 1;
	  for (unsigned int i=0; i<diag; ++i)
	    {
	      Assert (col_indices[i] >= this_cols[counter], ExcInternalError());

	      while (this_cols[counter] < col_indices[i])
		++counter;

	      Assert (this_cols[counter] == col_indices[i] || values[i] == 0,
		      ExcInvalidIndex(row,col_indices[i]));

	      val_ptr[counter] += values[i];
	    }

				   // add indices after diagonal
	  for (unsigned int i=post_diag; i<n_cols; ++i)
	    {
	      Assert (col_indices[i] >= this_cols[counter], ExcInternalError());

	      while (this_cols[counter] < col_indices[i])
		++counter;

	      Assert (this_cols[counter] == col_indices[i] || values[i] == 0,
		      ExcInvalidIndex(row,col_indices[i]));

	      val_ptr[counter] += values[i];
	    }

	  Assert (counter < cols->row_length(row), ExcInternalError());
	}
      else
	{
	  unsigned int counter = 0;
	  for (unsigned int i=0; i<n_cols; ++i)
	    {
	      Assert (col_indices[i] >= this_cols[counter], ExcInternalError());

	      while (this_cols[counter] < col_indices[i])
		++counter;

	      Assert (this_cols[counter] == col_indices[i] || values[i] == 0,
		      ExcInvalidIndex(row,col_indices[i]));

	      val_ptr[counter] += values[i];
	    }
	  Assert (counter < cols->row_length(row), ExcInternalError());
	}
      return;
    }

				   // unsorted case: first, search all the
				   // indices to find out which values we
				   // actually need to add.
  const unsigned int * const my_cols = cols->get_column_numbers();
  unsigned int index = cols->get_rowstart_indices()[row];
  const unsigned int next_row_index = cols->get_rowstart_indices()[row+1];

  for (unsigned int j=0; j<n_cols; ++j)
    {
      const number value = values[j];
      Assert (numbers::is_finite(value),
	      ExcMessage("The given value is not finite but either "
			 "infinite or Not A Number (NaN)"));

#ifdef DEBUG
      if (elide_zero_values==true && value == 0)
	continue;
#else
      if (value == 0)
	continue;
#endif

				   // check whether the next index to add is
				   // the next present index in the sparsity
				   // pattern (otherwise, do a binary
				   // search)
      if (index < next_row_index && my_cols[index] == col_indices[j])
	goto add_value;

      index = cols->operator()(row, col_indices[j]);

				   // it is allowed to add elements to
				   // the matrix that are not part of
				   // the sparsity pattern, if the
				   // value we add is zero
      if (index == SparsityPattern::invalid_entry)
	{
	  Assert (value == 0., ExcInvalidIndex(row,col_indices[j]));
	  continue;
	}

    add_value:
      val[index] += value;
      ++index;
    }
}



template <typename number>
template <typename number2>
void
SparseMatrix<number>::set (const unsigned int  row,
			   const unsigned int  n_cols,
			   const unsigned int *col_indices,
			   const number2      *values,
			   const bool          elide_zero_values)
{
  Assert (cols != 0, ExcNotInitialized());
  Assert (row < m(), ExcInvalidIndex1(row));

				   // First, search all the indices to find
				   // out which values we actually need to
				   // set.
  const unsigned int * my_cols = cols->get_column_numbers();
  std::size_t index = cols->get_rowstart_indices()[row], next_index = index;
  const std::size_t next_row_index = cols->get_rowstart_indices()[row+1];

  if (elide_zero_values == true)
    {
      for (unsigned int j=0; j<n_cols; ++j)
	{
	  const number value = values[j];
	  Assert (numbers::is_finite(value),
		  ExcMessage("The given value is not finite but either "
			     "infinite or Not A Number (NaN)"));

	  if (value == 0)
	    continue;

				   // check whether the next index to set is
				   // the next present index in the sparsity
				   // pattern (otherwise, do a binary
				   // search)
	  if (index != next_row_index && my_cols[index] == col_indices[j])
	    goto set_value;

	  next_index = cols->operator()(row, col_indices[j]);

				   // it is allowed to set elements in
				   // the matrix that are not part of
				   // the sparsity pattern, if the
				   // value to which we set it is zero
	  if (next_index == SparsityPattern::invalid_entry)
	    {
	      Assert (false, ExcInvalidIndex(row,col_indices[j]));
	      continue;
	    }
	  index = next_index;

	set_value:
	  val[index] = value;
	  ++index;
	}
    }
  else
    {
				// same code as above, but now check for zeros
      for (unsigned int j=0; j<n_cols; ++j)
	{
	  const number value = values[j];
	  Assert (numbers::is_finite(value),
		  ExcMessage("The given value is not finite but either "
			     "infinite or Not A Number (NaN)"));

	  if (index != next_row_index && my_cols[index] == col_indices[j])
	    goto set_value_checked;

	  next_index = cols->operator()(row, col_indices[j]);

	  if (next_index == SparsityPattern::invalid_entry)
	    {
	      Assert (value == 0., ExcInvalidIndex(row,col_indices[j]));
	      continue;
	    }
	  index = next_index;

	set_value_checked:
	  val[index] = value;
	  ++index;
	}
    }
}



template <typename number>
template <class OutVector, class InVector>
void
SparseMatrix<number>::vmult (OutVector& dst,
			     const InVector& src) const
{
  Assert (cols != 0, ExcNotInitialized());
  Assert (val != 0, ExcNotInitialized());
  Assert(m() == dst.size(), ExcDimensionMismatch(m(),dst.size()));
  Assert(n() == src.size(), ExcDimensionMismatch(n(),src.size()));

  Assert (!PointerComparison::equal(&src, &dst), ExcSourceEqualsDestination());

  parallel::apply_to_subranges (0U, m(),
				std_cxx1x::bind (internal::SparseMatrix::vmult_on_subrange
						 <number,InVector,OutVector>,
						 _1, _2,
						 val,
						 cols->rowstart,
						 cols->colnums,
						 std_cxx1x::cref(src),
						 std_cxx1x::ref(dst),
						 false),
				internal::SparseMatrix::minimum_parallel_grain_size);
}



template <typename number>
template <class OutVector, class InVector>
void
SparseMatrix<number>::Tvmult (OutVector& dst,
                              const InVector& src) const
{
  Assert (val != 0, ExcNotInitialized());
  Assert (cols != 0, ExcNotInitialized());
  Assert(n() == dst.size(), ExcDimensionMismatch(n(),dst.size()));
  Assert(m() == src.size(), ExcDimensionMismatch(m(),src.size()));

  Assert (!PointerComparison::equal(&src, &dst), ExcSourceEqualsDestination());

  dst = 0;

  for (unsigned int i=0;i<m();i++)
    {
      for (unsigned int j=cols->rowstart[i]; j<cols->rowstart[i+1] ;j++)
	{
	  const unsigned int p = cols->colnums[j];
	  dst(p) += val[j] * src(i);
	}
    }
}



template <typename number>
template <class OutVector, class InVector>
void
SparseMatrix<number>::vmult_add (OutVector& dst,
                                 const InVector& src) const
{
  Assert (cols != 0, ExcNotInitialized());
  Assert (val != 0, ExcNotInitialized());
  Assert(m() == dst.size(), ExcDimensionMismatch(m(),dst.size()));
  Assert(n() == src.size(), ExcDimensionMismatch(n(),src.size()));

  Assert (!PointerComparison::equal(&src, &dst), ExcSourceEqualsDestination());

  parallel::apply_to_subranges (0U, m(),
				std_cxx1x::bind (internal::SparseMatrix::vmult_on_subrange
						 <number,InVector,OutVector>,
						 _1, _2,
						 val,
						 cols->rowstart,
						 cols->colnums,
						 std_cxx1x::cref(src),
						 std_cxx1x::ref(dst),
						 true),
				internal::SparseMatrix::minimum_parallel_grain_size);
}



template <typename number>
template <class OutVector, class InVector>
void
SparseMatrix<number>::Tvmult_add (OutVector& dst,
                                  const InVector& src) const
{
  Assert (val != 0, ExcNotInitialized());
  Assert (cols != 0, ExcNotInitialized());
  Assert(n() == dst.size(), ExcDimensionMismatch(n(),dst.size()));
  Assert(m() == src.size(), ExcDimensionMismatch(m(),src.size()));

  Assert (!PointerComparison::equal(&src, &dst), ExcSourceEqualsDestination());

  for (unsigned int i=0;i<m();i++)
    for (unsigned int j=cols->rowstart[i]; j<cols->rowstart[i+1] ;j++)
      {
        const unsigned int p = cols->colnums[j];
        dst(p) += val[j] * src(i);
      }
}


namespace internal
{
  namespace SparseMatrix
  {
				     /**
				      * Perform a vmult using the SparseMatrix
				      * data structures, but only using a
				      * subinterval for the row indices.
				      *
				      * In the sequential case, this function
				      * is called on all rows, in the parallel
				      * case it may be called on a subrange,
				      * at the discretion of the task
				      * scheduler.
				      */
    template <typename number,
	      typename InVector>
    number matrix_norm_sqr_on_subrange (const unsigned int  begin_row,
					const unsigned int  end_row,
					const number       *values,
					const std::size_t  *rowstart,
					const unsigned int *colnums,
					const InVector     &v)
    {
      number norm_sqr=0.;

      for (unsigned int i=begin_row; i<end_row; ++i)
	{
	  number s = 0;
	  for (unsigned int j=rowstart[i]; j<rowstart[i+1] ;j++)
	    s += values[j] * v(colnums[j]);
	  norm_sqr += v(i)*numbers::NumberTraits<number>::conjugate(s);
	}
      return norm_sqr;
    }
  }
}



template <typename number>
template <typename somenumber>
somenumber
SparseMatrix<number>::matrix_norm_square (const Vector<somenumber>& v) const
{
  Assert (cols != 0, ExcNotInitialized());
  Assert (val != 0, ExcNotInitialized());
  Assert(m() == v.size(), ExcDimensionMismatch(m(),v.size()));
  Assert(n() == v.size(), ExcDimensionMismatch(n(),v.size()));

  return
    parallel::accumulate_from_subranges<number>
    (std_cxx1x::bind (internal::SparseMatrix::matrix_norm_sqr_on_subrange
		      <number,Vector<somenumber> >,
		      _1, _2,
		      val, cols->rowstart, cols->colnums,
		      std_cxx1x::cref(v)),
     0, m(),
     internal::SparseMatrix::minimum_parallel_grain_size);
}



namespace internal
{
  namespace SparseMatrix
  {
				     /**
				      * Perform a vmult using the SparseMatrix
				      * data structures, but only using a
				      * subinterval for the row indices.
				      *
				      * In the sequential case, this function
				      * is called on all rows, in the parallel
				      * case it may be called on a subrange,
				      * at the discretion of the task
				      * scheduler.
				      */
    template <typename number,
	      typename InVector>
    number matrix_scalar_product_on_subrange (const unsigned int  begin_row,
					      const unsigned int  end_row,
					      const number       *values,
					      const std::size_t  *rowstart,
					      const unsigned int *colnums,
					      const InVector     &u,
					      const InVector     &v)
    {
      number norm_sqr=0.;

      for (unsigned int i=begin_row; i<end_row; ++i)
	{
	  number s = 0;
	  for (unsigned int j=rowstart[i]; j<rowstart[i+1] ;j++)
	    s += values[j] * v(colnums[j]);
	  norm_sqr += u(i)*numbers::NumberTraits<number>::conjugate(s);
	}
      return norm_sqr;
    }
  }
}



template <typename number>
template <typename somenumber>
somenumber
SparseMatrix<number>::matrix_scalar_product (const Vector<somenumber>& u,
					     const Vector<somenumber>& v) const
{
  Assert (cols != 0, ExcNotInitialized());
  Assert (val != 0, ExcNotInitialized());
  Assert(m() == u.size(), ExcDimensionMismatch(m(),u.size()));
  Assert(n() == v.size(), ExcDimensionMismatch(n(),v.size()));

  return
    parallel::accumulate_from_subranges<number>
    (std_cxx1x::bind (internal::SparseMatrix::matrix_scalar_product_on_subrange
		      <number,Vector<somenumber> >,
		      _1, _2,
		      val, cols->rowstart, cols->colnums,
		      std_cxx1x::cref(u),
		      std_cxx1x::cref(v)),
     0, m(),
     internal::SparseMatrix::minimum_parallel_grain_size);
}



template <typename number>
template <typename numberB, typename numberC, typename numberV>
void
SparseMatrix<number>::mmult (SparseMatrix<numberC>       &C,
			     const SparseMatrix<numberB> &B,
			     const Vector<numberV>       &V,
			     const bool                   rebuild_sparsity_C) const
{
  const bool use_vector = V.size() == n() ? true : false;
  Assert (n() == B.m(), ExcDimensionMismatch(n(), B.m()));
  Assert (cols != 0, ExcNotInitialized());
  Assert (B.cols != 0, ExcNotInitialized());
  Assert (C.cols != 0, ExcNotInitialized());

  const SparsityPattern & sp_A = *cols;
  const SparsityPattern & sp_B = *B.cols;

				   // clear previous content of C
  if  (rebuild_sparsity_C == true)
    {
				   // need to change the sparsity pattern of
				   // C, so cast away const-ness.
      SparsityPattern & sp_C =
	*(const_cast<SparsityPattern *>(&C.get_sparsity_pattern()));
      C.clear();
      sp_C.reinit (0,0,0);

				   // create a sparsity pattern for the
				   // matrix. we will go through all the
				   // rows in the matrix A, and for each
				   // column in a row we add the whole row
				   // of matrix B with that row number. This
				   // means that we will insert a lot of
				   // entries to each row, which is best
				   // handled by the
				   // CompressedSimpleSparsityPattern class.
      {
	CompressedSimpleSparsityPattern csp (m(), B.n());
	for (unsigned int i = 0; i < csp.n_rows(); ++i)
	  {
	    const unsigned int * rows =
	      &sp_A.get_column_numbers()[sp_A.get_rowstart_indices()[i]];
	    const unsigned int *const end_rows =
	      &sp_A.get_column_numbers()[sp_A.get_rowstart_indices()[i+1]];
	    for (; rows != end_rows; ++rows)
	      {
		const unsigned int col = *rows;
		unsigned int * new_cols = const_cast<unsigned int*>
		  (&sp_B.get_column_numbers()
		   [sp_B.get_rowstart_indices()[col]]);
		unsigned int * end_new_cols = const_cast<unsigned int*>
		  (&sp_B.get_column_numbers()
		   [sp_B.get_rowstart_indices()[col+1]]);

				   // if B has a diagonal, need to add that
				   // manually. this way, we maintain
				   // sortedness.
		if (sp_B.optimize_diagonal() == true)
		  {
		    ++new_cols;
		    csp.add(i, col);
		  }

		csp.add_entries (i, new_cols, end_new_cols, true);
	      }
	  }
	sp_C.copy_from (csp);
      }

				   // reinit matrix C from that information
      C.reinit (sp_C);
    }

  Assert (C.m() == m(), ExcDimensionMismatch(C.m(), m()));
  Assert (C.n() == B.n(), ExcDimensionMismatch(C.n(), B.n()));

				   // create an array that caches some
				   // elements that are going to be written
				   // into the new matrix.
  unsigned int max_n_cols_B = 0;
  for (unsigned int i=0; i<B.m(); ++i)
    max_n_cols_B = std::max (max_n_cols_B, sp_B.row_length(i));
  std::vector<numberC> new_entries(max_n_cols_B);

				   // now compute the actual entries: a
				   // matrix-matrix product involves three
				   // nested loops. One over the rows of A,
				   // for each row we then loop over all the
				   // columns, and then we need to multiply
				   // each element with all the elements in
				   // that row in B.
  for (unsigned int i=0; i<C.m(); ++i)
    {
      const unsigned int * rows =
	&sp_A.get_column_numbers()[sp_A.get_rowstart_indices()[i]];
      const unsigned int *const end_rows =
	&sp_A.get_column_numbers()[sp_A.get_rowstart_indices()[i+1]];
      for (; rows != end_rows; ++rows)
	{
	  const double A_val = global_entry
	    (rows-&sp_A.get_column_numbers()[sp_A.get_rowstart_indices()[0]]);
	  const unsigned int col = *rows;
	  const unsigned int * new_cols =
	    (&sp_B.get_column_numbers()[sp_B.get_rowstart_indices()[col]]);

				   // special treatment for diagonal
	  if (sp_B.optimize_diagonal())
	    {
	      C.add (i, *new_cols, A_val *
		     B.global_entry(new_cols-&sp_B.get_column_numbers()
				    [sp_B.get_rowstart_indices()[0]]) *
		     (use_vector ? V(col) : 1));
	      ++new_cols;
	    }

				   // now the innermost loop that goes over
				   // all the elements in row 'col' of
				   // matrix B. Cache the elements, and then
				   // write them into C at once
	  numberC * new_ptr = &new_entries[0];
	  const numberB * B_val_ptr =
	    &B.val[new_cols-&sp_B.get_column_numbers()[sp_B.get_rowstart_indices()[0]]];
	  const numberB * const end_cols =
	    &B.val[&sp_B.get_column_numbers()[sp_B.get_rowstart_indices()[col+1]]-
		   &sp_B.get_column_numbers()[sp_B.get_rowstart_indices()[0]]];
	  for (; B_val_ptr != end_cols; ++B_val_ptr)
	    *new_ptr++ = A_val * *B_val_ptr * (use_vector ? V(col) : 1);

	  C.add (i, new_ptr-&new_entries[0], new_cols, &new_entries[0],
		 false, true);
	}
    }
}




template <typename number>
template <typename numberB, typename numberC, typename numberV>
void
SparseMatrix<number>::Tmmult (SparseMatrix<numberC>       &C,
			      const SparseMatrix<numberB> &B,
			      const Vector<numberV>       &V,
			      const bool                   rebuild_sparsity_C) const
{
  const bool use_vector = V.size() == m() ? true : false;
  Assert (m() == B.m(), ExcDimensionMismatch(m(), B.m()));
  Assert (cols != 0, ExcNotInitialized());
  Assert (B.cols != 0, ExcNotInitialized());
  Assert (C.cols != 0, ExcNotInitialized());

  const SparsityPattern & sp_A = *cols;
  const SparsityPattern & sp_B = *B.cols;

				   // clear previous content of C
  if  (rebuild_sparsity_C == true)
    {
				   // need to change the sparsity pattern of
				   // C, so cast away const-ness.
      SparsityPattern & sp_C =
	*(const_cast<SparsityPattern *>(&C.get_sparsity_pattern()));
      C.clear();
      sp_C.reinit (0,0,0);

				   // create a sparsity pattern for the
				   // matrix. we will go through all the
				   // rows in the matrix A, and for each
				   // column in a row we add the whole row
				   // of matrix B with that row number. This
				   // means that we will insert a lot of
				   // entries to each row, which is best
				   // handled by the
				   // CompressedSimpleSparsityPattern class.
      {
	CompressedSimpleSparsityPattern csp (n(), B.n());
	for (unsigned int i = 0; i < sp_A.n_rows(); ++i)
	  {
	    const unsigned int * rows =
	      &sp_A.get_column_numbers()[sp_A.get_rowstart_indices()[i]];
	    const unsigned int *const end_rows =
	      &sp_A.get_column_numbers()[sp_A.get_rowstart_indices()[i+1]];
	    unsigned int * new_cols = const_cast<unsigned int*>
	      (&sp_B.get_column_numbers()
	       [sp_B.get_rowstart_indices()[i]]);
	    unsigned int * end_new_cols = const_cast<unsigned int*>
	      (&sp_B.get_column_numbers()
	       [sp_B.get_rowstart_indices()[i+1]]);

	    if (sp_B.optimize_diagonal() == true)
	      ++new_cols;

	    for (; rows != end_rows; ++rows)
	      {
		const unsigned int row = *rows;

				   // if B has a diagonal, need to add that
				   // manually. this way, we maintain
				   // sortedness.
		if (sp_B.optimize_diagonal() == true)
		  csp.add(row, i);

		csp.add_entries (row, new_cols, end_new_cols, true);
	      }
	  }
	sp_C.copy_from (csp);
      }

				   // reinit matrix C from that information
      C.reinit (sp_C);
    }

  Assert (C.m() == n(), ExcDimensionMismatch(C.m(), n()));
  Assert (C.n() == B.n(), ExcDimensionMismatch(C.n(), B.n()));

				   // create an array that caches some
				   // elements that are going to be written
				   // into the new matrix.
  unsigned int max_n_cols_B = 0;
  for (unsigned int i=0; i<B.m(); ++i)
    max_n_cols_B = std::max (max_n_cols_B, sp_B.row_length(i));
  std::vector<numberC> new_entries(max_n_cols_B);

				   // now compute the actual entries: a
				   // matrix-matrix product involves three
				   // nested loops. One over the rows of A,
				   // for each row we then loop over all the
				   // columns, and then we need to multiply
				   // each element with all the elements in
				   // that row in B.
  for (unsigned int i=0; i<m(); ++i)
    {
      const unsigned int * rows =
	&sp_A.get_column_numbers()[sp_A.get_rowstart_indices()[i]];
      const unsigned int *const end_rows =
	&sp_A.get_column_numbers()[sp_A.get_rowstart_indices()[i+1]];
      const unsigned int * new_cols =
	(&sp_B.get_column_numbers()[sp_B.get_rowstart_indices()[i]]);
      if (sp_B.optimize_diagonal())
	++new_cols;

      const numberB * const end_cols =
	&B.val[&sp_B.get_column_numbers()[sp_B.get_rowstart_indices()[i+1]]-
	       &sp_B.get_column_numbers()[sp_B.get_rowstart_indices()[0]]];

      for (; rows != end_rows; ++rows)
	{
	  const unsigned int row = *rows;
	  const double A_val = global_entry
	    (rows-&sp_A.get_column_numbers()[sp_A.get_rowstart_indices()[0]]);

				   // special treatment for diagonal
	  if (sp_B.optimize_diagonal())
	    C.add (row, i, A_val *
		   B.global_entry(new_cols-1-&sp_B.get_column_numbers()
				  [sp_B.get_rowstart_indices()[0]]) *
		   (use_vector ? V(i) : 1));

				   // now the innermost loop that goes over
				   // all the elements in row 'col' of
				   // matrix B. Cache the elements, and then
				   // write them into C at once
	  numberC * new_ptr = &new_entries[0];
	  const numberB * B_val_ptr =
	    &B.val[new_cols-&sp_B.get_column_numbers()[sp_B.get_rowstart_indices()[0]]];
	  for (; B_val_ptr != end_cols; ++B_val_ptr)
	    *new_ptr++ = A_val * *B_val_ptr * (use_vector ? V(i) : 1);

	  C.add (row, new_ptr-&new_entries[0], new_cols, &new_entries[0],
		 false, true);
	}
    }
}



template <typename number>
typename SparseMatrix<number>::real_type
SparseMatrix<number>::l1_norm () const
{
  Assert (cols != 0, ExcNotInitialized());
  Assert (val != 0, ExcNotInitialized());

  Vector<real_type> column_sums(n());
  const unsigned int n_rows = m();
  for (unsigned int row=0; row<n_rows; ++row)
    for (unsigned int j=cols->rowstart[row]; j<cols->rowstart[row+1] ; ++j)
      column_sums(cols->colnums[j]) += numbers::NumberTraits<number>::abs(val[j]);

  return column_sums.linfty_norm();
}



template <typename number>
typename SparseMatrix<number>::real_type
SparseMatrix<number>::linfty_norm () const
{
  Assert (cols != 0, ExcNotInitialized());
  Assert (val != 0, ExcNotInitialized());

  const number *val_ptr = &val[cols->rowstart[0]];

  real_type max=0;
  const unsigned int n_rows = m();
  for (unsigned int row=0; row<n_rows; ++row)
    {
      real_type sum = 0;
      const number *const val_end_of_row = &val[cols->rowstart[row+1]];
      while (val_ptr != val_end_of_row)
	sum += numbers::NumberTraits<number>::abs(*val_ptr++);
      if (sum > max)
	max = sum;
    }
  return max;
}



template <typename number>
typename SparseMatrix<number>::real_type
SparseMatrix<number>::frobenius_norm () const
{
                                   // simply add up all entries in the
                                   // sparsity pattern, without taking any
                                   // reference to rows or columns
  real_type norm_sqr = 0;
  const unsigned int n_rows = m();
  for (const number *ptr = &val[0];
       ptr != &val[cols->rowstart[n_rows]]; ++ptr)
    norm_sqr +=  numbers::NumberTraits<number>::abs_square(*ptr);

  return std::sqrt (norm_sqr);
}



namespace internal
{
  namespace SparseMatrix
  {
				     /**
				      * Perform a vmult using the SparseMatrix
				      * data structures, but only using a
				      * subinterval for the row indices.
				      *
				      * In the sequential case, this function
				      * is called on all rows, in the parallel
				      * case it may be called on a subrange,
				      * at the discretion of the task
				      * scheduler.
				      */
    template <typename number,
	      typename InVector,
	      typename OutVector>
    number residual_sqr_on_subrange (const unsigned int  begin_row,
				     const unsigned int  end_row,
				     const number       *values,
				     const std::size_t  *rowstart,
				     const unsigned int *colnums,
				     const InVector     &u,
				     const InVector     &b,
				     OutVector          &dst)
    {
      number norm_sqr=0.;

      for (unsigned int i=begin_row; i<end_row; ++i)
	{
	  number s = b(i);
	  for (unsigned int j=rowstart[i]; j<rowstart[i+1] ;j++)
	    s -= values[j] * u(colnums[j]);
	  dst(i) = s;
	  norm_sqr += s*numbers::NumberTraits<number>::conjugate(s);
	}
      return norm_sqr;
    }
  }
}


template <typename number>
template <typename somenumber>
somenumber
SparseMatrix<number>::residual (Vector<somenumber>       &dst,
				const Vector<somenumber> &u,
				const Vector<somenumber> &b) const
{
  Assert (cols != 0, ExcNotInitialized());
  Assert (val != 0, ExcNotInitialized());
  Assert(m() == dst.size(), ExcDimensionMismatch(m(),dst.size()));
  Assert(m() == b.size(), ExcDimensionMismatch(m(),b.size()));
  Assert(n() == u.size(), ExcDimensionMismatch(n(),u.size()));

  Assert (&u != &dst, ExcSourceEqualsDestination());

  return
    std::sqrt (parallel::accumulate_from_subranges<number>
	       (std_cxx1x::bind (internal::SparseMatrix::residual_sqr_on_subrange
				 <number,Vector<somenumber>,Vector<somenumber> >,
				 _1, _2,
				 val, cols->rowstart, cols->colnums,
				 std_cxx1x::cref(u),
				 std_cxx1x::cref(b),
				 std_cxx1x::ref(dst)),
		0, m(),
		internal::SparseMatrix::minimum_parallel_grain_size));
}



template <typename number>
template <typename somenumber>
void
SparseMatrix<number>::precondition_Jacobi (Vector<somenumber>       &dst,
					   const Vector<somenumber> &src,
					   const number              om) const
{
  Assert (cols != 0, ExcNotInitialized());
  Assert (val != 0, ExcNotInitialized());
  Assert (cols->optimize_diagonal(),
	  typename SparsityPattern::ExcDiagonalNotOptimized());

  Assert (dst.size() == n(), ExcDimensionMismatch (dst.size(), n()));
  Assert (src.size() == n(), ExcDimensionMismatch (src.size(), n()));

  const unsigned int n = src.size();
  somenumber              *dst_ptr = dst.begin();
  const somenumber        *src_ptr = src.begin();
  const std::size_t  *rowstart_ptr = &cols->rowstart[0];

				   // optimize the following loop for
				   // the case that the relaxation
				   // factor is one. In that case, we
				   // can save one FP multiplication
				   // per row
				   //
				   // note that for square matrices,
				   // the diagonal entry is the first
				   // in each row, i.e. at index
				   // rowstart[i]. and we do have a
				   // square matrix by above assertion
  if (om != 1.)
    for (unsigned int i=0; i<n; ++i, ++dst_ptr, ++src_ptr, ++rowstart_ptr)
      *dst_ptr = om * *src_ptr / val[*rowstart_ptr];
  else
    for (unsigned int i=0; i<n; ++i, ++dst_ptr, ++src_ptr, ++rowstart_ptr)
      *dst_ptr = *src_ptr / val[*rowstart_ptr];
}



template <typename number>
template <typename somenumber>
void
SparseMatrix<number>::precondition_SSOR (Vector<somenumber>              &dst,
					 const Vector<somenumber>        &src,
					 const number                     om,
					 const std::vector<unsigned int> &pos_right_of_diagonal) const
{
				   // to understand how this function works
				   // you may want to take a look at the CVS
				   // archives to see the original version
				   // which is much clearer...
  Assert (cols != 0, ExcNotInitialized());
  Assert (val != 0, ExcNotInitialized());
  Assert (cols->optimize_diagonal(),
	  typename SparsityPattern::ExcDiagonalNotOptimized());

  Assert (dst.size() == n(), ExcDimensionMismatch (dst.size(), n()));
  Assert (src.size() == n(), ExcDimensionMismatch (src.size(), n()));

  const unsigned int  n            = src.size();
  const std::size_t  *rowstart_ptr = &cols->rowstart[0];
  somenumber         *dst_ptr      = &dst(0);

				   // case when we have stored the position
				   // just right of the diagonal (then we
				   // don't have to search for it).
  if (pos_right_of_diagonal.size() != 0)
    {
      Assert (pos_right_of_diagonal.size() == dst.size(),
	      ExcDimensionMismatch (pos_right_of_diagonal.size(), dst.size()));

				   // forward sweep
      for (unsigned int row=0; row<n; ++row, ++dst_ptr, ++rowstart_ptr)
	{
	  *dst_ptr = src(row);
	  const unsigned int first_right_of_diagonal_index =
	    pos_right_of_diagonal[row];
	  Assert (first_right_of_diagonal_index <= *(rowstart_ptr+1),
		  ExcInternalError());
	  number s = 0;
	  for (unsigned int j=(*rowstart_ptr)+1; j<first_right_of_diagonal_index; ++j)
	    s += val[j] * dst(cols->colnums[j]);

				   // divide by diagonal element
	  *dst_ptr -= s * om;
	  *dst_ptr /= val[*rowstart_ptr];
	};

      rowstart_ptr = &cols->rowstart[0];
      dst_ptr      = &dst(0);
      for ( ; rowstart_ptr!=&cols->rowstart[n]; ++rowstart_ptr, ++dst_ptr)
	*dst_ptr *= (2.-om)*val[*rowstart_ptr];

				   // backward sweep
      rowstart_ptr = &cols->rowstart[n-1];
      dst_ptr      = &dst(n-1);
      for (int row=n-1; row>=0; --row, --rowstart_ptr, --dst_ptr)
	{
	  const unsigned int end_row = *(rowstart_ptr+1);
	  const unsigned int first_right_of_diagonal_index
	    = pos_right_of_diagonal[row];
	  number s = 0;
	  for (unsigned int j=first_right_of_diagonal_index; j<end_row; ++j)
	    s += val[j] * dst(cols->colnums[j]);

	  *dst_ptr -= s * om;
	  *dst_ptr /= val[*rowstart_ptr];
	};
      return;
    }

				   // case when we need to get the position
				   // of the first element right of the
				   // diagonal manually for each sweep.
				   // forward sweep
  for (unsigned int row=0; row<n; ++row, ++dst_ptr, ++rowstart_ptr)
    {
      *dst_ptr = src(row);
				       // find the first element in this line
				       // which is on the right of the diagonal.
				       // we need to precondition with the
				       // elements on the left only.
				       // note: the first entry in each
				       // line denotes the diagonal element,
				       // which we need not check.
      const unsigned int first_right_of_diagonal_index
	= (internals::SparsityPatternTools::optimized_lower_bound (&cols->colnums[*rowstart_ptr+1],
								   &cols->colnums[*(rowstart_ptr+1)],
								   row)
	   -
	   &cols->colnums[0]);

      number s = 0;
      for (unsigned int j=(*rowstart_ptr)+1; j<first_right_of_diagonal_index; ++j)
	s += val[j] * dst(cols->colnums[j]);

				       // divide by diagonal element
      *dst_ptr -= s * om;
      *dst_ptr /= val[*rowstart_ptr];
    };

  rowstart_ptr = &cols->rowstart[0];
  dst_ptr      = &dst(0);
  for (unsigned int row=0; row<n; ++row, ++rowstart_ptr, ++dst_ptr)
    *dst_ptr *= (2.-om)*val[*rowstart_ptr];

				   // backward sweep
  rowstart_ptr = &cols->rowstart[n-1];
  dst_ptr      = &dst(n-1);
  for (int row=n-1; row>=0; --row, --rowstart_ptr, --dst_ptr)
    {
      const unsigned int end_row = *(rowstart_ptr+1);
      const unsigned int first_right_of_diagonal_index
	= (internals::SparsityPatternTools::optimized_lower_bound (&cols->colnums[*rowstart_ptr+1],
								   &cols->colnums[end_row],
								   static_cast<unsigned int>(row)) -
	   &cols->colnums[0]);
      number s = 0;
      for (unsigned int j=first_right_of_diagonal_index; j<end_row; ++j)
	s += val[j] * dst(cols->colnums[j]);
      *dst_ptr -= s * om;
      *dst_ptr /= val[*rowstart_ptr];
    };
}


template <typename number>
template <typename somenumber>
void
SparseMatrix<number>::precondition_SOR (Vector<somenumber>& dst,
                                        const Vector<somenumber>& src,
                                        const number om) const
{
  Assert (cols != 0, ExcNotInitialized());
  Assert (val != 0, ExcNotInitialized());
  Assert (cols->optimize_diagonal(),
	  typename SparsityPattern::ExcDiagonalNotOptimized());


  dst = src;
  SOR(dst,om);
}


template <typename number>
template <typename somenumber>
void
SparseMatrix<number>::precondition_TSOR (Vector<somenumber>& dst,
                                         const Vector<somenumber>& src,
                                         const number om) const
{
  Assert (cols != 0, ExcNotInitialized());
  Assert (val != 0, ExcNotInitialized());
  Assert (cols->optimize_diagonal(),
	  typename SparsityPattern::ExcDiagonalNotOptimized());


  dst = src;
  TSOR(dst,om);
}


template <typename number>
template <typename somenumber>
void
SparseMatrix<number>::SOR (Vector<somenumber>& dst,
                           const number om) const
{
  Assert (cols != 0, ExcNotInitialized());
  Assert (val != 0, ExcNotInitialized());
  Assert (cols->optimize_diagonal(),
	  typename SparsityPattern::ExcDiagonalNotOptimized());

  Assert (m() == dst.size(), ExcDimensionMismatch(m(),dst.size()));

  for (unsigned int row=0; row<m(); ++row)
    {
      somenumber s = dst(row);
      for (unsigned int j=cols->rowstart[row]; j<cols->rowstart[row+1]; ++j)
	{
	  const unsigned int col = cols->colnums[j];
	  if (col < row)
	    s -= val[j] * dst(col);
	}

      dst(row) = s * om / val[cols->rowstart[row]];
    }
}


template <typename number>
template <typename somenumber>
void
SparseMatrix<number>::TSOR (Vector<somenumber>& dst,
                            const number om) const
{
  Assert (cols != 0, ExcNotInitialized());
  Assert (val != 0, ExcNotInitialized());
  Assert (cols->optimize_diagonal(),
	  typename SparsityPattern::ExcDiagonalNotOptimized());

  Assert (m() == dst.size(), ExcDimensionMismatch(m(),dst.size()));

  unsigned int row=m()-1;
  while (true)
    {
      somenumber s = dst(row);
      for (unsigned int j=cols->rowstart[row]; j<cols->rowstart[row+1]; ++j)
	if (cols->colnums[j] > row)
	  s -= val[j] * dst(cols->colnums[j]);

      dst(row) = s * om / val[cols->rowstart[row]];

      if (row == 0)
	break;

      --row;
    }
}


template <typename number>
template <typename somenumber>
void
SparseMatrix<number>::PSOR (Vector<somenumber>& dst,
                            const std::vector<unsigned int>& permutation,
                            const std::vector<unsigned int>& inverse_permutation,
                            const number om) const
{
  Assert (cols != 0, ExcNotInitialized());
  Assert (val != 0, ExcNotInitialized());
  Assert (cols->optimize_diagonal(),
	  typename SparsityPattern::ExcDiagonalNotOptimized());

  Assert (m() == dst.size(), ExcDimensionMismatch(m(),dst.size()));
  Assert (m() == permutation.size(),
	  ExcDimensionMismatch(m(), permutation.size()));
  Assert (m() == inverse_permutation.size(),
	  ExcDimensionMismatch(m(), inverse_permutation.size()));

  for (unsigned int urow=0; urow<m(); ++urow)
    {
      const unsigned int row = permutation[urow];
      somenumber s = dst(row);

      for (unsigned int j=cols->rowstart[row]; j<cols->rowstart[row+1]; ++j)
	{
	  const unsigned int col = cols->colnums[j];
	  if (inverse_permutation[col] < urow)
	    {
	      s -= val[j] * dst(col);
	    }
	}

      dst(row) = s * om / val[cols->rowstart[row]];
    }
}


template <typename number>
template <typename somenumber>
void
SparseMatrix<number>::TPSOR (Vector<somenumber>& dst,
                             const std::vector<unsigned int>& permutation,
                             const std::vector<unsigned int>& inverse_permutation,
                             const number om) const
{
  Assert (cols != 0, ExcNotInitialized());
  Assert (val != 0, ExcNotInitialized());
  Assert (cols->optimize_diagonal(),
	  typename SparsityPattern::ExcDiagonalNotOptimized());

  Assert (m() == dst.size(), ExcDimensionMismatch(m(),dst.size()));
  Assert (m() == permutation.size(),
	  ExcDimensionMismatch(m(), permutation.size()));
  Assert (m() == inverse_permutation.size(),
	  ExcDimensionMismatch(m(), inverse_permutation.size()));

  for (unsigned int urow=m(); urow != 0;)
    {
      --urow;
      const unsigned int row = permutation[urow];
      somenumber s = dst(row);
      for (unsigned int j=cols->rowstart[row]; j<cols->rowstart[row+1]; ++j)
	{
	  const unsigned int col = cols->colnums[j];
	  if (inverse_permutation[col] > urow)
	    s -= val[j] * dst(col);
	}

      dst(row) = s * om / val[cols->rowstart[row]];
    }
}



template <typename number>
template <typename somenumber>
void
SparseMatrix<number>::Jacobi_step (Vector<somenumber> &v,
				   const Vector<somenumber> &b,
				   const number        om) const
{
  Assert (cols != 0, ExcNotInitialized());
  Assert (val != 0, ExcNotInitialized());
  Assert (cols->optimize_diagonal(),
	  typename SparsityPattern::ExcDiagonalNotOptimized());

  Assert (m() == v.size(), ExcDimensionMismatch(m(),v.size()));
  Assert (m() == b.size(), ExcDimensionMismatch(m(),b.size()));

  GrowingVectorMemory<Vector<somenumber> > mem;
  typename VectorMemory<Vector<somenumber> >::Pointer w(mem);
  w->reinit(v);

  if (!v.all_zero())
    {
      vmult (*w, v);
      *w -= b;
    }
  else
    w->equ (-1.,b);
  precondition_Jacobi (*w, *w, om);
  v -= *w;
}



template <typename number>
template <typename somenumber>
void
SparseMatrix<number>::SOR_step (Vector<somenumber> &v,
                                const Vector<somenumber> &b,
                                const number        om) const
{
  Assert (cols != 0, ExcNotInitialized());
  Assert (val != 0, ExcNotInitialized());
  Assert (cols->optimize_diagonal(),
	  typename SparsityPattern::ExcDiagonalNotOptimized());

  Assert (m() == v.size(), ExcDimensionMismatch(m(),v.size()));
  Assert (m() == b.size(), ExcDimensionMismatch(m(),b.size()));

  for (unsigned int row=0; row<m(); ++row)
    {
      somenumber s = b(row);
      for (unsigned int j=cols->rowstart[row]; j<cols->rowstart[row+1]; ++j)
	{
	  s -= val[j] * v(cols->colnums[j]);
	}
      v(row) += s * om / val[cols->rowstart[row]];
    }
}



template <typename number>
template <typename somenumber>
void
SparseMatrix<number>::TSOR_step (Vector<somenumber> &v,
                                 const Vector<somenumber> &b,
                                 const number        om) const
{
  Assert (cols != 0, ExcNotInitialized());
  Assert (val != 0, ExcNotInitialized());
  Assert (cols->optimize_diagonal(),
	  typename SparsityPattern::ExcDiagonalNotOptimized());

  Assert (m() == v.size(), ExcDimensionMismatch(m(),v.size()));
  Assert (m() == b.size(), ExcDimensionMismatch(m(),b.size()));

  for (int row=m()-1; row>=0; --row)
    {
      somenumber s = b(row);
      for (unsigned int j=cols->rowstart[row]; j<cols->rowstart[row+1]; ++j)
	{
	  s -= val[j] * v(cols->colnums[j]);
	}
      v(row) += s * om / val[cols->rowstart[row]];
    }
}



template <typename number>
template <typename somenumber>
void
SparseMatrix<number>::SSOR_step (Vector<somenumber> &v,
                                 const Vector<somenumber> &b,
                                 const number        om) const
{
  SOR_step(v,b,om);
  TSOR_step(v,b,om);
}



template <typename number>
template <typename somenumber>
void
SparseMatrix<number>::SSOR (Vector<somenumber>& dst,
                            const number om) const
{
  Assert (cols != 0, ExcNotInitialized());
  Assert (val != 0, ExcNotInitialized());
  Assert (cols->optimize_diagonal(),
	  typename SparsityPattern::ExcDiagonalNotOptimized());

  Assert (m() == dst.size(), ExcDimensionMismatch(m(),dst.size()));

  const unsigned int  n = dst.size();
  unsigned int  j;
  somenumber s;

  for (unsigned int i=0; i<n; i++)
    {
      s = 0.;
      for (j=cols->rowstart[i]; j<cols->rowstart[i+1] ;j++)
	{
	  const unsigned int p = cols->colnums[j];
	  if (p != SparsityPattern::invalid_entry)
	    {
	      if (i>j) s += val[j] * dst(p);
	    }
	}
      dst(i) -= s * om;
      dst(i) /= val[cols->rowstart[i]];
    }

  for (int i=n-1; i>=0; i--)  // this time, i is signed, but always positive!
    {
      s = 0.;
      for (j=cols->rowstart[i]; j<cols->rowstart[i+1] ;j++)
	{
	  const unsigned int p = cols->colnums[j];
	  if (p != SparsityPattern::invalid_entry)
	    {
	      if (static_cast<unsigned int>(i)<j) s += val[j] * dst(p);
	    }
	}
      dst(i) -= s * om / val[cols->rowstart[i]];
    }
}



template <typename number>
const SparsityPattern &
SparseMatrix<number>::get_sparsity_pattern () const
{
  Assert (cols != 0, ExcNotInitialized());
  return *cols;
}



template <typename number>
void SparseMatrix<number>::print_formatted (std::ostream &out,
					    const unsigned int precision,
					    const bool scientific,
					    const unsigned int width_,
					    const char* zero_string,
					    const double denominator) const
{
  Assert (cols != 0, ExcNotInitialized());
  Assert (val != 0, ExcNotInitialized());

  unsigned int width = width_;

  std::ios::fmtflags old_flags = out.flags();
  unsigned int old_precision = out.precision (precision);

  if (scientific)
    {
      out.setf (std::ios::scientific, std::ios::floatfield);
      if (!width)
	width = precision+7;
    } else {
      out.setf (std::ios::fixed, std::ios::floatfield);
      if (!width)
	width = precision+2;
    }

  for (unsigned int i=0; i<m(); ++i)
    {
      for (unsigned int j=0; j<n(); ++j)
	if ((*cols)(i,j) != SparsityPattern::invalid_entry)
	  out << std::setw(width)
	      << val[cols->operator()(i,j)] * denominator << ' ';
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
void SparseMatrix<number>::print_pattern (std::ostream &out,
					  const double threshold) const
{
  Assert (cols != 0, ExcNotInitialized());
  Assert (val != 0, ExcNotInitialized());

  for (unsigned int i=0; i<m(); ++i)
    {
      for (unsigned int j=0; j<n(); ++j)
	if ((*cols)(i,j) == SparsityPattern::invalid_entry)
	  out << '.';
	else
	  if (std::fabs(val[cols->operator()(i,j)]) > threshold)
	    out << '*';
	  else
	    out << ':';
      out << std::endl;
    };
  AssertThrow (out, ExcIO());
}



template <typename number>
void
SparseMatrix<number>::block_write (std::ostream &out) const
{
  AssertThrow (out, ExcIO());

                                   // first the simple objects,
                                   // bracketed in [...]
  out << '[' << max_len << "][";
                                   // then write out real data
  out.write (reinterpret_cast<const char*>(&val[0]),
	     reinterpret_cast<const char*>(&val[max_len])
	     - reinterpret_cast<const char*>(&val[0]));
  out << ']';

  AssertThrow (out, ExcIO());
}



template <typename number>
void
SparseMatrix<number>::block_read (std::istream &in)
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
  in.read (reinterpret_cast<char*>(&val[0]),
           reinterpret_cast<char*>(&val[max_len])
           - reinterpret_cast<char*>(&val[0]));
  in >> c;
  AssertThrow (c == ']', ExcIO());
}



template <typename number>
std::size_t
SparseMatrix<number>::memory_consumption () const
{
  return max_len*static_cast<std::size_t>(sizeof(number)) + sizeof(*this);
}


DEAL_II_NAMESPACE_CLOSE

#endif
