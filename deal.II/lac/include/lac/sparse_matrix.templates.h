//----------------------------  sparse_matrix.templates.h  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  sparse_matrix.templates.h  ---------------------------
#ifndef __deal2__sparse_matrix_templates_h
#define __deal2__sparse_matrix_templates_h


#include <lac/sparse_matrix.h>
#include <lac/vector.h>


#include <iostream>
#include <iomanip>
#include <algorithm>
#include <functional>
#include <cmath>

#ifdef DEAL_II_USE_MT
#  include <vector>
#  include <numeric>

#  include <base/thread_management.h>
#  include <base/multithread_info.h>
#endif




template <typename number>
SparseMatrix<number>::SparseMatrix () :
		cols(0),
		val(0),
		max_len(0)
{};



template <typename number>
SparseMatrix<number>::SparseMatrix (const SparseMatrix &m) :
		Subscriptor (m),
		cols(0),
		val(0),
		max_len(0)
{
  Assert (m.cols==0, ExcInvalidConstructorCall());
  Assert (m.val==0, ExcInvalidConstructorCall());
  Assert (m.max_len==0, ExcInvalidConstructorCall());
};



template <typename number>
SparseMatrix<number>&
SparseMatrix<number>::operator = (const SparseMatrix<number> &m)
{
  Assert (m.cols==0, ExcInvalidConstructorCall());
  Assert (m.val==0, ExcInvalidConstructorCall());
  Assert (m.max_len==0, ExcInvalidConstructorCall());

  return *this;
};



template <typename number>
SparseMatrix<number>::SparseMatrix (const SparsityPattern &c) :
		cols(&c),
		val(0),
		max_len(0)
{
  reinit();
};



template <typename number>
SparseMatrix<number>::~SparseMatrix ()
{
  cols = 0;
  
  if (val != 0)
    delete[] val;
};



template <typename number>
void
SparseMatrix<number>::reinit ()
{
  Assert (cols != 0, ExcMatrixNotInitialized());
  Assert (cols->compressed || cols->empty(), ExcNotCompressed());

  if (cols->empty()) 
    {
      if (val) delete[] val;
      val = 0;
      max_len = 0;
      return;
    };

  const unsigned int N = cols->n_nonzero_elements();
  if (N > max_len)
    {
      if (val) delete[] val;
      val = new number[N];
      max_len = N;
    };

  if (val)
    std::fill_n (&val[0], N, 0);
}



template <typename number>
void
SparseMatrix<number>::reinit (const SparsityPattern &sparsity)
{
  cols = &sparsity;
  reinit ();
};



template <typename number>
void
SparseMatrix<number>::clear ()
{
  cols = 0;
  if (val) delete[] val;
  val = 0;
  max_len = 0;
};



template <typename number>
bool
SparseMatrix<number>::empty () const
{
  if (cols == 0)
    return true;
  else
    return cols->empty();
};



template <typename number>
unsigned int
SparseMatrix<number>::n_nonzero_elements () const
{
  Assert (cols != 0, ExcMatrixNotInitialized());
  return cols->n_nonzero_elements ();
};



template <typename number>
unsigned int
SparseMatrix<number>::n_actually_nonzero_elements () const
{
  Assert (cols != 0, ExcMatrixNotInitialized());
  return std::count_if(&val[0], &val[n_nonzero_elements ()],
		       std::bind2nd(std::not_equal_to<double>(), 0));
};



template <typename number>
void
SparseMatrix<number>::symmetrize ()
{
  Assert (cols != 0, ExcMatrixNotInitialized());
  Assert (cols->rows == cols->cols, ExcMatrixNotSquare());
  
  const unsigned int n_rows = m();
  for (unsigned int row=0; row<n_rows; ++row)
    {
				       // first skip diagonal entry
      number             *val_ptr = &val[cols->rowstart[row]+1];
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
};



template <typename number>
template <typename somenumber>
SparseMatrix<number> &
SparseMatrix<number>::copy_from (const SparseMatrix<somenumber> &matrix)
{
  Assert (cols != 0, ExcMatrixNotInitialized());
  Assert (val != 0, ExcMatrixNotInitialized());
  Assert (cols == matrix.cols, ExcDifferentSparsityPatterns());

  std::copy (&matrix.val[0], &matrix.val[cols->n_nonzero_elements()],
	     &val[0]);
  
  return *this;
};


template <typename number>
template <typename somenumber>
void
SparseMatrix<number>::add_scaled (const number factor,
				  const SparseMatrix<somenumber> &matrix)
{
  Assert (cols != 0, ExcMatrixNotInitialized());
  Assert (val != 0, ExcMatrixNotInitialized());
  Assert (cols == matrix.cols, ExcDifferentSparsityPatterns());

  number             *val_ptr    = &val[0];
  const somenumber   *matrix_ptr = &matrix.val[0];
  const number *const end_ptr    = &val[cols->n_nonzero_elements()];

  while (val_ptr != end_ptr)
    *val_ptr++ += factor * *matrix_ptr++;
};


template <typename number>
template <typename somenumber>
void
SparseMatrix<number>::vmult (Vector<somenumber>& dst, const Vector<somenumber>& src) const
{
  Assert (cols != 0, ExcMatrixNotInitialized());
  Assert (val != 0, ExcMatrixNotInitialized());
  Assert(m() == dst.size(), ExcDimensionsDontMatch(m(),dst.size()));
  Assert(n() == src.size(), ExcDimensionsDontMatch(n(),src.size()));
  
  const unsigned int n_rows = m();

#ifdef DEAL_II_USE_MT
				   // in MT mode: start new threads only
				   // if the matrix is sufficiently large.
				   // the limit is mostly artificial
  if (n_rows/multithread_info.n_default_threads > 2000)
    {
      const unsigned int n_threads = multithread_info.n_default_threads;

      Threads::ThreadManager thread_manager;
      for (unsigned int i=0; i<n_threads; ++i)
	Threads::spawn (thread_manager,
			Threads::encapsulate (&SparseMatrix<number>::
					      template threaded_vmult<somenumber>)
			.collect_args (this, dst, src,
				       n_rows * i / n_threads,
				       n_rows * (i+1) / n_threads));
      thread_manager.wait ();

      return;
    };
#endif

				   // if not in MT mode or size<2000
				   // do it in an oldfashioned way
  const number       *val_ptr    = &val[cols->rowstart[0]];
  const unsigned int *colnum_ptr = &cols->colnums[cols->rowstart[0]];
  somenumber         *dst_ptr    = &dst(0);
  for (unsigned int row=0; row<n_rows; ++row)
    {
      somenumber s = 0.;
      const number *const val_end_of_row = &val[cols->rowstart[row+1]];
      while (val_ptr != val_end_of_row)
	s += *val_ptr++ * src(*colnum_ptr++);
      *dst_ptr++ = s;
    };
};


template <typename number>
template <typename somenumber>
void
SparseMatrix<number>::threaded_vmult (Vector<somenumber>       &dst,
				      const Vector<somenumber> &src,
				      const unsigned int        begin_row,
				      const unsigned int        end_row) const
{
#ifdef DEAL_II_USE_MT
  const number       *val_ptr    = &val[cols->rowstart[begin_row]];
  const unsigned int *colnum_ptr = &cols->colnums[cols->rowstart[begin_row]];
  somenumber         *dst_ptr    = &dst(begin_row);
  for (unsigned int row=begin_row; row<end_row; ++row)
    {
      somenumber s = 0.;
      const number *const val_end_of_row = &val[cols->rowstart[row+1]];
      while (val_ptr != val_end_of_row)
	s += *val_ptr++ * src(*colnum_ptr++);
      *dst_ptr++ = s;
    };
#else
				   // this function should not be called
				   // when not in parallel mode.
  Assert (false, ExcInternalError());
#endif
};


template <typename number>
template <typename somenumber>
void
SparseMatrix<number>::Tvmult (Vector<somenumber>& dst, const Vector<somenumber>& src) const
{
  Assert (val != 0, ExcMatrixNotInitialized());
  Assert (cols != 0, ExcMatrixNotInitialized());
  Assert(n() == dst.size(), ExcDimensionsDontMatch(n(),dst.size()));
  Assert(m() == src.size(), ExcDimensionsDontMatch(m(),src.size()));

  dst.clear ();

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
template <typename somenumber>
void
SparseMatrix<number>::vmult_add (Vector<somenumber>& dst, const Vector<somenumber>& src) const
{
  Assert (cols != 0, ExcMatrixNotInitialized());
  Assert (val != 0, ExcMatrixNotInitialized());
  Assert(m() == dst.size(), ExcDimensionsDontMatch(m(),dst.size()));
  Assert(n() == src.size(), ExcDimensionsDontMatch(n(),src.size()));

  const unsigned int  n_rows     = m();
  const number       *val_ptr    = &val[cols->rowstart[0]];
  const unsigned int *colnum_ptr = &cols->colnums[cols->rowstart[0]];
  somenumber         *dst_ptr    = &dst(0);
  for (unsigned int row=0; row<n_rows; ++row)
    {
      somenumber s = 0.;
      const number *const val_end_of_row = &val[cols->rowstart[row+1]];
      while (val_ptr != val_end_of_row)
	s += *val_ptr++ * src(*colnum_ptr++);
      *dst_ptr++ += s;
    };
};


template <typename number>
template <typename somenumber>
void
SparseMatrix<number>::Tvmult_add (Vector<somenumber>& dst, const Vector<somenumber>& src) const
{
  Assert (val != 0, ExcMatrixNotInitialized());
  Assert (cols != 0, ExcMatrixNotInitialized());
  Assert(n() == dst.size(), ExcDimensionsDontMatch(n(),dst.size()));
  Assert(m() == src.size(), ExcDimensionsDontMatch(m(),src.size()));

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
template <typename somenumber>
somenumber
SparseMatrix<number>::matrix_norm_square (const Vector<somenumber>& v) const
{
  Assert (cols != 0, ExcMatrixNotInitialized());
  Assert (val != 0, ExcMatrixNotInitialized());
  Assert(m() == v.size(), ExcDimensionsDontMatch(m(),v.size()));
  Assert(n() == v.size(), ExcDimensionsDontMatch(n(),v.size()));

  const unsigned int n_rows = m();
#ifdef DEAL_II_USE_MT
				   // if in MT mode and size sufficiently
				   // large: do it in parallel; the limit
				   // is mostly artificial
  if (n_rows/multithread_info.n_default_threads > 2000)
    {
      const unsigned int n_threads = multithread_info.n_default_threads;

				       // space for the norms of
				       // the different parts
      std::vector<somenumber> partial_sums (n_threads, 0);
      Threads::ThreadManager thread_manager;
				       // spawn some jobs...
      for (unsigned int i=0; i<n_threads; ++i)
	Threads::spawn (thread_manager,
			Threads::encapsulate (&SparseMatrix<number>::
					      template threaded_matrix_norm_square<somenumber>)
			.collect_args (this, v,
				       n_rows * i / n_threads,
				       n_rows * (i+1) / n_threads,
				       &partial_sums[i]));

				       // ... and wait until they're finished
      thread_manager.wait ();
				       // accumulate the partial results
      return accumulate (partial_sums.begin(),
			 partial_sums.end(),
			 0.);
    };
#endif
				   // if not in MT mode or the matrix is
				   // too small: do it one-by-one
  somenumber          sum        = 0.;
  const number       *val_ptr    = &val[cols->rowstart[0]];
  const unsigned int *colnum_ptr = &cols->colnums[cols->rowstart[0]];
  for (unsigned int row=0; row<n_rows; ++row)
    {
      somenumber s = 0.;
      const number *val_end_of_row = &val[cols->rowstart[row+1]];
      while (val_ptr != val_end_of_row)
	s += *val_ptr++ * v(*colnum_ptr++);

      sum += s* v(row);
    };

  return sum;
};



template <typename number>
template <typename somenumber>
void
SparseMatrix<number>::threaded_matrix_norm_square (const Vector<somenumber> &v,
						   const unsigned int        begin_row,
						   const unsigned int        end_row,
						   somenumber               *partial_sum) const
{
#ifdef DEAL_II_USE_MT
  somenumber sum = 0.;
  const number       *val_ptr    = &val[cols->rowstart[begin_row]];
  const unsigned int *colnum_ptr = &cols->colnums[cols->rowstart[begin_row]];
  for (unsigned int row=begin_row; row<end_row; ++row)
    {
      somenumber s = 0.;
      const number *val_end_of_row = &val[cols->rowstart[row+1]];
      while (val_ptr != val_end_of_row)
	s += *val_ptr++ * v(*colnum_ptr++);

      sum += s* v(row);
    };
  *partial_sum = sum;
  
#else
				   // function should not have been called
  Assert (false, ExcInternalError());
#endif
};



template <typename number>
template <typename somenumber>
somenumber
SparseMatrix<number>::matrix_scalar_product (const Vector<somenumber>& u,
					     const Vector<somenumber>& v) const
{
  Assert (cols != 0, ExcMatrixNotInitialized());
  Assert (val != 0, ExcMatrixNotInitialized());
  Assert(m() == u.size(), ExcDimensionsDontMatch(m(),u.size()));
  Assert(n() == v.size(), ExcDimensionsDontMatch(n(),v.size()));

  const unsigned int n_rows = m();
#ifdef DEAL_II_USE_MT
				   // if in MT mode and size sufficiently
				   // large: do it in parallel; the limit
				   // is mostly artificial
  if (n_rows/multithread_info.n_default_threads > 2000)
    {
      const unsigned int n_threads = multithread_info.n_default_threads;

				       // space for the norms of
				       // the different parts
      std::vector<somenumber> partial_sums (n_threads, 0);
      Threads::ThreadManager thread_manager;
				       // spawn some jobs...
      for (unsigned int i=0; i<n_threads; ++i)
	Threads::spawn (thread_manager,
			Threads::encapsulate (&SparseMatrix<number>::
					      template threaded_matrix_scalar_product<somenumber>)
			.collect_args (this, u, v,
				       n_rows * i / n_threads,
				       n_rows * (i+1) / n_threads,
				       &partial_sums[i]));

				       // ... and wait until they're finished
      thread_manager.wait ();
				       // accumulate the partial results
      return accumulate (partial_sums.begin(),
			 partial_sums.end(),
			 0.);
    };
#endif
				   // if not in MT mode or the matrix is
				   // too small: do it one-by-one
  somenumber          sum        = 0.;
  const number       *val_ptr    = &val[cols->rowstart[0]];
  const unsigned int *colnum_ptr = &cols->colnums[cols->rowstart[0]];
  for (unsigned int row=0; row<n_rows; ++row)
    {
      somenumber s = 0.;
      const number *val_end_of_row = &val[cols->rowstart[row+1]];
      while (val_ptr != val_end_of_row)
	s += *val_ptr++ * v(*colnum_ptr++);

      sum += s* u(row);
    };

  return sum;
};



template <typename number>
template <typename somenumber>
void
SparseMatrix<number>::threaded_matrix_scalar_product (const Vector<somenumber> &u,
						      const Vector<somenumber> &v,
						      const unsigned int        begin_row,
						      const unsigned int        end_row,
						      somenumber               *partial_sum) const
{
#ifdef DEAL_II_USE_MT
  somenumber sum = 0.;
  const number       *val_ptr    = &val[cols->rowstart[begin_row]];
  const unsigned int *colnum_ptr = &cols->colnums[cols->rowstart[begin_row]];
  for (unsigned int row=begin_row; row<end_row; ++row)
    {
      somenumber s = 0.;
      const number *val_end_of_row = &val[cols->rowstart[row+1]];
      while (val_ptr != val_end_of_row)
	s += *val_ptr++ * v(*colnum_ptr++);

      sum += s* u(row);
    };
  *partial_sum = sum;
  
#else
				   // function should not have been called
  Assert (false, ExcInternalError());
#endif
};





template <typename number>
number SparseMatrix<number>::l1_norm () const
{
  Assert (cols != 0, ExcMatrixNotInitialized());
  Assert (val != 0, ExcMatrixNotInitialized());

  Vector<number> column_sums(n());
  const unsigned int n_rows = m();
  for (unsigned int row=0; row<n_rows; ++row)
    for (unsigned int j=cols->rowstart[row]; j<cols->rowstart[row+1] ; ++j)
      column_sums(cols->colnums[j]) += std::fabs(val[j]);

  return column_sums.linfty_norm();
};


template <typename number>
number SparseMatrix<number>::linfty_norm () const
{
  Assert (cols != 0, ExcMatrixNotInitialized());
  Assert (val != 0, ExcMatrixNotInitialized());

  const number *val_ptr = &val[cols->rowstart[0]];

  number sum, max=0;
  const unsigned int n_rows = m();
  for (unsigned int row=0; row<n_rows; ++row)
    {
      sum=0;
      const number *const val_end_of_row = &val[cols->rowstart[row+1]];
      while (val_ptr != val_end_of_row)
	sum += std::fabs(*val_ptr++);
      if (sum > max)
	max = sum;
    }
  return max;
};


template <typename number>
template <typename somenumber>
somenumber
SparseMatrix<number>::residual (Vector<somenumber>       &dst,
				const Vector<somenumber> &u,
				const Vector<somenumber> &b) const
{
  Assert (cols != 0, ExcMatrixNotInitialized());
  Assert (val != 0, ExcMatrixNotInitialized());
  Assert(m() == dst.size(), ExcDimensionsDontMatch(m(),dst.size()));
  Assert(m() == b.size(), ExcDimensionsDontMatch(m(),b.size()));
  Assert(n() == u.size(), ExcDimensionsDontMatch(n(),u.size()));

  const unsigned int n_rows = m();
#ifdef DEAL_II_USE_MT
				   // if in MT mode and size sufficiently
				   // large: do it in parallel; the limit
				   // is mostly artificial
  if (n_rows/multithread_info.n_default_threads > 2000)
    {
      const unsigned int n_threads = multithread_info.n_default_threads;
 
				       // space for the square norms of
				       // the different parts
      std::vector<somenumber> partial_norms (n_threads, 0);
      Threads::ThreadManager thread_manager;
      for (unsigned int i=0; i<n_threads; ++i)
	Threads::spawn (thread_manager,
			Threads::encapsulate (&SparseMatrix<number>::
					      template threaded_residual<somenumber>)
			.collect_args (this, dst, u, b,
				       make_pair<unsigned int,unsigned int>
				       (n_rows * i / n_threads,
					n_rows * (i+1) / n_threads),
				       &partial_norms[i]));

				       // ... and wait until they're finished
      thread_manager.wait ();
				       // accumulate the partial results
      return sqrt(accumulate (partial_norms.begin(),
			      partial_norms.end(),
			      0.));
    };
#endif
  
  somenumber norm=0.;   
  
  for (unsigned int i=0; i<n_rows; ++i)
    {
      somenumber s = b(i);
      for (unsigned int j=cols->rowstart[i]; j<cols->rowstart[i+1] ;j++)
	{
	  const unsigned int p = cols->colnums[j];
	  s -= val[j] * u(p);
	}
      dst(i) = s;
      norm += dst(i)*dst(i);
    }
  return std::sqrt(norm);
}


template <typename number>
template <typename somenumber>
void
SparseMatrix<number>::threaded_residual (Vector<somenumber>       &dst,
					 const Vector<somenumber> &u,
					 const Vector<somenumber> &b,
					 const std::pair<unsigned int, unsigned int> interval,
					 somenumber               *partial_norm) const
{
  const unsigned int begin_row = interval.first,
		     end_row   = interval.second;
  
#ifdef DEAL_II_USE_MT
  somenumber norm=0.;   
  
  for (unsigned int i=begin_row; i<end_row; ++i)
    {
      somenumber s = b(i);
      for (unsigned int j=cols->rowstart[i]; j<cols->rowstart[i+1] ;j++)
	{
	  const unsigned int p = cols->colnums[j];
	  s -= val[j] * u(p);
	}
      dst(i) = s;
      norm += dst(i)*dst(i);
    };

  *partial_norm = norm;
#else
  Assert (false, ExcInternalError());
#endif
};



template <typename number>
template <typename somenumber>
void
SparseMatrix<number>::precondition_Jacobi (Vector<somenumber>& dst,
					   const Vector<somenumber>& src,
					   const number om) const
{
  Assert (cols != 0, ExcMatrixNotInitialized());
  Assert (val != 0, ExcMatrixNotInitialized());
  Assert (m() == n(), ExcMatrixNotSquare());

  const unsigned int n = src.size();
  somenumber              *dst_ptr = dst.begin();
  const somenumber        *src_ptr = src.begin();
  const unsigned int *rowstart_ptr = &cols->rowstart[0];
  
  for (unsigned int i=0; i<n; ++i, ++dst_ptr, ++src_ptr, ++rowstart_ptr)
				     // note that for square matrices,
				     // the diagonal entry is the first
				     // in each row, i.e. at index
				     // rowstart[i]
    *dst_ptr = om * *src_ptr / val[*rowstart_ptr];
};



template <typename number>
template <typename somenumber>
void
SparseMatrix<number>::precondition_SSOR (Vector<somenumber>& dst,
					 const Vector<somenumber>& src,
					 const number om) const
{
				   // to understand how this function works
				   // you may want to take a look at the CVS
				   // archives to see the original version
				   // which is much clearer...
  Assert (cols != 0, ExcMatrixNotInitialized());
  Assert (val != 0, ExcMatrixNotInitialized());
  Assert (m() == n(), ExcMatrixNotSquare());

  const unsigned int  n            = src.size();
  const unsigned int *rowstart_ptr = &cols->rowstart[0];
  somenumber         *dst_ptr      = &dst(0);

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
	= (SparsityPattern::optimized_lower_bound (&cols->colnums[*rowstart_ptr+1],
						   &cols->colnums[*(rowstart_ptr+1)],
						   row)
	   -
	   &cols->colnums[0]);
      
      for (unsigned int j=(*rowstart_ptr)+1; j<first_right_of_diagonal_index; ++j)
	*dst_ptr -= om* val[j] * dst(cols->colnums[j]);

				       // divide by diagonal element
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
      const unsigned int first_right_of_diagonal_index
	= (SparsityPattern::optimized_lower_bound (&cols->colnums[*rowstart_ptr+1],
						   &cols->colnums[*(rowstart_ptr+1)],
						   static_cast<unsigned int>(row)) -
	   &cols->colnums[0]);
      for (unsigned int j=first_right_of_diagonal_index; j<*(rowstart_ptr+1); ++j)
	if (cols->colnums[j] > static_cast<unsigned int>(row))
	  *dst_ptr -= om* val[j] * dst(cols->colnums[j]);
      
      *dst_ptr /= val[*rowstart_ptr];
    };
}


template <typename number>
template <typename somenumber>
void
SparseMatrix<number>::precondition_SOR (Vector<somenumber>& dst, const Vector<somenumber>& src,
			    const number om) const
{
  Assert (cols != 0, ExcMatrixNotInitialized());
  Assert (val != 0, ExcMatrixNotInitialized());
  Assert (m() == n(), ExcMatrixNotSquare());

  dst = src;
  SOR(dst,om);
};


template <typename number>
template <typename somenumber>
void
SparseMatrix<number>::precondition_TSOR (Vector<somenumber>& dst, const Vector<somenumber>& src,
			    const number om) const
{
  Assert (cols != 0, ExcMatrixNotInitialized());
  Assert (val != 0, ExcMatrixNotInitialized());
  Assert (m() == n(), ExcMatrixNotSquare());

  dst = src;
  TSOR(dst,om);
};


template <typename number>
template <typename somenumber>
void
SparseMatrix<number>::SOR (Vector<somenumber>& dst, const number om) const
{
  Assert (cols != 0, ExcMatrixNotInitialized());
  Assert (val != 0, ExcMatrixNotInitialized());
  Assert (m() == n(), ExcMatrixNotSquare());
  Assert (m() == dst.size(), ExcDimensionsDontMatch(m(),dst.size()));

  for (unsigned int row=0; row<m(); ++row)
    {
      somenumber s = dst(row);
      for (unsigned int j=cols->rowstart[row]; j<cols->rowstart[row+1]; ++j)
	if (cols->colnums[j] < row)
	  s -= val[j] * dst(cols->colnums[j]);

      dst(row) = s * om / val[cols->rowstart[row]];
    }
}


template <typename number>
template <typename somenumber>
void
SparseMatrix<number>::TSOR (Vector<somenumber>& dst, const number om) const
{
  Assert (cols != 0, ExcMatrixNotInitialized());
  Assert (val != 0, ExcMatrixNotInitialized());
  Assert (m() == n(), ExcMatrixNotSquare());
  Assert (m() == dst.size(), ExcDimensionsDontMatch(m(),dst.size()));

  for (unsigned int row=m(); row!=0;)
    {
      --row;
      somenumber s = dst(row);
      for (unsigned int j=cols->rowstart[row]; j<cols->rowstart[row+1]; ++j)
	if (cols->colnums[j] > row)
	  s -= val[j] * dst(cols->colnums[j]);

      dst(row) = s * om / val[cols->rowstart[row]];
    }
}


template <typename number>
template <typename somenumber>
void
SparseMatrix<number>::SOR_step (Vector<somenumber> &v,
	       const Vector<somenumber> &b,
	       const number        om) const
{
  Assert (cols != 0, ExcMatrixNotInitialized());
  Assert (val != 0, ExcMatrixNotInitialized());
  Assert (m() == n(), ExcMatrixNotSquare());
  Assert (m() == v.size(), ExcDimensionsDontMatch(m(),v.size()));
  Assert (m() == b.size(), ExcDimensionsDontMatch(m(),b.size()));

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
  Assert (cols != 0, ExcMatrixNotInitialized());
  Assert (val != 0, ExcMatrixNotInitialized());
  Assert (m() == n(), ExcMatrixNotSquare());
  Assert (m() == v.size(), ExcDimensionsDontMatch(m(),v.size()));
  Assert (m() == b.size(), ExcDimensionsDontMatch(m(),b.size()));

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
SparseMatrix<number>::SSOR (Vector<somenumber>& dst, const number om) const
{
  Assert (cols != 0, ExcMatrixNotInitialized());
  Assert (val != 0, ExcMatrixNotInitialized());
  Assert (m() == n(), ExcMatrixNotSquare());
  Assert (m() == dst.size(), ExcDimensionsDontMatch(m(),dst.size()));

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

  for (int i=n-1; i>=0; i--)  // this time, i is signed, but alsways positive!
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
  Assert (cols != 0, ExcMatrixNotInitialized());
  return *cols;
};



template <typename number>
void SparseMatrix<number>::print (std::ostream &out) const
{
  Assert (cols != 0, ExcMatrixNotInitialized());
  Assert (val != 0, ExcMatrixNotInitialized());

  for (unsigned int i=0; i<cols->rows; ++i)
    for (unsigned int j=cols->rowstart[i]; j<cols->rowstart[i+1]; ++j)
      out << "(" << i << "," << cols->colnums[j] << ") " << val[j] << std::endl;

  AssertThrow (out, ExcIO());
};


template <typename number>
void SparseMatrix<number>::print_formatted (std::ostream &out,
					    const unsigned int precision,
					    const bool scientific,
					    const unsigned int width_,
					    const char* zero_string,
					    const double denominator) const
{
  Assert (cols != 0, ExcMatrixNotInitialized());
  Assert (val != 0, ExcMatrixNotInitialized());

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
};



template <typename number>
unsigned int
SparseMatrix<number>::memory_consumption () const
{
  return sizeof(*this) + max_len*sizeof(number);
};


#endif
