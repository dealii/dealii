//----------------------------  filtered_matrix.templates.h  ---------------------------
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
//----------------------------  filtered_matrix.templates.h  ---------------------------
#ifndef __deal2__filtered_matrix_templates_h
#define __deal2__filtered_matrix_templates_h


#include <base/config.h>
#include <base/memory_consumption.h>
#include <lac/filtered_matrix.h>
#include <lac/sparse_matrix.h>
#include <lac/block_sparse_matrix.h>
#include <lac/vector.h>
#include <lac/block_vector.h>


template <class Matrix, class Vector>
FilteredMatrix<Matrix,Vector>::
FilteredMatrix () 
{};



template <class Matrix, class Vector>
FilteredMatrix<Matrix,Vector>::
FilteredMatrix (const FilteredMatrix &fm)
		:
		Subscriptor (),
		constraints (fm.constraints)
{
  set_referenced_matrix (*fm.matrix);
};



template <class Matrix, class Vector>
FilteredMatrix<Matrix,Vector>::
FilteredMatrix (const Matrix &m)
{
  set_referenced_matrix (m);
};



template <class Matrix, class Vector>
FilteredMatrix<Matrix,Vector> &
FilteredMatrix<Matrix,Vector>::operator = (const FilteredMatrix &fm)
{
  set_referenced_matrix (*fm.matrix);
  constraints = fm.constraints;
  return *this;
};



template <class Matrix, class Vector>
void
FilteredMatrix<Matrix,Vector>::
set_referenced_matrix (const Matrix &m)
{
  matrix = &m;
  allocate_tmp_vector ();
};



template <class Matrix, class Vector>
void
FilteredMatrix<Matrix,Vector>::clear_constraints () 
{
				   // swap vectors to release memory
  std::vector<IndexValuePair> empty;
  constraints.swap (empty);
};



template <class Matrix, class Vector>
void
FilteredMatrix<Matrix,Vector>::
apply_constraints (Vector     &v,
		   const bool  matrix_is_symmetric) const
{
				   // array that will hold the pairs
				   // of index/value of all nonzero
				   // entries in a given column
  std::vector<IndexValuePair> column_entries;
  
				   // iterate over all constraints and
				   // treat them one after the other
  const_index_value_iterator       i = constraints.begin();
  const const_index_value_iterator e = constraints.end();
  for (; i!=e; ++i)
    {
				       // define abbreviations
      const unsigned   index = i->first;
      const value_type value = i->second;
      
				       // check whether the value is
				       // zero, since in that case we do
				       // not have to modify other nodes
      if (value != 0)
	{
					   // first clear array of
					   // previous content
	  column_entries.clear ();
	  
					   // then get all entries in
					   // the present column
	  get_column_entries (index, column_entries, matrix_is_symmetric);
	  
					   // modify rhs for each entry
	  const_index_value_iterator       col     = column_entries.begin();
	  const const_index_value_iterator col_end = column_entries.end();
	  for (; col!=col_end; ++col)
	    v(col->first) -= col->second * value;
	};
    };

  
				       // finally set constrained
				       // entries themselves. we can't
				       // do it in the above loop
				       // since we might end up
				       // modifying an entry that we
				       // have already set if
				       // constrained dofs couple to
				       // each other
  for (i=constraints.begin(); i!=e; ++i)
    v(i->first) = i->second;
};



template <class Matrix, class Vector>
void
FilteredMatrix<Matrix,Vector>::pre_filter (Vector &v) const
{
				   // iterate over all constraints and
				   // zero out value
  const_index_value_iterator       i = constraints.begin();
  const const_index_value_iterator e = constraints.end();
  for (; i!=e; ++i)
    v(i->first) = 0;
};



template <class Matrix, class Vector>
void
FilteredMatrix<Matrix,Vector>::post_filter (const Vector &in,
					    Vector       &out) const
{
				   // iterate over all constraints and
				   // set value correctly
  const_index_value_iterator       i = constraints.begin();
  const const_index_value_iterator e = constraints.end();
  for (; i!=e; ++i)
    out(i->first) = in(i->first);
};



template <class Matrix, class Vector>
void
FilteredMatrix<Matrix,Vector>::vmult (Vector       &dst,
				      const Vector &src) const
{
  tmp_mutex.acquire ();
				   // first copy over src vector and
				   // pre-filter
  tmp_vector = src;
  pre_filter (tmp_vector);
				   // then let matrix do its work
  matrix->vmult (dst, tmp_vector);
				   // tmp_vector now no more needed
  tmp_mutex.release ();
				   // finally do post-filtering
  post_filter (src, dst);
};



template <class Matrix, class Vector>
typename FilteredMatrix<Matrix,Vector>::value_type
FilteredMatrix<Matrix,Vector>::residual (Vector       &dst,
					 const Vector &x,
					 const Vector &b) const
{
  tmp_mutex.acquire ();
				   // first copy over x vector and
				   // pre-filter
  tmp_vector = x;
  pre_filter (tmp_vector);
				   // then let matrix do its work
  value_type res  = matrix->residual (dst, tmp_vector, b);
  value_type res2 = res*res;
				   // tmp_vector now no more needed
  tmp_mutex.release ();
				   // finally do post-filtering. here,
				   // we set constrained indices to
				   // zero, but have to subtract their
				   // contributions to the residual
  const_index_value_iterator       i = constraints.begin();
  const const_index_value_iterator e = constraints.end();
  for (; i!=e; ++i)
    {
      const value_type v = dst(i->first);
      res2 -= v*v;
      dst(i->first) = 0;
    };
  
  Assert (res2>=0, ExcInternalError());
  return std::sqrt (res2);
};



template <class Matrix, class Vector>
void
FilteredMatrix<Matrix,Vector>::Tvmult (Vector       &dst,
				       const Vector &src) const
{
  tmp_mutex.acquire ();
				   // first copy over src vector and
				   // pre-filter
  tmp_vector = src;
  pre_filter (tmp_vector);
				   // then let matrix do its work
  matrix->Tvmult (dst, tmp_vector);
				   // tmp_vector now no more needed
  tmp_mutex.release ();
				   // finally do post-filtering
  post_filter (src, dst);
};

  

template <class Matrix, class Vector>
typename FilteredMatrix<Matrix,Vector>::value_type
FilteredMatrix<Matrix,Vector>::matrix_norm_square (const Vector &v) const
{
  tmp_mutex.acquire ();
  tmp_vector = v;

				   // zero out constrained entries and
				   // form matrix norm with original
				   // matrix. this is equivalent to
				   // forming the matrix norm of the
				   // original vector with the matrix
				   // where we have zeroed out rows
				   // and columns
  pre_filter (tmp_vector);
  const value_type ret = matrix->matrix_norm_square (tmp_vector);
  tmp_mutex.release ();
  return ret;
};



template <class Matrix, class Vector>
void
FilteredMatrix<Matrix,Vector>::
precondition_Jacobi (Vector           &dst,
		     const Vector     &src,
		     const value_type  omega) const
{
				   // first precondition as usual,
				   // using the fast algorithms of the
				   // matrix class
  matrix->precondition_Jacobi (dst, src, omega);

				   // then modify the constrained
				   // degree of freedom. as the
				   // diagonal entries of the filtered
				   // matrix would be 1.0, simply copy
				   // over old and new values
  const_index_value_iterator       i = constraints.begin();
  const const_index_value_iterator e = constraints.end();
  for (; i!=e; ++i)
    dst(i->first) = src(i->first);
};



template <class Matrix, class Vector>
unsigned int
FilteredMatrix<Matrix,Vector>::memory_consumption () const
{
  return (MemoryConsumption::memory_consumption (matrix) +
	  MemoryConsumption::memory_consumption (constraints) +
	  MemoryConsumption::memory_consumption (tmp_vector));
};



#endif
