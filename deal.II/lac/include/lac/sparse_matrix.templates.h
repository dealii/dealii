// $Id$

// This file was once part of the DEAL Library
// DEAL is Copyright(1995) by
// Roland Becker, Guido Kanschat, Franz-Theo Suttmeier
// Revised, modified and extended by Wolfgang Bangerth, 1998, 1999


#include <lac/sparsematrix.h>
#include <lac/vector.h>


#include <iostream>
#include <iomanip>
#include <algorithm>




template <typename number>
SparseMatrix<number>::SparseMatrix () :
		cols(0),
		val(0),
		max_len(0) {};



template <typename number>
SparseMatrix<number>::SparseMatrix (const SparseMatrixStruct &c)
		: cols(&c), val(0), max_len(0)
{
  reinit();
};



template <typename number>
SparseMatrix<number>::~SparseMatrix ()
{
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
        
  if (max_len<cols->vec_len)
    {
      if (val) delete[] val;
      val = new number[cols->vec_len];
      max_len = cols->vec_len;
    };

  if (val)
    fill_n (&val[0], cols->vec_len, 0);
}



template <typename number>
void
SparseMatrix<number>::reinit (const SparseMatrixStruct &sparsity) {
  cols = &sparsity;
  reinit ();
};



template <typename number>
void
SparseMatrix<number>::clear () {
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
SparseMatrix<number>::n_nonzero_elements () const {
  Assert (cols != 0, ExcMatrixNotInitialized());
  return cols->n_nonzero_elements ();
};



template <typename number>
template <typename somenumber>
SparseMatrix<number> &
SparseMatrix<number>::copy_from (const SparseMatrix<somenumber> &matrix) {
  Assert (cols != 0, ExcMatrixNotInitialized());
  Assert (val != 0, ExcMatrixNotInitialized());
  Assert (cols == matrix.cols, ExcDifferentSparsityPatterns());

  number                 *val_ptr = &val[0];
  const somenumber    *matrix_ptr = &matrix.val[0];
  const number     *const end_ptr = &val[cols->vec_len];

  while (val_ptr != end_ptr)
    *val_ptr++ = *matrix_ptr++;
  
  return *this;
};



template <typename number>
template <typename somenumber>
void
SparseMatrix<number>::add_scaled (const number factor,
				  const SparseMatrix<somenumber> &matrix) {
  Assert (cols != 0, ExcMatrixNotInitialized());
  Assert (val != 0, ExcMatrixNotInitialized());
  Assert (cols == matrix.cols, ExcDifferentSparsityPatterns());

  number             *val_ptr    = &val[0];
  const somenumber   *matrix_ptr = &matrix.val[0];
  const number *const end_ptr    = &val[cols->vec_len];

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
  const number *val_ptr = &val[cols->rowstart[0]];
  const int *colnum_ptr = &cols->colnums[cols->rowstart[0]];
  somenumber   *dst_ptr = &dst(0);
  for (unsigned int row=0; row<n_rows; ++row)
    {
      double s = 0.;
      const number *const val_end_of_row = &val[cols->rowstart[row+1]];
      while (val_ptr != val_end_of_row)
	s += *val_ptr++ * src(*colnum_ptr++);
      *dst_ptr++ = s;
    };
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
	  int p = cols->colnums[j];
	  dst(p) += val[j] * src(i);
	}
    }
}



template <typename number>
template <typename somenumber>
double
SparseMatrix<number>::matrix_norm (const Vector<somenumber>& v) const
{
  Assert (cols != 0, ExcMatrixNotInitialized());
  Assert (val != 0, ExcMatrixNotInitialized());
  Assert(m() == v.size(), ExcDimensionsDontMatch(m(),v.size()));
  Assert(n() == v.size(), ExcDimensionsDontMatch(n(),v.size()));

  double sum = 0.;
  const unsigned int n_rows = m();
  const number *val_ptr = &val[cols->rowstart[0]];
  const int *colnum_ptr = &cols->colnums[cols->rowstart[0]];
  for (unsigned int row=0; row<n_rows; ++row)
    {
      double s = 0.;
      const number *val_end_of_row = &val[cols->rowstart[row+1]];
      while (val_ptr != val_end_of_row)
	s += *val_ptr++ * v(*colnum_ptr++);

      sum += s* v(row);
    };

  return sum;
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
      column_sums(cols->colnums[j])+=fabs(val[j]);

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
	sum += fabs(*val_ptr++);
      if (sum > max)
	max = sum;
    }
  return max;
};


template <typename number>
template <typename somenumber>
double
SparseMatrix<number>::residual (Vector<somenumber>& dst, const Vector<somenumber>& u, const Vector<somenumber>& b) const
{
  Assert (cols != 0, ExcMatrixNotInitialized());
  Assert (val != 0, ExcMatrixNotInitialized());
  Assert(m() == dst.size(), ExcDimensionsDontMatch(m(),dst.size()));
  Assert(m() == b.size(), ExcDimensionsDontMatch(m(),b.size()));
  Assert(n() == u.size(), ExcDimensionsDontMatch(n(),u.size()));

  double s,norm=0.;   
  
  for (unsigned int i=0;i<m();i++)
    {
      s = b(i);
      for (unsigned int j=cols->rowstart[i]; j<cols->rowstart[i+1] ;j++)
	{
	  int p = cols->colnums[j];
	  s -= val[j] * u(p);
	}
      dst(i) = s;
      norm += dst(i)*dst(i);
    }
  return sqrt(norm);
}



template <typename number>
template <typename somenumber>
void
SparseMatrix<number>::precondition_Jacobi (Vector<somenumber>& dst, const Vector<somenumber>& src,
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
SparseMatrix<number>::precondition_SSOR (Vector<somenumber>& dst, const Vector<somenumber>& src,
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
	= (lower_bound (&cols->colnums[*rowstart_ptr+1],
			&cols->colnums[*(rowstart_ptr+1)],
			static_cast<signed int>(row)) -
	   &cols->colnums[0]);
				       
      for (unsigned int j=(*rowstart_ptr)+1; j<first_right_of_diagonal_index; ++j)
	*dst_ptr -= om* val[j] * dst(cols->colnums[j]);
      *dst_ptr /= val[*rowstart_ptr];
    };
  
  rowstart_ptr = &cols->rowstart[0];
  dst_ptr      = &dst(0);
  for (unsigned int row=0; row<n; ++row, ++rowstart_ptr, ++dst_ptr)
    *dst_ptr *= (2.-om)*val[*rowstart_ptr];

  rowstart_ptr = &cols->rowstart[n-1];
  dst_ptr      = &dst(n-1);
  for (int row=n-1; row>=0; --row, --rowstart_ptr, --dst_ptr)
    {
      const unsigned int first_right_of_diagonal_index
	= (lower_bound (&cols->colnums[*rowstart_ptr+1],
			&cols->colnums[*(rowstart_ptr+1)],
			static_cast<signed int>(row)) -
	   &cols->colnums[0]);
      for (unsigned int j=first_right_of_diagonal_index; j<*(rowstart_ptr+1); ++j)
	if (cols->colnums[j] > row)
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
void SparseMatrix<number>::precondition (Vector<somenumber> &dst, const Vector<somenumber> &src) const {
  Assert (m() == n(), ExcMatrixNotSquare());
  dst=src;
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
	if ((unsigned int)cols->colnums[j] < row)
	  s -= val[j] * dst(cols->colnums[j]);

      dst(row) = s * om / val[cols->rowstart[row]];
    }
}



template <typename number>
template <typename somenumber>
void
SparseMatrix<number>::SSOR (Vector<somenumber>& dst, const number om) const
{
  Assert (cols != 0, ExcMatrixNotInitialized());
  Assert (val != 0, ExcMatrixNotInitialized());

  int p;
  const unsigned int  n = dst.size();
  unsigned int  j;
  double s;
  
  for (unsigned int i=0; i<n; i++)
    {
      s = 0.;
      for (j=cols->rowstart[i]; j<cols->rowstart[i+1] ;j++)
	{
	  p = cols->colnums[j];
	  if (p>=0)
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
	  p = cols->colnums[j];
	  if (p>=0)
	    {
	      if ((unsigned int)i<j) s += val[j] * dst(p);
	    }
	}
      dst(i) -= s * om / val[cols->rowstart[i]];
    }
}



template <typename number>
const SparseMatrixStruct & SparseMatrix<number>::get_sparsity_pattern () const {
  Assert (cols != 0, ExcMatrixNotInitialized());
  return *cols;
};



template <typename number>
void SparseMatrix<number>::print (ostream &out) const {
  Assert (cols != 0, ExcMatrixNotInitialized());
  Assert (val != 0, ExcMatrixNotInitialized());

  for (unsigned int i=0; i<cols->rows; ++i)
    for (unsigned int j=cols->rowstart[i]; j<cols->rowstart[i+1]; ++j)
      out << "(" << i << "," << cols->colnums[j] << ") " << val[j] << endl;

  AssertThrow (out, ExcIO());
};



template <typename number>
void SparseMatrix<number>::print_formatted (ostream &out, const unsigned int precision) const {
  Assert (cols != 0, ExcMatrixNotInitialized());
  Assert (val != 0, ExcMatrixNotInitialized());
  out.precision (precision);
  out.setf (ios::scientific, ios::floatfield);   // set output format
  
  for (unsigned int i=0; i<m(); ++i) 
    {
      for (unsigned int j=0; j<n(); ++j)
	if ((*cols)(i,j) != -1)
	  out << setw(precision+7)
	      << val[cols->operator()(i,j)] << ' ';
	else
	  out << setw(precision+8) << " ";
      out << endl;
    };
  AssertThrow (out, ExcIO());

  out.setf (0, ios::floatfield);                 // reset output format
};

