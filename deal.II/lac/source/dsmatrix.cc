// $Id$

// This file is part of the DEAL Library
// DEAL is Copyright(1995) by
// Roland Becker, Guido Kanschat, Franz-Theo Suttmeier

#include <lac/dsmatrix.h>
#include <iostream>
#include <iomanip>
#include <algorithm>



dSMatrixStruct::dSMatrixStruct (const unsigned int m, const unsigned int n,
				const unsigned int max_per_row) 
		: max_dim(0),
		  max_vec_len(0),
		  rowstart(0),
		  colnums(0)
{
  reinit (m,n,max_per_row);
};



dSMatrixStruct::dSMatrixStruct (const unsigned int n,
				const unsigned int max_per_row)
		: max_dim(0),
		  max_vec_len(0),
		  rowstart(0),
		  colnums(0)
{
  reinit (n,n,max_per_row);
};



dSMatrixStruct::~dSMatrixStruct ()
{
  if (rowstart != 0)
    delete[] rowstart;
  if (colnums != 0)
    delete[] colnums;
}




void
dSMatrixStruct::reinit (const unsigned int m, const unsigned int n,
			const unsigned int max_per_row)
{
  Assert (m>0, ExcInvalidNumber(m));
  Assert (n>0, ExcInvalidNumber(n));
  Assert (max_per_row>0, ExcInvalidNumber(max_per_row));
  rows = m;
  cols = n;
  vec_len = m * max_per_row;
  max_row_len = max_per_row;

  if (rows > max_dim)
    {
      if (rowstart) delete[] rowstart;
      max_dim = rows;
      rowstart = new unsigned int[max_dim+1];
    }

  if (vec_len > max_vec_len)
    {
      if (colnums) delete[] colnums;
      max_vec_len = vec_len;
      colnums = new int[max_vec_len];
    }

  for (unsigned int i=0; i<=rows; i++)
    rowstart[i] = i * max_per_row;
  fill_n (&colnums[0], vec_len, -1);

  if (rows == cols)
    for (unsigned int i=0;i<rows;i++)
      colnums[rowstart[i]] = i;

  compressed = false;
}


void
dSMatrixStruct::compress ()
{
  if (compressed) return;
  unsigned int next_free_entry = 0,
		next_row_start = 0,
		    row_length = 0;

				   // reserve temporary storage to
				   // store the entries of one wor
  int *tmp_entries = new int[max_row_len];

				   // Traverse all rows
  for (unsigned int line=0; line<rows; ++line)
    {
				       // copy used entries, break if
				       // first unused entry is reached
      row_length = 0;
      for (unsigned int j=rowstart[line]; j<rowstart[line+1]; ++j,++row_length)
	if (colnums[j] != -1)
	  tmp_entries[row_length] = colnums[j];
	else
	  break;
				       // now #rowstart# is
				       // the number of entries in
				       // this line

				       // for square matrices, the
				       // first entry in each row
				       // is the diagonal one. In
				       // this case only sort the
				       // remaining entries, otherwise
				       // sort all
      sort ((rows==cols) ? &tmp_entries[1] : &tmp_entries[0],
	    &tmp_entries[row_length]);

				       // Re-insert column numbers
				       // into the field
      for (unsigned int j=0; j<row_length; ++j)
	colnums[next_free_entry++] = tmp_entries[j];

				       // note new start of this and
				       // the next row
      rowstart[line] = next_row_start;
      next_row_start = next_free_entry;

				       // some internal checks
      Assert ((rows!=cols) ||
	      (colnums[rowstart[line]] == static_cast<signed int>(line)),
	      ExcInternalError());
				       // assert that the first entry
				       // does not show up in
				       // the remaining ones and that
				       // the remaining ones are unique
				       // among themselves (this handles
				       // both cases, quadratic and
				       // rectangular matrices)
      Assert (find (&colnums[rowstart[line]+1],
		    &colnums[next_row_start],
		    colnums[rowstart[line]]) ==
	      &colnums[next_row_start],
	      ExcInternalError());
      Assert (adjacent_find(&colnums[rowstart[line]+1],
			    &colnums[next_row_start]) ==
	      &colnums[next_row_start],
	      ExcInternalError());
    };
  
  vec_len = rowstart[rows] = next_row_start;
  compressed = true;

  delete[] tmp_entries;
};



int
dSMatrixStruct::operator () (const unsigned int i, const unsigned int j) const
{
  Assert (i<rows, ExcInvalidIndex(i,rows));
  Assert (j<cols, ExcInvalidIndex(j,cols));
  Assert (compressed, ExcNotCompressed());

				   // check first entry separately, since
				   // for square matrices this is
				   // the diagonal entry (check only
				   // if a first entry exists)
  if (rowstart[i] != rowstart[i+1]) 
    {
      if (static_cast<signed int>(j) == colnums[rowstart[i]])
	return rowstart[i];
    }
  else
				     // no first entry exists for this
				     // line
    return -1;

				   // all other entries are sorted, so
				   // we can use a binary seach algorithm
  const int* p = lower_bound (&colnums[rowstart[i]+1],
			      &colnums[rowstart[i+1]],
			      static_cast<signed int>(j));
  if (*p == static_cast<signed int>(j))
    return (p - &colnums[0]);
  else
    return -1;
}


void
dSMatrixStruct::add (const unsigned int i, const unsigned int j)
{
  Assert (i<rows, ExcInvalidIndex(i,rows));
  Assert (j<cols, ExcInvalidIndex(j,cols));
  Assert (compressed==false, ExcMatrixIsCompressed());

  for (unsigned int k=rowstart[i]; k<rowstart[i+1]; k++)
    {
				       // entry already exists
      if (colnums[k] == (signed int)j) return;
				       // empty entry found, put new
				       // entry here
      if (colnums[k] == -1)
	{
	  colnums[k] = j;
	  return;
	};
    };

				   // if we came thus far, something went
				   // wrong: there was not enough space
				   // in this line
  Assert (false, ExcNotEnoughSpace(i, rowstart[i+1]-rowstart[i]));
}



void
dSMatrixStruct::add_matrix (const unsigned int n, const int* rowcols)
{
  for (unsigned int i=0; i<n; ++i)
    for (unsigned int j=0; j<n; ++j)
      add(rowcols[i], rowcols[j]);
}



void
dSMatrixStruct::add_matrix (const unsigned int m, const unsigned int n,
			    const int* rows, const int* cols)
{
  for (unsigned i=0; i<m; ++i)
    for (unsigned j=0; j<n; ++j)
      add(rows[i], cols[j]);
}



void
dSMatrixStruct::add_matrix (const iVector& rowcols)
{
  unsigned int i,j;
  for (i=0;i<rowcols.n();i++)
    for (j=0;j<rowcols.n();j++)
      add(rowcols(i), rowcols(j));
}



void
dSMatrixStruct::add_matrix (const iVector& rows, const iVector& cols)
{
  unsigned int i,j;
  for (i=0;i<rows.n();i++)
    for (j=0;j<cols.n();j++)
      add(rows(i), cols(j));
}



void
dSMatrixStruct::print_gnuplot (ostream &out) const
{
  for (unsigned int i=0; i<rows; ++i)
    for (unsigned int j=rowstart[i]; j<rowstart[i+1]; ++j)
      if (colnums[j]>=0)
	out << i << " " << -colnums[j] << endl;
}



unsigned int
dSMatrixStruct::bandwidth () const
{
  unsigned int b=0;
  for (unsigned int i=0; i<rows; ++i)
    for (unsigned int j=rowstart[i]; j<rowstart[i+1]; ++j)
      if (colnums[j]>=0) 
	{
	  if ((unsigned int)abs(i-colnums[j]) > b)
	    b = abs(i-colnums[j]);
	}
      else
					 // leave if at the end of
					 // the entries of this line
	break;
  return b;
}





/*-------------------------------------------------------------------------*/


dSMatrix::dSMatrix () :
		cols(0),
		val(0),
		max_len(0) {};



dSMatrix::dSMatrix (dSMatrixStruct &c)
		: cols(&c), val(0), max_len(0)
{
  reinit();
};



dSMatrix::~dSMatrix ()
{
  delete[] val;
};



void
dSMatrix::reinit ()
{
  Assert (cols != 0, ExcMatrixNotInitialized());
  Assert (cols->compressed, ExcNotCompressed());
  
  if (max_len<cols->vec_len)
  {
    if (val) delete[] val;
    val = new double[cols->vec_len];
    max_len = cols->vec_len;
  }
//  memset(val, 0, sizeof(*val) * cols->vec_len);
  for (int i = cols->vec_len-1 ; i>=0 ; i--) val[i] = 0;
}



void
dSMatrix::reinit (dSMatrixStruct &sparsity) {
  cols = &sparsity;
  reinit ();
};



void
dSMatrix::clear () {
  cols = 0;
  if (val) delete[] val;
  val = 0;
  max_len = 0;
};



dSMatrix &
dSMatrix::copy_from (const dSMatrix &matrix) {
  Assert (cols != 0, ExcMatrixNotInitialized());
  Assert (cols == matrix.cols, ExcDifferentSparsityPatterns());

  for (unsigned int i = 0 ; i<cols->vec_len; ++i)
    val[i] = matrix.val[i];

  return *this;
};



void
dSMatrix::add_scaled (const double factor, const dSMatrix &matrix) {
  Assert (cols != 0, ExcMatrixNotInitialized());
  Assert (cols == matrix.cols, ExcDifferentSparsityPatterns());

  for (unsigned int i = 0 ; i<cols->vec_len; ++i)
    val[i] += factor*matrix.val[i];
};



void
dSMatrix::vmult (dVector& dst, const dVector& src) const
{
  Assert (cols != 0, ExcMatrixNotInitialized());
  Assert(m() == dst.n(), ExcDimensionsDontMatch(m(),dst.n()));
  Assert(n() == src.n(), ExcDimensionsDontMatch(n(),src.n()));

  for (unsigned int i=0; i<m(); ++i)
    {
      double s = 0.;
      for (unsigned int j=cols->rowstart[i]; j<cols->rowstart[i+1]; ++j) 
	s += val[j] * src(cols->colnums[j]);
      dst(i) = s;
    }
}


void
dSMatrix::Tvmult (dVector& dst, const dVector& src) const
{
  Assert (cols != 0, ExcMatrixNotInitialized());
  Assert(n() == dst.n(), ExcDimensionsDontMatch(n(),dst.n()));
  Assert(m() == src.n(), ExcDimensionsDontMatch(m(),src.n()));

  unsigned int i;
  
  for (i=0;i<n();i++) dst(i) = 0.;
  
  for (i=0;i<m();i++)
    {
      for (unsigned int j=cols->rowstart[i]; j<cols->rowstart[i+1] ;j++)
	{
	  int p = cols->colnums[j];
	  dst(p) += val[j] * src(i);
	}
    }
}



double
dSMatrix::matrix_norm (const dVector& v) const
{
  Assert (cols != 0, ExcMatrixNotInitialized());
  Assert(m() == v.n(), ExcDimensionsDontMatch(m(),v.n()));
  Assert(n() == v.n(), ExcDimensionsDontMatch(n(),v.n()));

  double sum = 0.;
  for (unsigned int i=0;i<m();i++)
    {
      double s = 0.;
      for (unsigned int j=cols->rowstart[i]; j<cols->rowstart[i+1]; ++j) 
	s += val[j] * v(cols->colnums[j]);
      sum += s*v(i);
    };

  return sum;
}



double
dSMatrix::residual (dVector& dst, const dVector& u, const dVector& b)
{
  Assert (cols != 0, ExcMatrixNotInitialized());
  Assert(m() == dst.n(), ExcDimensionsDontMatch(m(),dst.n()));
  Assert(m() == b.n(), ExcDimensionsDontMatch(m(),b.n()));
  Assert(n() == u.n(), ExcDimensionsDontMatch(n(),u.n()));

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

void
dSMatrix::Jacobi_precond (dVector& dst, const dVector& src, const double om)
{
  Assert (cols != 0, ExcMatrixNotInitialized());

  const unsigned int n = src.n();

  for (unsigned int i=0;i<n;++i)
    {
      dst(i) = om * src(i) * val[cols->rowstart[i]];
    }
}

void
dSMatrix::SSOR_precond (dVector& dst, const dVector& src, const double om)
{
  Assert (cols != 0, ExcMatrixNotInitialized());

  int p;
  unsigned int  n = src.n();
  unsigned int  j;
  
  for (unsigned i=0; i<n; i++)
    {
      dst(i) = src(i);
      for (j=cols->rowstart[i]; j<cols->rowstart[i+1] ;j++)
	{
	  p = cols->colnums[j];
	  if (p<(signed int)i)
	    dst(i) -= om* val[j] * dst(p);
	}
      dst(i) /= val[cols->rowstart[i]];
    }
  for (unsigned int i=0; i<n; i++)
    dst(i) *= (2.-om)*val[cols->rowstart[i]];
  
  for (int i=n-1; i>=0; i--)
    {
      for (j=cols->rowstart[i];j<cols->rowstart[i+1];j++)
	{
	  p = cols->colnums[j];
	  if (p>i) dst(i) -= om* val[j] * dst(p);
	}
      dst(i) /= val[cols->rowstart[i]];
    }
}

void
dSMatrix::SOR_precond (dVector& dst, const dVector& src, const double om)
{
  Assert (cols != 0, ExcMatrixNotInitialized());
  dst = src;
  SOR(dst,om);
}

void
dSMatrix::SOR (dVector& dst, const double om)
{
  Assert (cols != 0, ExcMatrixNotInitialized());
  Assert(n() == m(), ExcDimensionsDontMatch(n(),m()));
  Assert(m() == dst.n(), ExcDimensionsDontMatch(m(),dst.n()));

  for (unsigned int i=0;i<m();i++)
    {
      double s = dst(i);
      for (unsigned int j=cols->rowstart[i]; j<cols->rowstart[i+1] ;j++)
	{
	  int p = cols->colnums[j];
	  if (p<(signed int)i)
	    s -= val[j] * dst(p);
	}
      dst(i) = s * om / val[cols->rowstart[i]];
    }
}

void
dSMatrix::SSOR (dVector& dst, const double om)
{
  Assert (cols != 0, ExcMatrixNotInitialized());

  int p;
  unsigned int  n = dst.n();
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



const dSMatrixStruct & dSMatrix::get_sparsity_pattern () const {
  Assert (cols != 0, ExcMatrixNotInitialized());
  return *cols;
};



void dSMatrix::print (ostream &out) const {
  Assert (cols != 0, ExcMatrixNotInitialized());

  for (unsigned int i=0; i<cols->rows; ++i)
    for (unsigned int j=cols->rowstart[i]; j<cols->rowstart[i+1]; ++j)
      out << "(" << i << "," << cols->colnums[j] << ") " << val[j] << endl;
};



void dSMatrix::print_formatted (ostream &out, const unsigned int precision) const {
  Assert (cols != 0, ExcMatrixNotInitialized());
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

  out.setf (0, ios::floatfield);                 // reset output format
};
