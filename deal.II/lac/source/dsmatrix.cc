// $Id$

// This file is part of the DEAL Library
// DEAL is Copyright(1995) by
// Roland Becker, Guido Kanschat, Franz-Theo Suttmeier

#include <lac/dsmatrix.h>
#include <ostream.h>


/*----------------- from sort.h -------------------------*/


//////////
template<class T>
inline void swap(T* a, T* b)
{
  T x = *a;
  *a = *b;
  *b = x;
}

//////////
template<class T>
inline void simple_sort(long n, T* field)
{
  long i,j;
  for (i=1;i<n;i++)
  {
    for (j=i+1;j<=n;j++)
    {
      if (field[j] < field[i])
      {
	swap(field+i,field+j);
      }
    }
  }
}

template<class T>
inline void heapsort_sift(T* a, long l, long r)
{
  long i = l;
  long j = 2*i;
  T    x = a[i];

  while (j<=r)
  {
    if (j<r) 
    {
      if (a[j] < a[j+1]) j++;
    }
    if (!(x < a[j])) break;
    a[i] = a[j];
    i = j;
    j = 2*i;
  }
  a[i] = x;
}


//////////
template<class T>
inline void heapsort(int n, T* field)
{
  field--;

  long l =(n/2)+1;
  long r = n;

  while (l>1)
  {
    l--;
    heapsort_sift(field,l,r);
  }
  while (r>1)
  {
    swap(field+l,field+r);
    r--;
    heapsort_sift(field,l,r);
  }
}

//////////
template<class T>
inline void _quicksort(long r, T* a, long l)
{
  long i = l;
  long j = r;
  T*   x = &a[(l+r)/2];
  do
  {
    while (a[i] < *x) i++;
    while (*x < a[j]) j--;
    if (i<=j)
    {
      swap(a+i,a+j);
      i++;
      j--;
    }
  }
  while (i<=j);
  if (l<j) _quicksort(j,a,l);
  if (i<r) _quicksort(r,a,i);
}

template<class T>
inline void quicksort(long r, T* a)
{
  _quicksort(r,a,1);
}


/*----------------- from sort.h -------------------------*/



void
dSMatrixStruct::reinit(int m, int n, int max_per_row)
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
      rowstart = new int[max_dim+1];
    }

  if (vec_len > max_vec_len)
    {
      if (colnums) delete[] colnums;
      max_vec_len = vec_len;
      colnums = new int[max_vec_len];
    }
  int i;
  for (i=0;i<=rows;i++) rowstart[i] = i * max_per_row;
  for (i = vec_len-1 ; i>=0 ; i--) colnums[i] = -1;

  if (rows == cols)
    for (i=0;i<rows;i++)
      colnums[rowstart[i]] = i;

  compressed = 0;
}

void
dSMatrixStruct::compress()
{
  if (compressed) return;
  int i,j,k,l,s;
  int* entries;

  entries = new int[max_row_len];

  // Traverse all rows

  for (i=0,k=0,s=0 ; i < rows ; i++)
    {

      // Sort entries in ascending order

      for (j=rowstart[i],l=0; j < rowstart[i+1] ;++j,++l)
	entries[l] = colnums[j];
      heapsort(max_row_len, entries);

      // Re-insert column numbers into the field

      // Ensure diagonal entry first in quadratic matrix

      if (cols == rows) colnums[k++] = i;

      for (l=0; l < max_row_len ; ++l)
	{
	  if (entries[l] >= 0)
	    if ((i!=entries[l]) || (rows!=cols))
	      colnums[k++] = entries[l];
	}
      rowstart[i] = s;
      s = k;
    }
  vec_len = rowstart[i] = s;
  compressed = 1;

  delete[] entries;
}

int
dSMatrixStruct::operator () (int i, int j)
{
  for (int k=rowstart[i] ; k<rowstart[i+1] ; k++)
    if (colnums[k] == j) return k;
  return -1;
}

void
dSMatrixStruct::add (int i, int j)
{
  Assert ((i>=0) && (i<rows), ExcInvalidIndex(i,rows));
  Assert ((j>=0) && (j<cols), ExcInvalidIndex(j,cols));

  int k;
  for (k=rowstart[i]; k<rowstart[i+1]; k++)
    {
      if (colnums[k] == j) return;
      if (colnums[k] == -1)
	{
	  colnums[k] = j;
	  return;
	}
    }
}

void
dSMatrixStruct::add_matrix(int n, int* rowcols)
{
  int i,j;
  for (i=0;i<n;i++)
    for (j=0;j<n;j++)
      add(rowcols[i], rowcols[j]);
}


void
dSMatrixStruct::add_matrix(int m, int n, int* rows, int* cols)
{
  int i,j;
  for (i=0;i<m;i++)
    for (j=0;j<n;j++)
      add(rows[i], cols[j]);
}

void
dSMatrixStruct::add_matrix(const iVector& rowcols)
{
  int i,j;
  for (i=0;i<rowcols.n();i++)
    for (j=0;j<rowcols.n();j++)
      add(rowcols(i), rowcols(j));
}


void
dSMatrixStruct::add_matrix(const iVector& rows, const iVector& cols)
{
  int i,j;
  for (i=0;i<rows.n();i++)
    for (j=0;j<cols.n();j++)
      add(rows(i), cols(j));
}

void
dSMatrix::reinit()
{
  if(max_len<cols.vec_len)
  {
    if (val) delete[] val;
    val = new double[cols.vec_len];
    max_len = cols.vec_len;
  }
//  memset(val, 0, sizeof(*val) * cols.vec_len);
  for (int i = cols.vec_len-1 ; i>=0 ; i--) val[i] = 0;
}

void
dSMatrixStruct::print_gnuplot (ostream &out) const
{
  for (int i=0; i<rows; ++i)
    for (int j=rowstart[i]; j<rowstart[i+1]; ++j)
      if (colnums[j]>=0)
	out << i << " " << -colnums[j] << endl;
}


int
dSMatrixStruct::bandwidth () const
{
  int b=0;
  for (int i=0; i<rows; ++i)
    for (int j=rowstart[i]; j<rowstart[i+1]; ++j)
      if (colnums[j]>=0) 
	{
	  if (abs(i-colnums[j]) > b)
	    b = abs(i-colnums[j]);
	}
      else
					 // leave if at the end of
					 // the entries of this line
	break;
  return b;
}

void
dSMatrix::vmult(dVector& dst,const dVector& src)
{
  // Assert(m() = dst.n())
  // Assert(n() = src.n())

  for (int i=0;i<m();i++)
    {
      double s = 0.;
      for (int j=cols.rowstart[i]; j < cols.rowstart[i+1] ;j++) 
	{
	  int p = cols.colnums[j];
	  s += val[j] * src(p);
	}
      dst(i) = s;
    }
}

void
dSMatrix::Tvmult(dVector& dst,const dVector& src) 
{
  // Assert(n() = dst.n())
  // Assert(m() = src.n())

  int i;
  
  for (i=0;i<n();i++) dst(i) = 0.;
  
  for (i=0;i<m();i++)
    {
      for (int j=cols.rowstart[i]; j<cols.rowstart[i+1] ;j++)
	{
	  int p = cols.colnums[j];
	  dst(p) += val[j] * src(i);
	}
    }
}

double
dSMatrix::residual(dVector& dst,const dVector& u,const dVector& b)
{
  // Assert(m() = dst.n())
  // Assert(n() = src.n())

  double s,norm=0.;   
  
  for (int i=0;i<m();i++)
    {
      s = b(i);
      for (int j=cols.rowstart[i]; j<cols.rowstart[i+1] ;j++)
	{
	  int p = cols.colnums[j];
	  s -= val[j] * u(p);
	}
      dst(i) = s;
      norm += dst(i)*dst(i);
    }
  return sqrt(norm);
}

void
dSMatrix::Jacobi_precond(dVector& dst,const dVector& src,double om)
{
  int n = src.n();

  for (int i=0;i<n;++i)
    {
      dst(i) = om * src(i) * val[cols.rowstart[i]];
    }
}

void
dSMatrix::SSOR_precond(dVector& dst,const dVector& src,double om)
{
  int  p,n = src.n();
  int  i,j;
  
  for (i=0;i<n;i++)
    {
      dst(i) = src(i);
      for (j=cols.rowstart[i]; j<cols.rowstart[i+1] ;j++)
	{
	  p = cols.colnums[j];
	  if (p<i) dst(i) -= om* val[j] * dst(p);
	}
      dst(i) /= val[cols.rowstart[i]];
    }
  for (i=0;i<n;i++) dst(i) *= (2.-om)*val[cols.rowstart[i]];
  
  for (i=n-1;i>=0;i--)
    {
      for (j=cols.rowstart[i];j<cols.rowstart[i+1];j++)
	{
	  p = cols.colnums[j];
	  if (p>i) dst(i) -= om* val[j] * dst(p);
	}
      dst(i) /= val[cols.rowstart[i]];
    }
}

void
dSMatrix::SOR_precond(dVector& dst,const dVector& src,double om)
{
  dst = src;
  SOR(dst,om);
}

void
dSMatrix::SOR(dVector& dst, double om)
{
  // Assert(m()==n())
  // Assert(m()==dst.n())

  for (int i=0;i<m();i++)
    {
      double s = dst(i);
      for (int j=cols.rowstart[i]; j<cols.rowstart[i+1] ;j++)
	{
	  int p = cols.colnums[j];
	  if (p<i)
	    s -= val[j] * dst(p);
	}
      dst(i) = s * om / val[cols.rowstart[i]];
    }
}

void
dSMatrix::SSOR(dVector& dst, double om)
{
  int  p,n = dst.n();
  int  i,j;
  double s;
  
  for (i=0;i<n;i++)
    {
      s = 0.;
      for (j=cols.rowstart[i]; j<cols.rowstart[i+1] ;j++)
	{
	  p = cols.colnums[j];
	  if (p>=0)
	    {
	      if (i>j) s += val[j] * dst(p);
	    }
	}
      dst(i) -= s * om;
      dst(i) /= val[cols.rowstart[i]];
    }

  for (i=n-1;i>=0;i--)
    {
      s = 0.;
      for (j=cols.rowstart[i]; j<cols.rowstart[i+1] ;j++)
	{
	  p = cols.colnums[j];
	  if (p>=0)
	    {
	      if (i<j) s += val[j] * dst(p);
	    }
	}
      dst(i) -= s * om / val[cols.rowstart[i]];
    }
}


void dSMatrix::print (ostream &out) const {
  for (int i=0; i<cols.rows; ++i)
    for (int j=cols.rowstart[i]; j<cols.rowstart[i+1]; ++j)
      out << "(" << i << "," << cols.colnums[j] << ") " << val[j] << endl;
};
