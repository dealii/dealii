// $Id$

#include <lac/vector.h>
#include <lac/full_matrix.h>

#include <vector>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <iomanip>
#include <algorithm>


template <typename number>
FullMatrix<number>::FullMatrix (const unsigned int n)
{
  init (n,n);
};


template <typename number>
FullMatrix<number>::FullMatrix (const unsigned int m, const unsigned int n)
{
  init (m,n);
};


template <typename number>
FullMatrix<number>::FullMatrix (const FullMatrix &m):
		Subscriptor()
{
  init (m.dim_image, m.dim_range);
  if (dim_range*dim_image != 0)
    copy (&m.val[0], &m.val[dim_image*dim_range],
	  &val[0]);
};


template <typename number>
void
FullMatrix<number>::init (const unsigned int mm, const unsigned int nn)
{
  val_size = nn*mm;
  val = new number[val_size];
  dim_range = nn;
  dim_image = mm;
  clear ();
};


template <typename number>
FullMatrix<number>::~FullMatrix ()
{
  delete[] val;
};


template <typename number>
bool
FullMatrix<number>::all_zero () const
{
  const number *p = &val[0],
	       *e = &val[n()*m()];
  while (p!=e)
    if (*p++ != 0.0)
      return false;

  return true;
};


template <typename number>
void
FullMatrix<number>::reinit (const unsigned int mm,
			    const unsigned int nn)
{
  if (val_size<nn*mm)
    {
      delete[] val;
      init(mm, nn);
    }
  else
    {
      dim_range = nn;
      dim_image = mm;
      clear ();
    }
}


template <typename number>
void
FullMatrix<number>::reinit (const unsigned int n)
{
  reinit (n, n);
};


template <typename number>
template <typename number2>
void
FullMatrix<number>::reinit (const FullMatrix<number2> &B)
{
  reinit (B.m(), B.n());
};



template <typename number>
template <typename number2>
void
FullMatrix<number>::vmult (Vector<number2>& dst,
			   const Vector<number2>& src,
			   const bool adding) const
{
  Assert(dst.size() == m(), ExcDimensionMismatch(dst.size(), m()));
  Assert(src.size() == n(), ExcDimensionMismatch(src.size(), n()));

  number2 s;
  if ((n()==3) && (m()==3))
  {
    number2 s0,s1,s2;
    s   = src(0);
    s0  = s*val[0]; s1  = s*val[3]; s2  = s*val[6]; 
    s   = src(1);
    s0 += s*val[1]; s1 += s*val[4]; s2 += s*val[7];
    s   = src(2);
    s0 += s*val[2]; s1 += s*val[5]; s2 += s*val[8];

    if (!adding)
    {
      dst(0) = s0;
      dst(1) = s1;
      dst(2) = s2;
    }
    else
    {
      dst(0) += s0;
      dst(1) += s1;
      dst(2) += s2;
    }
  }
  else if ((n()==4) && (m()==4))
  {
    number2 s0,s1,s2,s3;
    s = src(0);
    s0  = s*val[0]; s1  = s*val[4]; s2  = s*val[8];  s3  = s*val[12];
    s = src(1);
    s0 += s*val[1]; s1 += s*val[5]; s2 += s*val[9];  s3 += s*val[13];
    s = src(2);
    s0 += s*val[2]; s1 += s*val[6]; s2 += s*val[10]; s3 += s*val[14];
    s = src(3);
    s0 += s*val[3]; s1 += s*val[7]; s2 += s*val[11]; s3 += s*val[15];
    
    if (!adding)
    {
      dst(0) = s0;
      dst(1) = s1;
      dst(2) = s2;
      dst(3) = s3;
    }
    else
    {
      dst(0) += s0;
      dst(1) += s1;
      dst(2) += s2;
      dst(3) += s3;
    }
  }
  else if ((n()==8) && (m()==8))
  {
    number2 s0,s1,s2,s3,s4,s5,s6,s7;
    s = src(0);
    s0 = s*val[0]; s1 = s*val[8]; s2 = s*val[16]; s3 = s*val[24];
    s4 = s*val[32]; s5 = s*val[40]; s6 = s*val[48]; s7 = s*val[56];
    s = src(1);
    s0 += s*val[1]; s1 += s*val[9]; s2 += s*val[17]; s3 += s*val[25];
    s4 += s*val[33]; s5 += s*val[41]; s6 += s*val[49]; s7 += s*val[57];
    s = src(2);
    s0 += s*val[2]; s1 += s*val[10]; s2 += s*val[18]; s3 += s*val[26];
    s4 += s*val[34]; s5 += s*val[42]; s6 += s*val[50]; s7 += s*val[58];
    s = src(3);
    s0 += s*val[3]; s1 += s*val[11]; s2 += s*val[19]; s3 += s*val[27];
    s4 += s*val[35]; s5 += s*val[43]; s6 += s*val[51]; s7 += s*val[59];
    s = src(4);
    s0 += s*val[4]; s1 += s*val[12]; s2 += s*val[20]; s3 += s*val[28];
    s4 += s*val[36]; s5 += s*val[44]; s6 += s*val[52]; s7 += s*val[60];
    s = src(5);
    s0 += s*val[5]; s1 += s*val[13]; s2 += s*val[21]; s3 += s*val[29];
    s4 += s*val[37]; s5 += s*val[45]; s6 += s*val[53]; s7 += s*val[61];
    s = src(6);
    s0 += s*val[6]; s1 += s*val[14]; s2 += s*val[22]; s3 += s*val[30];
    s4 += s*val[38]; s5 += s*val[46]; s6 += s*val[54]; s7 += s*val[62];
    s = src(7);
    s0 += s*val[7]; s1 += s*val[15]; s2 += s*val[23]; s3 += s*val[31];
    s4 += s*val[39]; s5 += s*val[47]; s6 += s*val[55]; s7 += s*val[63];
    
    if (!adding)
    {
      dst(0) = s0;
      dst(1) = s1;
      dst(2) = s2;
      dst(3) = s3;
      dst(4) = s4;
      dst(5) = s5;
      dst(6) = s6;
      dst(7) = s7;
    }
    else
    {
      dst(0) += s0;
      dst(1) += s1;
      dst(2) += s2;
      dst(3) += s3;
      dst(4) += s4;
      dst(5) += s5;
      dst(6) += s6;
      dst(7) += s7;
    }
  }
  else
  {    
    number* e = val;
    const unsigned int size_m = m(),
		       size_n = n();
    for (unsigned int i=0; i<size_m; ++i)
    {
      s = 0.;
      for (unsigned int j=0; j<size_n; ++j)
	s += src(j) * *(e++);
      if (!adding) dst(i) = s;
      else dst(i) += s;
    }
  }
}



template <typename number>
template <typename number2>
void FullMatrix<number>::Tvmult (Vector<number2>       &dst,
				 const Vector<number2> &src,
				 const bool             adding) const
{
  Assert(dst.size() == n(), ExcDimensionMismatch(dst.size(), n()));
  Assert(src.size() == m(), ExcDimensionMismatch(src.size(), m()));

  unsigned int i,j;
  number2 s;
  const unsigned int size_m = m(),
		     size_n = n();
  for (i=0; i<size_m; ++i)
  {
    s = 0.;
    for (j=0; j<size_n; ++j)
      s += src(j) * el(j,i);
    if(!adding) dst(i) = s;
    else dst(i) += s;
  }
}



template <typename number>
template <typename number2, typename number3>
double FullMatrix<number>::residual (Vector<number2>& dst,
				     const Vector<number2>& src,
				     const Vector<number3>& right) const
{
  Assert(dst.size() == m(), ExcDimensionMismatch(dst.size(), m()));
  Assert(src.size() == n(), ExcDimensionMismatch(src.size(), n()));
  Assert(right.size() == m(), ExcDimensionMismatch(right.size(), m()));

  number2 s, res = 0.;
  const unsigned int size_m = m(),
		     size_n = n();
  for (unsigned int i=0; i<size_n; ++i)
    {
      s = right(i);
      for (unsigned int j=0; j<size_m; ++j)
	s -= src(j) * el(i,j);
      dst(i) = s;
      res += s*s;
    }
  return sqrt(res);
}



template <typename number>
template <typename number2>
void FullMatrix<number>::forward (Vector<number2>& dst,
				  const Vector<number2>& src) const
{
  Assert(dst.size() == m(), ExcDimensionMismatch(dst.size(), m()));
  Assert(src.size() == n(), ExcDimensionMismatch(src.size(), n()));

  unsigned int i,j;
  unsigned int nu = ( (m()<n()) ? m() : n());
  number2 s;
  for (i=0; i<nu; ++i)
    {
      s = src(i);
      for (j=0; j<i; ++j) s -= dst(j) * el(i,j);
      dst(i) = s/el(i,i);
    }
}



template <typename number>
template <typename number2>
void FullMatrix<number>::backward (Vector<number2>& dst,
				   const Vector<number2>& src) const
{
  unsigned int j;
  unsigned int nu = (m()<n() ? m() : n());
  number2 s;
  for (int i=nu-1; i>=0; --i)
    {
      s = src(i);
      for (j=i+1; j<nu; ++j) s -= dst(j) * el(i,j);
      dst(i) = s/el(i,i);
    }
}



template <typename number>
FullMatrix<number>&
FullMatrix<number>::operator = (const FullMatrix<number>& m) 
{
  reinit(m);
  if (dim_range*dim_image != 0)
    copy (&m.val[0], &m.val[dim_image*dim_range],
	  &val[0]);
  
  return *this;
}



template <typename number>
template <typename number2>
FullMatrix<number>&
FullMatrix<number>::operator = (const FullMatrix<number2>& m) 
{
  reinit(m);
  if (dim_range*dim_image != 0)
    copy (&m.val[0], &m.val[dim_image*dim_range],
	  &val[0]);

  return *this;
}



template <typename number>
template <typename number2>
void FullMatrix<number>::fill (const FullMatrix<number2>& src,
			       const unsigned int i,
			       const unsigned int j)
{
  Assert (n() >= src.n() + j, ExcInvalidDestination(n(), src.n(), j));
  Assert (m() >= src.m() + i, ExcInvalidDestination(m(), src.m(), i));

  for (unsigned int ii=0; ii<src.m() ; ++ii)
    for (unsigned int jj=0; jj<src.n() ; ++jj)
      el(ii+i,jj+j) = src.el(ii,jj);
}



template <typename number>
void FullMatrix<number>::add_row (const unsigned int i,
				  const number s,
				  const unsigned int j)
{
  for (unsigned int k=0; k<m(); ++k)
    el(i,k) += s*el(j,k);
}



template <typename number>
void FullMatrix<number>::add_row (const unsigned int i,
				  const number s,
				  const unsigned int j,
				  const number t,
				  const unsigned int k)
{
  const unsigned int size_m = m();
  for (unsigned l=0; l<size_m; ++l)
    el(i,l) += s*el(j,l) + t*el(k,l);
}



template <typename number>
void FullMatrix<number>::add_col (const unsigned int i, const number s,
			const unsigned int j)
{
  for (unsigned int k=0; k<n(); ++k)
    el(k,i) += s*el(k,j);
}



template <typename number>
void FullMatrix<number>::add_col (const unsigned int i, const number s,
		        const unsigned int j, const number t,
			const unsigned int k)
{
  for (unsigned int l=0; l<n(); ++l)
    el(l,i) += s*el(l,j) + t*el(l,k);
}



template <typename number>
void FullMatrix<number>::swap_row (const unsigned int i, const unsigned int j)
{
  number s;
  for (unsigned int k=0; k<m(); ++k)
  {
    s = el(i,k); el(i,k) = el(j,k); el(j,k) = s;
  }
}



template <typename number>
void FullMatrix<number>::swap_col (const unsigned int i, const unsigned int j)
{
  number s;
  for (unsigned int k=0; k<n(); ++k)
  {
    s = el(k,i); el(k,i) = el(k,j); el(k,j) = s;
  }
}



template <typename number>
void FullMatrix<number>::diagadd (const number src)
{
  Assert (m() == n(), ExcDimensionMismatch(m(),n()));
  for (unsigned int i=0; i<n(); ++i)
    el(i,i) += src;
}



template <typename number>
template <typename number2>
void FullMatrix<number>::mmult (FullMatrix<number2>& dst,
				const FullMatrix<number2>& src) const
{
  Assert (n() == src.m(), ExcDimensionMismatch(n(), src.m()));
  unsigned int i,j,k;
  number2 s = 1.;
  dst.reinit(m(), src.n());

  for (i=0;i<m();i++)
    for (j=0; j<src.n(); ++j)
      {
	s = 0.;
	for (k=0;k<n();k++) s+= el(i,k) * src.el(k,j);
	dst.el(i,j) = s;
      }
}



/*void FullMatrix<number>::mmult (FullMatrix& dst, const FullMatrix& src) const
{
  Assert (m() == src.n(), ExcDimensionMismatch(m(), src.n()));

  unsigned int i,j,k;
  number2 s = 1.;

  dst.reinit(n(), src.m());

  for (i=0;i<n();i++)
    for (j=0;j<src.m();j++)
      {
	s = 0.;
	for (k=0;k<m();k++) s+= el(i,k) * src.el(k,j);
	dst.el(i,j) = s;
      }
}*/



template <typename number>
template <typename number2>
void FullMatrix<number>::Tmmult (FullMatrix<number2>& dst, const FullMatrix<number2>& src) const
{
  Assert (n() == src.n(), ExcDimensionMismatch(n(), src.n()));

  unsigned int i,j,k;
  number2 s = 1.;
  dst.reinit(m(), src.m());

  for (i=0;i<m();i++)
    for (j=0;j<src.m();j++)
      {
	s = 0.;
	for (k=0;k<n();k++) s+= el(k,i) * src.el(k,j);
	dst.el(i,j) = s;
      }
}



/*void FullMatrix<number>::Tmmult(FullMatrix& dst, const FullMatrix& src) const
{
  Assert (m() == src.n(), ExcDimensionMismatch(m(), src.n()));

  unsigned int i,j,k;
  number2 s = 1.;
  
  dst.reinit(n(), src.m());

  for (i=0;i<n();i++)
    for (j=0;j<src.m();j++)
      {
	s = 0.;
	for (k=0;k<m();k++) s+= el(k,i) * src.el(k,j);
	dst.el(i,j) = s;
      }
}*/



template <typename number>
template <typename number2>
double FullMatrix<number>::matrix_norm (const Vector<number2> &v) const
{
  Assert(m() == v.size(), ExcDimensionMismatch(m(),v.size()));
  Assert(n() == v.size(), ExcDimensionMismatch(n(),v.size()));

  number2 sum = 0.;
  const unsigned int n_rows = m();
  const number *val_ptr = &val[0];
  const number2 *v_ptr;
  
  for (unsigned int row=0; row<n_rows; ++row)
    {
      number2 s = 0.;
      const number * const val_end_of_row = val_ptr+n_rows;
      v_ptr = v.begin();
      while (val_ptr != val_end_of_row)
	s += *val_ptr++ * *v_ptr++;

      sum += s* v(row);
    };

  return sum;
};



template <typename number>
template <typename number2>
double FullMatrix<number>::matrix_scalar_product (const Vector<number2> &u, const Vector<number2> &v) const
{
  Assert(m() == u.size(), ExcDimensionMismatch(m(),v.size()));
  Assert(n() == v.size(), ExcDimensionMismatch(n(),v.size()));

  number2 sum = 0.;
  const unsigned int n_rows = m();
  const unsigned int n_cols = n();
  const number *val_ptr = &val[0];
  const number2 *v_ptr;
  
  for (unsigned int row=0; row<n_rows; ++row)
    {
      number2 s = 0.;
      const number * const val_end_of_row = val_ptr+n_cols;
      v_ptr = v.begin();
      while (val_ptr != val_end_of_row)
	s += *val_ptr++ * *v_ptr++;

      sum += s* u(row);
    };

  return sum;
};



template <typename number>
number FullMatrix<number>::l1_norm () const
{
  number sum=0, max=0;
  const unsigned int n_rows = m(), n_cols = n();
  
  for (unsigned int col=0; col<n_cols; ++col)
    {
      sum=0;
      for (unsigned int row=0; row<n_rows; ++row)
	sum += fabs(el(row,col));
      if (sum > max)
	max = sum;
    }
  return max;
};



template <typename number>
number FullMatrix<number>::linfty_norm () const
{
  number sum=0, max=0;
  const unsigned int n_rows = m(), n_cols = n();

  for (unsigned int row=0; row<n_rows; ++row)
    {
      sum=0;
      for (unsigned int col=0; col<n_cols; ++col)
	sum += fabs(el(row,col));
      if (sum > max)
	max = sum;
    }
  return max;
};



template <typename number>
void
FullMatrix<number>::print (ostream& s, int w, int p) const
{
  unsigned int i,j;
  for (i=0;i<m();i++)
    {
      for (j=0;j<n();j++) s << setw(w) << setprecision(p) << el(i,j);
      s << endl;
    }
}



template <typename number>
template <typename number2>
void
FullMatrix<number>::add (const number s,const FullMatrix<number2>& src)
{
  Assert (m() == src.m(), ExcDimensionMismatch(m(), src.m()));
  Assert (n() == src.n(), ExcDimensionMismatch(n(), src.n()));
  if ((n()==3) && (m()==3))
  {
    val[0] += s * src.val[0];
    val[1] += s * src.val[1];
    val[2] += s * src.val[2];
    val[3] += s * src.val[3];
    val[4] += s * src.val[4];
    val[5] += s * src.val[5];
    val[6] += s * src.val[6];
    val[7] += s * src.val[7];
    val[8] += s * src.val[8];
  }
  else if ((n()==4) && (m()==4))
  {
    val[0] += s * src.val[0];
    val[1] += s * src.val[1];
    val[2] += s * src.val[2];
    val[3] += s * src.val[3];
    val[4] += s * src.val[4];
    val[5] += s * src.val[5];
    val[6] += s * src.val[6];
    val[7] += s * src.val[7];
    val[8] += s * src.val[8];
    val[9] += s * src.val[9];
    val[10] += s * src.val[10];
    val[11] += s * src.val[11];
    val[12] += s * src.val[12];
    val[13] += s * src.val[13];
    val[14] += s * src.val[14];
    val[15] += s * src.val[15];
  }
  else if ((n()==8) && (m()==8))
  {
    val[0] += s * src.val[0];
    val[1] += s * src.val[1];
    val[2] += s * src.val[2];
    val[3] += s * src.val[3];
    val[4] += s * src.val[4];
    val[5] += s * src.val[5];
    val[6] += s * src.val[6];
    val[7] += s * src.val[7];
    val[8] += s * src.val[8];
    val[9] += s * src.val[9];
    val[10] += s * src.val[10];
    val[11] += s * src.val[11];
    val[12] += s * src.val[12];
    val[13] += s * src.val[13];
    val[14] += s * src.val[14];
    val[15] += s * src.val[15];
    val[16] += s * src.val[16];
    val[17] += s * src.val[17];
    val[18] += s * src.val[18];
    val[19] += s * src.val[19];

    val[20] += s * src.val[20];
    val[21] += s * src.val[21];
    val[22] += s * src.val[22];
    val[23] += s * src.val[23];
    val[24] += s * src.val[24];
    val[25] += s * src.val[25];
    val[26] += s * src.val[26];
    val[27] += s * src.val[27];
    val[28] += s * src.val[28];
    val[29] += s * src.val[29];

    val[30] += s * src.val[30];
    val[31] += s * src.val[31];
    val[32] += s * src.val[32];
    val[33] += s * src.val[33];
    val[34] += s * src.val[34];
    val[35] += s * src.val[35];
    val[36] += s * src.val[36];
    val[37] += s * src.val[37];
    val[38] += s * src.val[38];
    val[39] += s * src.val[39];

    val[40] += s * src.val[40];
    val[41] += s * src.val[41];
    val[42] += s * src.val[42];
    val[43] += s * src.val[43];
    val[44] += s * src.val[44];
    val[45] += s * src.val[45];
    val[46] += s * src.val[46];
    val[47] += s * src.val[47];
    val[48] += s * src.val[48];
    val[49] += s * src.val[49];

    val[50] += s * src.val[50];
    val[51] += s * src.val[51];
    val[52] += s * src.val[52];
    val[53] += s * src.val[53];
    val[54] += s * src.val[54];
    val[55] += s * src.val[55];
    val[56] += s * src.val[56];
    val[57] += s * src.val[57];
    val[58] += s * src.val[58];
    val[59] += s * src.val[59];

    val[60] += s * src.val[60];
    val[61] += s * src.val[61];
    val[62] += s * src.val[62];
    val[63] += s * src.val[63];
  }
  else
  {
    const unsigned int size = n()*m();
    for (unsigned int i=0; i<size; i++)
      val[i] += s * src.val[i];
  }
}



template <typename number>
template <typename number2>
void
FullMatrix<number>::add_diag (const number s, const FullMatrix<number2>& src)
{
  Assert (m() == src.m(), ExcDimensionMismatch(m(), src.m()));
  Assert (n() == src.n(), ExcDimensionMismatch(n(), src.n()));

  if ((n()==3) && (m()==3))
  {
    val[0] += s * src.val[0];
    val[0] += s * src.val[1];
    val[0] += s * src.val[2];
    val[3] += s * src.val[3];
    val[3] += s * src.val[4];
    val[3] += s * src.val[5];
    val[6] += s * src.val[6];
    val[6] += s * src.val[7];
    val[6] += s * src.val[8];
  }
  else if ((n()==4) && (m()==4))
  {
    val[0] += s * src.val[0];
    val[0] += s * src.val[1];
    val[0] += s * src.val[2];
    val[0] += s * src.val[3];
    val[4] += s * src.val[4];
    val[4] += s * src.val[5];
    val[4] += s * src.val[6];
    val[4] += s * src.val[7];
    val[8] += s * src.val[8];
    val[8] += s * src.val[9];
    val[8] += s * src.val[10];
    val[8] += s * src.val[11];
    val[12] += s * src.val[12];
    val[12] += s * src.val[13];
    val[12] += s * src.val[14];
    val[12] += s * src.val[15];
  }
  else if ((n()==8) && (m()==8))
  {
    val[0] += s * src.val[0];
    val[0] += s * src.val[1];
    val[0] += s * src.val[2];
    val[0] += s * src.val[3];
    val[0] += s * src.val[4];
    val[0] += s * src.val[5];
    val[0] += s * src.val[6];
    val[0] += s * src.val[7];
    val[8] += s * src.val[8];
    val[8] += s * src.val[9];
    val[8] += s * src.val[10];
    val[8] += s * src.val[11];
    val[8] += s * src.val[12];
    val[8] += s * src.val[13];
    val[8] += s * src.val[14];
    val[8] += s * src.val[15];
    val[16] += s * src.val[16];
    val[16] += s * src.val[17];
    val[16] += s * src.val[18];
    val[16] += s * src.val[19];

    val[16] += s * src.val[20];
    val[16] += s * src.val[21];
    val[16] += s * src.val[22];
    val[16] += s * src.val[23];
    val[24] += s * src.val[24];
    val[24] += s * src.val[25];
    val[24] += s * src.val[26];
    val[24] += s * src.val[27];
    val[24] += s * src.val[28];
    val[24] += s * src.val[29];

    val[24] += s * src.val[30];
    val[24] += s * src.val[31];
    val[32] += s * src.val[32];
    val[32] += s * src.val[33];
    val[32] += s * src.val[34];
    val[32] += s * src.val[35];
    val[32] += s * src.val[36];
    val[32] += s * src.val[37];
    val[32] += s * src.val[38];
    val[32] += s * src.val[39];

    val[40] += s * src.val[40];
    val[40] += s * src.val[41];
    val[40] += s * src.val[42];
    val[40] += s * src.val[43];
    val[40] += s * src.val[44];
    val[40] += s * src.val[45];
    val[40] += s * src.val[46];
    val[40] += s * src.val[47];
    val[48] += s * src.val[48];
    val[48] += s * src.val[49];

    val[48] += s * src.val[50];
    val[48] += s * src.val[51];
    val[48] += s * src.val[52];
    val[48] += s * src.val[53];
    val[48] += s * src.val[54];
    val[48] += s * src.val[55];
    val[56] += s * src.val[56];
    val[56] += s * src.val[57];
    val[56] += s * src.val[58];
    val[56] += s * src.val[59];

    val[56] += s * src.val[60];
    val[56] += s * src.val[61];
    val[56] += s * src.val[62];
    val[56] += s * src.val[63];
  }
  else
  {
    const unsigned int size = n()*m();
    for (unsigned int i=0; i<size; i++)
      val[i] += s * src.val[i];
  }
}



template <typename number>
template <typename number2>
void
FullMatrix<number>::Tadd (const number s, const FullMatrix<number2>& src)
{
  Assert (m() == n(),     ExcNotQuadratic());
  Assert (m() == src.m(), ExcDimensionMismatch(m(), src.m()));
  Assert (n() == src.n(), ExcDimensionMismatch(n(), src.n()));

  if ((n()==3) && (m()==3))
  {
    val[0] += s * src.val[0];
    val[1] += s * src.val[3];
    val[2] += s * src.val[6];

    val[3] += s * src.val[1];
    val[4] += s * src.val[4];
    val[5] += s * src.val[7];

    val[6] += s * src.val[2];
    val[7] += s * src.val[5];
    val[8] += s * src.val[8];
  }
  else if ((n()==4) && (m()==4))
  {
    val[0] += s * src.val[0];
    val[1] += s * src.val[4];
    val[2] += s * src.val[8];
    val[3] += s * src.val[12];

    val[4] += s * src.val[1];
    val[5] += s * src.val[5];
    val[6] += s * src.val[9];
    val[7] += s * src.val[13];

    val[8] += s * src.val[2];
    val[9] += s * src.val[6];
    val[10] += s * src.val[10];
    val[11] += s * src.val[14];

    val[12] += s * src.val[3];
    val[13] += s * src.val[7];
    val[14] += s * src.val[11];
    val[15] += s * src.val[15];
  }
  else if ((n()==8) && (m()==8))
  {
    val[0] += s * src.val[0];
    val[1] += s * src.val[8];
    val[2] += s * src.val[16];
    val[3] += s * src.val[24];
    val[4] += s * src.val[32];
    val[5] += s * src.val[40];
    val[6] += s * src.val[48];
    val[7] += s * src.val[56];

    val[8] += s * src.val[1];
    val[9] += s * src.val[9];
    val[10] += s * src.val[17];
    val[11] += s * src.val[25];
    val[12] += s * src.val[33];
    val[13] += s * src.val[41];
    val[14] += s * src.val[49];
    val[15] += s * src.val[57];

    val[16] += s * src.val[2];
    val[17] += s * src.val[10];
    val[18] += s * src.val[18];
    val[19] += s * src.val[26];
    val[20] += s * src.val[34];
    val[21] += s * src.val[42];
    val[22] += s * src.val[50];
    val[23] += s * src.val[58];

    val[24] += s * src.val[3];
    val[25] += s * src.val[11];
    val[26] += s * src.val[19];
    val[27] += s * src.val[27];
    val[28] += s * src.val[35];
    val[29] += s * src.val[43];
    val[30] += s * src.val[51];
    val[31] += s * src.val[59];

    val[32] += s * src.val[4];
    val[33] += s * src.val[12];
    val[34] += s * src.val[20];
    val[35] += s * src.val[28];
    val[36] += s * src.val[36];
    val[37] += s * src.val[44];
    val[38] += s * src.val[52];
    val[39] += s * src.val[60];

    val[40] += s * src.val[5];
    val[41] += s * src.val[13];
    val[42] += s * src.val[21];
    val[43] += s * src.val[29];
    val[44] += s * src.val[37];
    val[45] += s * src.val[45];
    val[46] += s * src.val[53];
    val[47] += s * src.val[61];

    val[48] += s * src.val[6];
    val[49] += s * src.val[14];
    val[50] += s * src.val[22];
    val[51] += s * src.val[30];
    val[52] += s * src.val[38];
    val[53] += s * src.val[46];
    val[54] += s * src.val[54];
    val[55] += s * src.val[62];

    val[56] += s * src.val[7];
    val[57] += s * src.val[15];
    val[58] += s * src.val[23];
    val[59] += s * src.val[31];
    val[60] += s * src.val[39];
    val[61] += s * src.val[47];
    val[62] += s * src.val[55];
    val[63] += s * src.val[63];
  }
  else
    Assert (false, ExcNotImplemented(n()));
}



template <typename number>
bool
FullMatrix<number>::operator == (const FullMatrix<number> &m) const
{
  bool q = (dim_range==m.dim_range) && (dim_image==m.dim_image);
  if (!q) return false;

  return equal (&val[0], &val[dim_range*dim_image],
		&m.val[0]);
};



template <typename number>
double
FullMatrix<number>::determinant () const
{
  Assert (dim_range == dim_image,
	  ExcDimensionMismatch(dim_range, dim_image));
  Assert ((dim_range>=1) && (dim_range<=3), ExcNotImplemented(dim_range));
  
  switch (dim_range) 
    {
      case 1:
            return el(0,0);
      case 2:
            return el(0,0)*el(1,1) - el(1,0)*el(0,1);
      case 3:
            return  (el(0,0)*el(1,1)*el(2,2)
                     -el(0,0)*el(1,2)*el(2,1)
                     -el(1,0)*el(0,1)*el(2,2)
                     +el(1,0)*el(0,2)*el(2,1)
                     +el(2,0)*el(0,1)*el(1,2)
                     -el(2,0)*el(0,2)*el(1,1));
      default:
            return 0;
    };
};



template <typename number>
number
FullMatrix<number>::norm2 () const
{
  number s = 0.;
  for (unsigned int i=0;i<dim_image*dim_range;++i)
    s += val[i]*val[i];
  return s;
}



template <typename number>
void FullMatrix<number>::clear ()
{
  fill_n (&val[0], n()*m(), 0);
};



template <typename number>
void
FullMatrix<number>::invert (const FullMatrix<number> &M)
{
  Assert (dim_range == dim_image, ExcNotQuadratic());
  Assert (dim_range == M.dim_range,
          ExcDimensionMismatch(dim_range,M.dim_range));
  Assert (dim_image == M.dim_image,
	  ExcDimensionMismatch(dim_image,M.dim_image));
  
  switch (dim_range) 
    {
      case 1:
	    val[0] = 1.0/M.val[0];
	    return;
      case 2:
					     // this is Maple output,
					     // thus a bit unstructured
      {
	    const number t4 = 1.0/(M.el(0,0)*M.el(1,1)-M.el(0,1)*M.el(1,0));
	    el(0,0) = M.el(1,1)*t4;
	    el(0,1) = -M.el(0,1)*t4;
	    el(1,0) = -M.el(1,0)*t4;
	    el(1,1) = M.el(0,0)*t4;
	    return;
      };
      
      case 3:
      {
	    const number t4 = M.el(0,0)*M.el(1,1),
			 t6 = M.el(0,0)*M.el(1,2),
			 t8 = M.el(0,1)*M.el(1,0),
			t00 = M.el(0,2)*M.el(1,0),
			t01 = M.el(0,1)*M.el(2,0),
			t04 = M.el(0,2)*M.el(2,0),
			t07 = 1.0/(t4*M.el(2,2)-t6*M.el(2,1)-t8*M.el(2,2)+
 				   t00*M.el(2,1)+t01*M.el(1,2)-t04*M.el(1,1));
	    el(0,0) = (M.el(1,1)*M.el(2,2)-M.el(1,2)*M.el(2,1))*t07;
	    el(0,1) = -(M.el(0,1)*M.el(2,2)-M.el(0,2)*M.el(2,1))*t07;
	    el(0,2) = -(-M.el(0,1)*M.el(1,2)+M.el(0,2)*M.el(1,1))*t07;
	    el(1,0) = -(M.el(1,0)*M.el(2,2)-M.el(1,2)*M.el(2,0))*t07;
	    el(1,1) = (M.el(0,0)*M.el(2,2)-t04)*t07;
	    el(1,2) = -(t6-t00)*t07;
	    el(2,0) = -(-M.el(1,0)*M.el(2,1)+M.el(1,1)*M.el(2,0))*t07;
	    el(2,1) = -(M.el(0,0)*M.el(2,1)-t01)*t07;
	    el(2,2) = (t4-t8)*t07;
	    return;
      };

      case 4:
      {
					 // with (linalg);
					 // a:=matrix(4,4);
					 // evalm(a);
					 // ai:=inverse(a);
					 // readlib(C);
					 // C(ai,optimized,filename=x4);

	const number t14 = M.el(0,0)*M.el(1,1);
	const number t15 = M.el(2,2)*M.el(3,3);
	const number t17 = M.el(2,3)*M.el(3,2);
	const number t19 = M.el(0,0)*M.el(2,1);
	const number t20 = M.el(1,2)*M.el(3,3);
	const number t22 = M.el(1,3)*M.el(3,2);
	const number t24 = M.el(0,0)*M.el(3,1);
	const number t25 = M.el(1,2)*M.el(2,3);
	const number t27 = M.el(1,3)*M.el(2,2);
	const number t29 = M.el(1,0)*M.el(0,1);
	const number t32 = M.el(1,0)*M.el(2,1);
	const number t33 = M.el(0,2)*M.el(3,3);
	const number t35 = M.el(0,3)*M.el(3,2);
	const number t37 = M.el(1,0)*M.el(3,1);
	const number t38 = M.el(0,2)*M.el(2,3);
	const number t40 = M.el(0,3)*M.el(2,2);
	const number t42 = t14*t15-t14*t17-t19*t20+t19*t22+
			   t24*t25-t24*t27-t29*t15+t29*t17+
			   t32*t33-t32*t35-t37*t38+t37*t40;
	const number t43 = M.el(2,0)*M.el(0,1);
	const number t46 = M.el(2,0)*M.el(1,1);
	const number t49 = M.el(2,0)*M.el(3,1);
	const number t50 = M.el(0,2)*M.el(1,3);
	const number t52 = M.el(0,3)*M.el(1,2);
	const number t54 = M.el(3,0)*M.el(0,1);
	const number t57 = M.el(3,0)*M.el(1,1);
	const number t60 = M.el(3,0)*M.el(2,1);
	const number t63 = t43*t20-t43*t22-t46*t33+t46*t35+
			   t49*t50-t49*t52-t54*t25+t54*t27+
			   t57*t38-t57*t40-t60*t50+t60*t52;
	const number t65 = 1/(t42+t63);
	const number t71 = M.el(0,2)*M.el(2,1);
	const number t73 = M.el(0,3)*M.el(2,1);
	const number t75 = M.el(0,2)*M.el(3,1);
	const number t77 = M.el(0,3)*M.el(3,1);
	const number t81 = M.el(0,1)*M.el(1,2);
	const number t83 = M.el(0,1)*M.el(1,3);
	const number t85 = M.el(0,2)*M.el(1,1);
	const number t87 = M.el(0,3)*M.el(1,1);
	const number t101 = M.el(1,0)*M.el(2,2);
	const number t103 = M.el(1,0)*M.el(2,3);
	const number t105 = M.el(2,0)*M.el(1,2);
	const number t107 = M.el(2,0)*M.el(1,3);
	const number t109 = M.el(3,0)*M.el(1,2);
	const number t111 = M.el(3,0)*M.el(1,3);
	const number t115 = M.el(0,0)*M.el(2,2);
	const number t117 = M.el(0,0)*M.el(2,3);
	const number t119 = M.el(2,0)*M.el(0,2);
	const number t121 = M.el(2,0)*M.el(0,3);
	const number t123 = M.el(3,0)*M.el(0,2);
	const number t125 = M.el(3,0)*M.el(0,3);
	const number t129 = M.el(0,0)*M.el(1,2);
	const number t131 = M.el(0,0)*M.el(1,3);
	const number t133 = M.el(1,0)*M.el(0,2);
	const number t135 = M.el(1,0)*M.el(0,3);
	el(0,0) = (M.el(1,1)*M.el(2,2)*M.el(3,3)-M.el(1,1)*M.el(2,3)*M.el(3,2)-
		   M.el(2,1)*M.el(1,2)*M.el(3,3)+M.el(2,1)*M.el(1,3)*M.el(3,2)+
		   M.el(3,1)*M.el(1,2)*M.el(2,3)-M.el(3,1)*M.el(1,3)*M.el(2,2))*t65;
	el(0,1) = -(M.el(0,1)*M.el(2,2)*M.el(3,3)-M.el(0,1)*M.el(2,3)*M.el(3,2)-
		    t71*M.el(3,3)+t73*M.el(3,2)+t75*M.el(2,3)-t77*M.el(2,2))*t65;
	el(0,2) = (t81*M.el(3,3)-t83*M.el(3,2)-t85*M.el(3,3)+t87*M.el(3,2)+
		   t75*M.el(1,3)-t77*M.el(1,2))*t65;
	el(0,3) = -(t81*M.el(2,3)-t83*M.el(2,2)-t85*M.el(2,3)+t87*M.el(2,2)+
		    t71*M.el(1,3)-t73*M.el(1,2))*t65;
	el(1,0) = -(t101*M.el(3,3)-t103*M.el(3,2)-t105*M.el(3,3)+t107*M.el(3,2)+
		    t109*M.el(2,3)-t111*M.el(2,2))*t65;
	el(1,1) = (t115*M.el(3,3)-t117*M.el(3,2)-t119*M.el(3,3)+t121*M.el(3,2)+
		   t123*M.el(2,3)-t125*M.el(2,2))*t65;
	el(1,2) = -(t129*M.el(3,3)-t131*M.el(3,2)-t133*M.el(3,3)+t135*M.el(3,2)+
		    t123*M.el(1,3)-t125*M.el(1,2))*t65;
	el(1,3) = (t129*M.el(2,3)-t131*M.el(2,2)-t133*M.el(2,3)+t135*M.el(2,2)+
		   t119*M.el(1,3)-t121*M.el(1,2))*t65;
	el(2,0) = (t32*M.el(3,3)-t103*M.el(3,1)-t46*M.el(3,3)+t107*M.el(3,1)+
		   t57*M.el(2,3)-t111*M.el(2,1))*t65;
	el(2,1) = -(t19*M.el(3,3)-t117*M.el(3,1)-t43*M.el(3,3)+t121*M.el(3,1)+
		    t54*M.el(2,3)-t125*M.el(2,1))*t65;
	el(2,2) = (t14*M.el(3,3)-t131*M.el(3,1)-t29*M.el(3,3)+t135*M.el(3,1)+
		   t54*M.el(1,3)-t125*M.el(1,1))*t65;
	el(2,3) = -(t14*M.el(2,3)-t131*M.el(2,1)-t29*M.el(2,3)+t135*M.el(2,1)+
		    t43*M.el(1,3)-t121*M.el(1,1))*t65;
	el(3,0) = -(t32*M.el(3,2)-t101*M.el(3,1)-t46*M.el(3,2)+t105*M.el(3,1)+
		    t57*M.el(2,2)-t109*M.el(2,1))*t65;
	el(3,1) = (t19*M.el(3,2)-t115*M.el(3,1)-t43*M.el(3,2)+t119*M.el(3,1)+
		   t54*M.el(2,2)-t123*M.el(2,1))*t65;
	el(3,2) = -(t14*M.el(3,2)-t129*M.el(3,1)-t29*M.el(3,2)+t133*M.el(3,1)+
		    t54*M.el(1,2)-t123*M.el(1,1))*t65;
	el(3,3) = (t14*M.el(2,2)-t129*M.el(2,1)-t29*M.el(2,2)+t133*M.el(2,1)+
		   t43*M.el(1,2)-t119*M.el(1,1))*t65;

	break;
      }

      default:
	    AssertThrow (false, ExcNotImplemented(dim_range));
    };    
};
  


template <typename number>
void
FullMatrix<number>::print_formatted (ostream &out, const unsigned int precision) const
{
  out.precision (precision);
  out.setf (ios::scientific, ios::floatfield);   // set output format
  
  for (unsigned int i=0; i<m(); ++i) 
    {
      for (unsigned int j=0; j<n(); ++j)
	if (el(i,j) != 0)
	  out << setw(precision+7)
	      << el(i,j) << ' ';
	else
	  out << setw(precision+8) << " ";
      out << endl;
    };

  AssertThrow (out, ExcIO());
  
  out.setf (0, ios::floatfield);                 // reset output format
};


// Gauss-Jordan-Algorithmus
// cf. Stoer I (4th Edition) p. 153

template <typename number>
void
FullMatrix<number>::gauss_jordan()
{
  Assert (dim_range == dim_image, ExcNotQuadratic());
  vector<unsigned int> p(n());

  unsigned int i,j,k,r;
  number max, hr;

  for (i=0; i<n(); ++i)
    p[i] = i;

  for (j=0; j<n(); ++j)
    {
				       // pivotsearch
      max = fabs(el(j,j));
      r = j;
      for (i=j+1; i<n(); ++i)
	{
	  if (fabs(el(i,j)) > max)
	    {
	      max = fabs(el(i,j));
	      r = i;
	    }
	}
      Assert(max>1.e-16, ExcNotRegular());
				       // rowinterchange
      if (r>j)
	{
	  for (k=0; k<n(); ++k)
	    {
	      hr = el(j,k) ; el(j,k) = el(r,k) ; el(r,k) = hr;
	    }
	  i = p[j];
	  p[j] = p[r];
	  p[r] = i;
	}

				       // transformation
      hr = 1./el(j,j);
      el(j,j) = hr;
      for (k=0; k<n(); ++k)
	{
	  if (k==j) continue;
	  for (i=0; i<n(); ++i)
	    {
	      if (i==j) continue;
	      el(i,k) -= el(i,j)*el(j,k)*hr;
	    }
	}
      for (i=0; i<n(); ++i)
	{
	  el(i,j) *= hr;
	  el(j,i) *= -hr;
	}
      el(j,j) = hr;
    }
				   // columninterchange
  Vector<number> hv(n());
  for (i=0; i<n(); ++i)
    {
      for (k=0; k<n(); ++k) hv(p[k]) = el(i,k);
      for (k=0; k<n(); ++k) el(i,k) = hv(k);
    }
}

// QR-transformation cf. Stoer 1 4.8.2 (p. 191)

template <typename number>
template <typename number2>
void
FullMatrix<number>::householder(Vector<number2>& src)
{
  // m > n, src.n() = m
  Assert (dim_range <= dim_image, ExcDimensionMismatch(dim_range, dim_image));
  Assert (src.size() == dim_image, ExcDimensionMismatch(src.size(), dim_image));

  for (unsigned int j=0 ; j<n() ; ++j)
  {
    number2 sigma = 0;
    unsigned int i;
    for (i=j ; i<m() ; ++i) sigma += el(i,j)*el(i,j);
    if (fabs(sigma) < 1.e-15) return;
    number2 s = el(j,j);
    s = (s<0) ? sqrt(sigma) : -sqrt(sigma);
    number2 dj = s;

    number2 beta = 1./(s*el(j,j)-sigma);
    el(j,j) -= s;

    for (unsigned int k=j+1 ; k<n() ; ++k)
    {
      number2 sum = 0.;
      for (i=j ; i<m() ; ++i) sum += el(i,j)*el(i,k);
      sum *= beta;

      for (i=j ; i<m() ; ++i) el(i,k) += sum*el(i,j);
    }

    number2 sum = 0.;
    for (i=j ; i<m() ; ++i) sum += el(i,j)*src(i);
    sum *= beta;

    for (i=j ; i<m() ; ++i) src(i) += sum*el(i,j);
    el(j,j) = dj;
  }
}



template <typename number>
template <typename number2>
double
FullMatrix<number>::least_squares(Vector<number2>& dst, Vector<number2>& src)
{
  // m > n, m = src.n, n = dst.n

  householder(src);
  backward(dst, src);

  number2 sum = 0.;
  for (unsigned int i=n() ; i<m() ; ++i) sum += src(i) * src(i);
  return sqrt(sum);
}
