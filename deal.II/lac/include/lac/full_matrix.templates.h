// $Id$

#include <lac/vector.h>
#include <lac/ivector.h>
#include <lac/fullmatrix.h>

#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <iomanip>


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
FullMatrix<number>::FullMatrix (const FullMatrix &m) 
{
  init (m.dim_image, m.dim_range);
  number       *       p = &val[0];
  const number *      vp = &m.val[0];
  const number * const e = &val[dim_image*dim_range];

  while (p!=e)
    *p++ = *vp++;
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
FullMatrix<number>::reinit (const unsigned int mm, const unsigned int nn)
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

  double s;
  if ((n()==3) && (m()==3))
  {
    double s0,s1,s2;
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
    double s0,s1,s2,s3;
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
    double s0,s1,s2,s3,s4,s5,s6,s7;
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
void FullMatrix<number>::gsmult (Vector<number2>& dst, const Vector<number2>& src, const iVector& gl) const
{
  Assert(n() == m(), ExcNotQuadratic());
  Assert(dst.size() == n(), ExcDimensionMismatch(dst.size(), n()));
  Assert(src.size() == n(), ExcDimensionMismatch(src.size(), n()));
  Assert(gl.n() == n(), ExcDimensionMismatch(gl.n(), n()));

  double s;
  if ((n()==3) && (m()==3))
  {
    double s0=0.,s1=0.,s2=0.;
    s   = src(0);
    if(gl(1)<gl(0)) s1  = s*val[3]; if(gl(2)<gl(0))  s2  = s*val[6]; 
    s   = src(1);
    if(gl(0)<gl(1)) s0 += s*val[1]; if(gl(2)<gl(1))  s2 += s*val[7];
    s   = src(2);
    if(gl(0)<gl(2)) s0 += s*val[2]; if(gl(1)<gl(2))  s1 += s*val[5]; 

    dst(0) += s0;
    dst(1) += s1;
    dst(2) += s2;
  }
  else if ((n()==4) && (m()==4))
  {
    double s0=0.,s1=0.,s2=0.,s3=0.;
    s = src(0);
    if(gl(1)<gl(0)) s1 = s*val[4];  if(gl(2)<gl(0)) s2 = s*val[8]; if(gl(3)<gl(0)) s3 = s*val[12];
    s = src(1);
    if(gl(0)<gl(1)) s0 += s*val[1]; if(gl(2)<gl(1)) s2 += s*val[9]; if(gl(3)<gl(1)) s3 += s*val[13];
    s = src(2);
    if(gl(0)<gl(2)) s0 += s*val[2]; if(gl(1)<gl(2)) s1 += s*val[6]; if(gl(3)<gl(2)) s3 += s*val[14];
    s = src(3);
    if(gl(0)<gl(3)) s0 += s*val[3]; if(gl(1)<gl(3)) s1 += s*val[7]; if(gl(2)<gl(3)) s2 += s*val[11];

    dst(0) += s0;
    dst(1) += s1;
    dst(2) += s2;
    dst(3) += s3;
  }
  else if ((n()==8) && (m()==8))
  {
    double s0=0.,s1=0.,s2=0.,s3=0.,s4=0.,s5=0.,s6=0.,s7=0.;
    s = src(0);
    if(gl(1)<gl(0)) s1 = s*val[8]; 
    if(gl(2)<gl(0)) s2 = s*val[16]; 
    if(gl(3)<gl(0)) s3 = s*val[24];
    if(gl(4)<gl(0)) s4 = s*val[32]; 
    if(gl(5)<gl(0)) s5 = s*val[40]; 
    if(gl(6)<gl(0)) s6 = s*val[48]; 
    if(gl(7)<gl(0)) s7 = s*val[56];
    s = src(1);
    if(gl(0)<gl(1)) s0 += s*val[1]; 
    if(gl(2)<gl(1)) s2 += s*val[17]; 
    if(gl(3)<gl(1)) s3 += s*val[25];
    if(gl(4)<gl(1)) s4 += s*val[33]; 
    if(gl(5)<gl(1)) s5 += s*val[41]; 
    if(gl(6)<gl(1)) s6 += s*val[49]; 
    if(gl(7)<gl(1)) s7 += s*val[57];
    s = src(2);
    if(gl(0)<gl(2)) s0 += s*val[2]; 
    if(gl(1)<gl(2)) s1 += s*val[10]; 
    if(gl(3)<gl(2)) s3 += s*val[26];
    if(gl(4)<gl(2)) s4 += s*val[34]; 
    if(gl(5)<gl(2)) s5 += s*val[42]; 
    if(gl(6)<gl(2)) s6 += s*val[50]; 
    if(gl(7)<gl(2)) s7 += s*val[58];
    s = src(3);
    if(gl(0)<gl(3)) s0 += s*val[3]; 
    if(gl(1)<gl(3)) s1 += s*val[11]; 
    if(gl(2)<gl(3)) s2 += s*val[19]; 
    if(gl(4)<gl(3)) s4 += s*val[35]; 
    if(gl(5)<gl(3)) s5 += s*val[43]; 
    if(gl(6)<gl(3)) s6 += s*val[51]; 
    if(gl(7)<gl(3)) s7 += s*val[59];
    s = src(4);
    if(gl(0)<gl(4)) s0 += s*val[4]; 
    if(gl(1)<gl(4)) s1 += s*val[12]; 
    if(gl(2)<gl(4)) s2 += s*val[20]; 
    if(gl(3)<gl(4)) s3 += s*val[28];
    if(gl(5)<gl(4)) s5 += s*val[44]; 
    if(gl(6)<gl(4)) s6 += s*val[52]; 
    if(gl(7)<gl(4)) s7 += s*val[60];
    s = src(5);
    if(gl(0)<gl(5)) s0 += s*val[5]; 
    if(gl(1)<gl(5)) s1 += s*val[13]; 
    if(gl(2)<gl(5)) s2 += s*val[21]; 
    if(gl(3)<gl(5)) s3 += s*val[29];
    if(gl(4)<gl(5)) s4 += s*val[37]; 
    if(gl(6)<gl(5)) s6 += s*val[53]; 
    if(gl(7)<gl(5)) s7 += s*val[61];
    s = src(6);
    if(gl(0)<gl(6)) s0 += s*val[6]; 
    if(gl(1)<gl(6)) s1 += s*val[14]; 
    if(gl(2)<gl(6)) s2 += s*val[22]; 
    if(gl(3)<gl(6)) s3 += s*val[30];
    if(gl(4)<gl(6)) s4 += s*val[38]; 
    if(gl(5)<gl(6)) s5 += s*val[46]; 
    if(gl(7)<gl(6)) s7 += s*val[62];
    s = src(7);
    if(gl(0)<gl(7)) s0 += s*val[7]; 
    if(gl(1)<gl(7)) s1 += s*val[15]; 
    if(gl(2)<gl(7)) s2 += s*val[23]; 
    if(gl(3)<gl(7)) s3 += s*val[31];
    if(gl(4)<gl(7)) s4 += s*val[39]; 
    if(gl(5)<gl(7)) s5 += s*val[47]; 
    if(gl(6)<gl(7)) s6 += s*val[55]; 
    
    dst(0) += s0;
    dst(1) += s1;
    dst(2) += s2;
    dst(3) += s3;
    dst(4) += s4;
    dst(5) += s5;
    dst(6) += s6;
    dst(7) += s7;
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
	if(gl(i)<gl(j)) s += src(j) * *(e++);
      dst(i) += s;
    }
  }
}

template <typename number>
template <typename number2>
void FullMatrix<number>::Tvmult (Vector<number2>& dst, const Vector<number2>& src, const bool adding) const
{
  Assert(dst.size() == n(), ExcDimensionMismatch(dst.size(), n()));
  Assert(src.size() == m(), ExcDimensionMismatch(src.size(), m()));

  unsigned int i,j;
  double s;
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
double FullMatrix<number>::residual (Vector<number2>& dst, const Vector<number2>& src,
			   const Vector<number3>& right) const
{
  Assert(dst.size() == m(), ExcDimensionMismatch(dst.size(), m()));
  Assert(src.size() == n(), ExcDimensionMismatch(src.size(), n()));
  Assert(right.size() == m(), ExcDimensionMismatch(right.size(), m()));

  unsigned int i,j;
  double s, res = 0.;
  const unsigned int size_m = m(),
		     size_n = n();
  for (i=0; i<size_n; ++i)
    {
      s = right(i);
      for (j=0; j<size_m; ++j)
	s -= src(j) * el(i,j);
      dst(i) = s;
      res += s*s;
    }
  return sqrt(res);
}

template <typename number>
template <typename number2>
void FullMatrix<number>::forward (Vector<number2>& dst, const Vector<number2>& src) const
{
  Assert(n() == m(), ExcNotQuadratic());
  Assert(dst.size() == n(), ExcDimensionMismatch(dst.size(), n()));
  Assert(src.size() == n(), ExcDimensionMismatch(src.size(), n()));

  unsigned int i,j;
  unsigned int nu = (m()<n() ? m() : n());
  double s;
  for (i=0; i<nu; ++i)
    {
      s = src(i);
      for (j=0; j<i; ++j) s -= dst(j) * el(i,j);
      dst(i) = s/el(i,i);
    }
}

template <typename number>
template <typename number2>
void FullMatrix<number>::backward (Vector<number2>& dst, const Vector<number2>& src) const
{
  Assert(n() == m(), ExcNotQuadratic());
  Assert(dst.size() == n(), ExcDimensionMismatch(dst.size(), n()));
  Assert(src.size() == n(), ExcDimensionMismatch(src.size(), n()));

  unsigned int j;
  unsigned int nu = (m()<n() ? m() : n());
  double s;
  for (int i=nu-1; i>=0; --i)
    {
      s = src(i);
      for (j=i+1; j<nu; ++j) s -= dst(j) * el(i,j);
      dst(i) = s/el(i,i);
    }
}

template <typename number>
template <typename number2>
FullMatrix<number>&
FullMatrix<number>::operator = (const FullMatrix<number2>& m) 
{
  reinit(m);

  number *             p = &val[0];
  const number2 *     vp = &m.val[0];
  const number * const e = &val[dim_image*dim_range];

  while (p!=e)
    *p++ = *vp++;

  return *this;
}

template <typename number>
template <typename number2>
void FullMatrix<number>::fill (const FullMatrix<number2>& src,
		     const unsigned int i, const unsigned int j)
{
  Assert (n() >= src.n() + j, ExcInvalidDestination(n(), src.n(), j));
  Assert (m() >= src.m() + i, ExcInvalidDestination(m(), src.m(), i));

  for (unsigned int ii=0; ii<src.m() ; ++ii)
    for (unsigned int jj=0; jj<src.n() ; ++jj)
      el(ii+i,jj+j) = src.el(ii,jj);
}

template <typename number>
void FullMatrix<number>::add_row (const unsigned int i,
			const number s, const unsigned int j)
{
  for (unsigned int k=0; k<m(); ++k)
    el(i,k) += s*el(j,k);
}

template <typename number>
void FullMatrix<number>::add_row (const unsigned int i, const number s,
		        const unsigned int j, const number t,
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
void FullMatrix<number>::mmult (FullMatrix<number2>& dst, const FullMatrix<number2>& src) const
{
  Assert (n() == src.m(), ExcDimensionMismatch(n(), src.m()));
  unsigned int i,j,k;
  double s = 1.;
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
  double s = 1.;

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
  double s = 1.;
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
  double s = 1.;
  
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

  double sum = 0.;
  const unsigned int n_rows = m();
  const number *val_ptr = &val[0];
  const number2 *v_ptr;
  
  for (unsigned int row=0; row<n_rows; ++row)
    {
      double s = 0.;
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

  double sum = 0.;
  const unsigned int n_rows = m();
  const unsigned int n_cols = n();
  const number *val_ptr = &val[0];
  const number2 *v_ptr;
  
  for (unsigned int row=0; row<n_rows; ++row)
    {
      double s = 0.;
      const number * const val_end_of_row = val_ptr+n_cols;
      v_ptr = v.begin();
      while (val_ptr != val_end_of_row)
	s += *val_ptr++ * *v_ptr++;

      sum += s* u(row);
    };

  return sum;
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
    val[0] += s * src.el(0);
    val[1] += s * src.el(1);
    val[2] += s * src.el(2);
    val[3] += s * src.el(3);
    val[4] += s * src.el(4);
    val[5] += s * src.el(5);
    val[6] += s * src.el(6);
    val[7] += s * src.el(7);
    val[8] += s * src.el(8);
  }
  else if ((n()==4) && (m()==4))
  {
    val[0] += s * src.el(0);
    val[1] += s * src.el(1);
    val[2] += s * src.el(2);
    val[3] += s * src.el(3);
    val[4] += s * src.el(4);
    val[5] += s * src.el(5);
    val[6] += s * src.el(6);
    val[7] += s * src.el(7);
    val[8] += s * src.el(8);
    val[9] += s * src.el(9);
    val[10] += s * src.el(10);
    val[11] += s * src.el(11);
    val[12] += s * src.el(12);
    val[13] += s * src.el(13);
    val[14] += s * src.el(14);
    val[15] += s * src.el(15);
  }
  else if ((n()==8) && (m()==8))
  {
    val[0] += s * src.el(0);
    val[1] += s * src.el(1);
    val[2] += s * src.el(2);
    val[3] += s * src.el(3);
    val[4] += s * src.el(4);
    val[5] += s * src.el(5);
    val[6] += s * src.el(6);
    val[7] += s * src.el(7);
    val[8] += s * src.el(8);
    val[9] += s * src.el(9);
    val[10] += s * src.el(10);
    val[11] += s * src.el(11);
    val[12] += s * src.el(12);
    val[13] += s * src.el(13);
    val[14] += s * src.el(14);
    val[15] += s * src.el(15);
    val[16] += s * src.el(16);
    val[17] += s * src.el(17);
    val[18] += s * src.el(18);
    val[19] += s * src.el(19);

    val[20] += s * src.el(20);
    val[21] += s * src.el(21);
    val[22] += s * src.el(22);
    val[23] += s * src.el(23);
    val[24] += s * src.el(24);
    val[25] += s * src.el(25);
    val[26] += s * src.el(26);
    val[27] += s * src.el(27);
    val[28] += s * src.el(28);
    val[29] += s * src.el(29);

    val[30] += s * src.el(30);
    val[31] += s * src.el(31);
    val[32] += s * src.el(32);
    val[33] += s * src.el(33);
    val[34] += s * src.el(34);
    val[35] += s * src.el(35);
    val[36] += s * src.el(36);
    val[37] += s * src.el(37);
    val[38] += s * src.el(38);
    val[39] += s * src.el(39);

    val[40] += s * src.el(40);
    val[41] += s * src.el(41);
    val[42] += s * src.el(42);
    val[43] += s * src.el(43);
    val[44] += s * src.el(44);
    val[45] += s * src.el(45);
    val[46] += s * src.el(46);
    val[47] += s * src.el(47);
    val[48] += s * src.el(48);
    val[49] += s * src.el(49);

    val[50] += s * src.el(50);
    val[51] += s * src.el(51);
    val[52] += s * src.el(52);
    val[53] += s * src.el(53);
    val[54] += s * src.el(54);
    val[55] += s * src.el(55);
    val[56] += s * src.el(56);
    val[57] += s * src.el(57);
    val[58] += s * src.el(58);
    val[59] += s * src.el(59);

    val[60] += s * src.el(60);
    val[61] += s * src.el(61);
    val[62] += s * src.el(62);
    val[63] += s * src.el(63);
  }
  else
  {
    const unsigned int size = n()*m();
    for (unsigned int i=0; i<size; i++)
      val[i] += s * src.el(i);
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
    val[0] += s * src.el(0);
    val[0] += s * src.el(1);
    val[0] += s * src.el(2);
    val[3] += s * src.el(3);
    val[3] += s * src.el(4);
    val[3] += s * src.el(5);
    val[6] += s * src.el(6);
    val[6] += s * src.el(7);
    val[6] += s * src.el(8);
  }
  else if ((n()==4) && (m()==4))
  {
    val[0] += s * src.el(0);
    val[0] += s * src.el(1);
    val[0] += s * src.el(2);
    val[0] += s * src.el(3);
    val[4] += s * src.el(4);
    val[4] += s * src.el(5);
    val[4] += s * src.el(6);
    val[4] += s * src.el(7);
    val[8] += s * src.el(8);
    val[8] += s * src.el(9);
    val[8] += s * src.el(10);
    val[8] += s * src.el(11);
    val[12] += s * src.el(12);
    val[12] += s * src.el(13);
    val[12] += s * src.el(14);
    val[12] += s * src.el(15);
  }
  else if ((n()==8) && (m()==8))
  {
    val[0] += s * src.el(0);
    val[0] += s * src.el(1);
    val[0] += s * src.el(2);
    val[0] += s * src.el(3);
    val[0] += s * src.el(4);
    val[0] += s * src.el(5);
    val[0] += s * src.el(6);
    val[0] += s * src.el(7);
    val[8] += s * src.el(8);
    val[8] += s * src.el(9);
    val[8] += s * src.el(10);
    val[8] += s * src.el(11);
    val[8] += s * src.el(12);
    val[8] += s * src.el(13);
    val[8] += s * src.el(14);
    val[8] += s * src.el(15);
    val[16] += s * src.el(16);
    val[16] += s * src.el(17);
    val[16] += s * src.el(18);
    val[16] += s * src.el(19);

    val[16] += s * src.el(20);
    val[16] += s * src.el(21);
    val[16] += s * src.el(22);
    val[16] += s * src.el(23);
    val[24] += s * src.el(24);
    val[24] += s * src.el(25);
    val[24] += s * src.el(26);
    val[24] += s * src.el(27);
    val[24] += s * src.el(28);
    val[24] += s * src.el(29);

    val[24] += s * src.el(30);
    val[24] += s * src.el(31);
    val[32] += s * src.el(32);
    val[32] += s * src.el(33);
    val[32] += s * src.el(34);
    val[32] += s * src.el(35);
    val[32] += s * src.el(36);
    val[32] += s * src.el(37);
    val[32] += s * src.el(38);
    val[32] += s * src.el(39);

    val[40] += s * src.el(40);
    val[40] += s * src.el(41);
    val[40] += s * src.el(42);
    val[40] += s * src.el(43);
    val[40] += s * src.el(44);
    val[40] += s * src.el(45);
    val[40] += s * src.el(46);
    val[40] += s * src.el(47);
    val[48] += s * src.el(48);
    val[48] += s * src.el(49);

    val[48] += s * src.el(50);
    val[48] += s * src.el(51);
    val[48] += s * src.el(52);
    val[48] += s * src.el(53);
    val[48] += s * src.el(54);
    val[48] += s * src.el(55);
    val[56] += s * src.el(56);
    val[56] += s * src.el(57);
    val[56] += s * src.el(58);
    val[56] += s * src.el(59);

    val[56] += s * src.el(60);
    val[56] += s * src.el(61);
    val[56] += s * src.el(62);
    val[56] += s * src.el(63);
  }
  else
  {
    const unsigned int size = n()*m();
    for (unsigned int i=0; i<size; i++)
      val[i] += s * src.el(i);
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
    val[0] += s * src.el(0);
    val[1] += s * src.el(3);
    val[2] += s * src.el(6);

    val[3] += s * src.el(1);
    val[4] += s * src.el(4);
    val[5] += s * src.el(7);

    val[6] += s * src.el(2);
    val[7] += s * src.el(5);
    val[8] += s * src.el(8);
  }
  else if ((n()==4) && (m()==4))
  {
    val[0] += s * src.el(0);
    val[1] += s * src.el(4);
    val[2] += s * src.el(8);
    val[3] += s * src.el(12);

    val[4] += s * src.el(1);
    val[5] += s * src.el(5);
    val[6] += s * src.el(9);
    val[7] += s * src.el(13);

    val[8] += s * src.el(2);
    val[9] += s * src.el(6);
    val[10] += s * src.el(10);
    val[11] += s * src.el(14);

    val[12] += s * src.el(3);
    val[13] += s * src.el(7);
    val[14] += s * src.el(11);
    val[15] += s * src.el(15);
  }
  else if ((n()==8) && (m()==8))
  {
    val[0] += s * src.el(0);
    val[1] += s * src.el(8);
    val[2] += s * src.el(16);
    val[3] += s * src.el(24);
    val[4] += s * src.el(32);
    val[5] += s * src.el(40);
    val[6] += s * src.el(48);
    val[7] += s * src.el(56);

    val[8] += s * src.el(1);
    val[9] += s * src.el(9);
    val[10] += s * src.el(17);
    val[11] += s * src.el(25);
    val[12] += s * src.el(33);
    val[13] += s * src.el(41);
    val[14] += s * src.el(49);
    val[15] += s * src.el(57);

    val[16] += s * src.el(2);
    val[17] += s * src.el(10);
    val[18] += s * src.el(18);
    val[19] += s * src.el(26);
    val[20] += s * src.el(34);
    val[21] += s * src.el(42);
    val[22] += s * src.el(50);
    val[23] += s * src.el(58);

    val[24] += s * src.el(3);
    val[25] += s * src.el(11);
    val[26] += s * src.el(19);
    val[27] += s * src.el(27);
    val[28] += s * src.el(35);
    val[29] += s * src.el(43);
    val[30] += s * src.el(51);
    val[31] += s * src.el(59);

    val[32] += s * src.el(4);
    val[33] += s * src.el(12);
    val[34] += s * src.el(20);
    val[35] += s * src.el(28);
    val[36] += s * src.el(36);
    val[37] += s * src.el(44);
    val[38] += s * src.el(52);
    val[39] += s * src.el(60);

    val[40] += s * src.el(5);
    val[41] += s * src.el(13);
    val[42] += s * src.el(21);
    val[43] += s * src.el(29);
    val[44] += s * src.el(37);
    val[45] += s * src.el(45);
    val[46] += s * src.el(53);
    val[47] += s * src.el(61);

    val[48] += s * src.el(6);
    val[49] += s * src.el(14);
    val[50] += s * src.el(22);
    val[51] += s * src.el(30);
    val[52] += s * src.el(38);
    val[53] += s * src.el(46);
    val[54] += s * src.el(54);
    val[55] += s * src.el(62);

    val[56] += s * src.el(7);
    val[57] += s * src.el(15);
    val[58] += s * src.el(23);
    val[59] += s * src.el(31);
    val[60] += s * src.el(39);
    val[61] += s * src.el(47);
    val[62] += s * src.el(55);
    val[63] += s * src.el(63);
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

  for (unsigned int i=0; i<dim_image; ++i)
    for (unsigned int j=0; j<dim_range; ++j)
      if (el(i,j) != m.el(i,j)) return false;
  return true;
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
double
FullMatrix<number>::norm2 () const
{
  double s = 0.;
  for (unsigned int i=0;i<dim_image*dim_range;++i)
    s += val[i]*val[i];
  return s;
}


template <typename number>
void FullMatrix<number>::clear ()
{
  number       *val_ptr = &val[0];
  const number *end_ptr = &val[n()*m()];
  while (val_ptr != end_ptr)
    *val_ptr++ = 0.;
};



template <typename number>
void
FullMatrix<number>::invert (const FullMatrix<number> &M)
{
  Assert (dim_range == dim_image, ExcNotQuadratic());
  Assert ((dim_range>=1) && (dim_range<=4), ExcNotImplemented(dim_range));
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
	    const double t4 = 1.0/(M.el(0,0)*M.el(1,1)-M.el(0,1)*M.el(1,0));
	    el(0,0) = M.el(1,1)*t4;
	    el(0,1) = -M.el(0,1)*t4;
	    el(1,0) = -M.el(1,0)*t4;
	    el(1,1) = M.el(0,0)*t4;
	    return;
      };
      
      case 3:
      {
	    const double t4 = M.el(0,0)*M.el(1,1),
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

	const double t14 = M.el(0,0)*M.el(1,1);
	const double t15 = M.el(2,2)*M.el(3,3);
	const double t17 = M.el(2,3)*M.el(3,2);
	const double t19 = M.el(0,0)*M.el(2,1);
	const double t20 = M.el(1,2)*M.el(3,3);
	const double t22 = M.el(1,3)*M.el(3,2);
	const double t24 = M.el(0,0)*M.el(3,1);
	const double t25 = M.el(1,2)*M.el(2,3);
	const double t27 = M.el(1,3)*M.el(2,2);
	const double t29 = M.el(1,0)*M.el(0,1);
	const double t32 = M.el(1,0)*M.el(2,1);
	const double t33 = M.el(0,2)*M.el(3,3);
	const double t35 = M.el(0,3)*M.el(3,2);
	const double t37 = M.el(1,0)*M.el(3,1);
	const double t38 = M.el(0,2)*M.el(2,3);
	const double t40 = M.el(0,3)*M.el(2,2);
	const double t42 = t14*t15-t14*t17-t19*t20+t19*t22+
			   t24*t25-t24*t27-t29*t15+t29*t17+
			   t32*t33-t32*t35-t37*t38+t37*t40;
	const double t43 = M.el(2,0)*M.el(0,1);
	const double t46 = M.el(2,0)*M.el(1,1);
	const double t49 = M.el(2,0)*M.el(3,1);
	const double t50 = M.el(0,2)*M.el(1,3);
	const double t52 = M.el(0,3)*M.el(1,2);
	const double t54 = M.el(3,0)*M.el(0,1);
	const double t57 = M.el(3,0)*M.el(1,1);
	const double t60 = M.el(3,0)*M.el(2,1);
	const double t63 = t43*t20-t43*t22-t46*t33+t46*t35+
			   t49*t50-t49*t52-t54*t25+t54*t27+
			   t57*t38-t57*t40-t60*t50+t60*t52;
	const double t65 = 1/(t42+t63);
	const double t71 = M.el(0,2)*M.el(2,1);
	const double t73 = M.el(0,3)*M.el(2,1);
	const double t75 = M.el(0,2)*M.el(3,1);
	const double t77 = M.el(0,3)*M.el(3,1);
	const double t81 = M.el(0,1)*M.el(1,2);
	const double t83 = M.el(0,1)*M.el(1,3);
	const double t85 = M.el(0,2)*M.el(1,1);
	const double t87 = M.el(0,3)*M.el(1,1);
	const double t101 = M.el(1,0)*M.el(2,2);
	const double t103 = M.el(1,0)*M.el(2,3);
	const double t105 = M.el(2,0)*M.el(1,2);
	const double t107 = M.el(2,0)*M.el(1,3);
	const double t109 = M.el(3,0)*M.el(1,2);
	const double t111 = M.el(3,0)*M.el(1,3);
	const double t115 = M.el(0,0)*M.el(2,2);
	const double t117 = M.el(0,0)*M.el(2,3);
	const double t119 = M.el(2,0)*M.el(0,2);
	const double t121 = M.el(2,0)*M.el(0,3);
	const double t123 = M.el(3,0)*M.el(0,2);
	const double t125 = M.el(3,0)*M.el(0,3);
	const double t129 = M.el(0,0)*M.el(1,2);
	const double t131 = M.el(0,0)*M.el(1,3);
	const double t133 = M.el(1,0)*M.el(0,2);
	const double t135 = M.el(1,0)*M.el(0,3);
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
      }
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
  iVector p(n());

  unsigned int i,j,k,r;
  double max, hr;

  for (i=0; i<n(); ++i) p(i) = i;

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
	  i = p(j) ; p(j) = p(r) ; p(r) = i;
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
      for (k=0; k<n(); ++k) hv(p(k)) = el(i,k);
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
  Assert (src.size() == dim_range, ExcDimensionMismatch(src.size(), dim_range));

  for (unsigned int j=0 ; j<n() ; ++j)
  {
    double sigma = 0;
    unsigned int i;
    for (i=j ; i<m() ; ++i) sigma += el(i,j)*el(i,j);
    if (fabs(sigma) < 1.e-15) return;
    double s = el(j,j);
    s = (s<0) ? sqrt(sigma) : -sqrt(sigma);
    double dj = s;

    double beta = 1./(s*el(j,j)-sigma);
    el(j,j) -= s;

    for (unsigned int k=j+1 ; k<n() ; ++k)
    {
      double sum = 0.;
      for (i=j ; i<m() ; ++i) sum += el(i,j)*el(i,k);
      sum *= beta;

      for (i=j ; i<m() ; ++i) el(i,k) += sum*el(i,j);
    }

    double sum = 0.;
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

  double sum = 0.;
  for (unsigned int i=n() ; i<m() ; ++i) sum += src(i) * src(i);
  return sqrt(sum);
}
