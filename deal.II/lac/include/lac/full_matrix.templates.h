//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1999, 2000, 2001, 2002, 2003, 2004, 2005 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------
#ifndef __deal2__full_matrix_templates_h
#define __deal2__full_matrix_templates_h


#include <base/config.h>
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
                :
                Table<2,number> (n,n)
{}


template <typename number>
FullMatrix<number>::FullMatrix (const unsigned int m,
				const unsigned int n) :
		Table<2,number> (m, n)
{}


template <typename number>
FullMatrix<number>::FullMatrix (const unsigned int m,
                                const unsigned int n,
                                const number* entries)
                :
                Table<2,number> (m, n)
{
  this->fill (entries);
}


template <typename number>
FullMatrix<number>::FullMatrix (const FullMatrix &m)
                :
		Table<2,number> (m)
{}


template <typename number>
FullMatrix<number>&
FullMatrix<number>::operator = (const FullMatrix<number>& M)
{
  Table<2,number>::operator=(M);
  return *this;
}


template <typename number>
template <typename number2>
FullMatrix<number>&
FullMatrix<number>::operator = (const FullMatrix<number2>& M)
{
  TableBase<2,number>::operator=(M);
  return *this;
}


template <typename number>
bool
FullMatrix<number>::all_zero () const
{
  Assert (!this->empty(), ExcEmptyMatrix());
  
  const number* p = this->data();
  const number* const e = this->data() + this->n_elements();
  while (p!=e)
    if (*p++ != 0.0)
      return false;

  return true;
}



template <typename number>
FullMatrix<number> &
FullMatrix<number>::operator *= (const double factor)
{
  number       *p = &this->el(0,0);
  const number *e = &this->el(0,0) + n()*m();
  while (p != e)
    *p++ *= factor;

  return *this;
}



template <typename number>
FullMatrix<number> &
FullMatrix<number>::operator /= (const double factor)
{
  number       *p = &this->el(0,0);
  const number *e = &this->el(0,0) + n()*m();

  const number factor_inv = 1./factor;
  
  while (p != e)
    *p++ *= factor_inv;

  return *this;
}



template <typename number>
template <typename number2>
void
FullMatrix<number>::vmult (Vector<number2>& dst,
			   const Vector<number2>& src,
			   const bool adding) const
{
  Assert (!this->empty(), ExcEmptyMatrix());
  
  Assert(dst.size() == m(), ExcDimensionMismatch(dst.size(), m()));
  Assert(src.size() == n(), ExcDimensionMismatch(src.size(), n()));

  Assert (&src != &dst, ExcSourceEqualsDestination());

  if ((n()==3) && (m()==3))
  {
    number2 s;
    number2 s0,s1,s2;
    s   = src(0);
    s0  = s*this->data()[0];
    s1  = s*this->data()[3];
    s2  = s*this->data()[6];
    
    s   = src(1);
    s0 += s*this->data()[1];
    s1 += s*this->data()[4];
    s2 += s*this->data()[7];
    
    s   = src(2);
    s0 += s*this->data()[2];
    s1 += s*this->data()[5];
    s2 += s*this->data()[8];

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
    number2 s;
    number2 s0,s1,s2,s3;
    s = src(0);
    s0  = s*this->data()[0];
    s1  = s*this->data()[4];
    s2  = s*this->data()[8];
    s3  = s*this->data()[12];
    
    s = src(1);
    s0 += s*this->data()[1];
    s1 += s*this->data()[5];
    s2 += s*this->data()[9];
    s3 += s*this->data()[13];
    
    s = src(2);
    s0 += s*this->data()[2];
    s1 += s*this->data()[6];
    s2 += s*this->data()[10];
    s3 += s*this->data()[14];
    
    s = src(3);
    s0 += s*this->data()[3];
    s1 += s*this->data()[7];
    s2 += s*this->data()[11];
    s3 += s*this->data()[15];
    
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
    number2 s;
    number2 s0,s1,s2,s3,s4,s5,s6,s7;
    s = src(0);
    s0 = s*this->data()[0]; s1 = s*this->data()[8]; s2 = s*this->data()[16]; s3 = s*this->data()[24];
    s4 = s*this->data()[32]; s5 = s*this->data()[40]; s6 = s*this->data()[48]; s7 = s*this->data()[56];
    s = src(1);
    s0 += s*this->data()[1]; s1 += s*this->data()[9]; s2 += s*this->data()[17]; s3 += s*this->data()[25];
    s4 += s*this->data()[33]; s5 += s*this->data()[41]; s6 += s*this->data()[49]; s7 += s*this->data()[57];
    s = src(2);
    s0 += s*this->data()[2]; s1 += s*this->data()[10]; s2 += s*this->data()[18]; s3 += s*this->data()[26];
    s4 += s*this->data()[34]; s5 += s*this->data()[42]; s6 += s*this->data()[50]; s7 += s*this->data()[58];
    s = src(3);
    s0 += s*this->data()[3]; s1 += s*this->data()[11]; s2 += s*this->data()[19]; s3 += s*this->data()[27];
    s4 += s*this->data()[35]; s5 += s*this->data()[43]; s6 += s*this->data()[51]; s7 += s*this->data()[59];
    s = src(4);
    s0 += s*this->data()[4]; s1 += s*this->data()[12]; s2 += s*this->data()[20]; s3 += s*this->data()[28];
    s4 += s*this->data()[36]; s5 += s*this->data()[44]; s6 += s*this->data()[52]; s7 += s*this->data()[60];
    s = src(5);
    s0 += s*this->data()[5]; s1 += s*this->data()[13]; s2 += s*this->data()[21]; s3 += s*this->data()[29];
    s4 += s*this->data()[37]; s5 += s*this->data()[45]; s6 += s*this->data()[53]; s7 += s*this->data()[61];
    s = src(6);
    s0 += s*this->data()[6]; s1 += s*this->data()[14]; s2 += s*this->data()[22]; s3 += s*this->data()[30];
    s4 += s*this->data()[38]; s5 += s*this->data()[46]; s6 += s*this->data()[54]; s7 += s*this->data()[62];
    s = src(7);
    s0 += s*this->data()[7]; s1 += s*this->data()[15]; s2 += s*this->data()[23]; s3 += s*this->data()[31];
    s4 += s*this->data()[39]; s5 += s*this->data()[47]; s6 += s*this->data()[55]; s7 += s*this->data()[63];
    
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
    const number* e = this->data();
    const unsigned int size_m = m(),
		       size_n = n();
    if (!adding)
      {
	for (unsigned int i=0; i<size_m; ++i)
	  {
	    number2 s = 0.;
	    for (unsigned int j=0; j<size_n; ++j)
	      s += src(j) * *(e++);
	    dst(i) = s;
	  };
      }
    else
      {
	for (unsigned int i=0; i<size_m; ++i)
	  {
	    number2 s = 0.;
	    for (unsigned int j=0; j<size_n; ++j)
	      s += src(j) * *(e++);
	    dst(i) += s;
	  };
      };
  };
}



template <typename number>
template <typename number2>
void FullMatrix<number>::Tvmult (Vector<number2>       &dst,
				 const Vector<number2> &src,
				 const bool             adding) const
{
  Assert (!this->empty(), ExcEmptyMatrix());
  
  Assert(dst.size() == n(), ExcDimensionMismatch(dst.size(), n()));
  Assert(src.size() == m(), ExcDimensionMismatch(src.size(), m()));

  Assert (&src != &dst, ExcSourceEqualsDestination());

  const unsigned int size_m = m(),
		     size_n = n();

  if (!adding)
    {
      for (unsigned int i=0; i<size_n; ++i)
	{
	  number2 s = 0.;
	  for (unsigned int j=0; j<size_m; ++j)
	    s += src(j) * this->el(j,i);
	  dst(i) = s;
	};
    }
  else
    {
      for (unsigned int i=0; i<size_n; ++i)
	{
	  number2 s = 0.;
	  for (unsigned int j=0; j<size_m; ++j)
	    s += src(j) * this->el(j,i);
	  dst(i) += s;
	};
    };
}


template <typename number>
template <typename number2, typename number3>
double FullMatrix<number>::residual (Vector<number2>& dst,
				     const Vector<number2>& src,
				     const Vector<number3>& right) const
{
  Assert (!this->empty(), ExcEmptyMatrix());
  
  Assert(dst.size() == m(), ExcDimensionMismatch(dst.size(), m()));
  Assert(src.size() == n(), ExcDimensionMismatch(src.size(), n()));
  Assert(right.size() == m(), ExcDimensionMismatch(right.size(), m()));

  Assert (&src != &dst, ExcSourceEqualsDestination());

  number2 s, res = 0.;
  const unsigned int size_m = m(),
		     size_n = n();
  for (unsigned int i=0; i<size_n; ++i)
    {
      s = right(i);
      for (unsigned int j=0; j<size_m; ++j)
	s -= src(j) * this->el(i,j);
      dst(i) = s;
      res += s*s;
    }
  return std::sqrt(res);
}



template <typename number>
template <typename number2>
void FullMatrix<number>::forward (Vector<number2>       &dst,
				  const Vector<number2> &src) const
{
  Assert (!this->empty(), ExcEmptyMatrix());
  
  Assert (dst.size() == m(), ExcDimensionMismatch(dst.size(), m()));
  Assert (src.size() == n(), ExcDimensionMismatch(src.size(), n()));

  unsigned int i,j;
  unsigned int nu = ( (m()<n()) ? m() : n());
  number2 s;
  for (i=0; i<nu; ++i)
    {
      s = src(i);
      for (j=0; j<i; ++j) s -= dst(j) * this->el(i,j);
      dst(i) = s/this->el(i,i);
    }
}



template <typename number>
template <typename number2>
void FullMatrix<number>::backward (Vector<number2>       &dst,
				   const Vector<number2> &src) const
{
  Assert (!this->empty(), ExcEmptyMatrix());
  
  unsigned int j;
  unsigned int nu = (m()<n() ? m() : n());
  number2 s;
  for (int i=nu-1; i>=0; --i)
    {
      s = src(i);
      for (j=i+1; j<nu; ++j) s -= dst(j) * this->el(i,j);
      dst(i) = s/this->el(i,i);
    }
}



template <typename number>
template <typename number2>
void FullMatrix<number>::fill (const FullMatrix<number2> &src,
			       const unsigned int dst_offset_i,
			       const unsigned int dst_offset_j,
			       const unsigned int src_offset_i,
			       const unsigned int src_offset_j)
{
				   // Compute maximal size of copied block
  const unsigned int rows = (m() - dst_offset_i >= src.m() - src_offset_i)
			    ? src.m()
			    : m();
  const unsigned int cols = (n() - dst_offset_j >= src.n() - src_offset_j)
			    ? src.n()
			    : n();
  
  for (unsigned int i=0; i<rows ; ++i)
    for (unsigned int j=0; j<cols ; ++j)
      this->el(dst_offset_i+i,dst_offset_j+j)
	= src.el(src_offset_i+i,src_offset_j+j);
}


template <typename number>
template <typename number2>
void FullMatrix<number>::fill_permutation (const FullMatrix<number2> &src,
					   const std::vector<unsigned int>& p_rows,
					   const std::vector<unsigned int>& p_cols)
{
  Assert (p_rows.size() == this->n_rows(),
	  ExcDimensionMismatch (p_rows.size(), this->n_rows()));
  Assert (p_cols.size() == this->n_cols(),
	  ExcDimensionMismatch (p_cols.size(), this->n_cols()));

  for (unsigned int i=0;i<this->n_rows();++i)
    for (unsigned int j=0;j<this->n_cols();++j)
      this->el(i,j) = src(p_rows[i], p_cols[j]);
}



/*  template <typename number> */
/*  template <typename number2> */
/*  void FullMatrix<number>::fill (const number2* entries) */
/*  { */
/*      if (n_cols()*n_rows() != 0) */
/*        std::copy (entries, entries+n_rows()*n_cols(), this->data()); */
/*  } */



template <typename number>
void FullMatrix<number>::add_row (const unsigned int i,
				  const number s,
				  const unsigned int j)
{
  Assert (!this->empty(), ExcEmptyMatrix());
  
  for (unsigned int k=0; k<m(); ++k)
    this->el(i,k) += s*this->el(j,k);
}


template <typename number>
void FullMatrix<number>::add_row (const unsigned int i,
				  const number s,
				  const unsigned int j,
				  const number t,
				  const unsigned int k)
{
  Assert (!this->empty(), ExcEmptyMatrix());
  
  const unsigned int size_m = m();
  for (unsigned l=0; l<size_m; ++l)
    this->el(i,l) += s*this->el(j,l) + t*this->el(k,l);
}


template <typename number>
void FullMatrix<number>::add_col (const unsigned int i, const number s,
                                  const unsigned int j)
{
  Assert (!this->empty(), ExcEmptyMatrix());
  
  for (unsigned int k=0; k<n(); ++k)
    this->el(k,i) += s*this->el(k,j);
}


template <typename number>
void FullMatrix<number>::add_col (const unsigned int i, const number s,
                                  const unsigned int j, const number t,
                                  const unsigned int k)
{
  Assert (!this->empty(), ExcEmptyMatrix());
  
  for (unsigned int l=0; l<n(); ++l)
    this->el(l,i) += s*this->el(l,j) + t*this->el(l,k);
}



template <typename number>
void FullMatrix<number>::swap_row (const unsigned int i, const unsigned int j)
{
  Assert (!this->empty(), ExcEmptyMatrix());
  
  number s;
  for (unsigned int k=0; k<m(); ++k)
  {
    s = this->el(i,k); this->el(i,k) = this->el(j,k); this->el(j,k) = s;
  }
}


template <typename number>
void FullMatrix<number>::swap_col (const unsigned int i, const unsigned int j)
{
  Assert (!this->empty(), ExcEmptyMatrix());
  
  number s;
  for (unsigned int k=0; k<n(); ++k)
  {
    s = this->el(k,i); this->el(k,i) = this->el(k,j); this->el(k,j) = s;
  }
}


template <typename number>
void FullMatrix<number>::diagadd (const number src)
{
  Assert (!this->empty(), ExcEmptyMatrix());  
  Assert (m() == n(), ExcDimensionMismatch(m(),n()));
  
  for (unsigned int i=0; i<n(); ++i)
    this->el(i,i) += src;
}


template <typename number>
template <typename number2>
void FullMatrix<number>::mmult (FullMatrix<number2>       &dst,
				const FullMatrix<number2> &src,
				const bool                 adding) const
{
  Assert (!this->empty(), ExcEmptyMatrix());
  Assert (n() == src.m(), ExcDimensionMismatch(n(), src.m()));
  Assert (dst.n() == src.n(), ExcDimensionMismatch(dst.n(), src.n()));
  Assert (dst.m() == m(), ExcDimensionMismatch(m(), dst.m()));

  if (!adding)
    for (unsigned int i=0; i<m(); i++)
      for (unsigned int j=0; j<src.n(); j++)
	{
	  number2 s = 0.;
	  for (unsigned k=0; k<n(); k++)
	    s+= this->el(i,k) * src.el(k,j);
	  dst.el(i,j) = s;
	}
  else
    for (unsigned int i=0; i<m(); i++)
      for (unsigned int j=0; j<src.n(); j++)
	{
	  number2 s = 0.;
	  for (unsigned k=0; k<n(); k++)
	    s+= this->el(i,k) * src.el(k,j);
	  dst.el(i,j) += s;
	}
}



template <typename number>
template <typename number2>
void FullMatrix<number>::Tmmult (FullMatrix<number2>       &dst,
				 const FullMatrix<number2> &src,
				 const bool                 adding) const
{
  Assert (!this->empty(), ExcEmptyMatrix());  
  Assert (m() == src.m(), ExcDimensionMismatch(m(), src.m()));
  Assert (n() == dst.m(), ExcDimensionMismatch(n(), dst.m()));
  Assert (src.n() == dst.n(), ExcDimensionMismatch(src.n(), dst.n()));

  if (!adding)
    for (unsigned int i=0; i<n(); i++)
      for (unsigned int j=0; j<src.n(); j++)
	{
	  number2 s = 0;
	  for (unsigned int k=0; k<m(); k++)
	    s += this->el(k,i) * src.el(k,j);
	  dst.el(i,j) = s;
	}
  else
    for (unsigned int i=0; i<n(); i++)
      for (unsigned int j=0; j<src.n(); j++)
	{
	  number2 s = 0;
	  for (unsigned int k=0; k<m(); k++)
	    s += this->el(k,i) * src.el(k,j);
	  dst.el(i,j) += s;
	}
}



template <typename number>
template <typename number2>
number2 FullMatrix<number>::matrix_norm_square (const Vector<number2> &v) const
{
  Assert (!this->empty(), ExcEmptyMatrix());
  
  Assert(m() == v.size(), ExcDimensionMismatch(m(),v.size()));
  Assert(n() == v.size(), ExcDimensionMismatch(n(),v.size()));

  number2 sum = 0.;
  const unsigned int n_rows = m();
  const number *val_ptr = this->data();
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
}


template <typename number>
template <typename number2>
number2 FullMatrix<number>::matrix_scalar_product (const Vector<number2> &u,
						   const Vector<number2> &v) const
{
  Assert (!this->empty(), ExcEmptyMatrix());
  
  Assert(m() == u.size(), ExcDimensionMismatch(m(),v.size()));
  Assert(n() == v.size(), ExcDimensionMismatch(n(),v.size()));

  number2 sum = 0.;
  const unsigned int n_rows = m();
  const unsigned int n_cols = n();
  const number *val_ptr = this->data();
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
}



template <typename number>
void
FullMatrix<number>::symmetrize ()
{
  Assert (m() == n(), ExcNotQuadratic());
  
  const unsigned int N = m();
  for (unsigned int i=0; i<N; ++i)
    for (unsigned int j=i+1; j<N; ++j)
      {
	const number t = (this->el(i,j) + this->el(j,i)) / 2;
	this->el(i,j) = this->el(j,i) = t;
      };
}


template <typename number>
number FullMatrix<number>::l1_norm () const
{
  Assert (!this->empty(), ExcEmptyMatrix());
  
  number sum=0, max=0;
  const unsigned int n_rows = m(), n_cols = n();
  
  for (unsigned int col=0; col<n_cols; ++col)
    {
      sum=0;
      for (unsigned int row=0; row<n_rows; ++row)
	sum += std::fabs(this->el(row,col));
      if (sum > max)
	max = sum;
    }
  return max;
}


template <typename number>
number FullMatrix<number>::linfty_norm () const
{
  Assert (!this->empty(), ExcEmptyMatrix());
  
  number sum=0, max=0;
  const unsigned int n_rows = m(), n_cols = n();

  for (unsigned int row=0; row<n_rows; ++row)
    {
      sum=0;
      for (unsigned int col=0; col<n_cols; ++col)
	sum += std::fabs(this->el(row,col));
      if (sum > max)
	max = sum;
    }
  return max;
}


template <typename number>
void
FullMatrix<number>::print (std::ostream       &s,
			   const unsigned int  w,
			   const unsigned int  p) const
{
  Assert (!this->empty(), ExcEmptyMatrix());
  
  for (unsigned int i=0; i<this->m(); ++i)
    {
      for (unsigned int j=0; j<this->n(); ++j)
	s << std::setw(w) << std::setprecision(p) << this->el(i,j);
      s << std::endl;
    }
}


template <typename number>
template <typename number2>
void
FullMatrix<number>::add_scaled (const number               s,
                                const FullMatrix<number2> &src)
{
  Assert (!this->empty(), ExcEmptyMatrix());
  
  Assert (m() == src.m(), ExcDimensionMismatch(m(), src.m()));
  Assert (n() == src.n(), ExcDimensionMismatch(n(), src.n()));

  number* val = const_cast<number*> (this->data());
  const number2* srcval = src.data();
  
  if ((n()==3) && (m()==3))
    {
      val[0] += s * srcval[0];
      val[1] += s * srcval[1];
      val[2] += s * srcval[2];
      val[3] += s * srcval[3];
      val[4] += s * srcval[4];
      val[5] += s * srcval[5];
      val[6] += s * srcval[6];
      val[7] += s * srcval[7];
      val[8] += s * srcval[8];
    }
  else if ((n()==4) && (m()==4))
    {
      val[0] += s * srcval[0];
      val[1] += s * srcval[1];
      val[2] += s * srcval[2];
      val[3] += s * srcval[3];
      val[4] += s * srcval[4];
      val[5] += s * srcval[5];
      val[6] += s * srcval[6];
      val[7] += s * srcval[7];
      val[8] += s * srcval[8];
      val[9] += s * srcval[9];
      val[10] += s * srcval[10];
      val[11] += s * srcval[11];
      val[12] += s * srcval[12];
      val[13] += s * srcval[13];
      val[14] += s * srcval[14];
      val[15] += s * srcval[15];
    }
  else if ((n()==8) && (m()==8))
    {
      val[0] += s * srcval[0];
      val[1] += s * srcval[1];
      val[2] += s * srcval[2];
      val[3] += s * srcval[3];
      val[4] += s * srcval[4];
      val[5] += s * srcval[5];
      val[6] += s * srcval[6];
      val[7] += s * srcval[7];
      val[8] += s * srcval[8];
      val[9] += s * srcval[9];
      val[10] += s * srcval[10];
      val[11] += s * srcval[11];
      val[12] += s * srcval[12];
      val[13] += s * srcval[13];
      val[14] += s * srcval[14];
      val[15] += s * srcval[15];
      val[16] += s * srcval[16];
      val[17] += s * srcval[17];
      val[18] += s * srcval[18];
      val[19] += s * srcval[19];

      val[20] += s * srcval[20];
      val[21] += s * srcval[21];
      val[22] += s * srcval[22];
      val[23] += s * srcval[23];
      val[24] += s * srcval[24];
      val[25] += s * srcval[25];
      val[26] += s * srcval[26];
      val[27] += s * srcval[27];
      val[28] += s * srcval[28];
      val[29] += s * srcval[29];

      val[30] += s * srcval[30];
      val[31] += s * srcval[31];
      val[32] += s * srcval[32];
      val[33] += s * srcval[33];
      val[34] += s * srcval[34];
      val[35] += s * srcval[35];
      val[36] += s * srcval[36];
      val[37] += s * srcval[37];
      val[38] += s * srcval[38];
      val[39] += s * srcval[39];

      val[40] += s * srcval[40];
      val[41] += s * srcval[41];
      val[42] += s * srcval[42];
      val[43] += s * srcval[43];
      val[44] += s * srcval[44];
      val[45] += s * srcval[45];
      val[46] += s * srcval[46];
      val[47] += s * srcval[47];
      val[48] += s * srcval[48];
      val[49] += s * srcval[49];

      val[50] += s * srcval[50];
      val[51] += s * srcval[51];
      val[52] += s * srcval[52];
      val[53] += s * srcval[53];
      val[54] += s * srcval[54];
      val[55] += s * srcval[55];
      val[56] += s * srcval[56];
      val[57] += s * srcval[57];
      val[58] += s * srcval[58];
      val[59] += s * srcval[59];

      val[60] += s * srcval[60];
      val[61] += s * srcval[61];
      val[62] += s * srcval[62];
      val[63] += s * srcval[63];
    }
  else
    {
      const unsigned int size = n()*m();
      for (unsigned int i=0; i<size; i++)
        val[i] += s * srcval[i];
    }
}

template <typename number>
template <typename number2>
void FullMatrix<number>::add (const FullMatrix<number2> &src,
			      const double factor,
			      const unsigned int dst_offset_i,
			      const unsigned int dst_offset_j,
			      const unsigned int src_offset_i,
			      const unsigned int src_offset_j)
{
				   // Compute maximal size of copied block
  const unsigned int rows = (m() - dst_offset_i >= src.m() - src_offset_i)
			    ? src.m()
			    : m();
  const unsigned int cols = (n() - dst_offset_j >= src.n() - src_offset_j)
			    ? src.n()
			    : n();
  
  for (unsigned int i=0; i<rows ; ++i)
    for (unsigned int j=0; j<cols ; ++j)
      this->el(dst_offset_i+i,dst_offset_j+j)
	+= factor * src.el(src_offset_i+i,src_offset_j+j);
}



template <typename number>
template <typename number2>
void FullMatrix<number>::Tadd (const FullMatrix<number2> &src,
			      const double factor,
			      const unsigned int dst_offset_i,
			      const unsigned int dst_offset_j,
			      const unsigned int src_offset_i,
			      const unsigned int src_offset_j)
{
				   // Compute maximal size of copied block
  const unsigned int rows = (m() - dst_offset_i >= src.n() - src_offset_i)
			    ? src.n()
			    : m();
  const unsigned int cols = (n() - dst_offset_j >= src.m() - src_offset_j)
			    ? src.m()
			    : n();
  
  for (unsigned int i=0; i<rows ; ++i)
    for (unsigned int j=0; j<cols ; ++j)
      this->el(dst_offset_i+i,dst_offset_j+j)
	+= factor * src.el(src_offset_i+j,src_offset_j+i);
}



template <typename number>
template <typename number2>
void
FullMatrix<number>::add (const number               s,
			 const FullMatrix<number2> &src)
{
                                   // this function is deprecated. forward to
                                   // the other one
  add_scaled (s, src);
}



template <typename number>
template <typename number2>
void
FullMatrix<number>::add_diag (const number s, const FullMatrix<number2>& src)
{
  Assert (!this->empty(), ExcEmptyMatrix());
  
  Assert (m() == src.m(), ExcDimensionMismatch(m(), src.m()));
  Assert (n() == src.n(), ExcDimensionMismatch(n(), src.n()));

  number* val = const_cast<number*> (this->data());
  const number2* srcval = src.data();
  
  if ((n()==3) && (m()==3))
  {
    val[0] += s * srcval[0];
    val[0] += s * srcval[1];
    val[0] += s * srcval[2];
    val[3] += s * srcval[3];
    val[3] += s * srcval[4];
    val[3] += s * srcval[5];
    val[6] += s * srcval[6];
    val[6] += s * srcval[7];
    val[6] += s * srcval[8];
  }
  else if ((n()==4) && (m()==4))
  {
    val[0] += s * srcval[0];
    val[0] += s * srcval[1];
    val[0] += s * srcval[2];
    val[0] += s * srcval[3];
    val[4] += s * srcval[4];
    val[4] += s * srcval[5];
    val[4] += s * srcval[6];
    val[4] += s * srcval[7];
    val[8] += s * srcval[8];
    val[8] += s * srcval[9];
    val[8] += s * srcval[10];
    val[8] += s * srcval[11];
    val[12] += s * srcval[12];
    val[12] += s * srcval[13];
    val[12] += s * srcval[14];
    val[12] += s * srcval[15];
  }
  else if ((n()==8) && (m()==8))
  {
    val[0] += s * srcval[0];
    val[0] += s * srcval[1];
    val[0] += s * srcval[2];
    val[0] += s * srcval[3];
    val[0] += s * srcval[4];
    val[0] += s * srcval[5];
    val[0] += s * srcval[6];
    val[0] += s * srcval[7];
    val[8] += s * srcval[8];
    val[8] += s * srcval[9];
    val[8] += s * srcval[10];
    val[8] += s * srcval[11];
    val[8] += s * srcval[12];
    val[8] += s * srcval[13];
    val[8] += s * srcval[14];
    val[8] += s * srcval[15];
    val[16] += s * srcval[16];
    val[16] += s * srcval[17];
    val[16] += s * srcval[18];
    val[16] += s * srcval[19];

    val[16] += s * srcval[20];
    val[16] += s * srcval[21];
    val[16] += s * srcval[22];
    val[16] += s * srcval[23];
    val[24] += s * srcval[24];
    val[24] += s * srcval[25];
    val[24] += s * srcval[26];
    val[24] += s * srcval[27];
    val[24] += s * srcval[28];
    val[24] += s * srcval[29];

    val[24] += s * srcval[30];
    val[24] += s * srcval[31];
    val[32] += s * srcval[32];
    val[32] += s * srcval[33];
    val[32] += s * srcval[34];
    val[32] += s * srcval[35];
    val[32] += s * srcval[36];
    val[32] += s * srcval[37];
    val[32] += s * srcval[38];
    val[32] += s * srcval[39];

    val[40] += s * srcval[40];
    val[40] += s * srcval[41];
    val[40] += s * srcval[42];
    val[40] += s * srcval[43];
    val[40] += s * srcval[44];
    val[40] += s * srcval[45];
    val[40] += s * srcval[46];
    val[40] += s * srcval[47];
    val[48] += s * srcval[48];
    val[48] += s * srcval[49];

    val[48] += s * srcval[50];
    val[48] += s * srcval[51];
    val[48] += s * srcval[52];
    val[48] += s * srcval[53];
    val[48] += s * srcval[54];
    val[48] += s * srcval[55];
    val[56] += s * srcval[56];
    val[56] += s * srcval[57];
    val[56] += s * srcval[58];
    val[56] += s * srcval[59];

    val[56] += s * srcval[60];
    val[56] += s * srcval[61];
    val[56] += s * srcval[62];
    val[56] += s * srcval[63];
  }
  else
  {
    const unsigned int size = n()*m();
    for (unsigned int i=0; i<size; i++)
      val[i] += s * srcval[i];
  }
}


template <typename number>
template <typename number2>
void
FullMatrix<number>::Tadd (const number s, const FullMatrix<number2>& src)
{
  Assert (!this->empty(), ExcEmptyMatrix());
  
  Assert (m() == n(),     ExcNotQuadratic());
  Assert (m() == src.m(), ExcDimensionMismatch(m(), src.m()));
  Assert (n() == src.n(), ExcDimensionMismatch(n(), src.n()));

  number* val = const_cast<number*> (this->data());
  const number2* srcval = src.data();
  
  if ((n()==3) && (m()==3))
  {
    val[0] += s * srcval[0];
    val[1] += s * srcval[3];
    val[2] += s * srcval[6];

    val[3] += s * srcval[1];
    val[4] += s * srcval[4];
    val[5] += s * srcval[7];

    val[6] += s * srcval[2];
    val[7] += s * srcval[5];
    val[8] += s * srcval[8];
  }
  else if ((n()==4) && (m()==4))
  {
    val[0] += s * srcval[0];
    val[1] += s * srcval[4];
    val[2] += s * srcval[8];
    val[3] += s * srcval[12];

    val[4] += s * srcval[1];
    val[5] += s * srcval[5];
    val[6] += s * srcval[9];
    val[7] += s * srcval[13];

    val[8] += s * srcval[2];
    val[9] += s * srcval[6];
    val[10] += s * srcval[10];
    val[11] += s * srcval[14];

    val[12] += s * srcval[3];
    val[13] += s * srcval[7];
    val[14] += s * srcval[11];
    val[15] += s * srcval[15];
  }
  else if ((n()==8) && (m()==8))
  {
    val[0] += s * srcval[0];
    val[1] += s * srcval[8];
    val[2] += s * srcval[16];
    val[3] += s * srcval[24];
    val[4] += s * srcval[32];
    val[5] += s * srcval[40];
    val[6] += s * srcval[48];
    val[7] += s * srcval[56];

    val[8] += s * srcval[1];
    val[9] += s * srcval[9];
    val[10] += s * srcval[17];
    val[11] += s * srcval[25];
    val[12] += s * srcval[33];
    val[13] += s * srcval[41];
    val[14] += s * srcval[49];
    val[15] += s * srcval[57];

    val[16] += s * srcval[2];
    val[17] += s * srcval[10];
    val[18] += s * srcval[18];
    val[19] += s * srcval[26];
    val[20] += s * srcval[34];
    val[21] += s * srcval[42];
    val[22] += s * srcval[50];
    val[23] += s * srcval[58];

    val[24] += s * srcval[3];
    val[25] += s * srcval[11];
    val[26] += s * srcval[19];
    val[27] += s * srcval[27];
    val[28] += s * srcval[35];
    val[29] += s * srcval[43];
    val[30] += s * srcval[51];
    val[31] += s * srcval[59];

    val[32] += s * srcval[4];
    val[33] += s * srcval[12];
    val[34] += s * srcval[20];
    val[35] += s * srcval[28];
    val[36] += s * srcval[36];
    val[37] += s * srcval[44];
    val[38] += s * srcval[52];
    val[39] += s * srcval[60];

    val[40] += s * srcval[5];
    val[41] += s * srcval[13];
    val[42] += s * srcval[21];
    val[43] += s * srcval[29];
    val[44] += s * srcval[37];
    val[45] += s * srcval[45];
    val[46] += s * srcval[53];
    val[47] += s * srcval[61];

    val[48] += s * srcval[6];
    val[49] += s * srcval[14];
    val[50] += s * srcval[22];
    val[51] += s * srcval[30];
    val[52] += s * srcval[38];
    val[53] += s * srcval[46];
    val[54] += s * srcval[54];
    val[55] += s * srcval[62];

    val[56] += s * srcval[7];
    val[57] += s * srcval[15];
    val[58] += s * srcval[23];
    val[59] += s * srcval[31];
    val[60] += s * srcval[39];
    val[61] += s * srcval[47];
    val[62] += s * srcval[55];
    val[63] += s * srcval[63];
  }
  else
    Assert (false, ExcNotImplemented());
}


template <typename number>
bool
FullMatrix<number>::operator == (const FullMatrix<number> &M) const
{
				   // the matrices may either be both
				   // empty, or of same size and with
				   // same values, if they shall be
				   // equal
  bool result = (this->data()==0) && (M.data()==0);
  result = result || ((m()==M.m()) && (n()==M.n()) &&
	     std::equal (this->data(), this->data()+m()*n(),
			 M.data()));
  
  return result;
}


template <typename number>
double
FullMatrix<number>::determinant () const
{
  Assert (!this->empty(), ExcEmptyMatrix());
  
  Assert (this->n_cols() == this->n_rows(),
	  ExcDimensionMismatch(this->n_cols(), this->n_rows()));
  
  switch (this->n_cols()) 
    {
      case 1:
            return this->el(0,0);
      case 2:
            return this->el(0,0)*this->el(1,1) - this->el(1,0)*this->el(0,1);
      case 3:
            return  (this->el(0,0)*this->el(1,1)*this->el(2,2)
                     -this->el(0,0)*this->el(1,2)*this->el(2,1)
                     -this->el(1,0)*this->el(0,1)*this->el(2,2)
                     +this->el(1,0)*this->el(0,2)*this->el(2,1)
                     +this->el(2,0)*this->el(0,1)*this->el(1,2)
                     -this->el(2,0)*this->el(0,2)*this->el(1,1));
      default:
            Assert (false, ExcNotImplemented());
            return 0;
    };
}


template <typename number>
number
FullMatrix<number>::frobenius_norm () const
{
  Assert (!this->empty(), ExcEmptyMatrix());
  
  number s = 0.;
  for (unsigned int i=0; i<this->n_rows()*this->n_cols(); ++i)
    s += this->data()[i]*this->data()[i];
  return std::sqrt(s);
}



template <typename number>
number
FullMatrix<number>::norm2 () const
{
  return frobenius_norm();
}



template <typename number>
number
FullMatrix<number>::relative_symmetry_norm2 () const
{
  Assert (!this->empty(), ExcEmptyMatrix());
  
  number s = 0.;
  number a = 0.;
  for (unsigned int i=0; i<this->n_rows(); ++i)
    for (unsigned int j=0; j<this->n_cols(); ++j)
      {
	a += ((*this)(i,j)-(*this)(j,i))*((*this)(i,j)-(*this)(j,i));
	s += (*this)(i,j)*(*this)(i,j);
      }
  if (s!=0.)
    return std::sqrt(a)/std::sqrt(s);
  return 0;
}



template <typename number>
template <typename number2>
void
FullMatrix<number>::invert (const FullMatrix<number2> &M)
{
  Assert (!this->empty(), ExcEmptyMatrix());
  
  Assert (this->n_cols() == this->n_rows(),
	  ExcNotQuadratic());
  Assert (this->n_cols() == M.n_cols(),
          ExcDimensionMismatch(this->n_cols(), M.n_cols()));
  Assert (this->n_rows() == M.n_rows(),
	  ExcDimensionMismatch(this->n_rows(), M.n_rows()));
  
  switch (this->n_cols()) 
    {
      case 1:
	    this->el(0,0) = 1.0/M.el(0,0);
	    return;
      case 2:
					     // this is Maple output,
					     // thus a bit unstructured
      {
	    const number t4 = 1.0/(M.el(0,0)*M.el(1,1)-M.el(0,1)*M.el(1,0));
	    this->el(0,0) = M.el(1,1)*t4;
	    this->el(0,1) = -M.el(0,1)*t4;
	    this->el(1,0) = -M.el(1,0)*t4;
	    this->el(1,1) = M.el(0,0)*t4;
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
	    this->el(0,0) = (M.el(1,1)*M.el(2,2)-M.el(1,2)*M.el(2,1))*t07;
	    this->el(0,1) = -(M.el(0,1)*M.el(2,2)-M.el(0,2)*M.el(2,1))*t07;
	    this->el(0,2) = -(-M.el(0,1)*M.el(1,2)+M.el(0,2)*M.el(1,1))*t07;
	    this->el(1,0) = -(M.el(1,0)*M.el(2,2)-M.el(1,2)*M.el(2,0))*t07;
	    this->el(1,1) = (M.el(0,0)*M.el(2,2)-t04)*t07;
	    this->el(1,2) = -(t6-t00)*t07;
	    this->el(2,0) = -(-M.el(1,0)*M.el(2,1)+M.el(1,1)*M.el(2,0))*t07;
	    this->el(2,1) = -(M.el(0,0)*M.el(2,1)-t01)*t07;
	    this->el(2,2) = (t4-t8)*t07;
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
	this->el(0,0) = (M.el(1,1)*M.el(2,2)*M.el(3,3)-M.el(1,1)*M.el(2,3)*M.el(3,2)-
			 M.el(2,1)*M.el(1,2)*M.el(3,3)+M.el(2,1)*M.el(1,3)*M.el(3,2)+
			 M.el(3,1)*M.el(1,2)*M.el(2,3)-M.el(3,1)*M.el(1,3)*M.el(2,2))*t65;
	this->el(0,1) = -(M.el(0,1)*M.el(2,2)*M.el(3,3)-M.el(0,1)*M.el(2,3)*M.el(3,2)-
			  t71*M.el(3,3)+t73*M.el(3,2)+t75*M.el(2,3)-t77*M.el(2,2))*t65;
	this->el(0,2) = (t81*M.el(3,3)-t83*M.el(3,2)-t85*M.el(3,3)+t87*M.el(3,2)+
			 t75*M.el(1,3)-t77*M.el(1,2))*t65;
	this->el(0,3) = -(t81*M.el(2,3)-t83*M.el(2,2)-t85*M.el(2,3)+t87*M.el(2,2)+
			  t71*M.el(1,3)-t73*M.el(1,2))*t65;
	this->el(1,0) = -(t101*M.el(3,3)-t103*M.el(3,2)-t105*M.el(3,3)+t107*M.el(3,2)+
			  t109*M.el(2,3)-t111*M.el(2,2))*t65;
	this->el(1,1) = (t115*M.el(3,3)-t117*M.el(3,2)-t119*M.el(3,3)+t121*M.el(3,2)+
			 t123*M.el(2,3)-t125*M.el(2,2))*t65;
	this->el(1,2) = -(t129*M.el(3,3)-t131*M.el(3,2)-t133*M.el(3,3)+t135*M.el(3,2)+
			  t123*M.el(1,3)-t125*M.el(1,2))*t65;
	this->el(1,3) = (t129*M.el(2,3)-t131*M.el(2,2)-t133*M.el(2,3)+t135*M.el(2,2)+
			 t119*M.el(1,3)-t121*M.el(1,2))*t65;
	this->el(2,0) = (t32*M.el(3,3)-t103*M.el(3,1)-t46*M.el(3,3)+t107*M.el(3,1)+
			 t57*M.el(2,3)-t111*M.el(2,1))*t65;
	this->el(2,1) = -(t19*M.el(3,3)-t117*M.el(3,1)-t43*M.el(3,3)+t121*M.el(3,1)+
			  t54*M.el(2,3)-t125*M.el(2,1))*t65;
	this->el(2,2) = (t14*M.el(3,3)-t131*M.el(3,1)-t29*M.el(3,3)+t135*M.el(3,1)+
			 t54*M.el(1,3)-t125*M.el(1,1))*t65;
	this->el(2,3) = -(t14*M.el(2,3)-t131*M.el(2,1)-t29*M.el(2,3)+t135*M.el(2,1)+
			  t43*M.el(1,3)-t121*M.el(1,1))*t65;
	this->el(3,0) = -(t32*M.el(3,2)-t101*M.el(3,1)-t46*M.el(3,2)+t105*M.el(3,1)+
			  t57*M.el(2,2)-t109*M.el(2,1))*t65;
	this->el(3,1) = (t19*M.el(3,2)-t115*M.el(3,1)-t43*M.el(3,2)+t119*M.el(3,1)+
			 t54*M.el(2,2)-t123*M.el(2,1))*t65;
	this->el(3,2) = -(t14*M.el(3,2)-t129*M.el(3,1)-t29*M.el(3,2)+t133*M.el(3,1)+
			  t54*M.el(1,2)-t123*M.el(1,1))*t65;
	this->el(3,3) = (t14*M.el(2,2)-t129*M.el(2,1)-t29*M.el(2,2)+t133*M.el(2,1)+
			 t43*M.el(1,2)-t119*M.el(1,1))*t65;
	
	break;
      }

      default:
					     // if no inversion is
					     // hardcoded, fall back
					     // to use the
					     // Gauss-Jordan algorithm
	    *this = M;
	    gauss_jordan();
    };    
}


template <typename number>
template <typename somenumber>
void
FullMatrix<number>::precondition_Jacobi (Vector<somenumber>       &dst,
					 const Vector<somenumber> &src,
					 const number              om) const
{
  Assert (m() == n(), ExcNotQuadratic());
  Assert (dst.size() == n(), ExcDimensionMismatch (dst.size(), n()));
  Assert (src.size() == n(), ExcDimensionMismatch (src.size(), n()));

  const unsigned int n = src.size();
  somenumber              *dst_ptr = dst.begin();
  const somenumber        *src_ptr = src.begin();
  
  for (unsigned int i=0; i<n; ++i, ++dst_ptr, ++src_ptr)
    *dst_ptr = om * *src_ptr / this->el(i,i);
}



template <typename number>
void
FullMatrix<number>::print_formatted (std::ostream       &out,
				     const unsigned int  precision,
				     const bool          scientific,
				     const unsigned int  width_,
				     const char         *zero_string,
				     const double        denominator) const
{
  unsigned int width = width_;
  
  Assert ((!this->empty()) || (this->n_cols()+this->n_rows()==0),
	  ExcInternalError());
  
				   // set output format, but store old
				   // state
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
	if (this->el(i,j) != 0)
	  out << std::setw(width)
	      << this->el(i,j) * denominator << ' ';
	else
	  out << std::setw(width) << zero_string << ' ';
      out << std::endl;
    };

  AssertThrow (out, ExcIO());
				   // reset output format
  out.flags (old_flags);
  out.precision(old_precision);
}



template <typename number>
void
FullMatrix<number>::gauss_jordan ()
{
  Assert (!this->empty(), ExcEmptyMatrix());  
  Assert (this->n_cols() == this->n_rows(), ExcNotQuadratic());
  
				   // Gauss-Jordan-Algorithmus
				   // cf. Stoer I (4th Edition) p. 153
  const unsigned int N = n();

				   // first get an estimate of the
				   // size of the elements of this
				   // matrix, for later checks whether
				   // the pivot element is large
				   // enough, or whether we have to
				   // fear that the matrix is not
				   // regular
  double diagonal_sum = 0;
  for (unsigned int i=0; i<N; ++i)
    diagonal_sum += std::fabs(this->el(i,i));
  const double typical_diagonal_element = diagonal_sum/N;
  
  std::vector<unsigned int> p(N);
  for (unsigned int i=0; i<N; ++i)
    p[i] = i;

  for (unsigned int j=0; j<N; ++j)
    {
				       // pivot search: search that
				       // part of the line on and
				       // right of the diagonal for
				       // the largest element
      number       max = std::fabs(this->el(j,j));
      unsigned int r   = j;
      for (unsigned int i=j+1; i<N; ++i)
	{
	  if (std::fabs(this->el(i,j)) > max)
	    {
	      max = std::fabs(this->el(i,j));
	      r = i;
	    }
	}
				       // check whether the pivot is
				       // too small
      Assert(max > 1.e-16*typical_diagonal_element,
	     ExcNotRegular(max));
      
				       // row interchange
      if (r>j)
	{
	  for (unsigned int k=0; k<N; ++k)
	    std::swap (this->el(j,k), this->el(r,k));

	  std::swap (p[j], p[r]);
	}

				       // transformation
      const number hr = 1./this->el(j,j);
      this->el(j,j) = hr;
      for (unsigned int k=0; k<N; ++k)
	{
	  if (k==j) continue;
	  for (unsigned int i=0; i<N; ++i)
	    {
	      if (i==j) continue;
	      this->el(i,k) -= this->el(i,j)*this->el(j,k)*hr;
	    }
	}
      for (unsigned int i=0; i<N; ++i)
	{
	  this->el(i,j) *= hr;
	  this->el(j,i) *= -hr;
	}
      this->el(j,j) = hr;
    }
				   // column interchange
  std::vector<number> hv(N);
  for (unsigned int i=0; i<N; ++i)
    {
      for (unsigned int k=0; k<N; ++k)
	hv[p[k]] = this->el(i,k);
      for (unsigned int k=0; k<N; ++k)
	this->el(i,k) = hv[k];
    }
}

// QR-transformation cf. Stoer 1 4.8.2 (p. 191)

template <typename number>
template <typename number2>
void
FullMatrix<number>::householder(Vector<number2>& src)
{
  Assert (!this->empty(), ExcEmptyMatrix());
  
				   // m > n, src.n() = m
  Assert (this->n_cols() <= this->n_rows(),
	  ExcDimensionMismatch(this->n_cols(), this->n_rows()));
  Assert (src.size() == this->n_rows(),
	  ExcDimensionMismatch(src.size(), this->n_rows()));

  for (unsigned int j=0 ; j<n() ; ++j)
  {
    number2 sigma = 0;
    unsigned int i;
    for (i=j ; i<m() ; ++i) sigma += this->el(i,j)*this->el(i,j);
    if (std::fabs(sigma) < 1.e-15) return;
    number2 s = this->el(j,j);
    s = (s<0) ? std::sqrt(sigma) : -std::sqrt(sigma);
    number2 dj = s;

    number2 beta = 1./(s*this->el(j,j)-sigma);
    this->el(j,j) -= s;

    for (unsigned int k=j+1 ; k<n() ; ++k)
    {
      number2 sum = 0.;
      for (i=j ; i<m() ; ++i) sum += this->el(i,j)*this->el(i,k);
      sum *= beta;

      for (i=j ; i<m() ; ++i) this->el(i,k) += sum*this->el(i,j);
    }

    number2 sum = 0.;
    for (i=j ; i<m() ; ++i) sum += this->el(i,j)*src(i);
    sum *= beta;

    for (i=j ; i<m() ; ++i) src(i) += sum*this->el(i,j);
    this->el(j,j) = dj;
  }
}


template <typename number>
template <typename number2>
double
FullMatrix<number>::least_squares (Vector<number2>& dst,
                                   Vector<number2>& src)
{
  Assert (!this->empty(), ExcEmptyMatrix());
  
  // m > n, m = src.n, n = dst.n

  householder(src);
  backward(dst, src);

  number2 sum = 0.;
  for (unsigned int i=n() ; i<m() ; ++i) sum += src(i) * src(i);
  return std::sqrt(sum);
}



template <typename number>
unsigned int
FullMatrix<number>::memory_consumption () const
{
  return sizeof(*this) - sizeof (Table<2,number>)
    + Table<2,number>::memory_consumption();
}


#endif
