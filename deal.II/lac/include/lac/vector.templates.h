//----------------------------  vector.templates.h  ---------------------------
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
//----------------------------  vector.templates.h  ---------------------------
#ifndef __deal2__vector_templates_h
#define __deal2__vector_templates_h


#include <lac/vector.h>
#include <cmath>
#include <algorithm>



template <typename Number>
static inline Number sqr (const Number x)
{
  return x*x;
};


template <typename Number>
Vector<Number>::Vector () :
		dim(0),
		maxdim(0),
		val(0)
{}


template <typename Number>
Vector<Number>::Vector (const unsigned int n) :
		dim(0),
		maxdim(0),
		val(0)
{
  reinit (n, false);
}


template <typename Number>
Vector<Number>::Vector (const Vector<Number>& v) :
		dim(v.size()),
		maxdim(v.size()),
		val(0)
{
  if (dim)
    {
      val = new Number[maxdim];
      Assert (val != 0, ExcOutOfMemory());
      std::copy (v.begin(), v.end(), begin());
    }
}


// see the .h file for why this function was disabled
//
// template <typename Number>
// template <typename OtherNumber>
// Vector<Number>::Vector (const Vector<OtherNumber>& v) :
// 		dim(v.size()),
// 		maxdim(v.size()),
// 		val(0)
// {
//   if (dim)
//     {
//       val = new Number[maxdim];
//       Assert (val != 0, ExcOutOfMemory());
//       copy (v.begin(), v.end(), begin());
//     }
// }


template <typename Number>
Vector<Number>::~Vector ()
{
  if (val) delete[] val;
}



template <typename Number>
void Vector<Number>::reinit (const unsigned int n, const bool fast) {
  if (n==0) 
    {
      if (val) delete[] val;
      val = 0;
      maxdim = dim = 0;
      return;
    };
  
  if (n>maxdim)
    {
      if (val) delete[] val;
      val = new Number[n];
      Assert (val != 0, ExcOutOfMemory());
      maxdim = n;
    };
  dim = n;
  if (fast == false)
    clear ();
}


template <typename Number>
void Vector<Number>::reinit (const Vector<Number>& v, const bool fast)
{
  reinit (v.size(), fast);
};


template <typename Number>
void Vector<Number>::clear ()
{
  if (dim>0)
    std::fill (begin(), end(), 0.);
}



template <typename Number>
void
Vector<Number>::swap (Vector<Number> &v)
{
  std::swap (dim,    v.dim);
  std::swap (maxdim, v.maxdim);
  std::swap (val,    v.val);
};



template <typename Number>
bool Vector<Number>::all_zero () const
{
  Assert (dim!=0, ExcEmptyVector());
  
  const_iterator p = begin(),
		 e = end();
  while (p!=e)
    if (*p++ != 0.0)
      return false;
  return true;
};


template <typename Number>
Number Vector<Number>::operator * (const Vector<Number>& v) const
{
  Assert (dim!=0, ExcEmptyVector());

  if (&v == this)
    return norm_sqr();
  
  Assert (dim == v.dim, ExcDimensionMismatch(dim, v.dim));
  
  Number sum0 = 0,
	 sum1 = 0,
	 sum2 = 0,
	 sum3 = 0;

				   // use modern processors better by
				   // allowing pipelined commands to be
				   // executed in parallel
  const_iterator ptr  = begin(),
		 vptr = v.begin(),
		 eptr = ptr + (dim/4)*4;
  while (ptr!=eptr)
    {
      sum0 += (*ptr++ * *vptr++);
      sum1 += (*ptr++ * *vptr++);
      sum2 += (*ptr++ * *vptr++);
      sum3 += (*ptr++ * *vptr++);
    };
				   // add up remaining elements
  while (ptr != end())
    sum0 += *ptr++ * *vptr++;
    
  return sum0+sum1+sum2+sum3;
}


template <typename Number>
Number Vector<Number>::norm_sqr () const
{
  Assert (dim!=0, ExcEmptyVector());

  Number sum0 = 0,
	 sum1 = 0,
	 sum2 = 0,
	 sum3 = 0;

				   // use modern processors better by
				   // allowing pipelined commands to be
				   // executed in parallel
  const_iterator ptr  = begin(),
		 eptr = ptr + (dim/4)*4;
  while (ptr!=eptr)
    {
      sum0 += ::sqr(*ptr++);
      sum1 += ::sqr(*ptr++);
      sum2 += ::sqr(*ptr++);
      sum3 += ::sqr(*ptr++);
    };
				   // add up remaining elements
  while (ptr != end())
    sum0 += ::sqr(*ptr++);
  
  return sum0+sum1+sum2+sum3;
};


template <typename Number>
Number Vector<Number>::mean_value () const
{
  Assert (dim!=0, ExcEmptyVector());

  Number sum0 = 0,
	 sum1 = 0,
	 sum2 = 0,
	 sum3 = 0;

				   // use modern processors better by
				   // allowing pipelined commands to be
				   // executed in parallel
  const_iterator ptr  = begin(),
		 eptr = ptr + (dim/4)*4;
  while (ptr!=eptr)
    {
      sum0 += *ptr++;
      sum1 += *ptr++;
      sum2 += *ptr++;
      sum3 += *ptr++;
    };
				   // add up remaining elements
  while (ptr != end())
    sum0 += *ptr++;
  
  return (sum0+sum1+sum2+sum3)/size();
};


template <typename Number>
Number Vector<Number>::l1_norm () const
{
  Assert (dim!=0, ExcEmptyVector());

  Number sum0 = 0,
	 sum1 = 0,
	 sum2 = 0,
	 sum3 = 0;

				   // use modern processors better by
				   // allowing pipelined commands to be
				   // executed in parallel
  const_iterator ptr  = begin(),
		 eptr = ptr + (dim/4)*4;
  while (ptr!=eptr)
    {
      sum0 += fabs(*ptr++);
      sum1 += fabs(*ptr++);
      sum2 += fabs(*ptr++);
      sum3 += fabs(*ptr++);
    };
				   // add up remaining elements
  while (ptr != end())
    sum0 += fabs(*ptr++);
  
  return sum0+sum1+sum2+sum3;
};


template <typename Number>
Number Vector<Number>::l2_norm () const
{
  return sqrt(norm_sqr());
};


template <typename Number>
Number Vector<Number>::linfty_norm () const {
  Assert (dim!=0, ExcEmptyVector());

  Number max0=0.,
	 max1=0.,
	 max2=0.,
	 max3=0.;
  for (unsigned int i=0; i<(dim/4); ++i) 
    {
      if (max0<fabs(val[4*i]))   max0=fabs(val[4*i]);
      if (max1<fabs(val[4*i+1])) max1=fabs(val[4*i+1]);
      if (max2<fabs(val[4*i+2])) max2=fabs(val[4*i+2]);
      if (max3<fabs(val[4*i+3])) max3=fabs(val[4*i+3]);
    };
				   // add up remaining elements
  for (unsigned int i=(dim/4)*4; i<dim; ++i)
    if (max0<fabs(val[i]))
      max0 = fabs(val[i]);

  return std::max (std::max(max0, max1),
		   std::max(max2, max3));
};


template <typename Number>
Vector<Number>& Vector<Number>::operator += (const Vector<Number>& v)
{
  Assert (dim!=0, ExcEmptyVector());

  add (v);
  return *this;
}


template <typename Number>
Vector<Number>& Vector<Number>::operator -= (const Vector<Number>& v)
{
  Assert (dim!=0, ExcEmptyVector());
  Assert (dim == v.dim, ExcDimensionMismatch(dim, v.dim));

  iterator i_ptr = begin(),
	   i_end = end();
  const_iterator v_ptr = v.begin();
  while (i_ptr!=i_end)
    *i_ptr++ -= *v_ptr++;

  return *this;
}


template <typename Number>
void Vector<Number>::add (const Number v)
{
  Assert (dim!=0, ExcEmptyVector());

  iterator i_ptr = begin(),
	   i_end = end();
  while (i_ptr!=i_end)
    *i_ptr++ += v;
}


template <typename Number>
void Vector<Number>::add (const Vector<Number>& v)
{
  Assert (dim!=0, ExcEmptyVector());
  Assert (dim == v.dim, ExcDimensionMismatch(dim, v.dim));

  iterator i_ptr = begin(),
	   i_end = end();
  const_iterator v_ptr = v.begin();
  while (i_ptr!=i_end)
    *i_ptr++ += *v_ptr++;
}


template <typename Number>
void Vector<Number>::add (const Number a, const Vector<Number>& v)
{
  Assert (dim!=0, ExcEmptyVector());
  Assert (dim == v.dim, ExcDimensionMismatch(dim, v.dim));

  iterator i_ptr = begin(),
	   i_end = end();
  const_iterator v_ptr = v.begin();
  while (i_ptr!=i_end)
    *i_ptr++ += a * *v_ptr++;
}


template <typename Number>
void Vector<Number>::add (const Number a, const Vector<Number>& v,
			  const Number b, const Vector<Number>& w)
{
  Assert (dim!=0, ExcEmptyVector());
  Assert (dim == v.dim, ExcDimensionMismatch(dim, v.dim));
  Assert (dim == w.dim, ExcDimensionMismatch(dim, w.dim));
  iterator i_ptr = begin(),
	   i_end = end();
  const_iterator v_ptr = v.begin(),
		 w_ptr = w.begin();
  while (i_ptr!=i_end)
    *i_ptr++ += a * *v_ptr++ + b * *w_ptr++;
}


template <typename Number>
void Vector<Number>::sadd (const Number x, const Vector<Number>& v)
{
  Assert (dim!=0, ExcEmptyVector());
  Assert (dim == v.dim, ExcDimensionMismatch(dim, v.dim));
  iterator i_ptr = begin(),
	   i_end = end();
  const_iterator v_ptr = v.begin();
  for (; i_ptr!=i_end; ++i_ptr)
    *i_ptr = x * *i_ptr  + *v_ptr++;
}


template <typename Number>
void Vector<Number>::sadd (const Number x, const Number a, const Vector<Number>& v)
{
  Assert (dim!=0, ExcEmptyVector());
  Assert (dim == v.dim, ExcDimensionMismatch(dim, v.dim));
  iterator i_ptr = begin(),
	   i_end = end();
  const_iterator v_ptr = v.begin();
  for (; i_ptr!=i_end; ++i_ptr)
    *i_ptr = x * *i_ptr  +  a * *v_ptr++;
}


template <typename Number>
void Vector<Number>::sadd (const Number x, const Number a,
			   const Vector<Number>& v, const Number b, const Vector<Number>& w)
{
  Assert (dim!=0, ExcEmptyVector());
  Assert (dim == v.dim, ExcDimensionMismatch(dim, v.dim));
  Assert (dim == w.dim, ExcDimensionMismatch(dim, w.dim));
  iterator i_ptr = begin(),
	   i_end = end();
  const_iterator v_ptr = v.begin(),
		 w_ptr = w.begin();
  for (; i_ptr!=i_end; ++i_ptr)
    *i_ptr = x * *i_ptr  +  a * *v_ptr++  + b * *w_ptr++;
}


template <typename Number>
void Vector<Number>::sadd (const Number x, const Number a,
			   const Vector<Number>& v, const Number b,
			   const Vector<Number>& w, const Number c, const Vector<Number>& y)
{
  Assert (dim!=0, ExcEmptyVector());
  Assert (dim == v.dim, ExcDimensionMismatch(dim, v.dim));
  Assert (dim == w.dim, ExcDimensionMismatch(dim, w.dim));
  Assert (dim == y.dim, ExcDimensionMismatch(dim, y.dim));
  iterator i_ptr = begin(),
	   i_end = end();
  const_iterator v_ptr = v.begin(),
		 w_ptr = w.begin(),
		 y_ptr = y.begin();
  
  for (; i_ptr!=i_end; ++i_ptr)
    *i_ptr = (x * *i_ptr)  +  (a * *v_ptr++)  +  (b * *w_ptr++)  + (c * *y_ptr++);
}


template <typename Number>
void Vector<Number>::scale (const Number factor)
{
  Assert (dim!=0, ExcEmptyVector());

  iterator ptr=begin(), eptr=end();
  while (ptr!=eptr)
    *ptr++ *= factor;
}


template <typename Number>
void Vector<Number>::equ (const Number a, const Vector<Number>& u,
			  const Number b, const Vector<Number>& v)
{
  Assert (dim!=0, ExcEmptyVector());
  Assert (dim == u.dim, ExcDimensionMismatch(dim, u.dim));
  Assert (dim == v.dim, ExcDimensionMismatch(dim, v.dim));
  iterator i_ptr = begin(),
	   i_end = end();
  const_iterator u_ptr = u.begin(),
		 v_ptr = v.begin();
  while (i_ptr!=i_end)
    *i_ptr++ = a * *u_ptr++  + b * *v_ptr++;
}



template <typename Number>
void Vector<Number>::equ (const Number a, const Vector<Number>& u)
{
  Assert (dim!=0, ExcEmptyVector());
  Assert (dim == u.dim, ExcDimensionMismatch(dim, u.dim));
  iterator i_ptr = begin(),
	   i_end = end();
  const_iterator u_ptr = u.begin();
  while (i_ptr!=i_end)
    *i_ptr++ = a * *u_ptr++;
}



template <typename Number>
void Vector<Number>::ratio (const Vector<Number> &a, const Vector<Number> &b)
{
  Assert (dim!=0, ExcEmptyVector());
  Assert (a.dim == b.dim, ExcDimensionMismatch (a.dim, b.dim));

				   // no need to reinit with zeros, since
				   // we overwrite them anyway
  reinit (a.size(), true);
  iterator i_ptr = begin(),
	   i_end = end();
  const_iterator a_ptr = a.begin(),
		 b_ptr = b.begin();
  while (i_ptr!=i_end)
    *i_ptr++ = *a_ptr++ / *b_ptr++;
};



template <typename Number>
Vector<Number>& Vector<Number>::operator = (const Number s)
{
  Assert (dim!=0, ExcEmptyVector());
  std::fill (begin(), end(), s);
  return *this;
}



template <typename Number>
Vector<Number>&
Vector<Number>::operator = (const Vector<Number>& v)
{
  if (v.dim != dim)
    reinit (v.dim, true);
  if (dim!=0)
    std::copy (v.begin(), v.end(), begin());
  
  return *this;
}



template <typename Number>
template<typename Number2>
Vector<Number>&
Vector<Number>::operator = (const Vector<Number2>& v)
{
  if (v.size() != dim)
    reinit (v.size(), true);
  if (dim!=0)
    std::copy (v.begin(), v.end(), begin());
  
  return *this;
}



template <typename Number>
void Vector<Number>::print (FILE* f, const char* format) const
{
  Assert (dim!=0, ExcEmptyVector());
  if (!format) format = " %5.2f";
  for (unsigned int j=0;j<size();j++)
    fprintf(f, format, val[j]);
  fputc('\n',f);
}



template <typename Number>
void Vector<Number>::print (const char* format) const
{
  Assert (dim!=0, ExcEmptyVector());
  if (!format) format = " %5.2f";
  for (unsigned int j=0;j<size();j++)
    printf (format, val[j]);
  printf ("\n");
}



template <typename Number>
void Vector<Number>::print (std::ostream      &out,
			    const unsigned int precision,
			    const bool         scientific,
			    const bool         across) const
{
  Assert (dim!=0, ExcEmptyVector());

  out.precision (precision);
  if (scientific)
    out.setf (std::ios::scientific, std::ios::floatfield);
  else
    out.setf (std::ios::fixed, std::ios::floatfield);

  if (across)
    for (unsigned int i=0; i<size(); ++i)
      out << val[i] << ' ';
  else
    for (unsigned int i=0; i<size(); ++i)
      out << val[i] << std::endl;
  out << std::endl;
  
  AssertThrow (out, ExcIO());
};



template <typename Number>
void Vector<Number>::block_write (std::ostream &out) const
{
  AssertThrow (out, ExcIO());

				   // other version of the following
				   //  out << size() << std::endl << '[';
				   // reason: operator<< seems to use
				   // some resources that lead to
				   // problems in a multithreaded
				   // environment
  const unsigned int sz = size();
  char buf[16];
  
  sprintf(buf, "%d", sz);
  strcat(buf, "\n[");
  
  out.write(buf, strlen(buf));
  out.write (reinterpret_cast<const char*>(begin()),
	     reinterpret_cast<const char*>(end())
	     - reinterpret_cast<const char*>(begin()));
  
				   // out << ']';
  const char outro = ']';
  out.write (&outro, 1);
  
  AssertThrow (out, ExcIO());
};



template <typename Number>
void Vector<Number>::block_read (std::istream &in)
{
  AssertThrow (in, ExcIO());

  unsigned int sz;

  char buf[16];
  

  in.getline(buf,16,'\n');
  sz=atoi(buf);
  
				   // fast initialization, since the
				   // data elements are overwritten anyway
  reinit (sz, true);     

  char c;
				   //  in >> c;
  in.read (&c, 1);
  AssertThrow (c=='[', ExcIO());
  
  in.read (reinterpret_cast<char*>(begin()),
	   reinterpret_cast<const char*>(end())
	   - reinterpret_cast<const char*>(begin()));
  
				   //  in >> c;
  in.read (&c, 1);
  AssertThrow (c==']', ExcIO());
  AssertThrow (in, ExcIO());
}



template <typename Number>
unsigned int
Vector<Number>::memory_consumption () const
{
  return sizeof(*this) + (maxdim * sizeof(Number));
};


#endif
