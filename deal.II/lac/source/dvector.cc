// $Id$

// This file is part of the DEAL Library
// DEAL is Copyright(1995) by
// Roland Becker, Guido Kanschat, Franz-Theo Suttmeier

#include <lac/dvector.h>
#include <cmath>
#include <algorithm>


static inline double sqr (const double x) {
  return x*x;
};



dVector::dVector () :
		dim(0),
		maxdim(0),
		val(0)
{}


dVector::dVector (const unsigned int n) :
		dim(0),
		maxdim(0),
		val(0)
{
  reinit (n, false);
}


dVector::dVector (const dVector& v) :
		dim(v.size()),
		maxdim(v.size()),
		val(0)
{
  if (dim)
    {
      val = new double[maxdim];
      Assert (val != 0, ExcOutOfMemory());
      copy (v.begin(), v.end(), begin());
    }
}



void dVector::reinit (const unsigned int n, const bool fast) {
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
      val = new double[n];
      Assert (val != 0, ExcOutOfMemory());
      maxdim = n;
    };
  dim = n;
  if (fast == false)
    clear ();
}



void dVector::reinit (const dVector& v, const bool fast) {
  reinit (v.size(), fast);
};




dVector::~dVector ()
{
  if (val) delete[] val;
}



void dVector::clear () {
  if (dim>0)
    fill (begin(), end(), 0.);
}



bool dVector::all_zero () const {
  Assert (dim!=0, ExcEmptyVector());
  
  const_iterator p = begin(),
		 e = end();
  while (p!=e)
    if (*p++ != 0.0)
      return false;
  return true;
};



double dVector::operator * (const dVector& v) const
{
  Assert (dim!=0, ExcEmptyVector());

  if (&v == this)
    return norm_sqr();
  
  Assert (dim == v.dim, ExcDimensionsDontMatch(dim, v.dim));
  
  double sum0 = 0,
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



double dVector::norm_sqr () const
{
  Assert (dim!=0, ExcEmptyVector());

  double sum0 = 0,
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
      sum0 += sqr(*ptr++);
      sum1 += sqr(*ptr++);
      sum2 += sqr(*ptr++);
      sum3 += sqr(*ptr++);
    };
				   // add up remaining elements
  while (ptr != end())
    sum0 += sqr(*ptr++);
  
  return sum0+sum1+sum2+sum3;
};



double dVector::mean_value () const
{
  Assert (dim!=0, ExcEmptyVector());

  double sum0 = 0,
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



double dVector::l1_norm () const
{
  Assert (dim!=0, ExcEmptyVector());

  double sum0 = 0,
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



double dVector::l2_norm () const
{
  return sqrt(norm_sqr());
};



double dVector::linfty_norm () const {
  Assert (dim!=0, ExcEmptyVector());

  double max0=0.,
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
    if (max0<val[i])
      max0 = val[i];

  return max (max(max0, max1),
	      max(max2, max3));
};
  




dVector& dVector::operator += (const dVector& v)
{
  Assert (dim!=0, ExcEmptyVector());

  add (v);
  return *this;
}



dVector& dVector::operator -= (const dVector& v)
{
  Assert (dim!=0, ExcEmptyVector());
  Assert (dim == v.dim, ExcDimensionsDontMatch(dim, v.dim));

  iterator i_ptr = begin(),
	   i_end = end();
  const_iterator v_ptr = v.begin();
  while (i_ptr!=i_end)
    *i_ptr++ -= *v_ptr++;

  return *this;
}



void dVector::add (const double v)
{
  Assert (dim!=0, ExcEmptyVector());

  iterator i_ptr = begin(),
	   i_end = end();
  while (i_ptr!=i_end)
    *i_ptr++ += v;
}



void dVector::add (const dVector& v)
{
  Assert (dim!=0, ExcEmptyVector());
  Assert (dim == v.dim, ExcDimensionsDontMatch(dim, v.dim));

  iterator i_ptr = begin(),
	   i_end = end();
  const_iterator v_ptr = v.begin();
  while (i_ptr!=i_end)
    *i_ptr++ += *v_ptr++;
}



void dVector::add (const double a, const dVector& v)
{
  Assert (dim!=0, ExcEmptyVector());
  Assert (dim == v.dim, ExcDimensionsDontMatch(dim, v.dim));

  iterator i_ptr = begin(),
	   i_end = end();
  const_iterator v_ptr = v.begin();
  while (i_ptr!=i_end)
    *i_ptr++ += a * *v_ptr++;
}



void dVector::add (const double a, const dVector& v,
		   const double b, const dVector& w)
{
  Assert (dim!=0, ExcEmptyVector());
  Assert (dim == v.dim, ExcDimensionsDontMatch(dim, v.dim));
  Assert (dim == w.dim, ExcDimensionsDontMatch(dim, w.dim));
  iterator i_ptr = begin(),
	   i_end = end();
  const_iterator v_ptr = v.begin(),
		 w_ptr = w.begin();
  while (i_ptr!=i_end)
    *i_ptr++ += a * *v_ptr++ + b * *w_ptr++;
}



void dVector::sadd (const double x, const dVector& v)
{
  Assert (dim!=0, ExcEmptyVector());
  Assert (dim == v.dim, ExcDimensionsDontMatch(dim, v.dim));
  iterator i_ptr = begin(),
	   i_end = end();
  const_iterator v_ptr = v.begin();
  for (; i_ptr!=i_end; ++i_ptr)
    *i_ptr = x * *i_ptr  + *v_ptr++;
}



void dVector::sadd (const double x, const double a, const dVector& v)
{
  Assert (dim!=0, ExcEmptyVector());
  Assert (dim == v.dim, ExcDimensionsDontMatch(dim, v.dim));
  iterator i_ptr = begin(),
	   i_end = end();
  const_iterator v_ptr = v.begin();
  for (; i_ptr!=i_end; ++i_ptr)
    *i_ptr = x * *i_ptr  +  a * *v_ptr++;
}



void dVector::sadd (const double x, const double a,
		    const dVector& v, const double b, const dVector& w)
{
  Assert (dim!=0, ExcEmptyVector());
  Assert (dim == v.dim, ExcDimensionsDontMatch(dim, v.dim));
  Assert (dim == w.dim, ExcDimensionsDontMatch(dim, w.dim));
  iterator i_ptr = begin(),
	   i_end = end();
  const_iterator v_ptr = v.begin(),
		 w_ptr = w.begin();
  for (; i_ptr!=i_end; ++i_ptr)
    *i_ptr = x * *i_ptr  +  a * *v_ptr++  + b * *w_ptr++;
}



void dVector::sadd (const double x, const double a,
		    const dVector& v, const double b,
		    const dVector& w, const double c, const dVector& y)
{
  Assert (dim!=0, ExcEmptyVector());
  Assert (dim == v.dim, ExcDimensionsDontMatch(dim, v.dim));
  Assert (dim == w.dim, ExcDimensionsDontMatch(dim, w.dim));
  Assert (dim == y.dim, ExcDimensionsDontMatch(dim, y.dim));
  iterator i_ptr = begin(),
	   i_end = end();
  const_iterator v_ptr = v.begin(),
		 w_ptr = w.begin(),
		 y_ptr = y.begin();
  
  for (; i_ptr!=i_end; ++i_ptr)
    *i_ptr = (x * *i_ptr)  +  (a * *v_ptr++)  +  (b * *w_ptr++)  + (c * *y_ptr++);
}



void dVector::scale (const double factor)
{
  Assert (dim!=0, ExcEmptyVector());

  iterator ptr=begin(), eptr=end();
  while (ptr!=eptr)
    *ptr++ *= factor;
}



void dVector::equ (const double a, const dVector& u,
		   const double b, const dVector& v)
{
  Assert (dim!=0, ExcEmptyVector());
  Assert (dim == u.dim, ExcDimensionsDontMatch(dim, u.dim));
  Assert (dim == v.dim, ExcDimensionsDontMatch(dim, v.dim));
  iterator i_ptr = begin(),
	   i_end = end();
  const_iterator u_ptr = u.begin(),
		 v_ptr = v.begin();
  while (i_ptr!=i_end)
    *i_ptr++ = a * *u_ptr++  + b * *v_ptr++;
}



void dVector::equ (const double a, const dVector& u)
{
  Assert (dim!=0, ExcEmptyVector());
  Assert (dim == u.dim, ExcDimensionsDontMatch(dim, u.dim));
  iterator i_ptr = begin(),
	   i_end = end();
  const_iterator u_ptr = u.begin();
  while (i_ptr!=i_end)
    *i_ptr++ = a * *u_ptr++;
}



void dVector::ratio (const dVector &a, const dVector &b) {
  Assert (dim!=0, ExcEmptyVector());
  Assert (a.dim == b.dim, ExcDimensionsDontMatch (a.dim, b.dim));

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



dVector& dVector::operator = (const double s)
{
  Assert (dim!=0, ExcEmptyVector());
  fill (begin(), end(), s);
  return *this;
}



dVector& dVector::operator = (const dVector& v)
{
  if (v.dim != dim)
    reinit (v.dim, true);

  if (dim!=0)
    copy (v.begin(), v.end(), begin());
  return *this;
}



void dVector::print (FILE* f, const char* format) const
{
  Assert (dim!=0, ExcEmptyVector());
  if (!format) format = " %5.2f";
  for (unsigned int j=0;j<size();j++)
    fprintf(f, format, val[j]);
  fputc('\n',f);
}



void dVector::print (const char* format) const
{
  Assert (dim!=0, ExcEmptyVector());
  if (!format) format = " %5.2f";
  for (unsigned int j=0;j<size();j++)
    printf (format, val[j]);
  printf ("\n");
}



void dVector::print (ostream &out) const {
  Assert (dim!=0, ExcEmptyVector());
  for (unsigned int i=0; i<size(); ++i)
    out << val[i] << endl;

  if (!out)
    throw GlobalExcIO ();
};




