// $Id$

// This file is part of the DEAL Library
// DEAL is Copyright(1995) by
// Roland Becker, Guido Kanschat, Franz-Theo Suttmeier

#include <lac/dvector.h>
#include <math.h>

dVector::dVector () :
		dim(0),
		maxdim(0),
		val(0)
{}


dVector::dVector (const unsigned int n) :
		dim(n),
		maxdim(n),
		val(0)
{
  Assert (n>0, ExcInvalidNumber(n));

  if (n)
    {
      val = new double[maxdim];
      Assert (val != 0, ExcOutOfMemory());
      clear ();
    }
}


dVector::dVector (const dVector& v) :
		VectorBase(v),
		dim(v.n()),
		maxdim(v.n()),
		val(0)
{
  if (dim)
    {
      val = new double[maxdim];
      Assert (val != 0, ExcOutOfMemory());
      for (unsigned int i=0; i<dim; i++)
	val[i] = v.val[i];
    }
}



void dVector::reinit (const unsigned int n, const bool fast)
{
  Assert (n>0, ExcInvalidNumber(n));
  if (n>maxdim)
  {
    if (val) delete[] val;
    val = new double[n];
    Assert (val != 0, ExcOutOfMemory());
    maxdim = n;
  }
  dim = n;
  if (!fast)
    clear ();
}



void dVector::reinit (const dVector& v, const bool fast)
{
  const unsigned int n = v.n();
  if (n>maxdim)
  {
    if (val) delete[] val;
    val = new double[n];
    Assert (val != 0, ExcOutOfMemory());
    maxdim = n;
  }
  dim = n;
  if (!fast)
    clear ();
}



dVector::~dVector ()
{
  if (val) delete[] val;
}



void dVector::clear () {
  for (unsigned int i=0; i<dim; ++i)
    val[i] = 0.;
}


double dVector::operator * (const dVector& v) const
{
  Assert (dim == v.dim, ExcDimensionsDontMatch(dim, v.dim));
  
  double sum0 = 0,
	 sum1 = 0,
	 sum2 = 0,
	 sum3 = 0;

				   // use modern processors better by
				   // allowing pipelined commands to be
				   // executed in parallel
  for (unsigned int i=0; i<(dim/4); ++i) 
    {
      sum0 += val[4*i] * v.val[4*i];
      sum1 += val[4*i+1] * v.val[4*i+1];
      sum2 += val[4*i+2] * v.val[4*i+2];
      sum3 += val[4*i+3] * v.val[4*i+3];
    };
				   // add up remaining elements
  for (unsigned int i=(dim/4)*4; i<dim; ++i)
    sum0 += val[i] * v.val[i];
  
  return sum0+sum1+sum2+sum3;
}



double dVector::norm_sqr () const
{
  double sum0 = 0,
	 sum1 = 0,
	 sum2 = 0,
	 sum3 = 0;

				   // use modern processors better by
				   // allowing pipelined commands to be
				   // executed in parallel
  for (unsigned int i=0; i<(dim/4); ++i) 
    {
      sum0 += val[4*i] * val[4*i];
      sum1 += val[4*i+1] * val[4*i+1];
      sum2 += val[4*i+2] * val[4*i+2];
      sum3 += val[4*i+3] * val[4*i+3];
    };
				   // add up remaining elements
  for (unsigned int i=(dim/4)*4; i<dim; ++i)
    sum0 += val[i] * val[i];
  
  return sum0+sum1+sum2+sum3;
};



double dVector::mean_value () const
{
  double sum0 = 0,
	 sum1 = 0,
	 sum2 = 0,
	 sum3 = 0;

				   // use modern processors better by
				   // allowing pipelined commands to be
				   // executed in parallel
  for (unsigned int i=0; i<(dim/4); ++i) 
    {
      sum0 += val[4*i];
      sum1 += val[4*i+1];
      sum2 += val[4*i+2];
      sum3 += val[4*i+3];
    };
				   // add up remaining elements
  for (unsigned int i=(dim/4)*4; i<dim; ++i)
    sum0 += val[i];
  
  return (sum0+sum1+sum2+sum3)/n();
};



double dVector::l1_norm () const
{
  double sum0 = 0,
	 sum1 = 0,
	 sum2 = 0,
	 sum3 = 0;

				   // use modern processors better by
				   // allowing pipelined commands to be
				   // executed in parallel
  for (unsigned int i=0; i<(dim/4); ++i) 
    {
      sum0 += fabs(val[4*i]);
      sum1 += fabs(val[4*i+1]);
      sum2 += fabs(val[4*i+2]);
      sum3 += fabs(val[4*i+3]);
    };
				   // add up remaining elements
  for (unsigned int i=(dim/4)*4; i<dim; ++i)
    sum0 += fabs(val[i]);
  
  return sum0+sum1+sum2+sum3;
};



double dVector::l2_norm () const
{
  double sum0 = 0,
	 sum1 = 0,
	 sum2 = 0,
	 sum3 = 0;

				   // use modern processors better by
				   // allowing pipelined commands to be
				   // executed in parallel
  for (unsigned int i=0; i<(dim/4); ++i) 
    {
      sum0 += val[4*i] * val[4*i];
      sum1 += val[4*i+1] * val[4*i+1];
      sum2 += val[4*i+2] * val[4*i+2];
      sum3 += val[4*i+3] * val[4*i+3];
    };
				   // add up remaining elements
  for (unsigned int i=(dim/4)*4; i<dim; ++i)
    sum0 += val[i] * val[i];
  
  return sqrt(sum0+sum1+sum2+sum3);
};





void dVector::add (const double v)
{
  for (unsigned int i = 0; i < dim; i++) val[i] += v;
}



void dVector::add (const dVector& v)
{
  Assert (dim == v.dim, ExcDimensionsDontMatch(dim, v.dim));
  for (unsigned int i = 0; i < dim; i++) val[i] += v(i);
}



void dVector::add (const double a, const dVector& v)
{
  Assert (dim == v.dim, ExcDimensionsDontMatch(dim, v.dim));
  for (unsigned int i = 0; i < dim; i++) val[i] += a * v(i);
}



void dVector::add (const double a, const dVector& v,
		   const double b, const dVector& w)
{
  Assert (dim == v.dim, ExcDimensionsDontMatch(dim, v.dim));
  Assert (dim == w.dim, ExcDimensionsDontMatch(dim, w.dim));
  for (unsigned int i = 0; i < dim; i++)
    val[i] += a * v.val[i] + b * w.val[i];
}



void dVector::sadd (const double x, const dVector& v)
{
  Assert (dim == v.dim, ExcDimensionsDontMatch(dim, v.dim));
  for (unsigned int i = 0; i < dim; i++)
    val[i] = x * val[i] + v.val[i];
}



void dVector::sadd (const double x, const double a, const dVector& v)
{
  Assert (dim == v.dim, ExcDimensionsDontMatch(dim, v.dim));
  for (unsigned int i = 0; i < dim; i++)
    val[i] = x * val[i] + a * v.val[i];
}



void dVector::sadd (const double x, const double a,
		    const dVector& v, const double b, const dVector& w)
{
  Assert (dim == v.dim, ExcDimensionsDontMatch(dim, v.dim));
  Assert (dim == w.dim, ExcDimensionsDontMatch(dim, w.dim));
  for (unsigned int i = 0; i < dim; i++)
    val[i] = x * val[i] + a * v.val[i] + b * w.val[i];
}



void dVector::sadd (const double x, const double a,
		    const dVector& v, const double b,
		    const dVector& w, const double c, const dVector& y)
{
  Assert (dim == v.dim, ExcDimensionsDontMatch(dim, v.dim));
  Assert (dim == w.dim, ExcDimensionsDontMatch(dim, w.dim));
  Assert (dim == y.dim, ExcDimensionsDontMatch(dim, y.dim));
  for (unsigned int i = 0; i < dim; i++)
    val[i] = x * val[i] + a * v.val[i] + b * w.val[i] 
	     + c * y.val[i];
}



void dVector::equ (const double a, const dVector& u,
		   const double b, const dVector& v)
{
  Assert (dim == v.dim, ExcDimensionsDontMatch(dim, v.dim));
  for (unsigned int i=0; i<dim; i++)
    val[i] = a*u.val[i] + b*v.val[i];
}



void dVector::equ (const double a)
{
  for (unsigned int i=0; i<dim; i++)
    val[i] *= a;
}



void dVector::equ (const double a, const dVector& u)
{
  Assert (dim == u.dim, ExcDimensionsDontMatch(dim, u.dim));
  for (unsigned int i=0; i<dim; i++)
    val[i] = a*u.val[i];
}



dVector& dVector::operator = (const double s)
{
  for (unsigned int i=0; i<dim; i++)
    val[i] = s;
  return *this;
}



dVector& dVector::operator = (const dVector& v)
{
  if (v.dim != dim) reinit(v,1);

  for (unsigned int i=0; i<dim; i++)
    val[i] = v.val[i];
  return *this;
}



void dVector::cadd (const unsigned int i, const VectorBase& V,
		    const double s, const unsigned int j)
{
  const dVector& v = (const dVector&) V;

  Assert (i<dim, ExcInvalidIndex(i,dim));
  Assert (j<v.dim, ExcInvalidIndex(j,v.dim));

  val[i] += s*v.val[j];
}



void dVector::cadd(unsigned int i, const VectorBase& V,
		   double s, unsigned int j,
		   double t, unsigned int k)
{
  const dVector& v = (const dVector&) V;

  Assert (i<dim, ExcInvalidIndex(i,dim));
  Assert (j<v.dim, ExcInvalidIndex(j,v.dim));
  Assert (k<v.dim, ExcInvalidIndex(k,v.dim));

  val[i] += s*v.val[j] + t*v.val[k];
}



void dVector::cadd (const unsigned int i, const VectorBase& V,
		    const double s, const unsigned int j,
		    const double t, const unsigned int k,
		    const double q, const unsigned int l,
		    const double r, const unsigned int m)
{
  const dVector& v = (const dVector&) V;

  Assert (i<dim, ExcInvalidIndex(i,dim));
  Assert (j<v.dim, ExcInvalidIndex(j,v.dim));
  Assert (k<v.dim, ExcInvalidIndex(k,v.dim));
  Assert (l<v.dim, ExcInvalidIndex(l,v.dim));
  Assert (m<v.dim, ExcInvalidIndex(m,v.dim));

  val[i] += s*v.val[j] + t*v.val[k] + q*v.val[l] + r*v.val[m];
}



void dVector::czero (const unsigned int i)
{
  Assert (i<dim, ExcInvalidIndex(i,dim));
  val[i] = 0.;
}



void dVector::cequ (const unsigned int i, const VectorBase& V,
		    const double s, const unsigned int j)
{
  const dVector& v = (const dVector&) V;

  Assert (i<dim, ExcInvalidIndex(i,dim));
  Assert (j<v.dim, ExcInvalidIndex(j,v.dim));

  val[i] = s*v.val[j];
}



void dVector::cequ (const unsigned int i, const VectorBase& V,
		    const double s, const unsigned int j,
		    const double t, const unsigned int k)
{
  const dVector& v = (const dVector&) V;

  Assert (i<dim, ExcInvalidIndex(i,dim));
  Assert (j<v.dim, ExcInvalidIndex(j,v.dim));
  Assert (k<v.dim, ExcInvalidIndex(k,v.dim));

  val[i] = s*v.val[j] + t*v.val[k];
}



void dVector::cequ (const unsigned int i, const VectorBase& V,
		    const double s, const unsigned int j,
		    const double t, const unsigned int k,
		    const double q, const unsigned int l,
		    const double r, const unsigned int m)
{
  const dVector& v = (const dVector&) V;

  Assert (i<dim, ExcInvalidIndex(i,dim));
  Assert (j<v.dim, ExcInvalidIndex(j,v.dim));
  Assert (k<v.dim, ExcInvalidIndex(k,v.dim));
  Assert (l<v.dim, ExcInvalidIndex(l,v.dim));
  Assert (m<v.dim, ExcInvalidIndex(m,v.dim));

  val[i] = s*v.val[j] + t*v.val[k] + q*v.val[l] + r*v.val[m];
}



const char* dVector::name () const
{
  return "dVector";
}



void dVector::print (FILE* f, const char* format) const
{
  if (!format) format = " %5.2f";
  for (unsigned int j=0;j<n();j++)
    fprintf(f, format, val[j]);
  fputc('\n',f);
}



void dVector::print (const char* format) const
{
  if (!format) format = " %5.2f";
  for (unsigned int j=0;j<n();j++)
    printf (format, val[j]);
  printf ("\n");
}



void dVector::print (ostream &out) const {
  for (unsigned int i=0; i<n(); ++i)
    out << val[i] << endl;
};
