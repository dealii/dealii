// $Id$

// This file is part of the DEAL Library
// DEAL is Copyright(1995) by
// Roland Becker, Guido Kanschat, Franz-Theo Suttmeier

#include <lac/dvector.h>
#include <math.h>

dVector::dVector() :
		dim(0),
		maxdim(0),
		val(0)
{}


dVector::dVector(int n)
{
  Assert (n>0, ExcInvalidNumber(n));

  dim = n;
  maxdim = n;
  if (n)
    {
      val = new double[maxdim];
      Assert (val != 0, ExcOutOfMemory());
      clear ();
    }
  else
    val = 0;
}


dVector::dVector(const dVector& v)
{
  dim = v.n();
  maxdim = v.n();

  if (dim)
    {
      val = new double[maxdim];
      Assert (val != 0, ExcOutOfMemory());
      for (int i=0; i<dim; i++)
	val[i] = v.val[i];
    }
  else
    val = 0;
}



void dVector::reinit(int n, int fast)
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



void dVector::reinit(const dVector& v, int fast)
{
  int n = v.n();
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



dVector::~dVector()
{
  if (val) delete[] val;
}



void dVector::clear () {
  for (int i=0; i<dim; ++i)
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
  for (int i=0; i<(dim/4); ++i) 
    {
      sum0 += val[4*i] * v.val[4*i];
      sum1 += val[4*i+1] * v.val[4*i+1];
      sum2 += val[4*i+2] * v.val[4*i+2];
      sum3 += val[4*i+3] * v.val[4*i+3];
    };
				   // add up remaining elements
  for (int i=(dim/4)*4; i<dim; ++i)
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
  for (int i=0; i<(dim/4); ++i) 
    {
      sum0 += val[4*i] * val[4*i];
      sum1 += val[4*i+1] * val[4*i+1];
      sum2 += val[4*i+2] * val[4*i+2];
      sum3 += val[4*i+3] * val[4*i+3];
    };
				   // add up remaining elements
  for (int i=(dim/4)*4; i<dim; ++i)
    sum0 += val[i] * val[i];
  
  return sum0+sum1+sum2+sum3;
}



void dVector::add(const double v)
{
  for (int i = 0; i < dim; i++) val[i] += v;
}



void dVector::add(const dVector& v)
{
  Assert (dim == v.dim, ExcDimensionsDontMatch(dim, v.dim));
  for (int i = 0; i < dim; i++) val[i] += v(i);
}



void dVector::add(double a, const dVector& v)
{
  Assert (dim == v.dim, ExcDimensionsDontMatch(dim, v.dim));
  for (int i = 0; i < dim; i++) val[i] += a * v(i);
}



void dVector::add(double a, const dVector& v, double b, const dVector& w)
{
  Assert (dim == v.dim, ExcDimensionsDontMatch(dim, v.dim));
  Assert (dim == w.dim, ExcDimensionsDontMatch(dim, w.dim));
  for (int i = 0; i < dim; i++) val[i] += a * v.val[i] + b * w.val[i];
}



void dVector::sadd(double x, const dVector& v)
{
  Assert (dim == v.dim, ExcDimensionsDontMatch(dim, v.dim));
  for (int i = 0; i < dim; i++) val[i] = x * val[i] + v.val[i];
}



void dVector::sadd(double x, double a, const dVector& v)
{
  Assert (dim == v.dim, ExcDimensionsDontMatch(dim, v.dim));
  for (int i = 0; i < dim; i++) val[i] = x * val[i] + a * v.val[i];
}



void dVector::sadd(double x, double a, const dVector& v, double b, const dVector& w)
{
  Assert (dim == v.dim, ExcDimensionsDontMatch(dim, v.dim));
  Assert (dim == w.dim, ExcDimensionsDontMatch(dim, w.dim));
  for (int i = 0; i < dim; i++) val[i] = x * val[i] + a * v.val[i] + b * w.val[i];
}



void dVector::sadd(double x, double a, const dVector& v, 
		   double b, const dVector& w, double c, const dVector& y)
{
  Assert (dim == v.dim, ExcDimensionsDontMatch(dim, v.dim));
  Assert (dim == w.dim, ExcDimensionsDontMatch(dim, w.dim));
  Assert (dim == y.dim, ExcDimensionsDontMatch(dim, y.dim));
  for (int i = 0; i < dim; i++)
    val[i] = x * val[i] + a * v.val[i] + b * w.val[i] 
	     + c * y.val[i];
}



void dVector::equ(double a, const dVector& u, double b, const dVector& v)
{
  Assert (dim == v.dim, ExcDimensionsDontMatch(dim, v.dim));
  for (int i=0;i<dim;i++) val[i] = a*u.val[i] + b*v.val[i];
}



void dVector::equ(double a)
{
  for (int i=0;i<dim;i++) val[i] *= a;
}



void dVector::equ(double a, const dVector& u)
{
  Assert (dim == u.dim, ExcDimensionsDontMatch(dim, u.dim));
  for (int i=0;i<dim;i++) val[i] = a*u.val[i];
}



dVector& dVector::operator = (double s)
{
  for (int i=0;i<dim;i++) val[i] = s;
  return *this;
}



dVector& dVector::operator = (const dVector& v)
{
  if (v.dim != dim) reinit(v,1);

  for (int i=0;i<dim;i++) val[i] = v.val[i];
  return *this;
}



void dVector::cadd(int i, const VectorBase& V, double s, int j)
{
  const dVector& v = (const dVector&) V;

  Assert ((i>=0) || (i<dim), ExcInvalidIndex(i,dim));
  Assert ((j>=0) || (j<v.dim), ExcInvalidIndex(j,v.dim));

  val[i] += s*v.val[j];
}



void dVector::cadd(int i, const VectorBase& V, double s, int j, double t, int k)
{
  const dVector& v = (const dVector&) V;

  Assert ((i>=0) || (i<dim), ExcInvalidIndex(i,dim));
  Assert ((j>=0) || (j<v.dim), ExcInvalidIndex(j,v.dim));
  Assert ((k>=0) || (k<v.dim), ExcInvalidIndex(k,v.dim));

  val[i] += s*v.val[j] + t*v.val[k];
}



void dVector::cadd(int i, const VectorBase& V, double s, int j,
		   double t, int k, double q, int l, double r, int m)
{
  const dVector& v = (const dVector&) V;

  Assert ((i>=0) || (i<dim), ExcInvalidIndex(i,dim));
  Assert ((j>=0) || (j<v.dim), ExcInvalidIndex(j,v.dim));
  Assert ((k>=0) || (k<v.dim), ExcInvalidIndex(k,v.dim));
  Assert ((l>=0) || (l<v.dim), ExcInvalidIndex(l,v.dim));
  Assert ((m>=0) || (m<v.dim), ExcInvalidIndex(m,v.dim));

  val[i] += s*v.val[j] + t*v.val[k] + q*v.val[l] + r*v.val[m];
}



void dVector::czero(int i)
{
  Assert ((i>=0) && (i<dim), ExcInvalidIndex(i,dim));
  val[i] = 0.;
}



void dVector::cequ(int i, const VectorBase& V, double s, int j)
{
  const dVector& v = (const dVector&) V;

  Assert ((i>=0) && (i<dim), ExcInvalidIndex(i,dim));
  Assert ((j>=0) && (j<v.dim), ExcInvalidIndex(j,v.dim));

  val[i] = s*v.val[j];
}



void dVector::cequ(int i, const VectorBase& V, double s, int j, double t, int k)
{
  const dVector& v = (const dVector&) V;

  Assert ((i>=0) || (i<dim), ExcInvalidIndex(i,dim));
  Assert ((j>=0) || (j<v.dim), ExcInvalidIndex(j,v.dim));
  Assert ((k>=0) || (k<v.dim), ExcInvalidIndex(k,v.dim));

  val[i] = s*v.val[j] + t*v.val[k];
}



void dVector::cequ(int i, const VectorBase& V, double s, int j,
		   double t, int k, double q, int l, double r, int m)
{
  const dVector& v = (const dVector&) V;

  Assert ((i>=0) || (i<dim), ExcInvalidIndex(i,dim));
  Assert ((j>=0) || (j<v.dim), ExcInvalidIndex(j,v.dim));
  Assert ((k>=0) || (k<v.dim), ExcInvalidIndex(k,v.dim));
  Assert ((l>=0) || (l<v.dim), ExcInvalidIndex(l,v.dim));
  Assert ((m>=0) || (m<v.dim), ExcInvalidIndex(m,v.dim));

  val[i] = s*v.val[j] + t*v.val[k] + q*v.val[l] + r*v.val[m];
}



const char* dVector::name() const
{
  return "dVector";
}



void dVector::print(FILE* f, const char* format) const
{
  if (!format) format = " %5.2f";
  for (int j=0;j<n();j++)
    fprintf(f, format, val[j]);
  fputc('\n',f);
}



void dVector::print(const char* format) const
{
  if (!format) format = " %5.2f";
  for (int j=0;j<n();j++)
    printf (format, val[j]);
  printf ("\n");
}



void dVector::print (ostream &out) const {
  for (int i=0; i<n(); ++i)
    out << val[i] << endl;
};
