// $Id$

// This file is part of the DEAL Library
// DEAL is Copyright(1995) by
// Roland Becker, Guido Kanschat, Franz-Theo Suttmeier

#include <lac/dvector.h>
#include <math.h>

dVector::dVector() : val(0)
{
  dim=maxdim=0;
}

dVector::dVector(int n) : val(0)
{
  Assert (n>0, ExcInvalidNumber(n));

  dim = n;
  maxdim = n;
  if (n)
    {
      val = new double[maxdim];
      Assert (val != 0, ExcOutOfMemory());
      (*this) = 0.;
    }
}

dVector::dVector(const dVector& v) : val(0)
{
  dim = v.n();
  maxdim = v.n();

  if (dim)
    {
      val = new double[maxdim];
      Assert (val != 0, ExcOutOfMemory());
      for (int i=0;i<dim;i++) val[i] = v.val[i];
    }
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
  if (!fast) (*this) = 0.;
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
  if (!fast) (*this) = 0.;
}

dVector::~dVector()
{
  if (val) delete[] val;
}

double dVector::operator * (const dVector& v) const
{
  int i;
  double s;
  for (i = 0, s = 0.; i < dim; i++)
  {
	 s += val[i] * v.val[i];
  }
  return s;
}

void dVector::add(const double v)
{
  int i;
  for (i = 0; i < dim; i++) val[i] += v;
}

void dVector::add(const dVector& v)
{
  int i;
  for (i = 0; i < dim; i++) val[i] += v(i);
}

void dVector::add(double a, const dVector& v)
{
  int i;
  for (i = 0; i < dim; i++) val[i] += a * v(i);
}

void dVector::add(double a, const dVector& v, double b, const dVector& w)
{
  int i;
  for (i = 0; i < dim; i++) val[i] += a * v.val[i] + b * w.val[i];
}

void dVector::sadd(double x, const dVector& v)
{
  int i;
  for (i = 0; i < dim; i++) val[i] = x * val[i] + v.val[i];
}

void dVector::sadd(double x, double a, const dVector& v)
{
  int i;
  for (i = 0; i < dim; i++) val[i] = x * val[i] + a * v.val[i];
}

void dVector::sadd(double x, double a, const dVector& v, double b, const dVector& w)
{
  int i;
  for (i = 0; i < dim; i++) val[i] = x * val[i] + a * v.val[i] + b * w.val[i];
}

void dVector::sadd(double x, double a, const dVector& v, 
		   double b, const dVector& w, double c, const dVector& y)
{
  int i;
  for (i = 0; i < dim; i++) val[i] = x * val[i] + a * v.val[i] + b * w.val[i] 
			      + c * y.val[i];
}

void dVector::equ(double a, const dVector& u, double b, const dVector& v)
{
  int i;
  for (i=0;i<dim;i++) val[i] = a*u.val[i] + b*v.val[i];
}

void dVector::equ(double a)
{
  int i;
  for (i=0;i<dim;i++) val[i] *= a;
}

void dVector::equ(double a, const dVector& u)
{
  int i;
  for (i=0;i<dim;i++) val[i] = a*u.val[i];
}

dVector& dVector::operator = (double s)
{
  int i;
  for (i=0;i<dim;i++) val[i] = s;
  return *this;
}

dVector& dVector::operator = (const dVector& v)
{
  if (v.dim != dim) reinit(v,1);

  int i;
  for (i=0;i<dim;i++) val[i] = v.val[i];
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
  for (int j=0;j<n();j++) fprintf(f, format, val[j]);
  fputc('\n',f);
}

void dVector::print(const char* format) const
{
  if (!format) format = " %5.2f";
  for (int j=0;j<n();j++) printf (format, val[j]);
  printf ("\n");
}
