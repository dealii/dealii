// $Id$

// This file is part of the DEAL Library
// DEAL is Copyright(1995) by
// Roland Becker, Guido Kanschat, Franz-Theo Suttmeier

#include <lac/ivector.h>



iVector::iVector()
{
  dim=maxdim=1;
  val = new int[1];
  Assert (val != 0, ExcOutOfMemory());
  val[0] = 0;
}

iVector::iVector(int n)
{
  Assert (n>0, ExcInvalidNumber(n));

  dim = n;
  maxdim = n;
  val = new int[maxdim];
  Assert (val != 0, ExcOutOfMemory());
  (*this) = 0;
}

iVector::iVector(const iVector& v)
{
  reinit(v.dim,1);
  int i;
  for (i=0;i<dim;i++) val[i] = v.val[i];
}

void iVector::reinit(int n, int fast)
{
  Assert (n>0, ExcInvalidNumber(n));
  if (n>maxdim)
  {
    delete[] val;
    val = new int[n];
    Assert (val != 0, ExcOutOfMemory());
    maxdim = n;
  }
  dim = n;
  if (!fast) (*this) = 0;
}

void iVector::reinit(const iVector& v, int fast)
{
  int n = v.n();
  if (n>maxdim)
  {
    delete[] val;
    val = new int[n];
    Assert (val != 0, ExcOutOfMemory());
    maxdim = n;
  }
  dim = n;
  if (!fast) (*this) = 0;
}

iVector::~iVector()
{
  delete[] val;
}

void iVector::add(const iVector& v)
{
  int i;
  for (i = 0; i < dim; i++) val[i] += v(i);
}

void iVector::add(int a, const iVector& v)
{
  int i;
  for (i = 0; i < dim; i++) val[i] += a * v(i);
}

void iVector::equ(int a, const iVector& u)
{
  int i;
  for (i=0;i<dim;i++) val[i] = a*u.val[i];
}

iVector& iVector::operator = (int s)
{
  int i;
//  if (s==0.) memset(val,sizeof(*val) * dim, 0);
  for (i=0;i<dim;i++) val[i] = s;
  return *this;
}

iVector& iVector::operator = (const iVector& v)
{
  if (v.dim != dim) reinit(v,1);

  int i;
  for (i=0;i<dim;i++) val[i] = v.val[i];
  return *this;
}
