// $Id$

// This file is part of the DEAL Library
// DEAL is Copyright(1995) by
// Roland Becker, Guido Kanschat, Franz-Theo Suttmeier

#include <lac/ivector.h>



iVector::iVector() :
		dim(1),
		maxdim(1),
		val(new int[1])
{
  Assert (val != 0, ExcOutOfMemory());
  clear ();
}


iVector::iVector(int n) :
		dim(n),
		maxdim(n),
		val(new int[n])
{
  Assert (n>0, ExcInvalidNumber(n));

  Assert (val != 0, ExcOutOfMemory());
  clear ();
}


iVector::iVector(const iVector& v) :
		dim(v.dim),
		maxdim(v.maxdim),
		val(0)
{
  reinit(v.dim,1);
  int i;
  for (i=0; i<dim; i++)
    val[i] = v.val[i];
}



void iVector::reinit(int n, int fast)
{
  Assert (n>0, ExcInvalidNumber(n));
  if (n>maxdim)
  {
    if (val != 0) delete[] val;
    val = new int[n];
    Assert (val != 0, ExcOutOfMemory());
    maxdim = n;
  }
  dim = n;
  if (!fast)
    clear ();
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
  if (!fast)
    clear ();
}



iVector::~iVector()
{
  if (val)
    delete[] val;
}



void iVector::clear () {
  for (int i=0; i<dim; ++i)
    val[i] = 0;
}


void iVector::add(const iVector& v)
{
  for (int i=0; i<dim; ++i) val[i] += v(i);
}



void iVector::add(int a, const iVector& v)
{
  for (int i=0; i<dim; ++i) val[i] += a * v(i);
}



void iVector::equ(int a, const iVector& u)
{
  for (int i=0; i<dim; ++i) val[i] = a*u.val[i];
}



iVector& iVector::operator = (int s)
{
//  if (s==0.) memset(val,sizeof(*val) * dim, 0);
  for (int i=0; i<dim; ++i) val[i] = s;
  return *this;
}


iVector& iVector::operator = (const iVector& v)
{
  if (v.dim != dim) reinit(v,1);

  for (int i=0; i<dim; ++i) val[i] = v.val[i];
  return *this;
}
