//----------------------------  symmetric_tensor_07.cc  ---------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2005 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  symmetric_tensor_07.cc  ---------------------------

// in symmetric_tensor_06 we have established that contracting with a
// symmetric tensor by hand works as with a full tensor that is stored
// in non-symmetric form. here make sure that we can abbreviate the contraction 

#include "../tests.h"
#include <base/symmetric_tensor.h>
#include <base/tensor.h>
#include <base/logstream.h>
#include <fstream>
#include <iostream>


SymmetricTensor<2,1>
operator * (const SymmetricTensor<4,1> &t,
	    const SymmetricTensor<2,1> &s)
{
  const unsigned int dim = 1;
  SymmetricTensor<2,dim> tmp;
  tmp[0][0] = t[0][0][0][0] * s[0][0];
  return tmp;
}



SymmetricTensor<2,1>
operator * (const SymmetricTensor<2,1> &s,
	    const SymmetricTensor<4,1> &t)
{
  const unsigned int dim = 1;
  SymmetricTensor<2,dim> tmp;
  tmp[0][0] = t[0][0][0][0] * s[0][0];
  return tmp;
}



SymmetricTensor<2,2>
operator * (const SymmetricTensor<4,2> &t,
	    const SymmetricTensor<2,2> &s)
{
  const unsigned int dim = 2;
  SymmetricTensor<2,dim> tmp;

  for (unsigned int i=0; i<dim; ++i)
    for (unsigned int j=i; j<dim; ++j)
      tmp[i][j] = t[i][j][0][0] * s[0][0] +
		  t[i][j][1][1] * s[1][1] +
		  2 * t[i][j][0][1] * s[0][1];

  return tmp;
}



SymmetricTensor<2,2>
operator * (const SymmetricTensor<2,2> &s,
	    const SymmetricTensor<4,2> &t)
{
  const unsigned int dim = 2;
  SymmetricTensor<2,dim> tmp;

  for (unsigned int i=0; i<dim; ++i)
    for (unsigned int j=i; j<dim; ++j)
      tmp[i][j] = s[0][0] * t[0][0][i][j] * +
		  s[1][1] * t[1][1][i][j] +
		  2 * s[0][1] * t[0][1][i][j];

  return tmp;
}



SymmetricTensor<2,3>
operator * (const SymmetricTensor<4,3> &t,
	    const SymmetricTensor<2,3> &s)
{
  const unsigned int dim = 3;
  SymmetricTensor<2,dim> tmp;

  for (unsigned int i=0; i<dim; ++i)
    for (unsigned int j=i; j<dim; ++j)
      tmp[i][j] = t[i][j][0][0] * s[0][0] +
		  t[i][j][1][1] * s[1][1] +
		  t[i][j][2][2] * s[2][2] +
		  2 * t[i][j][0][1] * s[0][1] +
		  2 * t[i][j][0][2] * s[0][2] +
		  2 * t[i][j][1][2] * s[1][2];

  return tmp;
}



SymmetricTensor<2,3>
operator * (const SymmetricTensor<2,3> &s,
	    const SymmetricTensor<4,3> &t)
{
  const unsigned int dim = 3;
  SymmetricTensor<2,dim> tmp;

  for (unsigned int i=0; i<dim; ++i)
    for (unsigned int j=i; j<dim; ++j)
      tmp[i][j] = s[0][0] * t[0][0][i][j] +
		  s[1][1] * t[1][1][i][j] +
		  s[2][2] * t[2][2][i][j] +
		  2 * s[0][1] * t[0][1][i][j] +
		  2 * s[0][2] * t[0][2][i][j] +
		  2 * s[1][2] * t[1][2][i][j];

  return tmp;
}

  

template <int dim>
void test ()
{
  const double lambda = 5,
	       mu     = 7;
  SymmetricTensor<4,dim> ts;
  Tensor<4,dim>          ta;
  for (unsigned int i=0; i<dim; ++i)
    for (unsigned int j=0; j<dim; ++j)
      for (unsigned int k=0; k<dim; ++k)
	for (unsigned int l=0; l<dim; ++l)
	  ts[i][j][k][l] = ta[i][j][k][l] = (((i==k) && (j==l) ? mu : 0) +
					     ((i==l) && (j==k) ? mu : 0) +
					     ((i==j) && (k==l) ? lambda : 0));
	  
  SymmetricTensor<2,dim> as, bs;
  Tensor<2,dim>          aa, ba;

  for (unsigned int i=0; i<dim; ++i)
    for (unsigned int j=0; j<dim; ++j)
      as[i][j] = aa[i][j] = (1. + (i+1)*(j+1));

  bs = ts * as;
  double_contract (ba, ta, aa);
  
  for (unsigned int i=0; i<dim; ++i)
    for (unsigned int j=0; j<dim; ++j) 
      {
	Assert (as[i][j] == aa[i][j], ExcInternalError());
	Assert (bs[i][j] == ba[i][j], ExcInternalError());

	deallog << as[i][j] << ' ' << bs[i][i] << std::endl;
      }

				   // also test that the elasticity
				   // tensor is positive definite
  deallog << as * ts * as << std::endl;
  Assert (as * ts * as > 0, ExcInternalError());
}

  


int main ()
{
  std::ofstream logfile("symmetric_tensor_07.output");
  logfile.precision(3);
  deallog.attach(logfile);
  deallog.depth_console(0);

  test<2> ();
  test<3> ();
  
  deallog << "OK" << std::endl;
}
