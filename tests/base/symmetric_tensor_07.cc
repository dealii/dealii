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


  

template <int dim>
void test ()
{
  const double lambda = 7,
	       mu     = 5;
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

	deallog << as[i][j] << ' ' << bs[i][j] << std::endl;
      }

				   // test distributivity of
				   // multiplication
  Assert ((as*ts)*as == as*(ts*as), ExcInternalError());
  
  
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
