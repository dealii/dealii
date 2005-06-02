//----------------------------  symmetric_tensor_18.cc  ---------------------------
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
//----------------------------  symmetric_tensor_18.cc  ---------------------------

// compute double contraction between two rank-4 tensors

#include "../tests.h"
#include <base/symmetric_tensor.h>
#include <base/logstream.h>
#include <fstream>
#include <iostream>


template <int dim>
void test ()
{
  deallog << "dim=" << dim << std::endl;

  SymmetricTensor<4,dim> T, A, B;
  for (unsigned int i=0; i<dim; ++i)
    for (unsigned int j=0; j<dim; ++j)
      for (unsigned int k=0; k<dim; ++k)
	for (unsigned int l=0; l<dim; ++l)
	  {
					     // write some entries
					     // into the tensors. may
					     // be overwritten by
					     // subsequent writes, but
					     // who cares?
	    A[i][j][k][l] = (i+1)*(j+1)*(l+1)*(k+1);
	    B[i][j][k][l] = (i+2)*(j+3)*(l+4)*(k+5);
	  }
  
  T = A*B;

  for (unsigned int i=0; i<dim; ++i)
    for (unsigned int j=0; j<dim; ++j)
      for (unsigned int k=0; k<dim; ++k)
	for (unsigned int l=0; l<dim; ++l)
	  {
	    deallog << (int)T[i][j][k][l] << std::endl;

					     // calculate result by
					     // hand
	    double tmp = 0;
	    for (unsigned int a=0; a<dim; ++a)
	      for (unsigned int b=0; b<dim; ++b)
		tmp += A[i][j][a][b] * B[a][b][k][l];

	    Assert (std::fabs(T[i][j][k][l] - tmp) < 1e-14*tmp,
		    ExcInternalError());
	  }
}

  


int main ()
{
  std::ofstream logfile("symmetric_tensor_18.output");
  logfile.precision(3);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test<1> ();
  test<2> ();
  test<3> ();
  
  deallog << "OK" << std::endl;
}
