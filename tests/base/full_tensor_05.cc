//----------------------------  full_tensor_05.cc  ---------------------------
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
//----------------------------  full_tensor_05.cc  ---------------------------

// make sure the tensor t_ijkl=delta_ik delta_jl
// actually maps a rank-2 tensor onto twice itself

#include "../tests.h"
#include <base/tensor.h>
#include <base/logstream.h>
#include <fstream>
#include <iostream>


template <int dim>
void test ()
{
  Tensor<4,dim> t;
  for (unsigned int i=0; i<dim; ++i)
    for (unsigned int j=0; j<dim; ++j)
      for (unsigned int k=0; k<dim; ++k)
	for (unsigned int l=0; l<dim; ++l)
	  t[i][j][k][l] = ((i==k) && (j==l) ? 1 : 0);
	  
  Tensor<2,dim> a, b;
  a[0][0] = 1;
  a[1][1] = 2;
  a[0][1] = 3;

  for (unsigned int i=0; i<dim; ++i)
    for (unsigned int j=0; j<dim; ++j)
      {
	double tmp_ij = 0;
	for (unsigned int k=0; k<dim; ++k)
	  for (unsigned int l=0; l<dim; ++l)
	    {
	      deallog << i << ' ' << j << ' ' << k << ' ' << l << ": "
		      << t[i][j][k][l] << ' ' << a[k][l]
		      << std::endl;
	      tmp_ij += t[i][j][k][l] * a[k][l];
	    }
	b[i][j] = tmp_ij;
      }  

  Assert (a == b, ExcInternalError());

				   // try the same thing with scaled
				   // tensors etc
  t *= 2;
  b.clear ();
  for (unsigned int i=0; i<dim; ++i)
    for (unsigned int j=0; j<dim; ++j)
      {
	double tmp_ij = 0;
	for (unsigned int k=0; k<dim; ++k)
	  for (unsigned int l=0; l<dim; ++l)
	    tmp_ij += t[i][j][k][l] * a[k][l];
	b[i][j] = tmp_ij;
      }

  Assert (a == b/2, ExcInternalError());
}

  


int main ()
{
  std::ofstream logfile("full_tensor_05.output");
  logfile.precision(3);
  deallog.attach(logfile);
  deallog.depth_console(0);

  test<2> ();
  test<3> ();
  
  deallog << "OK" << std::endl;
}
