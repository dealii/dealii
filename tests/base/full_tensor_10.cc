//----------------------------  full_tensor_10.cc  ---------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2005, 2009 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  full_tensor_10.cc  ---------------------------

// test the determinant code for n>3

#include "../tests.h"
#include <deal.II/base/tensor.h>
#include <deal.II/base/logstream.h>
#include <fstream>
#include <iomanip>


template <int dim>
void test ()
{
  Tensor<2,dim> t;

				   // choose the same symmetric tensor
				   // as in symmetric_tensor_10
  for (unsigned int i=0; i<dim; ++i)
    for (unsigned int j=0; j<dim; ++j)
      t[i][j] = 1. * rand() / RAND_MAX;

  for (unsigned int i=0; i<dim; ++i)
    for (unsigned int j=0; j<dim; ++j)
      deallog << "A[" << i+1 << ',' << j+1 << "] := " << t[i][j] << ';'
	      << std::endl;
  
  deallog << determinant(t) << std::endl;
}

  


int main ()
{
  std::ofstream logfile("full_tensor_10/output");
  deallog << std::setprecision(3);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test<4> ();
  test<5> ();
  test<6> ();
  test<7> ();
  test<8> ();
  
  deallog << "OK" << std::endl;
}
