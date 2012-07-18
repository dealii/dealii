//----------------------------  tensor_06.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2006, 2012 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  tensor_06.cc  ---------------------------

// check Tensor<1,dim>::operator*(Tensor<1,dim>)

#include "../tests.h"
#include <deal.II/base/tensor.h>
#include <deal.II/base/logstream.h>
#include <fstream>
#include <iomanip>

template <int dim>
void test_tensor ()
{
  Tensor<1,dim> t1, t2;
  for (unsigned int i=0; i<dim; ++i)
    {
      t1[i] = i+1;
      t2[i] = 4.*i-10.;
    }
  double res = t1 * t2;
  deallog << "dim = " << dim << ": " << res << std::endl;
}

int main ()
{
  std::ofstream logfile("tensor_06/output");
  deallog << std::setprecision(3);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test_tensor<1>();
  test_tensor<2>();
  test_tensor<3>();
  test_tensor<4>();
  test_tensor<7>();
  deallog << "OK" << std::endl;
}
