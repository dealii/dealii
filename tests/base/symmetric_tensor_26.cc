//----------------------------  symmetric_tensor_26.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2006, 2008, 2009, 2011 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  symmetric_tensor_26.cc  ---------------------------

// test multiplication with a Tensor<1,dim>

#include "../tests.h"
#include <base/symmetric_tensor.h>
#include <base/logstream.h>
#include <lac/vector.h>
#include <fstream>
#include <iomanip>


template <int dim>
void check ()
{
  SymmetricTensor<2,dim> S;
  for (unsigned int i=0; i<S.n_independent_components; ++i)
    S[S.unrolled_to_component_indices (i)] = std::rand () % 10;

  Tensor<1,dim> x;
  for (unsigned int i=0; i<dim; ++i)
    x[i] = std::rand() % 10;
  
  deallog << "S = " << S << std::endl;
  deallog << "x = " << x << std::endl;
  deallog << "S*x = " << S*x << std::endl;
}


int main ()
{
  std::ofstream logfile("symmetric_tensor_26/output");
  deallog << std::setprecision(3);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  check<1> ();
  check<2> ();
  check<3> ();
  
  deallog << "OK" << std::endl;
}
