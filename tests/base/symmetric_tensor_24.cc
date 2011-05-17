//----------------------------  symmetric_tensor_24.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2006, 2008, 2009 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  symmetric_tensor_24.cc  ---------------------------

// check SymmetricTensor<2,dim>::component_to_unrolled_index and the
// other way round

#include "../tests.h"
#include <deal.II/base/symmetric_tensor.h>
#include <deal.II/base/logstream.h>
#include <deal.II/lac/vector.h>
#include <fstream>
#include <iomanip>


template <int dim>
void check ()
{
  typedef SymmetricTensor<2,dim> S;
  for (unsigned int i=0; i<S::n_independent_components; ++i)
    {
      deallog << i << "  --  "
	      << S::unrolled_to_component_indices (i)
	      << std::endl;
      Assert (S::component_to_unrolled_index
	      (S::unrolled_to_component_indices (i))
	      ==
	      i,
	      ExcInternalError());
    }
}


int main ()
{
  std::ofstream logfile("symmetric_tensor_24/output");
  deallog << std::setprecision(3);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  check<1> ();
  check<2> ();
  check<3> ();

  deallog << "OK" << std::endl;
}
