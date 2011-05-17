//----------------------------  symmetric_tensor_25.cc  ---------------------------
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
//----------------------------  symmetric_tensor_25.cc  ---------------------------

// Test the invariants of tensors using the Cayley-Hamilton theorem,
// see http://en.wikipedia.org/wiki/Invariants_of_tensors

#include "../tests.h"
#include <deal.II/base/symmetric_tensor.h>
#include <deal.II/base/logstream.h>
#include <deal.II/lac/vector.h>
#include <fstream>
#include <iomanip>


template <int dim>
void check ()
{
  for (unsigned int round=0; round<10; ++round)
    {
      SymmetricTensor<2,dim> S;
      for (unsigned int i=0; i<S.n_independent_components; ++i)
	S[S.unrolled_to_component_indices (i)] = std::rand () % 10;

      deallog << "S = " << S << std::endl;
      deallog << "first invariant  = " << first_invariant(S) << std::endl;
      deallog << "second invariant = " << second_invariant(S) << std::endl;
      deallog << "third invariant  = " << third_invariant(S) << std::endl;

      Tensor<2,dim> S_cubed;
      for (unsigned int d=0; d<dim; ++d)
	for (unsigned int e=0; e<dim; ++e)
	  for (unsigned int f=0; f<dim; ++f)
	    for (unsigned int g=0; g<dim; ++g)
	      S_cubed[d][e] += S[d][f] * S[f][g] * S[g][e];

      Tensor<2,dim> S_squared;
      for (unsigned int d=0; d<dim; ++d)
	for (unsigned int e=0; e<dim; ++e)
	  for (unsigned int f=0; f<dim; ++f)
	    S_squared[d][e] += S[d][f] * S[f][e];

      Tensor<2,dim> R = S_cubed
			- first_invariant(S) * S_squared
			+ second_invariant(S) * S
			- third_invariant(S) * unit_symmetric_tensor<dim> ();
      deallog << R << std::endl;

      Assert (R.norm() < 1e-10, ExcInternalError());
    }
}


int main ()
{
  std::ofstream logfile("symmetric_tensor_25/output");
  deallog << std::setprecision(3);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  check<3> ();

  deallog << "OK" << std::endl;
}
