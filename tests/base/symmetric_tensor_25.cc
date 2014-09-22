// ---------------------------------------------------------------------
//
// Copyright (C) 2006 - 2013 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------


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
        S[S.unrolled_to_component_indices (i)] = Testing::rand() % 10;

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
  std::ofstream logfile("output");
  deallog << std::setprecision(3);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  check<3> ();

  deallog << "OK" << std::endl;
}
