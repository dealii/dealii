// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2013 by the deal.II authors
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



#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/lac/vector.h>
#include "quickmatrix.h"
#include <time.h>

#define ITER 100


int main()
{
  Vector<double> u;
  Vector<double> v;

  clock_t start;
  clock_t diff;

  deallog << "Iterations: " << ITER << std::endl;

  for (unsigned int nx=32; nx<8192 ; nx*=2)
    {
      const unsigned int dim=(nx-1)*(nx-1);

      deallog << "size = " << nx << "  dim = " << dim << std::endl;

      start = clock();
      for (unsigned int i=0; i<ITER; i++)
        {
          u.reinit(dim);
          v.reinit(dim);
        }
      diff = clock()-start;
      deallog << "reinit: " << double(diff)/(2*ITER) << std::endl;

      start = clock();
      for (unsigned int i=0; i<ITER; i++)
        {
          u = (double) i;
        }
      diff = clock()-start;
      deallog << "operator=(double): " << double(diff)/ITER << std::endl;

      QuickMatrix<double> A(nx,nx);

      start = clock();
      for (unsigned int i=0; i<ITER; i++)
        {
          A.vmult(v,u);
        }
      diff = clock()-start;
      deallog << "vmult: " << double(diff)/ITER << std::endl;
    }
}


