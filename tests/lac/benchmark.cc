//----------------------------  benchmark.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  benchmark.cc  ---------------------------


#include <base/logstream.h>
#include <lac/vector.h>
#include "quickmatrix.h"
#include <time.h>

#define ITER 100
main()
{
  Vector<double> u;
  Vector<double> v;

  clock_t start;
  clock_t diff;

  deallog << "Iterations: " << ITER << endl;

  for (unsigned int nx=32; nx<8192 ; nx*=2)
    {
      const unsigned int dim=(nx-1)*(nx-1);

      deallog << "size = " << nx << "  dim = " << dim << endl;

      start = clock();
      for (unsigned int i=0;i<ITER;i++)
	{
	  u.reinit(dim);
	  v.reinit(dim);
	}
      diff = clock()-start;
      deallog << "reinit: " << double(diff)/(2*ITER) << endl;

      start = clock();
      for (unsigned int i=0;i<ITER;i++)
	{
	  u = (double) i;
	}
      diff = clock()-start;
      deallog << "operator=(double): " << double(diff)/ITER << endl;

      QuickMatrix<double> A(nx,nx);

      start = clock();
      for (unsigned int i=0;i<ITER;i++)
	{
	  A.vmult(v,u);
	}
      diff = clock()-start;
      deallog << "vmult: " << double(diff)/ITER << endl;
    }
}


