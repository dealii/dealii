// ---------------------------------------------------------------------
//
// Copyright (C) 2004 - 2014 by the deal.II authors
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


// test setting some elements using a non-const matrix iterator and operator-=,
// and reading them back through the matrix itself

#include "../tests.h"
#include <deal.II/lac/sparse_matrix.h>
#include <fstream>
#include <iomanip>


void test ()
{
  SparsityPattern sp (5,5,3);
  for (unsigned int i=0; i<5; ++i)
    for (unsigned int j=0; j<5; ++j)
      if (((i+2*j+1) % 3 == 0)
          ||
          (i==j))
        sp.add (i,j);
  sp.compress ();

  SparseMatrix<double> m(sp);
  for (unsigned int i=0; i<5; ++i)
    for (unsigned int j=0; j<5; ++j)
      if (((i+2*j+1) % 3 == 0)
          ||
          (i==j))
        m.set(i,j,1.);

  SparseMatrix<double>::iterator i = m.begin();
  for (; i!=m.end(); ++i)
    i->value() -= i->row()*i->column();

  for (unsigned int i=0; i<5; ++i)
    for (unsigned int j=0; j<5; ++j)
      if (((i+2*j+1) % 3 == 0)
          ||
          (i==j))
        {
          deallog << i << ' ' << j << ' ' << m.el(i,j)
                  << std::endl;
          Assert (std::fabs(m.el(i,j)-(1.-i*j)) < 1e-14,
                  ExcInternalError());
        }

  deallog << "OK" << std::endl;
}



int main ()
{
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  try
    {
      test ();
    }
  catch (std::exception &exc)
    {
      deallog << std::endl << std::endl
              << "----------------------------------------------------"
              << std::endl;
      deallog << "Exception on processing: " << std::endl
              << exc.what() << std::endl
              << "Aborting!" << std::endl
              << "----------------------------------------------------"
              << std::endl;

      return 1;
    }
  catch (...)
    {
      deallog << std::endl << std::endl
              << "----------------------------------------------------"
              << std::endl;
      deallog << "Unknown exception!" << std::endl
              << "Aborting!" << std::endl
              << "----------------------------------------------------"
              << std::endl;
      return 1;
    };
}
