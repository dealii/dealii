// ---------------------------------------------------------------------
//
// Copyright (C) 2004 - 2013 by the deal.II authors
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


/*
  test ::Vector::swap() (fixed in r 25668)
 */


#include "../tests.h"
#include <deal.II/base/utilities.h>
#include <deal.II/lac/trilinos_vector.h>
#include <fstream>
#include <iostream>
#include <vector>

void print(TrilinosWrappers::Vector &v)
{
  deallog << "size= " << v.size()
          << " el(0)= " << v(0)
          << " l2norm()= " << v.l2_norm() << std::endl;
}


void test ()
{
  TrilinosWrappers::Vector v(5);
  for (unsigned int i=0; i<v.size(); ++i)
    v(i) = 1;
  TrilinosWrappers::Vector w(9);
  for (unsigned int i=0; i<w.size(); ++i)
    w(i) = 2;


  deallog << "v: ";
  print(v);
  deallog << "w: ";
  print(w);

  deallog << "**swap**" << std::endl;

  swap(v,w);

  deallog << "v: ";
  print(v);
  deallog << "w: ";
  print(w);

  Assert (v.size()==9, ExcInternalError());

  deallog << "OK" << std::endl;
}



int main (int argc, char **argv)
{
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  Utilities::MPI::MPI_InitFinalize mpi_initialization (argc, argv);


  try
    {
      test ();
    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Exception on processing: " << std::endl
                << exc.what() << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;

      return 1;
    }
  catch (...)
    {
      std::cerr << std::endl << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Unknown exception!" << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    };
}
