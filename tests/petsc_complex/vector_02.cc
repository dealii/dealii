// ---------------------------------------------------------------------
//
// Copyright (C) 2004 - 2016 by the deal.II authors
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



// check assignment of elements in Vector

#include "../tests.h"
#include <deal.II/lac/petsc_vector.h>
#include <deal.II/lac/vector.h>
#include <fstream>
#include <iostream>
#include <vector>


void test ()
{
  const unsigned int s = 10;
  PETScWrappers::Vector v(s);
  for (unsigned int k=0; k<v.size(); ++k)
    v(k) = k;

  v.compress (VectorOperation::insert);

  PETScWrappers::Vector v2(s);
  for (int k=0; k<v2.size(); ++k)
    v2(k) = PetscScalar (k,-k);

  v2.compress(VectorOperation::insert);

  // we now print the vector v, it is all real. Then print after
  // adding v2 which is complex, and then subtract v2 again to get the
  // original vector back.

  deallog << "before: " << std::endl;
  for (unsigned int k=0; k<s; ++k)
    deallog << "(" << v(k).real () << "," << v(k).imag () << "i) ";
  deallog << std::endl;

  v.add(1.0,v2);

  deallog << "after: " << std::endl;
  for (unsigned int k=0; k<s; ++k)
    deallog << "(" << v(k).real () << "," << v(k).imag () << "i) ";
  deallog << std::endl;

  v.add(-1.0,v2);

  deallog << "back to original: " << std::endl;
  for (unsigned int k=0; k<s; ++k)
    deallog << "(" << v(k).real () << "," << v(k).imag () << "i) ";
  deallog << std::endl;

  deallog << "OK" << std::endl;
}

int main (int argc, char **argv)
{
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  try
    {
      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
      {
        test ();
      }

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
