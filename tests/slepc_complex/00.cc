// ---------------------------------------------------------------------
// $Id$
//
// Copyright (C) 2013 by the deal.II authors
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



// just check initialising SLEPc can be done and that it initialises
// PETSc in the way we expect, ie. *a* PETSc object exist. This is an
// approximate copy of slepc/00.cc..

#include "../tests.h"
#include <deal.II/lac/slepc_solver.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/numbers.h>
#include <fstream>
#include <iostream>


std::ofstream logfile ("output");

int main (int argc, char **argv)
{
  deallog.attach (logfile);
  deallog.depth_console (1);

  try
    {
      logfile << "Initializing SLEPc (PETSc): "
              << std::flush;

      SlepcInitialize (&argc, &argv, 0, 0);
      {
        logfile << "ok"
                << std::endl;

        // Do something simple with PETSc
        logfile << "Using PetscScalar:"
                << std::endl;

        const PetscScalar pi_1i  = numbers::PI + 1.*PETSC_i;
        const PetscScalar two    = 2.;

        logfile << "   pi+1i:           " << pi_1i
                << std::endl
                << "   two:             " << two
                << std::endl
                << "   two times pi+1i: " << two *pi_1i
                << std::endl;


        logfile << "Finalizing SLEPc (PETSc): "
                << std::flush;

      }
      SlepcFinalize ();

      logfile << "ok"
              << std::endl << std::endl;
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


