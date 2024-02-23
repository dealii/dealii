// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2013 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// just check initialising SLEPc can be done and that it initialises
// PETSc in the way we expect, ie. *a* PETSc object exist.

#include <deal.II/base/numbers.h>

#include <deal.II/lac/slepc_solver.h>

#include <iostream>

#include "../tests.h"

int
main(int argc, char **argv)
{
  initlog();
  try
    {
      deallog.get_file_stream() << "Initializing SLEPc (PETSc): " << std::flush;

      SlepcInitialize(&argc, &argv, nullptr, nullptr);
      {
        deallog.get_file_stream() << "ok" << std::endl;

        // Do something simple with PETSc
        deallog.get_file_stream() << "Using PetscScalar:" << std::endl;

        const PetscScalar pi  = numbers::PI;
        const PetscScalar two = 2.;

        deallog.get_file_stream()
          << "   pi:           " << pi << std::endl
          << "   two:          " << two << std::endl
          << "   two times pi: " << two * pi << std::endl;


        deallog.get_file_stream() << "Finalizing SLEPc (PETSc): " << std::flush;
      }
      SlepcFinalize();

      deallog.get_file_stream() << "ok" << std::endl << std::endl;
    }

  catch (const std::exception &exc)
    {
      std::cerr << std::endl
                << std::endl
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
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Unknown exception!" << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    };
}
