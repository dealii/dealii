// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2004 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#include <deal.II/fe/fe_raviart_thomas.h>

#include "../tests.h"

// STL
#include <fstream>
#include <iostream>

// My test headers
#include "fe_conformity_test.h"

#define PRECISION 4

int
main(int, char **)
{
  std::ofstream logfile("output");
  dealii::deallog << std::setprecision(PRECISION);
  dealii::deallog << std::fixed;
  logfile << std::setprecision(PRECISION);
  logfile << std::fixed;
  dealii::deallog.attach(logfile);

  try
    {
      using namespace FEConforimityTest;
      constexpr int dim = 2;

      for (unsigned int fe_degree = 0; fe_degree < 5; ++fe_degree)
        {
          // H(div) conformal
          FE_RaviartThomasNodal<dim> fe(fe_degree);

          for (unsigned int this_switch = 0; this_switch < (dim == 2 ? 4 : 8);
               ++this_switch)
            {
              deallog << std::endl
                      << "*******   degree " << fe_degree
                      << "   *******   orientation case " << this_switch
                      << "   *******" << std::endl;

              FEConformityTest<dim> fe_conformity_tester(fe, this_switch);
              fe_conformity_tester.run();
            }
        } // ++fe_degree
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
    }

  return 0;
}
