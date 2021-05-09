// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2021 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------

#include <deal.II/fe/fe_abf.h>
#include <deal.II/fe/fe_tools.h>

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

      // Test only up to order 2 since the FE constructor will throw an
      // exception for higher order
      for (unsigned int fe_degree = 0; fe_degree < 2; ++fe_degree)
        {
          // H(div) conformal
          FE_ABF<dim> fe(fe_degree);

          {
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
          }
        } // ++fe_degree
    }
  catch (std::exception &exc)
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
