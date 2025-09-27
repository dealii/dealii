// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2005 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// output the number of dofs per cell for the RTN element; that seemed wrong
// in the past. compare it also with the number of dofs per cell of the RT
// element

#include <deal.II/fe/fe_raviart_thomas.h>

#include <string>
#include <vector>

#include "../tests.h"

#define PRECISION 8



template <int dim>
void
test()
{
  for (unsigned int degree = 0; degree < 9 - 2 * dim; ++degree)
    {
      FE_RaviartThomasNodal<dim> fe_rtn(degree);
      deallog << fe_rtn.get_name() << ' ' << fe_rtn.dofs_per_cell << std::endl;

      if (dim != 3)
        {
          FE_RaviartThomas<dim> fe_rt(degree);
          deallog << fe_rt.get_name() << ' ' << fe_rt.dofs_per_cell
                  << std::endl;
        }
    }
}


int
main()
{
  std::ofstream logfile("output");
  deallog << std::setprecision(PRECISION);
  deallog << std::fixed;
  deallog.attach(logfile);

  test<2>();
  test<3>();

  return 0;
}
