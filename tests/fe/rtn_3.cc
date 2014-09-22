// ---------------------------------------------------------------------
//
// Copyright (C) 2003 - 2013 by the deal.II authors
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


// output the number of dofs per cell for the RTN element; that seemed wrong
// in the past. compare it also with the number of dofs per cell of the RT
// element

#include "../tests.h"
#include <deal.II/fe/fe_raviart_thomas.h>

#include <vector>
#include <fstream>
#include <string>

#define PRECISION 2



template<int dim>
void test ()
{
  for (unsigned int degree=0; degree<9-2*dim; ++degree)
    {
      FE_RaviartThomasNodal<dim> fe_rtn(degree);
      deallog << fe_rtn.get_name() << ' ' << fe_rtn.dofs_per_cell << std::endl;

      if (dim != 3)
        {
          FE_RaviartThomas<dim> fe_rt(degree);
          deallog << fe_rt.get_name() << ' ' << fe_rt.dofs_per_cell << std::endl;
        }
    }
}


int
main()
{
  std::ofstream logfile ("output");
  deallog << std::setprecision(PRECISION);
  deallog << std::fixed;
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test<2>();
  test<3>();

  return 0;
}



