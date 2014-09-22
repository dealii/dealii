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



// Just output the restriction matrices of the BDM element

#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/fe/fe_bdm.h>

#include <fstream>
#include <string>

#define PRECISION 8



template<int dim>
void
test(const unsigned int degree)
{
  deallog << "FE_BDM<" << dim << "> (" << degree << ")"
          << std::endl;

  FE_BDM<dim> fe_rt(degree);

  for (unsigned int c=0; c<GeometryInfo<dim>::max_children_per_cell; ++c)
    {
      const FullMatrix<double> &m = fe_rt.get_restriction_matrix(c);

      for (unsigned int i=0; i<m.m(); ++i)
        {
          for (unsigned int j=0; j<m.n(); ++j)
            deallog << m(i,j) << ' ';
          deallog << std::endl;
        }

      deallog << std::endl;
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

  for (unsigned int degree=1; degree<4; ++degree)
    {
      test<2>(degree);
//      test<3>(degree);
    }

  return 0;
}
