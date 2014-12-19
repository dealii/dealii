// ---------------------------------------------------------------------
//
// Copyright (C) 2006 - 2013 by the deal.II authors
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


// Test FETools::get_fe_by_name

#include "../tests.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>

#include <deal.II/base/logstream.h>
#include <deal.II/fe/fe.h>
#include <deal.II/fe/fe_tools.h>
#include <deal.II/base/quadrature.h>

template <int dim>
void test_fe(const char *name)
{
  FiniteElement<dim> *fe = FETools::get_fe_from_name<dim>(std::string(name));

  deallog << fe->get_name() << std::endl
          << '\t' << fe->dofs_per_cell
          << '\t' << fe->dofs_per_vertex
          << '\t' << fe->dofs_per_line
          << '\t' << fe->dofs_per_quad
          << '\t' << fe->dofs_per_hex
          << std::endl;
}


int
main()
{
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.depth_console(0);

  // These names are all correct.
  test_fe<1>("FE_Q(1)");
  test_fe<1>("FE_Q(2)");
  test_fe<1>("FE_Q<1>(2)");
  test_fe<1>("FE_Q(QGaussLobatto(3))");
  test_fe<1>("FE_Q(QGaussLobatto(4))");
  test_fe<2>("FE_Q(1)");
  test_fe<2>("FE_Q<2>(2)");
  test_fe<2>("FE_Q<dim>(2)");
  test_fe<2>("FE_DGQ(1)");
  test_fe<2>("FE_RaviartThomas(1)");
  test_fe<2>("FE_Q(QGaussLobatto(3))");
  test_fe<2>("FE_Q(QGaussLobatto(4))");
  test_fe<3>("FE_Q(1)");
  test_fe<1>("FESystem<1>[FE_Q<dim>(2)^dim-FE_DGQ<d>(1)]");
  test_fe<2>("FESystem<2>[FE_Q<2>(2)^dim-FE_DGQ<2>(1)]");
  test_fe<2>("FESystem[FESystem<2>[FE_Q<2>(2)^2-FE_DGQ<2>(1)]^2-FE_Q(1)]");
  test_fe<2>("FESystem[FESystem[FESystem[FE_Q(1)^2-FE_Q(1)]^2]-FESystem[FE_Q(1)^2]-FESystem[FE_Q(1)-FE_DGP(0)]]");

  // Now set up a list of malformed
  // names
  std::vector<const char *> names;
//  names.push_back("FE_Q[2]");

  for (unsigned int i=0; i<names.size(); ++i)
    {
      try
        {
          test_fe<2>(names[i]);
        }
      catch (ExceptionBase &e)
        {
          logfile << e.what();
        }
    }
}
