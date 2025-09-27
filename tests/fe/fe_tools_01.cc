// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2006 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Test FETools::get_fe_by_name

#include <deal.II/base/quadrature.h>

#include <deal.II/fe/fe.h>
#include <deal.II/fe/fe_tools.h>

#include <iostream>
#include <sstream>

#include "../tests.h"

template <int dim>
void
test_fe(const char *name)
{
  std::unique_ptr<FiniteElement<dim>> fe =
    FETools::get_fe_by_name<dim, dim>(std::string(name));

  deallog << fe->get_name() << std::endl
          << '\t' << fe->dofs_per_cell << '\t' << fe->dofs_per_vertex << '\t'
          << fe->dofs_per_line << '\t' << fe->dofs_per_quad << '\t'
          << fe->dofs_per_hex << std::endl;
}


int
main()
{
  initlog();

  // These names are all correct.
  test_fe<1>("FE_Q(1)");
  test_fe<1>("FE_Q(2)");
  test_fe<1>("FE_Q<1>(2)");
  test_fe<1>("FE_Q(QIterated(QTrapezoid(),2))");
  test_fe<1>("FE_Q(QGaussLobatto(4))");
  test_fe<1>("FE_Q(3)");
  test_fe<1>("FE_Q(QIterated(QTrapezoid(),3))");
  test_fe<2>("FE_Q(1)");
  test_fe<2>("FE_Q<2>(2)");
  test_fe<2>("FE_Q<dim>(2)");
  test_fe<2>("FE_DGQ(1)");
  test_fe<2>("FE_Bernstein(2)");
  test_fe<2>("FE_RaviartThomas(1)");
  test_fe<2>("FE_Q(QGaussLobatto(3))");
  test_fe<2>("FE_Q(QGaussLobatto(4))");
  test_fe<2>("FE_Q(QIterated(QTrapezoid(),3))");
  test_fe<3>("FE_Q(1)");
  test_fe<1>("FESystem<1>[FE_Q<dim>(2)^dim-FE_DGQ<d>(1)]");
  test_fe<2>("FESystem<2>[FE_Q<2>(2)^dim-FE_DGQ<2>(1)]");
  test_fe<2>("FESystem[FESystem<2>[FE_Q<2>(2)^2-FE_DGQ<2>(1)]^2-FE_Q(1)]");
  test_fe<2>(
    "FESystem[FESystem[FESystem[FE_Q(1)^2-FE_Q(1)]^2]-FESystem[FE_Q(1)^2]-FESystem[FE_Q(1)-FE_DGP(0)]]");

  // Now set up a list of malformed
  // names
  std::vector<const char *> names;
  //  names.push_back("FE_Q[2]");

  for (unsigned int i = 0; i < names.size(); ++i)
    {
      try
        {
          test_fe<2>(names[i]);
        }
      catch (ExceptionBase &e)
        {
          deallog << e.what();
        }
    }
}
