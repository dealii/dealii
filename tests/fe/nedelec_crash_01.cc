// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2010 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// A test extracted from integrators/cochain_01 that crashed with the
// FE_Nedelec at the time

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_nedelec.h>

#include "../tests.h"

#include "../test_grids.h"


int
main()
{
  const std::string logname = "output";
  std::ofstream     logfile(logname);
  deallog.attach(logfile);

  // generate a version of the
  // Nedelec element but force it to
  // use the old-style constraints
  struct MyFE : FE_Nedelec<3>
  {
    MyFE()
      : FE_Nedelec<3>(0){};
    virtual bool
    hp_constraints_are_implemented() const
    {
      return false;
    }
  } fe;

  Triangulation<3> tr(Triangulation<3>::limit_level_difference_at_vertices);
  TestGrids::hypercube(tr, 2, true);

  DoFHandler<3> dof(tr);
  dof.distribute_dofs(fe);

  AffineConstraints<double> constraints;
  DoFTools::make_hanging_node_constraints(dof, constraints);
  constraints.close();

  constraints.print(deallog.get_file_stream());
}
