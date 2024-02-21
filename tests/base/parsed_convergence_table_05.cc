// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2019 - 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Test the functionality of the ParsedConvergenceTable class for
// vector function + scalar function

#include <deal.II/base/function_lib.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/parsed_convergence_table.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/vector.h>

#include <deal.II/numerics/vector_tools.h>

#include <map>

#include "../tests.h"

int
main()
{
  initlog();

  ParsedConvergenceTable table({"u", "u", "p"},
                               {{VectorTools::H1_norm},
                                {VectorTools::L2_norm}});

  ParameterHandler prm;
  table.add_parameters(prm);

  Triangulation<2> tria;
  GridGenerator::hyper_cube(tria);

  FESystem<2>   fe(FE_Q<2>(1), 3);
  DoFHandler<2> dh(tria);

  Functions::CosineFunction<2> exact(3);

  for (unsigned int i = 0; i < 5; ++i)
    {
      tria.refine_global(1);
      dh.distribute_dofs(fe);
      Vector<double> sol(dh.n_dofs());
      VectorTools::interpolate(dh, exact, sol);
      table.error_from_exact(dh, sol, exact);
    }
  table.output_table(deallog.get_file_stream());
}
