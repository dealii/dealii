// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2018 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// check that we get a reasonable error message when trying to call
// MGConstrainedDofs::initialize() with distributing level dofs.

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/multigrid/mg_constrained_dofs.h>

#include <algorithm>

#include "../tests.h"

int
main(int argc, char *argv[])
{
  initlog();
  deal_II_exceptions::disable_abort_on_exception();

  constexpr int      dim = 2;
  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria);
  DoFHandler<dim> dh(tria);
  dh.distribute_dofs(FE_Q<dim>(1));
  MGConstrainedDoFs mg_dofs;
  try
    {
      mg_dofs.initialize(dh);
    }
  catch (ExceptionBase &e)
    {
      deallog << e.get_exc_name() << std::endl;
    }
}
