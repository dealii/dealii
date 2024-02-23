// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2020 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// Verify functionality of
// internal::hp::DoFHandlerImplementation::dominated_future_fe_on_children()

#include <deal.II/base/logstream.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/hp/fe_collection.h>

#include "../tests.h"



template <int dim>
void
test()
{
  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria);
  tria.refine_global(1);

  hp::FECollection<dim> fes;
  fes.push_back(FE_Q<dim>(1));
  fes.push_back(FE_Q<dim>(2));

  DoFHandler<dim> dofh(tria);
  dofh.begin_active()->set_active_fe_index(1);
  dofh.distribute_dofs(fes);

  const auto        &parent = dofh.begin(/*level=*/0);
  const unsigned int parent_future_fe =
    internal::hp::DoFHandlerImplementation::dominated_future_fe_on_children<
      dim>(parent);

  deallog << parent_future_fe << std::endl;
}


int
main()
{
  initlog();

  deallog.push("1d");
  test<1>();
  deallog.pop();
  deallog.push("2d");
  test<2>();
  deallog.pop();
  deallog.push("3d");
  test<3>();
  deallog.pop();
}
