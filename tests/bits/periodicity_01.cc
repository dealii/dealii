// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2013 - 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// check periodic boundary conditions for a simple enough case where we know
// the exact set of constraints
//
// this test simply uses two hypercubes and matches the faces at the far ends

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include "../tests.h"



template <int dim>
void
test()
{
  // create a 2x1 (or 2x1x1) mesh
  Triangulation<dim>        triangulation;
  std::vector<unsigned int> repetitions(dim, 1);
  repetitions[0] = 2;
  GridGenerator::subdivided_hyper_rectangle(
    triangulation,
    repetitions,
    Point<dim>(),
    (dim == 1 ? Point<dim>(2) :
                (dim == 2 ? Point<dim>(2, 1) : Point<dim>(2, 1, 1))));

  FE_Q<dim>       fe(1);
  DoFHandler<dim> dof_handler(triangulation);
  dof_handler.distribute_dofs(fe);

  AffineConstraints<double> cm;
  DoFTools::make_periodicity_constraints(
    dof_handler.begin(0)->face(0),
    (std::next(dof_handler.begin(0)))->face(1),
    cm);

  deallog << "dim " << std::to_string(dim) << ':' << std::endl;
  cm.print(deallog.get_file_stream());
}



int
main()
{
  initlog();
  test<1>();
  test<2>();
  test<3>();
  return 0;
}
