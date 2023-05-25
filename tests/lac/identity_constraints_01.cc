// ---------------------------------------------------------------------
//
// Copyright (C) 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------


// Check for identity constraints using different finite elements
// in h-refined grids with hanging nodes.
//   +---+---+-------+
//   |   |   |       |
//   +---+---+       |
//   |   |   |       |
//   +---+---+-------+


#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>

#include <deal.II/grid/tria.h>

#include <deal.II/lac/affine_constraints.h>

#include "../tests.h"

#include "../test_grids.h"


template <int dim>
void
test(const FiniteElement<dim> &fe)
{
  Triangulation<dim> tria;
  TestGrids::hyper_line(tria, 2);
  tria.begin_active()->set_refine_flag();
  tria.execute_coarsening_and_refinement();

  DoFHandler<dim> dofh(tria);
  dofh.distribute_dofs(fe);

  AffineConstraints<double> constraints;
  DoFTools::make_hanging_node_constraints(dofh, constraints);
  constraints.close();

  deallog << fe.get_name() << ":";
  // print dof identities
  for (const auto &line : constraints.get_lines())
    if ((line.entries.size() == 1) && (line.entries[0].second == 1.))
      deallog << " " << line.index << "=" << line.entries[0].first;
  deallog << std::endl;

#if false
  constraints.print(deallog.get_file_stream());
#endif
}


int
main()
{
  initlog();

  {
    constexpr int dim = 2;

    for (unsigned int degree = 1; degree <= 9; ++degree)
      test<dim>(FE_Q<dim>(degree));
  }

  {
    constexpr int dim = 3;

    for (unsigned int degree = 1; degree <= 9; ++degree)
      test<dim>(FE_Q<dim>(degree));
  }
}
