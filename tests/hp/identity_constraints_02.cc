// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2021 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// Create a 2x2x1 grid with FE_Q elements with degrees 1 - 4 and 2 - 5 assigned
// in two respective scenarios.
// Verify that we do not unify dofs on lines in 3D.
// On each of the four lines on the interface between the Q2 and Q4 element, the
// central dofs are identical and will be treated with constraints.
//
// We put special emphasis on the central line.
// - In scenario 1, the central dof of the central line will be constrained
//   against the vertex dofs of the Q1 element. Thus, we only have 3 identity
//   constraints in total on the remaining lines of the interface between the Q2
//   and Q4 element.
// - In scenario 2, the central dof of the central line belongs to the Q2
//   element and remains unconstrained. Thus, we end up with 4 identity
//   constraints.
//
// Scenario 1:    Scenario 2:
// +----+----+    +----+----+
// | Q3 | Q4 |    | Q4 | Q5 |
// +----+----+    +----+----+
// | Q1 | Q2 |    | Q2 | Q3 |
// +----+----+    +----+----+

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/hp/fe_collection.h>

#include <deal.II/lac/affine_constraints.h>

#include <array>

#include "../tests.h"


template <int dim>
void
test(std::array<unsigned int, 4> fe_degrees)
{
  Triangulation<dim> tria;
  {
    std::vector<unsigned int> repetitions(dim);
    Point<dim>                bottom_left, top_right;
    for (unsigned int d = 0; d < dim; ++d)
      if (d < 2)
        {
          repetitions[d] = 2;
          bottom_left[d] = -1.;
          top_right[d]   = 1.;
        }
      else
        {
          repetitions[d] = 1;
          bottom_left[d] = 0.;
          top_right[d]   = 1.;
        }
    GridGenerator::subdivided_hyper_rectangle(tria,
                                              repetitions,
                                              bottom_left,
                                              top_right);
  }

  hp::FECollection<dim> fe;
  for (const auto d : fe_degrees)
    fe.push_back(FE_Q<dim>(d));

  DoFHandler<dim> dh(tria);
  {
    unsigned int i = 0;
    for (const auto &cell : dh.active_cell_iterators())
      cell->set_active_fe_index(i++);
    Assert(i == 4, ExcInternalError());
  }
  dh.distribute_dofs(fe);

  AffineConstraints<double> constraints;
  DoFTools::make_hanging_node_constraints(dh, constraints);
  constraints.close();

  deallog << "Total constraints:          " << constraints.n_constraints()
          << std::endl
          << "  Inhomogeneous constraints: " << constraints.n_inhomogeneities()
          << std::endl
          << "  Identity constraints:     " << constraints.n_identities()
          << std::endl;
}


int
main()
{
  initlog();

  deallog << "FE degrees: 1 - 4" << std::endl;
  deallog.push("2d");
  test<2>({{1, 2, 3, 4}});
  deallog.pop();
  deallog.push("3d");
  test<3>({{1, 2, 3, 4}});
  deallog.pop();

  deallog << "FE degrees: 2 - 5" << std::endl;
  deallog.push("2d");
  test<2>({{2, 3, 4, 5}});
  deallog.pop();
  deallog.push("3d");
  test<3>({{2, 3, 4, 5}});
  deallog.pop();
}
