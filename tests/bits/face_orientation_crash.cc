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


// trip up the new code handling hanging node constraints with a face in 3d
// that has face_orientation==false



#include <deal.II/base/function_lib.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/vector.h>

#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"



template <int dim>
void
check()
{
  // create a mesh with at least one cell
  // that has a face with
  // face_orientation==false. refine each of
  // the 7 cells in turn, to make sure we
  // have a face with hanging nodes that has
  // face_orientation==false at least once
  for (unsigned int i = 0; i < 7; ++i)
    {
      deallog << "Check " << i << std::endl;

      Triangulation<dim> tria;
      GridGenerator::hyper_ball(tria);

      typename Triangulation<dim>::active_cell_iterator cell =
        tria.begin_active();
      std::advance(cell, i);
      cell->set_refine_flag();
      tria.execute_coarsening_and_refinement();

      // attach a DoFHandler
      FE_Q<dim>       element(1);
      DoFHandler<dim> dof(tria);
      dof.distribute_dofs(element);

      // then build hanging node
      // constraints. this should trip the
      // new code using the hp-constraints,
      // added in late July 2006
      AffineConstraints<double> constraints;
      DoFTools::make_hanging_node_constraints(dof, constraints);

      for (unsigned int j = 0; j < dof.n_dofs(); ++j)
        if (constraints.is_constrained(j))
          deallog << j << std::endl;
    }

  deallog << "OK" << std::endl;
}


int
main()
{
  initlog();
  deallog << std::setprecision(2) << std::fixed;

  check<3>();
}
