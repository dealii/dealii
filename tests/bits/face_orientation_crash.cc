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


// trip up the new code handling hanging node constraints with a face in 3d
// that has face_orientation==false



#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/base/function_lib.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/lac/vector.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/mapping_q.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/error_estimator.h>

#include <fstream>


template <int dim>
void
check ()
{
  // create a mesh with at least one cell
  // that has a face with
  // face_orientation==false. refine each of
  // the 7 cells in turn, to make sure we
  // have a face with hanging nodes that has
  // face_orientation==false at least once
  for (unsigned int i=0; i<7; ++i)
    {
      deallog << "Check " << i << std::endl;

      Triangulation<dim> tria;
      GridGenerator::hyper_ball(tria);

      typename Triangulation<dim>::active_cell_iterator
      cell = tria.begin_active();
      std::advance(cell,i);
      cell->set_refine_flag();
      tria.execute_coarsening_and_refinement ();

      // attach a DoFHandler
      FE_Q<dim> element(1);
      DoFHandler<dim> dof(tria);
      dof.distribute_dofs(element);

      // then build hanging node
      // constraints. this should trip the
      // new code using the hp constraints,
      // added in late July 2006
      ConstraintMatrix constraints;
      DoFTools::make_hanging_node_constraints (dof,
                                               constraints);

      for (unsigned int j=0; j<dof.n_dofs(); ++j)
        if (constraints.is_constrained (j))
          deallog << j << std::endl;
    }

  deallog << "OK" << std::endl;
}


int main ()
{
  std::ofstream logfile ("output");
  deallog << std::setprecision (2);
  deallog << std::fixed;
  deallog.attach(logfile);
  deallog.depth_console (0);

  check<3> ();
}
