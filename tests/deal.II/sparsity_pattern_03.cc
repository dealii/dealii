// ---------------------------------------------------------------------
//
// Copyright (C) 2000 - 2013 by the deal.II authors
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



// use the version of DoFTools::make_sparsity_pattern that takes two
// DoFHandler arguments for two DoFHandlers that are actually from different
// meshes (though with the same base)
//
// like sparsity_pattern_03, but use different finite elements


#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/lac/sparsity_pattern.h>
#include <deal.II/lac/block_sparsity_pattern.h>
#include <deal.II/lac/compressed_sparsity_pattern.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_system.h>

#include <fstream>




template <int dim>
void
check ()
{
  // create two different triangulations
  Triangulation<dim> triangulation_1;
  if (dim==2)
    GridGenerator::hyper_ball(triangulation_1, Point<dim>(), 1);
  else
    GridGenerator::hyper_cube(triangulation_1, -1,1);
  triangulation_1.refine_global (1);
  triangulation_1.begin_active()->set_refine_flag ();
  triangulation_1.execute_coarsening_and_refinement ();
  triangulation_1.begin_active(2)->set_refine_flag ();
  triangulation_1.execute_coarsening_and_refinement ();
  if (dim==1)
    triangulation_1.refine_global(2);


  Triangulation<dim> triangulation_2;
  if (dim==2)
    GridGenerator::hyper_ball(triangulation_2, Point<dim>(), 1);
  else
    GridGenerator::hyper_cube(triangulation_2, -1,1);
  triangulation_2.refine_global (1);
  (++triangulation_2.begin_active())->set_refine_flag ();
  triangulation_2.execute_coarsening_and_refinement ();
  (++triangulation_2.begin_active(2))->set_refine_flag ();
  triangulation_2.execute_coarsening_and_refinement ();
  if (dim==1)
    triangulation_2.refine_global(2);



  FESystem<dim> element_1(FE_Q<dim>(1), 2,
                          FE_Q<dim>(2), 1);
  FESystem<dim> element_2(FE_Q<dim>(3), 1,
                          FE_DGQ<dim>(0), 2);
  DoFHandler<dim> dof_1(triangulation_1);
  DoFHandler<dim> dof_2(triangulation_2);
  dof_1.distribute_dofs(element_1);
  dof_2.distribute_dofs(element_2);

  SparsityPattern sparsity (dof_1.n_dofs(), dof_2.n_dofs(),
                            std::max(dof_1.n_dofs(), dof_2.n_dofs()));
  DoFTools::make_sparsity_pattern (dof_1, dof_2, sparsity);
  sparsity.compress ();

  sparsity.print (deallog.get_file_stream());
}



int main ()
{
  std::ofstream logfile ("output");
  deallog << std::setprecision (2);
  deallog << std::fixed;
  deallog.attach(logfile);
  deallog.depth_console (0);

  deallog.push ("1d");
  check<1> ();
  deallog.pop ();
  deallog.push ("2d");
  check<2> ();
  deallog.pop ();
  deallog.push ("3d");
  check<3> ();
  deallog.pop ();
}
