// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2016 - 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Like the data_out_faces_postprocess_boundary_id test, except that
// it uses the class that is now actually part of the library, rather
// than its own implementation.


#include <deal.II/base/function.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_system.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/vector.h>

#include <deal.II/numerics/data_out_faces.h>
#include <deal.II/numerics/data_postprocessor.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"



std::ofstream logfile("output");


template <int dim>
void
test()
{
  Triangulation<dim> triangulation;
  FE_DGQ<dim>        fe(1);
  DoFHandler<dim>    dof_handler(triangulation);

  GridGenerator::hyper_cube(triangulation, 0, 1, true);
  triangulation.refine_global(1);
  triangulation.begin_active()->set_refine_flag();
  triangulation.execute_coarsening_and_refinement();

  dof_handler.distribute_dofs(fe);

  // Create a dummy vector. We will ignore its contents.
  Vector<double> solution(dof_handler.n_dofs());

  DataPostprocessors::BoundaryIds<dim> p;
  DataOutFaces<dim>                    data_out;
  data_out.attach_dof_handler(dof_handler);
  data_out.add_data_vector(solution, p);
  data_out.build_patches();
  data_out.write_gnuplot(logfile);
}


int
main()
{
  logfile << std::setprecision(2);
  deallog << std::setprecision(2);

  test<2>();
  test<3>();

  return 0;
}
