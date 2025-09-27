// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2014 - 2021 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Test by Alexander Grayver: The FE_Nedelec did not correctly compute face
// interpolation matrices from lower to higher order elements. This
// consequently led to a situation where
// DoFTools::make_hanging_node_constraints got into trouble because it could
// not find a primary DoF for one particular dependent DoF.


char logname[] = "output";


#include <deal.II/base/function.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_nedelec.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/vector.h>

#include <deal.II/numerics/vector_tools.h>

#include <vector>

#include "../tests.h"



template <int dim>
void
test()
{
  // create a mesh like this (viewed
  // from top, if in 3d):
  // *---*---*
  // | 0 | 1 |
  // *---*---*
  Triangulation<dim>        triangulation;
  std::vector<unsigned int> subdivisions(dim, 1);
  subdivisions[0] = 2;
  GridGenerator::subdivided_hyper_rectangle(triangulation,
                                            subdivisions,
                                            Point<dim>(),
                                            (dim == 3 ? Point<dim>(2, 1, 1) :
                                                        Point<dim>(2, 1)));
  (std::next(triangulation.begin_active()))->set_refine_flag();
  triangulation.execute_coarsening_and_refinement();

  hp::FECollection<dim> fe;
  fe.push_back(FE_Nedelec<dim>(1));
  fe.push_back(FE_Nedelec<dim>(2));

  DoFHandler<dim> dof_handler(triangulation);

  for (unsigned int i = 0; i < fe.size(); ++i)
    for (unsigned int j = 0; j < fe.size(); ++j)
      {
        deallog << "Testing " << fe[i].get_name() << " vs. " << fe[j].get_name()
                << std::endl;

        // set FE on coarse cell to 'i', on
        // all fine cells to 'j'
        typename DoFHandler<dim>::active_cell_iterator cell =
          dof_handler.begin_active();
        cell->set_active_fe_index(i);
        ++cell;

        for (; cell != dof_handler.end(); ++cell)
          cell->set_active_fe_index(j);

        dof_handler.distribute_dofs(fe);

        AffineConstraints<double> constraints;
        DoFTools::make_hanging_node_constraints(dof_handler, constraints);
        constraints.close();

        constraints.print(deallog.get_file_stream());
      }
}



int
main()
{
  std::ofstream logfile(logname);
  logfile.precision(7);

  deallog.attach(logfile);


  test<2>();
  test<3>();
}
