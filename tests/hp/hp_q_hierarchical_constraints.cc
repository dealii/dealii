// ---------------------------------------------------------------------
//
// Copyright (C) 2003 - 2020 by the deal.II authors
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


#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q_hierarchical.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/affine_constraints.h>

#include <deal.II/numerics/data_out.h>

#include <iostream>

#include "../tests.h"

/* A test to check that the resulting FE space is continuous
 * for different hp-refinement cases (p-only, "simple", "complex").
 * This was confirmed by visual examination (see #define below).
 * The correct constrains are written in a log to make a unit test.
 */

//#define FEQH_DEBUG_OUTPUT


template <int dim>
void
test(const bool apply_constrains, const unsigned int hp)
{
  Triangulation<dim> triangulation;
  {
    Triangulation<dim> triangulationL;
    Triangulation<dim> triangulationR;
    GridGenerator::hyper_cube(triangulationL,
                              -1,
                              0); // create a square [-1,0]^d domain
    GridGenerator::hyper_cube(triangulationR,
                              -1,
                              0); // create a square [-1,0]^d domain
    Point<dim> shift_vector;
    shift_vector[0] = 1.0;
    GridTools::shift(shift_vector, triangulationR);
    GridGenerator::merge_triangulations(triangulationL,
                                        triangulationR,
                                        triangulation);
  }

  hp::FECollection<dim>     fe;
  hp::DoFHandler<dim>       dof_handler(triangulation);
  AffineConstraints<double> constraints; // for boundary conditions


  // populate fe system:
  fe.push_back(FE_Q_Hierarchical<dim>(2));
  fe.push_back(FE_Q_Hierarchical<dim>(4));

  // set one cell to have different active_fe_index:
  typename hp::DoFHandler<dim>::active_cell_iterator cell =
    dof_handler.begin_active();
  cell->set_active_fe_index(1);

  // need to distribute dofs before refinement,
  // otherwise active fe index will not transfer to child cells
  dof_handler.distribute_dofs(fe);

  // refine first cell (simple case)
  std::string hp_string;
  if (hp == 1)
    {
      hp_string = "_hpSimple";
      cell->set_refine_flag();
    }
  // refine second cell
  else if (hp == 2)
    {
      hp_string = "_hpComplex";
      cell++;
      cell->set_refine_flag();
    }

  triangulation.execute_coarsening_and_refinement();

  dof_handler.distribute_dofs(fe);
  constraints.clear();
  DoFTools::make_hanging_node_constraints(dof_handler, constraints);
  constraints.close();

  const unsigned int n_dofs = dof_handler.n_dofs();

  Vector<double> v(n_dofs);

  if (apply_constrains)
    deallog << hp_string << std::endl;

  for (unsigned int i = 0; i < n_dofs; i++)
    {
      v    = 0.;
      v[i] = 1.;
      if (apply_constrains)
        {
          constraints.distribute(v);
          deallog << "i=" << i << std::endl;
          constraints.print(deallog.get_file_stream());
        }

#ifdef FEQH_DEBUG_OUTPUT
      DataOut<dim, hp::DoFHandler<dim>> data_out;
      data_out.attach_dof_handler(dof_handler);

      data_out.add_data_vector(v, "shape_function");

      // do so rather big number of subdivision to better see shape functions
      data_out.build_patches(20);

      std::ostringstream filename;
      filename << "shape_" << dim << "d"
               << (apply_constrains ? "_constrained" : "") << hp_string << "_"
               << i << ".vtk";

      std::ofstream output(filename.str().c_str());
      data_out.write_vtk(output);
#endif
    }
}

int
main(int argc, char *argv[])
{
  initlog();

  test<2>(true, 0);
  test<2>(true, 1);
  test<2>(true, 2);

#ifdef FEQH_DEBUG_OUTPUT
  test<2>(false, 0);
  test<2>(false, 1);
  test<2>(false, 2);
#endif
}
