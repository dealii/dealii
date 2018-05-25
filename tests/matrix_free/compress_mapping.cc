// ---------------------------------------------------------------------
//
// Copyright (C) 2013 - 2018 by the deal.II authors
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


// this function tests whether the compression of mapping (Jacobians) works
// properly. There should only be a few different Jacobians also when there
// are many cells as the weights should be identical

#include <deal.II/base/function.h>
#include <deal.II/base/utilities.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/constraint_matrix.h>

#include <deal.II/matrix_free/matrix_free.h>

#include "../tests.h"
#include "create_mesh.h"

std::ofstream logfile("output");



template <int dim>
void
test()
{
  deallog << "General mesh" << std::endl;
  Triangulation<dim> tria;
  create_mesh(tria);
  tria.begin_active()->set_refine_flag();
  tria.execute_coarsening_and_refinement();
  typename Triangulation<dim>::active_cell_iterator cell, endc;
  cell = tria.begin_active();
  endc = tria.end();
  for (; cell != endc; ++cell)
    if (cell->center().norm() < 0.5)
      cell->set_refine_flag();
  tria.execute_coarsening_and_refinement();
  tria.begin(tria.n_levels() - 1)->set_refine_flag();
  tria.execute_coarsening_and_refinement();
  for (unsigned int i = 0; i < 5 - dim; ++i)
    {
      cell                 = tria.begin_active();
      endc                 = tria.end();
      unsigned int counter = 0;
      for (; cell != endc; ++cell, ++counter)
        if (cell->center()[0] < 5.)
          cell->set_refine_flag();
      tria.execute_coarsening_and_refinement();
    }

  FE_Q<dim>       fe(1);
  DoFHandler<dim> dof(tria);
  dof.distribute_dofs(fe);
  ConstraintMatrix constraints;
  DoFTools::make_hanging_node_constraints(dof, constraints);
  constraints.close();

  const QGauss<1>                          quad(2);
  MatrixFree<dim>                          mf;
  typename MatrixFree<dim>::AdditionalData data;
  data.tasks_parallel_scheme = MatrixFree<dim>::AdditionalData::none;
  mf.reinit(dof, constraints, quad, data);

  const unsigned int        n_macro_cells = mf.n_macro_cells();
  std::vector<unsigned int> n_cell_types(4, 0);
  for (unsigned int i = 0; i < n_macro_cells; ++i)
    n_cell_types[mf.get_mapping_info().get_cell_type(i)]++;

  // should do at least some compression
  Assert(n_cell_types[0] + n_cell_types[1] > 0, ExcInternalError());
  Assert(mf.get_mapping_info().cell_data[0].jacobians[0].size() <
           (n_cell_types[3] * quad.size() + n_macro_cells - n_cell_types[3]),
         ExcInternalError());
  deallog << "OK" << std::endl;
}



template <int dim>
void
test_cube()
{
  deallog << "Hyper cube" << std::endl;
  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria);
  tria.refine_global(5 - dim);

  FE_Q<dim>       fe(1);
  DoFHandler<dim> dof(tria);
  dof.distribute_dofs(fe);
  ConstraintMatrix constraints;
  DoFTools::make_hanging_node_constraints(dof, constraints);
  constraints.close();

  const QGauss<1>                          quad(2);
  MatrixFree<dim>                          mf;
  typename MatrixFree<dim>::AdditionalData data;
  data.tasks_parallel_scheme = MatrixFree<dim>::AdditionalData::none;
  mf.reinit(dof, constraints, quad, data);

  const unsigned int        n_macro_cells = mf.n_macro_cells();
  std::vector<unsigned int> n_cell_types(4, 0);
  for (unsigned int i = 0; i < n_macro_cells; ++i)
    n_cell_types[mf.get_mapping_info().get_cell_type(i)]++;

  // should have one Cartesian cell and no other
  // cell type
  AssertDimension(n_cell_types[0], n_macro_cells);
  AssertDimension(mf.get_mapping_info().cell_data[0].jacobians[0].size(), 1);
  Assert(n_macro_cells > 1, ExcInternalError());
  deallog << "OK" << std::endl;
}



template <int dim>
void
test_parallelogram()
{
  deallog << "Parallelogram" << std::endl;
  Triangulation<dim> tria;
  Point<dim>         corners[dim];
  for (unsigned int d = 0; d < dim; ++d)
    {
      corners[d][d] = 1.;
      if (d > 0)
        corners[d][0] = 0.5 + 0.5 * d;
    }
  GridGenerator::parallelepiped(tria, corners);
  tria.refine_global(5 - dim);

  FE_Q<dim>       fe(1);
  DoFHandler<dim> dof(tria);
  dof.distribute_dofs(fe);
  ConstraintMatrix constraints;
  DoFTools::make_hanging_node_constraints(dof, constraints);
  constraints.close();

  const QGauss<1>                          quad(2);
  MatrixFree<dim>                          mf;
  typename MatrixFree<dim>::AdditionalData data;
  data.tasks_parallel_scheme = MatrixFree<dim>::AdditionalData::none;
  mf.reinit(dof, constraints, quad, data);

  const unsigned int        n_macro_cells = mf.n_macro_cells();
  std::vector<unsigned int> n_cell_types(4, 0);
  for (unsigned int i = 0; i < n_macro_cells; ++i)
    n_cell_types[mf.get_mapping_info().get_cell_type(i)]++;

  // should have one affine cell and no other
  // cell type
  AssertDimension(n_cell_types[1], n_macro_cells);
  AssertDimension(mf.get_mapping_info().cell_data[0].jacobians[0].size(), 1);
  Assert(n_macro_cells > 1, ExcInternalError());
  deallog << "OK" << std::endl;
}



int
main()
{
  deallog.attach(logfile);

  deallog << std::setprecision(3);

  {
    deallog.push("2d");
    test<2>();
    test_cube<2>();
    test_parallelogram<2>();
    deallog.pop();
    deallog.push("3d");
    test<3>();
    test_cube<3>();
    test_parallelogram<3>();
    deallog.pop();
  }
}
