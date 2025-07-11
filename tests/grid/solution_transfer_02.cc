// ---------------------------------------------------------------------
//
// Copyright (C) 2008 - 2024 by the deal.II authors
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



// Test SolutionTransfer::interpolate() for serial triangulations.

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_description.h>

#include <deal.II/lac/vector.h>

#include <deal.II/numerics/solution_transfer.h>
#include <deal.II/numerics/vector_tools.h>

#include "./tests.h"


template <int dim>
class InterpolationFunction : public Function<dim>
{
public:
  InterpolationFunction()
    : Function<dim>(1)
  {}

  virtual double
  value(const Point<dim> &p, const unsigned int component = 0) const
  {
    return p[0];
  }
};

template <int dim>
void
test(const unsigned int type)
{
  Triangulation<dim> triangulation;
  GridGenerator::hyper_cube(triangulation);
  triangulation.refine_global(3);

  const FE_Q<dim> fe(2);

  DoFHandler<dim> dof_handler(triangulation);
  dof_handler.distribute_dofs(fe);

  using VectorType = Vector<double>;

  VectorType vector(dof_handler.n_dofs());
  VectorTools::interpolate(dof_handler, InterpolationFunction<dim>(), vector);

  SolutionTransfer<dim, VectorType> solution_transfer(dof_handler);

  triangulation.prepare_coarsening_and_refinement();
  solution_transfer.prepare_for_coarsening_and_refinement(vector);

  if (type == 0)
    {
      triangulation.refine_global(1);
    }
  else if (type == 1)
    {
      for (const auto &cell : triangulation.active_cell_iterators())
        if (cell->center()[0] < 0.5)
          cell->set_refine_flag();
      triangulation.execute_coarsening_and_refinement();
    }
  else if (type == 2)
    {
      for (const auto &cell : triangulation.active_cell_iterators())
        if (cell->center()[0] < 0.5)
          cell->set_coarsen_flag();
      triangulation.execute_coarsening_and_refinement();
    }

  dof_handler.distribute_dofs(fe);

  vector.reinit(dof_handler.n_dofs());

  solution_transfer.interpolate(vector);

  VectorType error(dof_handler.n_dofs());
  VectorTools::interpolate(dof_handler, InterpolationFunction<dim>(), error);

  error -= vector;

  deallog << (error.linfty_norm() < 1e-16 ? "PASSED" : "FAILED") << std::endl;
}


int
main(int argc, char **argv)
{
  initlog();

  deallog.push("2d");
  {
    constexpr int dim = 2;

    test<dim>(0);
    test<dim>(1);
    test<dim>(2);
  }
  deallog.pop();

  deallog.push("3d");
  {
    constexpr int dim = 3;

    test<dim>(0);
    test<dim>(1);
    test<dim>(2);
  }
  deallog.pop();
}
