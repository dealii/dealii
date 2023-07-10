// ---------------------------------------------------------------------
//
// Copyright (C) 2008 - 2021 by the deal.II authors
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



// Test distributed solution transfer with fullydistributed triangulations.

#include <deal.II/distributed/solution_transfer.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_description.h>

#include <deal.II/lac/vector.h>

#include <deal.II/numerics/vector_tools.h>

#include "./tests.h"

using namespace dealii;

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
    return p.norm();
  }
};

template <int dim, typename TriangulationType>
void
test(TriangulationType &triangulation)
{
  DoFHandler<dim> dof_handler(triangulation);
  dof_handler.distribute_dofs(FE_Q<dim>(2));

  using VectorType = Vector<double>;

  VectorType vector(dof_handler.n_dofs());

  VectorTools::interpolate(dof_handler, InterpolationFunction<dim>(), vector);

  VectorType vector_loaded(dof_handler.n_dofs());

  const std::string filename =
    "solution_transfer_" + std::to_string(dim) + "d_out";

  {
    parallel::distributed::SolutionTransfer<dim, VectorType> solution_transfer(
      dof_handler);
    solution_transfer.prepare_for_serialization(vector);

    triangulation.save(filename);
  }

  triangulation.clear();

  {
    triangulation.load(filename);
    dof_handler.distribute_dofs(FE_Q<dim>(2));

    parallel::distributed::SolutionTransfer<dim, VectorType> solution_transfer(
      dof_handler);
    solution_transfer.deserialize(vector_loaded);
  }

  // Verify that error is 0.
  VectorType error(vector);
  error.add(-1, vector_loaded);

  deallog << (error.linfty_norm() < 1e-16 ? "PASSED" : "FAILED") << std::endl;
}


int
main(int argc, char **argv)
{
  initlog();

  deallog.push("2d");
  {
    constexpr int dim = 2;

    Triangulation<dim> triangulation;
    GridGenerator::hyper_cube(triangulation);
    triangulation.refine_global(3);

    test<dim>(triangulation);
  }
  deallog.pop();

  deallog.push("3d");
  {
    constexpr int dim = 3;

    Triangulation<dim> triangulation;
    GridGenerator::hyper_cube(triangulation);
    triangulation.refine_global(3);

    test<dim>(triangulation);
  }
  deallog.pop();
}
