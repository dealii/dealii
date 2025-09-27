// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2023 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// Test serialization with triangulations.

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/vector.h>

#include <deal.II/numerics/vector_tools.h>

#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>

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

  print_statistics(triangulation, false);

  std::stringstream stream;
  {
    boost::archive::text_oarchive oa(stream);
    oa << triangulation;
    oa << vector;
  }

  triangulation.clear();

  {
    boost::archive::text_iarchive ia(stream);
    ia >> triangulation;
    ia >> vector_loaded;
  }
  print_statistics(triangulation, false);

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
