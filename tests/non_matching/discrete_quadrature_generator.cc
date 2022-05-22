// ---------------------------------------------------------------------
//
// Copyright (C) 2021 - 2022 by the deal.II authors
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

#include <deal.II/base/exceptions.h>
#include <deal.II/base/function_level_set.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/non_matching/quadrature_generator.h>

#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"

#include "quadrature_printing.h"

using namespace dealii;

template <int dim>
class Test
{
public:
  Test();

  void
  run();

private:
  void
  setup_mesh();

  // Setup a discrete level set function corresponding to
  // $\psi(x) = (x_0 - 1.5) = 0$
  void
  setup_discrete_level_set();

  void
  test_discrete_quadrature_generator();

  void
  test_discrete_face_quadrature_generator();

  Triangulation<dim>    triangulation;
  hp::FECollection<dim> fe_collection;
  DoFHandler<dim>       dof_handler;

  hp::QCollection<1> q_collection1D;

  Vector<double> level_set;
};



template <int dim>
Test<dim>::Test()
  : dof_handler(triangulation)
{
  fe_collection.push_back(FE_Q<dim>(1));
  const unsigned int n_quadrature_points = 1;
  q_collection1D.push_back(QGauss<1>(n_quadrature_points));
}



template <int dim>
void
Test<dim>::run()
{
  setup_mesh();
  dof_handler.distribute_dofs(fe_collection);
  setup_discrete_level_set();
  test_discrete_quadrature_generator();
  test_discrete_face_quadrature_generator();
}



template <int dim>
void
Test<dim>::setup_mesh()
{
  const Point<dim> lower_left;
  Point<dim>       upper_right;
  upper_right[0] = 2.;

  for (unsigned int d = 1; d < dim; ++d)
    {
      upper_right[d] = 1.;
    }

  GridGenerator::hyper_rectangle(triangulation, lower_left, upper_right);
}



template <int dim>
void
Test<dim>::setup_discrete_level_set()
{
  Point<dim> point_on_zero_contour;
  point_on_zero_contour[0] = 1.5;
  const Functions::LevelSet::Plane<dim> analytical_levelset(
    point_on_zero_contour, Point<dim>::unit_vector(0));

  level_set.reinit(dof_handler.n_dofs());
  VectorTools::interpolate(dof_handler, analytical_levelset, level_set);
}



template <int dim>
void
Test<dim>::test_discrete_quadrature_generator()
{
  NonMatching::DiscreteQuadratureGenerator<dim> quadrature_generator(
    q_collection1D, dof_handler, level_set);
  quadrature_generator.generate(triangulation.begin_active());

  print_quadrature(quadrature_generator.get_inside_quadrature());
  print_quadrature(quadrature_generator.get_outside_quadrature());
  print_surface_quadrature(quadrature_generator.get_surface_quadrature());
}



template <int dim>
void
Test<dim>::test_discrete_face_quadrature_generator()
{
  NonMatching::DiscreteFaceQuadratureGenerator<dim> face_quadrature_generator(
    q_collection1D, dof_handler, level_set);

  const auto &cell = triangulation.begin_active();
  for (const auto f : cell->face_indices())
    {
      face_quadrature_generator.generate(cell, f);

      print_quadrature(face_quadrature_generator.get_inside_quadrature());
      print_quadrature(face_quadrature_generator.get_outside_quadrature());
      print_surface_quadrature(
        face_quadrature_generator.get_surface_quadrature());
    }
}



template <int dim>
void
run_test()
{
  deallog << "dim = " << dim << std::endl;
  Test<dim> test;
  test.run();
  deallog << std::endl;
}



int
main()
{
  initlog();

  run_test<1>();
  run_test<2>();
  run_test<3>();
}
