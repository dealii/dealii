// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2022 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#include <deal.II/base/exceptions.h>
#include <deal.II/base/function_signed_distance.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_q_iso_q1.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/non_matching/quadrature_generator.h>

#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"

#include "quadrature_printing.h"


template <int dim>
class Test
{
public:
  Test(const FE_Poly<dim> &fe);

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
Test<dim>::Test(const FE_Poly<dim> &fe)
  : dof_handler(triangulation)
{
  const FE_Q_iso_Q1<dim> *fe_q_iso_q1 =
    dynamic_cast<const FE_Q_iso_Q1<dim> *>(&fe);
  const bool is_fe_q_iso_q1 = fe_q_iso_q1 != nullptr;
  fe_collection.push_back(fe);

  const unsigned int n_quadrature_points = is_fe_q_iso_q1 ? 1 : fe.get_degree();
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
  const Functions::SignedDistance::Plane<dim> analytical_levelset(
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

  deallog << "inside" << std::endl;
  print_quadrature(quadrature_generator.get_inside_quadrature());
  deallog << "outside" << std::endl;
  print_quadrature(quadrature_generator.get_outside_quadrature());
  deallog << "surface" << std::endl;
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
      deallog << "face_index = " << f << std::endl;
      deallog << "inside" << std::endl;
      print_quadrature(face_quadrature_generator.get_inside_quadrature());
      deallog << "outside" << std::endl;
      print_quadrature(face_quadrature_generator.get_outside_quadrature());
      deallog << "surface" << std::endl;
      print_surface_quadrature(
        face_quadrature_generator.get_surface_quadrature());
    }
}



template <int dim>
void
run_test(const FE_Poly<dim> &fe)
{
  deallog << "dim = " << dim << std::endl;
  Test<dim> test(fe);
  test.run();
  deallog << std::endl;
}



int
main()
{
  initlog();

  std::shared_ptr<FE_Poly<1>> fe_1;
  std::shared_ptr<FE_Poly<2>> fe_2;
  std::shared_ptr<FE_Poly<3>> fe_3;
  fe_1 = std::make_shared<FE_Q<1>>(1);
  fe_2 = std::make_shared<FE_Q<2>>(1);
  fe_3 = std::make_shared<FE_Q<3>>(1);
  run_test<1>(*fe_1);
  run_test<2>(*fe_2);
  run_test<3>(*fe_3);
  fe_1 = std::make_shared<FE_Q_iso_Q1<1>>(2);
  fe_2 = std::make_shared<FE_Q_iso_Q1<2>>(2);
  fe_3 = std::make_shared<FE_Q_iso_Q1<3>>(2);
  run_test<1>(*fe_1);
  run_test<2>(*fe_2);
  run_test<3>(*fe_3);
}
