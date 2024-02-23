// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2018 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Make sure that periodic boundary conditions also work correctly if we have
// multiple periodic boundary pairs that meet at an edge.
// This tests the 2D case.


#include <deal.II/base/function.h>
#include <deal.II/base/quadrature.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/vector.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>

#include <fstream>
#include <iostream>

#include "../tests.h"


AffineConstraints<double>
make_constraint_matrix(const DoFHandler<2> &dof_handler, int version)
{
  constexpr int             dim = 2;
  AffineConstraints<double> constraints;
  constraints.clear();
  DoFTools::make_hanging_node_constraints(dof_handler, constraints);

  std::vector<
    GridTools::PeriodicFacePair<typename DoFHandler<dim>::cell_iterator>>
    periodicity_vectorDof;
  switch (version)
    {
      case 0:
        GridTools::collect_periodic_faces(
          dof_handler, 0, 1, 0, periodicity_vectorDof);
        GridTools::collect_periodic_faces(
          dof_handler, 2, 3, 1, periodicity_vectorDof);
      case 1:
        GridTools::collect_periodic_faces(
          dof_handler, 0, 1, 0, periodicity_vectorDof);
        GridTools::collect_periodic_faces(
          dof_handler, 3, 2, 1, periodicity_vectorDof);
      case 2:
        GridTools::collect_periodic_faces(
          dof_handler, 1, 0, 0, periodicity_vectorDof);
        GridTools::collect_periodic_faces(
          dof_handler, 2, 3, 1, periodicity_vectorDof);
      case 3:
        GridTools::collect_periodic_faces(
          dof_handler, 1, 0, 0, periodicity_vectorDof);
        GridTools::collect_periodic_faces(
          dof_handler, 3, 2, 1, periodicity_vectorDof);
      case 4:
        GridTools::collect_periodic_faces(
          dof_handler, 2, 3, 1, periodicity_vectorDof);
        GridTools::collect_periodic_faces(
          dof_handler, 0, 1, 0, periodicity_vectorDof);
      case 5:
        GridTools::collect_periodic_faces(
          dof_handler, 3, 2, 1, periodicity_vectorDof);
        GridTools::collect_periodic_faces(
          dof_handler, 0, 1, 0, periodicity_vectorDof);
      case 6:
        GridTools::collect_periodic_faces(
          dof_handler, 2, 3, 1, periodicity_vectorDof);
        GridTools::collect_periodic_faces(
          dof_handler, 1, 0, 0, periodicity_vectorDof);
      case 7:
        GridTools::collect_periodic_faces(
          dof_handler, 3, 2, 1, periodicity_vectorDof);
        GridTools::collect_periodic_faces(
          dof_handler, 1, 0, 0, periodicity_vectorDof);
    }

  DoFTools::make_periodicity_constraints<dim, dim>(periodicity_vectorDof,
                                                   constraints);

  constraints.close();
  std::map<types::global_dof_index, Point<dim>> support_points;
  DoFTools::map_dofs_to_support_points(MappingQ<dim, dim>(1),
                                       dof_handler,
                                       support_points);
  for (const auto &line : constraints.get_lines())
    for (const auto &entry : line.entries)
      deallog << "DoF " << line.index << " at " << support_points[line.index]
              << " is constrained to "
              << " DoF " << entry.first << " at " << support_points[entry.first]
              << " with value " << entry.second << std::endl;
  return constraints;
}

template <int dim>
class PeriodicReference : public Function<dim>
{
public:
  PeriodicReference()
    : Function<dim>()
  {}
  virtual double
  value(const Point<dim> &p, const unsigned int component = 0) const override
  {
    if (dim == 3)
      return std::sin(p[0] + 1.) * std::sin(p[1] + 2.) * std::sin(p[2] + 3.);
    return std::sin(p[0] + 1.) * std::sin(p[1] + 2.);
  }
};


template <int dim>
void
get_point_value(const DoFHandler<dim> &dof_handler,
                const Point<dim>      &point,
                const Vector<double>  &solution,
                Vector<double>        &value)
{
  VectorTools::point_value(dof_handler, solution, point, value);
}


void
check_periodicity(const DoFHandler<2> &dof_handler,
                  Vector<double>      &solution,
                  const unsigned int   cycle)
{
  unsigned int n_points = 2;
  for (unsigned int i = 0; i < cycle; ++i)
    n_points *= 2;

  // don't test exactly at the support points, since point_value is not stable
  // there
  const double eps = 1. / (16. * n_points);

  for (unsigned int i = 1; i < n_points; ++i)
    {
      Vector<double> value1(1);
      Vector<double> value2(1);

      Point<2> point1;
      point1[0] = -numbers::PI + 2. * i / n_points + eps;
      point1[1] = -numbers::PI;
      Point<2> point2;
      point2[0] = -numbers::PI + 2. * i / n_points + eps;
      point2[1] = numbers::PI;

      VectorTools::point_value(dof_handler, solution, point1, value1);
      VectorTools::point_value(dof_handler, solution, point2, value2);

      if (std::abs(value2[0] - value1[0]) > 1e-8)
        {
          std::cout << point1 << "\t"
                    << "fail" << std::endl;
          std::cout << point1 << "\t" << value1[0] << std::endl;
          std::cout << point2 << "\t" << value2[0] << std::endl;
          AssertThrow(false, ExcInternalError());
        }
      else
        {
          std::cout << point1 << "\t"
                    << "pass" << std::endl;
        }
    }
  for (unsigned int i = 1; i < n_points; ++i)
    {
      Vector<double> value1(1);
      Vector<double> value2(1);

      Point<2> point1;
      point1[1] = -numbers::PI + 2. * i / n_points + eps;
      point1[0] = -numbers::PI;
      Point<2> point2;
      point2[1] = -numbers::PI + 2. * i / n_points + eps;
      point2[0] = numbers::PI;

      VectorTools::point_value(dof_handler, solution, point1, value1);
      VectorTools::point_value(dof_handler, solution, point2, value2);

      if (std::abs(value2[0] - value1[0]) > 1e-8)
        {
          std::cout << point1 << "\t"
                    << "fail" << std::endl;
          std::cout << point1 << "\t" << value1[0] << std::endl;
          std::cout << point2 << "\t" << value2[0] << std::endl;
          Assert(false, ExcInternalError());
        }
      else
        {
          std::cout << point1 << "\t"
                    << "pass" << std::endl;
        }
    }
}


int
main(int argc, char *argv[])
{
  initlog();

  constexpr int      dim = 2;
  const double       L   = numbers::PI;
  Triangulation<dim> triangulation;
  GridGenerator::hyper_cube(triangulation, -L, L, true);

  triangulation.refine_global(1);
  typename Triangulation<dim>::active_cell_iterator cellBegin =
    triangulation.begin_active();
  cellBegin->set_refine_flag();
  triangulation.execute_coarsening_and_refinement();

  FE_Q<dim>       fe(1);
  DoFHandler<dim> dof_handler(triangulation);
  dof_handler.distribute_dofs(fe);

  std::vector<AffineConstraints<double>> constraints(8);

  PeriodicReference<dim> periodic_function;

  std::vector<Vector<double>> projection(8,
                                         Vector<double>(dof_handler.n_dofs()));

  for (unsigned int i = 0; i < 8; ++i)
    {
      deallog << "Testing version " << i << std::endl;
      constraints[i] = make_constraint_matrix(dof_handler, i);
      VectorTools::project(dof_handler,
                           constraints[i],
                           QGauss<dim>(3),
                           periodic_function,
                           projection[i]);
      check_periodicity(dof_handler, projection[i], i);
    }
}
