// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Evaluates the area of a complex geometry when mixing different manifolds on
// the surfaces than the rest of the cell. The mesh uses two layers of
// cylinders that finally are changed into a square on the outside of the
// domain.

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_accessor.h>

#include <deal.II/fe/fe_nothing.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/vector.h>

#include <sstream>
#include <string>
#include <vector>

#include "../tests.h"

const double D   = 0.1;
const double R   = D / 2.0;
const double R_1 = 1.2 * R;
const double R_2 = 1.7 * R;
const double H   = 0.41;
const double X_0 = 0.0;
const double X_1 = 0.3;
const double X_C = 0.5; // center
const double X_2 = 0.7;

const double Y_0 = 0.0;
const double Y_C = 0.2; // center

const unsigned int MANIFOLD_ID = 1;


void
create_triangulation(Triangulation<2> &tria)
{
  AssertThrow(std::abs((X_2 - X_1) - 2.0 * (X_C - X_1)) < 1.0e-12,
              ExcMessage("Geometry parameters X_1,X_2,X_C invalid!"));
  SphericalManifold<2> spherical_manifold(Point<2>(X_C, Y_C));

  Triangulation<2> circle_1, circle_2, circle_tmp, middle, middle_tmp,
    middle_tmp2, tmp_3D;
  std::vector<unsigned int> ref_1(2, 2);
  ref_1[1] = 2;

  // create middle part first as a hyper shell
  const double       outer_radius = (X_2 - X_1) / 2.0;
  const unsigned int n_cells      = 4;
  GridGenerator::hyper_shell(
    middle, Point<2>(X_C, Y_C), R_2, outer_radius, n_cells, true);
  middle.set_all_manifold_ids(MANIFOLD_ID);
  middle.set_manifold(MANIFOLD_ID, spherical_manifold);
  middle.refine_global(1);

  // two inner circles in order to refine towards the cylinder surface
  const unsigned int n_cells_circle = 8;
  GridGenerator::hyper_shell(
    circle_1, Point<2>(X_C, Y_C), R, R_1, n_cells_circle, true);
  circle_1.set_all_manifold_ids(MANIFOLD_ID);
  circle_1.set_manifold(MANIFOLD_ID, spherical_manifold);

  GridGenerator::hyper_shell(
    circle_2, Point<2>(X_C, Y_C), R_1, R_2, n_cells_circle, true);
  circle_2.set_all_manifold_ids(MANIFOLD_ID);
  circle_2.set_manifold(MANIFOLD_ID, spherical_manifold);

  // then move the vertices to the points where we want them to be to create a
  // slightly asymmetric cube with a hole
  for (Triangulation<2>::cell_iterator cell = middle.begin();
       cell != middle.end();
       ++cell)
    {
      for (const unsigned int v : GeometryInfo<2>::vertex_indices())
        {
          Point<2> &vertex = cell->vertex(v);
          if (std::abs(vertex[0] - X_2) < 1e-10 &&
              std::abs(vertex[1] - Y_C) < 1e-10)
            {
              vertex = Point<2>(X_2, H / 2.0);
            }
          else if (std::abs(vertex[0] -
                            (X_C + (X_2 - X_1) / 2.0 / std::sqrt(2))) < 1e-10 &&
                   std::abs(vertex[1] -
                            (Y_C + (X_2 - X_1) / 2.0 / std::sqrt(2))) < 1e-10)
            {
              vertex = Point<2>(X_2, H);
            }
          else if (std::abs(vertex[0] -
                            (X_C + (X_2 - X_1) / 2.0 / std::sqrt(2))) < 1e-10 &&
                   std::abs(vertex[1] -
                            (Y_C - (X_2 - X_1) / 2.0 / std::sqrt(2))) < 1e-10)
            {
              vertex = Point<2>(X_2, Y_0);
            }
          else if (std::abs(vertex[0] - X_C) < 1e-10 &&
                   std::abs(vertex[1] - (Y_C + (X_2 - X_1) / 2.0)) < 1e-10)
            {
              vertex = Point<2>(X_C, H);
            }
          else if (std::abs(vertex[0] - X_C) < 1e-10 &&
                   std::abs(vertex[1] - (Y_C - (X_2 - X_1) / 2.0)) < 1e-10)
            {
              vertex = Point<2>(X_C, Y_0);
            }
          else if (std::abs(vertex[0] -
                            (X_C - (X_2 - X_1) / 2.0 / std::sqrt(2))) < 1e-10 &&
                   std::abs(vertex[1] -
                            (Y_C + (X_2 - X_1) / 2.0 / std::sqrt(2))) < 1e-10)
            {
              vertex = Point<2>(X_1, H);
            }
          else if (std::abs(vertex[0] -
                            (X_C - (X_2 - X_1) / 2.0 / std::sqrt(2))) < 1e-10 &&
                   std::abs(vertex[1] -
                            (Y_C - (X_2 - X_1) / 2.0 / std::sqrt(2))) < 1e-10)
            {
              vertex = Point<2>(X_1, Y_0);
            }
          else if (std::abs(vertex[0] - X_1) < 1e-10 &&
                   std::abs(vertex[1] - Y_C) < 1e-10)
            {
              vertex = Point<2>(X_1, H / 2.0);
            }
        }
    }

  // must copy the triangulation because we cannot merge triangulations with
  // refinement...
  GridGenerator::flatten_triangulation(middle, middle_tmp);
  GridGenerator::merge_triangulations(circle_1, circle_2, circle_tmp);
  GridGenerator::merge_triangulations(middle_tmp, circle_tmp, tria);

  // Set the cylinder boundary  to 2, outflow to 1, the rest to 0.
  // tria.set_all_manifold_ids(0);
  for (Triangulation<2>::active_cell_iterator cell = tria.begin();
       cell != tria.end();
       ++cell)
    {
      if (Point<2>(X_C, Y_C).distance(cell->center()) <= R_2)
        cell->set_all_manifold_ids(MANIFOLD_ID);
    }
  tria.set_manifold(MANIFOLD_ID, FlatManifold<2>());
}

void
create_triangulation(Triangulation<3> &tria)
{
  Triangulation<2> tria_2d;
  create_triangulation(tria_2d);
  GridGenerator::extrude_triangulation(tria_2d, 3, H, tria);

  // Set the cylinder boundary  to 2, outflow to 1, the rest to 0.
  tria.set_all_manifold_ids(0);
  for (Triangulation<3>::active_cell_iterator cell = tria.begin();
       cell != tria.end();
       ++cell)
    {
      if (Point<3>(X_C, Y_C, cell->center()[2]).distance(cell->center()) <= R_2)
        cell->set_all_manifold_ids(MANIFOLD_ID);
    }
  tria.set_manifold(0, FlatManifold<3>());
  tria.set_manifold(MANIFOLD_ID, FlatManifold<3>());
}

template <int dim>
struct ManifoldWrapper
{
  Manifold<dim> *
  operator()(const Tensor<1, dim> &direction, const Point<dim> &center) const;
};

template <>
Manifold<2> *
ManifoldWrapper<2>::operator()(const Tensor<1, 2> & /*direction*/,
                               const Point<2> &center) const
{
  return new SphericalManifold<2>(center);
}

template <>
Manifold<3> *
ManifoldWrapper<3>::operator()(const Tensor<1, 3> &direction,
                               const Point<3>     &center) const
{
  return new CylindricalManifold<3>(direction, center);
}

template <int dim>
void
test()
{
  Point<dim> center;
  center[0] = X_C;
  center[1] = Y_C;
  Tensor<1, dim> direction;
  direction[dim - 1] = 1.;

  std::shared_ptr<Manifold<dim>> cylinder_manifold(
    ManifoldWrapper<dim>()(direction, center));
  Triangulation<dim> tria;
  create_triangulation(tria);
  tria.set_manifold(MANIFOLD_ID, *cylinder_manifold);

  FE_Nothing<dim> fe;
  for (unsigned int degree = 1; degree < 7; ++degree)
    {
      MappingQ<dim> mapping(degree);
      QGauss<dim>   quad(degree + 1);
      FEValues<dim> fe_values(mapping, fe, quad, update_JxW_values);
      double        sum = 0.;
      for (typename Triangulation<dim>::active_cell_iterator cell =
             tria.begin_active();
           cell != tria.end();
           ++cell)
        {
          fe_values.reinit(cell);
          double local_sum = 0;
          for (unsigned int q = 0; q < quad.size(); ++q)
            local_sum += fe_values.JxW(q);
          sum += local_sum;
        }
      const double reference =
        (dim == 2 ? 1. : H) * (H * (X_2 - X_1) - D * D * numbers::PI * 0.25);
      deallog << "Volume " << dim << "D mapping degree " << degree << ": "
              << sum << " error: " << (sum - reference) / reference
              << std::endl;
    }
}


int
main()
{
  initlog();

  test<2>();
  test<3>();
}
