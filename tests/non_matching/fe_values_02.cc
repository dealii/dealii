// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2021 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#include <deal.II/base/function_signed_distance.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/non_matching/fe_values.h>

#include <deal.II/numerics/vector_tools.h>

#include <optional>

#include "../tests.h"

/*
 * Set up a triangulation with 2 unit elements in a row: |-0-|-1-|, and a level
 * set function such that it cuts exactly through the face that is shared
 * between the elements. Test that the sum of the JxW-values we get from
 * inside/outside/surface_fe_values is the correct volume/surface area. That is,
 * all should be equal to 1.
 */
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

  /*
   * Setup a discrete level set function corresponding to
   * $\psi(x) = x_0$
   */
  void
  setup_discrete_level_set();

  void
  compute_volumes_and_surface_area(NonMatching::FEValues<dim> &fe_values) const;

  Triangulation<dim>    triangulation;
  hp::FECollection<dim> fe_collection;
  DoFHandler<dim>       dof_handler;

  hp::MappingCollection<dim> mapping_collection;
  hp::QCollection<dim>       q_collection;
  hp::QCollection<1>         q_collection1D;

  Vector<double>                   level_set;
  NonMatching::MeshClassifier<dim> mesh_classifier;
};



template <int dim>
Test<dim>::Test()
  : dof_handler(triangulation)
  , mesh_classifier(dof_handler, level_set)
{
  fe_collection.push_back(FE_Q<dim>(1));
  mapping_collection.push_back(MappingCartesian<dim>());
  const unsigned int n_quadrature_points = 1;
  q_collection.push_back(QGauss<dim>(n_quadrature_points));
  q_collection1D.push_back(QGauss<1>(n_quadrature_points));
}



template <int dim>
void
Test<dim>::run()
{
  setup_mesh();
  dof_handler.distribute_dofs(fe_collection);
  setup_discrete_level_set();
  mesh_classifier.reclassify();

  NonMatching::RegionUpdateFlags region_update_flags;
  region_update_flags.inside  = update_JxW_values;
  region_update_flags.outside = update_JxW_values;
  region_update_flags.surface = update_JxW_values;

  NonMatching::FEValues<dim> fe_values(fe_collection,
                                       q_collection1D[0],
                                       region_update_flags,
                                       mesh_classifier,
                                       dof_handler,
                                       level_set);
  compute_volumes_and_surface_area(fe_values);
}



template <int dim>
void
Test<dim>::setup_mesh()
{
  Point<dim> lower_left;
  lower_left[0] = -1;
  Point<dim> upper_right;

  std::vector<unsigned int> repetitions;
  upper_right[0] = 1;
  repetitions.push_back(2);
  for (unsigned int d = 1; d < dim; ++d)
    {
      upper_right[d] = 1;
      repetitions.push_back(1);
    }

  GridGenerator::subdivided_hyper_rectangle(triangulation,
                                            repetitions,
                                            lower_left,
                                            upper_right);
}



template <int dim>
void
Test<dim>::setup_discrete_level_set()
{
  const Point<dim> point_on_zero_contour;

  const Functions::SignedDistance::Plane<dim> analytical_levelset(
    point_on_zero_contour, Point<dim>::unit_vector(0));

  level_set.reinit(dof_handler.n_dofs());
  VectorTools::interpolate(dof_handler, analytical_levelset, level_set);
}



template <class FE_VALUES_TYPE>
double
sum_JxW_values(const FE_VALUES_TYPE &fe_values)
{
  double sum = 0;
  for (unsigned int q : fe_values.quadrature_point_indices())
    sum += fe_values.JxW(q);
  return sum;
}



template <int dim>
void
Test<dim>::compute_volumes_and_surface_area(
  NonMatching::FEValues<dim> &fe_values) const
{
  double surface_area   = 0;
  double inside_volume  = 0;
  double outside_volume = 0;

  for (auto &cell : dof_handler.active_cell_iterators())
    {
      fe_values.reinit(cell);

      const auto &inside_fe_values = fe_values.get_inside_fe_values();
      if (inside_fe_values)
        inside_volume += sum_JxW_values(*inside_fe_values);

      const auto &outside_fe_values = fe_values.get_outside_fe_values();
      if (outside_fe_values)
        outside_volume += sum_JxW_values(*outside_fe_values);

      const auto &fe_surface_values = fe_values.get_surface_fe_values();
      if (fe_surface_values)
        surface_area += sum_JxW_values(*fe_surface_values);
    }

  deallog << "inside volume = " << inside_volume << std::endl;
  deallog << "surface area = " << surface_area << std::endl;
  deallog << "outside volume = " << outside_volume << std::endl;
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
