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



// For all combinations of 2-cell triangulations that can be created by
// GridGenerator::non_standard_orientation_mesh, set up a discrete level set
// function corresponding to
// (y  - 0.25) + (z - 0.25) = 0
// and use NonMatching::FEInterfaceValues to generate the face quadratures on
// the shared face at x = 0. Check that the quadrature points on both sides of
// the face are the same.


#include <deal.II/base/exceptions.h>
#include <deal.II/base/function_signed_distance.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/non_matching/fe_values.h>

#include <deal.II/numerics/vector_tools.h>

#include <bitset>
#include <optional>

#include "../tests.h"


template <int dim>
class Test
{
public:
  Test(const Triangulation<dim> &triangulation);

  void
  run();

private:
  void
  setup_discrete_level_set();

  /*
   * For the incoming FEInterfaceValues object, print the distance between the
   * real space quadrature points of the underlying FEFaceValues objects, for
   * all quadrature points.
   */
  void
  print_quadrature_point_difference_on_both_sides(
    const dealii::FEInterfaceValues<dim> &fe_interface_values);

  hp::FECollection<dim> fe_collection;
  DoFHandler<dim>       dof_handler;
  Vector<double>        level_set;

  NonMatching::MeshClassifier<dim> mesh_classifier;

  QGauss<1> quadrature_1D;
};



template <int dim>
Test<dim>::Test(const Triangulation<dim> &triangulation)
  : dof_handler(triangulation)
  , mesh_classifier(dof_handler, level_set)
  , quadrature_1D(2)
{
  fe_collection.push_back(FE_Q<dim>(1));
}



template <int dim>
void
Test<dim>::setup_discrete_level_set()
{
  dof_handler.distribute_dofs(fe_collection);
  level_set.reinit(dof_handler.n_dofs());

  Point<dim>     point_on_zero_contour;
  Tensor<1, dim> plane_normal;
  for (unsigned int i = 1; i < dim; i++)
    {
      point_on_zero_contour[i] = 0.25;
      plane_normal[i]          = 1;
    }

  const Functions::SignedDistance::Plane<dim> analytical_levelset(
    point_on_zero_contour, plane_normal);

  VectorTools::interpolate(dof_handler, analytical_levelset, level_set);
}



template <int dim>
void
Test<dim>::print_quadrature_point_difference_on_both_sides(
  const dealii::FEInterfaceValues<dim> &fe_interface_values)
{
  for (const unsigned int q : fe_interface_values.quadrature_point_indices())
    {
      const Point<dim> &point_on_side_0 =
        fe_interface_values.get_fe_face_values(0).quadrature_point(q);
      const Point<dim> &point_on_side_1 =
        fe_interface_values.get_fe_face_values(1).quadrature_point(q);

      deallog << point_on_side_0.distance(point_on_side_1) << std::endl;
    }
}



template <int dim>
void
Test<dim>::run()
{
  setup_discrete_level_set();
  mesh_classifier.reclassify();

  NonMatching::RegionUpdateFlags region_update_flags;
  region_update_flags.inside  = update_quadrature_points;
  region_update_flags.outside = update_quadrature_points;

  NonMatching::FEInterfaceValues<dim> fe_values(fe_collection,
                                                quadrature_1D,
                                                region_update_flags,
                                                mesh_classifier,
                                                dof_handler,
                                                level_set);

  const auto cell = dof_handler.begin_active();

  // Test that the points agree only on the interface between the two cells,
  // which face this is depends on the incoming triangulation.
  for (const unsigned int face : cell->face_indices())
    if (!cell->face(face)->at_boundary())
      {
        const unsigned int invalid_subface =
          dealii::numbers::invalid_unsigned_int;

        fe_values.reinit(cell,
                         face,
                         invalid_subface,
                         cell->neighbor(face),
                         cell->neighbor_of_neighbor(face),
                         invalid_subface);

        Assert(fe_values.get_inside_fe_values(), ExcInternalError());
        Assert(fe_values.get_outside_fe_values(), ExcInternalError());

        print_quadrature_point_difference_on_both_sides(
          *fe_values.get_inside_fe_values());
        print_quadrature_point_difference_on_both_sides(
          *fe_values.get_outside_fe_values());
      }
}



int
main()
{
  initlog();
  const int dim = 3;

  // non_standard_orientation_mesh takes 4 boolean parameters. Test all
  // combinations.
  const unsigned int n_arguments    = 4;
  const unsigned int n_combinations = std::pow(2, 4);
  for (unsigned int i = 0; i < n_combinations; i++)
    {
      const std::bitset<n_arguments> bits(i);

      const bool face_orientation     = bits[0];
      const bool face_flip            = bits[1];
      const bool face_rotation        = bits[2];
      const bool manipulate_left_cube = bits[3];

      Triangulation<dim> triangulation;
      GridGenerator::non_standard_orientation_mesh(triangulation,
                                                   face_orientation,
                                                   face_flip,
                                                   face_rotation,
                                                   manipulate_left_cube);
      Test<dim> test(triangulation);
      test.run();
    }
}
