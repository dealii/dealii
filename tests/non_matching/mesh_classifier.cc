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

// Test the MeshClassifier class by classifying a triangulation with a single
// cell with a few different level set functions.

#include <deal.II/base/exceptions.h>
#include <deal.II/base/function.h>
#include <deal.II/base/function_signed_distance.h>
#include <deal.II/base/point.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/vector.h>

#include <deal.II/non_matching/mesh_classifier.h>

#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"


// Return a string with the same name as the incoming value of
// LocationToLevelSet.
std::string
location_to_string(const NonMatching::LocationToLevelSet location)
{
  std::string name;
  switch (location)
    {
      case NonMatching::LocationToLevelSet::inside:
        name = "inside";
        break;
      case NonMatching::LocationToLevelSet::outside:
        name = "outside";
        break;
      case NonMatching::LocationToLevelSet::intersected:
        name = "intersected";
        break;
      case NonMatching::LocationToLevelSet::unassigned:
        name = "unassigned";
        break;
      default:
        AssertThrow(false, ExcInternalError());
    }
  return name;
}



// Print LocationToLevelSet (as determined by the incoming MeshClassifier) to
// deallog, for the incoming cell and all of its faces.
template <int dim>
void
print_cell_and_face_locations(
  const NonMatching::MeshClassifier<dim>                  &classifier,
  const typename Triangulation<dim>::active_cell_iterator &cell)
{
  const NonMatching::LocationToLevelSet cell_location =
    classifier.location_to_level_set(cell);
  deallog << "cell " << location_to_string(cell_location) << std::endl;

  for (const unsigned int f : cell->face_indices())
    {
      const NonMatching::LocationToLevelSet face_location =
        classifier.location_to_level_set(cell, f);
      deallog << "face " << f << ' ' << location_to_string(face_location)
              << std::endl;
    }
}



// Test the version MeshClassifier that takes a Vector and a DoFHandler.
//
// Set up single cell triangulation over [-1, 1]^dim and a DoFHandler using
// FE_Q<1>(1). Interpolate the incoming level set function to the discrete
// space. Classify the cell and its faces and print the result to deallog.
template <int dim>
void
classify_with_discrete_level_set(const Function<dim> &level_set)
{
  deallog << "discrete:" << std::endl;
  Triangulation<dim> triangulation;
  GridGenerator::hyper_cube(triangulation, -1, 1);

  const FE_Q<dim> element(1);
  DoFHandler<dim> dof_handler(triangulation);
  dof_handler.distribute_dofs(element);

  Vector<double> discrete_level_set;
  discrete_level_set.reinit(dof_handler.n_dofs());
  VectorTools::interpolate(dof_handler, level_set, discrete_level_set);

  NonMatching::MeshClassifier<dim> classifier(dof_handler, discrete_level_set);
  classifier.reclassify();

  const typename Triangulation<dim>::active_cell_iterator cell =
    triangulation.begin_active();
  print_cell_and_face_locations(classifier, cell);
  deallog << std::endl;
}



// Test the version of MeshClassifier that takes a Function.
//
// Set up single cell triangulation over [-1, 1]^dim. Classify the cell and its
// faces and print the result to deallog.
template <int dim>
void
classify_with_analytic_level_set(const Function<dim> &level_set)
{
  deallog << "analytic:" << std::endl;
  Triangulation<dim> triangulation;
  GridGenerator::hyper_cube(triangulation, -1, 1);

  const FE_Q<dim> element(1);

  NonMatching::MeshClassifier<dim> classifier(triangulation,
                                              level_set,
                                              element);
  classifier.reclassify();

  const typename Triangulation<dim>::active_cell_iterator cell =
    triangulation.begin_active();
  print_cell_and_face_locations(classifier, cell);
  deallog << std::endl;
}



// Test both the version of MeshClassifier that takes a Vector and a DoFHandler,
// and the version that takes a Function.
template <int dim>
void
classify_with_discrete_and_analytic_level_set(const Function<dim> &level_set)
{
  deallog << std::endl;

  classify_with_discrete_level_set(level_set);

  classify_with_analytic_level_set(level_set);
}



// Test MeshClassifier with a level set function that is constant and positive.
template <int dim>
void
test_positive_function()
{
  deallog << "test_positive_function" << std::endl;

  const Functions::ConstantFunction<dim> level_set(1);

  classify_with_discrete_and_analytic_level_set(level_set);
}



// Test MeshClassifier with a level set function that is constant and negative.
template <int dim>
void
test_negative_function()
{
  deallog << "test_negative_function" << std::endl;

  const Functions::ConstantFunction<dim> level_set(-1);

  classify_with_discrete_and_analytic_level_set(level_set);
}



// Test MeshClassifier with a level set function corresponding to the plane
// (x = 0) intersecting the hypercube [-1, 1]^dim.
template <int dim>
void
test_intersection_x_eq_0_plane()
{
  deallog << "test_intersection_x_eq_0_plane" << std::endl;

  Tensor<1, dim> plane_normal;
  plane_normal[0] = 1;
  const Point<dim> origo;

  const Functions::SignedDistance::Plane<dim> level_set(origo, plane_normal);

  classify_with_discrete_and_analytic_level_set(level_set);
}



// Set up local level set coefficients for an Q2 element such that all
// coefficients are positive but the cell is still intersected.
template <int dim>
void
setup_intersected_Q2_positive_coefficients(Vector<double> &level_set){
  // This test case only makes sense in 2D and 3D, specialize for these below
  // and do nothing by default.
};



template <>
void
setup_intersected_Q2_positive_coefficients<2>(Vector<double> &level_set)
{
  const double delta = 1e-5;

  // Listing entries lexiographically.
  level_set(0) = 100 * delta;
  level_set(6) = delta;
  level_set(1) = delta;
  level_set(4) = .5;
  level_set(8) = .5;
  level_set(5) = .5;
  level_set(2) = 1;
  level_set(7) = 1;
  level_set(3) = 1;
}



template <>
void
setup_intersected_Q2_positive_coefficients<3>(Vector<double> &level_set)
{
  const double delta = 1e-5;

  level_set(0)  = 100 * delta;
  level_set(10) = 100 * delta;
  level_set(1)  = 100 * delta;

  level_set(8)  = 100 * delta;
  level_set(24) = delta;
  level_set(9)  = delta;

  level_set(2)  = 100 * delta;
  level_set(11) = 100 * delta;
  level_set(3)  = 100 * delta;

  level_set(16) = .5;
  level_set(22) = .5;
  level_set(17) = .5;

  level_set(20) = .5;
  level_set(26) = .5;
  level_set(21) = .5;

  level_set(18) = .5;
  level_set(23) = .5;
  level_set(19) = .5;

  level_set(4)  = 1;
  level_set(14) = 1;
  level_set(5)  = 1;

  level_set(12) = 1;
  level_set(25) = 1;
  level_set(13) = 1;

  level_set(6)  = 1;
  level_set(15) = 1;
  level_set(7)  = 1;
}



// Test the case that all the Lagrange coefficients of a Q2 element are positive
// but the cell is still intersected. This can happen because the Lagrange shape
// functions are negative between the support points.
template <int dim>
void
test_lagrange_coefficients_positive()
{
  deallog << "test_lagrange_coefficients_positive" << std::endl;
  Triangulation<dim> triangulation;
  GridGenerator::hyper_cube(triangulation);

  const FE_Q<dim> element(2);

  DoFHandler<dim> dof_handler(triangulation);
  dof_handler.distribute_dofs(element);

  Vector<double> level_set(element.dofs_per_cell);
  setup_intersected_Q2_positive_coefficients<dim>(level_set);

  NonMatching::MeshClassifier<dim> classifier(dof_handler, level_set);
  classifier.reclassify();

  const typename Triangulation<dim>::active_cell_iterator cell =
    triangulation.begin_active();
  print_cell_and_face_locations(classifier, cell);
  deallog << std::endl;
}



// Check that the values of LocationToLevelSet for the cells and faces get
// updated correctly when calling reclassify() multiple times.
//
// First, make the level set function all negative, call reclassify(), and check
// that the values of LocationToLevelSet for all cells and faces equals
// LocationToLevelSet::inside. Then, change the level set function to all
// positive, call reclassify() again, and check that all values have been
// changed to LocationToLevelSet::outside.
template <int dim>
void
test_reclassify_called_multiple_times()
{
  deallog << "test_reclassify_called_multiple_times" << std::endl;
  Triangulation<dim> triangulation;
  GridGenerator::hyper_cube(triangulation);

  const FE_Q<dim> element(1);

  DoFHandler<dim> dof_handler(triangulation);
  dof_handler.distribute_dofs(element);

  Vector<double>                   level_set(element.dofs_per_cell);
  NonMatching::MeshClassifier<dim> classifier(dof_handler, level_set);

  const typename Triangulation<dim>::active_cell_iterator cell =
    triangulation.begin_active();

  deallog << "Level set negative" << std::endl;
  level_set = -1;
  classifier.reclassify();
  print_cell_and_face_locations(classifier, cell);

  deallog << "Level set positive" << std::endl;
  level_set = 1;
  classifier.reclassify();
  print_cell_and_face_locations(classifier, cell);

  deallog << std::endl;
}



template <int dim>
void
run_test()
{
  deallog << "dim = " << dim << std::endl;

  test_negative_function<dim>();
  test_positive_function<dim>();
  test_intersection_x_eq_0_plane<dim>();
  // This test doesn't make sense in 1D.
  if (dim != 1)
    test_lagrange_coefficients_positive<dim>();

  test_reclassify_called_multiple_times<dim>();
}



int
main()
{
  initlog();
  run_test<1>();
  run_test<2>();
  run_test<3>();
}
