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

#include <deal.II/base/exceptions.h>
#include <deal.II/base/function_signed_distance.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_cartesian.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/non_matching/fe_values.h>

#include <deal.II/numerics/vector_tools.h>

#include <optional>

#include "../tests.h"


// Assert that the two incoming cells are the same. Throw an exception if not.
template <int dim>
void
assert_cells_are_the_same(
  const typename Triangulation<dim>::cell_iterator &expected,
  const typename Triangulation<dim>::cell_iterator &cell)
{
  AssertThrow(expected == cell, ExcInternalError());
  deallog << "OK" << std::endl;
}



// Assert that the incoming optional does not have a value.
template <class FEVALUES>
void
assert_doesnt_have_value(const std::optional<FEVALUES> &fe_values)
{
  AssertThrow(!fe_values, ExcInternalError());
  deallog << "OK" << std::endl;
}



// Set up a triangulation with 3 elements in a row |-0-|-1-|-2-|.
// and a level set function such that
//
// element  LocationToLevelSet
//   0           inside
//   1         intersected
//   2           outside
//
// Test the following:
// 1. That we can construct a NonMatching::FEValues object.
// 2. That when we call reinit on each cell in the triangulation,
// the std::optionals that we get from
// get_inside/outside/surface_fe_values are set up correctly. That is, the
// optionals that should contain a value do and the ones that should not contain
// a value do not.
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

  template <typename IteratorType>
  void
  test_fe_values_reinitializes_correctly(NonMatching::FEValues<dim> &fe_values,
                                         IteratorType cell) const;

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

  const NonMatching::RegionUpdateFlags region_update_flags;

  {
    // Test with the "simple" constructor.
    NonMatching::FEValues<dim> fe_values(fe_collection,
                                         q_collection1D[0],
                                         region_update_flags,
                                         mesh_classifier,
                                         dof_handler,
                                         level_set);
    test_fe_values_reinitializes_correctly(fe_values,
                                           triangulation.begin_active());
    test_fe_values_reinitializes_correctly(fe_values,
                                           dof_handler.begin_active());
  }
  {
    // Test with the "more advanced" constructor.
    NonMatching::FEValues<dim> fe_values(mapping_collection,
                                         fe_collection,
                                         q_collection,
                                         q_collection1D,
                                         region_update_flags,
                                         mesh_classifier,
                                         dof_handler,
                                         level_set);
    test_fe_values_reinitializes_correctly(fe_values,
                                           triangulation.begin_active());
    test_fe_values_reinitializes_correctly(fe_values,
                                           dof_handler.begin_active());
  }
}



template <int dim>
void
Test<dim>::setup_mesh()
{
  const Point<dim> lower_left;
  Point<dim>       upper_right;

  std::vector<unsigned int> repetitions;
  upper_right[0] = 3;
  repetitions.push_back(3);
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
  Point<dim> point_on_zero_contour;
  point_on_zero_contour[0] = 1.5;
  const Functions::SignedDistance::Plane<dim> analytical_levelset(
    point_on_zero_contour, Point<dim>::unit_vector(0));

  level_set.reinit(dof_handler.n_dofs());
  VectorTools::interpolate(dof_handler, analytical_levelset, level_set);
}



template <int dim>
template <typename IteratorType>
void
Test<dim>::test_fe_values_reinitializes_correctly(
  NonMatching::FEValues<dim> &fe_values,
  IteratorType                cell) const
{
  //  The first is inside so only the inside FEValues object should be
  //  initialized.
  fe_values.reinit(cell);

  assert_cells_are_the_same<dim>(cell,
                                 fe_values.get_inside_fe_values()->get_cell());
  assert_doesnt_have_value(fe_values.get_surface_fe_values());
  assert_doesnt_have_value(fe_values.get_outside_fe_values());


  //  The second is intersected so all FEValues object should be
  //  initialized.
  cell++;
  fe_values.reinit(cell);

  assert_cells_are_the_same<dim>(cell,
                                 fe_values.get_inside_fe_values()->get_cell());
  assert_cells_are_the_same<dim>(cell,
                                 fe_values.get_outside_fe_values()->get_cell());
  assert_cells_are_the_same<dim>(cell,
                                 fe_values.get_surface_fe_values()->get_cell());

  //  The third is outside so only the outside FEValues object should be
  //  initialized.
  cell++;
  fe_values.reinit(cell);

  assert_doesnt_have_value(fe_values.get_inside_fe_values());
  assert_doesnt_have_value(fe_values.get_surface_fe_values());
  assert_cells_are_the_same<dim>(cell,
                                 fe_values.get_outside_fe_values()->get_cell());
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
