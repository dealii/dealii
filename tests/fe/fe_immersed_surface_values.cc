// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2020 by the deal.II authors
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


#include <deal.II/base/point.h>
#include <deal.II/base/tensor.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_cartesian.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/non_matching/immersed_surface_quadrature.h>

#include "../tests.h"


using namespace dealii;

// Set up a triangulation with a single Cartesian cell and a finite element
// space with FE_Q-elements. Test that we can construct an
// FEImmersedSurfaceValues object and compute known values when using a
// Cartesian Mapping.
template <int dim>
class Test
{
public:
  Test();

  void
  run();

private:
  //  Construct a triangulation with a single cartesian cell, with extents.
  // [0,2]             in 1D,
  // [0,2]x[0,3]       in 2D, and
  // [0,2]x[0,3]x[0,4] in 3D.
  // Use extends that are different in all directions so that it's easier to
  // check that the various quantities are mapped correctly.
  void
  setup_single_cartesian_cell_triangulation();

  // Make the ImmersedSurfaceQuadraturewith contain a single point in the center
  // of the cell, with weight .5 and a normal different in all components.
  void
  setup_single_point_quadrature();

  // Return a FEImmersedSurfaceValues object created with the the incoming
  // UpdateFlags.
  std::unique_ptr<const FEImmersedSurfaceValues<dim>>
  get_fe_values_object(const UpdateFlags update_flags);

  // Test that the constant n_quadrature_points corresponds to what we sent in.
  void
  test_n_quadrature_points();

  // Test that we can use the get_quadrature function to get the quadrature
  // we passed to the FEImmersedSurfaceValues constructor.
  void
  test_get_quadrature();

  // Print the quadrature point in real space to make sure that it's mapped to
  // real space correctly.
  void
  test_point_mapped_correctly();

  // Print the normal in real space to make sure that it's mapped to real space
  // correctly.
  void
  test_normal();

  // Print JxW to check that the value is correct.
  void
  test_JxW();

  // Test that we can call the shape_surface_grad function.
  void
  test_shape_surface_grad();

  const FE_Q<dim>                             element;
  Triangulation<dim>                          triangulation;
  DoFHandler<dim>                             dof_handler;
  const MappingCartesian<dim>                 mapping;
  NonMatching::ImmersedSurfaceQuadrature<dim> quadrature;
};



template <int dim>
Test<dim>::Test()
  : element(1)
  , dof_handler(triangulation)
{}



template <int dim>
void
Test<dim>::run()
{
  setup_single_cartesian_cell_triangulation();
  dof_handler.distribute_dofs(element);
  setup_single_point_quadrature();

  test_n_quadrature_points();
  test_get_quadrature();
  test_point_mapped_correctly();
  test_normal();
  test_JxW();
  test_shape_surface_grad();
}



template <int dim>
void
Test<dim>::setup_single_cartesian_cell_triangulation()
{
  const Point<dim> lower_left;
  Point<dim>       upper_right;
  for (unsigned int d = 0; d < dim; ++d)
    {
      upper_right[d] = 2 + d;
    }

  GridGenerator::hyper_rectangle(triangulation, lower_left, upper_right);
}



template <int dim>
void
Test<dim>::setup_single_point_quadrature()
{
  Point<dim>     point;
  const double   weight = .5;
  Tensor<1, dim> normal;
  for (unsigned int d = 0; d < dim; ++d)
    {
      point[d]  = .5;
      normal[d] = d + 1;
    }
  normal /= normal.norm();
  quadrature.push_back(point, weight, normal);
}



template <int dim>
std::unique_ptr<const FEImmersedSurfaceValues<dim>>
Test<dim>::get_fe_values_object(const UpdateFlags update_flags)
{
  std::unique_ptr<FEImmersedSurfaceValues<dim>> fe_values(
    new FEImmersedSurfaceValues<dim>(
      mapping, element, quadrature, update_flags));

  const auto cell = dof_handler.begin_active();
  fe_values->reinit(cell);
  return std::move(fe_values);
}



template <int dim>
void
Test<dim>::test_n_quadrature_points()
{
  std::unique_ptr<const FEImmersedSurfaceValues<dim>> fe_values =
    get_fe_values_object(update_default);

  deallog << "n_quadrature_points = " << fe_values->n_quadrature_points
          << std::endl;
}



template <int dim>
void
Test<dim>::test_get_quadrature()
{
  std::unique_ptr<const FEImmersedSurfaceValues<dim>> fe_values =
    get_fe_values_object(update_default);

  const NonMatching::ImmersedSurfaceQuadrature<dim> &stored_quadrature =
    fe_values->get_quadrature();

  for (unsigned int q = 0; q < stored_quadrature.size(); q++)
    {
      deallog << "(point, weight, normal) = ([" << stored_quadrature.point(q)
              << "], " << stored_quadrature.weight(q) << ", ["
              << stored_quadrature.normal_vector(q) << "])" << std::endl;
    }
}



template <int dim>
void
Test<dim>::test_point_mapped_correctly()
{
  std::unique_ptr<const FEImmersedSurfaceValues<dim>> fe_values =
    get_fe_values_object(update_quadrature_points);

  deallog << "point = " << fe_values->quadrature_point(0) << std::endl;
}



template <int dim>
void
Test<dim>::test_normal()
{
  std::unique_ptr<const FEImmersedSurfaceValues<dim>> fe_values =
    get_fe_values_object(update_normal_vectors);

  deallog << "normal = " << fe_values->normal_vector(0) << std::endl;
}



template <int dim>
void
Test<dim>::test_JxW()
{
  std::unique_ptr<const FEImmersedSurfaceValues<dim>> fe_values =
    get_fe_values_object(update_JxW_values);

  deallog << "JxW = " << fe_values->JxW(0) << std::endl;
}



template <int dim>
void
Test<dim>::test_shape_surface_grad()
{
  std::unique_ptr<const FEImmersedSurfaceValues<dim>> fe_values =
    get_fe_values_object(update_gradients | update_normal_vectors);

  const unsigned int function_index = 0;
  const unsigned int q_index        = 0;
  const unsigned int component      = 0;

  deallog << "shape_surface_grad = "
          << fe_values->shape_surface_grad(function_index, q_index)
          << std::endl;

  deallog << "shape_surface_grad_component = "
          << fe_values->shape_surface_grad_component(function_index,
                                                     q_index,
                                                     component)
          << std::endl;
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
