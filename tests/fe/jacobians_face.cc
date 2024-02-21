// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2015 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Show the Jacobians and inverse Jacobians of FEFaceValues and
// FESubfaceValues on a hyperball mesh with one quadrature point

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/fe/fe_nothing.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"

template <int dim>
void
test()
{
  Triangulation<dim> tria;
  GridGenerator::hyper_ball(tria);
  static const SphericalManifold<dim> manifold;
  tria.set_all_manifold_ids_on_boundary(0);
  tria.set_manifold(0, manifold);
  tria.begin_active()->set_refine_flag();
  tria.execute_coarsening_and_refinement();

  MappingQ<dim>   mapping(5);
  FE_Nothing<dim> dummy;
  // choose a point that is not right in the middle of the cell so that the
  // Jacobian contains many nonzero entries
  Point<dim - 1> quad_p;
  for (int d = 0; d < dim - 1; ++d)
    quad_p[d] = 0.42 + 0.11 * d;
  Quadrature<dim - 1> quad(quad_p);

  {
    FEFaceValues<dim>    fe_val(mapping, dummy, quad, update_jacobians);
    FESubfaceValues<dim> fe_sub_val(mapping, dummy, quad, update_jacobians);

    deallog << dim << "D Jacobians:" << std::endl;
    typename Triangulation<dim>::active_cell_iterator cell =
                                                        tria.begin_active(),
                                                      endc = tria.end();
    for (; cell != endc; ++cell)
      for (const unsigned int f : GeometryInfo<dim>::face_indices())
        {
          fe_val.reinit(cell, f);

          for (unsigned int d = 0; d < dim; ++d)
            for (unsigned int e = 0; e < dim; ++e)
              deallog << fe_val.jacobian(0)[d][e] << ' ';
          deallog << std::endl;

          // Also check the Jacobian with FESubfaceValues
          if (cell->at_boundary(f) == false &&
              cell->neighbor(f)->level() < cell->level())
            {
              fe_sub_val.reinit(cell->neighbor(f),
                                cell->neighbor_face_no(f),
                                cell->neighbor_of_coarser_neighbor(f).second);

              for (unsigned int d = 0; d < dim; ++d)
                for (unsigned int e = 0; e < dim; ++e)
                  deallog << fe_sub_val.jacobian(0)[d][e] << ' ';
              deallog << std::endl;
            }
        }
    deallog << std::endl;
  }

  {
    FEFaceValues<dim>    fe_val(mapping, dummy, quad, update_inverse_jacobians);
    FESubfaceValues<dim> fe_sub_val(mapping,
                                    dummy,
                                    quad,
                                    update_inverse_jacobians);

    deallog << dim << "D inverse Jacobians:" << std::endl;
    typename Triangulation<dim>::active_cell_iterator cell =
                                                        tria.begin_active(),
                                                      endc = tria.end();
    for (; cell != endc; ++cell)
      for (const unsigned int f : GeometryInfo<dim>::face_indices())
        {
          fe_val.reinit(cell, f);

          for (unsigned int d = 0; d < dim; ++d)
            for (unsigned int e = 0; e < dim; ++e)
              deallog << fe_val.inverse_jacobian(0)[d][e] << ' ';
          deallog << std::endl;

          // Also check the inverse Jacobian with FESubfaceValues
          if (cell->at_boundary(f) == false &&
              cell->neighbor(f)->level() < cell->level())
            {
              fe_sub_val.reinit(cell->neighbor(f),
                                cell->neighbor_face_no(f),
                                cell->neighbor_of_coarser_neighbor(f).second);

              for (unsigned int d = 0; d < dim; ++d)
                for (unsigned int e = 0; e < dim; ++e)
                  deallog << fe_sub_val.inverse_jacobian(0)[d][e] << ' ';
              deallog << std::endl;
            }
        }
    deallog << std::endl;
  }

  {
    FEFaceValues<dim>    fe_val(mapping, dummy, quad, update_jacobian_grads);
    FESubfaceValues<dim> fe_sub_val(mapping,
                                    dummy,
                                    quad,
                                    update_jacobian_grads);

    deallog << dim << "D Jacobian gradients:" << std::endl;
    typename Triangulation<dim>::active_cell_iterator cell =
                                                        tria.begin_active(),
                                                      endc = tria.end();
    for (; cell != endc; ++cell)
      for (const unsigned int f : GeometryInfo<dim>::face_indices())
        {
          fe_val.reinit(cell, f);

          for (unsigned int d = 0; d < dim; ++d)
            for (unsigned int e = 0; e < dim; ++e)
              for (unsigned int f = 0; f < dim; ++f)
                deallog << fe_val.jacobian_grad(0)[d][e][f] << ' ';
          deallog << std::endl;

          // Also check the Jacobian with FESubfaceValues
          if (cell->at_boundary(f) == false &&
              cell->neighbor(f)->level() < cell->level())
            {
              fe_sub_val.reinit(cell->neighbor(f),
                                cell->neighbor_face_no(f),
                                cell->neighbor_of_coarser_neighbor(f).second);

              for (unsigned int d = 0; d < dim; ++d)
                for (unsigned int e = 0; e < dim; ++e)
                  for (unsigned int f = 0; f < dim; ++f)
                    deallog << fe_sub_val.jacobian_grad(0)[d][e][f] << ' ';
              deallog << std::endl;
            }
        }
    deallog << std::endl;
  }

  {
    FEFaceValues<dim>    fe_val(mapping,
                             dummy,
                             quad,
                             update_jacobian_pushed_forward_grads);
    FESubfaceValues<dim> fe_sub_val(mapping,
                                    dummy,
                                    quad,
                                    update_jacobian_pushed_forward_grads);

    deallog << dim << "D Jacobian pushed forward gradients:" << std::endl;
    typename Triangulation<dim>::active_cell_iterator cell =
                                                        tria.begin_active(),
                                                      endc = tria.end();
    for (; cell != endc; ++cell)
      for (const unsigned int f : GeometryInfo<dim>::face_indices())
        {
          fe_val.reinit(cell, f);

          for (unsigned int d = 0; d < dim; ++d)
            for (unsigned int e = 0; e < dim; ++e)
              for (unsigned int f = 0; f < dim; ++f)
                deallog << fe_val.jacobian_pushed_forward_grad(0)[d][e][f]
                        << ' ';
          deallog << std::endl;

          // Also check the Jacobian with FESubfaceValues
          if (cell->at_boundary(f) == false &&
              cell->neighbor(f)->level() < cell->level())
            {
              fe_sub_val.reinit(cell->neighbor(f),
                                cell->neighbor_face_no(f),
                                cell->neighbor_of_coarser_neighbor(f).second);

              for (unsigned int d = 0; d < dim; ++d)
                for (unsigned int e = 0; e < dim; ++e)
                  for (unsigned int f = 0; f < dim; ++f)
                    deallog
                      << fe_sub_val.jacobian_pushed_forward_grad(0)[d][e][f]
                      << ' ';
              deallog << std::endl;
            }
        }
    deallog << std::endl;
  }

  {
    FEFaceValues<dim>    fe_val(mapping,
                             dummy,
                             quad,
                             update_jacobian_2nd_derivatives);
    FESubfaceValues<dim> fe_sub_val(mapping,
                                    dummy,
                                    quad,
                                    update_jacobian_2nd_derivatives);

    deallog << dim << "D Jacobian hessians:" << std::endl;
    typename Triangulation<dim>::active_cell_iterator cell =
                                                        tria.begin_active(),
                                                      endc = tria.end();
    for (; cell != endc; ++cell)
      for (const unsigned int f : GeometryInfo<dim>::face_indices())
        {
          fe_val.reinit(cell, f);

          for (unsigned int d = 0; d < dim; ++d)
            for (unsigned int e = 0; e < dim; ++e)
              for (unsigned int f = 0; f < dim; ++f)
                for (unsigned int g = 0; g < dim; ++g)
                  deallog << fe_val.jacobian_2nd_derivative(0)[d][e][f][g]
                          << ' ';
          deallog << std::endl;

          // Also check the Jacobian with FESubfaceValues
          if (cell->at_boundary(f) == false &&
              cell->neighbor(f)->level() < cell->level())
            {
              fe_sub_val.reinit(cell->neighbor(f),
                                cell->neighbor_face_no(f),
                                cell->neighbor_of_coarser_neighbor(f).second);

              for (unsigned int d = 0; d < dim; ++d)
                for (unsigned int e = 0; e < dim; ++e)
                  for (unsigned int f = 0; f < dim; ++f)
                    for (unsigned int g = 0; g < dim; ++g)
                      deallog
                        << fe_sub_val.jacobian_2nd_derivative(0)[d][e][f][g]
                        << ' ';
              deallog << std::endl;
            }
        }
    deallog << std::endl;
  }

  {
    FEFaceValues<dim>    fe_val(mapping,
                             dummy,
                             quad,
                             update_jacobian_pushed_forward_2nd_derivatives);
    FESubfaceValues<dim> fe_sub_val(
      mapping, dummy, quad, update_jacobian_pushed_forward_2nd_derivatives);

    deallog << dim << "D Jacobian pushed forward hessians:" << std::endl;
    typename Triangulation<dim>::active_cell_iterator cell =
                                                        tria.begin_active(),
                                                      endc = tria.end();
    for (; cell != endc; ++cell)
      for (const unsigned int f : GeometryInfo<dim>::face_indices())
        {
          fe_val.reinit(cell, f);

          for (unsigned int d = 0; d < dim; ++d)
            for (unsigned int e = 0; e < dim; ++e)
              for (unsigned int f = 0; f < dim; ++f)
                for (unsigned int g = 0; g < dim; ++g)
                  deallog << fe_val.jacobian_pushed_forward_2nd_derivative(
                               0)[d][e][f][g]
                          << ' ';
          deallog << std::endl;

          // Also check the Jacobian with FESubfaceValues
          if (cell->at_boundary(f) == false &&
              cell->neighbor(f)->level() < cell->level())
            {
              fe_sub_val.reinit(cell->neighbor(f),
                                cell->neighbor_face_no(f),
                                cell->neighbor_of_coarser_neighbor(f).second);

              for (unsigned int d = 0; d < dim; ++d)
                for (unsigned int e = 0; e < dim; ++e)
                  for (unsigned int f = 0; f < dim; ++f)
                    for (unsigned int g = 0; g < dim; ++g)
                      deallog
                        << fe_sub_val.jacobian_pushed_forward_2nd_derivative(
                             0)[d][e][f][g]
                        << ' ';
              deallog << std::endl;
            }
        }
    deallog << std::endl;
  }

  {
    FEFaceValues<dim>    fe_val(mapping,
                             dummy,
                             quad,
                             update_jacobian_3rd_derivatives);
    FESubfaceValues<dim> fe_sub_val(mapping,
                                    dummy,
                                    quad,
                                    update_jacobian_3rd_derivatives);

    deallog << dim << "D Jacobian hessian gradients:" << std::endl;
    typename Triangulation<dim>::active_cell_iterator cell =
                                                        tria.begin_active(),
                                                      endc = tria.end();
    for (; cell != endc; ++cell)
      for (const unsigned int f : GeometryInfo<dim>::face_indices())
        {
          fe_val.reinit(cell, f);

          for (unsigned int d = 0; d < dim; ++d)
            for (unsigned int e = 0; e < dim; ++e)
              for (unsigned int f = 0; f < dim; ++f)
                for (unsigned int g = 0; g < dim; ++g)
                  for (unsigned int h = 0; h < dim; ++h)
                    deallog << fe_val.jacobian_3rd_derivative(0)[d][e][f][g][h]
                            << ' ';
          deallog << std::endl;

          // Also check the Jacobian with FESubfaceValues
          if (cell->at_boundary(f) == false &&
              cell->neighbor(f)->level() < cell->level())
            {
              fe_sub_val.reinit(cell->neighbor(f),
                                cell->neighbor_face_no(f),
                                cell->neighbor_of_coarser_neighbor(f).second);

              for (unsigned int d = 0; d < dim; ++d)
                for (unsigned int e = 0; e < dim; ++e)
                  for (unsigned int f = 0; f < dim; ++f)
                    for (unsigned int g = 0; g < dim; ++g)
                      for (unsigned int h = 0; h < dim; ++h)
                        deallog << fe_sub_val.jacobian_3rd_derivative(
                                     0)[d][e][f][g][h]
                                << ' ';
              deallog << std::endl;
            }
        }
    deallog << std::endl;
  }

  {
    FEFaceValues<dim>    fe_val(mapping,
                             dummy,
                             quad,
                             update_jacobian_pushed_forward_3rd_derivatives);
    FESubfaceValues<dim> fe_sub_val(
      mapping, dummy, quad, update_jacobian_pushed_forward_3rd_derivatives);

    deallog << dim
            << "D Jacobian pushed forward hessian gradients:" << std::endl;
    typename Triangulation<dim>::active_cell_iterator cell =
                                                        tria.begin_active(),
                                                      endc = tria.end();
    for (; cell != endc; ++cell)
      for (const unsigned int f : GeometryInfo<dim>::face_indices())
        {
          fe_val.reinit(cell, f);

          for (unsigned int d = 0; d < dim; ++d)
            for (unsigned int e = 0; e < dim; ++e)
              for (unsigned int f = 0; f < dim; ++f)
                for (unsigned int g = 0; g < dim; ++g)
                  for (unsigned int h = 0; h < dim; ++h)
                    deallog << fe_val.jacobian_pushed_forward_3rd_derivative(
                                 0)[d][e][f][g][h]
                            << ' ';
          deallog << std::endl;

          // Also check the Jacobian with FESubfaceValues
          if (cell->at_boundary(f) == false &&
              cell->neighbor(f)->level() < cell->level())
            {
              fe_sub_val.reinit(cell->neighbor(f),
                                cell->neighbor_face_no(f),
                                cell->neighbor_of_coarser_neighbor(f).second);

              for (unsigned int d = 0; d < dim; ++d)
                for (unsigned int e = 0; e < dim; ++e)
                  for (unsigned int f = 0; f < dim; ++f)
                    for (unsigned int g = 0; g < dim; ++g)
                      for (unsigned int h = 0; h < dim; ++h)
                        deallog
                          << fe_sub_val.jacobian_pushed_forward_3rd_derivative(
                               0)[d][e][f][g][h]
                          << ' ';
              deallog << std::endl;
            }
        }
    deallog << std::endl;
  }
}


int
main()
{
  initlog();
  deallog << std::setprecision(8) << std::fixed;

  test<2>();
  test<3>();

  return 0;
}
