// ---------------------------------------------------------------------
//
// Copyright (C) 2012 - 2020 by the deal.II authors
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


// Show the Jacobians, inverse Jacobians, Jacobian gradients,
// Jacobian Hessians, Jacobian Hessian gradients, and their
// pushed forward versions on hyperball with one quadrature
// point for MappingQ and MappingFEField

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/fe/fe_nothing.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_fe_field.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>

#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"


template <int dim>
void
do_test(const Triangulation<dim> &tria, const Mapping<dim> &mapping)
{
  FE_Nothing<dim> dummy;
  // choose a point that is not right in the
  // middle of the cell so that the Jacobian
  // contains many nonzero entries
  Point<dim> quad_p;
  for (int d = 0; d < dim; ++d)
    quad_p(d) = 0.42 + 0.11 * d;
  Quadrature<dim> quad(quad_p);

  {
    deallog << dim << "D Jacobians:" << std::endl;
    FEValues<dim> fe_val(mapping, dummy, quad, update_jacobians);
    typename Triangulation<dim>::active_cell_iterator cell =
                                                        tria.begin_active(),
                                                      endc = tria.end();
    for (; cell != endc; ++cell)
      {
        fe_val.reinit(cell);

        for (unsigned int d = 0; d < dim; ++d)
          for (unsigned int e = 0; e < dim; ++e)
            deallog << fe_val.jacobian(0)[d][e] << " ";
        deallog << std::endl;
      }
    deallog << std::endl;
  }
  {
    deallog << dim << "D inverse Jacobians:" << std::endl;
    FEValues<dim> fe_val(mapping, dummy, quad, update_inverse_jacobians);
    typename Triangulation<dim>::active_cell_iterator cell =
                                                        tria.begin_active(),
                                                      endc = tria.end();
    for (; cell != endc; ++cell)
      {
        fe_val.reinit(cell);

        for (unsigned int d = 0; d < dim; ++d)
          for (unsigned int e = 0; e < dim; ++e)
            deallog << fe_val.inverse_jacobian(0)[d][e] << " ";
        deallog << std::endl;
      }
    deallog << std::endl;
  }
  {
    deallog << dim << "D Jacobian gradients:" << std::endl;
    FEValues<dim> fe_val(mapping, dummy, quad, update_jacobian_grads);
    typename Triangulation<dim>::active_cell_iterator cell =
                                                        tria.begin_active(),
                                                      endc = tria.end();
    for (; cell != endc; ++cell)
      {
        fe_val.reinit(cell);

        for (unsigned int d = 0; d < dim; ++d)
          for (unsigned int e = 0; e < dim; ++e)
            for (unsigned int f = 0; f < dim; ++f)
              deallog << fe_val.jacobian_grad(0)[d][e][f] << " ";
        deallog << std::endl;
      }
    deallog << std::endl;
  }
  {
    deallog << dim << "D Jacobian pushed forward gradients:" << std::endl;
    FEValues<dim>                                     fe_val(mapping,
                         dummy,
                         quad,
                         update_jacobian_pushed_forward_grads);
    typename Triangulation<dim>::active_cell_iterator cell =
                                                        tria.begin_active(),
                                                      endc = tria.end();
    for (; cell != endc; ++cell)
      {
        fe_val.reinit(cell);

        for (unsigned int d = 0; d < dim; ++d)
          for (unsigned int e = 0; e < dim; ++e)
            for (unsigned int f = 0; f < dim; ++f)
              deallog << fe_val.jacobian_pushed_forward_grad(0)[d][e][f] << " ";
        deallog << std::endl;
      }
    deallog << std::endl;
  }
  {
    deallog << dim << "D Jacobian hessians:" << std::endl;
    FEValues<dim> fe_val(mapping, dummy, quad, update_jacobian_2nd_derivatives);
    typename Triangulation<dim>::active_cell_iterator cell =
                                                        tria.begin_active(),
                                                      endc = tria.end();
    for (; cell != endc; ++cell)
      {
        fe_val.reinit(cell);

        for (unsigned int d = 0; d < dim; ++d)
          for (unsigned int e = 0; e < dim; ++e)
            for (unsigned int f = 0; f < dim; ++f)
              for (unsigned int g = 0; g < dim; ++g)
                deallog << fe_val.jacobian_2nd_derivative(0)[d][e][f][g] << " ";
        deallog << std::endl;
      }
    deallog << std::endl;
  }
  {
    deallog << dim << "D Jacobian pushed forward hessians:" << std::endl;
    FEValues<dim>                                     fe_val(mapping,
                         dummy,
                         quad,
                         update_jacobian_pushed_forward_2nd_derivatives);
    typename Triangulation<dim>::active_cell_iterator cell =
                                                        tria.begin_active(),
                                                      endc = tria.end();
    for (; cell != endc; ++cell)
      {
        fe_val.reinit(cell);

        for (unsigned int d = 0; d < dim; ++d)
          for (unsigned int e = 0; e < dim; ++e)
            for (unsigned int f = 0; f < dim; ++f)
              for (unsigned int g = 0; g < dim; ++g)
                deallog << fe_val.jacobian_pushed_forward_2nd_derivative(
                             0)[d][e][f][g]
                        << " ";
        deallog << std::endl;
      }
    deallog << std::endl;
  }
  {
    deallog << dim << "D Jacobian hessian gradients:" << std::endl;
    FEValues<dim> fe_val(mapping, dummy, quad, update_jacobian_3rd_derivatives);
    typename Triangulation<dim>::active_cell_iterator cell =
                                                        tria.begin_active(),
                                                      endc = tria.end();
    for (; cell != endc; ++cell)
      {
        fe_val.reinit(cell);

        for (unsigned int d = 0; d < dim; ++d)
          for (unsigned int e = 0; e < dim; ++e)
            for (unsigned int f = 0; f < dim; ++f)
              for (unsigned int g = 0; g < dim; ++g)
                for (unsigned int h = 0; h < dim; ++h)
                  deallog << fe_val.jacobian_3rd_derivative(0)[d][e][f][g][h]
                          << " ";
        deallog << std::endl;
      }
    deallog << std::endl;
  }
  {
    deallog << dim
            << "D Jacobian pushed forward hessian gradients:" << std::endl;
    FEValues<dim>                                     fe_val(mapping,
                         dummy,
                         quad,
                         update_jacobian_pushed_forward_3rd_derivatives);
    typename Triangulation<dim>::active_cell_iterator cell =
                                                        tria.begin_active(),
                                                      endc = tria.end();
    for (; cell != endc; ++cell)
      {
        fe_val.reinit(cell);

        for (unsigned int d = 0; d < dim; ++d)
          for (unsigned int e = 0; e < dim; ++e)
            for (unsigned int f = 0; f < dim; ++f)
              for (unsigned int g = 0; g < dim; ++g)
                for (unsigned int h = 0; h < dim; ++h)
                  deallog << fe_val.jacobian_pushed_forward_3rd_derivative(
                               0)[d][e][f][g][h]
                          << " ";
        deallog << std::endl;
      }
    deallog << std::endl;
  }
}



template <int dim>
void
test()
{
  Triangulation<dim> tria;
  GridGenerator::hyper_ball(tria);
  static const SphericalManifold<dim> manifold;
  tria.set_all_manifold_ids_on_boundary(0);
  tria.set_manifold(0, manifold);

  {
    deallog << "========== MappingQ ==========" << std::endl;
    MappingQ<dim> mapping(5);
    do_test(tria, mapping);
  }

  {
    deallog << "========== MappingFEField ==========" << std::endl;
    FESystem<dim>   fe_euler(FE_Q<dim>(QGaussLobatto<1>(4)), dim);
    DoFHandler<dim> map_dh(tria);
    map_dh.distribute_dofs(fe_euler);

    Vector<double> euler_vec(map_dh.n_dofs());
    VectorTools::get_position_vector(map_dh, euler_vec);

    MappingFEField<dim> mapping(map_dh, euler_vec);
    do_test(tria, mapping);
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
