// ---------------------------------------------------------------------
//
// Copyright (C) 2007 - 2020 by the deal.II authors
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



// test FEInterfaceViews and extractor classes. these tests use a primitive
// finite element and vector extractors. To simplify the output, we skip writing
// results of high-order tensors and zeros.

#include <deal.II/base/function.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_interface_values.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/manifold_lib.h>

#include <deal.II/lac/vector.h>

#include <iostream>

#include "../tests.h"

template <int dim>
void
test(const Triangulation<dim> &tr, const FiniteElement<dim> &fe)
{
  DoFHandler<dim> dof(tr);
  dof.distribute_dofs(fe);

  deallog << "FE=" << fe.get_name() << std::endl;

  const QGauss<dim - 1>  face_quadrature(2);
  FEInterfaceValues<dim> fe_iv(fe,
                               face_quadrature,
                               update_values | update_gradients |
                                 update_quadrature_points | update_hessians |
                                 update_3rd_derivatives);
  auto                   cell = dof.begin_active();


  for (const auto f : GeometryInfo<dim>::face_indices())
    if (!cell->face(f)->at_boundary())
      {
        fe_iv.reinit(cell,
                     f,
                     numbers::invalid_unsigned_int,
                     cell->neighbor(f),
                     cell->neighbor_of_neighbor(f),
                     numbers::invalid_unsigned_int);
        const auto &       q_points    = fe_iv.get_quadrature_points();
        const unsigned int n_q_points  = q_points.size();
        const unsigned int n_dofs_face = fe_iv.n_current_interface_dofs();

        for (unsigned int c = 0; c < fe.n_components(); ++c)
          {
            // use a vector extractor if there
            // are sufficiently many components
            // left after the current component
            // 'c'
            if (c + dim <= fe.n_components())
              {
                FEValuesExtractors::Vector vec_components(c);
                for (unsigned int i = 0; i < n_dofs_face; ++i)
                  for (unsigned int q = 0; q < n_q_points; ++q)
                    {
                      if (fe_iv[vec_components].value(true, i, q).norm() != 0)
                        deallog << fe_iv[vec_components].value(true, i, q)
                                << "  ";
                      if (fe_iv[vec_components].value(false, i, q).norm() != 0)
                        deallog << fe_iv[vec_components].value(false, i, q)
                                << "  ";
                      if (fe_iv[vec_components].jump(i, q).norm() != 0)
                        deallog << fe_iv[vec_components].jump(i, q) << "  ";
                      if (fe_iv[vec_components].average(i, q).norm() != 0)
                        deallog << fe_iv[vec_components].average(i, q) << "  ";
                      if (fe_iv[vec_components].jump_gradient(i, q).norm() != 0)
                        deallog << fe_iv[vec_components].jump_gradient(i, q)
                                << "  ";
                      if (fe_iv[vec_components].average_gradient(i, q).norm() !=
                          0)
                        deallog << fe_iv[vec_components].average_gradient(i, q)
                                << std::endl;

                      for (unsigned int d = 0; d < dim; ++d)
                        {
                          Assert(fe_iv[vec_components].value(true, i, q)[d] ==
                                   fe_iv.shape_value(true, i, q, c + d),
                                 ExcInternalError());
                          Assert(fe_iv[vec_components].value(false, i, q)[d] ==
                                   fe_iv.shape_value(false, i, q, c + d),
                                 ExcInternalError());
                          Assert(fe_iv[vec_components].jump(i, q)[d] ==
                                   fe_iv.jump(i, q, c + d),
                                 ExcInternalError());
                          Assert(fe_iv[vec_components].average(i, q)[d] ==
                                   fe_iv.average(i, q, c + d),
                                 ExcInternalError());
                          Assert(fe_iv[vec_components].average_hessian(i,
                                                                       q)[d] ==
                                   fe_iv.average_hessian(i, q, c + d),
                                 ExcInternalError());
                          Assert(fe_iv[vec_components].jump_hessian(i, q)[d] ==
                                   fe_iv.jump_hessian(i, q, c + d),
                                 ExcInternalError());

                          Assert(fe_iv[vec_components].jump_3rd_derivative(
                                   i, q)[d] ==
                                   fe_iv.jump_3rd_derivative(i, q, c + d),
                                 ExcInternalError());
                        }
                    }
              }
          }
        return;
      } // at boundary
}



template <int dim>
void
test_hyper_sphere()
{
  Triangulation<dim> tr;
  GridGenerator::hyper_ball(tr);

  static const SphericalManifold<dim> boundary;
  tr.set_manifold(0, boundary);

  FESystem<dim> fe(FE_Q<dim>(1),
                   1,
                   FE_DGQArbitraryNodes<dim>(QGauss<1>(2)),
                   dim);
  test(tr, fe);
}


int
main()
{
  initlog();
  deallog << std::setprecision(2);

  test_hyper_sphere<2>();
  test_hyper_sphere<3>();
}
