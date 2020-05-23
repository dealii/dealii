// ---------------------------------------------------------------------
//
// Copyright (C) 2004 - 2020 by the deal.II authors
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


// check the correctness of fe_values.shape_gradient for FE_BernardiRaugel
// this test is almost a verbatim copy of the fe_rt.cc file

#include <deal.II/base/function.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_bernardi_raugel.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_nedelec.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/lac/vector.h>

#include <sstream>

#include "../tests.h"



template <int dim>
Tensor<1, dim>
ones()
{
  Tensor<1, dim> result;
  for (unsigned int i = 0; i < dim; ++i)
    result[i] = 1.0;
  return result;
}

template <int dim>
void
test(const Triangulation<dim> &tr,
     const FiniteElement<dim> &fe,
     const double              tolerance)
{
  DoFHandler<dim> dof(tr);
  dof.distribute_dofs(fe);

  std::stringstream ss;

  deallog << "FE=" << fe.get_name() << std::endl;

  const QGauss<dim> quadrature(4);
  FEValues<dim>     fe_values(fe,
                          quadrature,
                          update_gradients | update_quadrature_points |
                            update_JxW_values);

  for (typename DoFHandler<dim>::active_cell_iterator cell = dof.begin_active();
       cell != dof.end();
       ++cell)
    {
      fe_values.reinit(cell);

      deallog << "Cell nodes:" << std::endl;
      for (const unsigned int i : GeometryInfo<dim>::vertex_indices())
        {
          deallog << i << ": ( ";
          for (unsigned int d = 0; d < dim; ++d)
            deallog << cell->vertex(i)[d] << " ";
          deallog << ")" << std::endl;
        }

      for (unsigned int c = 0; c < fe.n_components(); ++c)
        {
          FEValuesExtractors::Scalar single_component(c);

          for (unsigned int i = 0; i < fe_values.dofs_per_cell; ++i)
            {
              ss << "component=" << c << ", dof=" << i << std::endl;

              Tensor<1, dim> bulk_integral;
              for (const auto q : fe_values.quadrature_point_indices())
                {
                  bulk_integral += fe_values[single_component].gradient(i, q) *
                                   fe_values.JxW(q);
                }

              deallog << "    bulk integral=" << bulk_integral << std::endl;
            }
        }
    }
}



template <int dim>
void
test_hyper_cube(const double tolerance)
{
  Triangulation<dim> tr;
  GridGenerator::hyper_cube(tr);

  typename Triangulation<dim>::active_cell_iterator cell = tr.begin_active();
  for (const unsigned int i : GeometryInfo<dim>::vertex_indices())
    {
      Point<dim> &point = cell->vertex(i);
      if (std::abs(point(dim - 1) - 1.0) < 1e-5)
        point(dim - 1) += 0.15;
    }

  FE_BernardiRaugel<dim> fe(1);
  test(tr, fe, tolerance);
}


int
main()
{
  initlog();
  deallog << std::setprecision(3);

  test_hyper_cube<2>(1e-6);
  test_hyper_cube<3>(1e-6);
}
