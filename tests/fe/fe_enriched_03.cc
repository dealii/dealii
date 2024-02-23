// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2016 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Check FEFaceValues::shape_value()

#include <deal.II/base/function.h>
#include <deal.II/base/utilities.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_enriched.h>
#include <deal.II/fe/fe_nothing.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>

#include <deal.II/hp/fe_collection.h>
#include <deal.II/hp/fe_values.h>
#include <deal.II/hp/q_collection.h>

#include <deal.II/lac/vector.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/data_postprocessor.h>

#include <iostream>

#include "../tests.h"


template <int dim>
class EnrichmentFunction : public Function<dim>
{
public:
  EnrichmentFunction()
    : Function<dim>(1)
  {}

  virtual double
  value(const Point<dim> &point, const unsigned int component = 0) const
  {
    return std::exp(-point.norm());
  }

  virtual Tensor<1, dim>
  gradient(const Point<dim> &point, const unsigned int component = 0) const
  {
    Tensor<1, dim> res = point;
    Assert(point.norm() > 0,
           dealii::ExcMessage("gradient is not defined at zero"));
    res *= -value(point) / point.norm();
    return res;
  }
};


template <int dim>
void
test2()
{
  deallog << "FEFaceValues" << std::endl;
  deallog << "for same underlying FEs: f(qp) * N_{fe}(qp) == N_{pou}(qp)"
          << std::endl;

  Triangulation<dim> triangulation;
  DoFHandler<dim>    dof_handler(triangulation);

  EnrichmentFunction<dim> function;
  FE_Enriched<dim>        fe(FE_Q<dim>(1), FE_Q<dim>(1), &function);

  GridGenerator::hyper_cube(triangulation);
  dof_handler.distribute_dofs(fe);

  QGauss<dim - 1>   quadrature(1);
  FEFaceValues<dim> fe_face_values(
    fe, quadrature, update_values | update_gradients | update_JxW_values);

  typename DoFHandler<dim>::active_cell_iterator cell =
                                                   dof_handler.begin_active(),
                                                 endc = dof_handler.end();
  for (; cell != endc; ++cell)
    for (const unsigned int face : GeometryInfo<dim>::face_indices())
      {
        fe_face_values.reinit(cell, face);
        const unsigned int                     n_q_points = quadrature.size();
        const unsigned int                     dofs_per_cell = fe.dofs_per_cell;
        const std::vector<dealii::Point<dim>> &q_points =
          fe_face_values.get_quadrature_points();

        for (unsigned int i = 0; i < dofs_per_cell; ++i)
          for (const auto q_point : fe_face_values.quadrature_point_indices())
            deallog << "dof=" << i << " qp=" << q_points[q_point]
                    << " f(qp)=" << function.value(q_points[q_point])
                    << " N(qp)=" << fe_face_values.shape_value(i, q_point)
                    << std::endl;
      }

  dof_handler.clear();
}


int
main(int argc, char **argv)
{
  initlog();
  deallog << std::setprecision(4) << std::fixed;
  deallog.depth_console(0);

  try
    {
      test2<3>();
    }
  catch (const std::exception &exc)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Exception on processing: " << std::endl
                << exc.what() << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;

      return 1;
    }
  catch (...)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Unknown exception!" << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    };
}
