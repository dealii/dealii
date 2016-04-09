// ---------------------------------------------------------------------
//
// Copyright (C) 2007 - 2015 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------



// test the FEValues views and extractor classes. these tests use a primitive
// finite element and scalar extractors

#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/base/function.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/lac/vector.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q1.h>

#include <fstream>



template<int dim>
void test (const Triangulation<dim> &tr,
           const FiniteElement<dim> &fe)
{
  DoFHandler<dim> dof(tr);
  dof.distribute_dofs(fe);

  deallog << "FE=" << fe.get_name()
          << std::endl;

  const QGauss<dim> quadrature(2);
  FEValues<dim> fe_values (fe, quadrature,
                           update_values | update_gradients | update_hessians);
  fe_values.reinit (dof.begin_active());

  for (unsigned int c=0; c<fe.n_components(); ++c)
    {
      FEValuesExtractors::Scalar single_component (c);

      for (unsigned int i=0; i<fe_values.dofs_per_cell; ++i)
        for (unsigned int q=0; q<fe_values.n_quadrature_points; ++q)
          {
            deallog << "i=" << i << ", q=" << q << std::endl;
            deallog << "   "
                    << fe_values[single_component].value (i,q) << ' ';
            for (unsigned int k=0; k<dim; ++k)
              deallog << fe_values[single_component].gradient (i,q) << ' ';
            deallog << std::endl;
            for (unsigned int k=0; k<dim; ++k)
              for (unsigned int l=0; l<dim; ++l)
                deallog << fe_values[single_component].hessian (i,q)[k][l] << std::endl;

            Assert (fe_values[single_component].value (i,q)
                    ==
                    fe_values.shape_value_component (i,q,c),
                    ExcInternalError());

            Assert (fe_values[single_component].gradient (i,q)
                    ==
                    fe_values.shape_grad_component (i,q,c),
                    ExcInternalError());

            Assert (fe_values[single_component].hessian (i,q)
                    ==
                    fe_values.shape_hessian_component (i,q,c),
                    ExcInternalError());
          }
    }
}



template<int dim>
void test_hyper_sphere()
{
  Triangulation<dim> tr;
  GridGenerator::hyper_ball(tr);

  static const HyperBallBoundary<dim> boundary;
  tr.set_boundary (0, boundary);

  FESystem<dim> fe (FE_Q<dim>(1), 1,
                    FE_Q<dim>(2), 2,
                    FE_DGQArbitraryNodes<dim>(QIterated<1>(QTrapez<1>(),3)), dim);
  test(tr, fe);
}


int main()
{
  std::ofstream logfile ("output");
  deallog << std::setprecision (2);

  deallog.attach(logfile);
  deallog.threshold_double(1.e-7);

  test_hyper_sphere<2>();
  test_hyper_sphere<3>();
}
