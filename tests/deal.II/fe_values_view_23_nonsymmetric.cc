// ---------------------------------------------------------------------
//
// Copyright (C) 2007 - 2013 by the deal.II authors
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



// there was a bug in getting the divergence of shape functions for the
// Tensor extractors. test that it is fixed by comparing with
// get_function_divergences

#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/base/function.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/lac/vector.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q1.h>

#include <fstream>



template<int dim>
void test (const Triangulation<dim> &tr,
           const FiniteElement<dim> &fe)
{
  deallog << "FE=" << fe.get_name()
          << std::endl;

  DoFHandler<dim> dof(tr);
  dof.distribute_dofs(fe);

  Vector<double> fe_function(dof.n_dofs());
  for (unsigned int i=0; i<dof.n_dofs(); ++i)
    fe_function(i) = (i+1)*(i+2);

  const QGauss<dim> quadrature(2);
  FEValues<dim> fe_values (fe, quadrature,
                           update_values | update_gradients);
  fe_values.reinit (dof.begin_active());

  // let the FEValues object compute the
  // divergences at quadrature points
  std::vector<Tensor<1,dim> > divergences (quadrature.size());
  FEValuesExtractors::Tensor<2> extractor(0);
  fe_values[extractor]
  .get_function_divergences (fe_function, divergences);

  // now do the same "by hand"
  std::vector<types::global_dof_index> local_dof_indices (fe.dofs_per_cell);
  dof.begin_active()->get_dof_indices (local_dof_indices);

  for (unsigned int q=0; q<quadrature.size(); ++q)
    {
      Tensor<1,dim> div_alt;
      for (unsigned int i=0; i<fe.dofs_per_cell; ++i)
        div_alt += fe_values[extractor].divergence (i,q) *
                   fe_function(local_dof_indices[i]);

      deallog << "q_point=" << q << std::endl
              << "   method 1: " << divergences[q] << std::endl
              << "   method 2: " << div_alt << std::endl
              << std::endl;
      Assert ((divergences[q] - div_alt).norm() <= divergences[q].norm(),
              ExcInternalError());
    }
}



template<int dim>
void test_hyper_sphere()
{
  Triangulation<dim> tr;
  GridGenerator::hyper_ball(tr);

  static const HyperBallBoundary<dim> boundary;
  tr.set_boundary (0, boundary);

  FESystem<dim> fe (FE_Q<dim>(1),
                    Tensor<2,dim>::n_independent_components);
  test(tr, fe);
}


int main()
{
  std::ofstream logfile ("output");
  deallog << std::setprecision (3);

  deallog.attach(logfile);
  deallog.depth_console (0);
  deallog.threshold_double(1.e-7);

  test_hyper_sphere<2>();
  test_hyper_sphere<3>();
}
