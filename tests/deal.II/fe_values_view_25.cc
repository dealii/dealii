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



// like _24, but for a simpler mesh for which the output has been
// verified to be correct

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
    fe_function(i) = i+1;

  const QGauss<dim> quadrature(2);
  FEValues<dim> fe_values (fe, quadrature,
                           update_values | update_gradients | update_q_points);
  fe_values.reinit (dof.begin_active());

  // let the FEValues object compute the
  // divergences at quadrature points
  std::vector<Tensor<1,dim> > divergences (quadrature.size());
  FEValuesExtractors::SymmetricTensor<2> extractor(0);
  fe_values[extractor]
  .get_function_divergences (fe_function, divergences);

  // now do the same "by hand"
  std::vector<types::global_dof_index> local_dof_indices (fe.dofs_per_cell);
  dof.begin_active()->get_dof_indices (local_dof_indices);

  for (unsigned int i=0; i<fe.dofs_per_cell; ++i)
    {
      deallog << "i=" << i << std::endl;

      for (unsigned int q=0; q<quadrature.size(); ++q)
        deallog << "  q_point=" << fe_values.quadrature_point(q) << std::endl
                << "    value= " << fe_values[extractor].value (i,q) << std::endl
                << "    div= " << fe_values[extractor].divergence (i,q) << std::endl;
    }
}



template<int dim>
void test_hyper_cube()
{
  Triangulation<dim> tr;
  GridGenerator::hyper_cube(tr);

  FESystem<dim> fe (FE_Q<dim>(1),
                    SymmetricTensor<2,dim>::n_independent_components);
  test(tr, fe);
}


int main()
{
  std::ofstream logfile ("output");
  deallog << std::setprecision (3);

  deallog.attach(logfile);
  deallog.depth_console (0);
  deallog.threshold_double(1.e-7);

  test_hyper_cube<1>();
  test_hyper_cube<2>();
  test_hyper_cube<3>();
}
