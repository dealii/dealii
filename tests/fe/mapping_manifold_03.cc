// ---------------------------------------------------------------------
//
// Copyright (C) 2016 by the deal.II authors
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


// Check that the quadrature points obtained through transform_unit_to_real
// are the same of those given by FEValues


#include "../tests.h"

#include <deal.II/base/utilities.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/mapping_q_generic.h>
#include <deal.II/fe/mapping_manifold.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>

template<int dim, int spacedim>
void test()
{
  deallog << "dim=" << dim << ", spacedim=" << spacedim << std::endl;

  Triangulation<dim, spacedim>   triangulation;

  Point<spacedim> center;
  for (unsigned int i=0; i<spacedim; ++i)
    center[i] = 5+i;

  const double inner_radius = 0.5,
               outer_radius = 1.0;
  GridGenerator::hyper_shell (triangulation,
                              center, inner_radius, outer_radius);


  static  const SphericalManifold<dim> manifold(center);

  triangulation.set_all_manifold_ids(0);
  triangulation.set_manifold (0, manifold);

  const QGauss<spacedim> quad(3);

  FE_Q<dim,spacedim> fe(1);

  DoFHandler<dim,spacedim> dof(triangulation);
  dof.distribute_dofs(fe);

  MappingManifold<dim,spacedim> mapping_manifold;
  MappingQGeneric<dim,spacedim> mapping_q(1);

  FEValues<dim,spacedim> fe_values_mapping (mapping_manifold,
                                            fe, quad,
                                            update_jacobians);

  FEValues<dim,spacedim> fe_values_q (mapping_q,
                                      fe, quad,
                                      update_jacobians);

  typename Triangulation<dim,spacedim>::active_cell_iterator cell= triangulation.begin_active();

  fe_values_mapping.reinit(cell);
  fe_values_q.reinit(cell);
  std::vector<DerivativeForm<1,dim,spacedim> > jac_from_mapping_manifold = fe_values_mapping.get_jacobians();

  std::vector<DerivativeForm<1,dim,spacedim> > jac_from_mapping_q = fe_values_q.get_jacobians();

  AssertThrow(jac_from_mapping_q.size() == jac_from_mapping_manifold.size(), ExcInternalError());

  for (unsigned int q=0; q<jac_from_mapping_q.size(); ++q)
    {

      deallog << "Jacobian from mapping manifold at point "<< q << std::endl;
      for (unsigned int d=0; d<dim; ++d)
        deallog << jac_from_mapping_manifold[q][d] << std::endl;

      deallog << "Jacobian from mapping q at point "<< q << std::endl;
      for (unsigned int d=0; d<dim; ++d)
        deallog << jac_from_mapping_q[q][d] << std::endl;

    }

  deallog << "OK" << std::endl;
}


int
main()
{
  std::ofstream logfile ("output");
  deallog.attach(logfile);
  deallog.threshold_double(1.e-10);

  test<2,2>();
  test<3,3>();

  return 0;
}



