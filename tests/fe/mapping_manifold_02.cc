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

  MappingManifold<dim,spacedim> mapping;

  FEValues<dim,spacedim> fe_values (mapping, fe, quad,
                                    update_quadrature_points);

  for (typename Triangulation<dim,spacedim>::active_cell_iterator cell= triangulation.begin_active();
       cell != triangulation.end(); ++cell)
    {
      fe_values.reinit(cell);
      std::vector<Point<spacedim> > fev_qp = fe_values.get_quadrature_points();
      for (unsigned int q=0; q<fev_qp.size(); ++q)
        {
          const Point<spacedim> pq = mapping.transform_unit_to_real_cell(cell,quad.point(q));

          if (pq.distance(fev_qp[q]) > 1e-10)
            {
              deallog << "Expected: " << pq << ", got: "
                      << fev_qp[q] << std::endl;
            }
        }
    }
  deallog << "OK" << std::endl;
}


int
main()
{
  initlog();

  test<2,2>();
  test<3,3>();

  return 0;
}



