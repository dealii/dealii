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

#include <fstream>

template<int dim, int spacedim>
void test()
{
  deallog << "dim=" << dim << ", spacedim=" << spacedim << std::endl;

  Triangulation<dim, spacedim>   triangulation;

  const Point<2> center (7,9);
  const Point<2> origin (0,0);
  const double inner_radius = 0.5,
               outer_radius = 1.0;
  GridGenerator::hyper_shell (triangulation,
                              center, inner_radius, outer_radius,
                              10);


  static  const SphericalManifold<dim> manifold(center);

  triangulation.set_all_manifold_ids_on_boundary(0);
  triangulation.set_manifold (0, manifold);
  typename  Triangulation<dim,spacedim>::active_cell_iterator
  cell = triangulation.begin_active(),
  endc = triangulation.end();
  for (; cell!=endc; ++cell)
    cell->set_all_manifold_ids (0);

  const QGauss<spacedim> quadrature(1);

  std::vector<Point<spacedim> > q_points_from_quadrature = quadrature.get_points();

  FE_Q<dim,spacedim> fe(1);

  DoFHandler<dim,spacedim> dof(triangulation);
  dof.distribute_dofs(fe);


  FEValues<dim,spacedim> fe_values (fe, quadrature,
                                    update_quadrature_points);


  MappingManifold<dim,spacedim> map_manifold;

  cell = triangulation.begin_active();

  for (; cell!=endc; ++cell)
    {
      fe_values.reinit(cell);
      std::vector<Point<spacedim> > q_points_from_fe_values = fe_values.get_quadrature_points();
      for (unsigned int q=0; q<q_points_from_fe_values.size(); ++q)
        {
          const Point<spacedim> pq = map_manifold.transform_unit_to_real_cell(cell,q_points_from_quadrature[q]);

          if (pq.distance(q_points_from_fe_values[q]) > 1e-10)
            {
              deallog << "Expected: " << pq << ", got: "
                      << q_points_from_fe_values[q] << std::endl;
            }
        }
    }
  deallog << "OK" << std::endl;
}


int
main()
{
  std::ofstream logfile ("output");
  deallog.attach(logfile);
  deallog.threshold_double(1.e-10);


  //  test<1,1>();
  test<2,2>();
  //  test<3,3>();

  // test<1,2>();
  // test<1,3>();
  // test<2,3>();

  return 0;
}



