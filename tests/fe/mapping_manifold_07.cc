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

// Check MappingManifold on faces. Create a sphere around a one cell
// Triangulation, and make sure all points are on the circle.

#include "../tests.h"

#include <deal.II/base/utilities.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/fe/mapping_q_generic.h>
#include <deal.II/fe/mapping_manifold.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/grid/manifold_lib.h>

template<int dim, int spacedim>
void test()
{
  deallog << "dim=" << dim << ", spacedim=" << spacedim << std::endl;

  Triangulation<dim, spacedim>   triangulation;

  GridGenerator::hyper_cube (triangulation, 2.0, 3.0);

  typename Triangulation<dim,spacedim>::active_cell_iterator
  cell = triangulation.begin_active();

  // Center and radius of the Ball
  Point<spacedim> center = cell->center();
  double radius = center.distance(cell->vertex(0));

  static  const SphericalManifold<dim,spacedim> manifold(cell->center());

  triangulation.set_all_manifold_ids_on_boundary(0);
  triangulation.set_manifold (0, manifold);

  MappingManifold<dim,spacedim> map_manifold;
  FE_Q<dim,spacedim> fe(1);
  const QGauss<dim-1> quad(5);

  FEFaceValues<dim,spacedim> fe_v(map_manifold, fe, quad,
                                  update_quadrature_points);



  for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
    {
      fe_v.reinit(cell, f);
      const std::vector<Point<spacedim> > &qps = fe_v.get_quadrature_points();
      for (unsigned int i=0; i<qps.size(); ++i)
        if ( std::abs(qps[i].distance(center) - radius) > 1e-10)
          {
            deallog << "Expected radius: " << radius << ", got: "
                    <<  qps[i].distance(center)
                    << ", on point " << qps[i] << std::endl;
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

  test<2,2>();
  test<2,3>();

  test<3,3>();

  return 0;
}



