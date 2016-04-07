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

// Check that MappingManifold and MappingQ1 are the same thing on a
// FlatManifold. Test on the quadrature points.

#include "../tests.h"

#include <deal.II/base/utilities.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/fe/mapping_q_generic.h>
#include <deal.II/fe/mapping_manifold.h>
#include <deal.II/grid/manifold_lib.h>

template<int dim, int spacedim>
void test()
{
  deallog << "dim=" << dim << ", spacedim=" << spacedim << std::endl;

  Triangulation<dim, spacedim>   triangulation;

  GridGenerator::hyper_cube (triangulation, 2.0, 3.0);


  const QGauss<dim> quadrature(5);

  std::vector<Point<dim> > q_points = quadrature.get_points();

  MappingManifold<dim,spacedim> map_manifold;
  MappingQGeneric<dim,spacedim> map_q1(1);

  typename  Triangulation<dim,spacedim>::active_cell_iterator
  cell = triangulation.begin_active(),
  endc = triangulation.end();
  for (; cell!=endc; ++cell)
    {
      for (unsigned int i=0; i<q_points.size(); ++i)
        {
          const Point<spacedim> pq = map_manifold.transform_unit_to_real_cell(cell,q_points[i]);
          const Point<spacedim> pq1 = map_q1.transform_unit_to_real_cell(cell, q_points[i]);
          if (pq.distance(pq1) > 1e-10)
            {
              deallog << "Expected: " << pq << ", got: "
                      << pq1 << std::endl;
            }
        }
    }
  deallog << "OK" << std::endl;
}


int
main()
{
  initlog();

  test<1,1>();
  test<2,2>();
  test<3,3>();

  test<1,2>();
  test<1,3>();
  test<2,3>();

  return 0;
}



