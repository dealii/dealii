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



#include "../tests.h"

#include <deal.II/base/utilities.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/fe/mapping_q_generic.h>
#include <deal.II/fe/mapping_manifold.h>
#include <deal.II/grid/manifold_lib.h>

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
  triangulation.refine_global (0);

  const QIterated<QTrapez<spacedim> > quadrature(7);

  std::vector<Point<spacedim> > q_points = quadrature.get_points();


  MappingManifold<dim,spacedim> map_manifold;
  MappingQGeneric<dim,spacedim> map_q1(1);


  cell = triangulation.begin_active();

  std::ofstream two_be_plotted ("points.txt");
  for (; cell!=endc; ++cell)
    {
      for (unsigned int i=0; i<q_points.size(); ++i)
        {
          // const Point<spacedim> p = map.transform_unit_to_real_cell(cell,q_points[i]);
          // deallog << "mapping manifold ---> " << p.distance() << std::endl;

          const Point<spacedim> pq = map_manifold.transform_unit_to_real_cell(cell,q_points[i]);
          deallog << "mapping q1       ---> " << pq.distance(origin) << std::endl;
          for (unsigned int d=0; d<spacedim; ++d)
            two_be_plotted << pq[d] << "  ";
          two_be_plotted << std::endl;
        }
    }
}


int
main()
{
  std::ofstream logfile ("output");
  deallog.attach(logfile);
  deallog.threshold_double(1.e-10);


  //  test<1,1>();
  test<2,2>();
  // test<3,3>();

  // test<1,2>();
  // test<1,3>();
  // test<2,3>();

  return 0;
}



