// ---------------------------------------------------------------------
//
// Copyright (C) 2011 - 2018 by the deal.II authors
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

// Test interior mapping of flat manifold

#include "../tests.h"
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>


// all include files you need here
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/grid_out.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/mapping_q.h>


// Helper function
template <int dim, int spacedim>
void
test(unsigned int ref, const MappingQ<dim> &mapping)
{
  deallog << "Testing dim=" << dim
          << ", spacedim="<< spacedim << std::endl;

  Triangulation<dim,spacedim> tria;
  GridGenerator::hyper_cube (tria);
  tria.refine_global(ref);

  typename Triangulation<dim,spacedim>::active_cell_iterator
  cell;

  QGauss<dim> quadrature(4);
  FE_Q<dim> fe(2);

  FEValues<dim> fe_values (mapping, fe, quadrature,
                           update_gradients | update_values | update_quadrature_points  |
                           update_JxW_values);

  for (cell=tria.begin_active(); cell!=tria.end(); ++cell)
    {

      // check that FlatManifold returns the middle of the cell.
      deallog << "Cell: " << cell << std::endl;

      fe_values.reinit(cell);

      deallog << "  center: " << cell->center() << std::endl;
      for (unsigned int q=0; q<quadrature.size(); ++q)
        {

          deallog << "  JxW(" << q << "): " << fe_values.JxW(q) << std::endl;
          deallog << "  p("<<q<<"): " << fe_values.quadrature_point(q) << std::endl;
          for (unsigned int i=0; i<fe_values.dofs_per_cell; ++i)
            deallog << "  shape " << i << "," <<q<<": " << fe_values.shape_value(i,q)
                    << " " << fe_values.shape_grad(i,q) << std::endl;

        }

      if (cell->get_manifold().get_new_point_on_cell(cell).distance(cell->center()) > 1e-6)
        {
          deallog << "Default manifold: " << cell->get_manifold().get_new_point_on_cell(cell) << std::endl;
          deallog << "Center of cell  : " << cell->center() << std::endl;
        }
    }

  deallog << "OK" <<std::endl<<std::endl;

}

int
main ()
{
  initlog();

  test<2,2>(2, MappingQ<2>(4, false));
  test<2,2>(2, MappingQ<2>(4, true));


  return 0;
}
