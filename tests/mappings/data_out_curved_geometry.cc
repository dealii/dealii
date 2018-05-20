// ---------------------------------------------------------------------
//
// Copyright (C) 2016 - 2018 by the deal.II authors
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

// MappingQ::clone forgot to copy the flag that describes whether to
// also use curved geometries in interior cells. This could lead to
// wrong output in DataOut, which puts the mapping into a
// hp::MappingCollection, which requires cloning the mapping.
//
// based on a testcase by Sebastian Gonzalez-Pintor

#include "../tests.h"
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/numerics/data_out.h>

int
main()
{
  initlog();

  const Point<2>             center(0, 0);
  const SphericalManifold<2> boundary_description(center);

  for(unsigned int degree = 1; degree <= 4; ++degree)
    {
      deallog << "===== Mapping degree " << degree << std::endl;

      Triangulation<2> triangulation;
      GridGenerator::hyper_ball(triangulation, Point<2>(0.0, 0.0), 1.0);

      const MappingQ<2> mapping(degree, true);
      const FE_Q<2>     dummy_fe(1);
      DoFHandler<2>     dof_handler(triangulation);

      for(Triangulation<2>::active_cell_iterator cell
          = triangulation.begin_active();
          cell != triangulation.end();
          ++cell)
        {
          if(cell->center().square() < 1.e-5)
            cell->set_material_id(1);
          else
            cell->set_material_id(0);
        }

      Triangulation<2>::active_cell_iterator cell
        = triangulation.begin_active(),
        endc = triangulation.end();

      for(; cell != endc; ++cell)
        if(cell->material_id() != 1)
          cell->set_all_manifold_ids(100);

      triangulation.set_manifold(100, boundary_description);

      dof_handler.distribute_dofs(dummy_fe);

      Vector<double> dummy(triangulation.n_active_cells());

      DataOut<2> data_out;
      data_out.attach_dof_handler(dof_handler);
      data_out.add_data_vector(dummy, "dummy");
      data_out.build_patches(mapping, 5, DataOut<2>::curved_inner_cells);
      data_out.write_vtk(deallog.get_file_stream());
    }
}
