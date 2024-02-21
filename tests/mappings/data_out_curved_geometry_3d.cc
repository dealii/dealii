// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2016 - 2021 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// Like data_out_curved_geometry.cc, but in 3d

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/numerics/data_out.h>

#include "../tests.h"



int
main()
{
  initlog();

  const CylindricalManifold<3> boundary_description;

  for (unsigned int degree = 1; degree <= 4; ++degree)
    {
      deallog << "===== Mapping degree " << degree << std::endl;

      Triangulation<3> triangulation;
      GridGenerator::cylinder(triangulation);

      const MappingQ<3> mapping(degree);
      const FE_Q<3>     dummy_fe(1);
      DoFHandler<3>     dof_handler(triangulation);

      for (Triangulation<3>::active_cell_iterator cell =
             triangulation.begin_active();
           cell != triangulation.end();
           ++cell)
        {
          if (cell->center()[1] * cell->center()[1] +
                cell->center()[2] * cell->center()[2] <
              1.e-5)
            cell->set_material_id(1);
          else
            cell->set_material_id(0);
        }

      Triangulation<3>::active_cell_iterator cell =
                                               triangulation.begin_active(),
                                             endc = triangulation.end();

      for (; cell != endc; ++cell)
        if (cell->material_id() != 1)
          cell->set_all_manifold_ids(100);

      triangulation.set_manifold(100, boundary_description);

      dof_handler.distribute_dofs(dummy_fe);

      Vector<double> dummy(triangulation.n_active_cells());
      for (Triangulation<3>::active_cell_iterator cell =
             triangulation.begin_active();
           cell != triangulation.end();
           ++cell)
        dummy(cell->active_cell_index()) = cell->material_id();

      DataOut<3> data_out;
      data_out.attach_dof_handler(dof_handler);
      data_out.add_data_vector(dummy, "dummy");
      data_out.build_patches(mapping, 2, DataOut<3>::curved_inner_cells);
      data_out.write_vtk(deallog.get_file_stream());
    }
}
