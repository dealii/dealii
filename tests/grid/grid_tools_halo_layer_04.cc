// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2015 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



#include <deal.II/dofs/dof_handler.h>

#include <deal.II/grid/filtered_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>

#include <deal.II/numerics/data_out.h>

#include "../tests.h"

template <int dim>
void
write_active_fe_index_to_file(const DoFHandler<dim> &dof_handler)
{
  int                                            count = 0;
  typename DoFHandler<dim>::active_cell_iterator cell =
                                                   dof_handler.begin_active(),
                                                 endc = dof_handler.end();
  for (; cell != endc; ++cell, ++count)
    {
      deallog << count << ' ' << cell->active_fe_index() << std::endl;
    }
  deallog << std::endl;
}

template <int dim>
void
write_vtk(const DoFHandler<dim> &dof_handler, const std::string filename)
{
  Vector<double> active_fe_index(
    dof_handler.get_triangulation().n_active_cells());
  int                                            count = 0;
  typename DoFHandler<dim>::active_cell_iterator cell =
                                                   dof_handler.begin_active(),
                                                 endc = dof_handler.end();
  for (; cell != endc; ++cell, ++count)
    {
      active_fe_index[count] = cell->active_fe_index();
    }

  const std::vector<DataComponentInterpretation::DataComponentInterpretation>
    data_component_interpretation(
      1, DataComponentInterpretation::component_is_scalar);
  const std::vector<std::string> data_names(1, "active_fe_index");

  DataOut<dim> data_out;
  data_out.attach_dof_handler(dof_handler);
  data_out.add_data_vector(active_fe_index,
                           data_names,
                           DataOut<dim>::type_cell_data,
                           data_component_interpretation);
  data_out.build_patches();

  std::ofstream output(filename);
  data_out.write_vtk(output);
}


template <int dim>
void
test()
{
  deallog << "dim = " << dim << std::endl;

  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria);
  tria.refine_global(2);
  DoFHandler<dim> dof_handler(tria);

  using cell_iterator = typename DoFHandler<dim>::active_cell_iterator;

  // Mark a small block at the corner of the hypercube
  cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
  for (; cell != endc; ++cell)
    {
      bool mark = true;
      for (unsigned int d = 0; d < dim; ++d)
        if (cell->center()[d] > 0.5)
          {
            mark = false;
            break;
          }

      if (mark == true)
        cell->set_active_fe_index(2);
      else
        cell->set_active_fe_index(1);
    }

  deallog << "Grid without halo:" << std::endl;
  write_active_fe_index_to_file(dof_handler);
  // Write to file to visually check result
  {
    const std::string filename =
      "grid_no_halo_" + Utilities::int_to_string(dim) + "d.vtk";
    write_vtk(dof_handler, filename.c_str());
  }

  // Compute a halo layer around active FE index 2 and set it to active FE index
  // 3
  std::function<bool(const cell_iterator &)> predicate =
    IteratorFilters::ActiveFEIndexEqualTo(2, true);
  std::vector<cell_iterator> active_halo_layer =
    GridTools::compute_active_cell_halo_layer(dof_handler, predicate);
  AssertThrow(active_halo_layer.size() > 0, ExcMessage("No halo layer found."));
  for (typename std::vector<cell_iterator>::iterator it =
         active_halo_layer.begin();
       it != active_halo_layer.end();
       ++it)
    {
      (*it)->set_active_fe_index(3);
    }

  deallog << "Grid with halo:" << std::endl;
  write_active_fe_index_to_file(dof_handler);
  // Write to file to visually check result
  {
    const std::string filename =
      "grid_with_halo_" + Utilities::int_to_string(dim) + "d.vtk";
    write_vtk(dof_handler, filename.c_str());
  }
}


int
main()
{
  initlog();

  test<2>();
  test<3>();

  return 0;
}
