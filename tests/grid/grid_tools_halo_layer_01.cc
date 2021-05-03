// ---------------------------------------------------------------------
//
// Copyright (C) 2001 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------



#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"

template <int dim>
bool
pred_mat_id(const typename Triangulation<dim>::active_cell_iterator &cell)
{
  return cell->material_id() == 2;
}

template <int dim>
void
write_mat_id_to_file(const Triangulation<dim> &tria)
{
  int                                               count = 0;
  typename Triangulation<dim>::active_cell_iterator cell  = tria.begin_active(),
                                                    endc  = tria.end();
  for (; cell != endc; ++cell, ++count)
    {
      deallog << count << " " << static_cast<int>(cell->material_id())
              << std::endl;
    }
  deallog << std::endl;
}


template <int dim>
void
test()
{
  deallog << "dim = " << dim << std::endl;

  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria);
  tria.refine_global(2);

  using cell_iterator = typename Triangulation<dim>::active_cell_iterator;

  // Mark a small block at the corner of the hypercube
  cell_iterator cell = tria.begin_active(), endc = tria.end();
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
        cell->set_material_id(2);
      else
        cell->set_material_id(1);
    }

  deallog << "Grid without halo:" << std::endl;
  write_mat_id_to_file(tria);
  // Write to file to visually check result
  {
    const std::string filename =
      "grid_no_halo_" + Utilities::int_to_string(dim) + "d.vtk";
    std::ofstream f(filename.c_str());
    GridOut().write_vtk(tria, f);
  }

  std::function<bool(const cell_iterator &)> predicate = pred_mat_id<dim>;

  // Compute a halo layer around material id 2 and set it to material id 3
  const std::vector<cell_iterator> active_halo_layer =
    GridTools::compute_active_cell_halo_layer(tria,
                                              predicate); // General predicate

  AssertThrow(active_halo_layer.size() > 0, ExcMessage("No halo layer found."));
  for (typename std::vector<cell_iterator>::const_iterator it =
         active_halo_layer.begin();
       it != active_halo_layer.end();
       ++it)
    {
      (*it)->set_material_id(3);
    }

  deallog << "Grid with halo:" << std::endl;
  write_mat_id_to_file(tria);
  // Write to file to visually check result
  {
    const std::string filename =
      "grid_with_halo_" + Utilities::int_to_string(dim) + "d.vtk";
    std::ofstream f(filename.c_str());
    GridOut().write_vtk(tria, f);
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
