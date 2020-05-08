// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2020 by the deal.II authors
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



// Test GridGenerator::plate_with_a_hole() with all padding non-equal

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/numerics/data_out.h>

#include "../tests.h"

template <int dim>
void
test()
{
  Point<dim> new_center;
  for (unsigned int d = 0; d < dim; ++d)
    new_center[d] = 10.;

  Triangulation<dim> triangulation;
  GridGenerator::plate_with_a_hole(
    triangulation, 0.4, 1., 1., 2., 0.8, 1.2, new_center, 2, 3, 1., 2, true);

  triangulation.refine_global(1);

  using Tuple = std::tuple<Point<dim>, types::boundary_id, types::manifold_id>;
  std::vector<Tuple> boundary_faces;

  for (const auto &cell : triangulation.active_cell_iterators())
    for (const unsigned int face_n : GeometryInfo<dim>::face_indices())
      if (cell->face(face_n)->at_boundary())
        boundary_faces.push_back(
          std::make_tuple(cell->face(face_n)->center(),
                          cell->face(face_n)->boundary_id(),
                          cell->face(face_n)->manifold_id()));

  // see dof_tools.cc internal::ComparisonHelper
  std::sort(begin(boundary_faces),
            end(boundary_faces),
            [](const Tuple &t1, const Tuple &t2) {
              const auto &b1 = std::get<1>(t1);
              const auto &b2 = std::get<1>(t2);
              if (b1 != b2)
                return b1 < b2;

              const auto &p1 = std::get<0>(t1);
              const auto &p2 = std::get<0>(t2);

              for (unsigned int d = 0; d < dim; ++d)
                if (std::abs(p1[d] - p2[d]) > 1e-8)
                  return p1[d] < p2[d];

              return std::get<2>(t1) < std::get<2>(t2);
            });

  for (const auto el : boundary_faces)
    deallog << "center: " << std::get<0>(el)
            << " boundary id: " << std::get<1>(el)
            << " manifold id: " << std::get<2>(el) << std::endl;

  deallog << std::endl << std::endl;
  Vector<float> manifold_id(triangulation.n_active_cells());
  unsigned int  index = 0;
  for (const auto &cell : triangulation.active_cell_iterators())
    {
      deallog << "center: " << cell->center()
              << " manifold id: " << cell->manifold_id() << std::endl;
      manifold_id[index] = cell->manifold_id();
      index++;
    }

  // write in vtk for visual inspection
  if (false)
    {
      DoFHandler<dim> dof_handler(triangulation);
      FE_Q<dim>       fe(1);
      dof_handler.distribute_dofs(fe);

      DataOut<dim> data_out;

      data_out.attach_dof_handler(dof_handler);
      data_out.add_data_vector(manifold_id, "manifold_id");
      data_out.build_patches();

      const std::string filename =
        "output_" + Utilities::int_to_string(dim) + ".vtu";
      std::ofstream output(filename.c_str());
      data_out.write_vtu(output);
    }
}


int
main()
{
  initlog();

  deallog.push("2d");
  test<2>();
  deallog.pop();

  deallog.push("3d");
  test<3>();
  deallog.pop();
}
