// ---------------------------------------------------------------------
//
// Copyright (C) 2002 - 2020 by the deal.II authors
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


// read a file in the MSH format used by the GMSH program.
// and test that  files in GMSH-2 format and in GMSH-4 format produce the same
// result.

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <string>

#include "../tests.h"

template <int dim>
struct PointComparator
{
  bool
  operator()(const Point<dim> &lhs, const Point<dim> &rhs) const
  {
    double downstream_size = 0;
    double weight          = 1.;
    for (unsigned int d = 0; d < dim; ++d)
      {
        downstream_size += (rhs[d] - lhs[d]) * weight;
        weight *= 1e-5;
      }
    if (downstream_size < 0)
      return false;
    else if (downstream_size > 0)
      return true;
    else
      {
        for (unsigned int d = 0; d < dim; ++d)
          {
            if (lhs[d] == rhs[d])
              continue;
            return lhs[d] < rhs[d];
          }
        return false;
      }
  }
};


template <int dim>
void
gmsh_grid(const char *name_v2, const char *name_v4)
{
  Triangulation<dim> tria_v2;
  {
    GridIn<dim> grid_in;
    grid_in.attach_triangulation(tria_v2);
    std::ifstream input_file(name_v2);
    grid_in.read_msh(input_file);
  }

  Triangulation<dim> tria_v4;
  {
    GridIn<dim> grid_in;
    grid_in.attach_triangulation(tria_v4);
    std::ifstream input_file(name_v4);
    grid_in.read_msh(input_file);
  }

  AssertThrow(tria_v2.n_active_cells() == tria_v4.n_active_cells(),
              ExcInternalError());
  deallog << "  " << tria_v2.n_active_cells() << " active cells" << std::endl;

  // The cells and faces in the two files are in reversed order
  // but the vertex index is the same.
  auto cell_v2 = tria_v2.begin_active();
  auto cell_v4 =
    std::next(tria_v4.begin_active(), tria_v4.n_active_cells() - 1);
  const auto end_v2 = tria_v2.end();
  for (; cell_v2 != end_v2; ++cell_v2, --cell_v4)
    {
      AssertThrow(cell_v2->material_id() == cell_v4->material_id(),
                  ExcInternalError());
      std::set<unsigned int> vertices_v2;
      std::set<unsigned int> vertices_v4;
      for (const unsigned int i : GeometryInfo<dim>::vertex_indices())
        {
          AssertThrow((tria_v2.get_vertices()[cell_v2->vertex_index(i)] -
                       tria_v4.get_vertices()[cell_v2->vertex_index(i)])
                          .norm() < 1.e-10,
                      ExcInternalError());
          AssertThrow((tria_v2.get_vertices()[cell_v4->vertex_index(i)] -
                       tria_v4.get_vertices()[cell_v4->vertex_index(i)])
                          .norm() < 1.e-10,
                      ExcInternalError());
          vertices_v2.insert(cell_v2->vertex_index(i));
          vertices_v4.insert(cell_v4->vertex_index(i));
        }
      AssertThrow(vertices_v2 == vertices_v4, ExcInternalError());
      std::map<Point<dim>, types::boundary_id, PointComparator<dim>> faces_v2;
      std::map<Point<dim>, types::boundary_id, PointComparator<dim>> faces_v4;
      for (const unsigned int i : GeometryInfo<dim>::face_indices())
        {
          faces_v2[cell_v2->face(i)->center()] =
            cell_v2->face(i)->boundary_id();
          faces_v4[cell_v4->face(i)->center()] =
            cell_v4->face(i)->boundary_id();
        }
      AssertThrow(faces_v2 == faces_v4, ExcInternalError());
      std::map<Point<dim>, types::boundary_id, PointComparator<dim>> lines_v2;
      std::map<Point<dim>, types::boundary_id, PointComparator<dim>> lines_v4;
      for (unsigned int i = 0; i < GeometryInfo<dim>::lines_per_cell; ++i)
        {
          lines_v2[cell_v2->line(i)->center()] =
            cell_v2->line(i)->boundary_id();
          lines_v4[cell_v4->line(i)->center()] =
            cell_v4->line(i)->boundary_id();
        }
      AssertThrow(lines_v2 == lines_v4, ExcInternalError());
    }
  deallog << "  OK" << std::endl;
}

void
filename_resolution()
{
  deallog << "grid_in_msh_version_2/hole81" << std::endl;
  gmsh_grid<2>(SOURCE_DIR "/grid_in_msh_version_2/hole81.msh",
               SOURCE_DIR "/grid_in_msh_version_4/hole81.msh");
  deallog << "grid_in_msh_version_2/hole8170" << std::endl;
  gmsh_grid<2>(SOURCE_DIR "/grid_in_msh_version_2/hole8170.msh",
               SOURCE_DIR "/grid_in_msh_version_4/hole8170.msh");
}


int
main()
{
  initlog();
  deallog.get_file_stream() << std::setprecision(5);

  filename_resolution();
}
