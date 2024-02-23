// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2010 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// Test interaction with p4est with a 3d mesh. this is a redux of the
// _03 test where p4est errors out for a mesh that has multiple
// disconnected components
//
// because it happens to be simple to do that here, we also test the
// corresponding 2d case here

#include <deal.II/base/tensor.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"

#include "coarse_grid_common.h"



// create a mesh that consists of two disconnected hypercubes
template <int dim>
void
create_disconnected_mesh(Triangulation<dim> &tria)
{
  std::vector<Point<dim>>    vertices(2 * GeometryInfo<dim>::vertices_per_cell);
  std::vector<CellData<dim>> cells(2);

  // cell 1
  {
    Point<dim> p1 = Point<dim>();
    Point<dim> p2;
    for (unsigned int i = 0; i < dim; ++i)
      p2[i] = 1;

    switch (dim)
      {
        case 2:
          vertices[0] = vertices[1] = p1;
          vertices[2] = vertices[3] = p2;

          vertices[1][0] = p2[0];
          vertices[2][0] = p1[0];
          break;
        case 3:
          vertices[0] = vertices[1] = vertices[2] = vertices[3] = p1;
          vertices[4] = vertices[5] = vertices[6] = vertices[7] = p2;

          vertices[1][0] = p2[0];
          vertices[2][1] = p2[1];
          vertices[3][0] = p2[0];
          vertices[3][1] = p2[1];

          vertices[4][0] = p1[0];
          vertices[4][1] = p1[1];
          vertices[5][1] = p1[1];
          vertices[6][0] = p1[0];

          break;
        default:
          DEAL_II_NOT_IMPLEMENTED();
      }

    // Prepare cell data
    for (const unsigned int i : GeometryInfo<dim>::vertex_indices())
      cells[0].vertices[i] = i;
    cells[0].material_id = 0;
  }

  // cell 2. shifted 2 units in x-direction
  {
    Point<dim> p1 = Point<dim>();
    Point<dim> p2;
    for (unsigned int i = 0; i < dim; ++i)
      p2[i] = 1;

    p1[0] += 2;
    p2[0] += 2;

    switch (dim)
      {
        case 2:
          vertices[GeometryInfo<dim>::vertices_per_cell + 0] =
            vertices[GeometryInfo<dim>::vertices_per_cell + 1] = p1;
          vertices[GeometryInfo<dim>::vertices_per_cell + 2] =
            vertices[GeometryInfo<dim>::vertices_per_cell + 3] = p2;

          vertices[GeometryInfo<dim>::vertices_per_cell + 1][0] = p2[0];
          vertices[GeometryInfo<dim>::vertices_per_cell + 2][0] = p1[0];
          break;
        case 3:
          vertices[GeometryInfo<dim>::vertices_per_cell + 0] =
            vertices[GeometryInfo<dim>::vertices_per_cell + 1] =
              vertices[GeometryInfo<dim>::vertices_per_cell + 2] =
                vertices[GeometryInfo<dim>::vertices_per_cell + 3] = p1;
          vertices[GeometryInfo<dim>::vertices_per_cell + 4] =
            vertices[GeometryInfo<dim>::vertices_per_cell + 5] =
              vertices[GeometryInfo<dim>::vertices_per_cell + 6] =
                vertices[GeometryInfo<dim>::vertices_per_cell + 7] = p2;

          vertices[GeometryInfo<dim>::vertices_per_cell + 1][0] = p2[0];
          vertices[GeometryInfo<dim>::vertices_per_cell + 2][1] = p2[1];
          vertices[GeometryInfo<dim>::vertices_per_cell + 3][0] = p2[0];
          vertices[GeometryInfo<dim>::vertices_per_cell + 3][1] = p2[1];

          vertices[GeometryInfo<dim>::vertices_per_cell + 4][0] = p1[0];
          vertices[GeometryInfo<dim>::vertices_per_cell + 4][1] = p1[1];
          vertices[GeometryInfo<dim>::vertices_per_cell + 5][1] = p1[1];
          vertices[GeometryInfo<dim>::vertices_per_cell + 6][0] = p1[0];

          break;
        default:
          DEAL_II_NOT_IMPLEMENTED();
      }

    // Prepare cell data
    for (const unsigned int i : GeometryInfo<dim>::vertex_indices())
      cells[1].vertices[i] = GeometryInfo<dim>::vertices_per_cell + i;
    cells[1].material_id = 0;
  }

  tria.create_triangulation(vertices, cells, SubCellData());
}



template <int dim>
void
test(std::ostream & /*out*/)
{
  parallel::distributed::Triangulation<dim> tr(
    MPI_COMM_WORLD,
    Triangulation<dim>::none,
    parallel::distributed::Triangulation<dim>::communicate_vertices_to_p4est);

  create_disconnected_mesh(tr);

  write_vtk(tr, "1");
}


int
main(int argc, char *argv[])
{
  initlog();
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  deallog.push("2d");
  test<2>(deallog.get_file_stream());
  deallog.pop();

  deallog.push("3d");
  test<3>(deallog.get_file_stream());
  deallog.pop();
}
