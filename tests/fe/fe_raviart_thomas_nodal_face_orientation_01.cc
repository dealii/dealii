// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2026 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------

// Check that the normal component of FE_RaviartThomasNodal is continuous
// across a face shared by two cells for all relative orientations of the
// cells in 3d: The second cell is rotated by each of the 24 rotations of the
// reference cube, such that each of its six faces is the shared face with
// each of the face orientations. Each global basis function is evaluated on
// both sides of the shared face.

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_raviart_thomas.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/tria.h>

#include <algorithm>
#include <set>

#include "../tests.h"

#include "../test_grids.h"


void
test(const unsigned int degree)
{
  constexpr int                 dim = 3;
  FE_RaviartThomasNodal<dim>    fe(degree);
  const MappingQ<dim>           mapping(1);
  const QGauss<dim - 1>         quad(degree + 2);
  const auto                    rots     = TestGrids::hypercube_rotations();
  double                        max_jump = 0;
  std::set<std::pair<int, int>> face_and_orientation;

  for (const auto &rotation : rots)
    {
      // two unit cubes sharing the face x=1
      std::vector<Point<dim>> vertices;
      for (unsigned int k = 0; k < 2; ++k)
        for (unsigned int j = 0; j < 2; ++j)
          for (unsigned int i = 0; i < 3; ++i)
            vertices.emplace_back(i, j, k);
      std::vector<CellData<dim>> cells(2);
      for (unsigned int c = 0; c < 2; ++c)
        for (unsigned int v = 0; v < 8; ++v)
          cells[c].vertices[v] =
            (c + (v & 1)) + 3 * ((v >> 1) & 1) + 6 * (v >> 2);
      TestGrids::rotate_cell(cells[1], rotation);

      Triangulation<dim> tria;
      tria.create_triangulation(vertices, cells, SubCellData());
      DoFHandler<dim> dof_handler(tria);
      dof_handler.distribute_dofs(fe);

      const auto   cell     = dof_handler.begin_active();
      const auto   neighbor = std::next(cell);
      unsigned int face     = 0;
      while (cell->at_boundary(face))
        ++face;
      const unsigned int neighbor_face = cell->neighbor_face_no(face);
      face_and_orientation.emplace(
        neighbor_face, neighbor->combined_face_orientation(neighbor_face));

      FEFaceValues<dim> fe_values(mapping,
                                  fe,
                                  quad,
                                  update_values | update_quadrature_points |
                                    update_normal_vectors);
      FEFaceValues<dim> fe_values_neighbor(
        mapping, fe, quad, update_values | update_quadrature_points);
      fe_values.reinit(cell, face);
      fe_values_neighbor.reinit(neighbor, neighbor_face);
      std::vector<types::global_dof_index> dof_indices(fe.n_dofs_per_cell()),
        dof_indices_neighbor(fe.n_dofs_per_cell());
      cell->get_dof_indices(dof_indices);
      neighbor->get_dof_indices(dof_indices_neighbor);

      const FEValuesExtractors::Vector velocity(0);
      for (types::global_dof_index i = 0; i < dof_handler.n_dofs(); ++i)
        for (const unsigned int q : fe_values.quadrature_point_indices())
          {
            // find the matching quadrature point on the neighbor
            unsigned int q_neighbor = 0;
            for (const unsigned int p :
                 fe_values_neighbor.quadrature_point_indices())
              if (fe_values_neighbor.quadrature_point(p).distance(
                    fe_values.quadrature_point(q)) < 1e-12)
                q_neighbor = p;
            Tensor<1, dim> value, value_neighbor;
            for (unsigned int j = 0; j < fe.n_dofs_per_cell(); ++j)
              {
                if (dof_indices[j] == i)
                  value += fe_values[velocity].value(j, q);
                if (dof_indices_neighbor[j] == i)
                  value_neighbor +=
                    fe_values_neighbor[velocity].value(j, q_neighbor);
              }
            max_jump = std::max(max_jump,
                                std::abs((value - value_neighbor) *
                                         fe_values.normal_vector(q)));
          }
    }

  deallog << fe.get_name() << ": " << face_and_orientation.size()
          << " combinations of face and orientation, normal component "
          << (max_jump < 1e-12 ? "continuous" : "NOT continuous") << std::endl;
}



int
main()
{
  initlog();
  for (unsigned int degree = 1; degree < 3; ++degree)
    test(degree);
}
