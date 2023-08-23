// ---------------------------------------------------------------------
//
// Copyright (C) 2023 by the deal.II authors
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

// Test a mesh with two tetrahedra for all possible orientations. Similar to
// orientation_02 but also checks that quadrature points on faces (computed via
// FEFaceValues) are correct.

#include <deal.II/fe/fe_nothing.h>
#include <deal.II/fe/fe_simplex_p.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/reference_cell.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"

void
test_face(const std::vector<Point<3>>    &vertices_,
          const std::vector<CellData<3>> &cell_data_,
          const unsigned int              face_n)
{
  const ReferenceCell ref_cell  = ReferenceCells::Tetrahedron;
  auto                vertices  = vertices_;
  auto                cell_data = cell_data_;

  Point<3> extra_vertex;
  for (unsigned int i = 0; i < 3; ++i)
    extra_vertex += ref_cell.template vertex<3>(ref_cell.face_to_cell_vertices(
      face_n, i, ReferenceCell::default_combined_face_orientation()));

  extra_vertex /= 3.0;
  extra_vertex += ref_cell.template unit_normal_vectors<3>(face_n);

  vertices.push_back(extra_vertex);

  cell_data.emplace_back();
  cell_data.back().vertices.resize(0);
  for (unsigned int i = 0; i < 3; ++i)
    cell_data.back().vertices.push_back(ref_cell.face_to_cell_vertices(
      face_n, i, ref_cell.default_combined_face_orientation()));
  cell_data.back().vertices.push_back(ref_cell.n_vertices());
  std::sort(cell_data.back().vertices.begin(), cell_data.back().vertices.end());

  unsigned int permutation_n = 0;
  do
    {
      Triangulation<3> tria;
      tria.create_triangulation(vertices, cell_data, SubCellData());
      deallog << "  p2 = " << permutation_n;

      FE_Nothing<3>     fe(ref_cell);
      const Mapping<3> &mapping =
        ref_cell.template get_default_linear_mapping<3>();
      Quadrature<2> face_quad(FE_SimplexP<3>(1).get_unit_face_support_points());

      FEFaceValues<3> cell_face_values(mapping,
                                       fe,
                                       face_quad,
                                       update_quadrature_points);
      FEFaceValues<3> neighbor_cell_face_values(mapping,
                                                fe,
                                                face_quad,
                                                update_quadrature_points);

      for (const auto &cell : tria.active_cell_iterators())
        for (unsigned int face_no : cell->face_indices())
          {
            auto neighbor_cell = cell->neighbor(face_no);
            if (neighbor_cell == tria.end())
              continue;

            auto face = cell->face(face_no);
            cell_face_values.reinit(cell, face);
            unsigned int neighbor_face_no = 0;
            for (; neighbor_face_no < ref_cell.n_faces(); ++neighbor_face_no)
              if (neighbor_cell->face(neighbor_face_no) == face)
                break;
            neighbor_cell_face_values.reinit(neighbor_cell, neighbor_face_no);

            deallog << " : " << int(cell->combined_face_orientation(face_no))
                    << ", "
                    << int(neighbor_cell->combined_face_orientation(
                         neighbor_face_no));

            double max_distance = 0.0;
            for (unsigned int q = 0; q < face_quad.size(); ++q)
              max_distance = std::max(
                max_distance,
                cell_face_values.get_quadrature_points()[q].distance(
                  neighbor_cell_face_values.get_quadrature_points()[q]));

            if (max_distance > 1e-12)
              {
                deallog << "!!!!! Found a wrong point permutation !!!!!"
                        << std::endl;

                deallog << "cell points =" << std::endl;
                for (unsigned int q = 0; q < face_quad.size(); ++q)
                  {
                    deallog << " "
                            << cell_face_values.get_quadrature_points()[q]
                            << std::endl;
                  }
                deallog << std::endl;

                deallog << "neighbor_cell points =" << std::endl;
                for (unsigned int q = 0; q < face_quad.size(); ++q)
                  {
                    deallog
                      << " "
                      << neighbor_cell_face_values.get_quadrature_points()[q]
                      << std::endl;
                  }
                deallog << std::endl;

                AssertThrow(false, ExcMessage("Should match"));
              }
          }
      deallog << std::endl;
      ++permutation_n;
    }
  while (std::next_permutation(cell_data.back().vertices.begin(),
                               cell_data.back().vertices.end()));
}

void
test()
{
  std::vector<Point<3>> vertices;
  vertices.emplace_back(0.0, 0.0, 0.0);
  vertices.emplace_back(1.0, 0.0, 0.0);
  vertices.emplace_back(0.0, 1.0, 0.0);
  vertices.emplace_back(0.0, 0.0, 1.0);

  std::vector<CellData<3>> cells;
  cells.emplace_back();
  cells.back().vertices = {0, 1, 2, 3};

  for (unsigned int face_n = 0; face_n < 4; ++face_n)
    {
      std::sort(cells.back().vertices.begin(), cells.back().vertices.end());
      unsigned int permutation_n = 0;
      deallog << "face_n = " << face_n << std::endl;
      do
        {
          deallog << "p1 = " << permutation_n << std::endl;
          test_face(vertices, cells, face_n);
          ++permutation_n;
        }
      while (std::next_permutation(cells.back().vertices.begin(),
                                   cells.back().vertices.end()));
    }
}

int
main()
{
  initlog();

  test();
}
