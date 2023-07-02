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

// Verify that the combination of QProjector, FEFaceValues, and various face
// orientations place face quadrature points in the same locations on abutting
// cells.

#include <deal.II/fe/fe_nothing.h>
#include <deal.II/fe/fe_simplex_p.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"

template <int dim>
void
test(const unsigned int orientation)
{
  Triangulation<dim> tria;

  if (false)
    {
      unsigned int ref = 2;
      GridGenerator::subdivided_hyper_cube_with_simplices(
        tria, Utilities::pow(2, ref), -1.0, +1.0);
    }
  else
    {
      const unsigned int face_no = 0;

      Triangulation<3> dummy;
      GridGenerator::reference_cell(dummy, ReferenceCells::Tetrahedron);

      auto vertices = dummy.get_vertices();

      std::vector<CellData<3>> cells;

      {
        CellData<3> cell;
        cell.vertices    = {0, 1, 2, 3};
        cell.material_id = 0;
        cells.push_back(cell);
      }

      if (false)
        {
          const auto &face = dummy.begin()->face(face_no);
          const auto  permuted =
            ReferenceCell(ReferenceCells::Triangle)
              .permute_according_orientation(
                std::array<unsigned int, 3>{{face->vertex_index(0),
                                             face->vertex_index(1),
                                             face->vertex_index(2)}},
                orientation);

          for (const auto o : permuted)
            std::cout << o << " ";
          std::cout << std::endl;

          auto direction =
            cross_product_3d(vertices[permuted[1]] - vertices[permuted[0]],
                             vertices[permuted[2]] - vertices[permuted[0]]);
          direction = direction / direction.norm();

          std::cout << direction << std::endl;

          vertices.emplace_back(0.0, 0.0, direction[2]);

          CellData<3> cell;
          cell.vertices = {permuted[0], permuted[1], permuted[2], 4u};

          cell.material_id = 1;
          cells.push_back(cell);
        }
      else
        {
          const auto &face = dummy.begin()->face(face_no);
          const auto  permuted =
            ReferenceCell(ReferenceCells::Triangle)
              .permute_according_orientation(
                std::array<unsigned int, 3>{{0, 1, 2}}, orientation);

          for (const auto o : permuted)
            std::cout << o << " ";
          std::cout << std::endl;

          auto direction =
            cross_product_3d(vertices[permuted[1]] - vertices[permuted[0]],
                             vertices[permuted[2]] - vertices[permuted[0]]);
          direction = direction / direction.norm();

          std::cout << direction << std::endl;

          vertices.emplace_back(0.0, 0.0, direction[2]);

          CellData<3> cell;
          cell.vertices.resize(4);

          cell.vertices[permuted[0]] = face->vertex_index(0);
          cell.vertices[permuted[1]] = face->vertex_index(1);
          cell.vertices[permuted[2]] = face->vertex_index(2);
          cell.vertices[3]           = 4;

          cell.material_id = 1;
          cells.push_back(cell);
        }

      tria.create_triangulation(vertices, cells, {});

      for (const auto &cell : tria.active_cell_iterators())
        {
          for (const auto l : cell->line_indices())
            std::cout << cell->line_orientation(l) << " ";
          std::cout << std::endl;
        }
      std::cout << std::endl;
    }

  const ReferenceCell ref_cell = ReferenceCells::get_simplex<dim>();
  FE_Nothing<dim>     fe(ref_cell);
  const Mapping<dim> &mapping =
    ref_cell.template get_default_linear_mapping<dim>();

  // Quadrature<dim - 1> face_quad =
  //  ref_cell.face_reference_cell(0).template get_gauss_type_quadrature<dim -
  //  1>(
  //    2u);

  Quadrature<dim - 1> face_quad(
    FE_SimplexP<dim>(1).get_unit_face_support_points());

  FEFaceValues<dim> cell_face_values(mapping,
                                     fe,
                                     face_quad,
                                     update_quadrature_points);
  FEFaceValues<dim> ncell_face_values(mapping,
                                      fe,
                                      face_quad,
                                      update_quadrature_points);

  for (const auto &cell : tria.active_cell_iterators())
    for (unsigned int f : cell->face_indices())
      {
        auto ncell = cell->neighbor(f);
        if (ncell == tria.end())
          continue;

        auto face = cell->face(f);
        cell_face_values.reinit(cell, face);
        unsigned int nf = 0;
        for (; nf < ref_cell.n_faces(); ++nf)
          if (ncell->face(nf) == face)
            break;
        ncell_face_values.reinit(ncell, nf);

        double max_distance = 0.0;
        for (unsigned int q = 0; q < face_quad.size(); ++q)
          max_distance =
            std::max(max_distance,
                     cell_face_values.get_quadrature_points()[q].distance(
                       ncell_face_values.get_quadrature_points()[q]));

        std::cout << orientation << " "
                  << int(ncell->combined_face_orientation(nf)) << std::endl;

        if (max_distance > 1e-12)
          {
            std::cout << "cell points =" << std::endl;
            for (unsigned int q = 0; q < face_quad.size(); ++q)
              {
                std::cout << " " << cell_face_values.get_quadrature_points()[q]
                          << std::endl;
              }
            std::cout << std::endl;
            std::cout << "face orientation = "
                      << int(cell->combined_face_orientation(f)) << std::endl;

            std::cout << "ncell points =" << std::endl;
            for (unsigned int q = 0; q < face_quad.size(); ++q)
              {
                std::cout << " " << ncell_face_values.get_quadrature_points()[q]
                          << std::endl;
              }
            std::cout << std::endl;
            std::cout << "nf orientation = "
                      << int(ncell->combined_face_orientation(nf)) << std::endl;

            AssertDimension(orientation,
                            int(ncell->combined_face_orientation(nf)));
            AssertThrow(false, ExcMessage("Should match"));
          }
      }
}

int
main()
{
  initlog();
  // test<2>();

  for (unsigned int o = 0; o < 6; ++o)
    test<3>(o);

  deallog << "OK!" << std::endl;
}
