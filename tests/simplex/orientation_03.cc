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
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"

template <int dim>
void
test()
{
  unsigned int       ref = 2;
  Triangulation<dim> tria;
  GridGenerator::subdivided_hyper_cube_with_simplices(tria,
                                                      Utilities::pow(2, ref),
                                                      -1.0,
                                                      +1.0);

  const ReferenceCell ref_cell = ReferenceCells::get_simplex<dim>();
  FE_Nothing<dim>     fe(ref_cell);
  const Mapping<dim> &mapping =
    ref_cell.template get_default_linear_mapping<dim>();

  Quadrature<dim - 1> face_quad =
    ref_cell.face_reference_cell(0).template get_gauss_type_quadrature<dim - 1>(
      2u);

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

        if (max_distance > 1e-12)
          {
            std::cout << "cell points =";
            for (unsigned int q = 0; q < face_quad.size(); ++q)
              {
                std::cout << " " << cell_face_values.get_quadrature_points()[q];
              }
            std::cout << std::endl;
            std::cout << "face orientation = "
                      << int(cell->combined_face_orientation(f)) << std::endl;

            std::cout << "ncell points =";
            for (unsigned int q = 0; q < face_quad.size(); ++q)
              {
                std::cout << " "
                          << ncell_face_values.get_quadrature_points()[q];
              }
            std::cout << std::endl;
            std::cout << "nf orientation = "
                      << int(ncell->combined_face_orientation(nf)) << std::endl;
            AssertThrow(false, ExcMessage("Should match"));
          }
      }
}

int
main()
{
  initlog();
  test<2>();
  test<3>();

  deallog << "OK!" << std::endl;
}
