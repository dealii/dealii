// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2022 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_tests_tensor_produce_matrix_h
#define dealii_tests_tensor_produce_matrix_h

#include <deal.II/fe/fe_nothing.h>
#include <deal.II/fe/fe_values.h>

namespace dealii
{
  namespace GridTools
  {
    // Compute harminic extent of all non-artificial cells.
    template <int dim>
    std::vector<std::array<double, dim>>
    compute_harmonic_cell_extent(const Mapping<dim>        &mapping,
                                 const Triangulation<dim>  &triangulation,
                                 const Quadrature<dim - 1> &quadrature)
    {
      std::vector<std::array<double, dim>> result(
        triangulation.n_active_cells());

      FE_Nothing<dim>   fe_nothing;
      FEFaceValues<dim> fe_face_values_0(mapping,
                                         fe_nothing,
                                         quadrature,
                                         update_quadrature_points);
      FEFaceValues<dim> fe_face_values_1(mapping,
                                         fe_nothing,
                                         quadrature,
                                         update_quadrature_points);

      for (const auto &cell : triangulation.active_cell_iterators())
        if (cell->is_artificial() == false)
          {
            for (unsigned int d = 0; d < dim; ++d)
              {
                fe_face_values_0.reinit(cell, 2 * d + 0);
                fe_face_values_1.reinit(cell, 2 * d + 1);

                double extent = 0.0;

                for (unsigned int q = 0; q < quadrature.size(); ++q)
                  extent += fe_face_values_0.quadrature_point(q).distance(
                              fe_face_values_1.quadrature_point(q)) *
                            quadrature.weight(q);

                result[cell->active_cell_index()][d] = extent;
              }
          }

      return result;
    }

    // Compute harmonic extent of each locally owned cell including of each
    // of its neighbors. If there is no neighbor, its extent is zero.
    template <int dim>
    std::vector<dealii::ndarray<double, dim, 3>>
    compute_harmonic_patch_extent(const Mapping<dim>        &mapping,
                                  const Triangulation<dim>  &triangulation,
                                  const Quadrature<dim - 1> &quadrature)
    {
      // 1) compute extent of each non-artificial cell
      const auto harmonic_cell_extents =
        GridTools::compute_harmonic_cell_extent(mapping,
                                                triangulation,
                                                quadrature);

      // 2) accumulate for each face the normal extent for the
      // neighboring cell(s); here we also consider periodicies
      std::vector<double> face_extent(triangulation.n_faces(), 0.0);

      for (const auto &cell : triangulation.active_cell_iterators())
        if (cell->is_artificial() == false)
          for (unsigned int d = 0; d < dim; ++d)
            {
              const auto extent =
                harmonic_cell_extents[cell->active_cell_index()][d];

              const auto add_extent_to_faces = [&](const unsigned int face_no) {
                face_extent[cell->face(face_no)->index()] += extent;

                if (cell->has_periodic_neighbor(face_no) &&
                    (cell->periodic_neighbor(face_no)->is_artificial() ==
                     false))
                  face_extent[cell->periodic_neighbor(face_no)
                                ->face(cell->periodic_neighbor_face_no(face_no))
                                ->index()] += extent;
              };

              add_extent_to_faces(2 * d + 0); // face 0
              add_extent_to_faces(2 * d + 1); // face 1
            }

      // 3) collect cell extent including those of the neighboring
      // cells, which corresponds to the difference of extent of the
      // current cell and the face extent
      std::vector<dealii::ndarray<double, dim, 3>> result(
        triangulation.n_active_cells());

      for (const auto &cell : triangulation.active_cell_iterators())
        if (cell->is_locally_owned())
          for (unsigned int d = 0; d < dim; ++d)
            {
              const auto cell_extent =
                harmonic_cell_extents[cell->active_cell_index()][d];

              const auto index = cell->active_cell_index();

              result[index][d][0] =
                face_extent[cell->face(2 * d + 0)->index()] - cell_extent;
              result[index][d][1] = cell_extent;
              result[index][d][2] =
                face_extent[cell->face(2 * d + 1)->index()] - cell_extent;
            }

      return result;
    }

  } // namespace GridTools
} // namespace dealii

#endif
