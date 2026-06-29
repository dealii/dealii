/* ------------------------------------------------------------------------
 *
 * SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
 * Copyright (C) 2023 - 2026 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * Part of the source code is dual licensed under Apache-2.0 WITH
 * LLVM-exception OR LGPL-2.1-or-later. Detailed license information
 * governing the source code and code contributions can be found in
 * LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
 *
 * ------------------------------------------------------------------------
 */

// Check that adjust_quad_dof_index_for_face_orientation gives consistent DoF
// indices. Test a mesh with two tetrahedra for all possible orientations.
// Similar to orientation_02 but checks that DoFs on the faces are correct.

#include <deal.II/base/function_lib.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_simplex_p.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/vector.h>

#include <deal.II/numerics/vector_tools_integrate_difference.h>
#include <deal.II/numerics/vector_tools_interpolate.h>

#include "../tests.h"

std::tuple<unsigned int, unsigned int>
create_triangulation(const std::vector<Point<3>>    &vertices_,
                     const std::vector<CellData<3>> &cell_data_,
                     const unsigned int              face_n,
                     const unsigned int              n_permutations,
                     Triangulation<3>               &tria)
{
  const ReferenceCell<3> ref_cell  = ReferenceCells::Tetrahedron;
  auto                   vertices  = vertices_;
  auto                   cell_data = cell_data_;

  Point<3> extra_vertex;
  for (unsigned int i = 0; i < 3; ++i)
    extra_vertex += ref_cell.vertex(ref_cell.face_to_cell_vertices(
      face_n, i, numbers::default_geometric_orientation));

  extra_vertex /= 3.0;
  extra_vertex += ref_cell.face_normal_vector(face_n);

  vertices.push_back(extra_vertex);

  cell_data.emplace_back();
  cell_data.back().vertices.resize(0);
  for (unsigned int i = 0; i < 3; ++i)
    cell_data.back().vertices.push_back(ref_cell.face_to_cell_vertices(
      face_n, i, numbers::default_geometric_orientation));
  cell_data.back().vertices.push_back(ref_cell.n_vertices());
  std::sort(cell_data.back().vertices.begin(), cell_data.back().vertices.end());

  unsigned int permutation_n = 0;
  do
    {
      tria.clear();
      tria.create_triangulation(vertices, cell_data, SubCellData());
      ++permutation_n;
    }
  while ((permutation_n < n_permutations) &&
         std::next_permutation(cell_data.back().vertices.begin(),
                               cell_data.back().vertices.end()));

  const auto cell = tria.begin();

  const auto face = cell->face(face_n);

  auto ncell = tria.begin();
  ncell++;

  unsigned int nf = 0;
  for (; nf < ref_cell.n_faces(); ++nf)
    if (ncell->face(nf) == face)
      break;

  return {nf, ncell->combined_face_orientation(nf)};
}


template <int dim>
void
test(const unsigned int degree)
{
  for (unsigned int f = 0; f < 4; ++f)
    {
      for (unsigned int r = 0; r < 24; ++r)
        {
          unsigned int orientation = r;
          unsigned int face_no     = f;

          Triangulation<dim> tria;

          Triangulation<dim> dummy;
          GridGenerator::reference_cell(dummy, ReferenceCells::Tetrahedron);

          auto vertices = dummy.get_vertices();

          std::vector<CellData<dim>> cells;

          {
            CellData<dim> cell;
            cell.vertices    = {0, 1, 2, 3};
            cell.material_id = 0;
            cells.push_back(cell);
          }

          std::tie(face_no, orientation) =
            create_triangulation(vertices, cells, face_no, r, tria);

          FE_SimplexP<dim> fe(degree);
          DoFHandler<dim>  dof_handler(tria);
          dof_handler.distribute_dofs(fe);

          const auto cell          = dof_handler.begin();
          auto       neighbor_cell = dof_handler.begin();
          neighbor_cell++;

          Assert(cell->neighbor(f) == neighbor_cell, ExcInternalError());
          Assert(neighbor_cell->neighbor(face_no) == cell, ExcInternalError());

          deallog << "Degree: " << degree << std::endl;
          deallog << "Face number, orientation: " << f << " "
                  << std::to_string(cell->combined_face_orientation(f))
                  << " neighbor face number, orientation: " << face_no << " "
                  << orientation << std::endl;


          std::vector<types::global_dof_index> dof_indices(
            fe.n_dofs_per_cell());
          std::vector<types::global_dof_index> neighbor_dof_indices(
            fe.n_dofs_per_cell());

          cell->get_dof_indices(dof_indices);
          neighbor_cell->get_dof_indices(neighbor_dof_indices);

          std::vector<types::global_dof_index> dof_indices_on_quad(
            fe.n_dofs_per_quad(f), numbers::invalid_dof_index);
          std::vector<types::global_dof_index> neighbor_dof_indices_on_quad(
            fe.n_dofs_per_quad(face_no), numbers::invalid_dof_index);

          for (unsigned int i = 0; i < fe.n_dofs_per_quad(f); ++i)
            dof_indices_on_quad[i] = dof_indices[4 + 6 * fe.n_dofs_per_line() +
                                                 f * fe.n_dofs_per_quad(f) + i];

          for (unsigned int i = 0; i < fe.n_dofs_per_quad(face_no); ++i)
            neighbor_dof_indices_on_quad[i] =
              neighbor_dof_indices[4 + 6 * fe.n_dofs_per_line() +
                                   face_no * fe.n_dofs_per_quad(face_no) + i];

          if (degree == 3)
            {
              // 1 DoF per face
              for (unsigned int i = 0; i < neighbor_dof_indices_on_quad.size();
                   ++i)
                {
                  if (neighbor_dof_indices_on_quad[i] == dof_indices_on_quad[i])
                    deallog << " ok ";
                  else
                    deallog << " error !!!!" << std::endl;
                }
            }
          else if (degree == 4)
            {
              for (unsigned int i = 0; i < neighbor_dof_indices_on_quad.size();
                   ++i)
                {
                  // 3 DoFs per face
                  // permute according to orientation
                  // same as applying the inverse permutation to
                  // neighbor_dof_indices_on_quad
                  if (dof_indices_on_quad[fe.reference_cell()
                                            .face_to_cell_vertices(
                                              0, i, orientation)] ==
                      neighbor_dof_indices_on_quad[i])
                    deallog << " ok ";
                  else
                    {
                      deallog << " error !!!!" << std::endl;
                    }
                }
            }
          deallog << std::endl;
        }
    }

  deallog << std::endl;
}

int
main()
{
  initlog();

  test<3>(3);
}
