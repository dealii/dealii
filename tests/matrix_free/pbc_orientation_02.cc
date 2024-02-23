// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2020 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// Like pbc_orientation_02 but with cell-centric loops.

#include <deal.II/base/function.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_dgq.h>

#include <deal.II/grid/grid_tools.h>

#include <deal.II/matrix_free/fe_evaluation.h>
#include <deal.II/matrix_free/matrix_free.h>
#include <deal.II/matrix_free/operators.h>
#include <deal.II/matrix_free/tools.h>

#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"


template <int dim>
class ExactSolution : public Function<dim>
{
public:
  ExactSolution()
  {}

  virtual double
  value(const Point<dim> &p, const unsigned int /*component*/ = 0) const
  {
    // independent of p[2] since we apply PBC in z-direction
    return p[0] + p[1];
  }
};



void
generate_grid(Triangulation<3> &triangulation, int orientation)
{
  Point<3>              vertices_1[] = {Point<3>(-0., -0., -0.),
                                        Point<3>(+1., -0., -0.),
                                        Point<3>(-0., +1., -0.),
                                        Point<3>(+1., +1., -0.),
                                        Point<3>(-0., -0., +0.5),
                                        Point<3>(+1., -0., +0.5),
                                        Point<3>(-0., +1., +0.5),
                                        Point<3>(+1., +1., +0.5),
                                        Point<3>(-0., -0., +1.),
                                        Point<3>(+1., -0., +1.),
                                        Point<3>(-0., +1., +1.),
                                        Point<3>(+1., +1., +1.)};
  std::vector<Point<3>> vertices(&vertices_1[0], &vertices_1[12]);

  std::vector<CellData<3>> cells(2, CellData<3>());

  /* cell 0 */
  int cell_vertices_0[GeometryInfo<3>::vertices_per_cell] = {
    0, 1, 2, 3, 4, 5, 6, 7};

  /* cell 1 */
  int cell_vertices_1[8][GeometryInfo<3>::vertices_per_cell] = {
    {4, 5, 6, 7, 8, 9, 10, 11},
    {5, 7, 4, 6, 9, 11, 8, 10},
    {7, 6, 5, 4, 11, 10, 9, 8},
    {6, 4, 7, 5, 10, 8, 11, 9},
    {9, 8, 11, 10, 5, 4, 7, 6},
    {8, 10, 9, 11, 4, 6, 5, 7},
    {10, 11, 8, 9, 6, 7, 4, 5},
    {11, 9, 10, 8, 7, 5, 6, 4}};

  for (const unsigned int j : GeometryInfo<3>::vertex_indices())
    {
      cells[orientation < 8 ? 0 : 1].vertices[j] = cell_vertices_0[j];
      cells[orientation < 8 ? 1 : 0].vertices[j] =
        cell_vertices_1[orientation % 8][j];
    }
  cells[0].material_id = 0;
  cells[1].material_id = 0;


  triangulation.create_triangulation(vertices, cells, SubCellData());
}



template <int dim, int fe_degree>
void
test()
{
  if (dim == 2)
    return;

  Triangulation<3> tria;
  for (unsigned int orientation = 0; orientation < 16; ++orientation)
    {
      deallog << "Testing orientation case " << orientation << std::endl;
      tria.clear();
      generate_grid(tria, orientation);

      FE_DGQ<3>     fe(fe_degree);
      DoFHandler<3> dof_handler(tria);
      dof_handler.distribute_dofs(fe);

      for (auto &cell : tria.cell_iterators())
        for (auto &face : cell->face_iterators())
          {
            if (std::abs(face->center()[2] - 0.0) < 10e-6)
              face->set_boundary_id(4);
            if (std::abs(face->center()[2] - 1.0) < 10e-6)
              face->set_boundary_id(5);
          }

      std::vector<dealii::GridTools::PeriodicFacePair<
        typename dealii::Triangulation<dim>::cell_iterator>>
        periodic_faces;

      dealii::GridTools::collect_periodic_faces(tria, 4, 5, 2, periodic_faces);

      tria.add_periodicity(periodic_faces);

      QGauss<1> quadrature(fe_degree + 1);
      typename MatrixFree<dim, double>::AdditionalData additional_data;
      additional_data.mapping_update_flags                = update_values;
      additional_data.mapping_update_flags_inner_faces    = update_values;
      additional_data.mapping_update_flags_boundary_faces = update_values;
      additional_data.mapping_update_flags_faces_by_cells = update_values;
      additional_data.tasks_parallel_scheme =
        MatrixFree<dim, double>::AdditionalData::none;

      MatrixFreeTools::categorize_by_boundary_ids(
        dof_handler.get_triangulation(), additional_data);

      MatrixFree<dim, double> data;

      AffineConstraints<double> dummy;
      dummy.close();
      data.reinit(
        MappingQ1<dim>{}, dof_handler, dummy, quadrature, additional_data);

      using VectorType = LinearAlgebra::distributed::Vector<double>;

      VectorType vec;
      data.initialize_dof_vector(vec);

      VectorTools::interpolate(dof_handler, ExactSolution<dim>(), vec);



      data.template loop_cell_centric<VectorType, VectorType>(
        [](const auto &data, auto &, const auto &src, const auto cell_range) {
          FEFaceEvaluation<dim, fe_degree, fe_degree + 1, 1, double> eval_minus(
            data, true);
          FEFaceEvaluation<dim, fe_degree, fe_degree + 1, 1, double> eval_plus(
            data, false);

          for (unsigned int cell = cell_range.first; cell < cell_range.second;
               ++cell)
            for (unsigned int face = 0;
                 face < GeometryInfo<dim>::faces_per_cell;
                 ++face)
              {
                const auto boundary_ids =
                  data.get_faces_by_cells_boundary_id(cell, face);
                const auto boundary_id = boundary_ids[0];

                if (boundary_id != numbers::internal_face_boundary_id)
                  continue;

                eval_minus.reinit(cell, face);
                eval_minus.gather_evaluate(src, EvaluationFlags::values);
                eval_plus.reinit(cell, face);
                eval_plus.gather_evaluate(src, EvaluationFlags::values);

                for (unsigned int q = 0; q < eval_minus.n_q_points; ++q)
                  {
                    const auto u_minus = eval_minus.get_value(q);
                    const auto u_plus  = eval_plus.get_value(q);

                    for (unsigned int v = 0;
                         v < VectorizedArray<double>::size();
                         ++v)
                      {
                        Assert(std::abs(u_minus[v] - u_plus[v]) < 1e-10,
                               ExcMessage("Entries do not match!"));
                      }
                  }
              }
        },
        vec,
        vec);
    }
}

int
main()
{
  initlog();
  test<3, 3>();
}
