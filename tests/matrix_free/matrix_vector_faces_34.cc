// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2021 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// similar to matrix_vector_faces_25 (matrix-free face evaluation,
// matrix-vector products as compared to the same implementation with
// MeshWorker, gather_evaluate and integrate_scatter, go through all different
// orientations), but for a system of equations and with local refinement

#include <deal.II/base/function.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_system.h>

#include "../tests.h"

#include "matrix_vector_faces_common.h"



void
generate_grid(Triangulation<3> &triangulation, int orientation)
{
  Point<3>              vertices_1[] = {Point<3>(-1., -1., -3.),
                                        Point<3>(+1., -1., -3.),
                                        Point<3>(-1., +1., -3.),
                                        Point<3>(+1., +1., -3.),
                                        Point<3>(-1., -1., -1.),
                                        Point<3>(+1., -1., -1.),
                                        Point<3>(-1., +1., -1.),
                                        Point<3>(+1., +1., -1.),
                                        Point<3>(-1., -1., +1.),
                                        Point<3>(+1., -1., +1.),
                                        Point<3>(-1., +1., +1.),
                                        Point<3>(+1., +1., +1.)};
  std::vector<Point<3>> vertices(&vertices_1[0], &vertices_1[12]);

  std::vector<CellData<3>> cells(2, CellData<3>());

  /* cell 0 */
  int cell_vertices_0[GeometryInfo<3>::vertices_per_cell] = {
    0, 1, 2, 3, 4, 5, 6, 7};

  /* cell 1 */
  int cell_vertices_1[8][GeometryInfo<3>::vertices_per_cell] = {
    {4, 5, 6, 7, 8, 9, 10, 11},
    {6, 4, 7, 5, 10, 8, 11, 9},
    {9, 8, 11, 10, 5, 4, 7, 6},
    {8, 10, 9, 11, 4, 6, 5, 7},
    {5, 7, 4, 6, 9, 11, 8, 10},
    {7, 6, 5, 4, 11, 10, 9, 8},
    {10, 11, 8, 9, 6, 7, 4, 5},
    {11, 9, 10, 8, 7, 5, 6, 4}};

  for (const unsigned int j : GeometryInfo<3>::vertex_indices())
    {
      cells[0].vertices[j] = cell_vertices_0[j];
      cells[1].vertices[j] = cell_vertices_1[orientation][j];
    }
  cells[0].material_id = 0;
  cells[1].material_id = 0;


  triangulation.create_triangulation(vertices, cells, SubCellData());

  const auto cell = ++(triangulation.begin());
  for (const auto face_no : cell->face_indices())
    if (!cell->face(face_no)->at_boundary())
      {
        deallog << "Orientation index within MatrixFree: "
                << (!cell->face_orientation(face_no) +
                    2 * cell->face_rotation(face_no) +
                    4 * cell->face_flip(face_no))
                << std::endl;
      }
}



void
run_test(const unsigned int fe_degree)
{
  FESystem<3> fe2(FE_DGQHermite<3>(fe_degree), 3);

  for (unsigned int orientation = 0; orientation < 8; ++orientation)
    {
      deallog << "Testing orientation case " << orientation << std::endl;
      for (unsigned int refine = 0; refine < 2; ++refine)
        {
          Triangulation<3> tria;
          generate_grid(tria, orientation);
          if (refine == 0)
            {
              tria.begin()->set_refine_flag();
              deallog << "Standard cell refined, oriented cell unrefined"
                      << std::endl;
            }
          else
            {
              (++tria.begin())->set_refine_flag();
              deallog << "Standard cell unrefined, oriented cell refined"
                      << std::endl;
            }
          tria.execute_coarsening_and_refinement();

          DoFHandler<3> dof(tria);
          dof.distribute_dofs(fe2);
          AffineConstraints<double> constraints;
          constraints.close();

          deallog << "Testing " << dof.get_fe().get_name();
          deallog << std::endl;

          MappingQ<3> mapping(dof.get_fe().degree + 1);

          Vector<double> in(dof.n_dofs()), out(dof.n_dofs());
          Vector<double> out_dist(out);

          // Set random seed for reproducibility
          Testing::srand(42);
          for (unsigned int i = 0; i < dof.n_dofs(); ++i)
            {
              if (constraints.is_constrained(i))
                continue;
              const double entry = Testing::rand() / (double)RAND_MAX;
              in(i)              = entry;
            }

          // assemble sparse matrix with MeshWorker for degrees less than 3
          // (the latter would be extraordinarily expensive already on 9 cells
          // as done here)
          if (fe_degree < 3)
            {
              SparsityPattern      sparsity;
              SparseMatrix<double> matrix;
              {
                DynamicSparsityPattern d_sparsity(dof.n_dofs());
                DoFTools::make_flux_sparsity_pattern(dof, d_sparsity);
                sparsity.copy_from(d_sparsity);
              }
              matrix.reinit(sparsity);
              MeshWorker::IntegrationInfoBox<3> info_box;
              UpdateFlags                       update_flags =
                update_values | update_gradients | update_jacobians;
              info_box.add_update_flags_all(update_flags);
              info_box.initialize_gauss_quadrature(dof.get_fe().degree + 1,
                                                   dof.get_fe().degree + 1,
                                                   dof.get_fe().degree + 1);
              info_box.initialize(dof.get_fe(), mapping);

              MeshWorker::DoFInfo<3> dof_info(dof);

              MeshWorker::Assembler::MatrixSimple<SparseMatrix<double>>
                assembler;
              assembler.initialize(matrix);

              MatrixIntegrator<3> integrator;
              MeshWorker::integration_loop<3, 3>(dof.begin_active(),
                                                 dof.end(),
                                                 dof_info,
                                                 info_box,
                                                 integrator,
                                                 assembler);

              matrix.vmult(out, in);
            }

          MatrixFree<3, double> mf_data;
          const QGauss<1>       quad(dof.get_fe().degree + 1);
          typename MatrixFree<3, double>::AdditionalData data;
          data.tasks_parallel_scheme =
            MatrixFree<3, double>::AdditionalData::none;
          data.tasks_block_size = 3;
          data.mapping_update_flags_inner_faces =
            (update_gradients | update_JxW_values);
          data.mapping_update_flags_boundary_faces =
            (update_gradients | update_JxW_values);

          mf_data.reinit(mapping, dof, constraints, quad, data);

          {
            MatrixFreeTest<3, -1, 0, double, Vector<double>, 3> mf(mf_data);
            mf.vmult(out_dist, in);

            if (fe_degree < 3)
              {
                out_dist -= out;
                const double diff_norm =
                  out_dist.linfty_norm() / out.linfty_norm();
                deallog << "Norm of difference basic:           " << diff_norm
                        << std::endl;
              }
            else
              {
                // Cannot compare in case we have no matrix reference. Then
                // use this result for the other variant instead
                out = out_dist;
              }
          }

          {
            MatrixFreeVariant<3, -1, 0, double, Vector<double>, 3> mf(mf_data,
                                                                      true);
            mf.vmult(out_dist, in);

            out_dist -= out;
            const double diff_norm = out_dist.linfty_norm() / out.linfty_norm();
            deallog << "Norm of difference gather_evaluate: " << diff_norm
                    << std::endl;
          }
        }
    }
}



template <int dim, int fe_degree>
void
test()
{
  if (dim == 2)
    return;
  else
    run_test(fe_degree + 1);
}
