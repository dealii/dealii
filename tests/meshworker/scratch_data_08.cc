/* ------------------------------------------------------------------------
 *
 * SPDX-License-Identifier: LGPL-2.1-or-later
 * Copyright (C) 2019 - 2023 by the deal.II authors
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

// Solve Laplacian using SIPG + mesh_loop + ScratchData + FEInterfaceData

#include <deal.II/base/function_lib.h>
#include <deal.II/base/patterns.h>
#include <deal.II/base/thread_management.h>

#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_dgq.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/grid_tools_cache.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/sparse_direct.h>

#include <deal.II/meshworker/copy_data.h>
#include <deal.II/meshworker/mesh_loop.h>
#include <deal.II/meshworker/scratch_data.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>

#include <cmath>
#include <fstream>
#include <iostream>
#include <unordered_map>

#include "../tests.h"

using namespace MeshWorker;

template <int dim, int spacedim>
void
test()
{
  Triangulation<dim, spacedim> tria;
  FE_DGQ<dim, spacedim>        fe(1);
  DoFHandler<dim, spacedim>    dh(tria);

  Functions::ConstantFunction<spacedim> rhs_function(1);
  Functions::ConstantFunction<spacedim> boundary_function(0);

  AffineConstraints<double> constraints;
  constraints.close();

  GridGenerator::hyper_cube(tria);
  tria.refine_global(3);
  tria.execute_coarsening_and_refinement();
  dh.distribute_dofs(fe);


  SparsityPattern sparsity;

  {
    DynamicSparsityPattern dsp(dh.n_dofs(), dh.n_dofs());
    DoFTools::make_flux_sparsity_pattern(dh, dsp);
    sparsity.copy_from(dsp);
  }

  SparseMatrix<double> matrix;
  matrix.reinit(sparsity);

  Vector<double> solution(dh.n_dofs());
  Vector<double> rhs(dh.n_dofs());

  QGauss<dim>     quad(3);
  QGauss<dim - 1> face_quad(3);

  UpdateFlags cell_flags = update_values | update_gradients |
                           update_quadrature_points | update_JxW_values;
  UpdateFlags face_flags = update_values | update_gradients |
                           update_quadrature_points | update_normal_vectors |
                           update_JxW_values;

  // Stabilization for SIPG
  double gamma = 100;

  using ScratchData = MeshWorker::ScratchData<dim, spacedim>;

  auto cell = dh.begin_active();
  auto endc = dh.end();

  using Iterator = decltype(cell);

  struct CopyDataFace
  {
    FullMatrix<double>                   cell_matrix;
    std::vector<types::global_dof_index> joint_dof_indices;
  };

  struct CopyData
  {
    FullMatrix<double>                   cell_matrix;
    Vector<double>                       cell_rhs;
    std::vector<types::global_dof_index> local_dof_indices;
    std::vector<CopyDataFace>            face_data;

    void
    reinit(const Iterator &cell, unsigned int dofs_per_cell)
    {
      cell_matrix.reinit(dofs_per_cell, dofs_per_cell);
      cell_rhs.reinit(dofs_per_cell);

      local_dof_indices.resize(dofs_per_cell);
      cell->get_dof_indices(local_dof_indices);
    }
  };

  ScratchData scratch(fe, quad, cell_flags, face_quad, face_flags);
  CopyData    copy;

  auto cell_worker =
    [&rhs_function](const Iterator &cell, ScratchData &s, CopyData &c) {
      const auto &fev = s.reinit(cell);
      const auto &JxW = s.get_JxW_values();
      const auto &p   = s.get_quadrature_points();

      c.reinit(cell, s.get_local_dof_indices().size());

      for (unsigned int q = 0; q < p.size(); ++q)
        for (unsigned int i = 0; i < fev.dofs_per_cell; ++i)
          {
            for (unsigned int j = 0; j < fev.dofs_per_cell; ++j)
              {
                c.cell_matrix(i, j) +=
                  fev.shape_grad(i, q) * fev.shape_grad(j, q) * JxW[q];
              }
            c.cell_rhs(i) +=
              fev.shape_value(i, q) * rhs_function.value(p[q]) * JxW[q];
          }
    };

  auto boundary_worker = [gamma, &boundary_function](const Iterator     &cell,
                                                     const unsigned int &f,
                                                     ScratchData        &s,
                                                     CopyData           &c) {
    const auto &fev = s.reinit(cell, f);
    const auto &JxW = s.get_JxW_values();
    const auto &p   = s.get_quadrature_points();
    const auto &n   = s.get_normal_vectors();

    const double gh = gamma / cell->diameter();

    for (unsigned int q = 0; q < p.size(); ++q)
      for (unsigned int i = 0; i < fev.dofs_per_cell; ++i)
        {
          for (unsigned int j = 0; j < fev.dofs_per_cell; ++j)
            {
              c.cell_matrix(i, j) +=
                (-fev.shape_grad(i, q) * n[q] * fev.shape_value(j, q) +
                 -fev.shape_grad(j, q) * n[q] * fev.shape_value(i, q) +
                 gh * fev.shape_value(i, q) * fev.shape_value(j, q)) *
                JxW[q];
            }
          c.cell_rhs(i) +=
            ((gh * fev.shape_value(i, q) - fev.shape_grad(i, q) * n[q]) *
             boundary_function.value(p[q])) *
            JxW[q];
        }
  };

  auto face_worker = [gamma](const Iterator     &cell,
                             const unsigned int &f,
                             const unsigned int &sf,
                             const Iterator     &ncell,
                             const unsigned int &nf,
                             const unsigned int &nsf,
                             ScratchData        &s,
                             CopyData           &c) {
    const auto &fev = s.reinit(cell, f, sf, ncell, nf, nsf);
    const auto &JxW = s.get_JxW_values();

    const auto &p = s.get_quadrature_points();
    const auto &n = s.get_normal_vectors();


    c.face_data.emplace_back();
    auto &copy_data_face             = c.face_data.back();
    auto &face_matrix                = copy_data_face.cell_matrix;
    copy_data_face.joint_dof_indices = fev.get_interface_dof_indices();
    const auto n_dofs                = fev.n_current_interface_dofs();
    face_matrix.reinit(n_dofs, n_dofs);

    const double gh = gamma / cell->diameter();

    for (unsigned int q = 0; q < p.size(); ++q)
      for (unsigned int i = 0; i < n_dofs; ++i)
        for (unsigned int j = 0; j < n_dofs; ++j)
          {
            face_matrix(i, j) += (-fev.jump_in_shape_gradients(i, q) * n[q] *
                                    fev.average_of_shape_values(j, q) -
                                  fev.average_of_shape_values(i, q) *
                                    fev.jump_in_shape_gradients(j, q) * n[q] +
                                  gh * fev.jump_in_shape_values(i, q) *
                                    fev.jump_in_shape_values(j, q)) *
                                 JxW[q];
          }
  };

  auto copier = [&constraints, &matrix, &rhs](const CopyData &c) {
    constraints.distribute_local_to_global(
      c.cell_matrix, c.cell_rhs, c.local_dof_indices, matrix, rhs);

    for (auto &cdf : c.face_data)
      constraints.distribute_local_to_global(cdf.cell_matrix,
                                             cdf.joint_dof_indices,
                                             matrix);
  };



  mesh_loop(cell,
            endc,
            cell_worker,
            copier,
            scratch,
            copy,
            assemble_own_cells | assemble_boundary_faces |
              assemble_own_interior_faces_once,
            boundary_worker,
            face_worker);

  SparseDirectUMFPACK inv;
  inv.initialize(matrix);

  inv.vmult(solution, rhs);

  deallog << "Linfty norm of solution " << solution.linfty_norm() << std::endl;
}


int
main()
{
  initlog();
  test<1, 1>();
  test<1, 2>();
  test<1, 3>();
  test<2, 2>();
  test<2, 3>();
  test<3, 3>();
}
