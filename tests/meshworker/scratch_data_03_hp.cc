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

// Solve Laplacian using SIPG + mesh_loop + ScratchData
// - hp variant

#include <deal.II/base/function_parser.h>
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

#include <deal.II/hp/fe_collection.h>
#include <deal.II/hp/q_collection.h>

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
  Triangulation<dim, spacedim>    tria;
  hp::FECollection<dim, spacedim> fe(FE_DGQ<dim, spacedim>(1),
                                     FE_DGQ<dim, spacedim>(2));
  DoFHandler<dim, spacedim>       dh(tria);

  FunctionParser<spacedim> rhs_function("1");
  FunctionParser<spacedim> boundary_function("0");

  AffineConstraints<double> constraints;
  constraints.close();

  GridGenerator::hyper_cube(tria);
  tria.refine_global(5);
  tria.execute_coarsening_and_refinement();

  const auto bounding_box = GridTools::compute_bounding_box(tria);
  for (auto &cell : dh.active_cell_iterators())
    {
      bool use_alt_index = true;
      for (unsigned int d = 0; d < spacedim; ++d)
        if (cell->center()[d] < bounding_box.center()[d])
          {
            use_alt_index = false;
            break;
          }

      if (use_alt_index)
        cell->set_active_fe_index(1);
      else
        cell->set_active_fe_index(0);
    }

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

  hp::QCollection<dim>     cell_quad(QGauss<dim>(3), QGauss<dim>(4));
  hp::QCollection<dim - 1> face_quad(QGauss<dim - 1>(3), QGauss<dim - 1>(4));

  UpdateFlags cell_flags = update_values | update_gradients |
                           update_quadrature_points | update_JxW_values;
  UpdateFlags face_flags = update_values | update_gradients |
                           update_quadrature_points | update_normal_vectors |
                           update_JxW_values;

  // Stabilization for SIPG
  double gamma = 100;

  using ScratchData = MeshWorker::ScratchData<dim, spacedim>;
  using CopyData = MeshWorker::CopyData<1 + GeometryInfo<dim>::faces_per_cell,
                                        1,
                                        1 + GeometryInfo<dim>::faces_per_cell>;

  ScratchData scratch(fe, cell_quad, cell_flags, face_quad, face_flags);
  CopyData    copy(0);

  auto cell = dh.begin_active();
  auto endc = dh.end();

  using Iterator = decltype(cell);

  auto cell_worker =
    [&rhs_function](const Iterator &cell, ScratchData &s, CopyData &c) {
      const auto &fev = s.reinit(cell);
      const auto &JxW = s.get_JxW_values();
      const auto &p   = s.get_quadrature_points();

      c.reinit(0, s.n_dofs_per_cell());
      c.local_dof_indices[0] = s.get_local_dof_indices();

      for (unsigned int q = 0; q < p.size(); ++q)
        for (unsigned int i = 0; i < fev.dofs_per_cell; ++i)
          {
            for (unsigned int j = 0; j < fev.dofs_per_cell; ++j)
              {
                c.matrices[0](i, j) +=
                  fev.shape_grad(i, q) * fev.shape_grad(j, q) * JxW[q];
              }
            c.vectors[0](i) +=
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

    for (unsigned int q = 0; q < p.size(); ++q)
      for (unsigned int i = 0; i < fev.dofs_per_cell; ++i)
        {
          for (unsigned int j = 0; j < fev.dofs_per_cell; ++j)
            {
              c.matrices[0](i, j) +=
                (-fev.shape_grad(i, q) * n[q] * fev.shape_value(j, q) +
                 -fev.shape_grad(j, q) * n[q] * fev.shape_value(i, q) +
                 gamma / cell->face(f)->diameter() * fev.shape_value(i, q) *
                   fev.shape_value(j, q)) *
                JxW[q];
            }
          c.vectors[0](i) +=
            ((gamma / cell->face(f)->diameter() * fev.shape_value(i, q) -
              fev.shape_grad(i, q) * n[q]) *
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
    const auto &fev  = s.reinit(cell, ncell, f, sf);
    const auto &JxW  = s.get_JxW_values();
    const auto &nfev = s.reinit_neighbor(cell, ncell, nf, nsf);

    c.reinit(f + 1, s.n_dofs_per_cell(), s.n_neighbor_dofs_per_cell());
    c.local_dof_indices[f + 1] = s.get_neighbor_dof_indices();

    const auto &p  = s.get_quadrature_points();
    const auto &n  = s.get_normal_vectors();
    const auto &nn = s.get_neighbor_normal_vectors();

    const double gh = gamma / cell->face(f)->diameter();

    for (unsigned int q = 0; q < p.size(); ++q)
      {
        for (unsigned int i = 0; i < fev.dofs_per_cell; ++i)
          {
            for (unsigned int j = 0; j < fev.dofs_per_cell; ++j)
              c.matrices[0](i, j) +=
                (-.5 * fev.shape_grad(i, q) * n[q] * fev.shape_value(j, q) +
                 -.5 * fev.shape_value(i, q) * n[q] * fev.shape_grad(j, q) +
                 gh * fev.shape_value(i, q) * fev.shape_value(j, q)) *
                JxW[q];

            for (unsigned int j = 0; j < nfev.dofs_per_cell; ++j)
              c.matrices[f + 1](i, j) +=
                (-.5 * fev.shape_grad(i, q) * nn[q] * nfev.shape_value(j, q) +
                 -.5 * fev.shape_value(i, q) * n[q] * nfev.shape_grad(j, q) -
                 gh * fev.shape_value(i, q) * nfev.shape_value(j, q)) *
                JxW[q];
          }
      }
  };

  auto copier = [&constraints, &matrix, &rhs](const CopyData &c) {
    constraints.distribute_local_to_global(
      c.matrices[0], c.vectors[0], c.local_dof_indices[0], matrix, rhs);

    // We don't necessarily visit all faces, so we have to filter
    // some of them out here.
    for (const unsigned int f : GeometryInfo<dim>::face_indices())
      if (c.local_dof_indices[f + 1].size() > 0)
        constraints.distribute_local_to_global(c.matrices[f + 1],
                                               c.local_dof_indices[0],
                                               c.local_dof_indices[f + 1],
                                               matrix);
  };

  mesh_loop(cell,
            endc,
            cell_worker,
            copier,
            scratch,
            copy,
            assemble_own_cells | assemble_boundary_faces |
              assemble_own_interior_faces_both,
            boundary_worker,
            face_worker);

  SparseDirectUMFPACK inv;
  inv.initialize(matrix);

  inv.vmult(solution, rhs);
  constraints.distribute(solution);

  deallog << "Linfty norm of solution " << solution.linfty_norm() << std::endl;
}


int
main()
{
  initlog();
  test<2, 2>();
}
