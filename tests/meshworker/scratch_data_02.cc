/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2018 - 2020 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * The deal.II library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE.md at
 * the top level directory of deal.II.
 *
 * ---------------------------------------------------------------------
 */

// Check mesh_loop + ScratchData for Laplacian (conforming)

#include <deal.II/base/function_parser.h>
#include <deal.II/base/patterns.h>
#include <deal.II/base/thread_management.h>

#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>

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
  FE_Q<dim, spacedim>          fe(1);
  DoFHandler<dim, spacedim>    dh(tria);

  FunctionParser<spacedim> rhs_function("1");
  FunctionParser<spacedim> boundary_function("0");

  AffineConstraints<double> constraints;

  GridGenerator::hyper_cube(tria);
  tria.refine_global(5);
  tria.execute_coarsening_and_refinement();
  dh.distribute_dofs(fe);


  DoFTools::make_hanging_node_constraints(dh, constraints);
  VectorTools::interpolate_boundary_values(dh,
                                           0,
                                           boundary_function,
                                           constraints);

  constraints.close();

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
  UpdateFlags face_flags = update_JxW_values;

  using ScratchData = MeshWorker::ScratchData<dim, spacedim>;
  using CopyData    = MeshWorker::CopyData<1, 1, 1>;

  ScratchData scratch(fe, quad, cell_flags, face_quad, face_flags);
  CopyData    copy(fe.dofs_per_cell);

  auto cell = dh.begin_active();
  auto endc = dh.end();

  typedef decltype(cell) Iterator;

  auto cell_worker =
    [&rhs_function](const Iterator &cell, ScratchData &s, CopyData &c) {
      const auto &fev = s.reinit(cell);
      const auto &JxW = s.get_JxW_values();
      const auto &p   = s.get_quadrature_points();

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

  auto copier = [&constraints, &matrix, &rhs](const CopyData &c) {
    constraints.distribute_local_to_global(
      c.matrices[0], c.vectors[0], c.local_dof_indices[0], matrix, rhs);
  };

  mesh_loop(cell, endc, cell_worker, copier, scratch, copy, assemble_own_cells);

  SparseDirectUMFPACK inv;
  inv.initialize(matrix);

  inv.vmult(solution, rhs);
  constraints.distribute(solution);

  deallog << "Linfty norm of solution: " << solution.linfty_norm() << std::endl;
}


int
main()
{
  initlog();

  test<2, 2>();
}
