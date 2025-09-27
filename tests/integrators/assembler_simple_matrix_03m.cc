// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2012 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


/**
 * @file Test whether Assembler::MatrixSimple writes the local blocks
 * into the right global positions of multiple matrices
 */

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_dgp.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_raviart_thomas.h>
#include <deal.II/fe/fe_system.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/block_indices.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>

#include <deal.II/meshworker/dof_info.h>
#include <deal.II/meshworker/local_results.h>
#include <deal.II/meshworker/simple.h>

#include "../tests.h"


template <typename number>
void
fill_matrices(MeshWorker::LocalResults<number> &results, bool face)
{
  for (unsigned int k = 0; k < results.n_matrices(); ++k)
    {
      FullMatrix<number> &M    = results.matrix(k, false).matrix;
      double              base = 1000 * (results.matrix(k).row + 1) +
                    100 * (results.matrix(k).column + 1);
      for (unsigned int i = 0; i < M.m(); ++i)
        for (unsigned int j = 0; j < M.n(); ++j)
          {
            number entry = base + 10 * i + j;
            //      if (k >= results.n_matrices()/2) entry *= -1;
            M(i, j) = entry;
            if (face)
              results.matrix(k, true).matrix(i, j) = entry;
          }
    }
}


template <int dim>
void
test(FiniteElement<dim> &fe)
{
  deallog << fe.get_name() << std::endl;

  Triangulation<dim> tr;
  GridGenerator::hyper_cube(tr);
  tr.refine_global(1);

  DoFHandler<dim> dof(tr);
  dof.distribute_dofs(fe);
  dof.initialize_local_block_info();
  DoFRenumbering::component_wise(dof);

  deallog << "DoFs " << dof.n_dofs() << std::endl;

  typename DoFHandler<dim>::active_cell_iterator cell     = dof.begin_active();
  typename DoFHandler<dim>::face_iterator        face     = cell->face(1);
  typename DoFHandler<dim>::active_cell_iterator neighbor = cell->neighbor(1);

  DynamicSparsityPattern csp(dof.n_dofs(), dof.n_dofs());
  DoFTools::make_flux_sparsity_pattern(dof, csp);
  SparsityPattern sparsity;
  sparsity.copy_from(csp);
  std::vector<SparseMatrix<double>> M(2);
  M[0].reinit(sparsity);
  M[1].reinit(sparsity);

  MeshWorker::Assembler::MatrixSimple<SparseMatrix<double>> ass;
  ass.initialize(M);
  MeshWorker::DoFInfo<dim> info(dof.block_info());
  ass.initialize_info(info, false);
  MeshWorker::DoFInfo<dim> infon(dof.block_info());
  ass.initialize_info(infon, true);

  deallog << "cell" << std::endl;
  info.reinit(cell);
  fill_matrices(info, false);
  ass.assemble(info);
  M[1].print_formatted(deallog.get_file_stream(), 0, false, 6);
  M[0] = 0.;
  M[1] = 0.;

  deallog << "face" << std::endl;
  ass.initialize_info(info, true);
  info.reinit(cell, face, 1);
  infon.reinit(neighbor, neighbor->face(0), 0);
  fill_matrices(info, true);
  fill_matrices(infon, true);
  ass.assemble(info, infon);
  M[1].print_formatted(deallog.get_file_stream(), 0, false, 6);
}

int
main()
{
  initlog();

  FE_DGP<2>           p0(0);
  FE_DGP<2>           p1(1);
  FE_RaviartThomas<2> rt0(0);
  FE_Q<2>             q2(2);

  FESystem<2> sys1(p0, 2, p1, 1);
  FESystem<2> sys2(p0, 2, rt0, 1);
  FESystem<2> sys3(rt0, 1, p0, 2);
  FESystem<2> sys4(p1, 2, q2, 2);
  FESystem<2> sys5(q2, 2, p1, 2);

  test(sys1);
  test(sys2);
  test(sys3);
  test(sys4);
  test(sys5);
}
