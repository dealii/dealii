// ---------------------------------------------------------------------
//
// Copyright (C) 2012 - 2013 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------


/**
 * @file Test whether Assembler::MatrixSimple writes the local blocks
 * into the right global positions
 */

#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/block_indices.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/meshworker/local_results.h>
#include <deal.II/meshworker/dof_info.h>
#include <deal.II/meshworker/simple.h>

#include <deal.II/fe/fe_dgp.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_raviart_thomas.h>
#include <deal.II/fe/fe_system.h>

using namespace dealii;

template <typename number>
void fill_matrices(MeshWorker::LocalResults<number> &results, bool face)
{
  for (unsigned int k=0; k<results.n_matrices(); ++k)
    {
      FullMatrix<number> &M = results.matrix(k, false).matrix;
      double base = 1000*(results.matrix(k).row+1) + 100*(results.matrix(k).column+1);
      for (unsigned int i=0; i<M.m(); ++i)
        for (unsigned int j=0; j<M.n(); ++j)
          {
            M(i,j) = base + 10*i+j;
            if (face)
              results.matrix(k, true).matrix(i,j) = base + 10*i+j;
          }
    }
}


template <int dim>
void test(FiniteElement<dim> &fe)
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

  typename DoFHandler<dim>::active_cell_iterator cell = dof.begin_active();
  typename DoFHandler<dim>::face_iterator face = cell->face(1);
  typename DoFHandler<dim>::active_cell_iterator neighbor = cell->neighbor(1);

  CompressedSparsityPattern csp(dof.n_dofs(),dof.n_dofs());
  DoFTools::make_flux_sparsity_pattern(dof, csp);
  TrilinosWrappers::SparsityPattern sparsity;
  sparsity.copy_from(csp);
  TrilinosWrappers::SparseMatrix M(sparsity);

  MeshWorker::Assembler::MatrixSimple<TrilinosWrappers::SparseMatrix> ass;
  ass.initialize(M);
  ass.initialize_local_blocks(dof.block_info().local());
  MeshWorker::DoFInfo<dim> info(dof.block_info());
  ass.initialize_info(info, false);
  MeshWorker::DoFInfo<dim> infon(dof.block_info());
  ass.initialize_info(infon, true);

  deallog << "cell" << std::endl;
  info.reinit(cell);
  fill_matrices(info, false);
  ass.assemble(info);
  M.print(deallog.get_file_stream());
  M = 0.;

  deallog << "face" << std::endl;
  ass.initialize_info(info, true);
  info.reinit(cell, face, 1);
  infon.reinit(neighbor, neighbor->face(0), 0);
  fill_matrices(info, true);
  fill_matrices(infon, true);
  ass.assemble(info, infon);
  M.print(deallog.get_file_stream());
}

int main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization (argc, argv);
  initlog();

  FE_DGP<2> p0(0);
  FE_DGP<2> p1(1);
  FE_RaviartThomas<2> rt0(0);
  FE_Q<2> q2(2);

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
