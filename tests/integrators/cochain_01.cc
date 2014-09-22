// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2013 by the deal.II authors
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


// Test, whether differential operators produce a cochain complex on
// the standard Hilbert space sequence

#include "../tests.h"
#include "../test_grids.h"

#include <deal.II/base/logstream.h>

#include <deal.II/lac/matrix_block.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/block_sparsity_pattern.h>
#include <deal.II/lac/compressed_sparsity_pattern.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_nedelec.h>
#include <deal.II/fe/fe_raviart_thomas.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/meshworker/integration_info.h>
#include <deal.II/meshworker/assembler.h>
#include <deal.II/meshworker/loop.h>

#include <deal.II/integrators/l2.h>
#include <deal.II/integrators/divergence.h>
#include <deal.II/integrators/laplace.h>
#include <deal.II/integrators/maxwell.h>
#include <deal.II/integrators/elasticity.h>

using namespace LocalIntegrators;

const bool debugging = false;

template <int dim>
void cell_matrix(
  MeshWorker::DoFInfo<dim> &dinfo,
  typename MeshWorker::IntegrationInfo<dim,dim> &info)
{
  unsigned int dm=0; // Matrix index
  unsigned int de=0; // Element index

  L2::mass_matrix(dinfo.matrix(dm++,false).matrix, info.fe_values(de));
  Laplace::cell_matrix(dinfo.matrix(dm++,false).matrix, info.fe_values(de));
  Divergence::gradient_matrix(dinfo.matrix(dm++,false).matrix, info.fe_values(de), info.fe_values(de+1));


  ++de;
  L2::mass_matrix(dinfo.matrix(dm++,false).matrix, info.fe_values(de));
  Maxwell::curl_curl_matrix(dinfo.matrix(dm++,false).matrix, info.fe_values(de));
  Maxwell::curl_matrix(dinfo.matrix(dm++,false).matrix, info.fe_values(de), info.fe_values(de+1));

  if (dim>2)
    {
      ++de;
      L2::mass_matrix(dinfo.matrix(dm++,false).matrix, info.fe_values(de));
      Divergence::grad_div_matrix(dinfo.matrix(dm++,false).matrix, info.fe_values(de));
      Divergence::cell_matrix(dinfo.matrix(dm++,false).matrix, info.fe_values(de), info.fe_values(de+1));
    }

  ++de;
  L2::mass_matrix(dinfo.matrix(dm++,false).matrix, info.fe_values(de));
}



template <int dim>
void
test_cochain(const Triangulation<dim> &tr, const FiniteElement<dim> &fe)
{
  MappingQ1<dim> mapping;
  // Initialize DofHandler for a
  // block system with local blocks
  DoFHandler<dim> dof(tr);
  dof.distribute_dofs(fe);
  DoFRenumbering::component_wise(dof);
  dof.initialize_local_block_info();

  // Initialize hanging node constraints
  ConstraintMatrix constraints;
//  DoFTools::make_zero_boundary_constraints(dof, constraints);
  DoFTools::make_hanging_node_constraints (dof, constraints);
  constraints.close();
  // Setup sparsity pattern for the
  // whole system
  BlockSparsityPattern sparsity;
  sparsity.reinit(dof.block_info().global().size(),
                  dof.block_info().global().size());
  BlockCompressedSparsityPattern c_sparsity(dof.block_info().global(),
                                            dof.block_info().global());
  DoFTools::make_flux_sparsity_pattern(dof, c_sparsity, constraints, false);
  sparsity.copy_from(c_sparsity);

  // Setup matrices
  MatrixBlockVector<SparseMatrix<double> > matrices;

  unsigned int d=0;
  matrices.add(d,d,"mass-H1");
  matrices.add(d,d,"Laplacian");
  matrices.add(d+1,d,"grad");

  ++d;
  matrices.add(d,d,"mass-Hcurl");
  matrices.add(d,d,"curl-curl");
  matrices.add(d+1,d,"curl");

  if (dim>2)
    {
      ++d;
      matrices.add(d,d,"mass-Hdiv");
      matrices.add(d,d,"grad-div");
      matrices.add(d+1,d,"div");
    }

  ++d;
  matrices.add(d,d,"mass-L2");

  matrices.reinit(sparsity);

  if (debugging)
    for (unsigned int i=0; i<matrices.size(); ++i)
      deallog << "Block " << '(' << matrices.block(i).row
              << ',' << matrices.block(i).column << ") "
              << std::setw(3) << std::right << matrices.matrix(i).m() << std::left
              << 'x' << std::setw(3) << matrices.matrix(i).n()
              << ' ' << matrices.name(i) << std::endl;

  // Build matrices

  MeshWorker::IntegrationInfoBox<dim> info_box;
  UpdateFlags update_flags = update_values | update_gradients;
  info_box.add_update_flags_cell(update_flags);
  info_box.initialize(fe, mapping, &dof.block_info());
  if (debugging)
    deallog << "Infobox ready" << std::endl;

  MeshWorker::DoFInfo<dim> dof_info(dof.block_info());

  MeshWorker::Assembler::MatrixLocalBlocksToGlobalBlocks<SparseMatrix<double> > assembler;
  assembler.initialize(&dof.block_info(), matrices);
  assembler.initialize(constraints);

  MeshWorker::loop<dim, dim, MeshWorker::DoFInfo<dim>, MeshWorker::IntegrationInfoBox<dim> >(
    dof.begin_active(), dof.end(), dof_info, info_box,
    cell_matrix<dim>, 0, 0, assembler);

  for (unsigned int b=0; b<matrices.size(); ++b)
    if (b%3 == 0)
      for (unsigned int i=0; i<matrices.block(b).matrix.m(); ++i)
        if (matrices.block(b).matrix.diag_element(i) == 0.)
          matrices.block(b).matrix.diag_element(i) = 1.;
  if (debugging)
    deallog << "Matrices ready" << std::endl;
  // Set up vectors
  BlockVector<double> source(dof.block_info().global());
  BlockVector<double> result1(dof.block_info().global());
  BlockVector<double> result2(dof.block_info().global());
  BlockVector<double> aux(dof.block_info().global());
  for (unsigned int i=0; i<source.size(); ++i)
    source(i) = i%5;

  // now check, whether d*d =
  // D^TM^{-1}D
  SolverControl control(100, 1.e-13, false, false);
  SolverCG<Vector<double> > solver(control);

  for (unsigned int d=0; d<dim; ++d)
    {
      deallog << "Form " << d << std::endl;
      const unsigned int m=3*d;

      if (d>0)
        {
          matrices.matrix(m+2).vmult(result2.block(d+1), aux.block(d));
          deallog << "d^2            " << result2.block(d+1).l2_norm() << std::endl;
          matrices.matrix(m+1).vmult(result2.block(d), aux.block(d));
          deallog << "d*d^2          " << result2.block(d).l2_norm() << std::endl;
        }

      matrices.matrix(m+2).vmult(result1.block(d+1), source.block(d));

      PreconditionSSOR<SparseMatrix<double> > precondition;
      precondition.initialize(matrices.matrix(m+3), 1.2);

      solver.solve(matrices.matrix(m+3), aux.block(d+1), result1.block(d+1), precondition);

      matrices.matrix(m+2).Tvmult(result1.block(d), aux.block(d+1));
      matrices.matrix(m+1).vmult(result2.block(d), source.block(d));

      if (debugging)
        deallog << "u " << source.block(d).l2_norm() << ' '
                << matrices.name(m+2) << ' ' << result1.block(d+1).l2_norm()
                << ' ' << matrices.name(m+3) << ' ' << aux.block(d+1).l2_norm()
                << ' ' << matrices.name(m+2) << "^T " << result1.block(d).l2_norm()
                << ' ' << matrices.name(m+1) << ' ' << result2.block(d).l2_norm() << std::endl;
      result2.block(d) -= result1.block(d);
      deallog << "Difference d*d " << result2.block(d).l2_norm() << std::endl;
    }
}

void run2d (unsigned int degree)
{
  std::ostringstream prefix;
  prefix << "d2-p" << degree;
  deallog.push(prefix.str());
  deallog << "Setup" << std::endl;

  FE_Q<2> h1(degree+1);
  FE_Nedelec<2> hdiv(degree);
  FE_DGQ<2> l2(degree);

  FESystem<2> fe(h1,1,hdiv,1,l2,1);

  if (true)
    {
      Triangulation<2> tr;
      TestGrids::hypercube(tr, 1);
      test_cochain(tr, fe);
    }

  if (true)
    {
      Triangulation<2> tr;
      TestGrids::hypercube(tr, 2, true);
      test_cochain(tr, fe);
    }

  if (true)
    {
      Triangulation<2> tr;
      TestGrids::hypercube(tr, 3, true);
      test_cochain(tr, fe);
    }

  deallog.pop();
}


void run3d (unsigned int degree)
{
  std::ostringstream prefix;
  prefix << "d3-p" << degree;
  deallog.push(prefix.str());
  deallog << "Setup" << std::endl;

  FE_Q<3> h1(degree+1);
  if (debugging) deallog << "H1" << std::endl;
  FE_Nedelec<3> hcurl(degree);
  if (debugging) deallog << "Hcurl" << std::endl;
  FE_RaviartThomas<3> hdiv(degree);
  if (debugging) deallog << "Hdiv" << std::endl;
  FE_DGQ<3> l2(degree);
  if (debugging) deallog << "L2" << std::endl;

  FESystem<3> fe(h1,1,hcurl,1,hdiv,1,l2,1);

  if (true)
    {
      Triangulation<3> tr;
      TestGrids::hypercube(tr, 1);
      test_cochain(tr, fe);
    }

  if (true)
    {
      Triangulation<3> tr(Triangulation<3>::limit_level_difference_at_vertices);
      TestGrids::hypercube(tr, 2, true);
      test_cochain(tr, fe);
    }
  deallog.pop();
}


int main()
{
  const std::string logname = "output";
  std::ofstream logfile(logname.c_str());
  deallog.attach(logfile);
  deallog.log_execution_time(false);
  if (!debugging)
    {
      deallog.depth_console(0);
      deallog.threshold_double(1.e-10);
    }

  run2d(0);
  run2d(1);
  run2d(2);

  run3d(0);
  run3d(1);
}
