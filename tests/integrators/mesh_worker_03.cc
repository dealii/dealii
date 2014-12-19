// ---------------------------------------------------------------------
//
// Copyright (C) 2000 - 2013 by the deal.II authors
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


// Test consistency of assemblers MatrixSimple with and without local blocks

#include "../tests.h"
#include <deal.II/meshworker/assembler.h>
#include <deal.II/meshworker/loop.h>

#include <deal.II/base/logstream.h>
#include <deal.II/lac/compressed_sparsity_pattern.h>
#include <deal.II/lac/sparsity_pattern.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/mapping_q1.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_dgp.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/multigrid/mg_tools.h>

#include <fstream>
#include <iomanip>

using namespace dealii;

// Local integrators fill every matrix entry with 10*block_row + block_col
template <int dim>
class MatrixIntegrator : public Subscriptor
{
public:
  static void cell(MeshWorker::DoFInfo<dim> &dinfo,
                   MeshWorker::IntegrationInfo<dim> &info);
  static void face(MeshWorker::DoFInfo<dim> &dinfo1,
                   MeshWorker::DoFInfo<dim> &dinfo2,
                   MeshWorker::IntegrationInfo<dim> &info1,
                   MeshWorker::IntegrationInfo<dim> &info2);
  static void block_cell(MeshWorker::DoFInfo<dim> &dinfo,
                         MeshWorker::IntegrationInfo<dim> &info);
  static void block_face(MeshWorker::DoFInfo<dim> &dinfo1,
                         MeshWorker::DoFInfo<dim> &dinfo2,
                         MeshWorker::IntegrationInfo<dim> &info1,
                         MeshWorker::IntegrationInfo<dim> &info2);
};


template <int dim>
void MatrixIntegrator<dim>::cell(MeshWorker::DoFInfo<dim> &dinfo,
                                 MeshWorker::IntegrationInfo<dim> &info)
{
  const FiniteElement<dim> &fe = info.fe_values().get_fe();
  FullMatrix<double> &local_matrix = dinfo.matrix(0).matrix;

  for (unsigned int i=0; i<fe.dofs_per_cell; ++i)
    for (unsigned int j=0; j<fe.dofs_per_cell; ++j)
      {
        local_matrix(i,j) = 10 * fe.system_to_block_index(i).first
                            + fe.system_to_block_index(j).first;
      }
}


template <int dim>
void MatrixIntegrator<dim>::face(MeshWorker::DoFInfo<dim> &dinfo1,
                                 MeshWorker::DoFInfo<dim> &dinfo2,
                                 MeshWorker::IntegrationInfo<dim> &info1,
                                 MeshWorker::IntegrationInfo<dim> &info2)
{
  const FiniteElement<dim> &fe1 = info1.fe_values().get_fe();
  const FiniteElement<dim> &fe2 = info2.fe_values().get_fe();
  FullMatrix<double> &matrix_v1u2 = dinfo1.matrix(0, true).matrix;
  FullMatrix<double> &matrix_v2u1 = dinfo2.matrix(0, true).matrix;

  for (unsigned int i=0; i<fe1.dofs_per_cell; ++i)
    for (unsigned int j=0; j<fe1.dofs_per_cell; ++j)
      {
        matrix_v1u2(i,j) = 10 * fe1.system_to_block_index(i).first
                           + fe2.system_to_block_index(j).first;
        matrix_v2u1(i,j) = 10 * fe2.system_to_block_index(i).first
                           + fe1.system_to_block_index(j).first;
      }
}


template <int dim>
void MatrixIntegrator<dim>::block_cell(
  MeshWorker::DoFInfo<dim> &dinfo,
  MeshWorker::IntegrationInfo<dim> &)
{
  for (unsigned int m=0; m<dinfo.n_matrices(); ++m)
    {
      MatrixBlock<FullMatrix<double> > &loc = dinfo.matrix(m);

      for (unsigned int i=0; i<loc.matrix.m(); ++i)
        for (unsigned int j=0; j<loc.matrix.n(); ++j)
          {
            loc.matrix(i,j) = 10 * loc.row + loc.column;
          }
    }
}


template <int dim>
void MatrixIntegrator<dim>::block_face(
  MeshWorker::DoFInfo<dim> &dinfo1,
  MeshWorker::DoFInfo<dim> &dinfo2,
  MeshWorker::IntegrationInfo<dim> &,
  MeshWorker::IntegrationInfo<dim> &)
{
  for (unsigned int m=0; m<dinfo1.n_matrices(); ++m)
    {
      MatrixBlock<FullMatrix<double> > &loc = dinfo1.matrix(m, true);

      for (unsigned int i=0; i<loc.matrix.m(); ++i)
        for (unsigned int j=0; j<loc.matrix.n(); ++j)
          {
            loc.matrix(i,j) = 10 * loc.row + loc.column;
          }
    }
  for (unsigned int m=0; m<dinfo2.n_matrices(); ++m)
    {
      MatrixBlock<FullMatrix<double> > &loc = dinfo2.matrix(m, true);

      for (unsigned int i=0; i<loc.matrix.m(); ++i)
        for (unsigned int j=0; j<loc.matrix.n(); ++j)
          {
            loc.matrix(i,j) = 10 * loc.row + loc.column;
          }
    }
}


template <int dim>
void
assemble(const DoFHandler<dim> &dof_handler, SparseMatrix<double> &matrix)
{
  const FiniteElement<dim> &fe = dof_handler.get_fe();
  MappingQ1<dim> mapping;

  MeshWorker::IntegrationInfoBox<dim> info_box;
  const unsigned int n_gauss_points = dof_handler.get_fe().tensor_degree()+1;
  info_box.initialize_gauss_quadrature(n_gauss_points, n_gauss_points, n_gauss_points);
  info_box.initialize_update_flags();
  UpdateFlags update_flags = update_values | update_gradients;
  info_box.add_update_flags(update_flags, true, true, true, true);
  info_box.initialize(fe, mapping);

  MeshWorker::DoFInfo<dim> dof_info(dof_handler);

  MeshWorker::Assembler::MatrixSimple<SparseMatrix<double> > assembler;
  assembler.initialize(matrix);

  MeshWorker::loop<dim, dim, MeshWorker::DoFInfo<dim>, MeshWorker::IntegrationInfoBox<dim> >
  (dof_handler.begin_active(),
   dof_handler.end(),
   dof_info,
   info_box,
   &MatrixIntegrator<dim>::cell,
   0,
   &MatrixIntegrator<dim>::face,
   assembler);
}


template <int dim>
void
assemble(const MGDoFHandler<dim> &dof_handler,
         MGLevelObject<SparseMatrix<double> > matrix,
         MGLevelObject<SparseMatrix<double> > dg_up,
         MGLevelObject<SparseMatrix<double> > dg_down)
{
  const FiniteElement<dim> &fe = dof_handler.get_fe();
  MappingQ1<dim> mapping;

  MeshWorker::IntegrationInfoBox<dim> info_box;
  const unsigned int n_gauss_points = dof_handler.get_fe().tensor_degree()+1;
  info_box.initialize_gauss_quadrature(n_gauss_points, n_gauss_points, n_gauss_points);
  info_box.initialize_update_flags();
  UpdateFlags update_flags = update_values | update_gradients;
  info_box.add_update_flags(update_flags, true, true, true, true);
  info_box.initialize(fe, mapping);

  MeshWorker::DoFInfo<dim> dof_info(dof_handler);

  MeshWorker::Assembler::MGMatrixSimple<SparseMatrix<double> > assembler;
  assembler.initialize(matrix);

  MeshWorker::loop<MeshWorker::DoFInfo<dim>, MeshWorker::IntegrationInfoBox<dim> >
  (dof_handler.begin(),
   dof_handler.end(),
   dof_info,
   info_box,
   &MatrixIntegrator<dim>::cell,
   &MatrixIntegrator<dim>::bdry,
   &MatrixIntegrator<dim>::face,
   assembler);
}


template <int dim>
void
test_simple(MGDoFHandler<dim> &mgdofs)
{
  SparsityPattern pattern;
  SparseMatrix<double> matrix;
  Vector<double> v;

  const DoFHandler<dim> &dofs = mgdofs;
  const FiniteElement<dim> &fe = dofs.get_fe();
  pattern.reinit (dofs.n_dofs(), dofs.n_dofs(),
                  (GeometryInfo<dim>::faces_per_cell
                   *GeometryInfo<dim>::max_children_per_face+1)*fe.dofs_per_cell);
  DoFTools::make_flux_sparsity_pattern (dofs, pattern);
  pattern.compress();
  matrix.reinit (pattern);

  assemble(dofs, matrix);
  deallog << std::setprecision(3);
  matrix.print_formatted(deallog.get_file_stream(),0,false,4);

  MGLevelObject<SparsityPattern> mg_sparsity;
  MGLevelObject<SparsityPattern> mg_sparsity_dg_interface;
  MGLevelObject<SparseMatrix<double> > mg_matrix;
  MGLevelObject<SparseMatrix<double> > mg_matrix_dg_up;
  MGLevelObject<SparseMatrix<double> > mg_matrix_dg_down;

  const unsigned int n_levels = mgdofs.get_tria().n_levels();

  mg_sparsity.resize(0, n_levels-1);
  mg_sparsity_dg_interface.resize(0, n_levels-1);
  mg_matrix.resize(0, n_levels-1);
  mg_matrix_dg_up.resize(0, n_levels-1);
  mg_matrix_dg_down.resize(0, n_levels-1);

  for (unsigned int level=mg_sparsity.min_level();
       level<=mg_sparsity.max_level(); ++level)
    {
      CompressedSparsityPattern c_sparsity(mgdofs.n_dofs(level));
      CompressedSparsityPattern ci_sparsity;
      if (level>0)
        ci_sparsity.reinit(mgdofs.n_dofs(level-1), mgdofs.n_dofs(level));

      MGTools::make_flux_sparsity_pattern(mgdofs, c_sparsity, level);
      if (level>0)
        MGTools::make_flux_sparsity_pattern_edge(mgdofs, ci_sparsity, level);

      mg_sparsity[level].copy_from(c_sparsity);
      mg_matrix[level].reinit(mg_sparsity[level]);
      if (level>0)
        {
          mg_sparsity_dg_interface[level].copy_from(ci_sparsity);
          mg_matrix_dg_up[level].reinit(mg_sparsity_dg_interface[level]);
          mg_matrix_dg_down[level].reinit(mg_sparsity_dg_interface[level]);
        }
    }

}


template<int dim>
void
test(const FiniteElement<dim> &fe)
{
  deallog << fe.get_name() << std::endl;
  Triangulation<dim> tr;
  GridGenerator::hyper_L(tr);
  tr.begin()->set_refine_flag();
  tr.execute_coarsening_and_refinement();
//  tr.begin(2)->set_refine_flag();
//  tr.execute_coarsening_and_refinement();
//  tr.refine_global(1);
  deallog << "Triangulation levels";
  for (unsigned int l=0; l<tr.n_levels(); ++l)
    deallog << ' ' << l << ':' << tr.n_cells(l);
  deallog << std::endl;

  unsigned int cn = 0;
  for (typename Triangulation<dim>::cell_iterator cell = tr.begin();
       cell != tr.end(); ++cell, ++cn)
    cell->set_user_index(cn);

  MGDoFHandler<dim> dofs(tr);
  dofs.distribute_dofs(fe);
  deallog << "DoFHandler " << dofs.n_dofs() << " levels";
  for (unsigned int l=0; l<tr.n_levels(); ++l)
    deallog << ' ' << l << ':' << dofs.n_dofs(l);
  deallog << std::endl;

  test_simple(dofs);
}


int main ()
{
  const std::string logname = "output";
  std::ofstream logfile(logname.c_str());
  deallog.attach(logfile);
  deallog.depth_console (0);

  FE_DGP<2> dgp0(0);
  FE_DGP<2> dgp1(1);
  FE_Q<2> q1(1);

  std::vector<std_cxx11::shared_ptr<FiniteElement<2> > > fe2;
  fe2.push_back(std_cxx11::shared_ptr<FiniteElement<2> >(new FE_DGP<2>(0)));
  fe2.push_back(std_cxx11::shared_ptr<FiniteElement<2> >(new FESystem<2>(dgp0,3)));
  fe2.push_back(std_cxx11::shared_ptr<FiniteElement<2> >(new FESystem<2>(dgp1,2)));
  fe2.push_back(std_cxx11::shared_ptr<FiniteElement<2> >(new FESystem<2>(q1,2)));
  fe2.push_back(std_cxx11::shared_ptr<FiniteElement<2> >(new FESystem<2>(dgp0,1,q1,1)));
//  fe2.push_back(std_cxx11::shared_ptr<FiniteElement<2> >(new  FE_Q<2>(1)));

  for (unsigned int i=0; i<fe2.size(); ++i)
    test(*fe2[i]);
}
