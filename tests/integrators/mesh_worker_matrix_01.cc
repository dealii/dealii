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


// Test whether the various assember classes put the right data in the
// right place.

#include "../tests.h"
#include <deal.II/meshworker/assembler.h>
#include <deal.II/meshworker/loop.h>

#include <deal.II/base/std_cxx11/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/lac/sparsity_pattern.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/mapping_q1.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_dgp.h>
#include <deal.II/fe/fe_system.h>

#include <fstream>
#include <iomanip>

using namespace dealii;


// Define a class that fills all available entries in the info objects
// with recognizable numbers.
template <int dim>
class Local : public Subscriptor
{
public:
  typedef MeshWorker::IntegrationInfo<dim> CellInfo;

  void cell(MeshWorker::DoFInfo<dim> &dinfo, CellInfo &info) const;
  void bdry(MeshWorker::DoFInfo<dim> &dinfo, CellInfo &info) const;
  void face(MeshWorker::DoFInfo<dim> &dinfo1, MeshWorker::DoFInfo<dim> &dinfo2,
            CellInfo &info1, CellInfo &info2) const;

  bool cells;
  bool faces;
};

template <int dim>
void
Local<dim>::cell(MeshWorker::DoFInfo<dim> &info, CellInfo &) const
{
  if (!cells) return;
  for (unsigned int k=0; k<info.n_matrices(); ++k)
    {
      const unsigned int block_row = info.matrix(k).row;
      const unsigned int block_col = info.matrix(k).column;
      FullMatrix<double> &M1 = info.matrix(k).matrix;
      if (block_row == block_col)
        for (unsigned int i=0; i<M1.m(); ++i)
          for (unsigned int j=0; j<M1.n(); ++j)
            {
              M1(i,j) = 10.;
            }
    }
}


template <int dim>
void
Local<dim>::bdry(MeshWorker::DoFInfo<dim> &info, CellInfo &) const
{
  if (!faces) return;
  for (unsigned int k=0; k<info.n_matrices(); ++k)
    {
      const unsigned int block_row = info.matrix(k).row;
      const unsigned int block_col = info.matrix(k).column;
      FullMatrix<double> &M1 = info.matrix(k).matrix;
      if (block_row == block_col)
        for (unsigned int i=0; i<M1.m(); ++i)
          for (unsigned int j=0; j<M1.n(); ++j)
            {
              M1(i,j) = 1.;
            }
    }
}


template <int dim>
void
Local<dim>::face(MeshWorker::DoFInfo<dim> &info1, MeshWorker::DoFInfo<dim> &info2,
                 CellInfo &, CellInfo &) const
{
  if (!faces) return;
  for (unsigned int k=0; k<info1.n_matrices(); ++k)
    {
      const unsigned int block_row = info1.matrix(k).row;
      const unsigned int block_col = info1.matrix(k).column;
      FullMatrix<double> &M1 = info1.matrix(k).matrix;
      if (block_row == block_col)
        for (unsigned int i=0; i<M1.m(); ++i)
          for (unsigned int j=0; j<M1.n(); ++j)
            {
              info1.matrix(k,false).matrix(i,j) = 1.;
              info2.matrix(k,false).matrix(i,j) = 1.;
              info1.matrix(k,true).matrix(i,j) = -1.;
              info2.matrix(k,true).matrix(i,j) = -1.;
            }
    }
}


template <int dim>
void
test_simple(DoFHandler<dim> &dofs, bool faces)
{
  SparsityPattern pattern;
  SparseMatrix<double> matrix;

  const FiniteElement<dim> &fe = dofs.get_fe();
  pattern.reinit (dofs.n_dofs(), dofs.n_dofs(),
                  (GeometryInfo<dim>::faces_per_cell
                   *GeometryInfo<dim>::max_children_per_face+1)*fe.dofs_per_cell);
  DoFTools::make_flux_sparsity_pattern (dofs, pattern);
  pattern.compress();
  matrix.reinit (pattern);

  Local<dim> local;
  local.cells = true;
  local.faces = faces;

  MappingQ1<dim> mapping;

  MeshWorker::IntegrationInfoBox<dim> info_box;
  info_box.initialize_gauss_quadrature(1, 1, 1);
  info_box.initialize_update_flags();
  info_box.initialize(fe, mapping);

  MeshWorker::DoFInfo<dim> dof_info(dofs.block_info());

  MeshWorker::Assembler::MatrixSimple<SparseMatrix<double> > assembler;
  assembler.initialize(matrix);

  MeshWorker::loop<dim, dim, MeshWorker::DoFInfo<dim>, MeshWorker::IntegrationInfoBox<dim> >
  (dofs.begin_active(), dofs.end(),
   dof_info, info_box,
   std_cxx11::bind (&Local<dim>::cell, local, std_cxx11::_1, std_cxx11::_2),
   std_cxx11::bind (&Local<dim>::bdry, local, std_cxx11::_1, std_cxx11::_2),
   std_cxx11::bind (&Local<dim>::face, local, std_cxx11::_1, std_cxx11::_2, std_cxx11::_3, std_cxx11::_4),
   assembler, true);

  //matrix.print_formatted(deallog.get_file_stream(), 0, false, 4);
  matrix.print(deallog.get_file_stream(), false, false);
}


template<int dim>
void
test(const FiniteElement<dim> &fe)
{
  Triangulation<dim> tr;
  GridGenerator::hyper_cube(tr);
  tr.refine_global(2);
  // tr.begin(1)->set_refine_flag();
  // tr.execute_coarsening_and_refinement();
  // tr.begin(2)->set_refine_flag();
  // tr.execute_coarsening_and_refinement();
//  tr.refine_global(1);
  deallog << "Triangulation levels";
  for (unsigned int l=0; l<tr.n_levels(); ++l)
    deallog << ' ' << l << ':' << tr.n_cells(l);
  deallog << std::endl;

  unsigned int cn = 0;
  for (typename Triangulation<dim>::cell_iterator cell = tr.begin();
       cell != tr.end(); ++cell, ++cn)
    cell->set_user_index(cn);

  DoFHandler<dim> dofs(tr);
  dofs.distribute_dofs(fe);
  dofs.initialize_local_block_info();
  deallog << "DoFHandler " << dofs.n_dofs() << std::endl;

  test_simple(dofs, false);
  deallog << "now with jump terms" << std::endl;
  test_simple(dofs, true);
}


int main ()
{
  const std::string logname = "output";
  std::ofstream logfile(logname.c_str());
  deallog.attach(logfile);
  deallog.depth_console (0);

  FE_DGP<2> p0(0);
  FE_Q<2>   q1(1);
  FESystem<2,2> sys1(p0,1,q1,1);
  std::vector<FiniteElement<2>*> fe2;
  fe2.push_back(&p0);
  fe2.push_back(&q1);
  fe2.push_back(&sys1);

  for (unsigned int i=0; i<fe2.size(); ++i)
    test(*fe2[i]);
}
