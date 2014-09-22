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

// Fill local structures with the following format:
// 1. Residuals: CC.DDD
// 2. Matrices:  CC.DDDEEE

//   CC : 2 digits number of cell
//
//   DDD: degree of freedom of test function in cell, 1 digit block, 2
//        digits number in block
//   EEE: degree of freedom of trial function in cell, same format as
//        DDD

template <int dim>
void
Local<dim>::cell(MeshWorker::DoFInfo<dim> &info, CellInfo &) const
{
  if (!cells) return;
  const unsigned int cell = info.cell->user_index();

  // Fill local residuals
  for (unsigned int k=0; k<info.n_vectors(); ++k)
    for (unsigned int b=0; b<info.vector(k).n_blocks(); ++b)
      for (unsigned int i=0; i<info.vector(k).block(b).size(); ++i)
        {
          const double x = cell + 0.1 * b + 0.001 * i;
          info.vector(k).block(b)(i) = x;
        }

  for (unsigned int k=0; k<info.n_matrices(); ++k)
    {
      const unsigned int block_row = info.matrix(k).row;
      const unsigned int block_col = info.matrix(k).column;
      FullMatrix<double> &M1 = info.matrix(k).matrix;
      for (unsigned int i=0; i<M1.m(); ++i)
        for (unsigned int j=0; j<M1.n(); ++j)
          {
            double x = .1 * block_row + .001 * i;
            x = .1 * block_col + .001 * j + .001 * x;
            M1(i,j) = cell + x;
          }
    }
}


template <int dim>
void
Local<dim>::bdry(MeshWorker::DoFInfo<dim> & , CellInfo &) const
{}


template <int dim>
void
Local<dim>::face(MeshWorker::DoFInfo<dim> & , MeshWorker::DoFInfo<dim> &,
                 CellInfo &, CellInfo &) const
{}


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
  v.reinit (dofs.n_dofs());

  Local<dim> local;
  local.cells = true;
  local.faces = false;

  MappingQ1<dim> mapping;

  MeshWorker::IntegrationInfoBox<dim> info_box;
  info_box.initialize_gauss_quadrature(1, 1, 1);
  info_box.initialize_update_flags();
  info_box.initialize(fe, mapping);

  MeshWorker::DoFInfo<dim> dof_info(dofs);

  MeshWorker::Assembler::SystemSimple<SparseMatrix<double>, Vector<double> >
  assembler;
  assembler.initialize(matrix, v);

  MeshWorker::loop<dim, dim, MeshWorker::DoFInfo<dim>, MeshWorker::IntegrationInfoBox<dim> >
  (dofs.begin_active(), dofs.end(),
   dof_info, info_box,
   std_cxx11::bind (&Local<dim>::cell, local, std_cxx11::_1, std_cxx11::_2),
   std_cxx11::bind (&Local<dim>::bdry, local, std_cxx11::_1, std_cxx11::_2),
   std_cxx11::bind (&Local<dim>::face, local, std_cxx11::_1, std_cxx11::_2, std_cxx11::_3, std_cxx11::_4),
   assembler, true);

  for (unsigned int i=0; i<v.size(); ++i)
    deallog << ' ' << std::setprecision(3) << v(i);
  deallog << std::endl;

  deallog << std::setprecision(6);
  matrix.print(deallog, true);
}


template<int dim>
void
test(const FiniteElement<dim> &fe)
{
  Triangulation<dim> tr;
  GridGenerator::hyper_cube(tr);
  tr.refine_global(1);
  tr.begin(1)->set_refine_flag();
  tr.execute_coarsening_and_refinement();
  tr.begin(2)->set_refine_flag();
  tr.execute_coarsening_and_refinement();
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

  std::vector<std_cxx11::shared_ptr<FiniteElement<2> > > fe2;
  fe2.push_back(std_cxx11::shared_ptr<FiniteElement<2> >(new  FE_DGP<2>(1)));
  fe2.push_back(std_cxx11::shared_ptr<FiniteElement<2> >(new  FE_Q<2>(1)));

  for (unsigned int i=0; i<fe2.size(); ++i)
    test(*fe2[i]);
}
