// ---------------------------------------------------------------------
//
// Copyright (C) 2000 - 2015 by the deal.II authors
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


// Test whether Assember::CellsAndFaces adds up correctly.

#include "../tests.h"
#include "empty_info.h"
#include <deal.II/meshworker/assembler.h>
#include <deal.II/meshworker/loop.h>

#include <deal.II/base/std_cxx11/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_dgp.h>
#include <deal.II/lac/block_vector.h>

#include <fstream>
#include <iomanip>

using namespace dealii;


// Define a class that fills all available entries in the info objects
// with recognizable numbers.

// Fills the following functionals:
// 0: number of cells
// 1: number of interior faces
// 2: number of boundary faces

const unsigned int n_functionals = 3;

template <int dim>
class Local : public Subscriptor
{
public:
  typedef EmptyInfo CellInfo;

  void cell(MeshWorker::DoFInfo<dim> &dinfo, CellInfo &info) const;
  void bdry(MeshWorker::DoFInfo<dim> &dinfo, CellInfo &info) const;
  void face(MeshWorker::DoFInfo<dim> &dinfo1, MeshWorker::DoFInfo<dim> &dinfo2,
            CellInfo &info1, CellInfo &info2) const;
};


template <int dim>
void
Local<dim>::cell(MeshWorker::DoFInfo<dim> &info, CellInfo &) const
{
  info.value(0) = 1.;
}


template <int dim>
void
Local<dim>::bdry(MeshWorker::DoFInfo<dim>  &info, CellInfo &) const
{
  info.value(2) = 1.;
}


template <int dim>
void
Local<dim>::face(MeshWorker::DoFInfo<dim>  &info1, MeshWorker::DoFInfo<dim> &info2,
                 CellInfo &, CellInfo &) const
{
  info1.value(1) = 1./2.;
  info2.value(1) = 1./2.;
}


template <int dim>
void
test_mesh(DoFHandler<dim> &mgdofs)
{
  const DoFHandler<dim> &dofs = mgdofs;

  BlockVector<double> cells(n_functionals);
  BlockVector<double> faces(n_functionals);
  for (unsigned int i=0; i<n_functionals; ++i)
    {
      cells.block(i).reinit(dofs.get_tria().n_cells());
      faces.block(i).reinit(dofs.get_tria().n_faces());
    }
  cells.collect_sizes();
  faces.collect_sizes();

  NamedData<BlockVector<double>* > out_data;
  BlockVector<double> *p = &cells;
  out_data.add(p, "cells");
  p = &faces;
  out_data.add(p, "faces");

  Local<dim> local;
  EmptyInfoBox info_box;
  MeshWorker::DoFInfo<dim> dof_info(dofs);

  MeshWorker::Assembler::CellsAndFaces<double> assembler;
  assembler.initialize(out_data, true);

  MeshWorker::LoopControl lctrl;
  lctrl.cells_first = true;
  lctrl.own_faces = MeshWorker::LoopControl::one;

  MeshWorker::loop<dim, dim, MeshWorker::DoFInfo<dim>, EmptyInfoBox>
  (dofs.begin_active(), dofs.end(),
   dof_info, info_box,
   std_cxx11::bind (&Local<dim>::cell, local, std_cxx11::_1, std_cxx11::_2),
   std_cxx11::bind (&Local<dim>::bdry, local, std_cxx11::_1, std_cxx11::_2),
   std_cxx11::bind (&Local<dim>::face, local, std_cxx11::_1, std_cxx11::_2, std_cxx11::_3, std_cxx11::_4),
   assembler, lctrl);

  deallog << "  Results cells";
  for (unsigned int i=0; i<n_functionals; ++i)
    deallog << '\t' << cells.block(i).l1_norm();
  deallog << std::endl;

  deallog << "  Results faces";
  for (unsigned int i=0; i<n_functionals; ++i)
    deallog << '\t' << faces.block(i).l1_norm();
  deallog << std::endl;

  cells = 0.;
  faces = 0.;

  MeshWorker::DoFInfo<dim> mg_dof_info(mgdofs);
  MeshWorker::loop<dim, dim, MeshWorker::DoFInfo<dim>, EmptyInfoBox>
  (mgdofs.begin_mg(), mgdofs.end_mg(),
   mg_dof_info, info_box,
   std_cxx11::bind (&Local<dim>::cell, local, std_cxx11::_1, std_cxx11::_2),
   std_cxx11::bind (&Local<dim>::bdry, local, std_cxx11::_1, std_cxx11::_2),
   std_cxx11::bind (&Local<dim>::face, local, std_cxx11::_1, std_cxx11::_2, std_cxx11::_3, std_cxx11::_4),
   assembler, lctrl);

  deallog << "MGResults cells";
  for (unsigned int i=0; i<n_functionals; ++i)
    deallog << '\t' << cells.block(i).l1_norm();
  deallog << std::endl;

  deallog << "MGResults faces";
  for (unsigned int i=0; i<n_functionals; ++i)
    deallog << '\t' << faces.block(i).l1_norm();
  deallog << std::endl;
}


template<int dim>
void
test(const FiniteElement<dim> &fe)
{
  Triangulation<dim> tr;
  DoFHandler<dim> dofs(tr);
  GridGenerator::hyper_cube(tr);
  tr.refine_global(1);
  deallog.push("1");
  dofs.distribute_dofs(fe);
  dofs.distribute_mg_dofs(fe);
  test_mesh(dofs);
  deallog.pop();
  tr.begin(1)->set_refine_flag();
  tr.execute_coarsening_and_refinement();
  deallog.push("2");
  dofs.distribute_dofs(fe);
  dofs.distribute_mg_dofs(fe);
  test_mesh(dofs);
  deallog.pop();
  tr.begin(2)->set_refine_flag();
  tr.execute_coarsening_and_refinement();
  deallog.push("3");
  dofs.distribute_dofs(fe);
  dofs.distribute_mg_dofs(fe);
  test_mesh(dofs);
  deallog.pop();
}


int main ()
{
  const std::string logname = "output";
  std::ofstream logfile(logname.c_str());
  deallog.attach(logfile);
  deallog.depth_console (0);

  FE_DGP<2> el2(0);
  FE_DGP<3> el3(0);

  deallog.push("2D");
  test(el2);
  deallog.pop();
  deallog.push("3D");
  test(el3);
  deallog.pop();
}
