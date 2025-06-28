// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2011 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Test whether Assembler::CellsAndFaces adds up correctly.

#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_dgp.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/lac/block_vector.h>

#include <deal.II/meshworker/assembler.h>
#include <deal.II/meshworker/loop.h>

#include <functional>

#include "../tests.h"

#include "empty_info.h"



// Define a class that fills all available entries in the info objects
// with recognizable numbers.

// Fills the following functionals:
// 0: number of cells
// 1: number of interior faces
// 2: number of boundary faces

const unsigned int n_functionals = 3;

template <int dim>
class Local : public EnableObserverPointer
{
public:
  using CellInfo = EmptyInfo;

  void
  cell(MeshWorker::DoFInfo<dim> &dinfo, CellInfo &info) const;
  void
  bdry(MeshWorker::DoFInfo<dim> &dinfo, CellInfo &info) const;
  void
  face(MeshWorker::DoFInfo<dim> &dinfo1,
       MeshWorker::DoFInfo<dim> &dinfo2,
       CellInfo                 &info1,
       CellInfo                 &info2) const;
};


template <int dim>
void
Local<dim>::cell(MeshWorker::DoFInfo<dim> &info, CellInfo &) const
{
  info.value(0) = 1.;
}


template <int dim>
void
Local<dim>::bdry(MeshWorker::DoFInfo<dim> &info, CellInfo &) const
{
  info.value(2) = 1.;
}


template <int dim>
void
Local<dim>::face(MeshWorker::DoFInfo<dim> &info1,
                 MeshWorker::DoFInfo<dim> &info2,
                 CellInfo &,
                 CellInfo &) const
{
  info1.value(1) = 1. / 2.;
  info2.value(1) = 1. / 2.;
}


template <int dim>
void
test_mesh(DoFHandler<dim> &mgdofs)
{
  const DoFHandler<dim> &dofs = mgdofs;

  BlockVector<double> cells(n_functionals);
  BlockVector<double> faces(n_functionals);
  for (unsigned int i = 0; i < n_functionals; ++i)
    {
      cells.block(i).reinit(dofs.get_triangulation().n_cells());
      faces.block(i).reinit(dofs.get_triangulation().n_faces());
    }
  cells.collect_sizes();
  faces.collect_sizes();

  AnyData out_data;
  out_data.add<BlockVector<double> *>(&cells, "cells");
  out_data.add<BlockVector<double> *>(&faces, "faces");

  Local<dim>               local;
  EmptyInfoBox             info_box;
  MeshWorker::DoFInfo<dim> dof_info(dofs);

  MeshWorker::Assembler::CellsAndFaces<double> assembler;
  assembler.initialize(out_data, true);

  MeshWorker::LoopControl lctrl;
  lctrl.cells_first = true;
  lctrl.own_faces   = MeshWorker::LoopControl::one;

  MeshWorker::loop<dim, dim, MeshWorker::DoFInfo<dim>, EmptyInfoBox>(
    dofs.begin_active(),
    dofs.end(),
    dof_info,
    info_box,
    std::bind(
      &Local<dim>::cell, local, std::placeholders::_1, std::placeholders::_2),
    std::bind(
      &Local<dim>::bdry, local, std::placeholders::_1, std::placeholders::_2),
    std::bind(&Local<dim>::face,
              local,
              std::placeholders::_1,
              std::placeholders::_2,
              std::placeholders::_3,
              std::placeholders::_4),
    assembler,
    lctrl);

  deallog << "  Results cells";
  for (unsigned int i = 0; i < n_functionals; ++i)
    deallog << '\t' << cells.block(i).l1_norm();
  deallog << std::endl;

  deallog << "  Results faces";
  for (unsigned int i = 0; i < n_functionals; ++i)
    deallog << '\t' << faces.block(i).l1_norm();
  deallog << std::endl;

  cells = 0.;
  faces = 0.;

  MeshWorker::DoFInfo<dim> mg_dof_info(mgdofs);
  MeshWorker::loop<dim, dim, MeshWorker::DoFInfo<dim>, EmptyInfoBox>(
    mgdofs.begin_mg(),
    mgdofs.end_mg(),
    mg_dof_info,
    info_box,
    std::bind(
      &Local<dim>::cell, local, std::placeholders::_1, std::placeholders::_2),
    std::bind(
      &Local<dim>::bdry, local, std::placeholders::_1, std::placeholders::_2),
    std::bind(&Local<dim>::face,
              local,
              std::placeholders::_1,
              std::placeholders::_2,
              std::placeholders::_3,
              std::placeholders::_4),
    assembler,
    lctrl);

  deallog << "MGResults cells";
  for (unsigned int i = 0; i < n_functionals; ++i)
    deallog << '\t' << cells.block(i).l1_norm();
  deallog << std::endl;

  deallog << "MGResults faces";
  for (unsigned int i = 0; i < n_functionals; ++i)
    deallog << '\t' << faces.block(i).l1_norm();
  deallog << std::endl;
}


template <int dim>
void
test(const FiniteElement<dim> &fe)
{
  Triangulation<dim> tr(Triangulation<dim>::limit_level_difference_at_vertices);
  DoFHandler<dim>    dofs(tr);
  GridGenerator::hyper_cube(tr);
  tr.refine_global(1);
  deallog.push("1");
  dofs.distribute_dofs(fe);
  dofs.distribute_mg_dofs();
  test_mesh(dofs);
  deallog.pop();
  tr.begin(1)->set_refine_flag();
  tr.execute_coarsening_and_refinement();
  deallog.push("2");
  dofs.distribute_dofs(fe);
  dofs.distribute_mg_dofs();
  test_mesh(dofs);
  deallog.pop();
  tr.begin(2)->set_refine_flag();
  tr.execute_coarsening_and_refinement();
  deallog.push("3");
  dofs.distribute_dofs(fe);
  dofs.distribute_mg_dofs();
  test_mesh(dofs);
  deallog.pop();
}


int
main()
{
  const std::string logname = "output";
  std::ofstream     logfile(logname);
  deallog.attach(logfile);

  FE_DGP<2> el2(0);
  FE_DGP<3> el3(0);

  deallog.push("2D");
  test(el2);
  deallog.pop();
  deallog.push("3D");
  test(el3);
  deallog.pop();
}
