// ---------------------------------------------------------------------
//
// Copyright (C) 2017 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------


// check that the MeshWorker framework is instantiated for all dim<=spacedim
// modified from a bug report submitted by Andrea Bonito

#include <deal.II/base/function_lib.h>

#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/mapping.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/meshworker/assembler.h>
#include <deal.II/meshworker/dof_info.h>
#include <deal.II/meshworker/integration_info.h>
#include <deal.II/meshworker/loop.h>

#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"



template <int dim, int spacedim>
class TestIntegrator : public MeshWorker::LocalIntegrator<dim, spacedim>
{
public:
  TestIntegrator(){};

  void
  cell(MeshWorker::DoFInfo<dim, spacedim> &                 dinfo,
       typename MeshWorker::IntegrationInfo<dim, spacedim> &info) const {};

  void
  boundary(MeshWorker::DoFInfo<dim, spacedim> &                 dinfo,
           typename MeshWorker::IntegrationInfo<dim, spacedim> &info) const {};

  void
  face(MeshWorker::DoFInfo<dim, spacedim> &                 dinfo1,
       MeshWorker::DoFInfo<dim, spacedim> &                 dinfo2,
       typename MeshWorker::IntegrationInfo<dim, spacedim> &info1,
       typename MeshWorker::IntegrationInfo<dim, spacedim> &info2) const {};
};

template <int dim, int spacedim>
void
test()
{
  MappingQ1<dim, spacedim> mapping;

  Triangulation<dim, spacedim> triangulation;
  GridGenerator::hyper_cube(triangulation);

  FE_Q<dim, spacedim>       fe(1);
  DoFHandler<dim, spacedim> dof_handler(triangulation);
  dof_handler.distribute_dofs(fe);

  SparseMatrix<double> matrix;

  MeshWorker::IntegrationInfoBox<dim, spacedim> info_box;
  info_box.initialize(fe, mapping);

  MeshWorker::DoFInfo<dim, spacedim> dof_info(dof_handler);

  MeshWorker::Assembler::MatrixSimple<SparseMatrix<double>> assembler;

  assembler.initialize(matrix);

  TestIntegrator<dim, spacedim> integrator;

  MeshWorker::integration_loop<dim, spacedim>(dof_handler.begin_active(),
                                              dof_handler.end(),
                                              dof_info,
                                              info_box,
                                              integrator,
                                              assembler);
  deallog << "dim = " << dim << ", spacedim = " << spacedim << ": OK"
          << std::endl;
}


int
main()
{
  initlog();

  test<1, 1>();
  test<1, 2>();
  test<2, 2>();
  test<1, 3>();
  test<2, 3>();
  test<3, 3>();

  return 0;
}
