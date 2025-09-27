// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2014 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// test meshworker LoopControl
// variation of mesh_worker_01 with more cpus and cells

#include <deal.II/distributed/tria.h>

#include <deal.II/fe/fe_dgp.h>

#include <deal.II/grid/filtered_iterator.h>
#include <deal.II/grid/grid_generator.h>

#include <deal.II/meshworker/assembler.h>
#include <deal.II/meshworker/loop.h>

#include "../tests.h"



template <int dim>
class myIntegrator : public dealii::MeshWorker::LocalIntegrator<dim>
{
public:
  using CellInfo = MeshWorker::IntegrationInfo<dim>;

  void
  cell(MeshWorker::DoFInfo<dim> &dinfo, CellInfo &info) const;
  void
  boundary(MeshWorker::DoFInfo<dim> &dinfo, CellInfo &info) const;
  void
  face(MeshWorker::DoFInfo<dim> &dinfo1,
       MeshWorker::DoFInfo<dim> &dinfo2,
       CellInfo                 &info1,
       CellInfo                 &info2) const;
};

template <int dim>
void
myIntegrator<dim>::cell(MeshWorker::DoFInfo<dim> &info, CellInfo &) const
{
  deallog << "C " << info.cell->id() << std::endl;
}


template <int dim>
void
myIntegrator<dim>::boundary(MeshWorker::DoFInfo<dim> &info, CellInfo &) const
{
  // deallog << "B cell = " << info.cell->id() << " face = " << info.face_number
  // << std::endl;
}


template <int dim>
void
myIntegrator<dim>::face(MeshWorker::DoFInfo<dim> &info1,
                        MeshWorker::DoFInfo<dim> &info2,
                        CellInfo &,
                        CellInfo &) const
{
  deallog << "F cell1 = " << info1.cell->id() << " face = " << info1.face_number
          << " cell2 = " << info2.cell->id() << " face2 = " << info2.face_number
          << std::endl;
}


class DoNothingAssembler
{
public:
  template <class DOFINFO>
  void
  initialize_info(DOFINFO &info, bool face) const
  {}
  template <class DOFINFO>
  void
  assemble(const DOFINFO &info)
  {}
  template <class DOFINFO>
  void
  assemble(const DOFINFO &info1, const DOFINFO &info2)
  {}
};

template <int dim>
void
test_simple(DoFHandler<dim> &dofs, MeshWorker::LoopControl &lctrl)
{
  myIntegrator<dim>                   local;
  DoNothingAssembler                  assembler;
  MeshWorker::IntegrationInfoBox<dim> info_box;

  MeshWorker::DoFInfo<dim> dof_info(dofs.block_info());

  //  integration_loop(ITERATOR begin,
  //                std_cxx20::type_identity_t<ITERATOR> end,
  //                DOFINFO &dinfo,
  //                INFOBOX &info,
  //                const std::function<void (DOFINFO &, typename
  //                INFOBOX::CellInfo &)> &cell_worker, const std::function<void
  //                (DOFINFO &, typename INFOBOX::CellInfo &)> &boundary_worker,
  //                const std::function<void (DOFINFO &, DOFINFO &,
  //                                                typename INFOBOX::CellInfo
  //                                                &, typename
  //                                                INFOBOX::CellInfo &)>
  //                                                &face_worker,
  //                ASSEMBLER &assembler,
  //                const LoopControl &lctrl)
  //


  MeshWorker::integration_loop<dim,
                               dim,
                               typename DoFHandler<dim>::active_cell_iterator,
                               DoNothingAssembler>(dofs.begin_active(),
                                                   dofs.end(),
                                                   dof_info,
                                                   info_box,
                                                   local,
                                                   assembler,
                                                   lctrl);

  //  MeshWorker::loop<dim, dim, MeshWorker::DoFInfo<dim>,
  //  MeshWorker::IntegrationInfoBox<dim> >
  //    (dofs.begin_active(), dofs.end(),
  //   dof_info, info_box,
  //       std::bind (&Integrator<dim>::cell, local, std::placeholders::_1,
  //       std::placeholders::_2),
  //   std::bind (&Integrator<dim>::bdry, local, std::placeholders::_1,
  //   std::placeholders::_2), std::bind (&Integrator<dim>::face, local,
  //   std::placeholders::_1, std::placeholders::_2, std::placeholders::_3,
  //   std::placeholders::_4),
  //     local,
  //     lctrl);
}

template <int dim>
void
test_loop(DoFHandler<dim> &dofs, MeshWorker::LoopControl &lctrl)
{
  deallog << "* own_cells=" << lctrl.own_cells
          << " ghost_cells=" << lctrl.ghost_cells
          << " own_faces=" << lctrl.own_faces
          << " faces_to_ghost=" << lctrl.faces_to_ghost << std::endl;
  test_simple(dofs, lctrl);
}

template <int dim>
void
test()
{
  parallel::distributed::Triangulation<dim> tr(MPI_COMM_WORLD,
                                               Triangulation<dim>::none/*,
                   parallel::distributed::Triangulation<dim>::construct_multigrid_hierarchy*/);
  GridGenerator::hyper_cube(tr);
  tr.refine_global(2);

  FE_DGP<dim> fe(0);

  DoFHandler<dim> dofs(tr);
  dofs.distribute_dofs(fe);

  dofs.initialize_local_block_info();
  deallog << "DoFHandler ndofs=" << dofs.n_dofs() << std::endl;

  MeshWorker::LoopControl lctrl;

  deallog << "*** 1. CELLS ***" << std::endl;
  /*
  lctrl.own_faces = MeshWorker::LoopControl::never;
  lctrl.faces_to_ghost = MeshWorker::LoopControl::never;

  lctrl.own_cells = false; lctrl.ghost_cells = false;
  test_loop(dofs, lctrl);

  lctrl.own_cells = true; lctrl.ghost_cells = false;
  test_loop(dofs, lctrl);

  lctrl.own_cells = false; lctrl.ghost_cells = true;
  test_loop(dofs, lctrl);

  lctrl.own_cells = true; lctrl.ghost_cells = true;
  test_loop(dofs, lctrl);
  */
  deallog << "*** 2. FACES ***" << std::endl;

  lctrl.own_cells   = false;
  lctrl.ghost_cells = false;

  lctrl.own_faces      = MeshWorker::LoopControl::one;
  lctrl.faces_to_ghost = MeshWorker::LoopControl::never;
  test_loop(dofs, lctrl);

  lctrl.own_faces      = MeshWorker::LoopControl::both;
  lctrl.faces_to_ghost = MeshWorker::LoopControl::never;
  test_loop(dofs, lctrl);

  lctrl.own_faces      = MeshWorker::LoopControl::never;
  lctrl.faces_to_ghost = MeshWorker::LoopControl::one;
  test_loop(dofs, lctrl);

  lctrl.own_faces      = MeshWorker::LoopControl::never;
  lctrl.faces_to_ghost = MeshWorker::LoopControl::both;
  test_loop(dofs, lctrl);

  //
  //
  //  for (int gc=0;gc<2;gc++)
  //  for (int oc=0;oc<2;oc++)
  //  for (int of=0;of<3;of++)
  //    {
  //      lctrl.own_cells = !!oc;
  //      lctrl.ghost_cells = !!gc;
  //
  //      lctrl.own_faces = (MeshWorker::LoopControl::FaceOption)of;
  //      test_loop(dofs, lctrl);
  //    }
}


int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    log;

  test<2>();
}
