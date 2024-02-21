// ------------------------------------------------------------------------
//
// Copyright (C) 2020 - 2021 by the deal.II Authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Test mesh_loop with anisotropic grids in 2 and 3 dimensions

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/meshworker/mesh_loop.h>

#include <fstream>

#include "../tests.h"



// Create 2d grid
void
create_grid(Triangulation<2> &tria)
{
  GridGenerator::hyper_cube(tria);

  tria.begin_active()->set_refine_flag(RefinementCase<2>::cut_x);
  tria.execute_coarsening_and_refinement();
  tria.begin_active()->set_refine_flag(RefinementCase<2>::cut_x);
  tria.execute_coarsening_and_refinement();

  {
    auto it = tria.begin();
    std::advance(it, 4);
    it->set_refine_flag(RefinementCase<2>::cut_x);
    tria.execute_coarsening_and_refinement();
  }

  {
    auto it = tria.begin();
    std::advance(it, 6);
    it->set_refine_flag(RefinementCase<2>::cut_y);
    tria.execute_coarsening_and_refinement();
  }

  {
    auto it = tria.begin();
    std::advance(it, 8);
    it->set_refine_flag(RefinementCase<2>::cut_x);
    tria.execute_coarsening_and_refinement();
  }

  {
    auto it = tria.begin();
    std::advance(it, 2);
    it->set_refine_flag(RefinementCase<2>::cut_y);
    tria.execute_coarsening_and_refinement();
  }
}


// Create 3d grid
void
create_grid(Triangulation<3> &tria)
{
  GridGenerator::hyper_cube<3>(tria);
  tria.begin_active()->set_refine_flag(dealii::RefinementCase<3>::cut_x);
  tria.execute_coarsening_and_refinement();
  tria.begin_active()->set_refine_flag(dealii::RefinementCase<3>::cut_y);
  tria.execute_coarsening_and_refinement();
  {
    auto it = tria.begin();
    std::advance(it, 2);
    it->set_refine_flag(dealii::RefinementCase<3>::cut_z);
    tria.execute_coarsening_and_refinement();
  }
}



template <int dim>
void
test()
{
  deallog << "Testing dim = " << dim << std::endl
          << "===============" << std::endl
          << std::endl;
  Triangulation<dim> tria;
  create_grid(tria);

  struct ScratchData
  {
    unsigned int foo;
  };
  struct CopyData
  {
    unsigned int foo;
  };

  ScratchData scratch;
  CopyData    copy;

  FE_Q<dim>       fe(1);
  DoFHandler<dim> dh(tria);
  dh.distribute_dofs(fe);

  auto workerFoo = [](const typename DoFHandler<dim>::active_cell_iterator &,
                      ScratchData &,
                      CopyData &) {};
  auto copierFoo = [](const CopyData &) {};
  auto faceFoo = [](const typename DoFHandler<dim>::active_cell_iterator &cell,
                    const unsigned int                                    f,
                    const unsigned int                                    sf,
                    const typename DoFHandler<dim>::active_cell_iterator &ncell,
                    const unsigned int                                    nf,
                    const unsigned int                                    nsf,
                    ScratchData &,
                    CopyData &) {
    deallog << "Cell (+): " << cell << " -- " << f << ", " << (int)sf
            << std::endl
            << "Cell (-): " << ncell << " -- " << nf << ", " << (int)nsf
            << std::endl
            << std::endl;
  };
  auto boundaryFoo = [](const typename DoFHandler<dim>::active_cell_iterator &,
                        const unsigned int,
                        ScratchData &,
                        CopyData &) {};
  std::function<void(const typename DoFHandler<dim>::active_cell_iterator &,
                     ScratchData &,
                     CopyData &)>
    emptyCellWorker(workerFoo);
  std::function<void(const typename DoFHandler<dim>::active_cell_iterator &,
                     const unsigned int,
                     const unsigned int,
                     const typename DoFHandler<dim>::active_cell_iterator &,
                     const unsigned int,
                     const unsigned int,
                     ScratchData &,
                     CopyData &)>
    emptyFaceWorker(faceFoo);
  std::function<void(const typename DoFHandler<dim>::active_cell_iterator &,
                     const unsigned int,
                     ScratchData &,
                     CopyData &)>
                                        emptyBoundaryWorker(boundaryFoo);
  std::function<void(const CopyData &)> emptyCopier(copierFoo);

  MeshWorker::AssembleFlags flags;
  flags = MeshWorker::AssembleFlags::assemble_own_cells |
          MeshWorker::AssembleFlags::assemble_own_interior_faces_once |
          MeshWorker::AssembleFlags::assemble_boundary_faces;

  MeshWorker::mesh_loop(dh.begin_active(),
                        dh.end(),
                        emptyCellWorker,
                        emptyCopier,
                        scratch,
                        copy,
                        flags,
                        emptyBoundaryWorker,
                        emptyFaceWorker);
}



int
main()
{
  initlog();
  MultithreadInfo::set_thread_limit(1);
  test<2>();
  test<3>();
}
