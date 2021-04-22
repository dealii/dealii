// ---------------------------------------------------------------------
//
// Copyright (C) 2007 - 2018 by the deal.II authors
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



// test mesh_loop on 1d meshes with local refinement

#include <deal.II/base/work_stream.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/meshworker/assemble_flags.h>
#include <deal.II/meshworker/mesh_loop.h>

#include <fstream>

#include "../tests.h"

struct ScratchData
{};
struct CopyData
{};

using namespace MeshWorker;

template <int dim, int spacedim>
void
test()
{
  Triangulation<dim, spacedim> tria;
  GridGenerator::hyper_cube(tria);
  tria.refine_global(1);
  tria.begin_active()->set_refine_flag();
  tria.execute_coarsening_and_refinement();

  ScratchData scratch;
  CopyData    copy;

  auto cell = tria.begin_active();
  auto endc = tria.end();

  using Iterator = decltype(cell);



  auto cell_worker = [](const Iterator &cell, ScratchData &s, CopyData &c) {
    deallog << "Cell worker on : " << cell << std::endl;
  };

  auto boundary_worker =
    [](const Iterator &cell, const unsigned int &f, ScratchData &, CopyData &) {
      deallog << "Boundary worker on : " << cell << ", Face : " << f
              << std::endl;
    };

  auto face_worker = [](const Iterator &    cell,
                        const unsigned int &f,
                        const unsigned int &sf,
                        const Iterator &    ncell,
                        const unsigned int &nf,
                        const unsigned int &nsf,
                        ScratchData &       s,
                        CopyData &          c) {
    deallog << "Face worker on : " << cell << ", Neighbor cell : " << ncell
            << ", Face : " << f << ", Neighbor Face : " << nf
            << ", Subface: " << sf << ", Neighbor Subface: " << nsf
            << std::endl;
  };

  auto copier = [](const CopyData &) { deallog << "copier" << std::endl; };

  deallog << "ONE DIMENSION" << std::endl << std::endl;

  mesh_loop(cell,
            endc,
            cell_worker,
            copier,
            scratch,
            copy,
            MeshWorker::assemble_own_cells |
              MeshWorker::assemble_own_interior_faces_once |
              MeshWorker::assemble_boundary_faces,
            boundary_worker,
            face_worker);
}


int
main()
{
  initlog();
  MultithreadInfo::set_thread_limit(1); // to make output deterministic

  test<1, 1>();
}
