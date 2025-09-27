// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2019 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Test assembly using a class and mesh_loop

#include <deal.II/base/work_stream.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/meshworker/assemble_flags.h>
#include <deal.II/meshworker/mesh_loop.h>

#include <fstream>

#include "../tests.h"

using namespace MeshWorker;

struct ScratchData
{};
struct CopyData
{};

template <int dim, int spacedim>
class TestClass
{
public:
  using CellIteratorType =
    typename Triangulation<dim, spacedim>::active_cell_iterator;

  void
  cell_worker(const CellIteratorType &cell, ScratchData &, CopyData &)
  {
    deallog << "Working on cell " << cell << std::endl;
  }

  void
  copier(const CopyData &)
  {}

  void
  boundary_worker(const CellIteratorType &cell,
                  const unsigned int      f,
                  ScratchData &,
                  CopyData &)
  {
    deallog << "Boundary worker on : " << cell << ", Face : " << f << std::endl;
  }

  void
  face_worker(const CellIteratorType &cell,
              const unsigned int      f,
              const unsigned int      sf,
              const CellIteratorType &ncell,
              const unsigned int      nf,
              const unsigned int      nsf,
              ScratchData &,
              CopyData &)
  {
    deallog << "Face worker on : " << cell << ", Neighbor cell : " << ncell
            << ", Face : " << f << ", Neighbor Face : " << nf
            << ", Subface: " << sf << ", Neighbor Subface: " << nsf
            << std::endl;
  }
};

template <int dim, int spacedim>
void
test()
{
  Triangulation<dim, spacedim> tria;
  GridGenerator::hyper_cube(tria);
  tria.refine_global(1);
  tria.begin_active()->set_refine_flag();
  tria.execute_coarsening_and_refinement();

  TestClass<dim, spacedim> test_class;
  ScratchData              scratch;
  CopyData                 copy;

  auto cell = tria.begin_active();
  auto endc = tria.end();

  deallog << "CELLS ONLY" << std::endl;
  mesh_loop(cell,
            endc,
            test_class,
            &TestClass<dim, spacedim>::cell_worker,
            &TestClass<dim, spacedim>::copier,
            scratch,
            copy,
            assemble_own_cells);

  deallog << "CELLS+BOUNDARY" << std::endl;

  mesh_loop(cell,
            endc,
            test_class,
            &TestClass<dim, spacedim>::cell_worker,
            &TestClass<dim, spacedim>::copier,
            scratch,
            copy,
            assemble_own_cells | assemble_boundary_faces,
            &TestClass<dim, spacedim>::boundary_worker);


  deallog << "CELLS+BOUNDARY+FACES" << std::endl;

  mesh_loop(cell,
            endc,
            test_class,
            &TestClass<dim, spacedim>::cell_worker,
            &TestClass<dim, spacedim>::copier,
            scratch,
            copy,
            assemble_own_cells | assemble_boundary_faces |
              assemble_own_interior_faces_once,
            &TestClass<dim, spacedim>::boundary_worker,
            &TestClass<dim, spacedim>::face_worker);
}


int
main()
{
  initlog();
  MultithreadInfo::set_thread_limit(1); // to make output deterministic

  test<2, 2>();
}
