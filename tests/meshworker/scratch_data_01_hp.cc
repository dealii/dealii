/* ------------------------------------------------------------------------
 *
 * SPDX-License-Identifier: LGPL-2.1-or-later
 * Copyright (C) 2019 - 2023 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * Part of the source code is dual licensed under Apache-2.0 WITH
 * LLVM-exception OR LGPL-2.1-or-later. Detailed license information
 * governing the source code and code contributions can be found in
 * LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
 *
 * ------------------------------------------------------------------------
 */

// Test ScratchData on simple volume computation
// - hp variant

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/hp/fe_collection.h>
#include <deal.II/hp/q_collection.h>

#include <deal.II/meshworker/copy_data.h>
#include <deal.II/meshworker/mesh_loop.h>
#include <deal.II/meshworker/scratch_data.h>

#include <cmath>
#include <fstream>
#include <iostream>
#include <unordered_map>

#include "../tests.h"

using namespace MeshWorker;

template <int dim, int spacedim>
void
test()
{
  Triangulation<dim, spacedim>    tria;
  hp::FECollection<dim, spacedim> fe(FE_Q<dim, spacedim>(1),
                                     FE_Q<dim, spacedim>(2));
  DoFHandler<dim, spacedim>       dh(tria);
  GridGenerator::hyper_cube(tria);
  tria.refine_global(1);
  tria.begin_active()->set_refine_flag();
  tria.execute_coarsening_and_refinement();

  const auto bounding_box = GridTools::compute_bounding_box(tria);
  for (auto &cell : dh.active_cell_iterators())
    {
      bool use_alt_index = true;
      for (unsigned int d = 0; d < spacedim; ++d)
        if (cell->center()[d] < bounding_box.center()[d])
          {
            use_alt_index = false;
            break;
          }

      if (use_alt_index)
        cell->set_active_fe_index(1);
      else
        cell->set_active_fe_index(0);
    }

  dh.distribute_dofs(fe);

  hp::QCollection<dim>     cell_quad(QGauss<dim>(3), QGauss<dim>(4));
  hp::QCollection<dim - 1> face_quad(QGauss<dim - 1>(3), QGauss<dim - 1>(4));

  UpdateFlags cell_flags = update_JxW_values;
  UpdateFlags face_flags = update_JxW_values;

  using ScratchData = MeshWorker::ScratchData<dim, spacedim>;
  using CopyData    = MeshWorker::CopyData<0, 1, 0>;

  ScratchData scratch(fe, cell_quad, cell_flags, face_quad, face_flags);
  CopyData    copy(3); // store measures

  auto cell = dh.begin_active();
  auto endc = dh.end();

  using Iterator = decltype(cell);

  auto measures = std::make_tuple(0.0, 0.0, 0.0);


  auto cell_worker = [](const Iterator &cell, ScratchData &s, CopyData &c) {
    const auto &fev = s.reinit(cell);
    const auto &JxW = s.get_JxW_values();
    for (auto w : JxW)
      c.vectors[0][0] += w;
  };

  auto boundary_worker = [](const Iterator     &cell,
                            const unsigned int &f,
                            ScratchData        &s,
                            CopyData           &c) {
    const auto &fev = s.reinit(cell, f);
    const auto &JxW = s.get_JxW_values();
    for (auto w : JxW)
      c.vectors[0][1] += w;
  };

  auto face_worker = [](const Iterator     &cell,
                        const unsigned int &f,
                        const unsigned int &sf,
                        const Iterator     &ncell,
                        const unsigned int &nf,
                        const unsigned int &nsf,
                        ScratchData        &s,
                        CopyData           &c) {
    const auto &fev  = s.reinit(cell, f, sf);
    const auto &nfev = s.reinit_neighbor(ncell, nf, nsf);

    const auto &JxW  = s.get_JxW_values();
    const auto &nJxW = s.get_neighbor_JxW_values();
    for (auto w : JxW)
      c.vectors[0][2] += w;
  };

  auto copier = [&measures](const CopyData &c) {
    std::get<0>(measures) += c.vectors[0][0];
    std::get<1>(measures) += c.vectors[0][1];
    std::get<2>(measures) += c.vectors[0][2];
  };

  mesh_loop(cell,
            endc,
            cell_worker,
            copier,
            scratch,
            copy,
            assemble_own_cells | assemble_boundary_faces |
              assemble_own_interior_faces_once,
            boundary_worker,
            face_worker);

  deallog << "Testing <" << dim << "," << spacedim << ">" << std::endl;

  deallog << "Volume: " << std::get<0>(measures) << std::endl;

  deallog << "Boundary surface: " << std::get<1>(measures) << std::endl;

  deallog << "Interior surface: " << std::get<2>(measures) << std::endl;
}



int
main()
{
  initlog();

  test<2, 2>();
  test<2, 3>();
  test<3, 3>();
}
