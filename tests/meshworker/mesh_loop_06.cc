/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2019 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * The deal.II library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE.md at
 * the top level directory of deal.II.
 *
 * ---------------------------------------------------------------------
 */

// test mesh_loop for filtered iterator ranges
// This test is based off base/work_stream_06.cc, but for filtered iterator
// ranges. In particular, we test the material ID predicate.

#include <deal.II/base/iterator_range.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/filtered_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

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
test(const types::material_id test_id)
{
  Assert(test_id == 0 || test_id == 1,
         ExcMessage("Only material IDs 0 and 1 are valid."));

  Triangulation<dim, spacedim> tria;
  FE_Q<dim, spacedim>          fe(1);
  DoFHandler<dim, spacedim>    dh(tria);
  GridGenerator::hyper_cube(tria);
  tria.refine_global(1);
  tria.begin_active()->set_material_id(1);
  tria.begin_active()->set_refine_flag();
  tria.execute_coarsening_and_refinement();
  dh.distribute_dofs(fe);

  QGauss<dim>     quad(3);
  QGauss<dim - 1> face_quad(3);

  UpdateFlags cell_flags = update_JxW_values;
  UpdateFlags face_flags = update_JxW_values;

  using ScratchData = MeshWorker::ScratchData<dim, spacedim>;
  using CopyData    = MeshWorker::CopyData<0, 1, 0>;

  ScratchData scratch(fe, quad, cell_flags, face_quad, face_flags);
  CopyData    copy(3); // store measures

  using Iterator = decltype(dh.begin_active());

  auto measures = std::make_tuple(0.0, 0.0, 0.0);


  auto cell_worker = [](const Iterator &cell, ScratchData &s, CopyData &c) {
    const auto &fev = s.reinit(cell);
    const auto &JxW = s.get_JxW_values();
    c               = 0;
    for (auto w : JxW)
      c.vectors[0][0] += w;
  };

  auto boundary_worker = [](const Iterator &    cell,
                            const unsigned int &f,
                            ScratchData &       s,
                            CopyData &          c) {
    const auto &fev = s.reinit(cell, f);
    const auto &JxW = s.get_JxW_values();
    for (auto w : JxW)
      c.vectors[0][1] += w;
  };

  auto face_worker = [](const Iterator &    cell,
                        const unsigned int &f,
                        const unsigned int &sf,
                        const Iterator &    ncell,
                        const unsigned int &nf,
                        const unsigned int &nsf,
                        ScratchData &       s,
                        CopyData &          c) {
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

  const auto filtered_iterator_range =
    filter_iterators(dh.active_cell_iterators(),
                     IteratorFilters::MaterialIdEqualTo(test_id));

  mesh_loop(filtered_iterator_range,
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

  for (auto mat_id : {0, 1})
    {
      deallog.push("Material ID " + Utilities::to_string(mat_id));
      {
        test<2, 2>(mat_id);
        test<2, 3>(mat_id);
        test<3, 3>(mat_id);
      }
      deallog.pop();
    }
}
