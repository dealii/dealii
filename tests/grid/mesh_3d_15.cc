// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2003 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// make sure TriaCellAccessor::neighbor_child_on_subface does what it
// is supposed to do. check it for dof accessors, since they simply
// call the tria accessors, and this way we catch both cases at the
// same time

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/vector.h>

#include "../tests.h"

#include "mesh_3d.h"



void
check_this(Triangulation<3> &tria)
{
  FE_Q<3> fe(1);

  DoFHandler<3> dof_handler(tria);
  dof_handler.distribute_dofs(fe);

  DoFHandler<3>::active_cell_iterator cell = dof_handler.begin_active();
  for (; cell != dof_handler.end(); ++cell)
    for (const unsigned int face_no : GeometryInfo<3>::face_indices())
      if (!cell->at_boundary(face_no) &&
          cell->neighbor(face_no)->has_children())
        for (unsigned int subface_no = 0;
             subface_no < GeometryInfo<3>::max_children_per_face;
             ++subface_no)
          {
            // get an iterator
            // pointing to the cell
            // behind the present
            // subface

            // way a) construct it ourselves,
            // considering orientation and
            // rotation of the face
            const DoFHandler<3>::cell_iterator neighbor =
              cell->neighbor(face_no);
            const unsigned int neighbor_neighbor =
              cell->neighbor_of_neighbor(face_no);
            const unsigned int neighbor_child_index =
              (GeometryInfo<3>::child_cell_on_face(
                RefinementCase<3>::isotropic_refinement,
                neighbor_neighbor,
                (subface_no),
                neighbor->face_orientation(neighbor_neighbor)));
            const DoFHandler<3>::active_cell_iterator neighbor_child =
              neighbor->child(neighbor_child_index);

            // way b) use the convenient
            // function
            // neighbor_child_on_subface

            // make sure, that both ways yield
            // the same result
            AssertThrow(neighbor_child ==
                          cell->neighbor_child_on_subface(face_no, subface_no),
                        ExcInternalError());
          }
}



void
check(Triangulation<3> &tria)
{
  (std::next(tria.begin_active()))->set_refine_flag();
  tria.execute_coarsening_and_refinement();

  deallog << "Initial check" << std::endl;
  check_this(tria);
  // TODO:[WB] Is there a reason to do this three times?
  // Changed to two. Guido
  for (unsigned int r = 0; r < 2; ++r)
    {
      tria.refine_global(1);
      deallog << "Check " << r << std::endl;
      check_this(tria);
    }

  coarsen_global(tria);
  deallog << "Check " << 1 << std::endl;
  check_this(tria);

  tria.refine_global(1);
  deallog << "Check " << 2 << std::endl;
  check_this(tria);
}


int
main()
{
  initlog();

  {
    Triangulation<3> coarse_grid;
    create_two_cubes(coarse_grid);
    check(coarse_grid);
  }

  {
    Triangulation<3> coarse_grid;
    create_L_shape(coarse_grid);
    check(coarse_grid);
  }

  {
    Triangulation<3> coarse_grid;
    GridGenerator::hyper_ball(coarse_grid);
    check(coarse_grid);
  }
}
