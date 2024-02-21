// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2021 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// Verify CellAccessor::neighbor_child_on_subface for a triangle mesh.


#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"


template <int dim>
void
test()
{
  // setup
  //   +---+---+
  //   |\  |  /|
  //   |  \|/  |
  //   |   +---+
  //   |    \  |
  //   |      \|
  //   +-------+
  Triangulation<dim> tria;
  GridGenerator::subdivided_hyper_cube_with_simplices(tria, 1);
  {
    auto cell = tria.begin_active();
    (++cell)->set_refine_flag();
    Assert(++cell == tria.end(), ExcInternalError());
  }
  tria.execute_coarsening_and_refinement();

  // find face index on unrefined cell to neighboring cells
  const auto  &unrefined_cell = tria.begin_active(0);
  unsigned int unrefined_f    = numbers::invalid_unsigned_int;
  for (unsigned int f = 0; f < unrefined_cell->n_faces(); ++f)
    if (!unrefined_cell->face(f)->at_boundary())
      unrefined_f = f;

  // verify whether unrefined cell and neighboring children have matching
  // vertices on their corresponding subface
  for (unsigned int sf = 0; sf < GeometryInfo<dim>::max_children_per_face; ++sf)
    {
      // unrefined vertex on subface
      const Point<dim> &vertex_unrefined =
        unrefined_cell->face(unrefined_f)->vertex(sf);

      // child on subface [! focus of this particular test !]
      const auto &neighboring_child =
        unrefined_cell->neighbor_child_on_subface(unrefined_f, sf);

      // face of child
      unsigned int neighboring_child_f = numbers::invalid_unsigned_int;
      for (unsigned int f = 0; f < neighboring_child->n_faces(); ++f)
        if (!neighboring_child->face(f)->at_boundary() &&
            neighboring_child->neighbor(f) == unrefined_cell)
          neighboring_child_f = f;
      const auto &neighboring_child_face =
        neighboring_child->face(neighboring_child_f);

      // find unrefined vertex among child face vertices
      bool vertex_found = false;
      for (unsigned int v = 0; v < neighboring_child_face->n_vertices(); ++v)
        if (neighboring_child_face->vertex(v) == vertex_unrefined)
          vertex_found = true;

      Assert(vertex_found,
             ExcMessage(
               "Wrong assignment of neighboring children to subfaces."));
    }

  deallog << "OK" << std::endl;
}


int
main()
{
  initlog();

  deallog.push("2d");
  test<2>();
  deallog.pop();
}
