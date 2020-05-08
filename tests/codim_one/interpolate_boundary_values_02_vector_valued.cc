// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2020 by the deal.II authors
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



// test VectorTools::interpolate_boundary_values for codim=1. like _02
// but for vector-valued elements

#include <deal.II/base/function_lib.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/numerics/vector_tools.h>

#include <string>

#include "../tests.h"

template <int dim>
class X : public Function<dim>
{
public:
  X()
    : Function<dim>(dim)
  {}

  double
  value(const Point<dim> &p, const unsigned int component) const
  {
    return p[component];
  }
};

void
test()
{
  const int dim      = 1;
  const int spacedim = 2;

  Triangulation<dim, spacedim> tria;
  GridGenerator::hyper_cube(tria);
  deallog << tria.n_active_cells() << " active cells" << std::endl;

  FESystem<dim, spacedim>   fe(FE_Q<dim, spacedim>(2), spacedim);
  DoFHandler<dim, spacedim> dof_handler(tria);
  dof_handler.distribute_dofs(fe);

  deallog << dof_handler.n_dofs() << " degrees of freedom" << std::endl;

  // test left and right boundary
  // separatel
  for (unsigned int boundary_id = 0; boundary_id < 2; ++boundary_id)
    {
      std::map<types::global_dof_index, double> bv;
      VectorTools::interpolate_boundary_values(dof_handler,
                                               boundary_id,
                                               X<spacedim>(),
                                               bv);
      deallog << bv.size() << " boundary degrees of freedom" << std::endl;

      for (std::map<types::global_dof_index, double>::const_iterator i =
             bv.begin();
           i != bv.end();
           ++i)
        deallog << i->first << ' ' << i->second << std::endl;

      for (DoFHandler<dim, spacedim>::active_cell_iterator cell =
             dof_handler.begin_active();
           cell != dof_handler.end();
           ++cell)
        for (const unsigned int f : GeometryInfo<dim>::face_indices())
          if (cell->at_boundary(f) &&
              (cell->face(f)->boundary_id() == boundary_id))
            for (unsigned int v = 0; v < GeometryInfo<dim>::vertices_per_face;
                 ++v)
              for (unsigned int i = 0; i < fe.dofs_per_vertex; ++i)
                {
                  AssertThrow(bv.find(cell->face(f)->vertex_dof_index(v, i)) !=
                                bv.end(),
                              ExcInternalError());
                  AssertThrow(bv[cell->face(f)->vertex_dof_index(v, i)] ==
                                X<spacedim>().value(cell->face(f)->vertex(v),
                                                    i),
                              ExcInternalError());
                }
    }
}



int
main()
{
  initlog();

  test();

  return 0;
}
