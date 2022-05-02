// ---------------------------------------------------------------------
//
// Copyright (C) 2010 - 2020 by the deal.II authors
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


// verify that we can indeed call extract_boundary_mesh with DoFHandler
// arguments. based on a test by Korosh Taebi

#include <deal.II/base/function.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q_eulerian.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/vector_tools.h>

#include <iostream>

#include "../tests.h"

namespace Step38
{
  template <int spacedim>
  class Extract_Mesh_Test
  {
  public:
    Extract_Mesh_Test();
    void
    run();

  private:
    static const unsigned int boundary_dim = spacedim - 1;

    Triangulation<spacedim>               volume_mesh_triangulation;
    Triangulation<boundary_dim, spacedim> boundary_triangulation;

    FE_Q<spacedim>               space_fe;
    FE_Q<boundary_dim, spacedim> boundary_fe;

    DoFHandler<spacedim>               space_dof_handler;
    DoFHandler<boundary_dim, spacedim> contact_dof_handler;
  };


  template <int spacedim>
  Extract_Mesh_Test<spacedim>::Extract_Mesh_Test()
    : space_fe(spacedim)
    , boundary_fe(1)
    , space_dof_handler(volume_mesh_triangulation)
    , contact_dof_handler(boundary_triangulation)
  {}

  template <int spacedim>
  void
  Extract_Mesh_Test<spacedim>::run()
  {
    GridGenerator::hyper_cube(volume_mesh_triangulation);
    std::set<types::boundary_id> boundary_ids;
    boundary_ids.insert(0);

    space_dof_handler.distribute_dofs(space_fe);

    std::map<typename DoFHandler<boundary_dim, spacedim>::cell_iterator,
             typename DoFHandler<spacedim>::face_iterator>
      element_assignment =
        GridGenerator::extract_boundary_mesh(space_dof_handler,
                                             contact_dof_handler,
                                             boundary_ids);

    contact_dof_handler.distribute_dofs(boundary_fe);

    typename std::map<
      typename DoFHandler<boundary_dim, spacedim>::cell_iterator,
      typename DoFHandler<spacedim>::face_iterator>::iterator Iterator;

    for (Iterator = element_assignment.begin();
         Iterator != element_assignment.end();
         ++Iterator)
      {
        deallog << "element_assignment maps " << Iterator->first << " onto "
                << Iterator->second << std::endl;
      }
  }
} // namespace Step38


int
main()
{
  initlog();

  {
    using namespace Step38;

    Extract_Mesh_Test<2> Test;
    Test.run();
  }
}
