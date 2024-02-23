// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// A test for Issue #14183. On a very coarse mesh, the
// vertices chosen to compute the normal vector may
// lead to two parallel tangent vectors, resulting
// in a zero normal vector. This can be fixed by
// choosing other combinations of vertices when the
// normal vector is too close to a zero vector.

#include <deal.II/base/exceptions.h>
#include <deal.II/base/function.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/vector_tools_constraints.h>

#include "../tests.h"

int
main()
{
  initlog();
  deallog.get_file_stream().precision(7);
  deallog.get_file_stream().setf(std::ios::fixed);

  {
    const unsigned int dim = 3;
    Triangulation<dim> triangulation(
      Triangulation<dim>::limit_level_difference_at_vertices);
    GridGenerator::half_hyper_shell(triangulation, Point<dim>(), 0.5, 1.);

    std::set<types::boundary_id> no_flux_boundary{0};
    MappingQ<dim>                mapping(2);

    FESystem<dim>   fe(FE_Q<dim>(2), dim);
    DoFHandler<dim> dofh(triangulation);

    dofh.distribute_dofs(fe);

    AffineConstraints<double> constraints;


    const unsigned int                 face_no = 0;
    const std::vector<Point<dim - 1>> &unit_support_points =
      fe.get_unit_face_support_points(face_no);

    Assert(unit_support_points.size() == fe.n_dofs_per_face(),
           ExcInternalError());


    Quadrature<dim - 1> face_quadrature(unit_support_points);

    FEFaceValues<dim> fe_face_values(mapping,
                                     fe,
                                     face_quadrature,
                                     update_quadrature_points |
                                       update_normal_vectors);


    std::set<types::boundary_id>::iterator b_id;

    for (const auto &cell : dofh.active_cell_iterators())
      if (!cell->is_artificial())
        for (const unsigned int face_no : cell->face_indices())
          if ((b_id = no_flux_boundary.find(
                 cell->face(face_no)->boundary_id())) != no_flux_boundary.end())
            {
              fe_face_values.reinit(cell, face_no);
              for (unsigned int i = 0; i < face_quadrature.size(); ++i)
                Tensor<1, dim> normal_vector =
                  (cell->face(face_no)->get_manifold().normal_vector(
                    cell->face(face_no), fe_face_values.quadrature_point(i)));
            }

    deallog << " OK! " << std::endl;
  }
}
