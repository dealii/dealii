// ---------------------------------------------------------------------
//
// Copyright (C) 2012 - 2018 by the deal.II authors
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



// we were computing normal vectors through the boundary object, but the
// boundary object doesn't know about orientation. this leads to trouble if we
// get normal vectors pointing in opposite directions on neighboring cells


#include <deal.II/base/exceptions.h>
#include <deal.II/base/function.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"



template <int dim>
void
run()
{
  Triangulation<dim> triangulation;

  GridIn<dim> gridin;
  gridin.attach_triangulation(triangulation);
  std::ifstream f(SOURCE_DIR "/no_flux_11.msh");
  Assert(f, ExcIO());
  gridin.read_msh(f);

  {
    typename Triangulation<dim>::active_cell_iterator cell = triangulation
                                                               .begin_active(),
                                                      endc =
                                                        triangulation.end();
    for (; cell != endc; ++cell)
      {
        if (cell->is_locally_owned())
          {
            for (const unsigned int face : GeometryInfo<dim>::face_indices())
              {
                if (cell->face(face)->at_boundary())
                  {
                    if ((std::fabs(cell->face(face)->center()(0)) < 0.1) &&
                        (std::fabs(cell->face(face)->center()(dim - 1)) <
                         1e-12))
                      {
                        cell->face(face)->set_boundary_id(1);
                      }

                    if ((std::fabs(cell->face(face)->center()(0)) < 1e-12) &&
                        (std::fabs(cell->face(face)->center()(dim - 1)) < 0.1))
                      {
                        cell->face(face)->set_boundary_id(1);
                      }

                    if ((std::fabs(1.0 - cell->face(face)->center()(0)) <
                         0.1) &&
                        (std::fabs(1.0 - cell->face(face)->center()(dim - 1)) <
                         1e-12))
                      {
                        cell->face(face)->set_boundary_id(2);
                      }

                    if ((std::fabs(1.0 - cell->face(face)->center()(0)) <
                         1e-12) &&
                        (std::fabs(1.0 - cell->face(face)->center()(dim - 1)) <
                         0.1))
                      {
                        cell->face(face)->set_boundary_id(2);
                      }

                    // no normal flux boundary

                    if ((std::fabs(cell->face(face)->center()(0)) >= 0.1 &&
                         std::fabs(cell->face(face)->center()(0)) <= 1.0) &&
                        (std::fabs(cell->face(face)->center()(dim - 1)) <
                         1e-12))
                      {
                        cell->face(face)->set_boundary_id(3);
                      }

                    if ((std::fabs(cell->face(face)->center()(0)) >= 0.0 &&
                         std::fabs(cell->face(face)->center()(0)) <= 0.9) &&
                        (std::fabs(1.0 - cell->face(face)->center()(dim - 1)) <
                         1e-12))
                      {
                        cell->face(face)->set_boundary_id(5);
                      }

                    if ((std::fabs(1.0 - cell->face(face)->center()(0)) <
                         1e-12) &&
                        (std::fabs(cell->face(face)->center()(dim - 1)) >=
                           0.0 &&
                         std::fabs(cell->face(face)->center()(dim - 1)) <= 0.9))
                      {
                        cell->face(face)->set_boundary_id(4);
                      }

                    if ((std::fabs(cell->face(face)->center()(0)) < 1e-12) &&
                        (std::fabs(cell->face(face)->center()(dim - 1)) >=
                           0.1 &&
                         std::fabs(cell->face(face)->center()(dim - 1)) <= 1.0))
                      {
                        cell->face(face)->set_boundary_id(6);
                      }
                  }
              }
          }
      }
  }

  FESystem<dim>   fe(FE_Q<dim>(1), dim);
  DoFHandler<dim> dof_handler(triangulation);
  dof_handler.distribute_dofs(fe);

  AffineConstraints<double>    constraints;
  std::set<types::boundary_id> no_normal_flux_boundaries;
  no_normal_flux_boundaries.insert(6);
  VectorTools::compute_no_normal_flux_constraints(dof_handler,
                                                  0,
                                                  no_normal_flux_boundaries,
                                                  constraints);

  constraints.print(deallog.get_file_stream());

  deallog.get_file_stream() << std::flush;
  constraints.close();

  deallog << "OK" << std::endl;
}


int
main()
{
  initlog();
  deallog.get_file_stream().precision(4);
  deallog.get_file_stream().setf(std::ios::fixed);

  run<2>();
}
