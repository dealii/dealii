// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2019 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Test FESubfaceValues<dim, spacedim>::reinit(
// const typename Triangulation<dim, spacedim>::cell_iterator &,
// const unsigned int, const unsigned int) for periodic faces with hanging
// nodes.

#include <deal.II/base/mpi.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>

#include "../tests.h"

#define PRECISION 4

template <int dim>
void
run(unsigned int degree, unsigned int n_q_points)
{
  Triangulation<dim> tria;

  // Create grid
  GridGenerator::hyper_cube(tria);

  tria.begin_active()->face(0)->set_boundary_id(10);
  tria.begin_active()->face(1)->set_boundary_id(11);
  tria.begin_active()->face(2)->set_boundary_id(12);
  tria.begin_active()->face(3)->set_boundary_id(13);

  // Add periodic boundary conds
  std::vector<
    GridTools::PeriodicFacePair<typename Triangulation<dim>::cell_iterator>>
    periodic_faces;
  GridTools::collect_periodic_faces(tria, 10, 11, 0, periodic_faces);
  GridTools::collect_periodic_faces(tria, 12, 13, 1, periodic_faces);
  tria.add_periodicity(periodic_faces);

  // refine grid
  tria.refine_global(1);
  tria.begin_active()->set_refine_flag();
  tria.execute_coarsening_and_refinement();

  // Output grid file for each processor and global
  GridOut grid_out;
  grid_out.write_mesh_per_processor_as_vtu(tria, "mesh");

  FE_DGQ<dim>            fe(degree);                       // dummy
  QGaussLobatto<dim - 1> quad(n_q_points);                 // dummy
  UpdateFlags            flags = update_quadrature_points; // dummy

  FESubfaceValues<dim> sub_face(fe, quad, flags);
  sub_face.reinit(tria.create_cell_iterator(CellId(0, {1})), 0, 0);
  for (auto p : sub_face.get_quadrature_points())
    deallog << p << std::endl;
  deallog << std::endl;

  sub_face.reinit(tria.create_cell_iterator(CellId(0, {1})), 1, 0);
  for (auto p : sub_face.get_quadrature_points())
    deallog << p << std::endl;
  deallog << std::endl;
}

int
main()
{
  std::ofstream logfile("output");
  deallog << std::setprecision(PRECISION);
  deallog.attach(logfile);

  run<2>(2 /*=degree*/, 3 /*qpoints*/);
}
