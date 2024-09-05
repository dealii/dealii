// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2001 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"


Point<3>
trans_func(const Point<3> &p)
{
  Point<3> r(p[0] + p[1] * p[1], p[1], p[2]);
  return r;
}


void
check_transform(Triangulation<3, 3> *tria, Triangulation<3, 3> *tria_ref)
{
  // run the test
  std::vector<bool> treated_vertices(tria_ref->n_vertices(), false);
  auto              active_cell  = tria->begin_active();
  auto              tria_ref_end = tria_ref->end();
  for (auto active_cell_ref = tria_ref->begin_active();
       active_cell_ref != tria_ref_end;
       ++active_cell_ref, ++active_cell)
    {
      // First, check hanging nodes
      for (unsigned int face = 0; face < 6; ++face)
        if (active_cell_ref->face(face)->has_children() &&
            !active_cell_ref->face(face)->at_boundary())
          {
            for (unsigned int line = 0; line < 4; ++line)
              if (active_cell_ref->face(face)->line(line)->has_children())
                {
                  auto transformed_vertex =
                    (trans_func(
                       active_cell_ref->face(face)->line(line)->vertex(0)) +
                     trans_func(
                       active_cell_ref->face(face)->line(line)->vertex(1))) /
                    2.0;
                  treated_vertices[active_cell->face(face)
                                     ->line(line)
                                     ->child(0)
                                     ->vertex_index(1)] = true;
                  if (active_cell->face(face)->line(line)->child(0)->vertex(
                        1) != transformed_vertex)
                    {
                      std::cout << active_cell->face(face)
                                     ->line(line)
                                     ->child(0)
                                     ->vertex(1)
                                << " should be " << transformed_vertex
                                << std::endl;
                      Assert(false, ExcInternalError());
                    }
                }

            if (static_cast<uint8_t>(
                  active_cell_ref->face(face)->refinement_case()) ==
                RefinementCase<2>::isotropic_refinement)
              {
                // The middle of the anisotropic face is adjusted differently
                auto transformed_vertex =
                  (trans_func(active_cell_ref->face(face)->vertex(0)) +
                   trans_func(active_cell_ref->face(face)->vertex(1)) +
                   trans_func(active_cell_ref->face(face)->vertex(2)) +
                   trans_func(active_cell_ref->face(face)->vertex(3))) /
                  4.0;
                treated_vertices
                  [active_cell->face(face)->child(0)->vertex_index(3)] = true;
                if (active_cell->face(face)->child(0)->vertex(3) !=
                    transformed_vertex)
                  {
                    std::cout
                      << active_cell->face(face)->child(0)->vertex_index(3)
                      << " should be " << transformed_vertex << std::endl;
                    Assert(false, ExcInternalError());
                  }
              }
          }

      // Then, check all the other vertices
      // std::cout << "Testing remaining vertices on " <<
      // active_cell_ref->level()
      //          << "." << active_cell_ref->index() << std::endl;
      for (unsigned int vertex = 0; vertex < 8; vertex++)
        if (!treated_vertices[active_cell->vertex_index(vertex)])
          {
            auto transformed_vertex =
              trans_func(active_cell_ref->vertex(vertex));
            treated_vertices[active_cell->vertex_index(vertex)] = true;
            if (active_cell->vertex(vertex) != transformed_vertex)
              {
                std::cout << active_cell->vertex(vertex) << " should be "
                          << transformed_vertex << std::endl;
                Assert(false, ExcInternalError());
              }
          }
    }
}

void
refine_grid(const bool anisotropic, Triangulation<3, 3> *tria)
{
  tria->refine_global(1);
  tria->execute_coarsening_and_refinement();

  const unsigned int n_refs = anisotropic ? 3 : 2;
  auto               cell   = tria->begin_active();
  for (unsigned int i = 0; i < n_refs; i++)
    {
      if (anisotropic)
        cell->set_refine_flag(RefinementCase<3>::cut_axis(i));
      else
        {
          cell->set_refine_flag();
          cell++;
        }
      cell++;
    }
  tria->execute_coarsening_and_refinement();
}

void
test(const bool anisotropic)
{
  Triangulation<3> tria;
  {
    GridGenerator::hyper_cube(tria);
    refine_grid(anisotropic, &tria);

    if (anisotropic)
      deallog << "Anisotropic test" << std::endl;
    else
      deallog << "Isotropic test" << std::endl;

    deallog << "Unchanged grid:" << std::endl;
    GridOut().write_gnuplot(tria, deallog.get_file_stream());
    // {
    //   std::ofstream f("grid1.vtk");
    //   GridOut().write_vtk(tria, f);
    // }

    GridTools::transform(trans_func, tria);
    deallog << "transformed grid:" << std::endl;
    GridOut().write_gnuplot(tria, deallog.get_file_stream());
    // {
    //   std::ofstream f("grid2.vtk");
    //   GridOut().write_vtk(tria, f);
    // }
  }

  Triangulation<3> tria_ref;
  GridGenerator::hyper_cube(tria_ref);
  refine_grid(anisotropic, &tria_ref);

  check_transform(&tria, &tria_ref);
}


int
main()
{
  initlog();
  test(/* isotropic   */ true);
  test(/* anisotropic */ false);

  return 0;
}
