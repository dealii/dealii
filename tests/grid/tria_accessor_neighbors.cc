// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2012 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// check the neighbors of cell iterators in 2D
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria.h>

#include <stdlib.h> // for rand()

#include "../tests.h"

void
print_log(const typename Triangulation<2, 2>::cell_iterator c1,
          const typename Triangulation<2, 2>::cell_iterator c2,
          const unsigned int                                face)
{
  std::cout << "Neighbor of " << c1->level() << "." << c1->index()
            << " at face " << face << " should be " << c2->level() << "."
            << c2->index() << std::endl;
  return;
}

void
print_log_either(const typename Triangulation<2, 2>::cell_iterator c1,
                 const typename Triangulation<2, 2>::cell_iterator c2,
                 const typename Triangulation<2, 2>::cell_iterator c3,
                 const unsigned int                                face)
{
  std::cout << "Neighbor of " << c1->level() << "." << c1->index()
            << " at face " << face << " should be either " << c2->level() << "."
            << c3->index() << " or " << c3->level() << "." << c3->index()
            << std::endl;
  return;
}

// Output grid
void
output_grid(const Triangulation<2, 2> *tria)
{
  GridOutFlags::Svg svg_flags;
  svg_flags.label_level_number      = true;
  svg_flags.label_cell_index        = true;
  svg_flags.coloring                = GridOutFlags::Svg::level_number;
  svg_flags.background              = GridOutFlags::Svg::transparent;
  svg_flags.line_thickness          = 1;
  svg_flags.boundary_line_thickness = 2;
  GridOut grid_out;
  grid_out.set_flags(svg_flags);
  grid_out.write_svg(*tria, deallog.get_file_stream());
  return;
}

void
check_neighbors(const typename Triangulation<2, 2>::cell_iterator cell)
{
  // cell corresponds to the refined cell
  for (unsigned int face = 0; face < 4; face++)
    {
      if (cell->face(face)->at_boundary()) // skip boundary faces
        continue;

      // Check if the cells neighbor is refined
      if (cell->neighbor(face)->has_children())
        {
          // Yes:
          for (unsigned int child = 0; child < cell->n_children(); child++)
            {
              // Is this cell anisotropically refined, s.t. the
              // child face matches the face of its parent?
              if (cell->face(face) == cell->child(child)->face(face))
                { // yes
                  if (cell->child(child)->neighbor(face) !=
                      cell->neighbor(face))
                    {
                      print_log(cell->child(child), cell->neighbor(face), face);
                      output_grid(&cell->get_triangulation());
                    }

                  Assert(cell->child(child)->neighbor(face) ==
                           cell->neighbor(face),
                         ExcInternalError()); // Neighbors must coincide
                }

              // If not, neighbor of the children must be either of the children
              // of the neighbor of the parent
              for (unsigned int neighbor_child = 0;
                   neighbor_child < cell->neighbor(face)->n_children();
                   neighbor_child++)
                {
                  // Get matching face
                  if (cell->child(child)->face(face) ==
                      cell->neighbor(face)
                        ->child(neighbor_child)
                        ->face(cell->neighbor_of_neighbor(face)))
                    {
                      // face matches. Is the child refined? This is important
                      // in the following setup
                      //
                      //   + --------- + --- + --- +
                      //   |           |     |     |
                      //   |           |     |     |
                      //   |           |     |     |
                      //   + --------- + --- + --- +
                      //   |           |           |
                      //   |           |           |
                      //   |           |           |
                      //   + --------- + --- + --- +
                      //
                      //   Where the parent of the left two anisotropically
                      //   refined cells is the cell we are checking. The
                      //   neighbor at face 1 then returns the coarsest cell on
                      //   the right. Iterating its children gives the unrefined
                      //   and refined cell. the face matches for the already
                      //   refined face. So this will throw an error. This only
                      //   happens in anisotropic refinement. If the cell were
                      //   to be refined isotropically, as below
                      //
                      //   + --------- + --- + --- +
                      //   |           |     |     |
                      //   |           + --- + --- +
                      //   |           |     |     |
                      //   + --------- + --- + --- +
                      //   |           |           |
                      //   |           |           |
                      //   |           |           |
                      //   + --------- + --- + --- +
                      //
                      //   Then the neighbor should be the parent of the small
                      //   four children.
                      //
                      if (cell->neighbor(face)
                            ->child(neighbor_child)
                            ->has_children() &&
                          cell->neighbor(face)
                              ->child(neighbor_child)
                              ->n_children() < 3)
                        { // yes! Then, it must be either of the children.
                          if (cell->child(child)->neighbor(face) !=
                                cell->neighbor(face)
                                  ->child(neighbor_child)
                                  ->child(0) &&
                              cell->child(child)->neighbor(face) !=
                                cell->neighbor(face)
                                  ->child(neighbor_child)
                                  ->child(1))
                            {
                              print_log_either(cell->child(child),
                                               cell->neighbor(face)
                                                 ->child(neighbor_child)
                                                 ->child(0),
                                               cell->neighbor(face)
                                                 ->child(neighbor_child)
                                                 ->child(1),
                                               face);
                              output_grid(&cell->get_triangulation());
                            }

                          Assert(cell->child(child)->neighbor(face) ==
                                     cell->neighbor(face)
                                       ->child(neighbor_child)
                                       ->child(0) ||
                                   cell->child(child)->neighbor(face) ==
                                     cell->neighbor(face)
                                       ->child(neighbor_child)
                                       ->child(1),
                                 ExcInternalError());
                        }
                      else
                        {
                          if (cell->child(child)->neighbor(face) !=
                              cell->neighbor(face)->child(neighbor_child))
                            {
                              print_log(cell->child(child),
                                        cell->neighbor(face)->child(
                                          neighbor_child),
                                        face);
                              output_grid(&cell->get_triangulation());
                            }

                          Assert(cell->child(child)->neighbor(face) ==
                                   cell->neighbor(face)->child(neighbor_child),
                                 ExcInternalError());
                        }
                    }
                }
            }
        }
      else
        {
          // No:
          for (unsigned int child = 0; child < cell->n_children(); child++)
            {
              if (cell->face(face) == cell->child(child)->face(face))
                {
                  if (cell->child(child)->neighbor(face) !=
                      cell->neighbor(face))
                    {
                      print_log(cell->child(child), cell->neighbor(face), face);
                      output_grid(&cell->get_triangulation());
                    }

                  Assert(cell->child(child)->neighbor(face) ==
                           cell->neighbor(face),
                         ExcInternalError());
                }
            }
        }
    }
}

std::vector<RefinementCase<2>>
get_refinement_cases()
{
  std::vector<RefinementCase<2>> out = {
    RefinementCase<2>::isotropic_refinement,
    RefinementCase<2>::isotropic_refinement,
    RefinementCase<2>::isotropic_refinement,
    RefinementCase<2>::cut_axis(1),
    RefinementCase<2>::isotropic_refinement,
    RefinementCase<2>::isotropic_refinement,
    RefinementCase<2>::cut_axis(0),
    RefinementCase<2>::isotropic_refinement,
    RefinementCase<2>::cut_axis(0),
    RefinementCase<2>::isotropic_refinement,
    RefinementCase<2>::cut_axis(1),
    RefinementCase<2>::isotropic_refinement,
    RefinementCase<2>::isotropic_refinement,
    RefinementCase<2>::isotropic_refinement,
    RefinementCase<2>::isotropic_refinement,
    RefinementCase<2>::isotropic_refinement,
    RefinementCase<2>::isotropic_refinement,
    RefinementCase<2>::isotropic_refinement,
    RefinementCase<2>::isotropic_refinement,
    RefinementCase<2>::cut_axis(1),
    RefinementCase<2>::isotropic_refinement,
    RefinementCase<2>::cut_axis(1),
    RefinementCase<2>::isotropic_refinement,
    RefinementCase<2>::isotropic_refinement,
    RefinementCase<2>::cut_axis(0),
    RefinementCase<2>::isotropic_refinement,
    RefinementCase<2>::cut_axis(1),
    RefinementCase<2>::cut_axis(1),
    RefinementCase<2>::isotropic_refinement,
    RefinementCase<2>::cut_axis(1),
    RefinementCase<2>::cut_axis(0),
    RefinementCase<2>::isotropic_refinement,
    RefinementCase<2>::cut_axis(1),
    RefinementCase<2>::cut_axis(1),
    RefinementCase<2>::isotropic_refinement,
    RefinementCase<2>::isotropic_refinement,
    RefinementCase<2>::isotropic_refinement,
    RefinementCase<2>::isotropic_refinement,
    RefinementCase<2>::isotropic_refinement,
    RefinementCase<2>::isotropic_refinement,
    RefinementCase<2>::isotropic_refinement,
    RefinementCase<2>::cut_axis(1),
    RefinementCase<2>::isotropic_refinement,
    RefinementCase<2>::isotropic_refinement,
    RefinementCase<2>::isotropic_refinement,
    RefinementCase<2>::isotropic_refinement,
    RefinementCase<2>::cut_axis(1),
    RefinementCase<2>::isotropic_refinement,
    RefinementCase<2>::cut_axis(1),
    RefinementCase<2>::isotropic_refinement,
    RefinementCase<2>::cut_axis(0),
    RefinementCase<2>::isotropic_refinement,
    RefinementCase<2>::isotropic_refinement,
    RefinementCase<2>::isotropic_refinement,
    RefinementCase<2>::isotropic_refinement,
    RefinementCase<2>::isotropic_refinement,
    RefinementCase<2>::cut_axis(1),
    RefinementCase<2>::cut_axis(1),
    RefinementCase<2>::cut_axis(0),
    RefinementCase<2>::isotropic_refinement,
    RefinementCase<2>::cut_axis(0),
    RefinementCase<2>::cut_axis(1),
    RefinementCase<2>::isotropic_refinement,
    RefinementCase<2>::cut_axis(1),
    RefinementCase<2>::isotropic_refinement,
    RefinementCase<2>::isotropic_refinement,
    RefinementCase<2>::isotropic_refinement,
    RefinementCase<2>::cut_axis(1),
    RefinementCase<2>::cut_axis(1),
    RefinementCase<2>::isotropic_refinement,
    RefinementCase<2>::isotropic_refinement,
    RefinementCase<2>::isotropic_refinement,
    RefinementCase<2>::isotropic_refinement,
    RefinementCase<2>::cut_axis(0),
    RefinementCase<2>::isotropic_refinement,
    RefinementCase<2>::cut_axis(0),
    RefinementCase<2>::cut_axis(1),
    RefinementCase<2>::isotropic_refinement,
    RefinementCase<2>::cut_axis(1),
    RefinementCase<2>::cut_axis(1),
    RefinementCase<2>::isotropic_refinement,
    RefinementCase<2>::isotropic_refinement,
    RefinementCase<2>::cut_axis(1),
    RefinementCase<2>::isotropic_refinement,
    RefinementCase<2>::isotropic_refinement,
    RefinementCase<2>::isotropic_refinement,
    RefinementCase<2>::isotropic_refinement,
    RefinementCase<2>::isotropic_refinement,
    RefinementCase<2>::isotropic_refinement,
    RefinementCase<2>::isotropic_refinement,
    RefinementCase<2>::isotropic_refinement,
    RefinementCase<2>::isotropic_refinement,
    RefinementCase<2>::cut_axis(0),
    RefinementCase<2>::isotropic_refinement,
    RefinementCase<2>::isotropic_refinement,
    RefinementCase<2>::isotropic_refinement,
    RefinementCase<2>::cut_axis(0),
    RefinementCase<2>::isotropic_refinement,
    RefinementCase<2>::isotropic_refinement,
    RefinementCase<2>::isotropic_refinement};
  return out;
}

std::vector<unsigned int>
get_cell_numbers()
{
  std::vector<unsigned int> out = {
    2,   5,   5,   12,  2,   16,  2,   4,   18,  15,  8,   21,  13,  19,  11,
    15,  10,  9,   56,  50,  63,  23,  50,  89,  82,  96,  38,  105, 4,   73,
    38,  81,  89,  29,  143, 65,  166, 135, 12,  130, 46,  213, 132, 71,  106,
    102, 0,   226, 258, 175, 278, 196, 57,  272, 288, 260, 218, 33,  337, 86,
    255, 267, 163, 179, 223, 175, 291, 347, 260, 272, 410, 369, 4,   90,  170,
    417, 361, 275, 413, 287, 218, 362, 375, 476, 69,  24,  138, 195, 149, 32,
    219, 252, 112, 127, 583, 415, 583, 522, 288, 125};
  return out;
}


void
test()
{
  // Get a grid in 2D
  Triangulation<2> triangulation;
  GridGenerator::hyper_cube(triangulation);
  triangulation.refine_global();

  const std::vector<RefinementCase<2>> refinements  = get_refinement_cases();
  const std::vector<unsigned int>      cell_numbers = get_cell_numbers();

  // Generate some random refinement
  // The idea is as follows
  // ref_type is a random number between 0 and 2.
  //  -> If ref_type is 0, 1 we refine that axis
  //  -> If ref_type is 2 we use isotropic refinement
  // cell_number refers to the cell we want to refine
  // obtained by advancing triangulation.begin_active()
  // cell_number times
  for (unsigned int r = 0; r < 100; r++)
    {
      auto cell = triangulation.begin_active();
      std::advance(cell, cell_numbers[r]);

      cell->set_refine_flag(refinements[r]);

      triangulation.execute_coarsening_and_refinement();
      check_neighbors(cell);
    }
  output_grid(&triangulation);
}


int
main()
{
  initlog();

  test();

  return 0;
}
