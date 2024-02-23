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



// Test a selection of common configurations of pipe junction geometries

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"


/**
 * Refines the provided triangulation globally by the specified amount of times.
 *
 * Also computes the volume of the mesh for debugging purposes. Writes the
 * finest mesh into the log stream.
 */
template <int dim, int spacedim>
void
refine_and_write(Triangulation<dim, spacedim> &tria,
                 const unsigned int            n_global_refinements = 0,
                 const std::string             name                 = "grid")
{
  deallog << "name: " << name << std::endl;

  const auto output = [&](const unsigned int n) {
    deallog << "  n: " << n << ", volume: " << GridTools::volume(tria)
            << std::endl;
  };

  output(0);
  for (unsigned int n = 1; n <= n_global_refinements; ++n)
    {
      tria.refine_global();
      output(n);
    }

  std::ostream &logfile = deallog.get_file_stream();
  GridOut       grid_out;
  grid_out.write_vtk(tria, logfile);
}



/**
 * Tests a selection of common configurations of pipe junction geometries.
 */
void
test_selection()
{
  constexpr unsigned int dim      = 3;
  constexpr unsigned int spacedim = 3;

  // y-pipe in plane
  {
    const std::vector<std::pair<Point<spacedim>, double>> openings = {
      {{{-2., 0., 0.}, 1.},
       {{1., std::sqrt(3), 0.}, 1.},
       {{1., -std::sqrt(3), 0.}, 1.}}};

    const std::pair<Point<spacedim>, double> bifurcation = {{0., 0., 0.}, 1.};

    Triangulation<dim, spacedim> tria;
    GridGenerator::pipe_junction(tria, openings, bifurcation);

    refine_and_write(tria, 1, "ypipe");
  }

  // t-pipe in plane
  {
    const std::vector<std::pair<Point<spacedim>, double>> openings = {
      {{{-2., 0., 0.}, 1.}, {{0., 2., 0.}, 1.}, {{2., 0., 0.}, 1.}}};

    const std::pair<Point<spacedim>, double> bifurcation = {{0., 0., 0.}, 1.};

    Triangulation<dim, spacedim> tria;
    GridGenerator::pipe_junction(tria, openings, bifurcation);

    refine_and_write(tria, 1, "tpipe");
  }

  // corner piece
  {
    const std::vector<std::pair<Point<spacedim>, double>> openings = {
      {{{2., 0., 0.}, 1.}, {{0., 2., 0.}, 1.}, {{0., 0., 2.}, 1.}}};

    const std::pair<Point<spacedim>, double> bifurcation = {{0., 0., 0.}, 1.};

    Triangulation<dim, spacedim> tria;
    GridGenerator::pipe_junction(tria, openings, bifurcation);

    refine_and_write(tria, 1, "corner");
  }

  // irregular configuration with arbitrary points
  {
    const std::vector<std::pair<Point<spacedim>, double>> openings = {
      {{{-4., 0., 0.}, 1.5}, {{4., -8., -0.4}, 0.75}, {{0.1, 0., -6.}, 0.5}}};

    const std::pair<Point<spacedim>, double> bifurcation = {{0., 0., 0.}, 1.};

    Triangulation<dim, spacedim> tria;
    GridGenerator::pipe_junction(tria, openings, bifurcation);

    refine_and_write(tria, 1, "irregular");
  }
}



int
main()
{
  initlog();

  test_selection();
}
