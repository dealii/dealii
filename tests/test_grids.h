// ---------------------------------------------------------------------
//
// Copyright (C) 2006 - 2020 by the deal.II authors
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

#include <deal.II/base/point.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <vector>

#include "tests.h"

/**
 * A set of test meshes for the deal.II test suite.
 *
 * These meshes exhibit certain key features for writing tests. If you
 * want to test certain properties of algorithms, the following table
 * might be of help.
 *
 * <table border=1>
 * <tr><th>Mesh</th><th>Feature tested</th></tr>
 * <tr><td>#hypercube(tr)</td><td>works at all on a single
 * cell</td></tr>
 * <tr><td>#hypercube(tr,2)</td><td>works on uniform meshes</td></tr>
 * <tr><td>#hypercube(tr,3,true)</td><td>works with local
 * refinement</td></tr>
 * <tr><td>#star_shaped(tr,1)</td><td>method is robust if more than
 * usual cells meet in one vertex</td></tr>
 * <tr><td>#star_shaped(tr,2,true)</td><td>method is robust if more than
 * usual cells meet in one vertex and local refinement exceeds one
 * level</td></tr>
 * <tr><td>#hyper_line(tr,n)</td><td>aligns n unit cells in
 * x-direction</td></tr>
 * </table>
 */
namespace TestGrids
{
  /**
   * Generate grids based on hypercube. These meshes have a regular geometry and
   * topology.
   *
   * @param <tt>refinement</tt> denotes the number of refinement steps of the
   * root cell.
   *
   * @param if <tt>local</tt> is <tt>true</tt>, refine only the cell containing
   * the corner with only negative coordinates.
   */
  template <int dim>
  void
  hypercube(Triangulation<dim> &tr,
            unsigned int        refinement = 0,
            bool                local      = false)
  {
    GridGenerator::hyper_cube(tr, -1., 1.);
    if (refinement && !local)
      tr.refine_global(refinement);
    if (refinement && local)
      {
        tr.refine_global(1);
        for (unsigned int i = 1; i < refinement; ++i)
          {
            for (typename Triangulation<dim>::active_cell_iterator cell =
                   tr.begin_active();
                 cell != tr.end();
                 ++cell)
              {
                const Point<dim> &p        = cell->center();
                bool              negative = true;
                for (unsigned int d = 0; d < dim; ++d)
                  if (p(d) >= 0.)
                    negative = false;
                if (negative)
                  cell->set_refine_flag();
              }
            tr.execute_coarsening_and_refinement();
          }
      }
    deallog << "Triangulation hypercube " << dim << "D refinement "
            << refinement;
    if (local)
      deallog << " local ";
    deallog << " steps " << tr.n_active_cells() << " active cells "
            << tr.n_cells() << " total cells " << std::endl;
  }

  /**
   * Create a star-shaped mesh, having more than the average
   * <tt>2<sup>dim</sup></tt> cells in the central vertex.
   *
   * @param <tt>refinement</tt> denotes the number of refinement
   * steps of the root mesh.
   *
   * @param if <tt>local</tt> is <tt>true</tt>, refine only one of
   * the coarse cells.
   */
  template <int dim>
  void
  star_shaped(Triangulation<dim> &tr,
              unsigned int        refinement = 0,
              bool                local      = false);

  /**
   * Local refinement of every other cell in a checkerboard fashion.
   */
  template <int dim>
  void
  checkers(Triangulation<dim> &tr);

  /**
   * Islands of local refinement.
   */
  template <int dim>
  void
  islands(Triangulation<dim> &tr);

  /**
   * Local refinement with an unrefined hole.
   */
  template <int dim>
  void
  laguna(Triangulation<dim> &tr);

  /**
   * Generates grid with @p n_cells unit cells $[0,1]^d$ aligned in the
   * x-coordinate direction.
   *
   * The domain this grid represents covers $[0,n_cells] \times [0,1]^{d-1}$.
   * Thus, integer division of the x-cordinate $x // 1$ on any point in this
   * triangulation gives you the index of the cell the point is located in.
   */
  template <int dim, int spacedim>
  void
  hyper_line(Triangulation<dim, spacedim> &tr, const int n_cells)
  {
    std::vector<unsigned int> repetitions(dim, 1);
    repetitions[0] = n_cells;
    Point<dim> p1, p2;
    for (unsigned int d = 0; d < dim; ++d)
      {
        p1[d] = 0;
        p2[d] = (d == 0) ? n_cells : 1;
      }
    GridGenerator::subdivided_hyper_rectangle(tr, repetitions, p1, p2);
  }
} // namespace TestGrids
