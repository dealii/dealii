/* ---------------------------------------------------------------------
 *
 * Copyright (C) 1999 - 2014 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * The deal.II library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE at
 * the top level of the deal.II distribution.
 *
 * ---------------------------------------------------------------------

 */

// @sect3{Include files}

// The most fundamental class in the library is the Triangulation class, which
// is declared here:
#include <deal.II/grid/tria.h>
// We need the following two includes for loops over cells and/or faces:
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
// Here are some functions to generate standard grids:
#include <deal.II/grid/grid_generator.h>
// We would like to use boundaries which are not straight lines, so we import
// some classes which predefine some boundary descriptions:
#include <deal.II/grid/tria_boundary_lib.h>
// Output of grids in various graphics formats:
#include <deal.II/grid/grid_out.h>

// This is needed for C++ output:
#include <fstream>
// And this for the declarations of the `sqrt' and `fabs' functions:
#include <cmath>

// The final step in importing deal.II is this: All deal.II functions and
// classes are in a namespace <code>dealii</code>, to make sure they don't
// clash with symbols from other libraries you may want to use in conjunction
// with deal.II. One could use these functions and classes by prefixing every
// use of these names by <code>dealii::</code>, but that would quickly become
// cumbersome and annoying. Rather, we simply import the entire deal.II
// namespace for general use:
using namespace dealii;

// @sect3{Creating the first mesh}

// In the following, first function, we simply use the unit square as domain
// and produce a globally refined grid from it.
void first_grid ()
{
  // The first thing to do is to define an object for a triangulation of a
  // two-dimensional domain:
  Triangulation<2> triangulation;
  // Here and in many following cases, the string "<2>" after a class name
  // indicates that this is an object that shall work in two space
  // dimensions. Likewise, there are versions of the triangulation class that
  // are working in one ("<1>") and three ("<3>") space dimensions. The way
  // this works is through some template magic that we will investigate in
  // some more detail in later example programs; there, we will also see how
  // to write programs in an essentially dimension independent way.

  // Next, we want to fill the triangulation with a single cell for a square
  // domain. The triangulation is the refined four times, to yield $4^4=256$
  // cells in total:
  GridGenerator::hyper_cube (triangulation);
  triangulation.refine_global (4);

  // Now we want to write a graphical representation of the mesh to an output
  // file. The GridOut class of deal.II can do that in a number of different
  // output formats; here, we choose encapsulated postscript (eps) format:
  std::ofstream out ("grid-1.eps");
  GridOut grid_out;
  grid_out.write_eps (triangulation, out);
  std::cout << "Grid written to grid-1.eps" << std::endl;
}



// @sect3{Creating the second mesh}

// The grid in the following, second function is slightly more complicated in
// that we use a ring domain and refine the result once globally.
void second_grid ()
{
  // We start again by defining an object for a triangulation of a
  // two-dimensional domain:
  Triangulation<2> triangulation;

  // We then fill it with a ring domain. The center of the ring shall be the
  // point (1,0), and inner and outer radius shall be 0.5 and 1. The number of
  // circumferential cells could be adjusted automatically by this function,
  // but we choose to set it explicitly to 10 as the last argument:
  const Point<2> center (1,0);
  const double inner_radius = 0.5,
               outer_radius = 1.0;
  GridGenerator::hyper_shell (triangulation,
                              center, inner_radius, outer_radius,
                              10);
  // By default, the triangulation assumes that all boundaries are straight
  // and given by the cells of the coarse grid (which we just created). It
  // uses this information when cells at the boundary are refined and new
  // points need to be introduced on the boundary; if the boundary is assumed
  // to be straight, then new points will simply be in the middle of the
  // surrounding ones.
  //
  // Here, however, we would like to have a curved boundary. Fortunately, some
  // good soul implemented an object which describes the boundary of a ring
  // domain; it only needs the center of the ring and automatically figures
  // out the inner and outer radius when needed. Note that we associate this
  // boundary object with that part of the boundary that has the "boundary
  // indicator" zero. By default (at least in 2d and 3d, the 1d case is
  // slightly different), all boundary parts have this number, but you can
  // change this number for some parts of the boundary. In that case, the
  // curved boundary thus associated with number zero will not apply on those
  // parts with a non-zero boundary indicator, but other boundary description
  // objects can be associated with those non-zero indicators. If no boundary
  // description is associated with a particular boundary indicator, a
  // straight boundary is implied. (Boundary indicators are a slightly
  // complicated topic; if you're confused about what exactly is happening
  // here, you may want to look at the
  // @ref GlossBoundaryIndicator "glossary entry on this topic".)
  const HyperShellBoundary<2> boundary_description(center);
  triangulation.set_boundary (0, boundary_description);

  // In order to demonstrate how to write a loop over all cells, we will
  // refine the grid in five steps towards the inner circle of the domain:
  for (unsigned int step=0; step<5; ++step)
    {
      // Next, we need an iterator that points to a cell and which we will
      // move over all active cells one by one. In a sense, you can think of a
      // triangulation as a collection of cells. If it was an array, you would
      // just get a pointer that you move from one to the next. In
      // triangulations, cells aren't stored as an array, so simple pointers
      // do not work, but one can generalize pointers to iterators (see <a
      // href="http://en.wikipedia.org/wiki/Iterator#C.2B.2B">this wikipedia
      // link</a> for more information). We will then get an iterator to the
      // first cell and iterate over all of the cells until we hit the last
      // one.
      //
      // The second important piece is that we only need the active cells.
      // Active cells are those that are not further refined, and the only
      // ones that can be marked for further refinement, obviously. deal.II
      // provides iterator categories that allow us to iterate over <i>all</i>
      // cells (including the parent cells of active ones) or only over the
      // active cells. Because we want the latter, we need to choose
      // Triangulation::active_cell_iterator as data type.
      //
      // Finally, by convention, we almost always use the names
      // <code>cell</code> and <code>endc</code> for the iterator pointing to
      // the present cell and to the "one-past-the-end" iterator. This is, in
      // a sense a misnomer, because the object is not really a "cell": it is
      // an iterator/pointer to a cell. We should really have started to call
      // these objects <code>cell_iterator</code> when deal.II started in
      // 1998, but it is what it is.
      //
      // After declaring the iterator variable, the loop over all cells is
      // then rather trivial, and looks like any loop involving pointers
      // instead of iterators:
      Triangulation<2>::active_cell_iterator
      cell = triangulation.begin_active(),
      endc = triangulation.end();
      for (; cell!=endc; ++cell)
        {
          // @note Writing a loop like this requires a lot of typing, but it
          // is the only way of doing it in C++98 and C++03. However, if you
          // have a C++11-compliant compiler, you can also use the C++11
          // range-based for loop style that requires significantly less
          // typing. Take a look at @ref CPP11 "the deal.II C++11 page" to see
          // how this works.
          //
          // Next, we want to loop over all vertices of the cells. Since we are
          // in 2d, we know that each cell has exactly four vertices. However,
          // instead of penning down a 4 in the loop bound, we make a first
          // attempt at writing it in a dimension-independent way by which we
          // find out about the number of vertices of a cell. Using the
          // GeometryInfo class, we will later have an easier time getting the
          // program to also run in 3d: we only have to change all occurrences
          // of <code>&lt;2&gt;</code> to <code>&lt;3&gt;</code>, and do not
          // have to audit our code for the hidden appearance of magic numbers
          // like a 4 that needs to be replaced by an 8:
          for (unsigned int v=0;
               v < GeometryInfo<2>::vertices_per_cell;
               ++v)
            {
              // If this cell is at the inner boundary, then at least one of its
              // vertices must sit on the inner ring and therefore have a radial
              // distance from the center of exactly 0.5, up to floating point
              // accuracy. Compute this distance, and if we have found a vertex
              // with this property flag this cell for later refinement. We can
              // then also break the loop over all vertices and move on to the
              // next cell.
              const double distance_from_center
                = center.distance (cell->vertex(v));

              if (std::fabs(distance_from_center - inner_radius) < 1e-10)
                {
                  cell->set_refine_flag ();
                  break;
                }
            }
        }

      // Now that we have marked all the cells that we want refined, we let
      // the triangulation actually do this refinement. The function that does
      // so owes its long name to the fact that one can also mark cells for
      // coarsening, and the function does coarsening and refinement all at
      // once:
      triangulation.execute_coarsening_and_refinement ();
    }


  // Finally, after these five iterations of refinement, we want to again
  // write the resulting mesh to a file, again in eps format. This works just
  // as above:
  std::ofstream out ("grid-2.eps");
  GridOut grid_out;
  grid_out.write_eps (triangulation, out);

  std::cout << "Grid written to grid-2.eps" << std::endl;

  // At this point, all objects created in this function will be destroyed in
  // reverse order. Unfortunately, we defined the boundary object after the
  // triangulation, which still has a pointer to it and the library will
  // produce an error if the boundary object is destroyed before the
  // triangulation. We therefore have to release it, which can be done as
  // follows. Note that this sets the boundary object used for part "0" of the
  // boundary back to a default object, over which the triangulation has full
  // control.
  triangulation.set_boundary (0);
  // An alternative to doing so, and one that is frequently more convenient,
  // would have been to declare the boundary object before the triangulation
  // object. In that case, the triangulation would have let lose of the
  // boundary object upon its destruction, and everything would have been
  // fine.
}



// @sect3{The main function}

// Finally, the main function. There isn't much to do here, only to call the
// two subfunctions, which produce the two grids.
int main ()
{
  first_grid ();
  second_grid ();
}
