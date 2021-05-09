/* ---------------------------------------------------------------------
 *
 * Copyright (C) 1999 - 2021 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * The deal.II library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE.md at
 * the top level directory of deal.II.
 *
 * ---------------------------------------------------------------------

 */

// @sect3{Include files}

// The most fundamental class in the library is the Triangulation class, which
// is declared here:
#include <deal.II/grid/tria.h>
// Here are some functions to generate standard grids:
#include <deal.II/grid/grid_generator.h>
// Output of grids in various graphics formats:
#include <deal.II/grid/grid_out.h>

// This is needed for C++ output:
#include <iostream>
#include <fstream>
// And this for the declarations of the `std::sqrt` and `std::fabs` functions:
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
void first_grid()
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
  GridGenerator::hyper_cube(triangulation);
  triangulation.refine_global(4);

  // Now we want to write a graphical representation of the mesh to an output
  // file. The GridOut class of deal.II can do that in a number of different
  // output formats; here, we choose scalable vector graphics (SVG) format
  // that you can visualize using the web browser of your choice:
  std::ofstream out("grid-1.svg");
  GridOut       grid_out;
  grid_out.write_svg(triangulation, out);
  std::cout << "Grid written to grid-1.svg" << std::endl;
}



// @sect3{Creating the second mesh}

// The grid in the following, second function is slightly more complicated in
// that we use a ring domain and refine the result once globally.
void second_grid()
{
  // We start again by defining an object for a triangulation of a
  // two-dimensional domain:
  Triangulation<2> triangulation;

  // We then fill it with a ring domain. The center of the ring shall be the
  // point (1,0), and inner and outer radius shall be 0.5 and 1. The number of
  // circumferential cells could be adjusted automatically by this function,
  // but we choose to set it explicitly to 10 as the last argument:
  const Point<2> center(1, 0);
  const double   inner_radius = 0.5, outer_radius = 1.0;
  GridGenerator::hyper_shell(
    triangulation, center, inner_radius, outer_radius, 10);
  // By default, the triangulation assumes that all boundaries are straight
  // lines, and all cells are bi-linear quads or tri-linear hexes, and that
  // they are defined by the cells of the coarse grid (which we just
  // created). Unless we do something special, when new points need to be
  // introduced the domain is assumed to be delineated by the straight
  // lines of the coarse mesh, and new points will simply be in the middle
  // of the surrounding ones. Here, however, we know that the domain is
  // curved, and we would like to have the Triangulation place new points
  // according to the underlying geometry. Fortunately, some good soul
  // implemented an object which describes a spherical domain, of which the
  // ring is a section; it only needs the center of the ring and
  // automatically figures out how to instruct the Triangulation where to
  // place the new points. The way this works in deal.II is that you tag
  // parts of the triangulation you want to be curved with a number that is
  // usually referred to as "manifold indicator" and then tell the
  // triangulation to use a particular "manifold object" for all places
  // with this manifold indicator. How exactly this works is not important
  // at this point (you can read up on it in step-53 and @ref manifold).
  // The functions in GridGenerator handle this for us in most
  // circumstances: they attach the correct manifold to a domain so that
  // when the triangulation is refined new cells are placed in the correct
  // places. In the present case GridGenerator::hyper_shell attaches a
  // SphericalManifold to all cells: this causes cells to be refined with
  // calculations in spherical coordinates (so new cells have edges that
  // are either radial or lie along concentric circles around the origin).
  //
  // By default (i.e., for a Triangulation created by hand or without a
  // call to a GridGenerator function like GridGenerator::hyper_shell or
  // GridGenerator::hyper_ball), all cells and faces of the Triangulation
  // have their manifold_id set to numbers::flat_manifold_id, which is
  // the default if you want a manifold that produces straight edges, but
  // you can change this number for individual cells and faces. In that
  // case, the curved manifold thus associated with number zero will not
  // apply to those parts with a non-zero manifold indicator, but other
  // manifold description objects can be associated with those non-zero
  // indicators. If no manifold description is associated with a particular
  // manifold indicator, a manifold that produces straight edges is
  // implied. (Manifold indicators are a slightly complicated topic; if
  // you're confused about what exactly is happening here, you may want to
  // look at the @ref GlossManifoldIndicator "glossary entry on this
  // topic".) Since the default chosen by GridGenerator::hyper_shell is
  // reasonable we leave things alone.
  //
  // In order to demonstrate how to write a loop over all cells, we will
  // refine the grid in five steps towards the inner circle of the domain:
  for (unsigned int step = 0; step < 5; ++step)
    {
      // Next, we need to loop over the active cells of the triangulation. You
      // can think of a triangulation as a collection of cells. If it were an
      // array, you would just get a pointer that you increment from one
      // element to the next using the operator `++`. The cells of a
      // triangulation aren't stored as a simple array, but the concept of an
      // <i>iterator</i> generalizes how pointers work to arbitrary collections
      // of objects (see <a href=
      // "http://en.wikipedia.org/wiki/Iterator#C.2B.2B">wikipedia</a> for more
      // information). Typically, any container type in C++ will return an
      // iterator pointing to the start of the collection with a method called
      // `begin`, and an iterator point to 1 past the end of the collection with
      // a method called `end`. We can increment an iterator `it` with the
      // operator `++it`, dereference it to get the underlying data with `*it`,
      // and check to see if we're done by comparing `it != collection.end()`.
      //
      // The second important piece is that we only need the active cells.
      // Active cells are those that are not further refined, and the only
      // ones that can be marked for further refinement. deal.II provides
      // iterator categories that allow us to iterate over <i>all</i> cells
      // (including the parent cells of active ones) or only over the active
      // cells. Because we want the latter, we need to call the method
      // Triangulation::active_cell_iterators().
      //
      // Putting all of this together, we can loop over all the active cells of
      // a triangulation with
      // @code{.cpp}
      //     for (auto it = triangulation.active_cell_iterators().begin();
      //          it != triangulation.active_cell_iterators().end();
      //          ++it)
      //       {
      //         auto cell = *it;
      //         // Then a miracle occurs...
      //       }
      // @endcode
      // In the initializer of this loop, we've used the `auto` keyword for the
      // type of the iterator `it`. The `auto` keyword means that the type of
      // the object being declared will be inferred from the context. This
      // keyword is useful when the actual type names are long or possibly even
      // redundant. If you're unsure of what the type is and want to look up
      // what operations the result supports, you can go to the documentation
      // for the method Triangulation::active_cell_iterators(). In this case,
      // the type of `it` is `Triangulation::active_cell_iterator`.
      //
      // While the `auto` keyword can save us from having to type out long names
      // of data types, we still have to type a lot of redundant declarations
      // about the start and end iterator and how to increment it. Instead of
      // doing that, we'll use
      // <a href="http://en.cppreference.com/w/cpp/language/range-for">range-
      // based for loops</a>, which wrap up all of the syntax shown above into a
      // much shorter form:
      for (auto &cell : triangulation.active_cell_iterators())
        {
          // @note See @ref Iterators for more information about the iterator
          // classes used in deal.II, and @ref CPP11 for more information about
          // range-based for loops and the `auto` keyword.
          //
          // Next, we loop over all vertices of the cells. For that purpose
          // we query an iterator over the vertex indices (in 2d, this is an
          // array that contains the elements `{0,1,2,3}`, but since
          // `cell->vertex_indices()` knows the dimension the cell lives in, the
          // array so returned is correct in all dimensions and this enables
          // this code to be correct whether we run it in 2d or 3d, i.e., it
          // enables "dimension-independent programming" -- a big part of what
          // we will discuss in step-4).
          for (const auto v : cell->vertex_indices())
            {
              // If this cell is at the inner boundary, then at least one of its
              // vertices must sit on the inner ring and therefore have a radial
              // distance from the center of exactly 0.5, up to floating point
              // accuracy. So we compute this distance, and if we find a vertex
              // with this property, we flag this cell for later refinement. We
              // can then also break the loop over all vertices and move on to
              // the next cell.
              //
              // Because the distance from the center is computed as a floating
              // point number, we have to expect that whatever we compute is
              // only accurate to within
              // [round-off](https://en.wikipedia.org/wiki/Round-off_error). As
              // a consequence, we can never expect to compare the distance
              // with the inner radius by equality: A statement such as
              // `if (distance_from_center == inner_radius)` will fail
              // unless we get exceptionally lucky. Rather, we need to do this
              // comparison with a certain tolerance, and the usual way to do
              // this is to write it as `if (std::abs(distance_from_center -
              // inner_radius) <= tolerance)`
              // where `tolerance` is some small number larger
              // than round-off. The question is how to choose it: We could just
              // pick, say, `1e-10`, but this is only appropriate if the objects
              // we compare are of size one. If we had created a mesh with cells
              // of size `1e+10`, then `1e-10` would be far lower than round-off
              // and, as before, the comparison will only succeed if we get
              // exceptionally lucky. Rather, it is almost always useful to make
              // the tolerance *relative* to a typical "scale" of the objects
              // being compared. Here, the "scale" would be the inner radius, or
              // maybe the diameter of cells. We choose the former and set the
              // tolerance equal to $10^{-6}$ times the inner radius of the
              // annulus.
              const double distance_from_center =
                center.distance(cell->vertex(v));

              if (std::fabs(distance_from_center - inner_radius) <=
                  1e-6 * inner_radius)
                {
                  cell->set_refine_flag();
                  break;
                }
            }
        }

      // Now that we have marked all the cells that we want refined, we let
      // the triangulation actually do this refinement. The function that does
      // so owes its long name to the fact that one can also mark cells for
      // coarsening, and the function does coarsening and refinement all at
      // once:
      triangulation.execute_coarsening_and_refinement();
    }


  // Finally, after these five iterations of refinement, we want to again
  // write the resulting mesh to a file, again in SVG format. This works just
  // as above:
  std::ofstream out("grid-2.svg");
  GridOut       grid_out;
  grid_out.write_svg(triangulation, out);

  std::cout << "Grid written to grid-2.svg" << std::endl;
}



// @sect3{The main function}

// Finally, the main function. There isn't much to do here, only to call the
// two subfunctions, which produce the two grids.
int main()
{
  first_grid();
  second_grid();
}
