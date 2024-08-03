/* ------------------------------------------------------------------------
 *
 * SPDX-License-Identifier: LGPL-2.1-or-later
 * Copyright (C) 1999 - 2024 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * Part of the source code is dual licensed under Apache-2.0 WITH
 * LLVM-exception OR LGPL-2.1-or-later. Detailed license information
 * governing the source code and code contributions can be found in
 * LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
 *
 * ------------------------------------------------------------------------
 */


// The first few includes are just like in the previous program, so do not
// require additional comments:
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>

// However, the next file is new. We need this include file for the
// association of degrees of freedom ("DoF"s) to vertices, lines, and cells:
#include <deal.II/dofs/dof_handler.h>

// The following include contains the description of the bilinear finite
// element, including the facts that it has one degree of freedom on each
// vertex of the triangulation, but none on faces and none in the interior of
// the cells.
//
// (In fact, the file contains the description of Lagrange elements in
// general, i.e. also the quadratic, cubic, etc versions, and not only for 2d
// but also 1d and 3d.)
#include <deal.II/fe/fe_q.h>
// In the following file, several tools for manipulating degrees of freedom
// can be found, and the one after it is necessary to call one of the
// functions imported from `dof_tools.h`:
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/mapping_q1.h>

// We will use a sparse matrix to visualize the pattern of nonzero entries
// resulting from the distribution of degrees of freedom on the grid. That
// class can be found here:
#include <deal.II/lac/sparse_matrix.h>
// We will also need to use an intermediate sparsity pattern structure, which
// is found in this file:
#include <deal.II/lac/dynamic_sparsity_pattern.h>

// We will want to use a special algorithm to renumber degrees of freedom. It
// is declared here:
#include <deal.II/dofs/dof_renumbering.h>

// And this is again needed for C++ output:
#include <fstream>

// Finally, as in step-1, we import the deal.II namespace into the global
// scope:
using namespace dealii;


// @sect3{Mesh generation}

// This is the function that produced the circular grid in the previous step-1
// example program with fewer refinements steps. The sole difference is that it
// returns the grid it produces via its argument.
//
// At the end of the function, we also output this mesh into a
// file. We will use this as one piece of information when visualizing
// the location of degrees of freedom. To output a mesh, we use the
// GridOut class that you have already seen in step-1; the difference
// is only that we use gnuplot rather than SVG format, because gnuplot
// is the program we will use to visualize DoF locations.
void make_grid(Triangulation<2> &triangulation)
{
  const Point<2> center(1, 0);
  const double   inner_radius = 0.5, outer_radius = 1.0;
  GridGenerator::hyper_shell(
    triangulation, center, inner_radius, outer_radius, 5);

  for (unsigned int step = 0; step < 3; ++step)
    {
      for (const auto &cell : triangulation.active_cell_iterators())
        for (const auto v : cell->vertex_indices())
          {
            const double distance_from_center =
              center.distance(cell->vertex(v));

            if (std::fabs(distance_from_center - inner_radius) <=
                1e-6 * inner_radius)
              {
                cell->set_refine_flag();
                break;
              }
          }

      triangulation.execute_coarsening_and_refinement();
    }

  std::ofstream mesh_file("mesh.gnuplot");
  GridOut().write_gnuplot(triangulation, mesh_file);
}


// @sect3{Outputting the location of degrees of freedom}

// The next function outputs the locations of degrees of freedom for
// later visualization. Where each DoF is located is something the
// DoFHandler object knows, so that is one of the arguments to this
// function. Since we want to do all of this twice (once for the
// original enumeration and once for the renumbered set of degrees of
// freedom), the function also takes as a second argument the name of
// the file into which we want the output to be written.
//
// In order to learn deal.II, it is probably not terribly important to
// understand exactly what this function does, and you can skip over
// it. But if you would like to know anyway: We want to call the
// function DoFTools::map_dofs_to_support_points() that returns a list
// of locations. It does so in the form of a map through which we can
// query (in a statement such as `dof_location_map[42]`) where the DoF
// is located (in the example, where the 42nd DoF is). It puts this
// information into the `dof_location_map` object.
//
// We then use the function
// DoFTools::write_gnuplot_dof_support_point_info() to write this
// information into a file in a format that is understandable to the
// gnuplot program that we will use for visualization in the results
// section.
void write_dof_locations(const DoFHandler<2> &dof_handler,
                         const std::string   &filename)
{
  const std::map<types::global_dof_index, Point<2>> dof_location_map =
    DoFTools::map_dofs_to_support_points(MappingQ1<2>(), dof_handler);

  std::ofstream dof_location_file(filename);
  DoFTools::write_gnuplot_dof_support_point_info(dof_location_file,
                                                 dof_location_map);
}


// @sect3{Creation of a DoFHandler}

// Up to now, we only have a grid, i.e. some geometrical (the position of the
// vertices) and some topological information (how vertices are connected to
// lines, and lines to cells, as well as which cells neighbor which other
// cells). To use numerical algorithms, one needs some logic information in
// addition to that: we would like to associate degree of freedom numbers to
// each vertex (or line, or cell, in case we were using higher order elements)
// to later generate matrices and vectors which describe a finite element
// field on the triangulation.
//
// This function shows how to do this. The object to consider is the
// <code>DoFHandler</code> class template.  Before we do so, however, we first
// need something that describes how many degrees of freedom are to be
// associated to each of these objects. Since this is one aspect of the
// definition of a finite element space, the finite element base class stores
// this information. In the present context, we therefore create an object of
// the derived class <code>FE_Q</code> that describes Lagrange elements. Its
// constructor takes one argument that states the polynomial degree of the
// element, which here is one (indicating a bi-linear element); this then
// corresponds to one degree of freedom for each vertex, while there are none
// on lines and inside the quadrilateral. A value of, say, three given to the
// constructor would instead give us a bi-cubic element with one degree of
// freedom per vertex, two per line, and four inside the cell. In general,
// <code>FE_Q</code> denotes the family of continuous elements with complete
// polynomials (i.e. tensor-product polynomials) up to the specified order.
//
// We first need to create an object of this class and then pass it on to the
// <code>DoFHandler</code> object to allocate storage for the degrees of
// freedom (in deal.II lingo: we <i>distribute degrees of
// freedom</i>).
void distribute_dofs(DoFHandler<2> &dof_handler)
{
  const FE_Q<2> finite_element(1);
  dof_handler.distribute_dofs(finite_element);

  // Now that we have associated a degree of freedom with a global
  // number to each vertex, Let us output this information using the
  // function above:
  write_dof_locations(dof_handler, "dof-locations-1.gnuplot");

  // In practice, we do not often care about where a degree of freedom
  // is geometrically located, and so other than seeing it once via
  // the call above is not practically useful. But where two degrees
  // of freedom are in relation to each other matters in other ways.
  //
  // Associated with each vertex of the triangulation is a shape
  // function. Assume we want to solve something like Laplace's equation, then
  // the different matrix entries will be the integrals over the gradient of
  // each pair of such shape functions. Obviously, since the shape functions
  // are nonzero only on the cells adjacent to the vertex they are associated
  // with, matrix entries will be nonzero only if the supports of the shape
  // functions associated to that column and row %numbers intersect. This is
  // only the case for adjacent shape functions, and therefore only for
  // adjacent vertices. Now, since the vertices are numbered more or less
  // randomly by the above function (DoFHandler::distribute_dofs), the pattern
  // of nonzero entries in the matrix will be somewhat ragged, and we will
  // take a look at it now.
  //
  // First we have to create a structure which we use to store the places of
  // nonzero elements. This can then later be used by one or more sparse
  // matrix objects that store the values of the entries in the locations
  // stored by this sparsity pattern. The class that stores the locations is
  // the SparsityPattern class. As it turns out, however, this class has some
  // drawbacks when we try to fill it right away: its data structures are set
  // up in such a way that we need to have an estimate for the maximal number
  // of entries we may wish to have in each row. In two space dimensions,
  // reasonable values for this estimate are available through the
  // DoFHandler::max_couplings_between_dofs() function, but in three
  // dimensions the function almost always severely overestimates the true
  // number, leading to a lot of wasted memory, sometimes too much for the
  // machine used, even if the unused memory can be released immediately after
  // computing the sparsity pattern. In order to avoid this, we use an
  // intermediate object of type DynamicSparsityPattern that uses a
  // different %internal data structure and that we can later copy into the
  // SparsityPattern object without much overhead. (Some more information on
  // these data structures can be found in the @ref Sparsity topic.) In order
  // to initialize this intermediate data structure, we have to give it the
  // size of the matrix, which in our case will be square with as many rows
  // and columns as there are degrees of freedom on the grid:
  DynamicSparsityPattern dynamic_sparsity_pattern(dof_handler.n_dofs(),
                                                  dof_handler.n_dofs());

  // We then fill this object with the places where nonzero elements will be
  // located given the present numbering of degrees of freedom:
  DoFTools::make_sparsity_pattern(dof_handler, dynamic_sparsity_pattern);

  // Now we are ready to create the actual sparsity pattern that we could
  // later use for our matrix. It will just contain the data already assembled
  // in the DynamicSparsityPattern.
  SparsityPattern sparsity_pattern;
  sparsity_pattern.copy_from(dynamic_sparsity_pattern);

  // With this, we can now write the results to a file:
  std::ofstream out("sparsity-pattern-1.svg");
  sparsity_pattern.print_svg(out);
  // The result is stored in an <code>.svg</code> file, where each nonzero entry
  // in the matrix corresponds with a red square in the image. The output will
  // be shown below.
  //
  // If you look at it, you will note that the sparsity pattern is
  // symmetric. This should not come as a surprise, since we have not given
  // the <code>DoFTools::make_sparsity_pattern</code> any information that
  // would indicate that our bilinear form may couple shape functions in a
  // non-symmetric way. You will also note that it has several distinct
  // region, which stem from the fact that the numbering starts from the
  // coarsest cells and moves on to the finer ones; since they are all
  // distributed symmetrically around the origin, this shows up again in the
  // sparsity pattern.
}


// @sect3{Renumbering of DoFs}

// In the sparsity pattern produced above, the nonzero entries extended quite
// far off from the diagonal. For some algorithms, for example for incomplete
// LU decompositions or Gauss-Seidel preconditioners, this is unfavorable, and
// we will show a simple way how to improve this situation.
//
// Remember that for an entry $(i,j)$ in the matrix to be nonzero, the
// supports of the shape functions i and j needed to intersect (otherwise in
// the integral, the integrand would be zero everywhere since either the one
// or the other shape function is zero at some point). However, the supports
// of shape functions intersected only if they were adjacent to each other, so
// in order to have the nonzero entries clustered around the diagonal (where
// $i$ equals $j$), we would like to have adjacent shape functions to be
// numbered with indices (DoF numbers) that differ not too much.
//
// This can be accomplished by a simple front marching algorithm, where one
// starts at a given vertex and gives it the index zero. Then, its neighbors
// are numbered successively, making their indices close to the original
// one. Then, their neighbors, if not yet numbered, are numbered, and so on.
//
// One algorithm that adds a little bit of sophistication along these lines is
// the one by Cuthill and McKee. We will use it in the following function to
// renumber the degrees of freedom such that the resulting sparsity pattern is
// more localized around the diagonal. The only interesting part of the
// function is the first call to <code>DoFRenumbering::Cuthill_McKee</code>,
// the rest is essentially as before:
void renumber_dofs(DoFHandler<2> &dof_handler)
{
  DoFRenumbering::Cuthill_McKee(dof_handler);

  write_dof_locations(dof_handler, "dof-locations-2.gnuplot");


  DynamicSparsityPattern dynamic_sparsity_pattern(dof_handler.n_dofs(),
                                                  dof_handler.n_dofs());
  DoFTools::make_sparsity_pattern(dof_handler, dynamic_sparsity_pattern);

  SparsityPattern sparsity_pattern;
  sparsity_pattern.copy_from(dynamic_sparsity_pattern);

  std::ofstream out("sparsity-pattern-2.svg");
  sparsity_pattern.print_svg(out);
}

// Again, the output is shown below. Note that the nonzero entries are
// clustered far better around the diagonal than before. This effect is even
// more distinguished for larger matrices (the present one has 1260 rows and
// columns, but large matrices often have several 100,000s).

// It is worth noting that the <code>DoFRenumbering</code> class offers a
// number of other algorithms as well to renumber degrees of freedom. For
// example, it would of course be ideal if all couplings were in the lower or
// upper triangular part of a matrix, since then solving the linear system
// would amount to only forward or backward substitution. This is of course
// unachievable for symmetric sparsity patterns, but in some special
// situations involving transport equations, this is possible by enumerating
// degrees of freedom from the inflow boundary along streamlines to the
// outflow boundary. Not surprisingly, <code>DoFRenumbering</code> also has
// algorithms for this.


// @sect3{The main function}

// Finally, this is the main program. The only thing it does is to allocate
// and create the triangulation, then create a <code>DoFHandler</code> object
// and associate it to the triangulation, and finally call above two functions
// on it:
int main()
{
  Triangulation<2> triangulation;
  make_grid(triangulation);

  DoFHandler<2> dof_handler(triangulation);

  distribute_dofs(dof_handler);
  renumber_dofs(dof_handler);
}
