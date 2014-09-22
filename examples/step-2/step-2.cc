/* ---------------------------------------------------------------------
 *
 * Copyright (C) 1999 - 2013 by the deal.II authors
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

 *
 * Author: Wolfgang Bangerth, University of Heidelberg, 1999
 */


// The first few includes are just like in the previous program, so do not
// require additional comments:
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_boundary_lib.h>

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
// can be found:
#include <deal.II/dofs/dof_tools.h>
// We will use a sparse matrix to visualize the pattern of nonzero entries
// resulting from the distribution of degrees of freedom on the grid. That
// class can be found here:
#include <deal.II/lac/sparse_matrix.h>
// We will also need to use an intermediate sparsity pattern structure, which
// is found in this file:
#include <deal.II/lac/compressed_sparsity_pattern.h>

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
// example program. The sole difference is that it returns the grid it
// produces via its argument.
//
// The details of what the function does are explained in step-1. The only
// thing we would like to comment on is this:
//
// Since we want to export the triangulation through this function's
// parameter, we need to make sure that the boundary object lives at least as
// long as the triangulation does. However, in step-1, the boundary object is
// a local variable, and it would be deleted at the end of the function, which
// is too early. We avoid the problem by declaring it 'static' which makes
// sure that the object is initialized the first time control the program
// passes this point, but at the same time assures that it lives until the end
// of the program.
void make_grid (Triangulation<2> &triangulation)
{
  const Point<2> center (1,0);
  const double inner_radius = 0.5,
               outer_radius = 1.0;
  GridGenerator::hyper_shell (triangulation,
                              center, inner_radius, outer_radius,
                              10);

  static const HyperShellBoundary<2> boundary_description(center);
  triangulation.set_boundary (0, boundary_description);

  for (unsigned int step=0; step<5; ++step)
    {
      Triangulation<2>::active_cell_iterator
      cell = triangulation.begin_active(),
      endc = triangulation.end();

      for (; cell!=endc; ++cell)
        for (unsigned int v=0;
             v < GeometryInfo<2>::vertices_per_cell;
             ++v)
          {
            const double distance_from_center
              = center.distance (cell->vertex(v));

            if (std::fabs(distance_from_center - inner_radius) < 1e-10)
              {
                cell->set_refine_flag ();
                break;
              }
          }

      triangulation.execute_coarsening_and_refinement ();
    }
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
// freedom (in deal.II lingo: we <code>distribute degrees of
// freedom</code>). Note that the DoFHandler object will store a reference to
// this finite element object, so we have to make sure its lifetime is at
// least as long as that of the <code>DoFHandler</code>; one way to make sure
// this is so is to make it static as well, in order to prevent its preemptive
// destruction. (However, the library would warn us if we forgot about this
// and abort the program if that occurred. You can check this, if you want, by
// removing the 'static' declaration.)
void distribute_dofs (DoFHandler<2> &dof_handler)
{
  // As described above, let us first create a finite element object, and then
  // use it to allocate degrees of freedom on the triangulation with which the
  // dof_handler object is associated:
  static const FE_Q<2> finite_element(1);
  dof_handler.distribute_dofs (finite_element);

  // Now that we have associated a degree of freedom with a global number to
  // each vertex, we wonder how to visualize this?  There is no simple way to
  // directly visualize the DoF number associated with each vertex. However,
  // such information would hardly ever be truly important, since the
  // numbering itself is more or less arbitrary. There are more important
  // factors, of which we will demonstrate one in the following.
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
  // intermediate object of type CompressedSparsityPattern that uses a
  // different %internal data structure and that we can later copy into the
  // SparsityPattern object without much overhead. (Some more information on
  // these data structures can be found in the @ref Sparsity module.) In order
  // to initialize this intermediate data structure, we have to give it the
  // size of the matrix, which in our case will be square with as many rows
  // and columns as there are degrees of freedom on the grid:
  CompressedSparsityPattern compressed_sparsity_pattern(dof_handler.n_dofs(),
                                                        dof_handler.n_dofs());

  // We then fill this object with the places where nonzero elements will be
  // located given the present numbering of degrees of freedom:
  DoFTools::make_sparsity_pattern (dof_handler, compressed_sparsity_pattern);

  // Now we are ready to create the actual sparsity pattern that we could
  // later use for our matrix. It will just contain the data already assembled
  // in the CompressedSparsityPattern.
  SparsityPattern sparsity_pattern;
  sparsity_pattern.copy_from (compressed_sparsity_pattern);

  // With this, we can now write the results to a file:
  std::ofstream out ("sparsity_pattern.1");
  sparsity_pattern.print_gnuplot (out);
  // The result is in GNUPLOT format, where in each line of the output file,
  // the coordinates of one nonzero entry are listed. The output will be shown
  // below.
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
void renumber_dofs (DoFHandler<2> &dof_handler)
{
  DoFRenumbering::Cuthill_McKee (dof_handler);

  CompressedSparsityPattern compressed_sparsity_pattern(dof_handler.n_dofs(),
                                                        dof_handler.n_dofs());
  DoFTools::make_sparsity_pattern (dof_handler, compressed_sparsity_pattern);

  SparsityPattern sparsity_pattern;
  sparsity_pattern.copy_from (compressed_sparsity_pattern);

  std::ofstream out ("sparsity_pattern.2");
  sparsity_pattern.print_gnuplot (out);
}

// Again, the output is shown below. Note that the nonzero entries are
// clustered far better around the diagonal than before. This effect is even
// more distinguished for larger matrices (the present one has 1260 rows and
// columns, but large matrices often have several 100,000s).

// It is worth noting that the <code>DoFRenumbering</code> class offers a
// number of other algorithms as well to renumber degrees of freedom. For
// example, it would of course be ideal if all couplings were in the lower or
// upper triangular part of a matrix, since then solving the linear system
// would among to only forward or backward substitution. This is of course
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
int main ()
{
  Triangulation<2> triangulation;
  make_grid (triangulation);

  DoFHandler<2> dof_handler (triangulation);

  distribute_dofs (dof_handler);
  renumber_dofs (dof_handler);
}
