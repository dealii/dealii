/* $Id$ */
/* Author: Wolfgang Bangerth, University of Heidelberg, 1999 */

				 // The following includes are just
				 // like for the previous program, so
				 // will not be commented further
#include <grid/tria.h>
#include <grid/tria_accessor.h>
#include <grid/tria_iterator.h>
#include <grid/grid_generator.h>
#include <grid/tria_boundary_lib.h>

				 // We need this include file for the
				 // association of degrees of freedom
				 // ("DoF"s) to vertices, lines, and
				 // cells.
#include <dofs/dof_handler.h>
				 // The following include contains the
				 // description of the bilinear finite
				 // element, including the facts that
				 // it has one degree of freedom on
				 // each vertex of the triangulation,
				 // but none on faces and none in the
				 // interior of the cells.
				 //
				 // In fact, the file contains the
				 // description of several more finite
				 // elements as well, such as
				 // biquadratic, bicubic and biquartic
				 // elements, but not only for two
				 // space dimensions, but also for one
				 // and three dimensions.
#include <fe/fe_lib.lagrange.h>
				 // In the following file, several
				 // tools for manipulating degrees of
				 // freedom can be found:
#include <dofs/dof_tools.h>
				 // We will use a sparse matrix to
				 // visualize the pattern of nonzero
				 // entries resulting from the
				 // distribution of degrees of freedom
				 // on the grid. That class can be
				 // found here:
#include <lac/sparse_matrix.h>
				 // We will want to use a special
				 // algorithm to renumber degrees of
				 // freedom. It is declared here:
#include <numerics/dof_renumbering.h>

				 // This is needed for C++ output:
#include <fstream>



				 // This is the function that produced
				 // the circular grid in the previous
				 // example. The sole difference is
				 // that it returns the grid it
				 // produces via its argument.
				 //
				 // We won't comment on the internals
				 // of this function, since this has
				 // been done in the previous
				 // example. If you don't understand
				 // what is happening here, look
				 // there.
void make_grid (Triangulation<2> &triangulation)
{
  const Point<2> center (1,0);
  const double inner_radius = 0.5,
	       outer_radius = 1.0;
  GridGenerator::hyper_shell (triangulation,
			      center, inner_radius, outer_radius);

				   // This is the single difference to
				   // the respetive function in the
				   // previous program: since we want
				   // to export the triangulation
				   // through this function's
				   // parameter, we need to make sure
				   // that the boundary object lives
				   // at least as long as the
				   // triangulation does. However,
				   // since the boundary object is a
				   // local variable, it would be
				   // deleted at the end of this
				   // function, which is too early; by
				   // declaring it 'static', we can
				   // assure that it lives until the
				   // end of the program.
  static const HyperShellBoundary<2> boundary_description(center);
  triangulation.set_boundary (0, boundary_description);
  
  for (unsigned int step=0; step<5; ++step)
    {
      Triangulation<2>::active_cell_iterator cell, endc;
      cell = triangulation.begin_active();
      endc = triangulation.end();

      for (; cell!=endc; ++cell)
	for (unsigned int vertex=0;
	     vertex < GeometryInfo<2>::vertices_per_cell;
	     ++vertex)
	  {
	    const Point<2> vector_to_center
	      = (cell->vertex(vertex) - center);
	    const double distance_from_center
	      = sqrt(vector_to_center.square());
	    
	    if (fabs(distance_from_center - inner_radius) < 1e-10)
	      {
		cell->set_refine_flag ();
		break;
	      };
	  };

      triangulation.execute_coarsening_and_refinement ();
    };
};


				 // Up to now, we only have a grid,
				 // i.e. some geometrical (the
				 // position of the vertices and which
				 // vertices make up which cell) and
				 // some topological information
				 // (neighborhoods of cells). To use
				 // numerical algorithms, one needs
				 // some logic information in addition
				 // to that: we would like to
				 // associate degree of freedom
				 // numbers to each vertex (or line,
				 // or cell, in case we were using
				 // higher order elements) to later
				 // generate matrices and vectors
				 // which describe a finite element
				 // field on the triangulation.
void distribute_dofs (DoFHandler<2> &dof_handler) 
{
				   // In order to associate degrees of
				   // freedom with features of a
				   // triangulation (vertices, lines,
				   // quadrilaterals), we need an
				   // object which describes how many
				   // degrees of freedom are to be
				   // associated to each of these
				   // objects. For (bi-, tri-)linear
				   // finite elements, this is done
				   // using the FEQ1 class, which
				   // states that one degree of
				   // freedom is to be assigned to
				   // each vertex, while there are
				   // none on lines and inside the
				   // quadrilateral. We first need to
				   // create an object of this class
				   // and use it to distribute the
				   // degrees of freedom. Note that
				   // the DoFHandler object will store
				   // a reference to this object, so
				   // we need to make it static as
				   // well, in order to prevent its
				   // preemptive
				   // destruction. (However, the
				   // library would warn us about this
				   // and exit the program if that
				   // occured. You can check this, if
				   // you want, by removing the
				   // 'static' declaration.)
  static const FEQ1<2> finite_element;
  dof_handler.distribute_dofs (finite_element);

				   // Now we have associated a number
				   // to each vertex, but how can we
				   // visualize this? Unfortunately,
				   // presently there is no way
				   // implemented to directly show the
				   // DoF number associated with each
				   // vertex. However, such
				   // information would hardly ever be
				   // truly important, since the
				   // numbering itself is more or less
				   // arbitrary. There are more
				   // important factors, of which we
				   // will visualize one in the
				   // following.
				   //
				   // Associated with each vertex of
				   // the triangulation is a shape
				   // function. Assume we want to
				   // solve something like Laplace's
				   // equation, then the different
				   // matrix entries will be the
				   // integrals over the gradient of
				   // each two such shape
				   // functions. Obviously, since the
				   // shape functions are not equal to
				   // zero only on the cells adjacent
				   // to the vertex they are
				   // associated to, matrix entries
				   // will be nonzero only of the
				   // supports of the shape functions
				   // associated to the column and row
				   // numbers intersect. This is only
				   // the case for adjacent shape
				   // functions, and therefore only
				   // for adjacent vertices. Now,
				   // since the vertices are numbered
				   // more or less randomly be the
				   // above function
				   // (distribute_dofs), the pattern
				   // of nonzero entries in the matrix
				   // will be somewhat ragged, and we
				   // will take a look at it now.
				   //
				   // First we have to create a
				   // structure which we use to store
				   // the places of nonzero
				   // elements. We have to give it the
				   // size of the matrix, which in our
				   // case will be square with that
				   // many rows and columns as there
				   // are degrees of freedom on the
				   // grid:
  SparsityPattern sparsity_pattern (dof_handler.n_dofs(),
				    dof_handler.n_dofs());
				   // We fill it with the places where
				   // nonzero elements will be located
				   // given the present numbering of
				   // degrees of freedom:
  DoFTools::make_sparsity_pattern (dof_handler, sparsity_pattern);
				   // Before further work can be done
				   // on the object, we have to allow
				   // for some internal
				   // reorganization:
  sparsity_pattern.compress ();

				   // Now write the results to a file
  ofstream out ("sparsity_pattern.1");
  sparsity_pattern.print_gnuplot (out);
				   // The result is in GNUPLOT format,
				   // where in each line of the output
				   // file, the coordinates of one
				   // nonzero entry are listed. The
				   // output will be shown below.
				   //
				   // If you look at it, you will note
				   // that the sparsity pattern is
				   // symmetric, which is quite often
				   // so, unless you have a rather
				   // special equation you want to
				   // solve. You will also note that
				   // it has several distinct region,
				   // which stem from the fact that
				   // the numbering starts from the
				   // coarsest cells and moves on to
				   // the finer ones; since they are
				   // all distributed symmetrically
				   // around the origin, this shows up
				   // again in the sparsity pattern.
};



				 // In the sparsity pattern produced
				 // above, the nonzero entries
				 // extended quite far off from the
				 // diagonal. For some algorithms,
				 // this is unfavorable, and we will
				 // show a simple way how to improve
				 // this situation.
				 //
				 // Remember that for an entry (i,j)
				 // in the matrix to be nonzero, the
				 // supports of the shape functions i
				 // and j needed to intersect
				 // (otherwise in the integral, the
				 // integrand would be zero everywhere
				 // since either the one or the other
				 // shape function is zero at some
				 // point). However, the supports of
				 // shape functions intersected only
				 // of they were adjacent to each
				 // other, so in order to have the
				 // nonzero entries clustered around
				 // the diagonal (where i equals j),
				 // we would like to have adjacent
				 // shape functions to be numbered
				 // with indices (DoF numbers) that
				 // differ not too much.
				 //
				 // This can be accomplished by a
				 // simple front marching algorithm,
				 // where one starts at a given vertex
				 // and gives it the index zero. Then,
				 // its neighbors are numbered
				 // successively, making their indices
				 // close to the original one. Then,
				 // their neighbors, if not yet
				 // numbered, are numbered, and so
				 // on. One such algorithm is the one
				 // by Cuthill and McKee, which is a
				 // little more complicated, but works
				 // along the same lines. We will use
				 // it to renumber the degrees of
				 // freedom such that the resulting
				 // sparsity pattern is more localized
				 // around the diagonal.
void renumber_dofs (DoFHandler<2> &dof_handler) 
{
				   // Renumber the degrees of freedom...
  DoFRenumbering::Cuthill_McKee (dof_handler);
				   // ...regenerate the sparsity pattern...
  SparsityPattern sparsity_pattern (dof_handler.n_dofs(),
				    dof_handler.n_dofs());
  DoFTools::make_sparsity_pattern (dof_handler, sparsity_pattern);
  sparsity_pattern.compress ();
				   // ...and output the result:
  ofstream out ("sparsity_pattern.2");
  sparsity_pattern.print_gnuplot (out);
				   // Again, the output is shown
				   // below. Note that the nonzero
				   // entries are clustered far better
				   // around the diagonal than
				   // before. This effect is even more
				   // distinguished for larger
				   // matrices (the present one has
				   // 1260 rows and columns, but large
				   // matrices often have several
				   // 100,000s).
};




				 // This is the main program, which
				 // only calls the other functions in
				 // their respective order.
int main () 
{
				   // Allocate space for a triangulation...
  Triangulation<2> triangulation;
				   // ...and create it
  make_grid (triangulation);

				   // A variable that will hold the
				   // information which vertex has
				   // which number. The geometric
				   // information is passed as
				   // parameter and a pointer to the
				   // triangulation will be stored
				   // inside the DoFHandler object.
  DoFHandler<2> dof_handler (triangulation);
				   // Associate vertices and degrees
				   // of freedom.
  distribute_dofs (dof_handler);

				   // Show the effect of renumbering
				   // of degrees of freedom to the
				   // sparsity pattern of the matrix.
  renumber_dofs (dof_handler);
};
