/* $Id$ */
/* Author: Wolfgang Bangerth, University of Heidelberg, 1999 */

/*    $Id$       */
/*                                                                */
/*    Copyright (C) 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007 by the deal.II authors */
/*                                                                */
/*    This file is subject to QPL and may not be  distributed     */
/*    without copyright and license information. Please refer     */
/*    to the file deal.II/doc/license.html for the  text  and     */
/*    further information on this license.                        */

                                 // @sect3{Include files}

				 // The first few (many?) include
				 // files have already been used in
				 // the previous example, so we will
				 // not explain their meaning here
				 // again.
#include <deal.II/grid/tria.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/numerics/vectors.h>
#include <deal.II/numerics/matrices.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>

#include <deal.II/numerics/data_out.h>
#include <fstream>
#include <iostream>

				 // This is new, however: in the previous
				 // example we got some unwanted output from
				 // the linear solvers. If we want to suppress
				 // it, we have to include this file and add a
				 // single line somewhere to the program (see
				 // the main() function below for that):
#include <deal.II/base/logstream.h>


#include <algorithm>
#include <numeric>
#include <deal.II/grid/tria_boundary.h>

				 // The last step is as in all
				 // previous programs:
using namespace dealii;

class PointCloudSurface : public StraightBoundary<3>
{
  public:
				     /**
				      * Constructor.
				      */
    PointCloudSurface (const std::string &filename);
    
				     /**
				      * Let the new point be the
				      * arithmetic mean of the two
				      * vertices of the line.
				      *
				      * Refer to the general
				      * documentation of this class
				      * and the documentation of the
				      * base class for more
				      * information.
				      */
    virtual Point<3>
    get_new_point_on_line (const Triangulation<3>::line_iterator &line) const;

				     /**
				      * Let the new point be the
				      * arithmetic mean of the four
				      * vertices of this quad and the
				      * four midpoints of the lines,
				      * which are already created at
				      * the time of calling this
				      * function.
				      *
				      * Refer to the general
				      * documentation of this class
				      * and the documentation of the
				      * base class for more
				      * information.
				      */
    virtual Point<3>
    get_new_point_on_quad (const Triangulation<3>::quad_iterator &quad) const;

				     /**
				      * Gives <tt>n=points.size()</tt>
				      * points that splits the
				      * StraightBoundary line into
				      * $n+1$ partitions of equal
				      * lengths.
				      *
				      * Refer to the general
				      * documentation of this class
				      * and the documentation of the
				      * base class.
				      */
    virtual void
    get_intermediate_points_on_line (const Triangulation<3>::line_iterator &line,
				     std::vector<Point<3> > &points) const;

				     /**
				      * Gives <tt>n=points.size()=m*m</tt>
				      * points that splits the
				      * p{StraightBoundary} quad into
				      * <tt>(m+1)(m+1)</tt> subquads of equal
				      * size.
				      *
				      * Refer to the general
				      * documentation of this class
				      * and the documentation of the
				      * base class.
				      */
    virtual void
    get_intermediate_points_on_quad (const Triangulation<3>::quad_iterator &quad,
				     std::vector<Point<3> > &points) const;

				     /**
				      * A function that, given a point
				      * <code>p</code>, returns the closest
				      * point on the surface defined by the
				      * input file. For the time being, we
				      * simply return the closest point in the
				      * point cloud, rather than doing any
				      * sort of interpolation.
				      */
    Point<3> closest_point (const Point<3> &p) const;    
  private:
    std::vector<Point<3> > point_list;
};


PointCloudSurface::PointCloudSurface (const std::string &filename)
{
				   // first read in all the points
  {
    std::ifstream in (filename.c_str());
    AssertThrow (in, ExcIO());
    
    while (in)
      {
	Point<3> p;
	in >> p;
	point_list.push_back (p);
      }

    AssertThrow (point_list.size() > 1, ExcIO());
  }
  
				   // next fit a linear model through the data
				   // cloud to rectify it in a local
				   // coordinate system
				   //
				   // the first step is to move the center of
				   // mass of the points to the origin
  {
    const Point<3> c_o_m = std::accumulate (point_list.begin(),
					    point_list.end(),
					    Point<3>()) /
			     point_list.size();
    for (unsigned int i=0; i<point_list.size(); ++i)
      point_list[i] -= c_o_m;
  }

				   // next do a least squares fit to the
				   // function ax+by. this leads to the
				   // following equations:
  
				   // min f(a,b) = sum_i (zi-a xi - b yi)^2 / 2
				   //
				   // f_a = sum_i (zi - a xi - b yi) xi = 0
				   // f_b = sum_i (zi - a xi - b yi) yi = 0
				   //
				   // f_a = (sum_i zi xi) - (sum xi^2) a - (sum xi yi) b = 0
				   // f_a = (sum_i zi yi) - (sum xi yi) a - (sum yi^2) b = 0
  {
    double A[2][2] = {{0,0},{0,0}};
    double B[2] = {0,0};

    for (unsigned int i=0; i<point_list.size(); ++i)
      {
	A[0][0] += point_list[i][0] * point_list[i][0];
	A[0][1] += point_list[i][0] * point_list[i][1];
	A[1][1] += point_list[i][1] * point_list[i][1];

	B[0] += point_list[i][0] * point_list[i][2];
	B[1] += point_list[i][1] * point_list[i][2];
      }

    const double det = A[0][0]*A[1][1]-2*A[0][1];
    const double a = (A[1][1] * B[0] - A[0][1] * B[1]) / det;
    const double b = (A[0][0] * B[1] - A[0][1] * B[0]) / det;


				     // with this information, we can rotate
				     // the points so that the corresponding
				     // least-squares fit would be the x-y
				     // plane
    const Point<2> gradient_direction
      = Point<2>(a,b) / std::sqrt(a*a+b*b);
    const Point<2> orthogonal_direction
      = Point<2>(-b,a) / std::sqrt(a*a+b*b);

    const double stretch_factor = std::sqrt(1.+a*a+b*b);
    
    for (unsigned int i=0; i<point_list.size(); ++i)
      {
					 // we can do that by, for each point,
					 // first subtract the points in the
					 // plane:
	point_list[i][2] -= a*point_list[i][0] + b*point_list[i][1];

				       // we made a mistake here, though:
				       // we've shrunk the plan in the
				       // direction parallel to the
				       // gradient. we will have to correct
				       // for this:
	const Point<2> xy (point_list[i][0],
			   point_list[i][1]);
	const double grad_distance = xy * gradient_direction;
	const double orth_distance = xy * orthogonal_direction;

					 // we then have to stretch the points
					 // in the gradient direction. the
					 // stretch factor is defined above
					 // (zero if the original plane was
					 // already the xy plane, infinity if
					 // it was vertical)
	const Point<2> new_xy
	  = (grad_distance * stretch_factor * gradient_direction +
	     orth_distance * orthogonal_direction);
	point_list[i][0] = new_xy[0];
	point_list[i][1] = new_xy[1];	
      }
  }
}


Point<3>
PointCloudSurface::closest_point (const Point<3> &p) const
{
  double distance = p.distance (point_list[0]);
  Point<3> point = point_list[0];
  
  for (std::vector<Point<3> >::const_iterator i=point_list.begin();
       i != point_list.end(); ++i)
    {
      const double d = p.distance (*i);
      if (d < distance)
	{
	  distance = d;
	  point = *i;
	}
    }

  return point;
}

  
Point<3>
PointCloudSurface::
get_new_point_on_line (const Triangulation<3>::line_iterator &line) const
{
  return closest_point (StraightBoundary<3>::get_new_point_on_line (line));
}



Point<3>
PointCloudSurface::
get_new_point_on_quad (const Triangulation<3>::quad_iterator &quad) const
{
  return closest_point (StraightBoundary<3>::get_new_point_on_quad (quad));
}



void
PointCloudSurface::
get_intermediate_points_on_line (const Triangulation<3>::line_iterator &line,
				 std::vector<Point<3> > &points) const
{
  StraightBoundary<3>::get_intermediate_points_on_line (line,
							points);
  for (unsigned int i=0; i<points.size(); ++i)
    points[i] = closest_point(points[i]);
}



void
PointCloudSurface::
get_intermediate_points_on_quad (const Triangulation<3>::quad_iterator &quad,
				 std::vector<Point<3> > &points) const
{
  StraightBoundary<3>::get_intermediate_points_on_quad (quad,
							points);
  for (unsigned int i=0; i<points.size(); ++i)
    points[i] = closest_point(points[i]);
}



PointCloudSurface pds("surface-points");








                                 // @sect3{The <code>LaplaceProblem</code> class template}

				 // This is again the same
				 // <code>LaplaceProblem</code> class as in the
				 // previous example. The only
				 // difference is that we have now
				 // declared it as a class with a
				 // template parameter, and the
				 // template parameter is of course
				 // the spatial dimension in which we
				 // would like to solve the Laplace
				 // equation. Of course, several of
				 // the member variables depend on
				 // this dimension as well, in
				 // particular the Triangulation
				 // class, which has to represent
				 // quadrilaterals or hexahedra,
				 // respectively. Apart from this,
				 // everything is as before.
template <int dim>
class LaplaceProblem 
{
  public:
    LaplaceProblem ();
    void run ();
    
  private:
    void make_grid_and_dofs ();
    void assemble_system ();
    void solve ();
    void output_results () const;

    Triangulation<dim>   triangulation;
    FE_Q<dim>            fe;
    DoFHandler<dim>      dof_handler;

    SparsityPattern      sparsity_pattern;
    SparseMatrix<double> system_matrix;

    Vector<double>       solution;
    Vector<double>       system_rhs;
};


                                 // @sect3{Right hand side and boundary values}




template <int dim>
class BoundaryValues : public Function<dim> 
{
  public:
    BoundaryValues () : Function<dim>() {}
    
    virtual double value (const Point<dim>   &p,
			  const unsigned int  component = 0) const;
};



template <int dim>
double BoundaryValues<dim>::value (const Point<dim> &p,
				   const unsigned int /*component*/) const 
{
  return std::max(p[dim-1], -5.);
}



                                 // @sect3{Implementation of the <code>LaplaceProblem</code> class}

                                 // Next for the implementation of the class
                                 // template that makes use of the functions
                                 // above. As before, we will write everything
                                 // as templates that have a formal parameter
                                 // <code>dim</code> that we assume unknown at the time
                                 // we define the template functions. Only
                                 // later, the compiler will find a
                                 // declaration of <code>LaplaceProblem@<2@></code> (in
                                 // the <code>main</code> function, actually) and
                                 // compile the entire class with <code>dim</code>
                                 // replaced by 2, a process referred to as
                                 // `instantiation of a template'. When doing
                                 // so, it will also replace instances of
                                 // <code>RightHandSide@<dim@></code> by
                                 // <code>RightHandSide@<2@></code> and instantiate the
                                 // latter class from the class template.
                                 //
                                 // In fact, the compiler will also find a
                                 // declaration <code>LaplaceProblem@<3@></code> in
                                 // <code>main()</code>. This will cause it to again go
                                 // back to the general
                                 // <code>LaplaceProblem@<dim@></code> template, replace
                                 // all occurrences of <code>dim</code>, this time by
                                 // 3, and compile the class a second
                                 // time. Note that the two instantiations
                                 // <code>LaplaceProblem@<2@></code> and
                                 // <code>LaplaceProblem@<3@></code> are completely
                                 // independent classes; their only common
                                 // feature is that they are both instantiated
                                 // from the same general template, but they
                                 // are not convertible into each other, for
                                 // example, and share no code (both
                                 // instantiations are compiled completely
                                 // independently).


                                 // @sect4{LaplaceProblem::LaplaceProblem}

				 // After this introduction, here is the
				 // constructor of the <code>LaplaceProblem</code>
				 // class. It specifies the desired polynomial
				 // degree of the finite elements and
				 // associates the DoFHandler to the
				 // triangulation just as in the previous
				 // example program, step-3:
template <int dim>
LaplaceProblem<dim>::LaplaceProblem () :
                fe (1),
		dof_handler (triangulation)
{}


                                 // @sect4{LaplaceProblem::make_grid_and_dofs}

				 // Grid creation is something
				 // inherently dimension
				 // dependent. However, as long as the
				 // domains are sufficiently similar
				 // in 2D or 3D, the library can
				 // abstract for you. In our case, we
				 // would like to again solve on the
				 // square [-1,1]x[-1,1] in 2D, or on
				 // the cube [-1,1]x[-1,1]x[-1,1] in
				 // 3D; both can be termed
				 // <code>hyper_cube</code>, so we may use the
				 // same function in whatever
				 // dimension we are. Of course, the
				 // functions that create a hypercube
				 // in two and three dimensions are
				 // very much different, but that is
				 // something you need not care
				 // about. Let the library handle the
				 // difficult things.
				 //
				 // Likewise, associating a degree of freedom
				 // with each vertex is something which
				 // certainly looks different in 2D and 3D,
				 // but that does not need to bother you
				 // either. This function therefore looks
				 // exactly like in the previous example,
				 // although it performs actions that in their
				 // details are quite different if <code>dim</code>
				 // happens to be 3. The only significant
				 // difference from a user's perspective is
				 // the number of cells resulting, which is
				 // much higher in three than in two space
				 // dimensions!
template <int dim>
void LaplaceProblem<dim>::make_grid_and_dofs ()
{
  GridGenerator::hyper_cube (triangulation, -30, 30);

  for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
    if (triangulation.begin()->face(f)->center()[2] > 15)
      {
	triangulation.begin()->face(f)->set_boundary_indicator (1);
	for (unsigned int i=0; i<GeometryInfo<dim>::lines_per_face; ++i)
	  triangulation.begin()->face(f)->line(i)->set_boundary_indicator (1);
	break;
      }
  triangulation.set_boundary (1, pds);
  
  
  for (unsigned int v=0; v<GeometryInfo<dim>::vertices_per_cell; ++v)
    if (triangulation.begin()->vertex(v)[2] > 0)
      triangulation.begin()->vertex(v)
	= pds.closest_point (Point<3>(triangulation.begin()->vertex(v)[0],
				      triangulation.begin()->vertex(v)[1],
				      0));
	
  for (unsigned int i=0; i<4; ++i)
    {
      for (typename Triangulation<dim>::active_cell_iterator
	     cell = triangulation.begin_active();
	   cell != triangulation.end(); ++cell)
	for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
	  if (cell->face(f)->boundary_indicator() == 1)
	    cell->set_refine_flag ();
      
      triangulation.execute_coarsening_and_refinement ();

      std::cout << "Refinement cycle " << i << std::endl
		<< "   Number of active cells: "
		<< triangulation.n_active_cells()
		<< std::endl
		<< "   Total number of cells: "
		<< triangulation.n_cells()
		<< std::endl;

    }
  
  
  dof_handler.distribute_dofs (fe);

  std::cout << "   Number of degrees of freedom: "
	    << dof_handler.n_dofs()
	    << std::endl;

  sparsity_pattern.reinit (dof_handler.n_dofs(),
			   dof_handler.n_dofs(),
			   dof_handler.max_couplings_between_dofs());
  DoFTools::make_sparsity_pattern (dof_handler, sparsity_pattern);
  sparsity_pattern.compress();

  system_matrix.reinit (sparsity_pattern);

  solution.reinit (dof_handler.n_dofs());
  system_rhs.reinit (dof_handler.n_dofs());
}


                                 // @sect4{LaplaceProblem::assemble_system}

				 // Unlike in the previous example, we
				 // would now like to use a
				 // non-constant right hand side
				 // function and non-zero boundary
				 // values. Both are tasks that are
				 // readily achieved with a only a few
				 // new lines of code in the
				 // assemblage of the matrix and right
				 // hand side.
				 //
				 // More interesting, though, is the
				 // way we assemble matrix and right
				 // hand side vector dimension
				 // independently: there is simply no
				 // difference to the 
				 // two-dimensional case. Since the
				 // important objects used in this
				 // function (quadrature formula,
				 // FEValues) depend on the dimension
				 // by way of a template parameter as
				 // well, they can take care of
				 // setting up properly everything for
				 // the dimension for which this
				 // function is compiled. By declaring
				 // all classes which might depend on
				 // the dimension using a template
				 // parameter, the library can make
				 // nearly all work for you and you
				 // don't have to care about most
				 // things.
template <int dim>
void LaplaceProblem<dim>::assemble_system () 
{  
  MatrixTools::create_laplace_matrix (dof_handler,
				      QGauss<dim>(2),
				      system_matrix);
  system_rhs = 0;
  
  std::map<unsigned int,double> boundary_values;
  VectorTools::interpolate_boundary_values (dof_handler,
					    0,
					    BoundaryValues<dim>(),
					    boundary_values);
  MatrixTools::apply_boundary_values (boundary_values,
				      system_matrix,
				      solution,
				      system_rhs);
}


                                 // @sect4{LaplaceProblem::solve}

				 // Solving the linear system of
				 // equations is something that looks
				 // almost identical in most
				 // programs. In particular, it is
				 // dimension independent, so this
				 // function is copied verbatim from the
				 // previous example.
template <int dim>
void LaplaceProblem<dim>::solve () 
{
				   // NEW
  SolverControl           solver_control (dof_handler.n_dofs(),
					  1e-12*system_rhs.l2_norm());
  SolverCG<>              cg (solver_control);

  PreconditionSSOR<> preconditioner;
  preconditioner.initialize(system_matrix, 1.2);

  cg.solve (system_matrix, solution, system_rhs,
	    preconditioner);
}


                                 // @sect4{LaplaceProblem::output_results}

				 // This function also does what the
				 // respective one did in step-3. No changes
				 // here for dimension independence either.
                                 //
                                 // The only difference to the previous
                                 // example is that we want to write output in
                                 // GMV format, rather than for gnuplot (GMV
                                 // is another graphics program that, contrary
                                 // to gnuplot, shows data in nice colors,
                                 // allows rotation of geometries with the
                                 // mouse, and generates reasonable
                                 // representations of 3d data; for ways to
                                 // obtain it see the ReadMe file of
                                 // deal.II). To write data in this format, we
                                 // simply replace the
                                 // <code>data_out.write_gnuplot</code> call by
                                 // <code>data_out.write_gmv</code>.
                                 //
                                 // Since the program will run both 2d and 3d
                                 // versions of the laplace solver, we use the
                                 // dimension in the filename to generate
                                 // distinct filenames for each run (in a
                                 // better program, one would check whether
                                 // `dim' can have other values than 2 or 3,
                                 // but we neglect this here for the sake of
                                 // brevity).
template <int dim>
void LaplaceProblem<dim>::output_results () const
{
  DataOut<dim> data_out;

  data_out.attach_dof_handler (dof_handler);
  data_out.add_data_vector (solution, "solution");

  data_out.build_patches ();

  std::ofstream output (dim == 2 ?
			"solution-2d.gmv" :
			"solution-3d.gmv");
  data_out.write_gmv (output);
}



                                 // @sect4{LaplaceProblem::run}

                                 // This is the function which has the
				 // top-level control over
				 // everything. Apart from one line of
				 // additional output, it is the same
				 // as for the previous example.
template <int dim>
void LaplaceProblem<dim>::run () 
{
  std::cout << "Solving problem in " << dim << " space dimensions." << std::endl;
  
  make_grid_and_dofs();
  assemble_system ();
  solve ();
  output_results ();
}


                                 // @sect3{The <code>main</code> function}

				 // And this is the main function. It also
				 // looks mostly like in step-3, but if you
				 // look at the code below, note how we first
				 // create a variable of type
				 // <code>LaplaceProblem@<2@></code> (forcing the
				 // compiler to compile the class template
				 // with <code>dim</code> replaced by <code>2</code>) and run a
				 // 2d simulation, and then we do the whole
				 // thing over in 3d.
				 //
				 // In practice, this is probably not what you
				 // would do very frequently (you probably
				 // either want to solve a 2d problem, or one
				 // in 3d, but not both at the same
				 // time). However, it demonstrates the
				 // mechanism by which we can simply change
				 // which dimension we want in a single place,
				 // and thereby force the compiler to
				 // recompile the dimension independent class
				 // templates for the dimension we
				 // request. The emphasis here lies on the
				 // fact that we only need to change a single
				 // place. This makes it rather trivial to
				 // debug the program in 2d where computations
				 // are fast, and then switch a single place
				 // to a 3 to run the much more computing
				 // intensive program in 3d for `real'
				 // computations.
				 //
				 // Each of the two blocks is enclosed in
				 // braces to make sure that the
				 // <code>laplace_problem_2d</code> variable goes out
				 // of scope (and releases the memory it
				 // holds) before we move on to allocate
				 // memory for the 3d case. Without the
				 // additional braces, the
				 // <code>laplace_problem_2d</code> variable would only
				 // be destroyed at the end of the function,
				 // i.e. after running the 3d problem, and
				 // would needlessly hog memory while the 3d
				 // run could actually use it.
                                 //
                                 // Finally, the first line of the function is
                                 // used to suppress some output.  Remember
                                 // that in the previous example, we had the
                                 // output from the linear solvers about the
                                 // starting residual and the number of the
                                 // iteration where convergence was
                                 // detected. This can be suppressed through
                                 // the <code>deallog.depth_console(0)</code> call.
                                 //
                                 // The rationale here is the following: the
                                 // deallog (i.e. deal-log, not de-allog)
                                 // variable represents a stream to which some
                                 // parts of the library write output. It
                                 // redirects this output to the console and
                                 // if required to a file. The output is
                                 // nested in a way so that each function can
                                 // use a prefix string (separated by colons)
                                 // for each line of output; if it calls
                                 // another function, that may also use its
                                 // prefix which is then printed after the one
                                 // of the calling function. Since output from
                                 // functions which are nested deep below is
                                 // usually not as important as top-level
                                 // output, you can give the deallog variable
                                 // a maximal depth of nested output for
                                 // output to console and file. The depth zero
                                 // which we gave here means that no output is
                                 // written. By changing it you can get more
                                 // information about the innards of the
                                 // library.
int main () 
{
  deallog.depth_console (0);
  {
    LaplaceProblem<3> laplace_problem_3d;
    laplace_problem_3d.run ();
  }
  
  return 0;
}
