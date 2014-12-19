// ---------------------------------------------------------------------
//
// Copyright (C) 2003 - 2013 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------



// a testcase by Joa, see mailing list 2007/02/24. checks that we can find a
// point inside a cell, where the point is definitely inside the domain. this
// testcase is a simple modification of step-6, only a few lines are added

#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/grid/tria.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/data_out.h>

#include <fstream>
#include <iomanip>

// From the following include file we
// will import the declaration of
// H1-conforming finite element shape
// functions. This family of finite
// elements is called <code>FE_Q</code>, and
// was used in all examples before
// already to define the usual bi- or
// tri-linear elements, but we will
// now use it for bi-quadratic
// elements:
#include <deal.II/fe/fe_q.h>
// We will not read the grid from a
// file as in the previous example,
// but generate it using a function
// of the library. However, we will
// want to write out the locally
// refined grids (just the grid, not
// the solution) in each step, so we
// need the following include file
// instead of <code>grid_in.h</code>:
#include <deal.II/grid/grid_out.h>


// When using locally refined grids,
// we will get so-called <code>hanging
// nodes</code>. However, the standard
// finite element methods assumes
// that the discrete solution spaces
// be continuous, so we need to make
// sure that the degrees of freedom
// on hanging nodes conform to some
// constraints such that the global
// solution is continuous. The
// following file contains a class
// which is used to handle these
// constraints:
#include <deal.II/lac/constraint_matrix.h>

// In order to refine our grids
// locally, we need a function from
// the library that decides which
// cells to flag for refinement or
// coarsening based on the error
// indicators we have computed. This
// function is defined here:
#include <deal.II/grid/grid_refinement.h>

// Finally, we need a simple way to
// actually compute the refinement
// indicators based on some error
// estimat. While in general,
// adaptivity is very
// problem-specific, the error
// indicator in the following file
// often yields quite nicely adapted
// grids for a wide class of
// problems.
#include <deal.II/numerics/error_estimator.h>

// Finally, this is as in previous
// programs:
using namespace dealii;


// @sect3{The <code>LaplaceProblem</code> class template}

// The main class is again almost
// unchanged. Two additions, however,
// are made: we have added the
// <code>refine_grid</code> function, which is
// used to adaptively refine the grid
// (instead of the global refinement
// in the previous examples), and a
// variable which will hold the
// constraints associated to the
// hanging nodes. In addition, we
// have added a destructor to the
// class for reasons that will become
// clear when we discuss its
// implementation.
template <int dim>
class LaplaceProblem
{
public:
  LaplaceProblem ();
  ~LaplaceProblem ();

  void run ();

private:
  void setup_system ();
  void assemble_system ();
  void solve ();
  void refine_grid ();
  void output_results (const unsigned int cycle) const;

  Triangulation<dim>   triangulation;

  DoFHandler<dim>      dof_handler;
  FE_Q<dim>            fe;

  // This is the new variable in
  // the main class. We need an
  // object which holds a list of
  // constraints originating from
  // the hanging nodes:
  ConstraintMatrix     hanging_node_constraints;

  SparsityPattern      sparsity_pattern;
  SparseMatrix<double> system_matrix;

  Vector<double>       solution;
  Vector<double>       system_rhs;
};


// @sect3{Nonconstant coefficients}

// The implementation of nonconstant
// coefficients is copied verbatim
// from step-5:

template <int dim>
class Coefficient : public Function<dim>
{
public:
  Coefficient () : Function<dim>() {}

  virtual double value (const Point<dim>   &p,
                        const unsigned int  component = 0) const;

  virtual void value_list (const std::vector<Point<dim> > &points,
                           std::vector<double>            &values,
                           const unsigned int              component = 0) const;
};



template <int dim>
double Coefficient<dim>::value (const Point<dim> &p,
                                const unsigned int) const
{
  if (p.square() < 0.5*0.5)
    return 20;
  else
    return 1;
}



template <int dim>
void Coefficient<dim>::value_list (const std::vector<Point<dim> > &points,
                                   std::vector<double>            &values,
                                   const unsigned int              component) const
{
  const unsigned int n_points = points.size();

  Assert (values.size() == n_points,
          ExcDimensionMismatch (values.size(), n_points));

  Assert (component == 0,
          ExcIndexRange (component, 0, 1));

  for (unsigned int i=0; i<n_points; ++i)
    {
      if (points[i].square() < 0.5*0.5)
        values[i] = 20;
      else
        values[i] = 1;
    }
}


// @sect3{The <code>LaplaceProblem</code> class implementation}

// @sect4{LaplaceProblem::LaplaceProblem}

// The constructor of this class is
// mostly the same as before, but
// this time we want to use the
// quadratic element. To do so, we
// only have to replace the
// constructor argument (which was
// <code>1</code> in all previous examples) by
// the desired polynomial degree
// (here <code>2</code>):
template <int dim>
LaplaceProblem<dim>::LaplaceProblem () :
  dof_handler (triangulation),
  fe (2)
{}


// @sect4{LaplaceProblem::~LaplaceProblem}

// Here comes the added destructor of
// the class. The reason why we want
// to add it is a subtle change in
// the order of data elements in the
// class as compared to all previous
// examples: the <code>dof_handler</code>
// object was defined before and not
// after the <code>fe</code> object. Of course
// we could have left this order
// unchanged, but we would like to
// show what happens if the order is
// reversed since this produces a
// rather nasty side-effect and
// results in an error which is
// difficult to track down if one
// does not know what happens.
//
// Basically what happens is the
// following: when we distribute the
// degrees of freedom using the
// function call
// <code>dof_handler.distribute_dofs()</code>,
// the <code>dof_handler</code> also stores a
// pointer to the finite element in
// use. Since this pointer is used
// every now and then until either
// the degrees of freedom are
// re-distributed using another
// finite element object or until the
// <code>dof_handler</code> object is
// destroyed, it would be unwise if
// we would allow the finite element
// object to be deleted before the
// <code>dof_handler</code> object. To
// disallow this, the DoF handler
// increases a counter inside the
// finite element object which counts
// how many objects use that finite
// element (this is what the
// <code>Subscriptor</code>/<code>SmartPointer</code>
// class pair is used for, in case
// you want something like this for
// your own programs; see step-7 for
// a more complete discussion
// of this topic). The finite
// element object will refuse its
// destruction if that counter is
// larger than zero, since then some
// other objects might rely on the
// persistence of the finite element
// object. An exception will then be
// thrown and the program will
// usually abort upon the attempt to
// destroy the finite element.
//
// To be fair, such exceptions about
// still used objects are not
// particularly popular among
// programmers using deal.II, since
// they only tell us that something
// is wrong, namely that some other
// object is still using the object
// that is presently being
// destructed, but most of the time
// not who this user is. It is
// therefore often rather
// time-consuming to find out where
// the problem exactly is, although
// it is then usually straightforward
// to remedy the situation. However,
// we believe that the effort to find
// invalid references to objects that
// do no longer exist is less if the
// problem is detected once the
// reference becomes invalid, rather
// than when non-existent objects are
// actually accessed again, since
// then usually only invalid data is
// accessed, but no error is
// immediately raised.
//
// Coming back to the present
// situation, if we did not write
// this destructor, the compiler will
// generate code that triggers
// exactly the behavior sketched
// above. The reason is that member
// variables of the
// <code>LaplaceProblem</code> class are
// destructed bottom-up (i.e. in
// reverse order of their declaration
// in the class), as always in
// C++. Thus, the finite element
// object will be destructed before
// the DoF handler object, since its
// declaration is below the one of
// the DoF handler. This triggers the
// situation above, and an exception
// will be raised when the <code>fe</code>
// object is destructed. What needs
// to be done is to tell the
// <code>dof_handler</code> object to release
// its lock to the finite element. Of
// course, the <code>dof_handler</code> will
// only release its lock if it really
// does not need the finite element
// any more, i.e. when all finite
// element related data is deleted
// from it. For this purpose, the
// <code>DoFHandler</code> class has a
// function <code>clear</code> which deletes
// all degrees of freedom, and
// releases its lock to the finite
// element. After this, you can
// safely destruct the finite element
// object since its internal counter
// is then zero.
//
// For completeness, we add the
// output of the exception that would
// have been triggered without this
// destructor, to the end of the
// results section of this example.
template <int dim>
LaplaceProblem<dim>::~LaplaceProblem ()
{
  dof_handler.clear ();
}


// @sect4{LaplaceProblem::setup_system}

// The next function is setting up
// all the variables that describe
// the linear finite element problem,
// such as the DoF handler, the
// matrices, and vectors. The
// difference to what we did in
// step-5 is only that we now also
// have to take care of handing node
// constraints. These constraints are
// handled almost transparently by
// the library, i.e. you only need to
// know that they exist and how to
// get them, but you do not have to
// know how they are formed or what
// exactly is done with them.
//
// At the beginning of the function,
// you find all the things that are
// the same as in step-5: setting up
// the degrees of freedom (this time
// we have quadratic elements, but
// there is no difference from a user
// code perspective to the linear --
// or cubic, for that matter --
// case), generating the sparsity
// pattern, and initializing the
// solution and right hand side
// vectors. Note that the sparsity
// pattern will have significantly
// more entries per row now, since
// there are now 9 degrees of freedom
// per cell, not only four, that can
// couple with each other. The
// <code>dof_Handler.max_couplings_between_dofs()</code>
// call will take care of this,
// however:
template <int dim>
void LaplaceProblem<dim>::setup_system ()
{
  dof_handler.distribute_dofs (fe);

  sparsity_pattern.reinit (dof_handler.n_dofs(),
                           dof_handler.n_dofs(),
                           dof_handler.max_couplings_between_dofs());
  DoFTools::make_sparsity_pattern (dof_handler, sparsity_pattern);

  solution.reinit (dof_handler.n_dofs());
  system_rhs.reinit (dof_handler.n_dofs());


  // After setting up all the degrees
  // of freedoms, here are now the
  // differences compared to step-5,
  // all of which are related to
  // constraints associated with the
  // hanging nodes. In the class
  // desclaration, we have already
  // allocated space for an object
  // <code>hanging_node_constraints</code>
  // that will hold a list of these
  // constraints (they form a matrix,
  // which is reflected in the name
  // of the class, but that is
  // immaterial for the moment). Now
  // we have to fill this
  // object. This is done using the
  // following function calls (the
  // first clears the contents of the
  // object that may still be left
  // over from computations on the
  // previous mesh before the last
  // adaptive refinement):
  hanging_node_constraints.clear ();
  DoFTools::make_hanging_node_constraints (dof_handler,
                                           hanging_node_constraints);

  // The next step is <code>closing</code>
  // this object. For this note that,
  // in principle, the
  // <code>ConstraintMatrix</code> class can
  // hold other constraints as well,
  // i.e. constraints that do not
  // stem from hanging
  // nodes. Sometimes, it is useful
  // to use such constraints, in
  // which case they may be added to
  // the <code>ConstraintMatrix</code> object
  // after the hanging node
  // constraints were computed. After
  // all constraints have been added,
  // they need to be sorted and
  // rearranged to perform some
  // actions more efficiently. This
  // postprocessing is done using the
  // <code>close()</code> function, after which
  // no further constraints may be
  // added any more:
  hanging_node_constraints.close ();

  // The constrained hanging nodes
  // will later be eliminated from
  // the linear system of
  // equations. When doing so, some
  // additional entries in the global
  // matrix will be set to non-zero
  // values, so we have to reserve
  // some space for them here. Since
  // the process of elimination of
  // these constrained nodes is
  // called <code>condensation</code>, the
  // functions that eliminate them
  // are called <code>condense</code> for both
  // the system matrix and right hand
  // side, as well as for the
  // sparsity pattern.
  hanging_node_constraints.condense (sparsity_pattern);

  // Now all non-zero entries of the
  // matrix are known (i.e. those
  // from regularly assembling the
  // matrix and those that were
  // introduced by eliminating
  // constraints). We can thus close
  // the sparsity pattern and remove
  // unneeded space:
  sparsity_pattern.compress();

  // Finally, the so-constructed
  // sparsity pattern serves as the
  // basis on top of which we will
  // create the sparse matrix:
  system_matrix.reinit (sparsity_pattern);
}

// @sect4{LaplaceProblem::assemble_system}

// Next, we have to assemble the
// matrix again. There are no code
// changes compared to step-5 except
// for a single place: We have to use
// a higher-order quadrature formula
// to account for the higher
// polynomial degree in the finite
// element shape functions. This is
// easy to change: the constructor of
// the <code>QGauss</code> class takes the
// number of quadrature points in
// each space direction. Previously,
// we had two points for bilinear
// elements. Now we should use three
// points for biquadratic elements.
//
// The rest of the code that forms
// the local contributions and
// transfers them into the global
// objects remains unchanged. It is
// worth noting, however, that under
// the hood several things are
// different than before. First, the
// variables <code>dofs_per_cell</code> and
// <code>n_q_points</code> now are 9 each,
// where they were 4
// before. Introducing such variables
// as abbreviations is a good
// strategy to make code work with
// different elements without having
// to change too much code. Secondly,
// the <code>fe_values</code> object of course
// needs to do other things as well,
// since the shape functions are now
// quadratic, rather than linear, in
// each coordinate variable. Again,
// however, this is something that is
// completely transparent to user
// code and nothing that you have to
// worry about.
template <int dim>
void LaplaceProblem<dim>::assemble_system ()
{
  const QGauss<dim>  quadrature_formula(3);

  FEValues<dim> fe_values (fe, quadrature_formula,
                           update_values    |  update_gradients |
                           update_q_points  |  update_JxW_values);

  const unsigned int   dofs_per_cell = fe.dofs_per_cell;
  const unsigned int   n_q_points    = quadrature_formula.size();

  FullMatrix<double>   cell_matrix (dofs_per_cell, dofs_per_cell);
  Vector<double>       cell_rhs (dofs_per_cell);

  std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);

  const Coefficient<dim> coefficient;
  std::vector<double>    coefficient_values (n_q_points);

  typename DoFHandler<dim>::active_cell_iterator
  cell = dof_handler.begin_active(),
  endc = dof_handler.end();
  for (; cell!=endc; ++cell)
    {
      cell_matrix = 0;
      cell_rhs = 0;

      fe_values.reinit (cell);

      coefficient.value_list (fe_values.get_quadrature_points(),
                              coefficient_values);

      for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
        for (unsigned int i=0; i<dofs_per_cell; ++i)
          {
            for (unsigned int j=0; j<dofs_per_cell; ++j)
              cell_matrix(i,j) += (coefficient_values[q_point] *
                                   fe_values.shape_grad(i,q_point) *
                                   fe_values.shape_grad(j,q_point) *
                                   fe_values.JxW(q_point));

            cell_rhs(i) += (fe_values.shape_value(i,q_point) *
                            1.0 *
                            fe_values.JxW(q_point));
          }

      cell->get_dof_indices (local_dof_indices);
      for (unsigned int i=0; i<dofs_per_cell; ++i)
        {
          for (unsigned int j=0; j<dofs_per_cell; ++j)
            system_matrix.add (local_dof_indices[i],
                               local_dof_indices[j],
                               cell_matrix(i,j));

          system_rhs(local_dof_indices[i]) += cell_rhs(i);
        }
    }

  // After the system of equations
  // has been assembled just as for
  // the previous examples, we still
  // have to eliminate the
  // constraints due to hanging
  // nodes. This is done using the
  // following two function calls:
  hanging_node_constraints.condense (system_matrix);
  hanging_node_constraints.condense (system_rhs);
  // Using them, degrees of freedom
  // associated to hanging nodes have
  // been removed from the linear
  // system and the independent
  // variables are only the regular
  // nodes. The constrained nodes are
  // still in the linear system
  // (there is a one on the diagonal
  // of the matrix and all other
  // entries for this line are set to
  // zero) but the computed values
  // are invalid (the <code>condense</code>
  // function modifies the system so
  // that the values in the solution
  // corresponding to constrained
  // nodes are invalid, but that the
  // system still has a well-defined
  // solution; we compute the correct
  // values for these nodes at the
  // end of the <code>solve</code> function).

  // As almost all the stuff before,
  // the interpolation of boundary
  // values works also for higher
  // order elements without the need
  // to change your code for that. We
  // note that for proper results, it
  // is important that the
  // elimination of boundary nodes
  // from the system of equations
  // happens *after* the elimination
  // of hanging nodes.
  std::map<types::global_dof_index,double> boundary_values;
  VectorTools::interpolate_boundary_values (dof_handler,
                                            0,
                                            ZeroFunction<dim>(),
                                            boundary_values);
  MatrixTools::apply_boundary_values (boundary_values,
                                      system_matrix,
                                      solution,
                                      system_rhs);
}



// @sect4{LaplaceProblem::solve}

// We continue with gradual
// improvements. The function that
// solves the linear system again
// uses the SSOR preconditioner, and
// is again unchanged except that we
// have to incorporate hanging node
// constraints. As mentioned above,
// the degrees of freedom
// corresponding to hanging node
// constraints have been removed from
// the linear system by giving the
// rows and columns of the matrix a
// special treatment. This way, the
// values for these degrees of
// freedom have wrong, but
// well-defined values after solving
// the linear system. What we then
// have to do is to use the
// constraints to assign to them the
// values that they should have. This
// process, called <code>distributing</code>
// hanging nodes, computes the values
// of constrained nodes from the
// values of the unconstrained ones,
// and requires only a single
// additional function call that you
// find at the end of this function:

template <int dim>
void LaplaceProblem<dim>::solve ()
{
  SolverControl           solver_control (1000, 1e-12);
  SolverCG<>              cg (solver_control);

  PreconditionSSOR<> preconditioner;
  preconditioner.initialize(system_matrix, 1.2);

  cg.solve (system_matrix, solution, system_rhs,
            preconditioner);

  hanging_node_constraints.distribute (solution);
}


// @sect4{LaplaceProblem::refine_grid}

// Instead of global refinement, we
// now use a slightly more elaborate
// scheme. We will use the
// <code>KellyErrorEstimator</code> class
// which implements an error
// estimator for the Laplace
// equation; it can in principle
// handle variable coefficients, but
// we will not use these advanced
// features, but rather use its most
// simple form since we are not
// interested in quantitative results
// but only in a quick way to
// generate locally refined grids.
//
// Although the error estimator
// derived by Kelly et al. was
// originally developed for the Laplace
// equation, we have found that it is
// also well suited to quickly
// generate locally refined grids for
// a wide class of
// problems. Basically, it looks at
// the jumps of the gradients of the
// solution over the faces of cells
// (which is a measure for the second
// derivatives) and scales it by the
// size of the cell. It is therefore
// a measure for the local smoothness
// of the solution at the place of
// each cell and it is thus
// understandable that it yields
// reasonable grids also for
// hyperbolic transport problems or
// the wave equation as well,
// although these grids are certainly
// suboptimal compared to approaches
// specially tailored to the
// problem. This error estimator may
// therefore be understood as a quick
// way to test an adaptive program.
//
// The way the estimator works is to
// take a <code>DoFHandler</code> object
// describing the degrees of freedom
// and a vector of values for each
// degree of freedom as input and
// compute a single indicator value
// for each active cell of the
// triangulation (i.e. one value for
// each of the
// <code>triangulation.n_active_cells()</code>
// cells). To do so, it needs two
// additional pieces of information:
// a quadrature formula on the faces
// (i.e. quadrature formula on
// <code>dim-1</code> dimensional objects. We
// use a 3-point Gauss rule again, a
// pick that is consistent and
// appropriate with the choice
// bi-quadratic finite element shape
// functions in this program.
// (What constitutes a suitable
// quadrature rule here of course
// depends on knowledge of the way
// the error estimator evaluates
// the solution field. As said
// above, the jump of the gradient
// is integrated over each face,
// which would be a quadratic
// function on each face for the
// quadratic elements in use in
// this example. In fact, however,
// it is the square of the jump of
// the gradient, as explained in
// the documentation of that class,
// and that is a quartic function,
// for which a 3 point Gauss
// formula is sufficient since it
// integrates polynomials up to
// order 5 exactly.)
//
// Secondly, the function wants a
// list of boundaries where we have
// imposed Neumann value, and the
// corresponding Neumann values. This
// information is represented by an
// object of type
// <code>FunctionMap@<dim@>::type</code> that is
// essentially a map from boundary
// indicators to function objects
// describing Neumann boundary values
// (in the present example program,
// we do not use Neumann boundary
// values, so this map is empty, and
// in fact constructed using the
// default constructor of the map in
// the place where the function call
// expects the respective function
// argument).
//
// The output, as mentioned is a
// vector of values for all
// cells. While it may make sense to
// compute the *value* of a degree of
// freedom very accurately, it is
// usually not helpful to compute the
// *error indicator* corresponding to
// a cell particularly accurately. We
// therefore typically use a vector
// of floats instead of a vector of
// doubles to represent error
// indicators.
template <int dim>
void LaplaceProblem<dim>::refine_grid ()
{
  Vector<float> estimated_error_per_cell (triangulation.n_active_cells());

  KellyErrorEstimator<dim>::estimate (dof_handler,
                                      QGauss<dim-1>(3),
                                      typename FunctionMap<dim>::type(),
                                      solution,
                                      estimated_error_per_cell);

  // The above function returned one
  // error indicator value for each
  // cell in the
  // <code>estimated_error_per_cell</code>
  // array. Refinement is now done as
  // follows: refine those 30 per
  // cent of the cells with the
  // highest error values, and
  // coarsen the 3 per cent of cells
  // with the lowest values.
  //
  // One can easily verify that if
  // the second number were zero,
  // this would approximately result
  // in a doubling of cells in each
  // step in two space dimensions,
  // since for each of the 30 per
  // cent of cells, four new would be
  // replaced, while the remaining 70
  // per cent of cells remain
  // untouched. In practice, some
  // more cells are usually produced
  // since it is disallowed that a
  // cell is refined twice while the
  // neighbor cell is not refined; in
  // that case, the neighbor cell
  // would be refined as well.
  //
  // In many applications, the number
  // of cells to be coarsened would
  // be set to something larger than
  // only three per cent. A non-zero
  // value is useful especially if
  // for some reason the initial
  // (coarse) grid is already rather
  // refined. In that case, it might
  // be necessary to refine it in
  // some regions, while coarsening
  // in some other regions is
  // useful. In our case here, the
  // initial grid is very coarse, so
  // coarsening is only necessary in
  // a few regions where
  // over-refinement may have taken
  // place. Thus a small, non-zero
  // value is appropriate here.
  //
  // The following function now takes
  // these refinement indicators and
  // flags some cells of the
  // triangulation for refinement or
  // coarsening using the method
  // described above. It is from a
  // class that implements
  // several different algorithms to
  // refine a triangulation based on
  // cell-wise error indicators.
  GridRefinement::refine_and_coarsen_fixed_number (triangulation,
                                                   estimated_error_per_cell,
                                                   0.3, 0.03);

  // After the previous function has
  // exited, some cells are flagged
  // for refinement, and some other
  // for coarsening. The refinement
  // or coarsening itself is not
  // performed by now, however, since
  // there are cases where further
  // modifications of these flags is
  // useful. Here, we don't want to
  // do any such thing, so we can
  // tell the triangulation to
  // perform the actions for which
  // the cells are flagged:
  triangulation.execute_coarsening_and_refinement ();
}


// @sect4{LaplaceProblem::output_results}

// At the end of computations on each
// grid, and just before we continue
// the next cycle with mesh
// refinement, we want to output the
// results from this cycle.
//
// In the present program, we will
// not write the solution (except for
// in the last step, see the next
// function), but only the meshes
// that we generated, as a
// two-dimensional Encapsulated
// Postscript (EPS) file.
//
// We have already seen in step-1 how
// this can be achieved. The only
// thing we have to change is the
// generation of the file name, since
// it should contain the number of
// the present refinement cycle
// provided to this function as an
// argument. The most general way is
// to use the std::stringstream class
// as shown in step-5, but here's a
// little hack that makes it simpler
// if we know that we have less than
// 10 iterations: assume that the
// numbers `0' through `9' are
// represented consecutively in the
// character set used on your machine
// (this is in fact the case in all
// known character sets), then
// '0'+cycle gives the character
// corresponding to the present cycle
// number. Of course, this will only
// work if the number of cycles is
// actually less than 10, and rather
// than waiting for the disaster to
// happen, we safeguard our little
// hack with an explicit assertion at
// the beginning of the function. If
// this assertion is triggered,
// i.e. when <code>cycle</code> is larger than
// or equal to 10, an exception of
// type <code>ExcNotImplemented</code> is
// raised, indicating that some
// functionality is not implemented
// for this case (the functionality
// that is missing, of course, is the
// generation of file names for that
// case):
template <int dim>
void LaplaceProblem<dim>::output_results (const unsigned int cycle) const
{
  Assert (cycle < 10, ExcNotImplemented());

  std::string filename = "grid-";
  filename += ('0' + cycle);
  filename += ".eps";

  std::ofstream output (filename.c_str());

  GridOut grid_out;
  grid_out.write_eps (triangulation, output);
}



// @sect4{LaplaceProblem::run}

// The final function before
// <code>main()</code> is again the main
// driver of the class, <code>run()</code>. It
// is similar to the one of step-5,
// except that we generate a file in
// the program again instead of
// reading it from disk, in that we
// adaptively instead of globally
// refine the mesh, and that we
// output the solution on the final
// mesh in the present function.
//
// The first block in the main loop
// of the function deals with mesh
// generation. If this is the first
// cycle of the program, instead of
// reading the grid from a file on
// disk as in the previous example,
// we now again create it using a
// library function. The domain is
// again a circle, which is why we
// have to provide a suitable
// boundary object as well. We place
// the center of the circle at the
// origin and have the radius be one
// (these are the two hidden
// arguments to the function, which
// have default values).
//
// You will notice by looking at the
// coarse grid that it is of inferior
// quality than the one which we read
// from the file in the previous
// example: the cells are less
// equally formed. However, using the
// library function this program
// works in any space dimension,
// which was not the case before.
//
// In case we find that this is not
// the first cycle, we want to refine
// the grid. Unlike the global
// refinement employed in the last
// example program, we now use the
// adaptive procedure described
// above.
//
// The rest of the loop looks as
// before:
template <int dim>
void LaplaceProblem<dim>::run ()
{
  for (unsigned int cycle=0; cycle<8; ++cycle)
    {
      deallog << "Cycle " << cycle << ':' << std::endl;

      if (cycle == 0)
        {
          GridGenerator::hyper_ball (triangulation);

          static const HyperBallBoundary<dim> boundary;
          triangulation.set_boundary (0, boundary);

          triangulation.refine_global (1);
        }
      else
        refine_grid ();


      deallog << "   Number of active cells:       "
              << triangulation.n_active_cells()
              << std::endl;

      setup_system ();

      deallog << "   Number of degrees of freedom: "
              << dof_handler.n_dofs()
              << std::endl;

      assemble_system ();
      solve ();
    }

  // try to find a bunch of points that are
  // definitely inside the domain (we here
  // have a circle of radius 1, so find
  // points inside a radius of 0.9)
  for (int i=0; i<1000; i++)
    {
      double r = sqrt((0.9*Testing::rand()/RAND_MAX));
      double phi = 2*3.14*(1.0*Testing::rand()/RAND_MAX);
      double x = r*cos(phi);
      double y = r*sin(phi);
      Point<2> p(x,y);
      VectorTools::point_value(dof_handler,solution,p);
    }
}


// @sect3{The <code>main</code> function}

// The main function is unaltered in
// its functionality from the
// previous example, but we have
// taken a step of additional
// caution. Sometimes, something goes
// wrong (such as insufficient disk
// space upon writing an output file,
// not enough memory when trying to
// allocate a vector or a matrix, or
// if we can't read from or write to
// a file for whatever reason), and
// in these cases the library will
// throw exceptions. Since these are
// run-time problems, not programming
// errors that can be fixed once and
// for all, this kind of exceptions
// is not switched off in optimized
// mode, in contrast to the
// <code>Assert</code> macro which we have
// used to test against programming
// errors. If uncaught, these
// exceptions propagate the call tree
// up to the <code>main</code> function, and
// if they are not caught there
// either, the program is aborted. In
// many cases, like if there is not
// enough memory or disk space, we
// can't do anything but we can at
// least print some text trying to
// explain the reason why the program
// failed. A way to do so is shown in
// the following. It is certainly
// useful to write any larger program
// in this way, and you can do so by
// more or less copying this function
// except for the <code>try</code> block that
// actually encodes the functionality
// particular to the present
// application.
int main ()
{
  std::ofstream logfile ("output");
  deallog << std::setprecision (3);
  deallog << std::fixed;
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  // The general idea behind the
  // layout of this function is as
  // follows: let's try to run the
  // program as we did before...
  try
    {
      LaplaceProblem<2> laplace_problem_2d;
      laplace_problem_2d.run ();
    }
  // ...and if this should fail, try
  // to gather as much information as
  // possible. Specifically, if the
  // exception that was thrown is an
  // object of a class that is
  // derived from the C++ standard
  // class <code>exception</code>, then we can
  // use the <code>what</code> member function
  // to get a string which describes
  // the reason why the exception was
  // thrown.
  //
  // The deal.II exception classes
  // are all derived from the
  // standard class, and in
  // particular, the <code>exc.what()</code>
  // function will return
  // approximately the same string as
  // would be generated if the
  // exception was thrown using the
  // <code>Assert</code> macro. You have seen
  // the output of such an exception
  // in the previous example, and you
  // then know that it contains the
  // file and line number of where
  // the exception occured, and some
  // other information. This is also
  // what the following statements
  // would print.
  //
  // Apart from this, there isn't
  // much that we can do except
  // exiting the program with an
  // error code (this is what the
  // <code>return 1;</code> does):
  catch (std::exception &exc)
    {
      deallog << std::endl << std::endl
              << "----------------------------------------------------"
              << std::endl;
      deallog << "Exception on processing: " << std::endl
              << exc.what() << std::endl
              << "Aborting!" << std::endl
              << "----------------------------------------------------"
              << std::endl;

      return 1;
    }
  // If the exception that was thrown
  // somewhere was not an object of a
  // class derived from the standard
  // <code>exception</code> class, then we
  // can't do anything at all. We
  // then simply print an error
  // message and exit.
  catch (...)
    {
      deallog << std::endl << std::endl
              << "----------------------------------------------------"
              << std::endl;
      deallog << "Unknown exception!" << std::endl
              << "Aborting!" << std::endl
              << "----------------------------------------------------"
              << std::endl;
      return 1;
    }

  // If we got to this point, there
  // was no exception which
  // propagated up to the main
  // function (there may have been
  // exceptions, but they were caught
  // somewhere in the program or the
  // library). Therefore, the program
  // performed as was expected and we
  // can return without error.
  return 0;
}
