/* $Id$ */
/* Author: Wolfgang Bangerth, University of Heidelberg, 2002 */

/*    $Id$       */
/*    Version: $Name$                                          */
/*                                                                */
/*    Copyright (C) 2002, 2003, 2004 by the deal.II authors                   */
/*                                                                */
/*    This file is subject to QPL and may not be  distributed     */
/*    without copyright and license information. Please refer     */
/*    to the file deal.II/doc/license.html for the  text  and     */
/*    further information on this license.                        */

                                 // As usual, most of the headers here have
                                 // already been used and discussed in
                                 // previous examples:
#include <base/quadrature_lib.h>
#include <base/function.h>
#include <base/logstream.h>
#include <lac/vector.h>
#include <lac/full_matrix.h>
#include <lac/sparse_matrix.h>
#include <lac/solver_cg.h>
#include <lac/vector_memory.h>
#include <lac/precondition.h>
#include <grid/tria.h>
#include <grid/grid_generator.h>
#include <grid/grid_refinement.h>
#include <grid/tria_accessor.h>
#include <grid/tria_iterator.h>
#include <grid/tria_boundary_lib.h>
#include <dofs/dof_handler.h>
#include <dofs/dof_accessor.h>
#include <dofs/dof_constraints.h>
#include <dofs/dof_tools.h>
#include <fe/fe_q.h>
#include <fe/fe_values.h>
#include <numerics/vectors.h>
#include <numerics/matrices.h>
#include <numerics/data_out.h>

                                 // This is probably the only new one - it
                                 // declares the class that we use to transfer
                                 // a solution on one grid to the one we
                                 // obtain after refining/coarsening it:
#include <numerics/solution_transfer.h>

                                 // And here comes the usual assortment of C++
                                 // header files:
#include <fstream>
#include <iostream>

#ifdef HAVE_STD_STRINGSTREAM
#  include <sstream>
#else
#  include <strstream>
#endif


                                 // The first thing we have here is a helper
                                 // function that computes an even power |v|^n
                                 // of a vector ``v'', by evaluating
                                 // (v*v)^(n/2). We need this in the
                                 // computations below where we do not want to
                                 // dwell on the fact that the gradient of the
                                 // solution is actually a scalar in the 1d
                                 // situation we consider in this program (in
                                 // 1d, the gradient is a vector with a single
                                 // element, which is easily extracted). Small
                                 // tricks like this make it significantly
                                 // simpler to later extend a program so that
                                 // it also runs in higher space dimensions.
                                 //
                                 // While the implementation of the function
                                 // is obvious, note the assertion at the
                                 // beginning of the function body, which
                                 // makes sure that the exponent is indeed an
                                 // even number (here, we use that ``n/2'' is
                                 // computed in integer arithmetic, i.e. any
                                 // remainder of the division is
                                 // lost). ``ExcMessage'' is a pre-defined
                                 // exception class that takes a string
                                 // argument explaining what goes wrong. It is
                                 // a simpler way to declare exceptions than
                                 // the ones shown in step-9 and step-13/14
                                 // where we explicitly declared exception
                                 // classes. However, by using a generic
                                 // exception class, we lose the ability to
                                 // attach additional information at run-time
                                 // to the exception message, such as the
                                 // value of the variable ``n''. By following
                                 // the way explained in above example
                                 // programs, adding this feature is simple,
                                 // though.
template <int dim>
inline
double gradient_power (const Tensor<1,dim> &v,
                       const unsigned int n)
{
  Assert ((n/2)*2 == n, ExcMessage ("Value of 'n' must be even"));
  double p = 1;
  for (unsigned int k=0; k<n; k+=2)
    p *= (v*v);
  return p;
}



                                 // Secondly, we declare a class that defines
                                 // our initial values for the nonlinear
                                 // iteration. It is a function object,
                                 // i.e. it has a member operator that returns
                                 // for a given point the value of the
                                 // function. The value we return is a random
                                 // perturbation of the ``x^1/3'' function
                                 // which we know is the optimal solution in a
                                 // larger function space. To make things a
                                 // little simpler on the optimizer, we return
                                 // zero if the proposed random value is
                                 // negative.
                                 //
                                 // Note that this class works strictly only
                                 // for 1d. If the program is to be extended
                                 // to higher space dimensions, so has to be
                                 // this class.
class InitializationValues : public Function<1> 
{
  public:
    InitializationValues () : Function<1>() {};
    
    virtual double value (const Point<1>     &p,
			  const unsigned int  component = 0) const;
};



                                 // So here comes the function that implements
                                 // the function object. The ``base'' value is
                                 // ``x^1/3'', while ``random'' is a random
                                 // number between -1 and 1 (note that
                                 // ``rand()'' returns a random integer value
                                 // between zero and ``RAND_MAX''; to convert
                                 // it to a floating point value between 0 and
                                 // 2, we have to divide by ``RAND_MAX'' and
                                 // multiply by two -- note that the first
                                 // multiplication has to happen in floating
                                 // point arithmetic, so that the division is
                                 // done in non-truncating floating point mode
                                 // as well; the final step is then to shift
                                 // the interval [0,2] to [-1,1]).
                                 //
                                 // In a second step, we add the base value
                                 // and a random value in [-0.1,0.1] together
                                 // and return it, unless it is less than
                                 // zero, in which case we take zero.
double InitializationValues::value (const Point<1> &p,
                                    const unsigned int) const 
{
  const double base = std::pow(p(0), 1./3.);
  const double random = 2.*rand()/RAND_MAX-1;
  return std::max (base+.1*random, 0.);
}



                                 // Next is the declaration of the main
                                 // class. As in most of the previous example
                                 // programs, the public interface of the
                                 // class consists only of a constructor and a
                                 // ``run'' function that does the actual
                                 // work. The constructor takes an additional
                                 // argument that indicates the number of the
                                 // run we are presently performing. This
                                 // value is only used at the very end when we
                                 // generate graphical output with a filename
                                 // that matches this number.
                                 //
                                 // The private section of the class has the
                                 // usual assortment of functions setting up
                                 // the computations, doing one nonlinear
                                 // step, refineming the mesh, doing a line
                                 // search for step length computations,
                                 // etc. The ``energy'' function computes the
                                 // value of the optimization functional on an
                                 // arbitrary finite element function with
                                 // nodal values given on the ``DoFHandler''
                                 // given as an argument. Since it does not
                                 // depend on the state of this object, we
                                 // declare this function as ``static''.
                                 //
                                 // The member variables of this class are
                                 // what we have seen before, and the
                                 // variables that characterize the linear
                                 // system to be solved in the next nonlinear
                                 // step, as well as the present approximation
                                 // of the solution.
template <int dim>
class MinimizationProblem 
{
  public:
    MinimizationProblem  (const unsigned int run_number);
    void run ();
    
  private:
    void initialize_solution ();
    void setup_system_on_mesh ();
    void assemble_step ();
    double line_search (const Vector<double> & update) const;
    void do_step ();
    void output_results () const;
    void refine_grid ();

    static double energy (const DoFHandler<dim> &dof_handler,
                          const Vector<double>  &function);
 

    const unsigned int run_number;
    
    Triangulation<dim>   triangulation;

    FE_Q<dim>            fe;
    DoFHandler<dim>      dof_handler;

    ConstraintMatrix     hanging_node_constraints;
    
    SparsityPattern      sparsity_pattern;
    SparseMatrix<double> matrix;

    Vector<double>       present_solution;
    Vector<double>       residual;
};



                                 // The constructor of this class is actually
                                 // somewhat boring:
template <int dim>
MinimizationProblem<dim>::MinimizationProblem (const unsigned int run_number)
                :
                run_number (run_number),
                fe (1),
		dof_handler (triangulation)
{}


                                 // Then, here is the function that
                                 // initializes the solution before the first
                                 // non-linear iteration, by setting the
                                 // initial values to the random function
                                 // described above and making sure that the
                                 // boundary values are set correctly. We will
                                 // then only seek updates to this function
                                 // with zero boundary values, so that the
                                 // boundary values are always correct.
                                 //
                                 // Note how we have specialized this function
                                 // to 1d only. We do this since the second
                                 // part of the function, where we deal with
                                 // boundary values, is only correct if we are
                                 // in 1d. Not generating a general template
                                 // for this function prevents the compiler
                                 // from erroneously compiling this function
                                 // for other space dimensions, then.
template <>
void MinimizationProblem<1>::initialize_solution () 
{
                                   // The first part is to assign the correct
                                   // size to the vector, and use library
                                   // function that takes a function object,
                                   // and interpolates the given vector living
                                   // on a ``DoFHandler'' to this function
                                   // object:
  present_solution.reinit (dof_handler.n_dofs());
  VectorTools::interpolate (dof_handler,
                            InitializationValues(),
                            present_solution);

                                   // Then we still have to make sure that we
                                   // get the boundary values right. This
                                   // could have been done inside the
                                   // ``InitializationValues'' class, but it
                                   // is instructive to see how it can also be
                                   // done, in particular since it is so
                                   // simple in 1d. First, start out with an
                                   // arbitrary cell on level 0, i.e. the
                                   // coarse mesh:
  DoFHandler<1>::cell_iterator cell;
  cell = dof_handler.begin(0);
                                   // Then move as far to the left as
                                   // possible. Note that while in two or more
                                   // space dimensions, there is is no
                                   // guarantee as to the coordinate
                                   // directions of a given face number of a
                                   // cell, in 1d the zeroth face (and
                                   // neighbor) is always the one to the left,
                                   // and the first one the one to the
                                   // right. Similarly, the zeroth child is
                                   // the left one, the first child is the
                                   // right one.
  while (cell->at_boundary(0) == false)
    cell = cell->neighbor(0);
                                   // Now that we are at the leftmost coarse
                                   // grid cell, go recursively through its
                                   // left children until we find a terminal
                                   // one:
  while (cell->has_children() == true)
    cell = cell->child(0);
                                   // Then set the value of the solution
                                   // corresponding to the zeroth degree of
                                   // freedom and the zeroth vertex of the
                                   // cell to zero. Note that the zeroth
                                   // vertex is the left one, and that zero is
                                   // the only valid second argument to the
                                   // call to ``vertex_dof_index'', since we
                                   // have a scalar finite element; thus,
                                   // there is only a single component.
  present_solution(cell->vertex_dof_index(0,0)) = 0;

                                   // Now do all the same with the right
                                   // boundary value, and set it to one:
  cell = dof_handler.begin(0);
  while (cell->at_boundary(1) == false)
    cell = cell->neighbor(1);
  while (cell->has_children())
    cell = cell->child(1);
  present_solution(cell->vertex_dof_index(1,0)) = 1;
}


                                 // The function that prepares the member
                                 // variables of this class for assembling the
                                 // linear system in each nonlinear step is
                                 // also not very interesting. This has all
                                 // been shown before in previous example
                                 // programs. Note, however, that all this
                                 // works in 1d just as in any other space
                                 // dimension, and would not require any
                                 // changes if we were to use the program in
                                 // another space dimension.
                                 //
                                 // Note that this function is only called
                                 // when the mesh has been changed (or before
                                 // the first nonlinear step). It only
                                 // initializes the variables to their right
                                 // sizes, but since these sizes don't change
                                 // as long as we don't change the mesh, we
                                 // can use them for more than just one
                                 // nonlinear iteration without reinitializing
                                 // them.
template <int dim>
void MinimizationProblem<dim>::setup_system_on_mesh ()
{
  hanging_node_constraints.clear ();
  DoFTools::make_hanging_node_constraints (dof_handler,
					   hanging_node_constraints);
  hanging_node_constraints.close ();

  sparsity_pattern.reinit (dof_handler.n_dofs(),
			   dof_handler.n_dofs(),
			   dof_handler.max_couplings_between_dofs());
  DoFTools::make_sparsity_pattern (dof_handler, sparsity_pattern);

  hanging_node_constraints.condense (sparsity_pattern);

  sparsity_pattern.compress();
}



                                 // Next is the function that assembles the
                                 // linear system. The first part,
                                 // initializing various local variables is
                                 // what we have been doing previously
                                 // already.
template <int dim>
void MinimizationProblem<dim>::assemble_step ()
{
                                   // The first two lines of the function
                                   // clear the matrix and right hand side
                                   // values of their prior content, which
                                   // could possibly still be there from the
                                   // previous nonlinear step.
  matrix.reinit (sparsity_pattern);
  residual.reinit (dof_handler.n_dofs());

                                   // Then we initialize a ``FEValues'' object
                                   // with a 4-point Gauss quadrature
                                   // formula. This object will be used to
                                   // compute the values and gradients of the
                                   // shape functions at the quadrature
                                   // points, which we need to assemble the
                                   // matrix and right hand side of the
                                   // nonlinear step as outlined in the
                                   // introduction to this example program. In
                                   // order to compute values and gradients,
                                   // we need to pass the ``update_values''
                                   // and ``update_gradients'' flags to the
                                   // constructor, and the
                                   // ``update_JxW_values'' flag for the
                                   // Jacobian times the weight at a
                                   // quadrature point. In addition, we need
                                   // to have the coordinate values of each
                                   // quadrature point in real space for the
                                   // ``x-u^3'' terms; to get these from the
                                   // ``FEValues'' object, we need to pass it
                                   // the ``update_q_points'' flag.
                                   //
                                   // It is a simple calculation to figure out
                                   // that for linear elements, the integrals
                                   // in the right hand side semilinear form
                                   // is a polynomial of sixth order. Thus,
                                   // the appropriate quadrature formula is
                                   // the one we have chosen here.
  QGauss4<dim>  quadrature_formula;
  FEValues<dim> fe_values (fe, quadrature_formula, 
			   UpdateFlags(update_values    |
				       update_gradients |
				       update_q_points  |
				       update_JxW_values));

                                   // Next, here are the usual two convenience
                                   // variables, followed by declarations for
                                   // the local contributions to matrix and
                                   // right hand side, as well as an array to
                                   // hold the indices of the local degrees of
                                   // freedom on each cell:
  const unsigned int dofs_per_cell = fe.dofs_per_cell;
  const unsigned int n_q_points    = quadrature_formula.n_quadrature_points;

  FullMatrix<double>   cell_matrix (dofs_per_cell, dofs_per_cell);
  Vector<double>       cell_rhs (dofs_per_cell);

  std::vector<unsigned int> local_dof_indices (dofs_per_cell);

                                   // The next two variables are needed since
                                   // the problem we consider is nonlinear,
                                   // and thus the right hand side depends on
                                   // the previous solution (in a Newton
                                   // method, for example, the left hand side
                                   // matrix would also depend on the previous
                                   // solution, but as explained in the
                                   // introduction, we only use a simple
                                   // gradient-type method in which the matrix
                                   // is a scaled Laplace-type matrix). In
                                   // order to compute the values of the
                                   // integrand for the right hand side, we
                                   // therefore need to have the values and
                                   // gradients of the previous solution at
                                   // the quadrature points. We will get them
                                   // from the ``FEValues'' object above, and
                                   // will put them into the following two
                                   // variables:
  std::vector<double>         local_solution_values (n_q_points);
  std::vector<Tensor<1,dim> > local_solution_grads (n_q_points);

                                   // Now, here comes the main loop over all
                                   // the cells of the mesh:
  typename DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();
  for (; cell!=endc; ++cell)
    {
                                       // First, clear the objects that hold
                                       // the local matrix and right hand side
                                       // contributions for this cell:
      cell_matrix = 0;
      cell_rhs = 0;

                                       // Then initialize the values and
                                       // gradients of the shape functions at
                                       // the quadrature points of this cell:
      fe_values.reinit (cell);

                                       // And get the values and gradients of
                                       // the previous solution at the
                                       // quadrature points. To get them, we
                                       // don't actually have to do much,
                                       // except for giving the ``FEValues''
                                       // object the global node vector from
                                       // which to compute this data, and a
                                       // reference to the objects into which
                                       // to put them. After the calls, the
                                       // ``local_solution_values'' and
                                       // ``local_solution_values'' variables
                                       // will contain values and gradients
                                       // for each of the quadrature points on
                                       // this cell.
      fe_values.get_function_values (present_solution,
                                     local_solution_values);
      fe_values.get_function_grads (present_solution,
                                    local_solution_grads);

                                       // Then loop over all quadrature
                                       // points:
      for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
        {
                                           // Have convenience variables for
                                           // the values and gradient of the
                                           // solution at the present
                                           // quadrature point, as well as the
                                           // location in real space of this
                                           // quadrature point, and of the
                                           // expression (x-u^3), since it
                                           // appears so often:
          const double u = local_solution_values[q_point],
                       x = fe_values.quadrature_point(q_point)(0);
          const double x_minus_u3 = (x-std::pow(u,3));
          const Tensor<1,dim> u_prime = local_solution_grads[q_point];

                                           // Then do the double loop over all
                                           // shape functions to compute the
                                           // local contribution to the
                                           // matrix. The terms are simple
                                           // equivalents of the formula
                                           // stated in the introduction. Note
                                           // how we extract the size of an
                                           // element from the iterator to the
                                           // present cell:
          for (unsigned int i=0; i<dofs_per_cell; ++i)
            for (unsigned int j=0; j<dofs_per_cell; ++j)
              cell_matrix(i,j)
                += (fe_values.shape_grad(i,q_point) *
                    fe_values.shape_grad(j,q_point) *
                    cell->diameter() *
                    cell->diameter()
                    +
                    fe_values.shape_value(i,q_point) *
                    fe_values.shape_value(j,q_point)) *
                fe_values.JxW(q_point);

                                           // And here comes the loop over all
                                           // local degrees of freedom to form
                                           // the right hand side. The formula
                                           // looks a little convoluted, but
                                           // is again a simple image of what
                                           // was given in the introduction:
          for (unsigned int i=0; i<dofs_per_cell; ++i)
            cell_rhs(i) += -((6. * x_minus_u3 *
                              gradient_power (u_prime, 4) *
                              fe_values.shape_value(i,q_point)
                              *
                              (x_minus_u3 *
                               (u_prime * 
                                fe_values.shape_grad(i,q_point))
                               -
                               (u_prime*u_prime) * u * u *
                               fe_values.shape_value(i,q_point))
                              )
                             *
                             fe_values.JxW(q_point));
        }
      
                                       // After summing up all the
                                       // contributions, we have to transfer
                                       // them to the global objects. This is
                                       // done in the same way as always
                                       // before:
      cell->get_dof_indices (local_dof_indices);
      for (unsigned int i=0; i<dofs_per_cell; ++i)
	{
	  for (unsigned int j=0; j<dofs_per_cell; ++j)
	    matrix.add (local_dof_indices[i],
			       local_dof_indices[j],
			       cell_matrix(i,j));
	  
	  residual(local_dof_indices[i]) += cell_rhs(i);
	}
    }

                                   // Now that we have all the local
                                   // contributions summed up, we have to
                                   // eliminate hanging node constraints and
                                   // boundary values. Hanging nodes are
                                   // simple:
  hanging_node_constraints.condense (matrix);
  hanging_node_constraints.condense (residual);

                                   // Boundary values are, too, but with a
                                   // twist this time: in all previous example
                                   // programs, we have used that by default
                                   // (i.e. unless something else is set), all
                                   // boundaries have indicator zero. To
                                   // figure out what boundary indicator a
                                   // face of a cell had, the library
                                   // functions would query an iterator
                                   // designating this face, which would in
                                   // turn pluck out this value from some of
                                   // the data structures in the
                                   // library. Unfortunately, in 1d cells have
                                   // no faces: these would only be points,
                                   // and we don't associated anything in the
                                   // library with points except for their
                                   // coordinates. Thus there are no face
                                   // iterators, and no way to figure out
                                   // which boundary indicator it may have. On
                                   // the other hand, in 1d, there can only be
                                   // two boundaries anyway for a connected
                                   // domain: the left end point and the right
                                   // end point. And in contrast to the case
                                   // in higher dimensions, where the
                                   // (changeable) default is zero for all
                                   // boundary parts, in 1d the convention is
                                   // that the left boundary point has
                                   // indicator zero, while the right boundary
                                   // point has indicator one. Since there are
                                   // no face iterators, it is also not
                                   // possible to change this, but you will
                                   // hardly ever have to. So in order to
                                   // assign zero boundary values on both
                                   // sides, in 1d we not only need to
                                   // evaluate boundary values for indicator
                                   // zero, but also for indicator one. If
                                   // this program is ever going to be run in
                                   // higher dimensions, then we should only
                                   // evaluate for indicator zero, which is
                                   // why we have placed the ``if'' statement
                                   // in front of the second function call.
                                   //
                                   // Note that we need zero boundary
                                   // conditions on both ends, since the space
                                   // in which search for the solution has
                                   // fixed boundary conditions zero and one,
                                   // and we have set the initial values to
                                   // already satisfy them. Thus, the updates
                                   // computed in each nonlinear step must
                                   // have zero boundary values.
  std::map<unsigned int,double> boundary_values;
  VectorTools::interpolate_boundary_values (dof_handler,
					    0,
					    ZeroFunction<dim>(),
					    boundary_values);
  if (dim == 1)
    VectorTools::interpolate_boundary_values (dof_handler,
                                              1,
                                              ZeroFunction<dim>(),
                                              boundary_values);
  Vector<double> dummy (residual.size());
  MatrixTools::apply_boundary_values (boundary_values,
				      matrix,
				      dummy,
				      residual);
}


                                 // Once we have a search (update) direction,
                                 // we need to figure out how far to go in
                                 // this direction. This is what line search
                                 // is good for, and this function does
                                 // exactly this: compute and return the
                                 // length of the update step.
                                 //
                                 // Since we already know the direction, we
                                 // only have to solve the one-dimensional
                                 // problem of minimizing the energy along
                                 // this direction. Note, however, that in
                                 // general we do not have the gradient of the
                                 // energy functional in this direction, so we
                                 // have to approximate it (and the second
                                 // derivatives) using finite differences.
                                 //
                                 // In most applications, it is sufficient to
                                 // find an approximate minimizer of this
                                 // one-dimensional problem, or even just a
                                 // point that may not be a minimizer but
                                 // instead just satisfies a few conditions
                                 // like those of Armijo and Goldstein. The
                                 // rational for this is generally that
                                 // evaluating the objective function too
                                 // often is too expensive. However, here, we
                                 // are a little more lenient, since the
                                 // overall run-time is dominated by inverting
                                 // the system matrix in each nonlinear
                                 // step. Thus, we will do this minimization
                                 // by using a fixed number of five Newton
                                 // steps in this one-dimensional problem, and
                                 // using a bisection algorithm as a substep
                                 // in it.
                                 //
                                 // As is quite common in step length
                                 // procedures, this function contains a fair
                                 // number of heuristics and strategies that
                                 // might not be obvious at first. Step length
                                 // determination is notorious for its
                                 // complications, and this implementation is
                                 // not an exception. Note that if one tries
                                 // to omit the special-casing, then one
                                 // oftentimes encounters situations where the
                                 // found step length is really not very good.
template <int dim>
double
MinimizationProblem<dim>::line_search (const Vector<double> &update) const
{
                                   // Start out with a zero step length:
  double alpha = 0.;
  Vector<double> tmp (present_solution.size());

                                   // Then do at most five Newton steps:
  for (unsigned int step=0; step<5; ++step)
    {
                                       // At the present location, which is
                                       // ``present_solution+alpha*update'',
                                       // evaluate the energy
      tmp = present_solution;
      tmp.add (alpha, update);
      const double f_a = energy (dof_handler, tmp);

                                       // Then determine a finite difference
                                       // step length ``dalpha'', and also
                                       // evaluate the energy functional at
                                       // positions ``alpha+dalpha'' and
                                       // ``alpha-dalpha'' along the search
                                       // direction:
      const double dalpha = (alpha != 0 ? alpha/100 : 0.01);
      
      tmp = present_solution;
      tmp.add (alpha+dalpha, update);
      const double f_a_plus = energy (dof_handler, tmp);

      tmp = present_solution;
      tmp.add (alpha-dalpha, update);
      const double f_a_minus = energy (dof_handler, tmp);

                                       // From these three data points, we can
                                       // compute a finite difference
                                       // approximation of the first and
                                       // second derivatives:
      const double f_a_prime       = (f_a_plus-f_a_minus) / (2*dalpha);
      const double f_a_doubleprime = ((f_a_plus-2*f_a+f_a_minus) /
                                      (dalpha*dalpha));

                                       // If the gradient is (relative to the
                                       // energy value) too small, then this
                                       // means that we have found a minimum
                                       // of the energy functional along the
                                       // search direction. In this case,
                                       // abort here and return the found step
                                       // length value:
      if (std::fabs(f_a_prime) < 1e-7*std::fabs(f_a))
        break;

                                       // Alternatively, also abort if the
                                       // curvature is too small, because we
                                       // can't compute a Newton step
                                       // then. This is somewhat
                                       // unsatisfactory, since we are not at
                                       // a minimum, and can certainly be
                                       // improved. There are a number of
                                       // other strategies for this case,
                                       // which we leave for interested
                                       // readers:
      if (std::fabs(f_a_doubleprime) < 1e-7*std::fabs(f_a_prime))
        break;

                                       // Then compute the Newton step as the
                                       // negative of the inverse Hessian
                                       // applied to the gradient.
      double step_length = -f_a_prime / f_a_doubleprime;

                                       // And do a number of correcting steps:
                                       // if the energy at the predicted new
                                       // position would be larger than at the
                                       // present position, then halve the
                                       // step length and try again. If this
                                       // does not help after three such
                                       // cycles, then simply give up and use
                                       // the value we have.
      for (unsigned int i=0; i<3; ++i)
        {
          tmp = present_solution;
          tmp.add (alpha+step_length, update);
          const double e = energy (dof_handler, tmp);
          
          if (e >= f_a)
            step_length /= 2;
          else
            break;
        }

                                       // After all this, update alpha and go
                                       // on to the next Newton step.
      alpha += step_length;
    }

                                   // Finally, return with the computed step length.
  return alpha;
}



                                 // The next function is again a rather boring
                                 // one: it does one nonlinear step, by
                                 // calling the function that assembles the
                                 // linear system, then solving it, computing
                                 // a step length, and finally updating the
                                 // solution vector. This should all be mostly
                                 // self-explanatory, given that we have shown
                                 // the solution of a linear system before.
template <int dim>
void MinimizationProblem<dim>::do_step ()
{          
  assemble_step ();

  Vector<double> update (present_solution.size());
  {
    SolverControl           solver_control (residual.size(),
                                            1e-2*residual.l2_norm());
    PrimitiveVectorMemory<> vector_memory;
    SolverCG<>              solver (solver_control, vector_memory);
    
    PreconditionSSOR<> preconditioner;
    preconditioner.initialize(matrix);

    solver.solve (matrix, update, residual,
                  preconditioner);
    hanging_node_constraints.distribute (update);
  }

  const double step_length = line_search (update);
  present_solution.add (step_length, update);
}



                                 // The same holds for the function that
                                 // outputs the solution in gnuplot format
                                 // into a file with a name that includes the
                                 // number of the run we are presently
                                 // performing.
template <int dim>
void
MinimizationProblem<dim>::output_results () const
{
  DataOut<dim> data_out;
  data_out.attach_dof_handler (dof_handler);
  data_out.add_data_vector (present_solution, "solution");
  data_out.build_patches ();

#ifdef HAVE_STD_STRINGSTREAM
  std::ostringstream filename;
#else
  std::ostrstream filename;
#endif
  filename << "solution-"
           << run_number
           << ".gnuplot"
           << std::ends;
#ifdef HAVE_STD_STRINGSTREAM
  std::ofstream out (filename.str().c_str());
#else
  std::ofstream out (filename.str());
#endif

  data_out.write_gnuplot (out);
}



                                 // The function to compute error indicator
                                 // and refine the mesh accordingly is a
                                 // little more interesting. In particular, it
                                 // shows some more of the techniques usually
                                 // used in 1d applications. First, note that
                                 // this again is a specialization that only
                                 // works in 1d. However, to make later
                                 // extension to higher space dimensions
                                 // simpler, we define a constant integer
                                 // ``dim'' at the beginning of the function;
                                 // by using this constant as template
                                 // argument in all places, we are actually
                                 // able to write most of the code as if it
                                 // were dimension independent, thus
                                 // minimizing the amount of later changes.
template <>
void MinimizationProblem<1>::refine_grid ()
{
  const unsigned int dim = 1;
  
  Vector<float> error_indicators (triangulation.n_active_cells());

                                   // Then define the quadrature formula, and
                                   // what values we will want to extract from
                                   // the solution. Here, we use the two-point
                                   // trapezoidal rule, i.e. we evaluate the
                                   // residual only at the end points of the
                                   // cells. Incidentally, this also makes
                                   // evaluating the jump terms between cells
                                   // simpler. Note that for the error
                                   // indicators, we not only need values and
                                   // gradients of the solution, but also its
                                   // second derivatives, as well as the
                                   // physical location of quadrature points.
  QTrapez<dim> quadrature;
  FEValues<dim> fe_values (fe, quadrature,
                           update_values   | update_gradients |
                           update_second_derivatives |
                           update_q_points | update_JxW_values);

                                   // The error indicator formula presented in
                                   // the introduction requires us to compute
                                   // jumps of the solution and gradient
                                   // across cell boundaries. Since the
                                   // solution itself is continuous, we only
                                   // need to evaluate the gradient on the
                                   // neighbor cells. To avoid some of the
                                   // work needed to reinitialize a
                                   // ``FEValues'' object on a cell, we define
                                   // another such object here that we will
                                   // only use for the neighbor cells. The
                                   // data we need from the side of the
                                   // present cell is provided by above
                                   // object.
  FEValues<dim> neighbor_fe_values (fe, quadrature,
                                    update_gradients);

                                   // Then, as before, we need objects holding
                                   // values and derivatives of the solution
                                   // at quadrature points. Here, we also need
                                   // second derivatives, which is simple,
                                   // however:
  std::vector<double> local_values (quadrature.n_quadrature_points);
  std::vector<Tensor<1,dim> > local_gradients (quadrature.n_quadrature_points);
  std::vector<Tensor<2,dim> > local_2nd_derivs (quadrature.n_quadrature_points);

                                   // With all this, we can start the loop
                                   // over all cells. Since we need to write
                                   // the result for each cell into
                                   // consecutive elements of a vector, we
                                   // also keep a running index ``cell_index''
                                   // that we increase with each cell treated.
  DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active (),
    endc = dof_handler.end ();
  for (unsigned int cell_index = 0; cell!=endc; ++cell, ++cell_index)
    {
                                       // After initializing the ``FEValues''
                                       // object on each cell, use it to
                                       // evaluate solution and first and
                                       // second derivatives of it at the
                                       // quadrature points:
      fe_values.reinit (cell);
      fe_values.get_function_values (present_solution, local_values);
      fe_values.get_function_grads (present_solution, local_gradients);
      fe_values.get_function_2nd_derivatives (present_solution, local_2nd_derivs);

                                       // Given the formula in the
                                       // introduction, the computation of the
                                       // cell residuals should actually be
                                       // relatively obvious. The result,
                                       // multiplied by the appropriate power
                                       // of the cell's size is then written
                                       // into the vector of error indicators.
                                       //
                                       // Note that in the following
                                       // computations, we have already made
                                       // use of the fact that we are in 1d,
                                       // since we extract the gradient as a
                                       // scalar value.
      double cell_residual_norm = 0;
      for (unsigned int q=0; q<quadrature.n_quadrature_points; ++q)
        {
          const double x             = fe_values.quadrature_point(q)[0];
          const double u             = local_values[q];
          const double u_prime       = local_gradients[q][0];
          const double u_doubleprime = local_2nd_derivs[q][0][0];
          const double local_residual_value
            = ((x-u*u*u) * std::pow(u_prime, 4) *
               (u*u*u_prime*u_prime
                +
                5*(x-u*u*u)*u_doubleprime
                +
                2*u_prime*(1-3*u*u*u_prime)));
          
          cell_residual_norm += (local_residual_value * local_residual_value *
                                 fe_values.JxW(q));
        }
      error_indicators(cell_index) = cell_residual_norm *
                                     cell->diameter() * cell->diameter();

                                       // The next step is to evaluate the
                                       // jump terms. To make computations
                                       // somewhat simpler (and to free up the
                                       // ``local_*'' variables for use on
                                       // neighboring elements), we define
                                       // some convenience variables for the
                                       // positions of the left and right cell
                                       // boundary point, as well as the
                                       // values and gradients at these
                                       // points.
                                       //
                                       // To be cautious, we don't blindly
                                       // trust that the trapezoidal rule has
                                       // its evaluation points as the left
                                       // and right end point of the cell (it
                                       // could in principle have them in the
                                       // reverse order, i.e. the zeroth point
                                       // is at x=1, and the first one at
                                       // x=0), and use an assertion to
                                       // actually check for this. If this
                                       // would not be the case, an exception
                                       // of the (predefined) class
                                       // ``ExcInternalError'' would be
                                       // thrown. Of course, this does not
                                       // happen in this program, but it shows
                                       // a way of defensive coding: if you
                                       // are not sure of an assumption, guard
                                       // it by a test. This also guards us
                                       // against possible future changes in
                                       // the library: the quadrature classes
                                       // do not promise any particular order
                                       // of their quadrature points, so the
                                       // ``QTrapez'' class could in principle
                                       // change the order of its two
                                       // evaluation points. In that case,
                                       // your code would tell you that
                                       // something changed, rather than
                                       // computing a wrong result when you
                                       // upgrade to a new version of the
                                       // library. (The point made here is
                                       // theoretical: we are not going to
                                       // change the order of evaluation
                                       // points; the intent is simply how to
                                       // add some defensive touches to a
                                       // program that make sure that it
                                       // really does what it is hoped to do.)
                                       //
                                       // Given that we are now sure that
                                       // ``x_left'' and ``x_right'',
                                       // extracted from the zeroth and first
                                       // quadrature point, are indeed the
                                       // left and right vertex of the cell,
                                       // we can also be sure that the values
                                       // we extract for ``u_left'' et al. are
                                       // the ones we expect them to be, since
                                       // the order of these values must of
                                       // course match the order of the
                                       // quadrature points.
      const double x_left  = fe_values.quadrature_point(0)[0];
      const double x_right = fe_values.quadrature_point(1)[0];

      Assert (x_left  == cell->vertex(0)[0], ExcInternalError());
      Assert (x_right == cell->vertex(1)[0], ExcInternalError());

      const double u_left  = local_values[0];
      const double u_right = local_values[1];

      const double u_prime_left  = local_gradients[0][0];
      const double u_prime_right = local_gradients[1][0];

                                       // Next, we have to check whether this
                                       // cell has a left neighbor:
      if (cell->at_boundary(0) == false)
        {
                                           // If so, find its left
                                           // neighbor. We do so by asking for
                                           // the cell that is immediately
                                           // adjacent to the left (the zeroth
                                           // neighbor in 1d). However, this
                                           // may be a cell that in itself has
                                           // children, so to get to the
                                           // active left neighbor, we have to
                                           // recursively check whether that
                                           // cell has children, and if so
                                           // take its right child, since that
                                           // is adjacent to the left of the
                                           // present cell. Note that unless
                                           // you are in 1d, there is no safe
                                           // way to assume that the first
                                           // child of the zeroth neighbor is
                                           // indeed adjacent to the present
                                           // cell. Rather, more than one of
                                           // the children of a neighbor may
                                           // be adjacent to the present
                                           // cell. Also note that in two or
                                           // higher space dimensions, a
                                           // neighbor of an active cell may
                                           // only be at most once refined,
                                           // since we have the rule that
                                           // there can only be one hanging
                                           // node per face. This rule does
                                           // not exist in 1d: neighboring
                                           // cells may have totally
                                           // independent refinement
                                           // levels. Thus, we really need the
                                           // ``while'' loop, not only an
                                           // ``if'' clause.
          DoFHandler<dim>::cell_iterator left_neighbor = cell->neighbor(0);
          while (left_neighbor->has_children())
            left_neighbor = left_neighbor->child(1);

                                           // With the so-found neighbor,
                                           // initialize the second
                                           // ``FEValues'' object to it,
                                           // extract the gradients of the
                                           // solution there, and from this
                                           // get the gradient at the
                                           // interface (this is the first
                                           // element of ``local_gradients'',
                                           // since the right end point of the
                                           // neighbor cell has index 1) as a
                                           // scalar value (this is the zeroth
                                           // component of
                                           // ``local_gradients[1]''.
          neighbor_fe_values.reinit (left_neighbor);
          neighbor_fe_values.get_function_grads (present_solution, local_gradients);

          const double neighbor_u_prime_left = local_gradients[1][0];

                                           // Then compute the jump, and add a
                                           // suitable multiple to the error
                                           // indicator for this cell:
          const double left_jump = std::pow(x_left-std::pow(u_left,3), 2) *
                                   (std::pow(neighbor_u_prime_left,5) -
                                    std::pow(u_prime_left,5));
          error_indicators(cell_index) += left_jump * left_jump *
                                          cell->diameter();
        }

                                       // Once we have done the left neighbor,
                                       // we can play exactly the same game
                                       // with the right neighbor:
      if (cell->at_boundary(1) == false)
        {
          DoFHandler<dim>::cell_iterator right_neighbor = cell->neighbor(1);
          while (right_neighbor->has_children())
            right_neighbor = right_neighbor->child(0);
          
          neighbor_fe_values.reinit (right_neighbor);
          neighbor_fe_values.get_function_grads (present_solution, local_gradients);

          const double neighbor_u_prime_right = local_gradients[0][0];

          const double right_jump = std::pow(x_right-std::pow(u_right,3), 2) *
                                   (std::pow(neighbor_u_prime_right,5) -
                                    std::pow(u_prime_right,5));
          error_indicators(cell_index) += right_jump * right_jump *
                                          cell->diameter();
        }      
    } 

                                   // Now we have all the refinement
                                   // indicators computed, and want to refine
                                   // the grid. In contrast to previous
                                   // examples, however, we would like to
                                   // transfer the solution vector from the
                                   // old to the new grid. This is what the
                                   // ``SolutionTransfer'' class is good for,
                                   // but it requires some preliminary
                                   // work. First, we need to tag the cells
                                   // that we want to refine or coarsen, as
                                   // usual:
  GridRefinement::refine_and_coarsen_fixed_number (triangulation,
						   error_indicators,
						   0.3, 0.03);
                                   // Then, however, we need an additional
                                   // step: if, for example, you flag a cell
                                   // that is once more refined than its
                                   // neighbor, and that neighbor is not
                                   // flagged for refinement, we would end up
                                   // with a jump of two refinement levels
                                   // across a cell interface. In 1d, this
                                   // would in general be allowed, but not in
                                   // higher space dimensions, and some mesh
                                   // smoothing algorithms in 1d may also
                                   // disallow this. To avoid these
                                   // situations, the library will silently
                                   // also have to refine the neighbor cell
                                   // once. It does so by calling the
                                   // ``Triangulation<dim>::prepare_coarsening_and_refinement''
                                   // function before actually doing the
                                   // refinement and coarsening. This function
                                   // flags a set of additional cells for
                                   // refinement or coarsening, to enforce
                                   // rules like the one-hanging-node
                                   // rule. The cells that are flagged for
                                   // refinement and coarsening after calling
                                   // this function are exactly the ones that
                                   // will actually be refined or
                                   // coarsened. Since the
                                   // ``SolutionTransfer'' class needs this
                                   // information in order to store the data
                                   // from the old mesh and transfer to the
                                   // new one.
  triangulation.prepare_coarsening_and_refinement();

                                   // With this out of the way, we initialize
                                   // a ``SolutionTransfer'' object with the
                                   // present ``DoFHandler'' and attach the
                                   // solution vector to it:
  SolutionTransfer<dim,double> solution_transfer(dof_handler);
  solution_transfer.prepare_for_coarsening_and_refinement (present_solution);

                                   // Then we do the actual refinement, and
                                   // distribute degrees of freedom on the new
                                   // mesh:
  triangulation.execute_coarsening_and_refinement ();
  dof_handler.distribute_dofs (fe);

                                   // Finally, we retrieve the old solution
                                   // interpolated to the new mesh. Since the
                                   // ``SolutionTransfer'' function does not
                                   // actually store the values of the old
                                   // solution, but rather indices, we need to
                                   // preserve the old solution vector until
                                   // we have gotten the new interpolated
                                   // values. Thus, we have the new values
                                   // written into a temporary vector, and
                                   // only afterwards write them into the
                                   // solution vector object:
  Vector<double> tmp (dof_handler.n_dofs());
  solution_transfer.interpolate (present_solution, tmp);
  present_solution = tmp;

                                   // Here is some final thing, that is
                                   // actually unnecessary in 1d, but
                                   // necessary for higher space dimensions,
                                   // so we show it anyway: the result of what
                                   // the ``SolutionTransfer'' class provides
                                   // is a vector that is interpolated from
                                   // the old to the new mesh. Unfortunately,
                                   // it does not necessarily have the right
                                   // values at constrained (hanging) nodes,
                                   // so we have to fix this up to make the
                                   // solution conforming again. The simplest
                                   // way to do this is this:
  hanging_node_constraints.clear ();
  DoFTools::make_hanging_node_constraints (dof_handler,
					   hanging_node_constraints);
  hanging_node_constraints.close ();  
  hanging_node_constraints.distribute (present_solution);
                                   // This is wasteful, since we create a
                                   // ``ConstraintMatrix'' object that will be
                                   // recreated again in the next call to
                                   // ``setup_system_on_mesh'' immediately
                                   // afterwards. A more efficient
                                   // implementation would make sure that it
                                   // is created only once. We don't care so
                                   // much here, since in 1d there are no
                                   // constraints, so all of these operations
                                   // are really cheap, but we do not
                                   // recommend this as general programming
                                   // strategy.
}



                                 // Before going over to the framework
                                 // functions, we still need to look at the
                                 // implementation of the function that
                                 // computes the energy of a nodal vector in
                                 // the functional considered in this example
                                 // program. Its idea is simple: take a nodal
                                 // vector and the ``DoFHandler'' object it is
                                 // living on, then loop over all cells and
                                 // add up the local contributions to the
                                 // energy:
template <int dim>
double
MinimizationProblem<dim>::energy (const DoFHandler<dim> &dof_handler,
                                  const Vector<double>  &function)
{
                                   // First define the quadrature formula and
                                   // a ``FEValues'' object with which to
                                   // compute the values of the input function
                                   // at the quadrature points. Note again
                                   // that the integrand is a polynomial of
                                   // degree six, so a 4-point Gauss formula
                                   // is appropriate:
  QGauss4<dim>  quadrature_formula;
  FEValues<dim> fe_values (dof_handler.get_fe(), quadrature_formula, 
			   UpdateFlags(update_values    |
				       update_gradients |
				       update_q_points  |
				       update_JxW_values));

  const unsigned int   n_q_points    = quadrature_formula.n_quadrature_points;

                                   // Then, just as when we integrated the
                                   // linear system, we need two variables
                                   // that will hold the values and gradients
                                   // of the given function at the quadrature
                                   // points:
  std::vector<double>         local_solution_values (n_q_points);
  std::vector<Tensor<1,dim> > local_solution_grads (n_q_points);

                                   // With this, define an energy variable,
                                   // and loop over all the cells:
  double energy = 0.;

  typename DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();
  for (; cell!=endc; ++cell)
    {
                                       // On each cell, initialize the
                                       // ``FEValues'' object, and extract
                                       // values and gradients of the given
                                       // function:
      fe_values.reinit (cell);
      fe_values.get_function_values (function,
                                     local_solution_values);
      fe_values.get_function_grads (function,
                                    local_solution_grads);

                                       // Then loop over all quadrature points
                                       // on this cell, and add up the
                                       // contribution of each to the global
                                       // energy:
      for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
        energy += (std::pow (fe_values.quadrature_point(q_point)(0)
                             -
                             std::pow (local_solution_values[q_point], 3),
                             2) *
                   gradient_power (local_solution_grads[q_point], 6) *
                   fe_values.JxW (q_point));
    }

                                   // Once we have done this, return the
                                   // integrated value.
  return energy;
}


                                 // So here is the driver function,
                                 // ``run()''. It generate a coarse mesh,
                                 // refines it a couple of times, and
                                 // initializes the starting values. It then
                                 // goes into a loop in which we first set up
                                 // the member variables for the new mesh, and
                                 // then do a fixed number of five gradient
                                 // steps. If after this the energy has not
                                 // significantly decreased compares to the
                                 // last time we checked, we assume that we
                                 // have converged and exit, otherwise we
                                 // refine the mesh and start over. Once we
                                 // have determined that the computations have
                                 // converged somewhere, we output the
                                 // results.
template <int dim>
void MinimizationProblem<dim>::run () 
{
  GridGenerator::hyper_cube (triangulation, 0., 1.);
  triangulation.refine_global (4);
  dof_handler.distribute_dofs (fe);
  initialize_solution ();

  double last_energy = energy (dof_handler, present_solution);
  
  while (true)
    {
      setup_system_on_mesh ();
      
      for (unsigned int iteration=0; iteration<5; ++iteration)
        do_step ();

      const double this_energy = energy (dof_handler, present_solution);
      std::cout << "   Energy: " << this_energy << std::endl;

      if ((last_energy-this_energy) < 1e-5*last_energy)
        break;

      last_energy = this_energy;

      refine_grid ();
    }
  
  output_results ();
  
  std::cout << std::endl;
}


                                 // Finally: ``main()''. This function does
                                 // what all its counterparts in previous
                                 // examples already did, i.e. create an
                                 // object of the main class, and hand off
                                 // work to them, only retaining its role as a
                                 // guard to catch exceptions and print some
                                 // information if we get one. The only
                                 // difference is that it generates objects
                                 // multiple times and runs them. Since the
                                 // initial value for each run is a random
                                 // function, where the random number
                                 // generator returns different values each
                                 // time, all these runs are actually
                                 // different, although it may seem that they
                                 // are independent of each other.
int main () 
{
  try
    {
      deallog.depth_console (0);

      const unsigned int n_realizations = 10;
      for (unsigned int realization=0; realization<n_realizations; ++realization)
        {
          std::cout << "Realization " << realization << ":" << std::endl;
  
          MinimizationProblem<1> minimization_problem_1d (realization);
          minimization_problem_1d.run ();
        }
    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl << std::endl
		<< "----------------------------------------------------"
		<< std::endl;
      std::cerr << "Exception on processing: " << std::endl
		<< exc.what() << std::endl
		<< "Aborting!" << std::endl
		<< "----------------------------------------------------"
		<< std::endl;
      return 1;
    }
  catch (...) 
    {
      std::cerr << std::endl << std::endl
		<< "----------------------------------------------------"
		<< std::endl;
      std::cerr << "Unknown exception!" << std::endl
		<< "Aborting!" << std::endl
		<< "----------------------------------------------------"
		<< std::endl;
      return 1;
    }
  return 0;
}
