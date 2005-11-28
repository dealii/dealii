/* $Id$ */
/* Author: Wolfgang Bangerth, University of Heidelberg, 1999 */

/*    $Id$       */
/*    Version: $Name$                                          */
/*                                                                */
/*    Copyright (C) 1999, 2000, 2001, 2002, 2003, 2004, 2005 by the deal.II authors */
/*                                                                */
/*    This file is subject to QPL and may not be  distributed     */
/*    without copyright and license information. Please refer     */
/*    to the file deal.II/doc/license.html for the  text  and     */
/*    further information on this license.                        */


const unsigned int degree = 2;

				 // The first few (many?) include
				 // files have already been used in
				 // the previous example, so we will
				 // not explain their meaning here
				 // again.
#include <grid/tria.h>
#include <dofs/dof_handler.h>
#include <dofs/dof_renumbering.h>
#include <grid/grid_generator.h>
#include <grid/tria_accessor.h>
#include <grid/tria_iterator.h>
#include <dofs/dof_accessor.h>
#include <fe/fe_raviart_thomas.h>
#include <fe/fe_dgq.h>
#include <fe/fe_system.h>
#include <dofs/dof_tools.h>
#include <fe/fe_values.h>
#include <base/quadrature_lib.h>
#include <base/function.h>
#include <numerics/vectors.h>
#include <numerics/matrices.h>
#include <lac/block_vector.h>
#include <lac/full_matrix.h>
#include <lac/block_sparse_matrix.h>
#include <lac/solver_cg.h>
#include <lac/solver_gmres.h>
#include <lac/precondition.h>

#include <numerics/data_out.h>
#include <fstream>
#include <iostream>

				 // This is new, however: in the
				 // previous example we got some
				 // unwanted output from the linear
				 // solvers. If we want to suppress
				 // it, we have to include this file
				 // and add a line somewhere to the
				 // program; in this program, it was
				 // added to the main function.
#include <base/logstream.h>



				 // This is again the same
				 // LaplaceProblem class as in the
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
    FESystem<dim>            fe;
    DoFHandler<dim>      dof_handler;

    BlockSparsityPattern      sparsity_pattern;
    BlockSparseMatrix<double> system_matrix;

    BlockVector<double>       solution;
    BlockVector<double>       system_rhs;
};


				 // In the following, we declare two
				 // more classes, which will represent
				 // the functions of the
				 // dim-dimensional space denoting the
				 // right hand side and the
				 // non-homogeneous Dirichlet boundary
				 // values.
				 //
				 // Each of these classes is derived
				 // from a common, abstract base class
				 // Function, which declares the
				 // common interface which all
				 // functions have to follow. In
				 // particular, concrete classes have
				 // to overload the `value' function,
				 // which takes a point in
				 // dim-dimensional space as
				 // parameters and shall return the
				 // value at that point as a `double'
				 // variable.
				 //
				 // The `value' function takes a
				 // second argument, which we have
				 // here named `component': This is
				 // only meant for vector valued
				 // functions, where you may want to
				 // access a certain component of the
				 // vector at the point `p'. However,
				 // our functions are scalar, so we
				 // need not worry about this
				 // parameter and we will not use it
				 // in the implementation of the
				 // functions. Note that in the base
				 // class (Function), the declaration
				 // of the `value' function has a
				 // default value of zero for the
				 // component, so we will access the
				 // `value' function of the right hand
				 // side with only one parameter,
				 // namely the point where we want to
				 // evaluate the function.
				 //
				 // Note that the C++ language forces
				 // us to declare and define a
				 // constructor to the following
				 // classes even though they are
				 // empty. This is due to the fact
				 // that the base class has no default
				 // constructor (i.e. one without
				 // arguments), even though it has a
				 // constructor which has default
				 // values for all arguments.
template <int dim>
class RightHandSide : public Function<dim> 
{
  public:
    RightHandSide () : Function<dim>() {};
    
    virtual double value (const Point<dim>   &p,
			  const unsigned int  component = 0) const;
};



template <int dim>
class BoundaryValues : public Function<dim> 
{
  public:
    BoundaryValues () : Function<dim>() {};
    
    virtual double value (const Point<dim>   &p,
			  const unsigned int  component = 0) const;
};




				 // We wanted the right hand side
				 // function to be 4*(x**4+y**4) in
				 // 2D, or 4*(x**4+y**4+z**4) in
				 // 3D. Unfortunately, this is not as
				 // elegantly feasible dimension
				 // independently as much of the rest
				 // of this program, so we have to do
				 // it using a small
				 // loop. Fortunately, the compiler
				 // knows the size of the loop at
				 // compile time, i.e. the number of
				 // times the body will be executed,
				 // so it can optimize away the
				 // overhead needed for the loop and
				 // the result will be as fast as if
				 // we had used the formulas above
				 // right away.
				 //
				 // Note that the different
				 // coordinates (i.e. `x', `y', ...)
				 // of the point are accessed using
				 // the () operator.
template <int dim>
double RightHandSide<dim>::value (const Point<dim> &p,
				  const unsigned int) const 
{
  double return_value = deal_II_numbers::PI * deal_II_numbers::PI * dim;
  for (unsigned int i=0; i<dim; ++i)
    return_value *= std::sin(deal_II_numbers::PI*p(i));

  return return_value;
}


				 // The boundary values were to be
				 // chosen to be x*x+y*y in 2D, and
				 // x*x+y*y+z*z in 3D. This happens to
				 // be equal to the square of the
				 // vector from the origin to the
				 // point at which we would like to
				 // evaluate the function,
				 // irrespective of the dimension. So
				 // that is what we return:
template <int dim>
double BoundaryValues<dim>::value (const Point<dim> &p,
				   const unsigned int) const 
{
  return p.square();
}




				 // This is the constructor of the
				 // LaplaceProblem class. It specifies
				 // the desired polynomial degree of
				 // the finite elements and associates
				 // the DoFHandler to the
				 // triangulation just as in the
				 // previous example.
template <int dim>
LaplaceProblem<dim>::LaplaceProblem () :
                fe (FE_RaviartThomas<dim>(degree),1,FE_DGQ<dim>(degree),1),
		dof_handler (triangulation)
{}



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
				 // ``hyper_cube'', so we may use the
				 // same function in whatever
				 // dimension we are. Of course, the
				 // functions that create a hypercube
				 // in two and three dimensions are
				 // very much different, but that is
				 // something you need not care
				 // about. Let the library handle the
				 // difficult things.
				 //
				 // Likewise, associating a degree of
				 // freedom with each vertex is
				 // something which certainly looks
				 // different in 2D and 3D, but that
				 // does not need to bother you. This
				 // function therefore looks exactly
				 // like in the previous example,
				 // although it performs actions that
				 // in their details are quite
				 // different. The only significant
				 // difference is the number of cells
				 // resulting, which is much higher in
				 // three than in two space
				 // dimensions!
template <int dim>
void LaplaceProblem<dim>::make_grid_and_dofs ()
{
  GridGenerator::hyper_cube (triangulation, 0, 1);
  triangulation.refine_global (4);
  
  std::cout << "   Number of active cells: "
	    << triangulation.n_active_cells()
	    << std::endl
	    << "   Total number of cells: "
	    << triangulation.n_cells()
	    << std::endl;

  dof_handler.distribute_dofs (fe);
  DoFRenumbering::component_wise (dof_handler);
  
  std::vector<unsigned int> dofs_per_component (dim+1);
  DoFTools::count_dofs_per_component (dof_handler, dofs_per_component);
  const unsigned int n_u = dofs_per_component[0],
                     n_p = dofs_per_component[dim];

  std::cout << "   Number of degrees of freedom: "
	    << dof_handler.n_dofs()
            << " (" << n_u << '+' << n_p << ')'
	    << std::endl;
  
  sparsity_pattern.reinit (2,2);
  sparsity_pattern.block(0,0).reinit (n_u, n_u,
                                      dof_handler.max_couplings_between_dofs());
  sparsity_pattern.block(1,0).reinit (n_p, n_u,
                                      dof_handler.max_couplings_between_dofs());
  sparsity_pattern.block(0,1).reinit (n_u, n_p,
                                      dof_handler.max_couplings_between_dofs());
  sparsity_pattern.block(1,1).reinit (n_p, n_p,
                                dof_handler.max_couplings_between_dofs());
  sparsity_pattern.collect_sizes();
  DoFTools::make_sparsity_pattern (dof_handler, sparsity_pattern);
  sparsity_pattern.compress();

  system_matrix.reinit (sparsity_pattern);

  std::vector<unsigned int> block_components (2);
  block_components[0] = n_u;
  block_components[1] = n_p;
  solution.reinit (block_components);
  system_rhs.reinit (block_components);
}


Tensor<1,2> extract_u (const FEValues<2> &fe_values,
                       const unsigned int j,
                       const unsigned int q)
{
  Tensor<1,2> tmp;
  tmp[0] = fe_values.shape_value_component (j,q,0);
  tmp[1] = fe_values.shape_value_component (j,q,1);
  return tmp;
}



Tensor<1,3> extract_u (const FEValues<3> &fe_values,
                       const unsigned int j,
                       const unsigned int q)
{
  Tensor<1,3> tmp;
  tmp[0] = fe_values.shape_value_component (j,q,0);
  tmp[1] = fe_values.shape_value_component (j,q,1);
  tmp[2] = fe_values.shape_value_component (j,q,2);
  return tmp;
}





double extract_div_u (const FEValues<2> &fe_values,
                      const unsigned int j,
                      const unsigned int q)
{
  return fe_values.shape_grad_component (j,q,0)[0] +
    fe_values.shape_grad_component (j,q,1)[1];
}


double extract_div_u (const FEValues<3> &fe_values,
                      const unsigned int j,
                      const unsigned int q)
{
  return fe_values.shape_grad_component (j,q,0)[0] +
    fe_values.shape_grad_component (j,q,1)[1] +
    fe_values.shape_grad_component (j,q,2)[2];
}

  
template <int dim>
double extract_p (const FEValues<dim> &fe_values,
                  const unsigned int j,
                  const unsigned int q)
{
  return fe_values.shape_value_component (j,q,dim);
}



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
				 // difference to the pure
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
  QGauss<dim> quadrature_formula(degree+2);

				   // We wanted to have a non-constant
				   // right hand side, so we use an
				   // object of the class declared
				   // above to generate the necessary
				   // data. Since this right hand side
				   // object is only used in this
				   // function, we only declare it
				   // here, rather than as a member
				   // variable of the LaplaceProblem
				   // class, or somewhere else.
  const RightHandSide<dim> right_hand_side;

				   // Compared to the previous
				   // example, in order to evaluate
				   // the non-constant right hand side
				   // function we now also need the
				   // quadrature points on the cell we
				   // are presently on (previously,
				   // they were only needed on the
				   // unit cell, in order to compute
				   // the values and gradients of the
				   // shape function, which are
				   // defined on the unit cell
				   // however). We can tell the
				   // FEValues object to do for us by
				   // giving it the update_q_points
				   // flag:
  FEValues<dim> fe_values (fe, quadrature_formula, 
			   UpdateFlags(update_values    |
				       update_gradients |
				       update_q_points  |
				       update_JxW_values));

				   // Note that the following numbers
				   // depend on the dimension which we
				   // are presently using. However,
				   // the FE and Quadrature classes do
				   // all the necessary work for you
				   // and you don't have to care about
				   // the dimension dependent parts:
  const unsigned int   dofs_per_cell = fe.dofs_per_cell;
  const unsigned int   n_q_points    = quadrature_formula.n_quadrature_points;

  FullMatrix<double>   local_matrix (dofs_per_cell, dofs_per_cell);
  Vector<double>       local_rhs (dofs_per_cell);
  std::vector<double> rhs_values (n_q_points);

  std::vector<unsigned int> local_dof_indices (dofs_per_cell);

				   // Note here, that a cell is a
				   // quadrilateral in two space
				   // dimensions, but a hexahedron in
				   // 3D. In fact, the
				   // active_cell_iterator data type
				   // is something different,
				   // depending on the dimension we
				   // are in, but to the outside world
				   // they look alike and you will
				   // probably never see a difference
				   // although they are totally
				   // unrelated.
  typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(),
						 endc = dof_handler.end();
  for (; cell!=endc; ++cell)
    {
      fe_values.reinit (cell);
      local_matrix = 0;
      local_rhs = 0;

      right_hand_side.value_list (fe_values.get_quadrature_points(),
                                  rhs_values);
      
      for (unsigned int q=0; q<n_q_points; ++q) 
        for (unsigned int i=0; i<dofs_per_cell; ++i)
          {
            const Tensor<1,dim> phi_i_u = extract_u (fe_values, i, q);
            const double div_phi_i_u = extract_div_u (fe_values, i, q);
            const double phi_i_p = extract_p (fe_values, i, q);
            
            for (unsigned int j=0; j<dofs_per_cell; ++j)
              {
                const Tensor<1,dim> phi_j_u = extract_u (fe_values, j, q);
                const double div_phi_j_u = extract_div_u (fe_values, j, q);
                const double phi_j_p = extract_p (fe_values, j, q);
                
                local_matrix(i,j) += (phi_i_u * phi_j_u
                                      - div_phi_i_u * phi_j_p
                                      + phi_i_p * div_phi_j_u)
                                     * fe_values.JxW(q);
              }

            local_rhs(i) += phi_i_p *
                            rhs_values[q] *
                            fe_values.JxW(q);
          }
      
      cell->get_dof_indices (local_dof_indices);
      for (unsigned int i=0; i<dofs_per_cell; ++i)
        for (unsigned int j=0; j<dofs_per_cell; ++j)
          system_matrix.add (local_dof_indices[i],
                             local_dof_indices[j],
                             local_matrix(i,j));
      for (unsigned int i=0; i<dofs_per_cell; ++i)
        system_rhs(local_dof_indices[i]) += local_rhs(i);
    }
}


class SchurComplement 
{
  public:
    SchurComplement (const BlockSparseMatrix<double> &A)
                    :
                    A (A),
                    tmp1 (A.block(0,0).m()),
                    tmp2 (A.block(0,0).m())
      {}

    void vmult (Vector<double>       &dst,
                const Vector<double> &src) const
      {
        A.block(0,1).vmult (tmp1, src);

        SolverControl           solver_control (tmp1.size(),
                                                1e-8*tmp1.l2_norm());
        SolverCG<>              cg (solver_control);

        PreconditionSSOR<> precondition;
        precondition.initialize(A.block(0,0));
        cg.solve (A.block(0,0), tmp2, tmp1, precondition);

        std::cout << "     " << solver_control.last_step()
                  << " inner iterations needed to obtain convergence."
                  << std::endl;
        
        A.block(1,0).vmult (dst, tmp2);

        dst *= -1;
      }

  private:
    const BlockSparseMatrix<double> &A;

    mutable Vector<double> tmp1, tmp2;
};

        

				 // Solving the linear system of
				 // equation is something that looks
				 // almost identical in most
				 // programs. In particular, it is
				 // dimension independent, so this
				 // function is mostly copied from the
				 // previous example.
template <int dim>
void LaplaceProblem<dim>::solve () 
{
  {
    SolverControl           solver_control (system_matrix.block(0,0).m(),
                                            1e-6*system_rhs.block(1).l2_norm());
    SolverCG<>              cg (solver_control);

    cg.solve (SchurComplement(system_matrix), solution.block(1),
              system_rhs.block(1),
              PreconditionIdentity());
  
                                     // We have made one addition,
                                     // though: since we suppress output
                                     // from the linear solvers, we have
                                     // to print the number of
                                     // iterations by hand.
    std::cout << "   " << solver_control.last_step()
              << " CG mass matrix iterations needed to obtain convergence."
              << std::endl;
  }
  {
    Vector<double> tmp (system_matrix.block(0,0).m());
    system_matrix.block(0,1).vmult (tmp, solution.block(1));
    
    SolverControl           solver_control (system_matrix.block(0,0).m(),
                                            1e-6*tmp.l2_norm());
    SolverGMRES<>              cg (solver_control);

    cg.solve (system_matrix.block(0,0), solution.block(0),
              tmp, PreconditionIdentity());
  
                                     // We have made one addition,
                                     // though: since we suppress output
                                     // from the linear solvers, we have
                                     // to print the number of
                                     // iterations by hand.
    std::cout << "   " << solver_control.last_step()
              << " CG Schur complement iterations needed to obtain convergence."
              << std::endl;
  }
}



				 // This function also does what the
				 // respective one did in the previous
				 // example. No changes here for
				 // dimension independence either.
template <int dim>
void LaplaceProblem<dim>::output_results () const
{
  DataOut<dim> data_out;

  data_out.attach_dof_handler (dof_handler);
  data_out.add_data_vector (solution, "solution");

  data_out.build_patches (3);

				   // Only difference to the previous
				   // example: write output in GMV
				   // format, rather than for
				   // gnuplot. We use the dimension in
				   // the filename to generate
				   // distinct filenames for each run
				   // (in a better program, one would
				   // check whether `dim' can have
				   // other values than 2 or 3, but we
				   // neglect this here for the sake
				   // of brevity).
  std::ofstream output (dim == 2 ?
			"solution-2d" :
			"solution-3d");
  data_out.write_gnuplot (output);
}



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

    

				 // And this is the main function. It
				 // also looks mostly like in the
				 // previous example:
int main () 
{
				   // In the previous example, we had
				   // the output from the linear
				   // solvers about the starting
				   // residual and the number of the
				   // iteration where convergence was
				   // detected. This can be suppressed
				   // like this:
  deallog.depth_console (0);
				   // The rationale here is the
				   // following: the deallog
				   // (i.e. deal-log, not de-allog)
				   // variable represents a stream to
				   // which some parts of the library
				   // write output. It redirects this
				   // output to the console and if
				   // required to a file. The output
				   // is nested in a way that each
				   // function can use a prefix string
				   // (separated by colons) for each
				   // line of output; if it calls
				   // another function, that may also
				   // use its prefix which is then
				   // printed after the one of the
				   // calling function. Since output
				   // from functions which are nested
				   // deep below is usually not as
				   // important as top-level output,
				   // you can give the deallog
				   // variable a maximal depth of
				   // nested output for output to
				   // console and file. The depth zero
				   // which we gave here means that no
				   // output is written.

				   // After having done this
				   // administrative stuff, we can go
				   // on just as before: define one of
				   // these top-level objects and
				   // transfer control to
				   // it. Actually, now is the point
				   // where we have to tell the
				   // compiler which dimension we
				   // would like to use; all functions
				   // up to now including the classes
				   // were only templates and nothing
				   // has been compiled by now, but by
				   // declaring the following objects,
				   // the compiler will start to
				   // compile all the functions at the
				   // top using the template parameter
				   // replaced with a concrete value.
				   //
				   // For demonstration, we will first
				   // let the whole thing run in 2D
				   // and then in 3D:
  LaplaceProblem<2> laplace_problem_2d;
  laplace_problem_2d.run ();

//   LaplaceProblem<3> laplace_problem_3d;
//   laplace_problem_3d.run ();
  
  return 0;
}
