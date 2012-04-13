/* $Id$ */
/* Author: Wolfgang Bangerth, University of Heidelberg, 2000 */

/*    $Id$       */
/*                                                                */
/*    Copyright (C) 2000, 2001, 2002, 2003, 2004, 2006, 2007, 2008, 2010, 2011 by the deal.II authors */
/*                                                                */
/*    This file is subject to QPL and may not be  distributed     */
/*    without copyright and license information. Please refer     */
/*    to the file deal.II/doc/license.html for the  text  and     */
/*    further information on this license.                        */

                                 // @sect3{Include files}

				 // The first few files have already
				 // been covered in previous examples
				 // and will thus not be further
				 // commented on.
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/compressed_sparsity_pattern.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/grid/tria.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/numerics/vectors.h>
#include <deal.II/numerics/matrices.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/error_estimator.h>


#include <fstream>
#include <iostream>

				// We will use adaptive mesh refinement between Newton
				// interations. To do so, we need to be able to work
				// with a solution on the new mesh, although it was
				// computed on the old one. The SolutionTransfer
				// class transfers the solution to the new mesh.

#include <deal.II/numerics/solution_transfer.h>

				// In this tutorial, we can't use the CG-method as a solver, as
				// described above, but we use the minimal residual method, which
				// is included with this file.

#include <deal.II/lac/solver_minres.h>

				// As in previous programs:

using namespace dealii;


								// @sect3{The <code>Step15</code> class template}

				// The class template is basically the same as in step 6.
				// Four additions are made: There are two solution vectors,
				// one for the Newton update, and one for the solution of
				// the original pde. Also we need a double for the residual
				// of the Newton method, an integer, which counts the mesh
				// refinements and a bool for the boundary condition in the first
				// Newton step.

template <int dim>
class Step15
{
  public:
    Step15 ();
    ~Step15 ();

    void run ();

  private:
    void setup_system ();
    void assemble_system ();
    void solve ();
    void refine_grid ();

    Triangulation<dim>   triangulation;

    DoFHandler<dim>      dof_handler;
    FE_Q<dim>            fe;

    ConstraintMatrix     hanging_node_constraints;

    SparsityPattern      sparsity_pattern;
    SparseMatrix<double> system_matrix;

    Vector<double>       present_solution;
    Vector<double>       newton_update;
    Vector<double>       system_rhs;



    double 				 res;
    unsigned int         refinement;

    			// As described in the Introduction, the first Newton iteration
    			// is special, because of the boundary condition. To implement
    			// these correctly, there is a bool, which is true in the first
    			// step and false ever after.
    bool 				 first_step;
};

								// @sect3{Boundary condition}

				// The boundary condition is implemented just like in step 4.
				// It was chosen as $g(x,y)=sin(2 \pi (x+y))$ in this example.

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
  double return_value=sin(2*M_PI*(p[0]+p[1]));
  return return_value;
}

								// @sect3{The <code>Step15</code> class implementation}

								// @sect4{Step15::Step15}

				// The constructor and destructor of the class are the same
				// as in the first few tutorials.

template <int dim>
Step15<dim>::Step15 ()
		:
		dof_handler (triangulation),
                fe (2)
{}



				//
template <int dim>
Step15<dim>::~Step15 ()
{
  dof_handler.clear ();
}

								// @sect4{Step15::setup_system}

				// As always in the setup-system function, we setup the variables
				// of the finite element method. There are same differences to
				// step 6, because we don't have to solve one pde over all,
				// but one in every Newton step. Also the starting function
				// has to be setup in the first step.

template <int dim>
void Step15<dim>::setup_system ()
{

				// This function will be called, every time we refine the mesh
				// to resize the system matrix,  Newton update - and right hand
				// side vector and to set the right values of hanging nodes to
				// get a continuous solution.
				// But only the first time, the starting solution has to be
				// initialized. Also the vector of the solution will be
				// resized in the <code>refine_grid</code> function, while the
				// vector is transfered to the new mesh.

  if(first_step)
  {
	dof_handler.distribute_dofs (fe);
	present_solution.reinit (dof_handler.n_dofs());
	for(unsigned int i=0; i<dof_handler.n_dofs();++i)
	{
	  present_solution(i)=0;
	}
				// The constraint matrix, holding a list of the hanging nodes,
				// will be setup in the <code>refine_grid</code> function
				// after refining the mesh.

	hanging_node_constraints.clear ();
	DoFTools::make_hanging_node_constraints (dof_handler,
				   hanging_node_constraints);
	hanging_node_constraints.close ();
  }


  	  	  	    // The remaining parts of the function are the same as in step 6.

  newton_update.reinit (dof_handler.n_dofs());
  system_rhs.reinit (dof_handler.n_dofs());

  CompressedSparsityPattern c_sparsity(dof_handler.n_dofs());
  DoFTools::make_sparsity_pattern (dof_handler, c_sparsity);

  hanging_node_constraints.condense (c_sparsity);

  sparsity_pattern.copy_from(c_sparsity);
  system_matrix.reinit (sparsity_pattern);
}

								// @sect4{Step15::assemble_system}

				// This function does the same as in the previous tutorials.
				// The only additional step is the correct implementation of
				// the boundary condition and the usage of the gradients of
				// the old solution.

template <int dim>
void Step15<dim>::assemble_system ()
{
  const QGauss<dim>  quadrature_formula(3);

  system_matrix = 0;
  system_rhs = 0;

  FEValues<dim> fe_values (fe, quadrature_formula,
			   update_values    |  update_gradients |
			   update_quadrature_points  |  update_JxW_values);

  const unsigned int   dofs_per_cell = fe.dofs_per_cell;
  const unsigned int   n_q_points    = quadrature_formula.size();

  FullMatrix<double>   cell_matrix (dofs_per_cell, dofs_per_cell);
  Vector<double>       cell_rhs (dofs_per_cell);

  std::vector<unsigned int> local_dof_indices (dofs_per_cell);

  typename DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();
  for (; cell!=endc; ++cell)
    {
      cell_matrix = 0;
      cell_rhs = 0;

      fe_values.reinit (cell);


      for (unsigned int q_point = 0; q_point < n_q_points; ++q_point) {

    	  	  // To setup up the linear system, the gradient of the old solution
    	  	  // in the quadrature points is needed. For this purpose there is
    	  	  // is a function, which will write these gradients in a vector,
    	  	  // where every component of the vector is a vector itself:

			std::vector<Tensor<1, dim> > gradients(n_q_points);
			fe_values.get_function_gradients(present_solution, gradients);

			  // Having the gradients of the old solution in the quadrature
			  // points, we are able to compute the coefficients $a_{n}$
			  // in these points.

			const double coeff = 1/sqrt(1 + gradients[q_point] * gradients[q_point]);

			  // The assembly of the system then is the same as always, except
			  // of the damping parameter of the Newton method, which we set on
			  // 0.1 in this case.

			for (unsigned int i = 0; i < dofs_per_cell; ++i) {
				for (unsigned int j = 0; j < dofs_per_cell; ++j) {
					cell_matrix(i, j) += (fe_values.shape_grad(i, q_point)
							* coeff
							* (fe_values.shape_grad(j, q_point)
								- coeff * coeff
								* (fe_values.shape_grad(j, q_point)
								* gradients[q_point])
								* gradients[q_point])
							* fe_values.JxW(q_point));
				}

				cell_rhs(i) -=0.1 *
						(fe_values.shape_grad(i, q_point) * coeff
								* gradients[q_point] * fe_values.JxW(q_point));
			}
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
  hanging_node_constraints.condense (system_matrix);
  hanging_node_constraints.condense (system_rhs);
  std::map<unsigned int,double> boundary_values;

    			// As described above, there is a different boundary condition
  	  	  	  	// in the first Newton step than in the later ones. This is
  	  	  	  	// implemented with the help of the bool first_step, which
  	  	  	  	// will later be false for all times. Starting with the zero-
  	  	  	  	// function in the first step, we have to set the boundary
  	  	  	    // condition $\delta u^{0}=g$ on $\partial \Omega $:

  if(first_step)
  {
	VectorTools::interpolate_boundary_values (dof_handler,
					      0,
				          BoundaryValues<dim>(),
						  boundary_values);
  }
  	  	  	    // In later steps, the Newton update has to have homogeneous
  	  	  	    // boundary conditions, in order for the solution to have the
  	  	  	    // right ones.

  else{
    VectorTools::interpolate_boundary_values (dof_handler,
	  				      0,
	  				      ZeroFunction<dim>(),
	  				      boundary_values);
  }

  MatrixTools::apply_boundary_values (boundary_values,
				      system_matrix,
				      newton_update,
				      system_rhs);
}

								// @sect4{Step15::solve}

				// The solve function is the same as always, we just have to
				// implement the minimal residual method as a solver and
				// apply the Newton update to the solution.

template <int dim>
void Step15<dim>::solve ()
{
  res=system_rhs.l2_norm();
  SolverControl           solver_control (1000, res*1e-6);
  SolverMinRes<>              solver (solver_control);

  PreconditionSSOR<> preconditioner;
  preconditioner.initialize(system_matrix, 1.2);

  solver.solve (system_matrix, newton_update, system_rhs,
		preconditioner);

  hanging_node_constraints.distribute (newton_update);

  	  	  	  	  // In this step, the old solution is updated to the new one:

  present_solution += newton_update;
}

								// @sect4{Step15::refine_grid}

				// The first part of this function is the same as in step 6.
				// But after refining the mesh we have to transfer the old
				// solution to the new one, which is done with the help of
				// the SolutionTransfer class.


template <int dim>
void Step15<dim>::refine_grid ()
{
  Vector<float> estimated_error_per_cell (triangulation.n_active_cells());

  KellyErrorEstimator<dim>::estimate (dof_handler,
				      QGauss<dim-1>(3),
				      typename FunctionMap<dim>::type(),
				      present_solution,
				      estimated_error_per_cell);

  GridRefinement::refine_and_coarsen_fixed_number (triangulation,
						   estimated_error_per_cell,
						   0.3, 0.03);

  	  	    	// Then we need an additional step: if, for example,
  	  	  	  	// you flag a cell that is once more refined than its neighbor,
  	  	  	  	// and that neighbor is not flagged for refinement, we would end
  	  	  	  	// up with a jump of two refinement levels across a cell interface.
  	  	  	  	// To avoid these situations, the library will
  	  	  	  	// silently also have to refine the neighbor cell once. It does so
  	  	  	  	// by calling the Triangulation::prepare_coarsening_and_refinement
  	  	  	  	// function before actually doing the refinement and coarsening.
  	  	  	  	// This function flags a set of additional cells for refinement or
  	  	  	  	// coarsening, to enforce rules like the one-hanging-node rule.
  	  	  	  	// The cells that are flagged for refinement and coarsening after
  	  	  	  	// calling this function are exactly the ones that will actually
  	  	  	  	// be refined or coarsened. Since the SolutionTransfer class needs
  	  	  	  	// this information in order to store the data from the old mesh
  	  	  	  	// and transfer to the new one.

  triangulation.prepare_coarsening_and_refinement ();

  	  	  	  	// With this out of the way, we initialize a SolutionTransfer
  	  	  	  	// object with the present DoFHandler and attach the solution
  	  	  	  	// vector to it:

  SolutionTransfer<dim> solution_transfer(dof_handler);
  solution_transfer.prepare_for_coarsening_and_refinement(present_solution);

  	  	  	  	// Then we do the actual refinement, and distribute degrees
  	  	  	  	// of freedom on the new mesh:

  triangulation.execute_coarsening_and_refinement();
  dof_handler.distribute_dofs(fe);

  	  	  	  	// Finally, we retrieve the old solution interpolated to the new
  	  	  	  	// mesh. Since the SolutionTransfer function does not actually
  	  	  	  	// store the values of the old solution, but rather indices, we
  	  	  	  	// need to preserve the old solution vector until we have gotten
  	  	  	  	// the new interpolated values. Thus, we have the new values
  	  	  	  	// written into a temporary vector, and only afterwards write
  	  	  	  	// them into the solution vector object:

  Vector<double> tmp(dof_handler.n_dofs());
  solution_transfer.interpolate(present_solution,tmp);
  present_solution=tmp;

  	  	  	  	// Having refined the mesh, there might be new nodal points on
  	  	  	  	// the boundary. These have just interpolated values, but
  	  	  	  	// not the right boundary values. This is fixed up, by
  	  	  	  	// setting all boundary nodals explicit to the right value:

  std::map<unsigned int, double> boundary_values2;
  VectorTools::interpolate_boundary_values(dof_handler,
		  	  	  	  	  	  	  	  	    0,
		  	  	  	  	  	  	  	  	    BoundaryValues<dim>(),
		  	  	  	  	  	  	  	  	    boundary_values2);
  for (std::map<unsigned int, double>::const_iterator
		  p = boundary_values2.begin();
		  p != boundary_values2.end();
		  ++p)
	  present_solution(p->first) = p->second;

  	  	  	  	// On the new mesh, there are different hanging nodes, which shall
  	  	  	  	// be enlisted in a matrix like before. To ensure there are no
  	  	  	  	// hanging nodes of the old mesh in the matrix, it's first cleared:
  hanging_node_constraints.clear();

  	  	  	  	// After doing so, the hanging nodes of the new mesh can be
  	  	  	  	// enlisted in the matrix, like before. Calling the
  	  	  	  	// <code>setup_system</code> function in the <code>run</code>
  	  	  	  	// function again after this, the hanging nodes don't have to
  	  	  	  	// be enlisted there once more.

  DoFTools::make_hanging_node_constraints(dof_handler, hanging_node_constraints);
  hanging_node_constraints.close();
  hanging_node_constraints.distribute(present_solution);
}

								// @sect4{Step15::run}

				// In the run function, the first grid is build. Also in this
				// function, the Newton iteration is implemented.

template <int dim>
void Step15<dim>::run ()
{

				// The integer refinement counts the mesh refinements. Obviously
				// starting the program, it should be zero.
	refinement=0;
	first_step=true;

				// As described in the introduction, the domain is a unitball around
				// the origin. The Mesh is globally refined two times, not to start
				// on the coarse mesh, which consists only of five cells.

	GridGenerator::hyper_ball (triangulation);
	static const HyperBallBoundary<dim> boundary;
	triangulation.set_boundary (0, boundary);
	triangulation.refine_global(2);

				// The Newton iteration starts here. During the first step, there is
				// no residual computed, so the bool is needed here to enter the
				// iteration scheme. Later the Newton method will continue until the
				// residual is less than $10^{-3}$.

	while(first_step || (res>1e-3))
	{

				// In the first step, we compute the solution on the two times globally 
				// refined mesh. After that the mesh will be refined 
				// adaptively, in order to not get too many cells. The refinement 
				// is the first thing done every time we restart the process in the while-loop.
		if(!first_step)
		{
			refine_grid();

			std::cout<<"********mesh-refinement:"<<refinement+1<<" ********"<<std::endl;
			refinement++;
		}


				// First thing to do after refining the mesh, is to setup the vectors,
				// matrices, etc., which is done in the <code>setup_system</code>
				// function.

		setup_system();

				// On every mesh there are done five Newton steps, in order to get a
				// better solution, before the mesh gets too fine and the computations
				// take more time.

		for(unsigned int i=0; i<5;++i)
		{

				// In every Newton step the system matrix and the right hand side
				// have to be computed.

			assemble_system ();
			solve ();
			first_step=false;
			std::cout<<"residual:"<<res<<std::endl;
		}

			    // The fifth solution, as well as the Newton update,
				// on every mesh will be written in a vtk-file,
				// in order to show the convergence of the solution.

		Assert (refinement < 100, ExcNotImplemented());

		DataOutBase::EpsFlags vtk_flags;

		DataOut<dim> data_out;
		data_out.set_flags (vtk_flags);

		data_out.attach_dof_handler (dof_handler);
		data_out.add_data_vector (newton_update, "update");
		data_out.add_data_vector (present_solution, "solution");
		data_out.build_patches (6);
		std::string filename = "solution-";
		if(refinement<10)
		{
		  filename += ('0' + refinement);
		}
		else{
		  filename += ('0' + refinement/10);
		  filename += ('0' + refinement%10);
		}
		filename += ".vtk";
		std::ofstream output (filename.c_str());
		data_out.write_vtk (output);

	}
}

								// @ sect4{The main function}

				// Finally the main function, this follows the scheme of all other main
				// functions:

int main ()
{

  try
    {
      deallog.depth_console (0);

      Step15<2> laplace_problem_2d;
      laplace_problem_2d.run ();
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

