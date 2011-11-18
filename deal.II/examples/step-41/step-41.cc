/* $Id: step-4.cc 24093 2011-08-16 13:58:12Z bangerth $ */
/* Author: Wolfgang Bangerth, University of Heidelberg, 1999 */

/*    $Id: step-4.cc 24093 2011-08-16 13:58:12Z bangerth $       */
/*                                                                */
/*    Copyright (C) 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008, 2009, 2010, 2011 by the deal.II authors */
/*                                                                */
/*    This file is subject to QPL and may not be  distributed     */
/*    without copyright and license information. Please refer     */
/*    to the file deal.II/doc/license.html for the  text  and     */
/*    further information on this license.                        */

                                 // @sect3{Include files}

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
#include <deal.II/lac/compressed_sparsity_pattern.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_bicgstab.h>
#include <deal.II/lac/precondition.h>

#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/trilinos_vector.h>
#include <deal.II/lac/trilinos_precondition.h>
#include <deal.II/lac/trilinos_solver.h>

#include <deal.II/numerics/data_out.h>
#include <fstream>
#include <iostream>
#include <list>

using namespace dealii;

                                 // @sect3{The <code>Step41</code> class template}

				 // This class supply all function and variables
                                 // to an obstacle problem. The projection_active_set
                                 // function and the ConstaintMatrix are important
                                 // for the handling of the active set as we see
                                 // later.

template <int dim>
class Step41 
{
  public:
    Step41 ();
    void run ();
    
  private:
    void make_grid ();
    void setup_system();
    void assemble_system ();
    void assemble_mass_matrix ();
    void projection_active_set ();
    void solve ();
    void output_results (const std::string& title) const;

    Triangulation<dim>   triangulation;
    FE_Q<dim>            fe;
    DoFHandler<dim>      dof_handler;
    unsigned int         n_refinements;

    ConstraintMatrix     constraints;
    
    SparsityPattern                sparsity_pattern;
    TrilinosWrappers::SparseMatrix system_matrix;
    TrilinosWrappers::SparseMatrix system_matrix_complete;
    TrilinosWrappers::SparseMatrix mass_matrix;

    TrilinosWrappers::Vector       solution;
    TrilinosWrappers::Vector       tmp_solution;
    TrilinosWrappers::Vector       system_rhs;
    TrilinosWrappers::Vector       system_rhs_complete;
    TrilinosWrappers::Vector       resid_vector;
    TrilinosWrappers::Vector       active_set;
    TrilinosWrappers::Vector       diag_mass_matrix_vector;
};


                                 // @sect3{Right hand side and boundary values}

template <int dim>
class RightHandSide : public Function<dim> 
{
  public:
    RightHandSide () : Function<dim>() {}
    
    virtual double value (const Point<dim>   &p,
			  const unsigned int  component = 0) const;
};

template <int dim>
class BoundaryValues : public Function<dim> 
{
  public:
    BoundaryValues () : Function<dim>() {}
    
    virtual double value (const Point<dim>   &p,
			  const unsigned int  component = 0) const;
};

template <int dim>
class Obstacle : public Function<dim> 
{
  public:
    Obstacle () : Function<dim>() {}
    
    virtual double value (const Point<dim>   &p,
			  const unsigned int  component = 0) const;
};



				 // For this example, we choose as right hand
				 // side function a constant force density
                                 // like the gravitation attraction.
template <int dim>
double RightHandSide<dim>::value (const Point<dim> &p,
				  const unsigned int /*component*/) const 
{
  double return_value = -10;

  return return_value;
}


				 // As boundary values, we choose the zero.
template <int dim>
double BoundaryValues<dim>::value (const Point<dim> &p,
				   const unsigned int /*component*/) const 
{
  double return_value = 0;

  return return_value;
}


				 // The obstacle function describes a cascaded
                                 // barrier. So if the gravitation attraction
                                 // pulls the membrane down it blows over the
                                 // steps.
template <int dim>
double Obstacle<dim>::value (const Point<dim> &p,
			     const unsigned int /*component*/) const 
{
  double return_value = 0;

  if (p (0) < -0.5)
    return_value = -0.2;
  else if (p (0) >= -0.5 && p (0) < 0.0)
    return_value = -0.4;
  else if (p (0) >= 0.0 && p (0) < 0.5)
    return_value = -0.6;
  else
    return_value = -0.8;

  return return_value;
}



                                 // @sect3{Implementation of the <code>Step41</code> class}
            

                                 // @sect4{Step41::Step41}

template <int dim>
Step41<dim>::Step41 ()
		:
                fe (1),
		dof_handler (triangulation)
{}


                                 // @sect4{Step41::make_grid}

                                 // We solve our obstacle problem on the square
                                 // $[-1,1]\times [-1,1]$ in 2D.
template <int dim>
void Step41<dim>::make_grid ()
{
  GridGenerator::hyper_cube (triangulation, -1, 1);
  n_refinements = 6;
  triangulation.refine_global (n_refinements);
  
  std::cout << "   Number of active cells: "
	    << triangulation.n_active_cells()
	    << std::endl
	    << "   Total number of cells: "
	    << triangulation.n_cells()
	    << std::endl;
}

                                 // @sect4{Step41::setup_system}

template <int dim>
void Step41<dim>::setup_system ()
{
  dof_handler.distribute_dofs (fe);

  std::cout << "   Number of degrees of freedom: "
	    << dof_handler.n_dofs()
	    << std::endl;

  CompressedSparsityPattern c_sparsity(dof_handler.n_dofs());
  DoFTools::make_sparsity_pattern (dof_handler, c_sparsity, constraints, false);
  sparsity_pattern.copy_from(c_sparsity);
  
  system_matrix.reinit (sparsity_pattern);
  system_matrix_complete.reinit (sparsity_pattern);
  mass_matrix.reinit (sparsity_pattern);
  
  solution.reinit (dof_handler.n_dofs());
  tmp_solution.reinit (dof_handler.n_dofs());
  system_rhs.reinit (dof_handler.n_dofs());
  system_rhs_complete.reinit (dof_handler.n_dofs());
  resid_vector.reinit (dof_handler.n_dofs());
  active_set.reinit (dof_handler.n_dofs());
  diag_mass_matrix_vector.reinit (dof_handler.n_dofs());
}


                                 // @sect4{Step41::assemble_system}


				 // At once with assembling the system matrix and
                                 // right-hand-side we apply the constraints
                                 // to our system. The constraint consists not
                                 // only of the zero Dirichlet boundary values,
                                 // in addition they contain the obstacle values.
                                 // The projection_active_set function are used
                                 // to fill the ConstraintMatrix.
template <int dim>
void Step41<dim>::assemble_system () 
{  
  QGauss<dim>  quadrature_formula(2);

  const RightHandSide<dim> right_hand_side;

  FEValues<dim> fe_values (fe, quadrature_formula, 
			   update_values   | update_gradients |
                           update_quadrature_points | update_JxW_values);

  const unsigned int   dofs_per_cell = fe.dofs_per_cell;
  const unsigned int   n_q_points    = quadrature_formula.size();

  FullMatrix<double>   cell_matrix (dofs_per_cell, dofs_per_cell);
  TrilinosWrappers::Vector       cell_rhs (dofs_per_cell);

  std::vector<unsigned int> local_dof_indices (dofs_per_cell);

  typename DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();
  
  for (; cell!=endc; ++cell)
    {
      fe_values.reinit (cell);
      cell_matrix = 0;
      cell_rhs = 0;

      for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
	for (unsigned int i=0; i<dofs_per_cell; ++i)
	  {
	    for (unsigned int j=0; j<dofs_per_cell; ++j)
	      cell_matrix(i,j) += (fe_values.shape_grad (i, q_point) *
				   fe_values.shape_grad (j, q_point) *
				   fe_values.JxW (q_point));

	    cell_rhs(i) += (fe_values.shape_value (i, q_point) *
			    right_hand_side.value (fe_values.quadrature_point (q_point)) *
			    fe_values.JxW (q_point));
	  }

      cell->get_dof_indices (local_dof_indices);

                                       // This function apply the constraints
                                       // to the system matrix and system rhs.
                                       // The true parameter is set to make sure
                                       // that the system rhs contains correct
                                       // values in the rows with inhomogeneity
                                       // constraints.
      constraints.distribute_local_to_global (cell_matrix, cell_rhs,
                                              local_dof_indices,
                                              system_matrix, system_rhs, true);
    }
}

template <int dim>
void Step41<dim>::assemble_mass_matrix () 
{  
  QTrapez<dim>  quadrature_formula;

  FEValues<dim> fe_values (fe, quadrature_formula, 
			   update_values   | update_quadrature_points | update_JxW_values);

  const unsigned int   dofs_per_cell = fe.dofs_per_cell;
  const unsigned int   n_q_points    = quadrature_formula.size();

  FullMatrix<double>   cell_matrix (dofs_per_cell, dofs_per_cell);

  std::vector<unsigned int> local_dof_indices (dofs_per_cell);

  typename DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();
  
  for (; cell!=endc; ++cell)
    {
      fe_values.reinit (cell);
      cell_matrix = 0;

      for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
	for (unsigned int i=0; i<dofs_per_cell; ++i)
	  for (unsigned int j=0; j<dofs_per_cell; ++j)
	    cell_matrix(i,j) += (fe_values.shape_value (i, q_point) *
				 fe_values.shape_value (j, q_point) *
				 fe_values.JxW (q_point));

      cell->get_dof_indices (local_dof_indices);

                                       // This function apply the constraints
                                       // to the system matrix and system rhs.
                                       // The true parameter is set to make sure
                                       // that the system rhs contains correct
                                       // values in the rows with inhomogeneity
                                       // constraints.
      constraints.distribute_local_to_global (cell_matrix,
                                              local_dof_indices,
                                              mass_matrix);
    }
}

                                 // @sect4{Step41::projection_active_set}

				 // Updating of the active set which means to
                                 // set a inhomogeneity constraint in the
                                 // ConstraintMatrix. At the same time we set
                                 // the solution to the correct value - the obstacle value.
                                 // To control the active set we use the vector
                                 // active_set which contains a zero in a component
                                 // that is not in the active set and elsewise a
                                 // one. With the output file you can visualize it.
template <int dim>
void Step41<dim>::projection_active_set ()
{
  const Obstacle<dim>     obstacle;
  std::vector<bool>       vertex_touched (triangulation.n_vertices(),
				    false);
  unsigned int            counter_contact_constraints = 0; 
  
  typename DoFHandler<dim>::active_cell_iterator
  cell = dof_handler.begin_active(),
  endc = dof_handler.end();

  constraints.clear();

                                       // to find and supply the constraints for the
                                       // obstacle condition
  active_set = 0.0;
  const double c = 100.0;
  for (; cell!=endc; ++cell)
    for (unsigned int v=0; v<GeometryInfo<2>::vertices_per_cell; ++v)
      {
	unsigned int index_x = cell->vertex_dof_index (v,0);

				       // the local row where
	Point<dim> point (cell->vertex (v)[0], cell->vertex (v)[1]);
	double obstacle_value = obstacle.value (point);
	double solution_index_x = solution (index_x);

                                       // To decide which dof belongs to the
                                       // active-set. For that we scale the
                                       // residual-vector with the cell-size and
                                       // the diag-entry of the mass-matrix.

                                       // TODO: I have to check the condition
	if (resid_vector (index_x) +
	   diag_mass_matrix_vector (index_x)*c*(obstacle_value - solution_index_x) > 0)
	{
	  constraints.add_line (index_x);
	  constraints.set_inhomogeneity (index_x, obstacle_value);
	  solution (index_x) = obstacle_value;
	  active_set (index_x) = 1.0;
	  
	  if (vertex_touched[cell->vertex_index(v)] == false)
	  {
	    vertex_touched[cell->vertex_index(v)] = true;
	    counter_contact_constraints += 1;
	  }
	}
      }
  std::cout<< "Number of Contact-Constaints: " << counter_contact_constraints <<std::endl;

                                       // To supply the boundary values of the
                                       // dirichlet-boundary in constraints
  VectorTools::interpolate_boundary_values (dof_handler,
					    0,
					    BoundaryValues<dim>(),
					    constraints);
  constraints.close ();
}

                                 // @sect4{Step41::solve}

template <int dim>
void Step41<dim>::solve () 
{
  ReductionControl        reduction_control (100, 1e-12, 1e-3);
  SolverCG<TrilinosWrappers::Vector>   solver (reduction_control); 
  TrilinosWrappers::PreconditionAMG precondition;
  precondition.initialize (system_matrix);

  solver.solve (system_matrix, solution, system_rhs, precondition);

  std::cout << "Initial error: " << reduction_control.initial_value() <<std::endl;
  std::cout << "   " << reduction_control.last_step()
	    << " CG iterations needed to obtain convergence with an error: "
	    <<  reduction_control.last_value()
	    << std::endl;

  constraints.distribute (solution);
}

                                 // @sect4{Step41::output_results}

				 // We use the vtk-format for the output.
                                 // The file contains the displacement,
                                 // the residual and active set vectors.
template <int dim>
void Step41<dim>::output_results (const std::string& title) const
{
  DataOut<dim> data_out;
  
  data_out.attach_dof_handler (dof_handler);
  // data_out.add_data_vector (tmp_solution, "Displacement");
  // data_out.add_data_vector (resid_vector, "Residual");
  data_out.add_data_vector (active_set, "ActiveSet");

  data_out.build_patches ();

  std::ofstream output_vtk ((title + ".vtk").c_str ());
  data_out.write_gnuplot (output_vtk);
}



                                 // @sect4{Step41::run}

                                 // This is the function which has the
				 // top-level control over everything.
                                 // Here the active set method is implemented.

                                 // TODO: I have to compare it with the algorithm
                                 // in the Wohlmuth-paper
template <int dim>
void Step41<dim>::run () 
{
  std::cout << "Solving problem in " << dim << " space dimensions." << std::endl;

  make_grid();
  setup_system ();

  constraints.clear ();
  VectorTools::interpolate_boundary_values (dof_handler,
					    0,
					    BoundaryValues<dim>(),
					    constraints);
  constraints.close ();
  ConstraintMatrix constraints_complete (constraints);
  assemble_system ();
  solve ();

                                       // to save the system_matrix and the
                                       // rhs to compute the residual in every
                                       // step of the active-set-iteration
  system_matrix_complete.copy_from (system_matrix);
  system_rhs_complete = system_rhs;
 
                                       // to compute the factor which is used
                                       // to scale the residual. You can consider
                                       // this diagonal matrix as the discretization
                                       // of a lagrange multiplier for the
                                       // contact force
  assemble_mass_matrix ();
  for (unsigned int j=0; j<solution.size (); j++)
    diag_mass_matrix_vector (j) = mass_matrix.diag_element (j);

  resid_vector = 0;
  resid_vector -= system_rhs_complete;
  system_matrix_complete.vmult_add  (resid_vector, solution);

                                      // to compute a start active set
  std::cout<< "Update Active Set:" <<std::endl;
  projection_active_set ();
  TrilinosWrappers::Vector       active_set_old (active_set);
  for (unsigned int i=0; i<solution.size (); i++)
    {
      std::cout<< "Assemble System:" <<std::endl;
      system_matrix = 0;
      system_rhs = 0;
      assemble_system ();

      std::cout<< "Solve System:" <<std::endl;
      solve ();
      tmp_solution = solution;

      resid_vector = 0;
      resid_vector -= system_rhs_complete;
      system_matrix_complete.vmult_add  (resid_vector, solution);

      std::cout<< "Update Active Set:"<<std::endl;
      projection_active_set ();

      for (unsigned int k = 0; k<solution.size (); k++)
      	if (active_set (k) == 1)
      	  resid_vector (k) = 0;

      std::cout<< "Create Output:" <<std::endl;
      std::ostringstream filename_output;
      filename_output << "output_";
      filename_output << i;
      output_results (filename_output.str ());

                                     // the residual of the non-contact part
                                     // of the system serves as an additional
                                     // control which is not necassary for
                                     // for the primal-dual active set strategy
      double resid = resid_vector.l2_norm ();
      std::cout<< i << ". Residual of the non-contact part of the system = " << resid <<std::endl;

                                      // if both the old and the new
                                      // active set are identical the
                                      // computation stops
      if (active_set == active_set_old)
	break;
      active_set_old = active_set;
    }
}


                                 // @sect3{The <code>main</code> function}

				 // And this is the main function. It also
				 // looks mostly like in step-3, but if you
				 // look at the code below, note how we first
				 // create a variable of type
				 // <code>Step41@<2@></code> (forcing
				 // the compiler to compile the class template
				 // with <code>dim</code> replaced by
				 // <code>2</code>) and run a 2d simulation,
				 // and then we do the whole thing over in 3d.
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
				 // <code>laplace_problem_2d</code> variable
				 // goes out of scope (and releases the memory
				 // it holds) before we move on to allocate
				 // memory for the 3d case. Without the
				 // additional braces, the
				 // <code>laplace_problem_2d</code> variable
				 // would only be destroyed at the end of the
				 // function, i.e. after running the 3d
				 // problem, and would needlessly hog memory
				 // while the 3d run could actually use it.
                                 //
                                 // Finally, the first line of the function is
                                 // used to suppress some output.  Remember
                                 // that in the previous example, we had the
                                 // output from the linear solvers about the
                                 // starting residual and the number of the
                                 // iteration where convergence was
                                 // detected. This can be suppressed through
                                 // the <code>deallog.depth_console(0)</code>
                                 // call.
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
int main (int argc, char *argv[]) 
{
  deallog.depth_console (0);

  Utilities::MPI::MPI_InitFinalize mpi_initialization (argc, argv);

  Step41<2> laplace_problem_2d;
  laplace_problem_2d.run ();
  
  return 0;
}
