/* $Id$ */
/* Author: Wolfgang Bangerth, University of Heidelberg, 2002 */

/*    $Id$       */
/*    Version: $Name$                                          */
/*                                                                */
/*    Copyright (C) 2002, 2003 by the deal.II authors                   */
/*                                                                */
/*    This file is subject to QPL and may not be  distributed     */
/*    without copyright and license information. Please refer     */
/*    to the file deal.II/doc/license.html for the  text  and     */
/*    further information on this license.                        */

#include <base/quadrature_lib.h>
#include <base/function.h>
#include <base/logstream.h>
#include <lac/vector.h>
#include <lac/full_matrix.h>
#include <lac/sparse_matrix.h>
                                 //
#include <lac/solver_gmres.h>
#include <lac/solver_cg.h>
#include <lac/vector_memory.h>
#include <lac/precondition.h>
#include <grid/tria.h>
#include <dofs/dof_handler.h>
#include <grid/grid_generator.h>
#include <grid/tria_accessor.h>
#include <grid/tria_iterator.h>
#include <grid/tria_boundary_lib.h>
#include <dofs/dof_accessor.h>
#include <dofs/dof_tools.h>
#include <fe/fe_values.h>
#include <numerics/vectors.h>
#include <numerics/matrices.h>
#include <numerics/data_out.h>
#include <numerics/solution_transfer.h>
#include <fe/fe_q.h>
#include <grid/grid_out.h>

#include <grid/grid_refinement.h>

#include <dofs/dof_constraints.h>

#include <numerics/error_estimator.h>

#include <fstream>
#include <iostream>

#ifdef HAVE_STD_STRINGSTREAM
#  include <sstream>
#else
#  include <strstream>
#endif



template <int dim>
class MinimizationProblem 
{
  public:
    MinimizationProblem  ();
    ~MinimizationProblem  ();
    void run ();
    
  private:
    void setup_system ();
    void assemble_newton_step (const bool p);
    void solve ();
    void initialize ();
    void refine_grid ();
    double energy () const;
    void output_results (const unsigned int cycle) const;

    Triangulation<dim>   triangulation;

    FE_Q<dim>            fe;
    DoFHandler<dim>      dof_handler;

    ConstraintMatrix     hanging_node_constraints;
    
    SparsityPattern      sparsity_pattern;
    SparseMatrix<double> newton_matrix;

    Vector<double>       present_solution;
    Vector<double>       newton_residual;
    Vector<double>       newton_update;
};



class InitializationValues : public Function<1> 
{
  public:
    InitializationValues () : Function<1>() {};
    
    virtual double value (const Point<1>     &p,
			  const unsigned int  component = 0) const;
};



double InitializationValues::value (const Point<1> &p,
                                    const unsigned int) const 
{
  return std::pow(p(0), 1./3.);
};



template <int dim>
MinimizationProblem<dim>::MinimizationProblem () :
                fe (1),
		dof_handler (triangulation)
{};


template <int dim>
MinimizationProblem<dim>::~MinimizationProblem () 
{
  dof_handler.clear ();
};



template <int dim>
void MinimizationProblem<dim>::setup_system ()
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
};


template <int dim>
double gradient_power (const Tensor<1,dim> &v,
                       const unsigned int n)
{
  Assert ((n/2)*2 == n, ExcMessage ("Value of 'n' must be even"));
  double p = 1;
  for (unsigned int k=0; k<n; k+=2)
    p += (v*v);
  return p;
};


template <int dim>
void MinimizationProblem<dim>::assemble_newton_step (const bool p) 
{
  if (p)
    newton_matrix.reinit (sparsity_pattern);
  newton_update.reinit (dof_handler.n_dofs());
  newton_residual.reinit (dof_handler.n_dofs());
  
  QGauss3<dim>  quadrature_formula;

  FEValues<dim> fe_values (fe, quadrature_formula, 
			   UpdateFlags(update_values    |
				       update_gradients |
				       update_q_points  |
				       update_JxW_values));

  const unsigned int   dofs_per_cell = fe.dofs_per_cell;
  const unsigned int   n_q_points    = quadrature_formula.n_quadrature_points;

  FullMatrix<double>   cell_matrix (dofs_per_cell, dofs_per_cell);
  Vector<double>       cell_rhs (dofs_per_cell);

  std::vector<unsigned int> local_dof_indices (dofs_per_cell);
  
  std::vector<double>         local_solution_values (n_q_points);
  std::vector<Tensor<1,dim> > local_solution_grads (n_q_points);
  
  typename DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();
  for (; cell!=endc; ++cell)
    {
      cell_matrix.clear ();
      cell_rhs.clear ();

      fe_values.reinit (cell);

      fe_values.get_function_values (present_solution,
                                     local_solution_values);
      fe_values.get_function_grads (present_solution,
                                    local_solution_grads);
      
      for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
        {
          const double u = local_solution_values[q_point],
                       x = fe_values.quadrature_point(q_point)(0);
          const double x_minus_u3 = (x-std::pow(u,3));

          const Tensor<1,dim> u_prime = local_solution_grads[q_point];
          
          for (unsigned int i=0; i<dofs_per_cell; ++i)
            {
              if (p)
              for (unsigned int j=0; j<dofs_per_cell; ++j)
                {
                  cell_matrix(i,j)
                    += (30.* x_minus_u3 * x_minus_u3 *
                        gradient_power (u_prime, 4) *
                        (fe_values.shape_grad(i,q_point)    *
                         fe_values.shape_grad(j,q_point))   *
                        fe_values.JxW(q_point));

                  cell_matrix(i,j)
                    += (-36. * x_minus_u3 * u * u
                        *
                        gradient_power(u_prime, 4)
                        *                        
                        (fe_values.shape_value(i,q_point)    *
                         (fe_values.shape_grad(j,q_point)   *
                          u_prime)
                         +
                         fe_values.shape_value(j,q_point)    *
                         (fe_values.shape_grad(i,q_point)   *
                          u_prime)
                         )*
                        fe_values.JxW(q_point));

                  cell_matrix(i,j)
                    += ((30.* std::pow(u,4.) - 12*x*u)
                         *
                         gradient_power(u_prime, 6)
                        *                        
                        (fe_values.shape_value(i,q_point)    *
                         fe_values.shape_value(j,q_point))   *
                        fe_values.JxW(q_point));
                };
              
              cell_rhs(i) += -((6. * x_minus_u3 *
                                gradient_power (local_solution_grads[q_point],
                                                4) *
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
            };
        };
      

      cell->get_dof_indices (local_dof_indices);
      for (unsigned int i=0; i<dofs_per_cell; ++i)
	{
          if (p)
	  for (unsigned int j=0; j<dofs_per_cell; ++j)
	    newton_matrix.add (local_dof_indices[i],
			       local_dof_indices[j],
			       cell_matrix(i,j));
	  
	  newton_residual(local_dof_indices[i]) += cell_rhs(i);
	};
    };

  hanging_node_constraints.condense (newton_matrix);
  hanging_node_constraints.condense (newton_residual);

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
  MatrixTools::apply_boundary_values (boundary_values,
				      newton_matrix,
				      newton_update,
				      newton_residual);
//  std::cout << "  res=" << newton_residual.l2_norm() << std::endl;
};



template <int dim>
void MinimizationProblem<dim>::solve () 
{
  SolverControl           solver_control (1000,
                                          1e-3*newton_residual.l2_norm());
  PrimitiveVectorMemory<> vector_memory;
  SolverGMRES<>           gmres (solver_control, vector_memory);

  PreconditionJacobi<> preconditioner;
  preconditioner.initialize(newton_matrix, 1.2);

//   std::cout << "------------------------------" << std::endl;
//   newton_matrix.print_formatted(cout);
//   std::cout << std::endl;
//   for (unsigned int i=0; i<newton_residual.size(); ++i)
//     std::cout << newton_residual(i) << " ";
//   std::cout << std::endl;
  
  
  gmres.solve (newton_matrix, newton_update, newton_residual,
               preconditioner);

  present_solution += newton_update;
  hanging_node_constraints.distribute (present_solution);
};



template <int dim>
void MinimizationProblem<dim>::initialize () 
{
  dof_handler.distribute_dofs (fe);
  present_solution.reinit (dof_handler.n_dofs());
  VectorTools::interpolate (dof_handler,
                            InitializationValues(),
                            present_solution);
};



template <int dim>
void MinimizationProblem<dim>::refine_grid ()
{
  Vector<float> estimated_error_per_cell (triangulation.n_active_cells());

  typename FunctionMap<dim>::type neumann_boundary;

  KellyErrorEstimator<dim>::estimate (dof_handler,
				      Quadrature<dim-1>(1),
				      neumann_boundary,
				      present_solution,
				      estimated_error_per_cell);

  GridRefinement::refine_and_coarsen_fixed_number (triangulation,
						   estimated_error_per_cell,
						   0.3, 0.03);

  SolutionTransfer<dim,double> solution_transfer(dof_handler);
  triangulation.prepare_coarsening_and_refinement();
  solution_transfer.prepare_for_coarsening_and_refinement (present_solution);
  triangulation.execute_coarsening_and_refinement ();
  dof_handler.distribute_dofs (fe);

  Vector<double> tmp (dof_handler.n_dofs());
  solution_transfer.interpolate (present_solution, tmp);
  present_solution = tmp;
};



template <int dim>
double MinimizationProblem<dim>::energy () const
{
  double energy = 0.;

  QGauss3<dim>  quadrature_formula;
  FEValues<dim> fe_values (fe, quadrature_formula, 
			   UpdateFlags(update_values    |
				       update_gradients |
				       update_q_points  |
				       update_JxW_values));

  const unsigned int   n_q_points    = quadrature_formula.n_quadrature_points;

  std::vector<double>         local_solution_values (n_q_points);
  std::vector<Tensor<1,dim> > local_solution_grads (n_q_points);
  
  typename DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();
  for (; cell!=endc; ++cell)
    {
      fe_values.reinit (cell);
      fe_values.get_function_values (present_solution,
                                     local_solution_values);
      fe_values.get_function_grads (present_solution,
                                    local_solution_grads);
      
      for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
        energy += (std::pow (fe_values.quadrature_point(q_point)(0)
                             -
                             std::pow (local_solution_values[q_point],
                                       3),
                             2) *
                   gradient_power (local_solution_grads[q_point],
                                   6) *
                   fe_values.JxW (q_point));
    };
  
  return energy;
};



template <int dim>
void MinimizationProblem<dim>::output_results (const unsigned int cycle) const
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
           << cycle
           << ".gnuplot"
           << std::ends;
#ifdef HAVE_STD_STRINGSTREAM
  std::ofstream out (filename.str().c_str());
#else
  std::ofstream out (filename.str());
#endif

  data_out.write_gnuplot (out);
};



template <int dim>
void MinimizationProblem<dim>::run () 
{
  GridGenerator::hyper_cube (triangulation, 0., 1.);
  triangulation.refine_global (1);
  initialize ();

  for (unsigned int refinement_cycle=0; refinement_cycle<50;
       ++refinement_cycle)
    {
      std::cout << "Cycle " << refinement_cycle << ':' << std::endl;

      std::cout << "   Number of active cells:       "
		<< triangulation.n_active_cells()
		<< std::endl;

      setup_system ();

      unsigned int iteration=0;
      for (; iteration<1000; ++iteration)
        {
          assemble_newton_step (true);
          solve ();

          if (newton_residual.l2_norm() < 1.e-12)
            break;
        };
      output_results (refinement_cycle);

      std::cout << "   Iterations            :       "
		<< iteration
                << std::endl;
      std::cout << "   Energy                :       "
		<< energy ()
		<< std::endl;
      
      refine_grid ();
    };
};

    
int main () 
{
  try
    {
      deallog.depth_console (0);

      MinimizationProblem<1> minimization_problem_1d;
      minimization_problem_1d.run ();
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
    };
  return 0;
};
