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

#include <base/quadrature_lib.h>
#include <base/function.h>
#include <base/logstream.h>
#include <lac/vector.h>
#include <lac/full_matrix.h>
#include <lac/sparse_matrix.h>
                                 //
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
    void output_results (const unsigned int cycle) const;
    
  private:
    void setup_system ();
    void assemble_step ();
    double line_search (const Vector<double> & update) const;
    void do_step ();
    void initialize ();
    void refine_grid ();

    static double energy (const DoFHandler<dim> &dof_handler,
                          const Vector<double>  &function);
 
    
    Triangulation<dim>   triangulation;

    FE_Q<dim>            fe;
    DoFHandler<dim>      dof_handler;

    ConstraintMatrix     hanging_node_constraints;
    
    SparsityPattern      sparsity_pattern;
    SparseMatrix<double> matrix;

    Vector<double>       present_solution;
    Vector<double>       residual;
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
  const double base = std::pow(p(0), 1./3.);
  const double random = 2.*rand()/RAND_MAX-1;
  if (base+.1*random < 0 )
    return 0;
  else
    return base+.1*random;
}



template <int dim>
MinimizationProblem<dim>::MinimizationProblem ()
                :
                fe (1),
		dof_handler (triangulation)
{}


template <int dim>
MinimizationProblem<dim>::~MinimizationProblem () 
{
  dof_handler.clear ();
}



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
}


template <int dim>
double gradient_power (const Tensor<1,dim> &v,
                       const unsigned int n)
{
  Assert ((n/2)*2 == n, ExcMessage ("Value of 'n' must be even"));
  double p = 1;
  for (unsigned int k=0; k<n; k+=2)
    p += (v*v);
  return p;
}


template <int dim>
void MinimizationProblem<dim>::assemble_step ()
{
  matrix.reinit (sparsity_pattern);
  residual.reinit (dof_handler.n_dofs());
  
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
              for (unsigned int j=0; j<dofs_per_cell; ++j)
                cell_matrix(i,j)
                  += (fe_values.shape_grad(i,q_point) *
                      fe_values.shape_grad(j,q_point) * cell->diameter() * cell->diameter() +
                      fe_values.shape_value(i,q_point) *
                      fe_values.shape_value(j,q_point)) *
                  fe_values.JxW(q_point);
              
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
            }
        }
      

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

  hanging_node_constraints.condense (matrix);
  hanging_node_constraints.condense (residual);

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



template <int dim>
double
MinimizationProblem<dim>::line_search (const Vector<double> &update) const
{
  double alpha = 0.;
  Vector<double> tmp (present_solution.size());
  
  for (unsigned int step=0; step<5; ++step)
    {
      tmp = present_solution;
      tmp.add (alpha, update);
      const double f_s = energy (dof_handler, tmp);
      
      const double dalpha = (alpha != 0 ? alpha/100 : 0.01);
      
      tmp = present_solution;
      tmp.add (alpha+dalpha, update);
      const double f_s_plus = energy (dof_handler, tmp);

      tmp = present_solution;
      tmp.add (alpha-dalpha, update);
      const double f_s_minus = energy (dof_handler, tmp);

      const double f_s_prime       = (f_s_plus-f_s_minus) / (2*dalpha);
      const double f_s_doubleprime = ((f_s_plus-2*f_s+f_s_minus) /
                                      (dalpha*dalpha));

      if (std::fabs(f_s_prime) < 1e-7*std::fabs(f_s))
        break;

      if (std::fabs(f_s_doubleprime) < 1e-7*std::fabs(f_s_prime))
        break;

      double step_length = -f_s_prime / f_s_doubleprime;
      for (unsigned int i=0; i<3; ++i)
        {
          tmp = present_solution;
          tmp.add (alpha+step_length, update);
          const double e = energy (dof_handler, tmp);
          
          if (e >= f_s)
            step_length /= 2;
          else
            break;
        }
      alpha += step_length;
    }

  return alpha;
}




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


template <>
void MinimizationProblem<1>::initialize () 
{
  dof_handler.distribute_dofs (fe);
  present_solution.reinit (dof_handler.n_dofs());
  VectorTools::interpolate (dof_handler,
                            InitializationValues(),
                            present_solution);
  DoFHandler<1>::cell_iterator cell;
  cell = dof_handler.begin(0);
  while (cell->has_children())
    cell = cell->child(0);
  present_solution(cell->vertex_dof_index(0,0)) = 0;
  
  cell = dof_handler.begin(0);
  while (cell->has_children())
    cell = cell->child(1);
  present_solution(cell->vertex_dof_index(1,0)) = 1;
}



template <>
void MinimizationProblem<1>::refine_grid ()
{
  const unsigned int dim = 1;
  
  Vector<float> estimated_error_per_cell (triangulation.n_active_cells());

  QTrapez<dim> quadrature;
  FEValues<dim> fe_values (fe, quadrature,
                           update_values   | update_gradients |
                           update_second_derivatives |
                           update_q_points | update_JxW_values);

  FEValues<dim> neighbor_fe_values (fe, quadrature,
                                    update_gradients);

  std::vector<double> local_values (quadrature.n_quadrature_points);
  std::vector<Tensor<1,dim> > local_gradients (quadrature.n_quadrature_points);
  std::vector<Tensor<2,dim> > local_2nd_derivs (quadrature.n_quadrature_points);

  DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active (),
    endc = dof_handler.end ();
  for (unsigned int index = 0; cell!=endc; ++cell, ++index)
    {
      fe_values.reinit (cell);
      fe_values.get_function_values (present_solution, local_values);
      fe_values.get_function_grads (present_solution, local_gradients);
      fe_values.get_function_2nd_derivatives (present_solution, local_2nd_derivs);

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

      estimated_error_per_cell(index) = cell_residual_norm *
                                        cell->diameter() * cell->diameter();

      const double u_left  = local_values[0];
      const double u_right = local_values[1];

      const double u_prime_left  = local_gradients[0][0];
      const double u_prime_right = local_gradients[1][0];

      const double x_left  = fe_values.quadrature_point(0)[0];
      const double x_right = fe_values.quadrature_point(1)[0];

      Assert (x_left  == cell->vertex(0)[0], ExcInternalError());
      Assert (x_right == cell->vertex(1)[0], ExcInternalError());

      if (cell->at_boundary(0) == false)
        {
          DoFHandler<dim>::cell_iterator left_neighbor = cell->neighbor(0);
          while (left_neighbor->has_children())
            left_neighbor = left_neighbor->child(1);
          
          neighbor_fe_values.reinit (left_neighbor);
          neighbor_fe_values.get_function_grads (present_solution, local_gradients);

          const double neighbor_u_prime_left = local_gradients[1][0];

          const double left_jump = std::pow(x_left-std::pow(u_left,3), 2) *
                                   (std::pow(neighbor_u_prime_left,5) -
                                    std::pow(u_prime_left,5));
          estimated_error_per_cell(index) += left_jump * left_jump *
                                             cell->diameter();
        }

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
          estimated_error_per_cell(index) += right_jump * right_jump *
                                             cell->diameter();
        }
      
    } 
  
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

  hanging_node_constraints.clear ();
  DoFTools::make_hanging_node_constraints (dof_handler,
					   hanging_node_constraints);
  hanging_node_constraints.close ();  
  hanging_node_constraints.distribute (present_solution);
}



template <int dim>
double
MinimizationProblem<dim>::energy (const DoFHandler<dim> &dof_handler,
                                  const Vector<double>  &function)
{
  double energy = 0.;

  QGauss3<dim>  quadrature_formula;
  FEValues<dim> fe_values (dof_handler.get_fe(), quadrature_formula, 
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
      fe_values.get_function_values (function,
                                     local_solution_values);
      fe_values.get_function_grads (function,
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
    }
  
  return energy;
}



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
}



template <int dim>
void MinimizationProblem<dim>::run () 
{
  GridGenerator::hyper_cube (triangulation, 0., 1.);
  triangulation.refine_global (4);
  initialize ();

  double last_energy = energy (dof_handler, present_solution);
  
  while (true)
    {
      setup_system ();

      unsigned int iteration=0;
      for (; iteration<5; ++iteration)
        do_step ();

      const double this_energy = energy (dof_handler, present_solution);
      std::cout << "   Energy: " << this_energy << std::endl;

      if ((last_energy-this_energy) < 1e-5*last_energy)
        break;

      last_energy = this_energy;

      refine_grid ();
    }
  std::cout << "   Final Energy: " << energy (dof_handler, present_solution) << std::endl;
  std::cout << std::endl;
}

    
int main () 
{
  try
    {
      deallog.depth_console (0);

      for (unsigned int realization=0; realization<100; ++realization)
        {
          std::cout << "Realization " << realization << ":" << std::endl;
  
          MinimizationProblem<1> minimization_problem_1d;
          minimization_problem_1d.run ();
          minimization_problem_1d.output_results (realization);
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
