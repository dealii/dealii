// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2013 by the deal.II authors
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



// a un-hp-ified version of hp/step-8


#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <fstream>
std::ofstream logfile("output");


#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/hp/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/hp/fe_values.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/numerics/error_estimator.h>

#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_q.h>

#include <fstream>
#include <iomanip>



template <int dim>
class ElasticProblem
{
public:
  ElasticProblem ();
  ~ElasticProblem ();
  void run ();

private:
  void setup_system ();
  void assemble_system ();
  void solve ();
  void refine_grid ();
  void output_results (const unsigned int cycle) const;

  Triangulation<dim>   triangulation;
  DoFHandler<dim>      dof_handler;

  FESystem<dim>        fe;

  ConstraintMatrix     hanging_node_constraints;

  SparsityPattern      sparsity_pattern;
  SparseMatrix<double> system_matrix;

  Vector<double>       solution;
  Vector<double>       system_rhs;
};



template <int dim>
class RightHandSide :  public Function<dim>
{
public:
  RightHandSide ();

  virtual void vector_value (const Point<dim> &p,
                             Vector<double>   &values) const;

  virtual void vector_value_list (const std::vector<Point<dim> > &points,
                                  std::vector<Vector<double> >   &value_list) const;
};


template <int dim>
RightHandSide<dim>::RightHandSide ()
  :
  Function<dim> (dim)
{}


template <int dim>
inline
void RightHandSide<dim>::vector_value (const Point<dim> &p,
                                       Vector<double>   &values) const
{
  Assert (values.size() == dim,
          ExcDimensionMismatch (values.size(), dim));
  Assert (dim >= 2, ExcNotImplemented());

  Point<dim> point_1, point_2;
  point_1(0) = 0.5;
  point_2(0) = -0.5;

  if (((p-point_1).square() < 0.2*0.2) ||
      ((p-point_2).square() < 0.2*0.2))
    values(0) = 1;
  else
    values(0) = 0;

  if (p.square() < 0.2*0.2)
    values(1) = 1;
  else
    values(1) = 0;
}



template <int dim>
void RightHandSide<dim>::vector_value_list (const std::vector<Point<dim> > &points,
                                            std::vector<Vector<double> >   &value_list) const
{
  Assert (value_list.size() == points.size(),
          ExcDimensionMismatch (value_list.size(), points.size()));

  const unsigned int n_points = points.size();

  for (unsigned int p=0; p<n_points; ++p)
    RightHandSide<dim>::vector_value (points[p],
                                      value_list[p]);
}






template <int dim>
ElasticProblem<dim>::ElasticProblem ()
  :
  dof_handler (triangulation),
  fe (FE_Q<dim>(1), dim)
{}




template <int dim>
ElasticProblem<dim>::~ElasticProblem ()
{
  dof_handler.clear ();
}



template <int dim>
void ElasticProblem<dim>::setup_system ()
{
  dof_handler.distribute_dofs (fe);
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

  system_matrix.reinit (sparsity_pattern);

  solution.reinit (dof_handler.n_dofs());
  system_rhs.reinit (dof_handler.n_dofs());
}



template <int dim>
void ElasticProblem<dim>::assemble_system ()
{
  QGauss<dim>  quadrature_formula (2);

  FEValues<dim> x_fe_values (fe, quadrature_formula,
                             update_values   | update_gradients |
                             update_q_points | update_JxW_values);

  const unsigned int   dofs_per_cell = fe.dofs_per_cell;
  const unsigned int   n_q_points    = quadrature_formula.size();

  FullMatrix<double>   cell_matrix (dofs_per_cell, dofs_per_cell);
  Vector<double>       cell_rhs (dofs_per_cell);

  std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);

  std::vector<double>     lambda_values (n_q_points);
  std::vector<double>     mu_values (n_q_points);

  ConstantFunction<dim> lambda(1.), mu(1.);

  RightHandSide<dim>      right_hand_side;
  std::vector<Vector<double> > rhs_values (n_q_points,
                                           Vector<double>(dim));


  typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(),
                                                 endc = dof_handler.end();
  for (; cell!=endc; ++cell)
    {
      cell_matrix = 0;
      cell_rhs = 0;

      x_fe_values.reinit (cell);
      const FEValues<dim> &fe_values = x_fe_values.get_present_fe_values();

      lambda.value_list (fe_values.get_quadrature_points(), lambda_values);
      mu.value_list     (fe_values.get_quadrature_points(), mu_values);

      right_hand_side.vector_value_list (fe_values.get_quadrature_points(),
                                         rhs_values);

      for (unsigned int i=0; i<dofs_per_cell; ++i)
        {
          const unsigned int
          component_i = fe.system_to_component_index(i).first;

          for (unsigned int j=0; j<dofs_per_cell; ++j)
            {
              const unsigned int
              component_j = fe.system_to_component_index(j).first;

              for (unsigned int q_point=0; q_point<n_q_points;
                   ++q_point)
                {
                  cell_matrix(i,j)
                  +=
                    (
                      (fe_values.shape_grad(i,q_point)[component_i] *
                       fe_values.shape_grad(j,q_point)[component_j] *
                       lambda_values[q_point])
                      +
                      (fe_values.shape_grad(i,q_point)[component_j] *
                       fe_values.shape_grad(j,q_point)[component_i] *
                       mu_values[q_point])
                      +
                      ((component_i == component_j) ?
                       (fe_values.shape_grad(i,q_point) *
                        fe_values.shape_grad(j,q_point) *
                        mu_values[q_point])  :
                       0)
                    )
                    *
                    fe_values.JxW(q_point);
                }
            }
        }

      for (unsigned int i=0; i<dofs_per_cell; ++i)
        {
          const unsigned int
          component_i = fe.system_to_component_index(i).first;

          for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
            cell_rhs(i) += fe_values.shape_value(i,q_point) *
                           rhs_values[q_point](component_i) *
                           fe_values.JxW(q_point);
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

  std::map<types::global_dof_index,double> boundary_values;
  VectorTools::interpolate_boundary_values (dof_handler,
                                            0,
                                            ZeroFunction<dim>(dim),
                                            boundary_values);
  MatrixTools::apply_boundary_values (boundary_values,
                                      system_matrix,
                                      solution,
                                      system_rhs);
}




template <int dim>
void ElasticProblem<dim>::solve ()
{
  SolverControl           solver_control (1000, 1e-12);
  SolverCG<>              cg (solver_control);

  PreconditionSSOR<> preconditioner;
  preconditioner.initialize(system_matrix, 1.2);

  cg.solve (system_matrix, solution, system_rhs,
            preconditioner);

  hanging_node_constraints.distribute (solution);
}



template <int dim>
void ElasticProblem<dim>::refine_grid ()
{
  Vector<float> estimated_error_per_cell (triangulation.n_active_cells());

  typename FunctionMap<dim>::type neumann_boundary;
  KellyErrorEstimator<dim>::estimate (dof_handler,
                                      QGauss<dim-1>(2),
                                      neumann_boundary,
                                      solution,
                                      estimated_error_per_cell);

  GridRefinement::refine_and_coarsen_fixed_number (triangulation,
                                                   estimated_error_per_cell,
                                                   0.3, 0.03);

  triangulation.execute_coarsening_and_refinement ();
}



template <int dim>
void ElasticProblem<dim>::output_results (const unsigned int cycle) const
{
  std::string filename = "solution-";
  filename += ('0' + cycle);
  Assert (cycle < 10, ExcInternalError());

  filename += ".gmv";

  DataOut<dim,DoFHandler<dim> > data_out;
  data_out.attach_dof_handler (dof_handler);



  std::vector<std::string> solution_names;
  switch (dim)
    {
    case 1:
      solution_names.push_back ("displacement");
      break;
    case 2:
      solution_names.push_back ("x_displacement");
      solution_names.push_back ("y_displacement");
      break;
    case 3:
      solution_names.push_back ("x_displacement");
      solution_names.push_back ("y_displacement");
      solution_names.push_back ("z_displacement");
      break;
    default:
      Assert (false, ExcNotImplemented());
    }

  data_out.add_data_vector (solution, solution_names);
  data_out.build_patches ();
  data_out.write_gmv (deallog.get_file_stream());
}




template <int dim>
void ElasticProblem<dim>::run ()
{
  for (unsigned int cycle=0; cycle<4; ++cycle)
    {
      deallog << "Cycle " << cycle << ':' << std::endl;

      if (cycle == 0)
        {
          GridGenerator::hyper_cube (triangulation, -1, 1);
          triangulation.refine_global (2);
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
      output_results (cycle);
    }
}


int main ()
{
  deallog << std::setprecision(2);
  logfile << std::setprecision(2);

  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  try
    {
      deallog.depth_console (0);

      ElasticProblem<2> elastic_problem_2d;
      elastic_problem_2d.run ();
    }
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

  return 0;
}
