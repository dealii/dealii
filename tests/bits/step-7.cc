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



// a un-hp-ified version of hp/step-7


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
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/data_out.h>

#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/base/smartpointer.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/base/convergence_table.h>
#include <deal.II/hp/fe_values.h>

#include <typeinfo>
#include <fstream>
#include <iomanip>



template <int dim>
class SolutionBase
{
protected:
  static const unsigned int n_source_centers = 3;
  static const Point<dim>   source_centers[n_source_centers];
  static const double       width;
};


template <>
const Point<1>
SolutionBase<1>::source_centers[SolutionBase<1>::n_source_centers]
  = { Point<1>(-1.0 / 3.0),
      Point<1>(0.0),
      Point<1>(+1.0 / 3.0)
    };

template <>
const Point<2>
SolutionBase<2>::source_centers[SolutionBase<2>::n_source_centers]
  = { Point<2>(-0.5, +0.5),
      Point<2>(-0.5, -0.5),
      Point<2>(+0.5, -0.5)
    };

template <int dim>
const double SolutionBase<dim>::width = 1./3.;



template <int dim>
class Solution : public Function<dim>,
  protected SolutionBase<dim>
{
public:
  Solution () : Function<dim>() {}

  virtual double value (const Point<dim>   &p,
                        const unsigned int  component = 0) const;

  virtual Tensor<1,dim> gradient (const Point<dim>   &p,
                                  const unsigned int  component = 0) const;
};


template <int dim>
double Solution<dim>::value (const Point<dim>   &p,
                             const unsigned int) const
{
  double return_value = 0;
  for (unsigned int i=0; i<this->n_source_centers; ++i)
    {
      const Point<dim> x_minus_xi = p - this->source_centers[i];
      return_value += std::exp(-x_minus_xi.square() /
                               (this->width * this->width));
    }

  return return_value;
}


template <int dim>
Tensor<1,dim> Solution<dim>::gradient (const Point<dim>   &p,
                                       const unsigned int) const
{
  Tensor<1,dim> return_value;

  for (unsigned int i=0; i<this->n_source_centers; ++i)
    {
      const Point<dim> x_minus_xi = p - this->source_centers[i];

      return_value += (-2 / (this->width * this->width) *
                       std::exp(-x_minus_xi.square() /
                                (this->width * this->width)) *
                       x_minus_xi);
    }

  return return_value;
}



template <int dim>
class RightHandSide : public Function<dim>,
  protected SolutionBase<dim>
{
public:
  RightHandSide () : Function<dim>() {}

  virtual double value (const Point<dim>   &p,
                        const unsigned int  component = 0) const;
};


template <int dim>
double RightHandSide<dim>::value (const Point<dim>   &p,
                                  const unsigned int) const
{
  double return_value = 0;
  for (unsigned int i=0; i<this->n_source_centers; ++i)
    {
      const Point<dim> x_minus_xi = p - this->source_centers[i];

      return_value += ((2*dim - 4*x_minus_xi.square()/
                        (this->width * this->width)) /
                       (this->width * this->width) *
                       std::exp(-x_minus_xi.square() /
                                (this->width * this->width)));
      return_value += std::exp(-x_minus_xi.square() /
                               (this->width * this->width));
    }

  return return_value;
}



template <int dim>
class HelmholtzProblem
{
public:
  enum RefinementMode
  {
    global_refinement, adaptive_refinement
  };

  HelmholtzProblem (const FiniteElement<dim> &fe,
                    const RefinementMode      refinement_mode);

  ~HelmholtzProblem ();

  void run ();

private:
  void setup_system ();
  void assemble_system ();
  void solve ();
  void refine_grid ();
  void process_solution (const unsigned int cycle);

  Triangulation<dim>                      triangulation;
  DoFHandler<dim>                         dof_handler;

  SmartPointer<const FiniteElement<dim> > fe;

  ConstraintMatrix                        hanging_node_constraints;

  SparsityPattern                         sparsity_pattern;
  SparseMatrix<double>                    system_matrix;

  Vector<double>                          solution;
  Vector<double>                          system_rhs;

  const RefinementMode                    refinement_mode;

  ConvergenceTable                        convergence_table;
};




template <int dim>
HelmholtzProblem<dim>::HelmholtzProblem (const FiniteElement<dim> &fe,
                                         const RefinementMode refinement_mode) :
  dof_handler (triangulation),
  fe (&fe),
  refinement_mode (refinement_mode)
{}



template <int dim>
HelmholtzProblem<dim>::~HelmholtzProblem ()
{
  dof_handler.clear ();
}



template <int dim>
void HelmholtzProblem<dim>::setup_system ()
{
  dof_handler.distribute_dofs (*fe);
  DoFRenumbering::Cuthill_McKee (dof_handler);

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
void HelmholtzProblem<dim>::assemble_system ()
{
  QGauss<dim>   quadrature_formula (3);
  QGauss<dim-1> face_quadrature_formula (3);

  const unsigned int n_q_points    = quadrature_formula.size();
  const unsigned int n_face_q_points = face_quadrature_formula.size();

  const unsigned int dofs_per_cell = (*fe).dofs_per_cell;

  FullMatrix<double>   cell_matrix (dofs_per_cell, dofs_per_cell);
  Vector<double>       cell_rhs (dofs_per_cell);

  std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);

  FEValues<dim>  x_fe_values (*fe, quadrature_formula,
                              update_values   | update_gradients |
                              update_q_points | update_JxW_values);

  FEFaceValues<dim> x_fe_face_values (*fe, face_quadrature_formula,
                                      update_values         | update_q_points  |
                                      update_normal_vectors | update_JxW_values);

  const RightHandSide<dim> right_hand_side;
  std::vector<double>  rhs_values (n_q_points);

  const Solution<dim> exact_solution;

  typename DoFHandler<dim>::active_cell_iterator
  cell = dof_handler.begin_active(),
  endc = dof_handler.end();
  for (; cell!=endc; ++cell)
    {
      cell_matrix = 0;
      cell_rhs = 0;

      x_fe_values.reinit (cell);
      const FEValues<dim> &fe_values = x_fe_values.get_present_fe_values();

      right_hand_side.value_list (fe_values.get_quadrature_points(),
                                  rhs_values);

      for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
        for (unsigned int i=0; i<dofs_per_cell; ++i)
          {
            for (unsigned int j=0; j<dofs_per_cell; ++j)
              cell_matrix(i,j) += ((fe_values.shape_grad(i,q_point) *
                                    fe_values.shape_grad(j,q_point)
                                    +
                                    fe_values.shape_value(i,q_point) *
                                    fe_values.shape_value(j,q_point)) *
                                   fe_values.JxW(q_point));

            cell_rhs(i) += (fe_values.shape_value(i,q_point) *
                            rhs_values [q_point] *
                            fe_values.JxW(q_point));
          }

      for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
        if (cell->face(face)->at_boundary()
            &&
            (cell->face(face)->boundary_indicator() == 1))
          {
            x_fe_face_values.reinit (cell, face);
            const FEFaceValues<dim> &fe_face_values
              = x_fe_face_values.get_present_fe_values ();

            for (unsigned int q_point=0; q_point<n_face_q_points; ++q_point)
              {
                const double neumann_value
                  = (exact_solution.gradient (fe_face_values.quadrature_point(q_point)) *
                     fe_face_values.normal_vector(q_point));

                for (unsigned int i=0; i<dofs_per_cell; ++i)
                  cell_rhs(i) += (neumann_value *
                                  fe_face_values.shape_value(i,q_point) *
                                  fe_face_values.JxW(q_point));
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

  std::map<types::global_dof_index,double> boundary_values;
  VectorTools::interpolate_boundary_values (dof_handler,
                                            0,
                                            Solution<dim>(),
                                            boundary_values);
  MatrixTools::apply_boundary_values (boundary_values,
                                      system_matrix,
                                      solution,
                                      system_rhs);
}



template <int dim>
void HelmholtzProblem<dim>::solve ()
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
void HelmholtzProblem<dim>::refine_grid ()
{
  switch (refinement_mode)
    {
    case global_refinement:
    {
      triangulation.refine_global (1);
      break;
    }

    case adaptive_refinement:
    {
      Vector<float> estimated_error_per_cell (triangulation.n_active_cells());

      typename FunctionMap<dim>::type neumann_boundary;
      KellyErrorEstimator<dim>::estimate (dof_handler,
                                          QGauss<dim-1>(3),
                                          neumann_boundary,
                                          solution,
                                          estimated_error_per_cell);

      GridRefinement::refine_and_coarsen_fixed_number (triangulation,
                                                       estimated_error_per_cell,
                                                       0.3, 0.03);

      triangulation.execute_coarsening_and_refinement ();

      break;
    }

    default:
    {
      Assert (false, ExcNotImplemented());
    }
    }
}



template <int dim>
void HelmholtzProblem<dim>::process_solution (const unsigned int cycle)
{
  Vector<float> difference_per_cell (triangulation.n_active_cells());
  VectorTools::integrate_difference (dof_handler,
                                     solution,
                                     Solution<dim>(),
                                     difference_per_cell,
                                     QGauss<dim>(3),
                                     VectorTools::L2_norm);
  const double L2_error = difference_per_cell.l2_norm();

  VectorTools::integrate_difference (dof_handler,
                                     solution,
                                     Solution<dim>(),
                                     difference_per_cell,
                                     QGauss<dim>(3),
                                     VectorTools::H1_seminorm);
  const double H1_error = difference_per_cell.l2_norm();

  const QTrapez<1>     q_trapez;
  const QIterated<dim> q_iterated (q_trapez, 5);
  VectorTools::integrate_difference (dof_handler,
                                     solution,
                                     Solution<dim>(),
                                     difference_per_cell,
                                     q_iterated,
                                     VectorTools::Linfty_norm);
  const double Linfty_error = difference_per_cell.linfty_norm();

  const unsigned int n_active_cells=triangulation.n_active_cells();
  const unsigned int n_dofs=dof_handler.n_dofs();

  deallog << "Cycle " << cycle << ':'
          << std::endl
          << "   Number of active cells:       "
          << n_active_cells
          << std::endl
          << "   Number of degrees of freedom: "
          << n_dofs
          << std::endl;

  convergence_table.add_value("cycle", cycle);
  convergence_table.add_value("cells", n_active_cells);
  convergence_table.add_value("dofs", n_dofs);
  convergence_table.add_value("L2", L2_error);
  convergence_table.add_value("H1", H1_error);
  convergence_table.add_value("Linfty", Linfty_error);
}



template <int dim>
void HelmholtzProblem<dim>::run ()
{
  for (unsigned int cycle=0; cycle<4; ++cycle)
    {
      if (cycle == 0)
        {
          GridGenerator::hyper_cube (triangulation, -1, 1);
          triangulation.refine_global (1);

          typename Triangulation<dim>::cell_iterator
          cell = triangulation.begin (),
          endc = triangulation.end();
          for (; cell!=endc; ++cell)
            for (unsigned int face=0;
                 face<GeometryInfo<dim>::faces_per_cell;
                 ++face)
              if ((cell->face(face)->center()(0) == -1)
                  ||
                  (cell->face(face)->center()(1) == -1))
                cell->face(face)->set_boundary_indicator (1);
        }
      else
        refine_grid ();


      setup_system ();

      assemble_system ();
      solve ();

      process_solution (cycle);
    }


  if (refinement_mode == adaptive_refinement)
    {
      std::string gmv_filename;
      switch (refinement_mode)
        {
        case global_refinement:
          gmv_filename = "solution-global";
          break;
        case adaptive_refinement:
          gmv_filename = "solution-adaptive";
          break;
        default:
          Assert (false, ExcNotImplemented());
        }

      switch ((*fe).degree)
        {
        case 1:
          gmv_filename += "-q1";
          break;
        case 2:
          gmv_filename += "-q2";
          break;

        default:
          Assert (false, ExcNotImplemented());
        }

      gmv_filename += ".gmv";

      DataOut<dim,DoFHandler<dim> > data_out;
      data_out.attach_dof_handler (dof_handler);
      data_out.add_data_vector (solution, "solution");

      data_out.build_patches ((*fe).degree);
      data_out.write_gmv (deallog.get_file_stream());
    }



  convergence_table.set_precision("L2", 3);
  convergence_table.set_precision("H1", 3);
  convergence_table.set_precision("Linfty", 3);

  convergence_table.set_scientific("L2", true);
  convergence_table.set_scientific("H1", true);
  convergence_table.set_scientific("Linfty", true);

  convergence_table.set_tex_caption("cells", "\\# cells");
  convergence_table.set_tex_caption("dofs", "\\# dofs");
  convergence_table.set_tex_caption("L2", "$L^2$-error");
  convergence_table.set_tex_caption("H1", "$H^1$-error");
  convergence_table.set_tex_caption("Linfty", "$L^\\infty$-error");

  convergence_table.set_tex_format("cells", "r");
  convergence_table.set_tex_format("dofs", "r");

  deallog << std::endl;
  convergence_table.write_text(deallog.get_file_stream());

  std::string error_filename = "error";
  switch (refinement_mode)
    {
    case global_refinement:
      error_filename += "-global";
      break;
    case adaptive_refinement:
      error_filename += "-adaptive";
      break;
    default:
      Assert (false, ExcNotImplemented());
    }

  switch ((*fe).degree)
    {
    case 1:
      error_filename += "-q1";
      break;
    case 2:
      error_filename += "-q2";
      break;
    default:
      Assert (false, ExcNotImplemented());
    }

  error_filename += ".tex";

  convergence_table.write_tex(deallog.get_file_stream());



  if (refinement_mode==global_refinement)
    {
      convergence_table.add_column_to_supercolumn("cycle", "n cells");
      convergence_table.add_column_to_supercolumn("cells", "n cells");

      std::vector<std::string> new_order;
      new_order.push_back("n cells");
      new_order.push_back("H1");
      new_order.push_back("L2");
      convergence_table.set_column_order (new_order);

      convergence_table
      .evaluate_convergence_rates("L2", ConvergenceTable::reduction_rate);
      convergence_table
      .evaluate_convergence_rates("L2", ConvergenceTable::reduction_rate_log2);
      convergence_table
      .evaluate_convergence_rates("H1", ConvergenceTable::reduction_rate_log2);

      deallog << std::endl;
      convergence_table.write_text(deallog.get_file_stream());

      std::string conv_filename = "convergence";
      switch (refinement_mode)
        {
        case global_refinement:
          conv_filename += "-global";
          break;
        case adaptive_refinement:
          conv_filename += "-adaptive";
          break;
        default:
          Assert (false, ExcNotImplemented());
        }
      switch ((*fe).degree)
        {
        case 1:
          conv_filename += "-q1";
          break;
        case 2:
          conv_filename += "-q2";
          break;
        default:
          Assert (false, ExcNotImplemented());
        }
      conv_filename += ".tex";

      convergence_table.write_tex(deallog.get_file_stream());
    }
}


int main ()
{
  deallog << std::setprecision(2);
  logfile << std::setprecision(2);

  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  const unsigned int dim = 2;

  try
    {
      deallog.depth_console (0);


      {
        deallog << "Solving with Q1 elements, adaptive refinement" << std::endl
                << "=============================================" << std::endl
                << std::endl;

        FE_Q<dim> fe(1);
        HelmholtzProblem<dim>
        helmholtz_problem_2d (fe, HelmholtzProblem<dim>::adaptive_refinement);

        helmholtz_problem_2d.run ();

        deallog << std::endl;
      }

      {
        deallog << "Solving with Q1 elements, global refinement" << std::endl
                << "===========================================" << std::endl
                << std::endl;

        FE_Q<dim> fe(1);
        HelmholtzProblem<dim>
        helmholtz_problem_2d (fe, HelmholtzProblem<dim>::global_refinement);

        helmholtz_problem_2d.run ();

        deallog << std::endl;
      }

      {
        deallog << "Solving with Q2 elements, global refinement" << std::endl
                << "===========================================" << std::endl
                << std::endl;

        FE_Q<dim> fe(2);
        HelmholtzProblem<dim>
        helmholtz_problem_2d (fe, HelmholtzProblem<dim>::global_refinement);

        helmholtz_problem_2d.run ();

        deallog << std::endl;
      }

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


template const double SolutionBase<2>::width;
