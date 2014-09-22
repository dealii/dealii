// ---------------------------------------------------------------------
//
// Copyright (C) 2012 - 2014 by the deal.II authors
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


// simplified version of step-48


#include "../tests.h"

#include <deal.II/base/logstream.h>
#include <deal.II/base/utilities.h>
#include <deal.II/base/function.h>
#include <deal.II/lac/vector.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/numerics/vector_tools.h>

#include <deal.II/lac/parallel_vector.h>
#include <deal.II/matrix_free/matrix_free.h>
#include <deal.II/matrix_free/fe_evaluation.h>

#include <fstream>
#include <iostream>
#include <iomanip>


namespace Step48
{
  using namespace dealii;

  const unsigned int fe_degree = 4;



  template <int dim, int fe_degree>
  class SineGordonOperation
  {
  public:
    SineGordonOperation(const MatrixFree<dim,double> &data_in,
                        const double                      time_step);

    void apply (parallel::distributed::Vector<double>                        &dst,
                const std::vector<parallel::distributed::Vector<double>*> &src) const;

  private:
    const MatrixFree<dim,double>         &data;
    const VectorizedArray<double>         delta_t_sqr;
    parallel::distributed::Vector<double> inv_mass_matrix;

    void local_apply (const MatrixFree<dim,double>               &data,
                      parallel::distributed::Vector<double>      &dst,
                      const std::vector<parallel::distributed::Vector<double>*> &src,
                      const std::pair<unsigned int,unsigned int> &cell_range) const;
  };




  template <int dim, int fe_degree>
  SineGordonOperation<dim,fe_degree>::
  SineGordonOperation(const MatrixFree<dim,double> &data_in,
                      const double                  time_step)
    :
    data(data_in),
    delta_t_sqr(make_vectorized_array(time_step *time_step))
  {
    VectorizedArray<double> one = make_vectorized_array (1.);

    data.initialize_dof_vector (inv_mass_matrix);

    FEEvaluation<dim,fe_degree> fe_eval(data);
    const unsigned int          n_q_points = fe_eval.n_q_points;

    for (unsigned int cell=0; cell<data.get_size_info().n_macro_cells; ++cell)
      {
        fe_eval.reinit(cell);
        for (unsigned int q=0; q<n_q_points; ++q)
          fe_eval.submit_value(one,q);
        fe_eval.integrate (true,false);
        fe_eval.distribute_local_to_global (inv_mass_matrix);
      }

    inv_mass_matrix.compress(VectorOperation::add);
    for (unsigned int k=0; k<inv_mass_matrix.local_size(); ++k)
      if (inv_mass_matrix.local_element(k)>1e-15)
        inv_mass_matrix.local_element(k) = 1./inv_mass_matrix.local_element(k);
      else
        inv_mass_matrix.local_element(k) = 0;
  }





  template <int dim, int fe_degree>
  void SineGordonOperation<dim, fe_degree>::
  local_apply (const MatrixFree<dim>                 &data,
               parallel::distributed::Vector<double>     &dst,
               const std::vector<parallel::distributed::Vector<double>*> &src,
               const std::pair<unsigned int,unsigned int> &cell_range) const
  {
    AssertDimension (src.size(), 2);
    FEEvaluation<dim,fe_degree> current (data), old (data);
    deallog << "submit / sine values: ";
    for (unsigned int cell=cell_range.first; cell<cell_range.second; ++cell)
      {
        current.reinit (cell);
        old.reinit (cell);

        current.read_dof_values (*src[0]);
        old.read_dof_values     (*src[1]);

        current.evaluate (true, true, false);
        old.evaluate (true, false, false);

        for (unsigned int q=0; q<current.n_q_points; ++q)
          {
            const VectorizedArray<double> current_value = current.get_value(q);
            const VectorizedArray<double> old_value     = old.get_value(q);

            const VectorizedArray<double> submit_value =
              2.*current_value - old_value - delta_t_sqr * std::sin(current_value);
            const VectorizedArray<double> simple_value =
              2.*current_value - old_value;
            const VectorizedArray<double> sin_value = delta_t_sqr * std::sin(current_value);
            current.submit_value (simple_value-sin_value,q);
            current.submit_gradient (- delta_t_sqr *
                                     current.get_gradient(q), q);

            // output first value on quadrature point for all vector
            // components (should be stable irrespective of vectorization
            // width)
            if (q==0)
              for (unsigned int v=0; v<VectorizedArray<double>::n_array_elements; ++v)
                deallog << submit_value[v] << " = " << simple_value[v] << " - "
                        << sin_value[v] << std::endl;
          }

        current.integrate (true,true);
        current.distribute_local_to_global (dst);
      }
  }




  template <int dim, int fe_degree>
  void SineGordonOperation<dim, fe_degree>::
  apply (parallel::distributed::Vector<double>                        &dst,
         const std::vector<parallel::distributed::Vector<double>*> &src) const
  {
    dst = 0;
    data.cell_loop (&SineGordonOperation<dim,fe_degree>::local_apply,
                    this, dst, src);
    dst.scale(inv_mass_matrix);
  }



  template <int dim>
  class ExactSolution : public Function<dim>
  {
  public:
    ExactSolution (const unsigned int n_components = 1,
                   const double time = 0.) : Function<dim>(n_components, time) {}
    virtual double value (const Point<dim> &p,
                          const unsigned int component = 0) const;
  };

  template <int dim>
  double ExactSolution<dim>::value (const Point<dim> &p,
                                    const unsigned int /* component */) const
  {
    return std::exp(-p.distance(Point<2>(0.03, -0.2))*p.distance(Point<2>(0.03,-0.2))*0.1);
  }




  template <int dim>
  class SineGordonProblem
  {
  public:
    SineGordonProblem ();
    void run ();

  private:
    void make_grid_and_dofs ();
    void output_norm ();

    Triangulation<dim>   triangulation;
    FE_Q<dim>            fe;
    DoFHandler<dim>      dof_handler;
    ConstraintMatrix     constraints;
    IndexSet             locally_relevant_dofs;

    MatrixFree<dim,double> matrix_free_data;

    parallel::distributed::Vector<double> solution, old_solution, old_old_solution;

    const unsigned int n_global_refinements;
    double time, time_step;
    const double final_time;
    const double cfl_number;
    const unsigned int output_timestep_skip;
  };



  template <int dim>
  SineGordonProblem<dim>::SineGordonProblem ()
    :
    fe (QGaussLobatto<1>(fe_degree+1)),
    dof_handler (triangulation),
    n_global_refinements (2),
    time (-10),
    final_time (-9.75),
    cfl_number (.1/fe_degree),
    output_timestep_skip (1)
  {}


  template <int dim>
  void SineGordonProblem<dim>::make_grid_and_dofs ()
  {
    GridGenerator::hyper_cube (triangulation, -15, 15);
    triangulation.refine_global (n_global_refinements);

    deallog << "   Number of global active cells: "
            << triangulation.n_active_cells()
            << std::endl;

    dof_handler.distribute_dofs (fe);

    deallog << "   Number of degrees of freedom: "
            << dof_handler.n_dofs()
            << std::endl;


    DoFTools::extract_locally_relevant_dofs (dof_handler,
                                             locally_relevant_dofs);
    constraints.close();

    QGaussLobatto<1> quadrature (fe_degree+1);
    typename MatrixFree<dim>::AdditionalData additional_data;
    additional_data.tasks_parallel_scheme =
      MatrixFree<dim>::AdditionalData::none;

    matrix_free_data.reinit (dof_handler, constraints,
                             quadrature, additional_data);

    matrix_free_data.initialize_dof_vector (solution);
    old_solution.reinit (solution);
    old_old_solution.reinit (solution);
  }




  template <int dim>
  void
  SineGordonProblem<dim>::output_norm ()
  {
    constraints.distribute (solution);
    solution.update_ghost_values();

    Vector<float> norm_per_cell (triangulation.n_active_cells());
    VectorTools::integrate_difference (dof_handler,
                                       solution,
                                       ZeroFunction<dim>(),
                                       norm_per_cell,
                                       QGauss<dim>(fe_degree+1),
                                       VectorTools::L2_norm);
    const double solution_norm = norm_per_cell.norm_sqr();

    deallog << "   Time: "
            << std::setw(6) << std::setprecision(3) << time
            << ", solution norm: "
            << std::setprecision(9) << std::setw(13) << solution_norm
            << std::endl;
  }



  template <int dim>
  void
  SineGordonProblem<dim>::run ()
  {
    make_grid_and_dofs();

    const double min_cell_diameter =
      triangulation.last()->diameter()/std::sqrt(dim);
    time_step = cfl_number * min_cell_diameter;
    time_step = (final_time-time)/(int((final_time-time)/time_step));
    deallog << "   Time step size: " << time_step << ", finest cell: "
            << min_cell_diameter << std::endl << std::endl;


    VectorTools::interpolate (dof_handler,
                              ExactSolution<dim> (1, time),
                              solution);
    VectorTools::interpolate (dof_handler,
                              ExactSolution<dim> (1, time-time_step),
                              old_solution);
    output_norm ();

    std::vector<parallel::distributed::Vector<double>*> previous_solutions;
    previous_solutions.push_back(&old_solution);
    previous_solutions.push_back(&old_old_solution);

    SineGordonOperation<dim,fe_degree> sine_gordon_op (matrix_free_data,
                                                       time_step);

    unsigned int timestep_number = 1;

    for (time+=time_step; time<=final_time; time+=time_step, ++timestep_number)
      {
        old_old_solution.swap (old_solution);
        old_solution.swap (solution);
        sine_gordon_op.apply (solution, previous_solutions);

        output_norm();
      }

    deallog << std::endl
            << "   Performed " << timestep_number << " time steps."
            << std::endl << std::endl;
  }
}



int main (int argc, char **argv)
{
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog << std::setprecision(4);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  Step48::SineGordonProblem<2> sg_problem;
  sg_problem.run ();
}

