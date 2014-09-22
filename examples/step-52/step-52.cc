/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2014 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * The deal.II library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE at
 * the top level of the deal.II distribution.
 *
 * ---------------------------------------------------------------------

 *
 * Authors: Damien Lebrun-Grandie, Bruno Turcksin, 2014
 */

// @sect3{Include files}

// The first task as usual is to include the functionality of these well-known
// deal.II library files and some C++ header files.
#include <deal.II/base/function.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_out.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/lac/sparse_direct.h>

#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/data_out.h>

#include <fstream>
#include <iostream>
#include <cmath>
#include <map>

// This is the only include file that is new: It includes all the Runge-Kutta
// methods.
#include <deal.II/base/time_stepping.h>


// The next step is like in all previous tutorial programs: We put everything
// into a namespace of its own and then import the deal.II classes and functions
// into it.
namespace Step52
{
  using namespace dealii;

  // @sect3{Diffusion}

  // The next piece is the declaration of the main class. Most of the functions in
  // this class are not new and have been explained in previous tutorials. The
  // only interesting functions are <code>evaluate_diffusion</code> and
  // <code>id_minus_tau_J_inverse</code>. <code>evaluate_diffusion</code> evaluates the
  // diffusion equation, $M^{-1}(f(t,y))$, at a given time, for a given $\tau$
  // and a given $y$. <code>id_minus_tau_J_inverse</code> evaluates
  // $\left(I-\tau M^{-1} \frac{\partial f(t,y)}{\partial y}\right)^{-1}$ or
  // equivalently $\left(M-\tau \frac{\partial f}{\partial y}\right)^{-1} M$ at
  // a given time, for a given $\tau$ and $y$. This function is needed when an
  // implicit method is used.
  class Diffusion
  {
  public:
    Diffusion();

    void run();

  private:
    void setup_system();

    void assemble_system();

    double get_source(double time,const Point<2> &point) const;

    Vector<double> evaluate_diffusion(const double time, const Vector<double> &y) const;

    Vector<double> id_minus_tau_J_inverse(const double time,
                                          const double tau,
                                          const Vector<double> &y);

    void output_results(unsigned int time_step,TimeStepping::runge_kutta_method method) const;

    // The next three functions are the driver for the explicit methods, the
    // implicit methods, and the embedded explicit methods respectively. The
    // driver function for embedded explicit methods returns the number of
    // steps executed.
    void explicit_method(TimeStepping::runge_kutta_method method,
                         const unsigned int n_time_steps,
                         const double       initial_time,
                         const double       final_time);

    void implicit_method(TimeStepping::runge_kutta_method method,
                         const unsigned int n_time_steps,
                         const double       initial_time,
                         const double       final_time);

    unsigned int embedded_explicit_method(TimeStepping::runge_kutta_method method,
                                          const unsigned int n_time_steps,
                                          const double initial_time,
                                          const double final_time);


    unsigned int                 fe_degree;

    double                       diffusion_coefficient;
    double                       absorption_xs;

    Triangulation<2>             triangulation;

    FE_Q<2>                      fe;

    DoFHandler<2>                dof_handler;

    ConstraintMatrix             constraint_matrix;

    SparsityPattern              sparsity_pattern;

    SparseMatrix<double>         system_matrix;
    SparseMatrix<double>         mass_matrix;
    SparseMatrix<double>         mass_minus_tau_Jacobian;

    SparseDirectUMFPACK          inverse_mass_matrix;

    Vector<double>               solution;
  };



  // We choose quadratic finite elements and we initialize the parameters.
  Diffusion::Diffusion()
    :
    fe_degree(2),
    diffusion_coefficient(1./30.),
    absorption_xs(1.),
    fe(fe_degree),
    dof_handler(triangulation)
  {}



  // @sect5{<code>Diffusion::setup_system</code>}
  // Now, we create the constraint matrix and the sparsity pattern. Then, we
  // initialize the matrices and the solution vector.
  void Diffusion::setup_system()
  {
    dof_handler.distribute_dofs(fe);

    VectorTools::interpolate_boundary_values(dof_handler,1,ZeroFunction<2>(),constraint_matrix);
    constraint_matrix.close();

    CompressedSparsityPattern c_sparsity(dof_handler.n_dofs());
    DoFTools::make_sparsity_pattern(dof_handler,c_sparsity,constraint_matrix);
    sparsity_pattern.copy_from(c_sparsity);

    system_matrix.reinit(sparsity_pattern);
    mass_matrix.reinit(sparsity_pattern);
    mass_minus_tau_Jacobian.reinit(sparsity_pattern);
    solution.reinit(dof_handler.n_dofs());
  }



  // @sect5{<code>Diffusion::assemble_system</code>}
  // In this function, we compute
  // $-\int D \nabla b_i \cdot \nabla b_j d\boldsymbol{r} - \int \Sigma_a b_i b_j d\boldsymbol{r}$
  // and the mass matrix $\int b_i b_j d\boldsymbol{r}$. The mass matrix is then
  // inverted using a direct solver.
  void Diffusion::assemble_system()
  {
    system_matrix = 0.;
    mass_matrix = 0.;

    const QGauss<2> quadrature_formula(fe_degree+1);

    FEValues<2> fe_values(fe, quadrature_formula,
                          update_values | update_gradients | update_JxW_values);


    const unsigned int dofs_per_cell = fe.dofs_per_cell;
    const unsigned int n_q_points    = quadrature_formula.size();

    FullMatrix<double> cell_matrix (dofs_per_cell, dofs_per_cell);
    FullMatrix<double> cell_mass_matrix (dofs_per_cell, dofs_per_cell);

    std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);

    typename DoFHandler<2>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();

    for (; cell!=endc; ++cell)
      {
        cell_matrix = 0.;
        cell_mass_matrix = 0.;

        fe_values.reinit (cell);

        for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
          for (unsigned int i=0; i<dofs_per_cell; ++i)
            for (unsigned int j=0; j<dofs_per_cell; ++j)
              {
                cell_matrix(i,j) += ((-diffusion_coefficient * fe_values.shape_grad(i,q_point) *
                                      fe_values.shape_grad(j,q_point) - absorption_xs *
                                      fe_values.shape_value(i,q_point) * fe_values.shape_value(j,q_point)) *
                                     fe_values.JxW(q_point));
                cell_mass_matrix(i,j) +=  fe_values.shape_value(i,q_point) *
                                          fe_values.shape_value(j,q_point) *
                                          fe_values.JxW(q_point);
              }

        cell->get_dof_indices(local_dof_indices);

        constraint_matrix.distribute_local_to_global(cell_matrix,local_dof_indices,system_matrix);
        constraint_matrix.distribute_local_to_global(cell_mass_matrix,local_dof_indices,mass_matrix);
      }

    inverse_mass_matrix.initialize(mass_matrix);
  }



  // @sect5{<code>Diffusion::get_source</code>}
  //
  // In this function, the source for a given time and a given point is
  // computed.
  double Diffusion::get_source(double time,const Point<2> &point) const
  {
    const double pi = 3.14159265358979323846;
    const double intensity = 10.;
    const double frequency = pi/10.;
    const double b = 5.;
    const double x = point(0);
    double source = 0.;

    source = intensity*(frequency*std::cos(frequency*time)*(b*x-x*x) + std::sin(frequency*time) *
                        (absorption_xs*(b*x-x*x)+2.*diffusion_coefficient));

    return source;
  }



  // @sect5{<code>Diffusion:evaluate_diffusion</code>}
  //
  // Now, the weak form of the diffusion equation is evaluated at a given
  // time $t$ and for a given vector $y$.
  Vector<double> Diffusion::evaluate_diffusion(const double time, const Vector<double> &y) const
  {
    Vector<double> tmp(dof_handler.n_dofs());
    tmp = 0.;
    system_matrix.vmult(tmp,y);

    const QGauss<2> quadrature_formula(fe_degree+1);

    FEValues<2> fe_values(fe, quadrature_formula,
                          update_values | update_quadrature_points | update_JxW_values);


    const unsigned int dofs_per_cell = fe.dofs_per_cell;
    const unsigned int n_q_points    = quadrature_formula.size();

    Vector<double>  cell_source(dofs_per_cell);

    std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);

    typename DoFHandler<2>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();

    for (; cell!=endc; ++cell)
      {
        cell_source = 0.;

        fe_values.reinit (cell);

        for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
          {
            double source = get_source(time,fe_values.quadrature_point(q_point)) ;
            for (unsigned int i=0; i<dofs_per_cell; ++i)
              cell_source(i) += source * fe_values.shape_value(i,q_point) *
                                fe_values.JxW(q_point);
          }

        cell->get_dof_indices(local_dof_indices);

        constraint_matrix.distribute_local_to_global(cell_source,local_dof_indices,tmp);
      }


    Vector<double> value(dof_handler.n_dofs());
    inverse_mass_matrix.vmult(value,tmp);

    return value;
  }


  // @sect5{<code>Diffusion::id_minus_tau_J_inverse</code>}
  //
  // We compute $\left(M-\tau \frac{\partial f}{\partial y}\right)^{-1} M$. This
  // is done in several steps:
  //   - compute $M-\tau \frac{\partial f}{\partial y}$
  //   - inverse the matrix to get $\left(M-\tau \frac{\partial f}{\partial y}\right)^{-1}$
  //   - compute $tmp=My$
  //   - compute $\left(M-\tau \frac{\partial f}{\partial y}\right)^{-1} tmp = \left(M-\tau \frac{\partial f}{\partial y}\right)^{-1} My$.
  Vector<double> Diffusion::id_minus_tau_J_inverse(const double time, const double tau,
                                                   const Vector<double> &y)
  {
    Vector<double> tmp(dof_handler.n_dofs());
    Vector<double> result(y);
    SparseDirectUMFPACK inverse_mass_minus_tau_Jacobian;

    mass_minus_tau_Jacobian.copy_from(mass_matrix);
    mass_minus_tau_Jacobian.add(-tau,system_matrix);

    inverse_mass_minus_tau_Jacobian.initialize(mass_minus_tau_Jacobian);

    mass_matrix.vmult(tmp,y);

    inverse_mass_minus_tau_Jacobian.vmult(result,tmp);

    return result;
  }



  // @sect5{<code>Diffusion::output_results</code>}
  //
  // We output the solution in vtu files.
  void Diffusion::output_results(unsigned int time_step,TimeStepping::runge_kutta_method method) const
  {
    std::string method_name;

    switch (method)
      {
      case TimeStepping::FORWARD_EULER :
      {
        method_name = "forward_euler";
        break;
      }
      case TimeStepping::RK_THIRD_ORDER :
      {
        method_name = "rk3";
        break;
      }
      case TimeStepping::RK_CLASSIC_FOURTH_ORDER :
      {
        method_name = "rk4";
        break;
      }
      case TimeStepping::BACKWARD_EULER :
      {
        method_name = "backward_euler";
        break;
      }
      case TimeStepping::IMPLICIT_MIDPOINT :
      {
        method_name = "implicit_midpoint";
        break;
      }
      case TimeStepping::SDIRK_TWO_STAGES :
      {
        method_name = "sdirk";
        break;
      }
      case TimeStepping::HEUN_EULER :
      {
        method_name = "heun_euler";
        break;
      }
      case TimeStepping::BOGACKI_SHAMPINE :
      {
        method_name = "bocacki_shampine";
        break;
      }
      case TimeStepping::DOPRI :
      {
        method_name = "dopri";
        break;
      }
      case TimeStepping::FEHLBERG :
      {
        method_name = "fehlberg";
        break;
      }
      case TimeStepping::CASH_KARP :
      {
        method_name = "cash_karp";
        break;
      }
      default :
      {
        break;
      }
      }

    DataOut<2> data_out;

    data_out.attach_dof_handler(dof_handler);
    data_out.add_data_vector(solution, "solution");

    data_out.build_patches();

    const std::string filename = "solution-" + method_name + "-" +
                                 Utilities::int_to_string (time_step, 3) +
                                 ".vtu";
    std::ofstream output(filename.c_str());
    data_out.write_vtu(output);
  }


  // @sect5{<code>Diffusion::explicit_method</code>}
  // This function is the driver for all the explicit method. It calls
  // <code>evolve_one_time_step</code> which performs one time step.
  // <code>evolve_one_time_step</code> needs to evaluate $M^{-1}(f(t,y))$,
  // i.e, it needs <code>evaluate_diffusion</code>. Because
  // <code>evaluate_diffusion</code> is a member function, it needs to be bind to
  // $this$. Finally, the solution is output every 10 time steps.
  void Diffusion::explicit_method(TimeStepping::runge_kutta_method method,
                                  const unsigned int                n_time_steps,
                                  const double                      initial_time,
                                  const double                      final_time)
  {
    const double time_step = (final_time-initial_time)/static_cast<double> (n_time_steps);
    double time = initial_time;
    solution = 0.;

    TimeStepping::ExplicitRungeKutta<Vector<double> > explicit_runge_kutta(method);
    output_results(0,method);
    for (unsigned int i=0; i<n_time_steps; ++i)
      {
        time = explicit_runge_kutta.evolve_one_time_step(
                 std_cxx11::bind(&Diffusion::evaluate_diffusion,this,std_cxx11::_1,std_cxx11::_2),
                 time,time_step,solution);

        if ((i+1)%10==0)
          output_results(i+1,method);
      }
  }



  // @sect5{<code>Diffusion::implicit_method</code>}
  // This function is equivalent to <code>explicit_method</code> but for implicit
  // methods. When using implicit methods, we need to evaluate $M^{-1}(f(t,y))$
  // and $\left(I-\tau M^{-1} \frac{\partial f(t,y)}{\partial y}\right)^{-1}$.
  void Diffusion::implicit_method(TimeStepping::runge_kutta_method method,
                                  const unsigned int               n_time_steps,
                                  const double                     initial_time,
                                  const double                     final_time)
  {
    const double time_step = (final_time-initial_time)/static_cast<double> (n_time_steps);
    double time = initial_time;
    solution = 0.;

    TimeStepping::ImplicitRungeKutta<Vector<double> > implicit_runge_kutta(method);
    output_results(0,method);
    for (unsigned int i=0; i<n_time_steps; ++i)
      {
        time = implicit_runge_kutta.evolve_one_time_step(
                 std_cxx11::bind(&Diffusion::evaluate_diffusion,this,std_cxx11::_1,std_cxx11::_2),
                 std_cxx11::bind(&Diffusion::id_minus_tau_J_inverse,this,std_cxx11::_1,std_cxx11::_2,
                                 std_cxx11::_3),
                 time,time_step,solution);

        if ((i+1)%10==0)
          output_results(i+1,method);
      }
  }



  // @sect5{<code>Diffusion::embedded_explicit_method</code>}
  // This function is the driver for the embedded explict methods. It requires
  // more parameters:
  //   - coarsen_param: factor multiplying the current time step when the error
  //   is below the threshold.
  //   - refine_param: factor multiplying the current time step when the error
  //   is above the threshold.
  //   - min_delta: smallest time step acceptable.
  //   - max_delta: largest time step acceptable.
  //   - refine_tol: threshold above which the time step is refined.
  //   - coarsen_tol: threshold below which the time step is coarsen.
  // Embedded methods use a guessed time step. If the error using this time step
  // is too large, the time step will be reduced. If the error is below the
  // threshold, a larger time step will be tried for the next time step.
  // </code>delta_t_guess</code> is the guessed time step produced by the embedded method.
  unsigned int Diffusion::embedded_explicit_method(TimeStepping::runge_kutta_method method,
                                                   const unsigned int n_time_steps,
                                                   const double initial_time,
                                                   const double final_time)
  {
    double time_step = (final_time-initial_time)/static_cast<double> (n_time_steps);
    double time = initial_time;
    const double coarsen_param = 1.2;
    const double refine_param = 0.8;
    const double min_delta = 1e-8;
    const double max_delta = 10*time_step;
    const double refine_tol = 1e-1;
    const double coarsen_tol = 1e-5;
    solution = 0.;

    TimeStepping::EmbeddedExplicitRungeKutta<Vector<double> > embedded_explicit_runge_kutta(method,
        coarsen_param,refine_param,min_delta,max_delta,refine_tol,coarsen_tol);
    output_results(0,method);
    unsigned int n_steps=0;
    while (time<final_time)
      {
        // The last time step is chosend such that the final time is exactly
        // reached.
        if (time+time_step>final_time)
          time_step = final_time-time;

        time = embedded_explicit_runge_kutta.evolve_one_time_step(
                 std_cxx11::bind(&Diffusion::evaluate_diffusion,this,std_cxx11::_1,std_cxx11::_2),
                 time,time_step,solution);

        if ((n_steps+1)%10==0)
          output_results(n_steps+1,method);

        time_step = embedded_explicit_runge_kutta.get_status().delta_t_guess;
        ++n_steps;
      }

    return n_steps;
  }



  // @sect5{<code>Diffusion::run</code>}
  void Diffusion::run()
  {
    // We create the grid (a [0,5]x[0,5] square) and refine the mesh four times.
    // The final grid has 16 by 16 cells, for a total of 256.
    GridGenerator::hyper_cube(triangulation, 0., 5.);
    triangulation.refine_global(4);

    // Set the boundary indicator for x=0 and x=5 to 1.
    typename Triangulation<2>::active_cell_iterator
    cell = triangulation.begin_active(),
    endc = triangulation.end();

    for (; cell!=endc; ++cell)
      for (unsigned int f=0; f<GeometryInfo<2>::faces_per_cell; ++f)
        if (cell->face(f)->at_boundary())
          {
            if ((cell->face(f)->center()[0]==0.) || (cell->face(f)->center()[0]==5.))
              cell->face(f)->set_boundary_indicator(1);
            else
              cell->face(f)->set_boundary_indicator(0);
          }

    setup_system();

    assemble_system();

    unsigned int n_steps = 0;
    const unsigned int n_time_steps = 200;
    const double initial_time = 0.;
    const double final_time = 10.;

    // Next, we solve the diffusion problem using different Runge-Kutta methods.
    std::cout << "Explicit methods:" << std::endl;
    explicit_method (TimeStepping::FORWARD_EULER,
                     n_time_steps,
                     initial_time,
                     final_time);
    std::cout << "Forward Euler:            error=" << solution.l2_norm() << std::endl;

    explicit_method (TimeStepping::RK_THIRD_ORDER,
                     n_time_steps,
                     initial_time,
                     final_time);
    std::cout << "Third order Runge-Kutta:  error=" << solution.l2_norm() << std::endl;

    explicit_method (TimeStepping::RK_CLASSIC_FOURTH_ORDER,
                     n_time_steps,
                     initial_time,
                     final_time);
    std::cout << "Fourth order Runge-Kutta: error=" << solution.l2_norm() << std::endl;
    std::cout << std::endl;


    std::cout << "Implicit methods:" << std::endl;
    implicit_method (TimeStepping::BACKWARD_EULER,
                     n_time_steps,
                     initial_time,
                     final_time);
    std::cout << "Backward Euler:           error=" << solution.l2_norm() << std::endl;

    implicit_method (TimeStepping::IMPLICIT_MIDPOINT,
                     n_time_steps,
                     initial_time,
                     final_time);
    std::cout << "Implicit Midpoint:        error=" << solution.l2_norm() << std::endl;

    implicit_method (TimeStepping::CRANK_NICOLSON,
                     n_time_steps,
                     initial_time,
                     final_time);
    std::cout << "Crank-Nicolson:           error=" << solution.l2_norm() << std::endl;

    implicit_method (TimeStepping::SDIRK_TWO_STAGES,
                     n_time_steps,
                     initial_time,
                     final_time);
    std::cout << "SDIRK:                    error=" << solution.l2_norm() << std::endl;
    std::cout << std::endl;


    std::cout << "Embedded explicit methods:" << std::endl;
    n_steps = embedded_explicit_method (TimeStepping::HEUN_EULER,
                                        n_time_steps,
                                        initial_time,
                                        final_time);
    std::cout << "Heun-Euler:               error=" << solution.l2_norm() << std::endl;
    std::cout << "                steps performed=" << n_steps << std::endl;

    n_steps = embedded_explicit_method (TimeStepping::BOGACKI_SHAMPINE,
                                        n_time_steps,
                                        initial_time,
                                        final_time);
    std::cout << "Bogacki-Shampine:         error=" << solution.l2_norm() << std::endl;
    std::cout << "                steps performed=" << n_steps << std::endl;

    n_steps = embedded_explicit_method (TimeStepping::DOPRI,
                                        n_time_steps,
                                        initial_time,
                                        final_time);
    std::cout << "Dopri:                    error=" << solution.l2_norm() << std::endl;
    std::cout << "                steps performed=" << n_steps << std::endl;

    n_steps = embedded_explicit_method (TimeStepping::FEHLBERG,
                                        n_time_steps,
                                        initial_time,
                                        final_time);
    std::cout << "Fehlberg:                 error=" << solution.l2_norm() << std::endl;
    std::cout << "                steps performed=" << n_steps << std::endl;

    n_steps = embedded_explicit_method (TimeStepping::CASH_KARP,
                                        n_time_steps,
                                        initial_time,
                                        final_time);
    std::cout << "Cash-Karp:                error=" << solution.l2_norm() << std::endl;
    std::cout << "                steps performed=" << n_steps << std::endl;
  }
}



// @sect3{The <code>main()</code> function}
//
// The following <code>main</code> function is similar to previous examples
// and need not be commented on.
int main ()
{
  try
    {
      Step52::Diffusion diffusion;
      diffusion.run();
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
}
