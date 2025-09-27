/* ------------------------------------------------------------------------
 *
 * SPDX-License-Identifier: LGPL-2.1-or-later
 * Copyright (C) 2021 - 2024 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * Part of the source code is dual licensed under Apache-2.0 WITH
 * LLVM-exception OR LGPL-2.1-or-later. Detailed license information
 * governing the source code and code contributions can be found in
 * LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
 *
 * ------------------------------------------------------------------------
 *
 * Author: Tyler Anderson, Colorado State University, 2021
 */


// @sect3{Include files}

// The program starts with the usual include files, all of which you should have
// seen before by now:
#include <deal.II/base/convergence_table.h>
#include <deal.II/base/function.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/utilities.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/vector.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/data_out_stack.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/matrix_creator.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/solution_transfer.h>
#include <deal.II/numerics/vector_tools.h>

#include <fstream>
#include <iostream>

// Then the usual placing of all content of this program into a namespace and
// the importation of the deal.II namespace into the one we will work in. We
// also define an identifier to allow for the MMS code to be run when
// <code>MMS</code> is defined. Otherwise, the program solves the original
// problem:
namespace BlackScholesSolver
{
  using namespace dealii;

#define MMS

  // @sect3{Solution Class}

  // This section creates a class for the known solution when testing using the
  // MMS. Here we are using $v(\tau,S) = -\tau^2 -S^2 + 6$ for the solution. We
  // need to include the solution equation and the gradient for the H1 seminorm
  // calculation.
  template <int dim>
  class Solution : public Function<dim>
  {
  public:
    Solution(const double maturity_time);

    virtual double value(const Point<dim>  &p,
                         const unsigned int component = 0) const override;

    virtual Tensor<1, dim>
    gradient(const Point<dim>  &p,
             const unsigned int component = 0) const override;

  private:
    const double maturity_time;
  };


  template <int dim>
  Solution<dim>::Solution(const double maturity_time)
    : maturity_time(maturity_time)
  {
    Assert(dim == 1, ExcNotImplemented());
  }


  template <int dim>
  double Solution<dim>::value(const Point<dim>  &p,
                              const unsigned int component) const
  {
    return -Utilities::fixed_power<2, double>(p[component]) -
           Utilities::fixed_power<2, double>(this->get_time()) + 6;
  }


  template <int dim>
  Tensor<1, dim> Solution<dim>::gradient(const Point<dim>  &p,
                                         const unsigned int component) const
  {
    return Point<dim>(-2 * p[component]);
  }



  // @sect3{Equation Data}

  // In the following classes and functions, we implement the right hand side
  // and boundary values that define this problem and for which we need function
  // objects. The right hand side is chosen as discussed at the end of the
  // introduction.
  //
  // First, we handle the initial condition.
  template <int dim>
  class InitialConditions : public Function<dim>
  {
  public:
    InitialConditions(const double strike_price);

    virtual double value(const Point<dim>  &p,
                         const unsigned int component = 0) const override;

  private:
    const double strike_price;
  };


  template <int dim>
  InitialConditions<dim>::InitialConditions(const double strike_price)
    : strike_price(strike_price)
  {}


  template <int dim>
  double InitialConditions<dim>::value(const Point<dim>  &p,
                                       const unsigned int component) const
  {
#ifdef MMS
    return -Utilities::fixed_power<2, double>(p[component]) + 6;
#else
    return std::max(p[component] - strike_price, 0.);
#endif
  }



  // Next, we handle the left boundary condition.
  template <int dim>
  class LeftBoundaryValues : public Function<dim>
  {
  public:
    virtual double value(const Point<dim>  &p,
                         const unsigned int component = 0) const override;
  };


  template <int dim>
  double LeftBoundaryValues<dim>::value(const Point<dim> &,
                                        const unsigned int /*component*/) const
  {
#ifdef MMS
    return -Utilities::fixed_power<2, double>(this->get_time()) + 6;
#else
    return 0.;
#endif
  }



  // Then, we handle the right boundary condition.
  template <int dim>
  class RightBoundaryValues : public Function<dim>
  {
  public:
    RightBoundaryValues(const double strike_price, const double interest_rate);

    virtual double value(const Point<dim>  &p,
                         const unsigned int component = 0) const override;

  private:
    const double strike_price;
    const double interest_rate;
  };


  template <int dim>
  RightBoundaryValues<dim>::RightBoundaryValues(const double strike_price,
                                                const double interest_rate)
    : strike_price(strike_price)
    , interest_rate(interest_rate)
  {}


  template <int dim>
  double RightBoundaryValues<dim>::value(const Point<dim>  &p,
                                         const unsigned int component) const
  {
#ifdef MMS
    return -Utilities::fixed_power<2, double>(p[component]) -
           Utilities::fixed_power<2, double>(this->get_time()) + 6;
#else
    return (p[component] - strike_price) *
           exp((-interest_rate) * (this->get_time()));
#endif
  }



  // Finally, we handle the right hand side.
  template <int dim>
  class RightHandSide : public Function<dim>
  {
  public:
    RightHandSide(const double asset_volatility, const double interest_rate);

    virtual double value(const Point<dim>  &p,
                         const unsigned int component = 0) const override;

  private:
    const double asset_volatility;
    const double interest_rate;
  };


  template <int dim>
  RightHandSide<dim>::RightHandSide(const double asset_volatility,
                                    const double interest_rate)
    : asset_volatility(asset_volatility)
    , interest_rate(interest_rate)
  {}


  template <int dim>
  double RightHandSide<dim>::value(const Point<dim>  &p,
                                   const unsigned int component) const
  {
#ifdef MMS
    return 2 * (this->get_time()) -
           Utilities::fixed_power<2, double>(asset_volatility * p[component]) -
           2 * interest_rate * Utilities::fixed_power<2, double>(p[component]) -
           interest_rate *
             (-Utilities::fixed_power<2, double>(p[component]) -
              Utilities::fixed_power<2, double>(this->get_time()) + 6);
#else
    (void)p;
    (void)component;
    return 0.0;
#endif
  }



  // @sect3{The <code>BlackScholes</code> Class}

  // The next piece is the declaration of the main class of this program. This
  // is very similar to the Step-26 tutorial, with some modifications. New
  // matrices had to be added to calculate the A and B matrices, as well as the
  // $V_{diff}$ vector mentioned in the introduction. We also define the
  // parameters used in the problem.
  //
  // - <code>maximum_stock_price</code>: The imposed upper bound on the spatial
  // domain. This is the maximum allowed stock price.
  // - <code>maturity_time</code>: The upper bound on the time domain. This is
  // when the option expires.\n
  // - <code>asset_volatility</code>: The volatility of the stock price.\n
  // - <code>interest_rate</code>: The risk free interest rate.\n
  // - <code>strike_price</code>: The agreed upon price that the buyer will
  // have the option of purchasing  the stocks at the expiration time.
  //
  // Some slight differences between this program and step-26 are the creation
  // of the <code>a_matrix</code> and the <code>b_matrix</code>, which is
  // described in the introduction. We then also need to store the current time,
  // the size of the time step, and the number of the current time step.
  // Next, we will store the output into a <code>DataOutStack</code>
  // variable because we will be layering the solution at each time on top of
  // one another to create the solution manifold. Then, we have a variable that
  // stores the current cycle and number of cycles that we will run when
  // calculating the solution. The cycle is one full solution calculation given
  // a mesh. We refine the mesh once in between each cycle to exhibit the
  // convergence properties of our program. Finally, we store the convergence
  // data into a convergence table.
  //
  // As far as member functions are concerned, we have a function that
  // calculates the convergence information for each cycle, called
  // <code>process_solution</code>. This is just like what is done in step-7.
  template <int dim>
  class BlackScholes
  {
  public:
    BlackScholes();

    void run();

  private:
    void setup_system();
    void solve_time_step();
    void refine_grid();
    void process_solution();
    void add_results_for_output();
    void write_convergence_table();

    const double maximum_stock_price;
    const double maturity_time;
    const double asset_volatility;
    const double interest_rate;
    const double strike_price;

    Triangulation<dim> triangulation;
    const FE_Q<dim>    fe;
    DoFHandler<dim>    dof_handler;

    AffineConstraints<double> constraints;

    SparsityPattern      sparsity_pattern;
    SparseMatrix<double> mass_matrix;
    SparseMatrix<double> laplace_matrix;
    SparseMatrix<double> a_matrix;
    SparseMatrix<double> b_matrix;
    SparseMatrix<double> system_matrix;

    Vector<double> solution;
    Vector<double> system_rhs;

    double time;
    double time_step;

    const double       theta;
    const unsigned int n_cycles;
    const unsigned int n_time_steps;

    DataOutStack<dim>        data_out_stack;
    std::vector<std::string> solution_names;

    ConvergenceTable convergence_table;
  };

  // @sect3{The <code>BlackScholes</code> Implementation}

  // Now, we get to the implementation of the main class. We will set the values
  // for the various parameters used in the problem. These were chosen because
  // they are fairly normal values for these parameters. Although the stock
  // price has no upper bound in reality (it is in fact infinite), we impose
  // an upper bound that is twice the strike price. This is a somewhat arbitrary
  // choice to be twice the strike price, but it is large enough to see the
  // interesting parts of the solution.
  template <int dim>
  BlackScholes<dim>::BlackScholes()
    : maximum_stock_price(1.)
    , maturity_time(1.)
    , asset_volatility(.2)
    , interest_rate(0.05)
    , strike_price(0.5)
    , fe(1)
    , dof_handler(triangulation)
    , time(0.0)
    , theta(0.5)
    , n_cycles(4)
    , n_time_steps(5000)
  {
    Assert(dim == 1, ExcNotImplemented());
  }

  // @sect4{<code>BlackScholes::setup_system</code>}

  // The next function sets up the DoFHandler object, computes
  // the constraints, and sets the linear algebra objects to their correct
  // sizes. We also compute the @ref GlossMassMatrix "mass matrix" here by calling a function from the
  // library. We will compute the other 3 matrices next, because these need to
  // be computed 'by hand'.
  //
  // Note, that the time step is initialized here because the maturity time was
  // needed to compute the time step.
  template <int dim>
  void BlackScholes<dim>::setup_system()
  {
    dof_handler.distribute_dofs(fe);

    time_step = maturity_time / n_time_steps;

    constraints.clear();
    DoFTools::make_hanging_node_constraints(dof_handler, constraints);
    constraints.close();
    DynamicSparsityPattern dsp(dof_handler.n_dofs());
    DoFTools::make_sparsity_pattern(dof_handler,
                                    dsp,
                                    constraints,
                                    /*keep_constrained_dofs = */ true);
    sparsity_pattern.copy_from(dsp);

    mass_matrix.reinit(sparsity_pattern);
    laplace_matrix.reinit(sparsity_pattern);
    a_matrix.reinit(sparsity_pattern);
    b_matrix.reinit(sparsity_pattern);
    system_matrix.reinit(sparsity_pattern);

    MatrixCreator::create_mass_matrix(dof_handler,
                                      QGauss<dim>(fe.degree + 1),
                                      mass_matrix);

    // Below is the code to create the Laplace matrix with non-constant
    // coefficients. This corresponds to the matrix D in the introduction. This
    // non-constant coefficient is represented in the
    // <code>current_coefficient</code> variable.
    const unsigned int dofs_per_cell = fe.dofs_per_cell;
    FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
    const QGauss<dim>  quadrature_formula(fe.degree + 1);
    FEValues<dim>      fe_values(fe,
                            quadrature_formula,
                            update_values | update_gradients |
                              update_quadrature_points | update_JxW_values);
    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
    for (const auto &cell : dof_handler.active_cell_iterators())
      {
        cell_matrix = 0.;
        fe_values.reinit(cell);
        for (const unsigned int q_index : fe_values.quadrature_point_indices())
          {
            const double current_coefficient =
              fe_values.quadrature_point(q_index).square();
            for (const unsigned int i : fe_values.dof_indices())
              {
                for (const unsigned int j : fe_values.dof_indices())
                  cell_matrix(i, j) +=
                    (current_coefficient *              // (x_q)^2
                     fe_values.shape_grad(i, q_index) * // grad phi_i(x_q)
                     fe_values.shape_grad(j, q_index) * // grad phi_j(x_q)
                     fe_values.JxW(q_index));           // dx
              }
          }
        cell->get_dof_indices(local_dof_indices);
        for (const unsigned int i : fe_values.dof_indices())
          {
            for (const unsigned int j : fe_values.dof_indices())
              laplace_matrix.add(local_dof_indices[i],
                                 local_dof_indices[j],
                                 cell_matrix(i, j));
          }
      }

    // Now we will create the A matrix. Below is the code to create the matrix A
    // as discussed in the introduction. The non constant coefficient is again
    // represented in  the <code>current_coefficient</code> variable.
    for (const auto &cell : dof_handler.active_cell_iterators())
      {
        cell_matrix = 0.;
        fe_values.reinit(cell);
        for (const unsigned int q_index : fe_values.quadrature_point_indices())
          {
            const Tensor<1, dim> current_coefficient =
              fe_values.quadrature_point(q_index);
            for (const unsigned int i : fe_values.dof_indices())
              {
                for (const unsigned int j : fe_values.dof_indices())
                  {
                    cell_matrix(i, j) +=
                      (current_coefficient *               // x_q
                       fe_values.shape_grad(i, q_index) *  // grad phi_i(x_q)
                       fe_values.shape_value(j, q_index) * // phi_j(x_q)
                       fe_values.JxW(q_index));            // dx
                  }
              }
          }
        cell->get_dof_indices(local_dof_indices);
        for (const unsigned int i : fe_values.dof_indices())
          {
            for (const unsigned int j : fe_values.dof_indices())
              a_matrix.add(local_dof_indices[i],
                           local_dof_indices[j],
                           cell_matrix(i, j));
          }
      }

    // Finally we will create the matrix B. Below is the code to create the
    // matrix B as discussed in the introduction. The non constant coefficient
    // is again represented in the <code>current_coefficient</code> variable.
    for (const auto &cell : dof_handler.active_cell_iterators())
      {
        cell_matrix = 0.;
        fe_values.reinit(cell);
        for (const unsigned int q_index : fe_values.quadrature_point_indices())
          {
            const Tensor<1, dim> current_coefficient =
              fe_values.quadrature_point(q_index);
            for (const unsigned int i : fe_values.dof_indices())
              {
                for (const unsigned int j : fe_values.dof_indices())
                  cell_matrix(i, j) +=
                    (current_coefficient *               // x_q
                     fe_values.shape_value(i, q_index) * // phi_i(x_q)
                     fe_values.shape_grad(j, q_index) *  // grad phi_j(x_q)
                     fe_values.JxW(q_index));            // dx
              }
          }
        cell->get_dof_indices(local_dof_indices);
        for (const unsigned int i : fe_values.dof_indices())
          {
            for (const unsigned int j : fe_values.dof_indices())
              b_matrix.add(local_dof_indices[i],
                           local_dof_indices[j],
                           cell_matrix(i, j));
          }
      }

    solution.reinit(dof_handler.n_dofs());
    system_rhs.reinit(dof_handler.n_dofs());
  }

  // @sect4{<code>BlackScholes::solve_time_step</code>}

  // The next function is the one that solves the actual linear system for a
  // single time step. The only interesting thing here is that the matrices
  // we have built are symmetric positive definite, so we can use the
  // conjugate gradient method.
  template <int dim>
  void BlackScholes<dim>::solve_time_step()
  {
    SolverControl                          solver_control(1000, 1e-12);
    SolverCG<Vector<double>>               cg(solver_control);
    PreconditionSSOR<SparseMatrix<double>> preconditioner;
    preconditioner.initialize(system_matrix, 1.0);
    cg.solve(system_matrix, solution, system_rhs, preconditioner);
    constraints.distribute(solution);
  }

  // @sect4{<code>BlackScholes::add_results_for_output</code>}

  // This is simply the function to stitch the solution pieces together. For
  // this, we create a new layer at each time, and then add the solution vector
  // for that timestep. The function then stitches this together with the old
  // solutions using 'build_patches'.
  template <int dim>
  void BlackScholes<dim>::add_results_for_output()
  {
    data_out_stack.new_parameter_value(time, time_step);
    data_out_stack.attach_dof_handler(dof_handler);
    data_out_stack.add_data_vector(solution, solution_names);
    data_out_stack.build_patches(2);
    data_out_stack.finish_parameter_value();
  }

  // @sect4{<code>BlackScholes::refine_grid</code>}

  // It is somewhat unnecessary to have a function for the global refinement
  // that we do. The reason for the function is to allow for the possibility of
  // an adaptive refinement later.
  template <int dim>
  void BlackScholes<dim>::refine_grid()
  {
    triangulation.refine_global(1);
  }

  // @sect4{<code>BlackScholes::process_solution</code>}

  // This is where we calculate the convergence and error data to evaluate the
  // effectiveness of the program. Here, we calculate the $L^2$, $H^1$ and
  // $L^{\infty}$ norms.
  template <int dim>
  void BlackScholes<dim>::process_solution()
  {
    Solution<dim> sol(maturity_time);
    sol.set_time(time);
    Vector<float> difference_per_cell(triangulation.n_active_cells());
    VectorTools::integrate_difference(dof_handler,
                                      solution,
                                      sol,
                                      difference_per_cell,
                                      QGauss<dim>(fe.degree + 1),
                                      VectorTools::L2_norm);
    const double L2_error =
      VectorTools::compute_global_error(triangulation,
                                        difference_per_cell,
                                        VectorTools::L2_norm);
    VectorTools::integrate_difference(dof_handler,
                                      solution,
                                      sol,
                                      difference_per_cell,
                                      QGauss<dim>(fe.degree + 1),
                                      VectorTools::H1_seminorm);
    const double H1_error =
      VectorTools::compute_global_error(triangulation,
                                        difference_per_cell,
                                        VectorTools::H1_seminorm);
    const QTrapezoid<1>  q_trapezoid;
    const QIterated<dim> q_iterated(q_trapezoid, fe.degree * 2 + 1);
    VectorTools::integrate_difference(dof_handler,
                                      solution,
                                      sol,
                                      difference_per_cell,
                                      q_iterated,
                                      VectorTools::Linfty_norm);
    const double Linfty_error =
      VectorTools::compute_global_error(triangulation,
                                        difference_per_cell,
                                        VectorTools::Linfty_norm);
    const unsigned int n_active_cells = triangulation.n_active_cells();
    const unsigned int n_dofs         = dof_handler.n_dofs();
    convergence_table.add_value("cells", n_active_cells);
    convergence_table.add_value("dofs", n_dofs);
    convergence_table.add_value("L2", L2_error);
    convergence_table.add_value("H1", H1_error);
    convergence_table.add_value("Linfty", Linfty_error);
  }

  //@sect4{<code>BlackScholes::write_convergence_table</code> }

  // This next part is building the convergence and error tables. By this, we
  // need to set the settings for how to output the data that was calculated
  // during <code>BlackScholes::process_solution</code>. First, we will create
  // the headings and set up the cells properly. During this, we will also
  // prescribe the precision of our results. Then we will write the calculated
  // errors based on the $L^2$, $H^1$, and $L^{\infty}$ norms to the console and
  // to the error LaTeX file.
  template <int dim>
  void BlackScholes<dim>::write_convergence_table()
  {
    convergence_table.set_precision("L2", 3);
    convergence_table.set_precision("H1", 3);
    convergence_table.set_precision("Linfty", 3);
    convergence_table.set_scientific("L2", true);
    convergence_table.set_scientific("H1", true);
    convergence_table.set_scientific("Linfty", true);
    convergence_table.set_tex_caption("cells", "\\# cells");
    convergence_table.set_tex_caption("dofs", "\\# dofs");
    convergence_table.set_tex_caption("L2", "@f$L^2@f$-error");
    convergence_table.set_tex_caption("H1", "@f$H^1@f$-error");
    convergence_table.set_tex_caption("Linfty", "@f$L^\\infty@f$-error");
    convergence_table.set_tex_format("cells", "r");
    convergence_table.set_tex_format("dofs", "r");
    std::cout << std::endl;
    convergence_table.write_text(std::cout);
    std::string error_filename = "error";
    error_filename += "-global";
    error_filename += ".tex";
    std::ofstream error_table_file(error_filename);
    convergence_table.write_tex(error_table_file);

    // Next, we will make the convergence table. We will again write this to
    // the console and to the convergence LaTeX file.
    convergence_table.add_column_to_supercolumn("cells", "n cells");
    std::vector<std::string> new_order;
    new_order.emplace_back("n cells");
    new_order.emplace_back("H1");
    new_order.emplace_back("L2");
    convergence_table.set_column_order(new_order);
    convergence_table.evaluate_convergence_rates(
      "L2", ConvergenceTable::reduction_rate);
    convergence_table.evaluate_convergence_rates(
      "L2", ConvergenceTable::reduction_rate_log2);
    convergence_table.evaluate_convergence_rates(
      "H1", ConvergenceTable::reduction_rate);
    convergence_table.evaluate_convergence_rates(
      "H1", ConvergenceTable::reduction_rate_log2);
    std::cout << std::endl;
    convergence_table.write_text(std::cout);
    std::string conv_filename = "convergence";
    conv_filename += "-global";
    switch (fe.degree)
      {
        case 1:
          conv_filename += "-q1";
          break;
        case 2:
          conv_filename += "-q2";
          break;
        default:
          DEAL_II_NOT_IMPLEMENTED();
      }
    conv_filename += ".tex";
    std::ofstream table_file(conv_filename);
    convergence_table.write_tex(table_file);
  }

  // @sect4{<code>BlackScholes::run</code>}

  // Now we get to the main driver of the program. This is where we do all the
  // work of looping through the time steps and calculating the solution vector
  // each time. Here at the top, we set the initial refinement value and then
  // create a mesh. Then we refine this mesh once. Next, we set up the
  // data_out_stack object to store our solution. Finally, we start a for loop
  // to loop through the cycles. This lets us recalculate a solution for each
  // successive mesh refinement. At the beginning of each iteration, we need to
  // reset the time and time step number. We introduce an if statement to
  // accomplish this because we don't want to do this on the first iteration.
  template <int dim>
  void BlackScholes<dim>::run()
  {
    GridGenerator::hyper_cube(triangulation, 0.0, maximum_stock_price, true);
    triangulation.refine_global(0);

    solution_names.emplace_back("u");
    data_out_stack.declare_data_vector(solution_names,
                                       DataOutStack<dim>::dof_vector);

    Vector<double> vmult_result;
    Vector<double> forcing_terms;

    for (unsigned int cycle = 0; cycle < n_cycles; ++cycle)
      {
        if (cycle != 0)
          {
            refine_grid();
            time = 0.0;
          }

        setup_system();

        std::cout << std::endl
                  << "===========================================" << std::endl
                  << "Cycle " << cycle << ':' << std::endl
                  << "Number of active cells: "
                  << triangulation.n_active_cells() << std::endl
                  << "Number of degrees of freedom: " << dof_handler.n_dofs()
                  << std::endl
                  << std::endl;

        VectorTools::interpolate(dof_handler,
                                 InitialConditions<dim>(strike_price),
                                 solution);

        if (cycle == (n_cycles - 1))
          {
            add_results_for_output();
          }

        // Next, we run the main loop which runs until we exceed the maturity
        // time. We first compute the right hand side of the equation, which is
        // described in the introduction. Recall that it contains the term
        // $\left[-\frac{1}{4}k_n\sigma^2\mathbf{D}-k_nr\mathbf{M}+k_n\sigma^2
        // \mathbf{B}-k_nr\mathbf{A}+\mathbf{M}\right]V^{n-1}$. We put these
        // terms into the variable system_rhs, with the help of a temporary
        // vector:
        vmult_result.reinit(dof_handler.n_dofs());
        forcing_terms.reinit(dof_handler.n_dofs());
        for (unsigned int timestep_number = 0; timestep_number < n_time_steps;
             ++timestep_number)
          {
            time += time_step;

            if (timestep_number % 1000 == 0)
              std::cout << "Time step " << timestep_number << " at t=" << time
                        << std::endl;

            mass_matrix.vmult(system_rhs, solution);

            laplace_matrix.vmult(vmult_result, solution);
            system_rhs.add(
              (-1) * (1 - theta) * time_step *
                Utilities::fixed_power<2, double>(asset_volatility) * 0.5,
              vmult_result);
            mass_matrix.vmult(vmult_result, solution);

            system_rhs.add((-1) * (1 - theta) * time_step * interest_rate * 2,
                           vmult_result);

            a_matrix.vmult(vmult_result, solution);
            system_rhs.add((-1) * time_step * interest_rate, vmult_result);

            b_matrix.vmult(vmult_result, solution);
            system_rhs.add(
              (-1) * Utilities::fixed_power<2, double>(asset_volatility) *
                time_step * 1,
              vmult_result);

            // The second piece is to compute the contributions of the source
            // terms. This corresponds to the term $-k_n\left[\frac{1}{2}F^{n-1}
            // +\frac{1}{2}F^n\right]$. The following code calls
            // VectorTools::create_right_hand_side to compute the vectors $F$,
            // where we set the time of the right hand side (source) function
            // before we evaluate it. The result of this all ends up in the
            // forcing_terms variable:
            RightHandSide<dim> rhs_function(asset_volatility, interest_rate);
            rhs_function.set_time(time);
            VectorTools::create_right_hand_side(dof_handler,
                                                QGauss<dim>(fe.degree + 1),
                                                rhs_function,
                                                forcing_terms);
            forcing_terms *= time_step * theta;
            system_rhs -= forcing_terms;

            rhs_function.set_time(time - time_step);
            VectorTools::create_right_hand_side(dof_handler,
                                                QGauss<dim>(fe.degree + 1),
                                                rhs_function,
                                                forcing_terms);
            forcing_terms *= time_step * (1 - theta);
            system_rhs -= forcing_terms;

            // Next, we add the forcing terms to the ones that come from the
            // time stepping, and also build the matrix $\left[\mathbf{M}+
            // \frac{1}{4}k_n\sigma^2\mathbf{D}+k_nr\mathbf{M}\right]$ that we
            // have to invert in each time step. The final piece of these
            // operations is to eliminate hanging node constrained degrees of
            // freedom from the linear system:
            system_matrix.copy_from(mass_matrix);
            system_matrix.add(
              (theta)*time_step *
                Utilities::fixed_power<2, double>(asset_volatility) * 0.5,
              laplace_matrix);
            system_matrix.add((time_step)*interest_rate * theta * (1 + 1),
                              mass_matrix);

            constraints.condense(system_matrix, system_rhs);

            // There is one more operation we need to do before we can solve it:
            // boundary values. To this end, we create a boundary value object,
            // set the proper time to the one of the current time step, and
            // evaluate it as we have done many times before. The result is used
            // to also set the correct boundary values in the linear system:
            {
              RightBoundaryValues<dim> right_boundary_function(strike_price,
                                                               interest_rate);
              LeftBoundaryValues<dim>  left_boundary_function;
              right_boundary_function.set_time(time);
              left_boundary_function.set_time(time);
              std::map<types::global_dof_index, double> boundary_values;
              VectorTools::interpolate_boundary_values(dof_handler,
                                                       0,
                                                       left_boundary_function,
                                                       boundary_values);
              VectorTools::interpolate_boundary_values(dof_handler,
                                                       1,
                                                       right_boundary_function,
                                                       boundary_values);
              MatrixTools::apply_boundary_values(boundary_values,
                                                 system_matrix,
                                                 solution,
                                                 system_rhs);
            }

            // With this out of the way, all we have to do is solve the system,
            // generate graphical data on the last cycle, and create the
            // convergence table data.
            solve_time_step();

            if (cycle == (n_cycles - 1))
              {
                add_results_for_output();
              }
          }
#ifdef MMS
        process_solution();
#endif
      }

    const std::string filename = "solution.vtk";
    std::ofstream     output(filename);
    data_out_stack.write_vtk(output);

#ifdef MMS
    write_convergence_table();
#endif
  }

} // namespace BlackScholesSolver

// @sect3{The <code>main</code> Function}

// Having made it this far, there is, again, nothing much to discuss for the
// main function of this program: it looks like all such functions since step-6.
int main()
{
  try
    {
      using namespace BlackScholesSolver;

      BlackScholes<1> black_scholes_solver;
      black_scholes_solver.run();
    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl
                << std::endl
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
      std::cerr << std::endl
                << std::endl
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
