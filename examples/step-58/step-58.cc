/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2018 - 2020 by the deal.II authors
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
 * Author: Wolfgang Bangerth, Colorado State University
 *         Yong-Yong Cai, Beijing Computational Science Research Center
 */

// @sect3{Include files}
// The program starts with the usual include files, all of which you should
// have seen before by now:
#include <deal.II/base/logstream.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/block_sparse_matrix.h>
#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/sparse_direct.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/matrix_tools.h>

#include <fstream>
#include <iostream>


// Then the usual placing of all content of this program into a namespace and
// the importation of the deal.II namespace into the one we will work in:
namespace Step58
{
  using namespace dealii;

  // @sect3{The <code>NonlinearSchroedingerEquation</code> class}
  //
  // Then the main class. It looks very much like the corresponding
  // classes in step-4 or step-6, with the only exception that the
  // matrices and vectors and everything else related to the
  // linear system are now storing elements of type `std::complex<double>`
  // instead of just `double`.
  template <int dim>
  class NonlinearSchroedingerEquation
  {
  public:
    NonlinearSchroedingerEquation();
    void run();

  private:
    void setup_system();
    void assemble_matrices();
    void do_half_phase_step();
    void do_full_spatial_step();
    void output_results() const;


    Triangulation<dim> triangulation;
    FE_Q<dim>          fe;
    DoFHandler<dim>    dof_handler;

    AffineConstraints<std::complex<double>> constraints;

    SparsityPattern                    sparsity_pattern;
    SparseMatrix<std::complex<double>> system_matrix;
    SparseMatrix<std::complex<double>> rhs_matrix;

    Vector<std::complex<double>> solution;
    Vector<std::complex<double>> system_rhs;

    double       time;
    double       time_step;
    unsigned int timestep_number;

    double kappa;
  };



  // @sect3{Equation data}

  // Before we go on filling in the details of the main class, let us define
  // the equation data corresponding to the problem, i.e. initial values, as
  // well as a right hand side class. (We will reuse the initial conditions
  // also for the boundary values, which we simply keep constant.) We do so
  // using classes derived
  // from the Function class template that has been used many times before, so
  // the following should not look surprising. The only point of interest is
  // that we here have a complex-valued problem, so we have to provide the
  // second template argument of the Function class (which would otherwise
  // default to `double`). Furthermore, the return type of the `value()`
  // functions is then of course also complex.
  //
  // What precisely these functions return has been discussed at the end of
  // the Introduction section.
  template <int dim>
  class InitialValues : public Function<dim, std::complex<double>>
  {
  public:
    InitialValues()
      : Function<dim, std::complex<double>>(1)
    {}

    virtual std::complex<double>
    value(const Point<dim> &p, const unsigned int component = 0) const override;
  };



  template <int dim>
  std::complex<double>
  InitialValues<dim>::value(const Point<dim> & p,
                            const unsigned int component) const
  {
    static_assert(dim == 2, "This initial condition only works in 2d.");

    (void)component;
    Assert(component == 0, ExcIndexRange(component, 0, 1));

    const std::vector<Point<dim>> vortex_centers = {{0, -0.3},
                                                    {0, +0.3},
                                                    {+0.3, 0},
                                                    {-0.3, 0}};

    const double R = 0.1;
    const double alpha =
      1. / (std::pow(R, dim) * std::pow(numbers::PI, dim / 2.));

    double sum = 0;
    for (const auto &vortex_center : vortex_centers)
      {
        const Tensor<1, dim> distance = p - vortex_center;
        const double         r        = distance.norm();

        sum += alpha * std::exp(-(r * r) / (R * R));
      }

    return {std::sqrt(sum), 0.};
  }



  template <int dim>
  class Potential : public Function<dim>
  {
  public:
    Potential() = default;
    virtual double value(const Point<dim> & p,
                         const unsigned int component = 0) const override;
  };



  template <int dim>
  double Potential<dim>::value(const Point<dim> & p,
                               const unsigned int component) const
  {
    (void)component;
    Assert(component == 0, ExcIndexRange(component, 0, 1));

    return (Point<dim>().distance(p) > 0.7 ? 1000 : 0);
  }



  // @sect3{Implementation of the <code>NonlinearSchroedingerEquation</code> class}

  // We start by specifying the implementation of the constructor
  // of the class. There is nothing of surprise to see here except
  // perhaps that we choose quadratic ($Q_2$) Lagrange elements --
  // the solution is expected to be smooth, so we choose a higher
  // polynomial degree than the bare minimum.
  template <int dim>
  NonlinearSchroedingerEquation<dim>::NonlinearSchroedingerEquation()
    : fe(2)
    , dof_handler(triangulation)
    , time(0)
    , time_step(1. / 128)
    , timestep_number(0)
    , kappa(1)
  {}


  // @sect4{Setting up data structures and assembling matrices}

  // The next function is the one that sets up the mesh, DoFHandler, and
  // matrices and vectors at the beginning of the program, i.e. before the
  // first time step. The first few lines are pretty much standard if you've
  // read through the tutorial programs at least up to step-6:
  template <int dim>
  void NonlinearSchroedingerEquation<dim>::setup_system()
  {
    GridGenerator::hyper_cube(triangulation, -1, 1);
    triangulation.refine_global(6);

    std::cout << "Number of active cells: " << triangulation.n_active_cells()
              << std::endl;

    dof_handler.distribute_dofs(fe);

    std::cout << "Number of degrees of freedom: " << dof_handler.n_dofs()
              << std::endl
              << std::endl;

    DynamicSparsityPattern dsp(dof_handler.n_dofs(), dof_handler.n_dofs());
    DoFTools::make_sparsity_pattern(dof_handler, dsp);
    sparsity_pattern.copy_from(dsp);

    system_matrix.reinit(sparsity_pattern);
    rhs_matrix.reinit(sparsity_pattern);

    solution.reinit(dof_handler.n_dofs());
    system_rhs.reinit(dof_handler.n_dofs());

    constraints.close();
  }



  // Next, we assemble the relevant matrices. The way we have written
  // the Crank-Nicolson discretization of the spatial step of the Strang
  // splitting (i.e., the second of the three partial steps in each time
  // step), we were led to the linear system
  // $\left[ -iM  +  \frac 14 k_{n+1} A + \frac 12 k_{n+1} W \right]
  //   \Psi^{(n,2)}
  //  =
  //  \left[ -iM  -  \frac 14 k_{n+1} A - \frac 12 k_{n+1} W \right]
  //   \Psi^{(n,1)}$.
  // In other words, there are two matrices in play here -- one for the
  // left and one for the right hand side. We build these matrices
  // separately. (One could avoid building the right hand side matrix
  // and instead just form the *action* of the matrix on $\Psi^{(n,1)}$
  // in each time step. This may or may not be more efficient, but
  // efficiency is not foremost on our minds for this program.)
  template <int dim>
  void NonlinearSchroedingerEquation<dim>::assemble_matrices()
  {
    const QGauss<dim> quadrature_formula(fe.degree + 1);

    FEValues<dim> fe_values(fe,
                            quadrature_formula,
                            update_values | update_gradients |
                              update_quadrature_points | update_JxW_values);

    const unsigned int dofs_per_cell = fe.dofs_per_cell;
    const unsigned int n_q_points    = quadrature_formula.size();

    FullMatrix<std::complex<double>> cell_matrix_lhs(dofs_per_cell,
                                                     dofs_per_cell);
    FullMatrix<std::complex<double>> cell_matrix_rhs(dofs_per_cell,
                                                     dofs_per_cell);

    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
    std::vector<double>                  potential_values(n_q_points);
    const Potential<dim>                 potential;

    for (const auto &cell : dof_handler.active_cell_iterators())
      {
        cell_matrix_lhs = std::complex<double>(0.);
        cell_matrix_rhs = std::complex<double>(0.);

        fe_values.reinit(cell);

        potential.value_list(fe_values.get_quadrature_points(),
                             potential_values);

        for (unsigned int q_index = 0; q_index < n_q_points; ++q_index)
          {
            for (unsigned int k = 0; k < dofs_per_cell; ++k)
              {
                for (unsigned int l = 0; l < dofs_per_cell; ++l)
                  {
                    const std::complex<double> i = {0, 1};

                    cell_matrix_lhs(k, l) +=
                      (-i * fe_values.shape_value(k, q_index) *
                         fe_values.shape_value(l, q_index) +
                       time_step / 4 * fe_values.shape_grad(k, q_index) *
                         fe_values.shape_grad(l, q_index) +
                       time_step / 2 * potential_values[q_index] *
                         fe_values.shape_value(k, q_index) *
                         fe_values.shape_value(l, q_index)) *
                      fe_values.JxW(q_index);

                    cell_matrix_rhs(k, l) +=
                      (-i * fe_values.shape_value(k, q_index) *
                         fe_values.shape_value(l, q_index) -
                       time_step / 4 * fe_values.shape_grad(k, q_index) *
                         fe_values.shape_grad(l, q_index) -
                       time_step / 2 * potential_values[q_index] *
                         fe_values.shape_value(k, q_index) *
                         fe_values.shape_value(l, q_index)) *
                      fe_values.JxW(q_index);
                  }
              }
          }

        cell->get_dof_indices(local_dof_indices);
        constraints.distribute_local_to_global(cell_matrix_lhs,
                                               local_dof_indices,
                                               system_matrix);
        constraints.distribute_local_to_global(cell_matrix_rhs,
                                               local_dof_indices,
                                               rhs_matrix);
      }
  }


  // @sect4{Implementing the Strang splitting steps}

  // Having set up all data structures above, we are now in a position to
  // implement the partial steps that form the Strang splitting scheme. We
  // start with the half-step to advance the phase, and that is used as the
  // first and last part of each time step.
  //
  // To this end, recall that for the first half step, we needed to
  // compute
  // $\psi^{(n,1)} = e^{-i\kappa|\psi^{(n,0)}|^2 \tfrac
  //  12\Delta t} \; \psi^{(n,0)}$. Here, $\psi^{(n,0)}=\psi^{(n)}$ and
  //  $\psi^{(n,1)}$
  // are functions of space and correspond to the output of the previous
  // complete time step and the result of the first of the three part steps,
  // respectively. A corresponding solution must be computed for the third
  // of the part steps, i.e.
  // $\psi^{(n,3)} = e^{-i\kappa|\psi^{(n,2)}|^2 \tfrac
  //  12\Delta t} \; \psi^{(n,2)}$, where $\psi^{(n,3)}=\psi^{(n+1)}$ is
  // the result of the time step as a whole, and its input $\psi^{(n,2)}$ is
  // the result of the spatial step of the Strang splitting.
  //
  // An important realization is that while $\psi^{(n,0)}(\mathbf x)$ may be a
  // finite element function (i.e., is piecewise polynomial), this may not
  // necessarily be the case for the "rotated" function in which we have updated
  // the phase using the exponential factor (recall that the amplitude of that
  // function remains constant as part of that step). In other words, we could
  // *compute* $\psi^{(n,1)}(\mathbf x)$ at every point $\mathbf x\in\Omega$,
  // but we can't represent it on a mesh because it is not a piecewise
  // polynomial function. The best we can do in a discrete setting is to compute
  // a projection or interpolation. In other words, we can compute
  // $\psi_h^{(n,1)}(\mathbf x) = \Pi_h
  //     \left(e^{-i\kappa|\psi_h^{(n,0)}(\mathbf x)|^2 \tfrac 12\Delta t}
  //     \; \psi_h^{(n,0)}(\mathbf x) \right)$ where $\Pi_h$ is a projection or
  // interpolation operator. The situation is particularly simple if we
  // choose the interpolation: Then, all we need to compute is the value of
  // the right hand side *at the node points* and use these as nodal
  // values for the vector $\Psi^{(n,1)}$ of degrees of freedom. This is
  // easily done because evaluating the right hand side at node points
  // for a Lagrange finite element as used here requires us to only
  // look at a single (complex-valued) entry of the node vector. In other
  // words, what we need to do is to compute
  // $\Psi^{(n,1)}_j = e^{-i\kappa|\Psi^{(n,0)}_j|^2 \tfrac
  //  12\Delta t} \; \Psi^{(n,0)}_j$ where $j$ loops over all of the entries
  // of our solution vector. This is what the function below does -- in fact,
  // it doesn't even use separate vectors for $\Psi^{(n,0)}$ and $\Psi^{(n,1)}$,
  // but just updates the same vector as appropriate.
  template <int dim>
  void NonlinearSchroedingerEquation<dim>::do_half_phase_step()
  {
    for (auto &value : solution)
      {
        const std::complex<double> i         = {0, 1};
        const double               magnitude = std::abs(value);

        value = std::exp(-i * kappa * magnitude * magnitude * (time_step / 2)) *
                value;
      }
  }



  // The next step is to solve for the linear system in each time step, i.e.,
  // the second half step of the Strang splitting we use. Recall that it had the
  // form $C\Psi^{(n,2)} = R\Psi^{(n,1)}$ where $C$ and $R$ are the matrices we
  // assembled earlier.
  //
  // The way we solve this here is using a direct solver. We first form the
  // right hand side $r=R\Psi^{(n,1)}$ using the SparseMatrix::vmult() function
  // and put the result into the `system_rhs` variable. We then call
  // SparseDirectUMFPACK::solver() which takes as argument the matrix $C$
  // and the right hand side vector and returns the solution in the same
  // vector `system_rhs`. The final step is then to put the solution so computed
  // back into the `solution` variable.
  template <int dim>
  void NonlinearSchroedingerEquation<dim>::do_full_spatial_step()
  {
    rhs_matrix.vmult(system_rhs, solution);

    SparseDirectUMFPACK direct_solver;
    direct_solver.solve(system_matrix, system_rhs);

    solution = system_rhs;
  }



  // @sect4{Creating graphical output}

  // The last of the helper functions and classes we ought to discuss are the
  // ones that create graphical output. The result of running the half and full
  // steps for the local and spatial parts of the Strang splitting is that we
  // have updated the `solution` vector $\Psi^n$ to the correct value at the end
  // of each time step. Its entries contain complex numbers for the solution at
  // the nodes of the finite element mesh.
  //
  // Complex numbers are not easily visualized. We can output their real and
  // imaginary parts, i.e., the fields $\text{Re}(\psi_h^{(n)}(\mathbf x))$ and
  // $\text{Im}(\psi_h^{(n)}(\mathbf x))$, and that is exactly what the DataOut
  // class does when one attaches as complex-valued vector via
  // DataOut::add_data_vector() and then calls DataOut::build_patches(). That is
  // indeed what we do below.

  // But oftentimes we are not particularly interested in real and imaginary
  // parts of the solution vector, but instead in derived quantities such as the
  // magnitude $|\psi|$ and phase angle $\text{arg}(\psi)$ of the solution. In
  // the context of quantum systems such as here, the magnitude itself is not so
  // interesting, but instead it is the "amplitude", $|\psi|^2$ that is a
  // physical property: it corresponds to the probability density of finding a
  // particle in a particular place of state. The way to put computed quantities
  // into output files for visualization -- as used in numerous previous
  // tutorial programs -- is to use the facilities of the DataPostprocessor and
  // derived classes. Specifically, both the amplitude of a complex number and
  // its phase angles are scalar quantities, and so the DataPostprocessorScalar
  // class is the right tool to base what we want to do on.
  //
  // Consequently, what we do here is to implement two classes
  // `ComplexAmplitude` and `ComplexPhase` that compute for each point at which
  // DataOut decides to generate output, the amplitudes $|\psi_h|^2$ and phases
  // $\text{arg}(\psi_h)$ of the solution for visualization. There is a fair
  // amount of boiler-plate code below, with the only interesting parts of
  // the first of these two classes being how its `evaluate_vector_field()`
  // function computes the `computed_quantities` object.
  //
  // (There is also the rather awkward fact that the <a
  // href="https://en.cppreference.com/w/cpp/numeric/complex/norm">std::norm()</a>
  // function does not compute what one would naively imagine, namely $|\psi|$,
  // but returns $|\psi|^2$ instead. It's certainly quite confusing to have a
  // standard function mis-named in such a way...)
  namespace DataPostprocessors
  {
    template <int dim>
    class ComplexAmplitude : public DataPostprocessorScalar<dim>
    {
    public:
      ComplexAmplitude();

      virtual void evaluate_vector_field(
        const DataPostprocessorInputs::Vector<dim> &inputs,
        std::vector<Vector<double>> &computed_quantities) const override;
    };


    template <int dim>
    ComplexAmplitude<dim>::ComplexAmplitude()
      : DataPostprocessorScalar<dim>("Amplitude", update_values)
    {}


    template <int dim>
    void ComplexAmplitude<dim>::evaluate_vector_field(
      const DataPostprocessorInputs::Vector<dim> &inputs,
      std::vector<Vector<double>> &               computed_quantities) const
    {
      Assert(computed_quantities.size() == inputs.solution_values.size(),
             ExcDimensionMismatch(computed_quantities.size(),
                                  inputs.solution_values.size()));

      for (unsigned int q = 0; q < computed_quantities.size(); ++q)
        {
          Assert(computed_quantities[q].size() == 1,
                 ExcDimensionMismatch(computed_quantities[q].size(), 1));
          Assert(inputs.solution_values[q].size() == 2,
                 ExcDimensionMismatch(inputs.solution_values[q].size(), 2));

          const std::complex<double> psi(inputs.solution_values[q](0),
                                         inputs.solution_values[q](1));
          computed_quantities[q](0) = std::norm(psi);
        }
    }



    // The second of these postprocessor classes computes the phase angle
    // of the complex-valued solution at each point. In other words, if we
    // represent $\psi(\mathbf x,t)=r(\mathbf x,t) e^{i\varphi(\mathbf x,t)}$,
    // then this class computes $\varphi(\mathbf x,t)$. The function
    // <a
    // href="https://en.cppreference.com/w/cpp/numeric/complex/arg">std::arg</a>
    // does this for us, and returns the angle as a real number between $-\pi$
    // and $+\pi$.
    //
    // For reasons that we will explain in detail in the results section, we
    // do not actually output this value at each location where output is
    // generated. Rather, we take the maximum over all evaluation points of the
    // phase and then fill each evaluation point's output field with this
    // maximum -- in essence, we output the phase angle as a piecewise constant
    // field, where each cell has its own constant value. The reasons for this
    // will become clear once you read through the discussion further down
    // below.
    template <int dim>
    class ComplexPhase : public DataPostprocessorScalar<dim>
    {
    public:
      ComplexPhase();

      virtual void evaluate_vector_field(
        const DataPostprocessorInputs::Vector<dim> &inputs,
        std::vector<Vector<double>> &computed_quantities) const override;
    };


    template <int dim>
    ComplexPhase<dim>::ComplexPhase()
      : DataPostprocessorScalar<dim>("Phase", update_values)
    {}


    template <int dim>
    void ComplexPhase<dim>::evaluate_vector_field(
      const DataPostprocessorInputs::Vector<dim> &inputs,
      std::vector<Vector<double>> &               computed_quantities) const
    {
      Assert(computed_quantities.size() == inputs.solution_values.size(),
             ExcDimensionMismatch(computed_quantities.size(),
                                  inputs.solution_values.size()));

      double max_phase = -numbers::PI;
      for (unsigned int q = 0; q < computed_quantities.size(); ++q)
        {
          Assert(computed_quantities[q].size() == 1,
                 ExcDimensionMismatch(computed_quantities[q].size(), 1));
          Assert(inputs.solution_values[q].size() == 2,
                 ExcDimensionMismatch(inputs.solution_values[q].size(), 2));

          max_phase =
            std::max(max_phase,
                     std::arg(
                       std::complex<double>(inputs.solution_values[q](0),
                                            inputs.solution_values[q](1))));
        }

      for (auto &output : computed_quantities)
        output(0) = max_phase;
    }

  } // namespace DataPostprocessors


  // Having so implemented these post-processors, we create output as we always
  // do. As in many other time-dependent tutorial programs, we attach flags to
  // DataOut that indicate the number of the time step and the current
  // simulation time.
  template <int dim>
  void NonlinearSchroedingerEquation<dim>::output_results() const
  {
    const DataPostprocessors::ComplexAmplitude<dim> complex_magnitude;
    const DataPostprocessors::ComplexPhase<dim>     complex_phase;

    DataOut<dim> data_out;

    data_out.attach_dof_handler(dof_handler);
    data_out.add_data_vector(solution, "Psi");
    data_out.add_data_vector(solution, complex_magnitude);
    data_out.add_data_vector(solution, complex_phase);
    data_out.build_patches();

    data_out.set_flags(DataOutBase::VtkFlags(time, timestep_number));

    const std::string filename =
      "solution-" + Utilities::int_to_string(timestep_number, 3) + ".vtu";
    std::ofstream output(filename);
    data_out.write_vtu(output);
  }



  // @sect4{Running the simulation}

  // The remaining step is how we set up the overall logic for this program.
  // It's really relatively simple: Set up the data structures; interpolate the
  // initial conditions onto finite element space; then iterate over all time
  // steps, and on each time step perform the three parts of the Strang
  // splitting method. Every tenth time step, we generate graphical output.
  // That's it.
  template <int dim>
  void NonlinearSchroedingerEquation<dim>::run()
  {
    setup_system();
    assemble_matrices();

    time = 0;
    VectorTools::interpolate(dof_handler, InitialValues<dim>(), solution);
    output_results();

    const double end_time = 1;
    for (; time <= end_time; time += time_step)
      {
        ++timestep_number;

        std::cout << "Time step " << timestep_number << " at t=" << time
                  << std::endl;

        do_half_phase_step();
        do_full_spatial_step();
        do_half_phase_step();

        if (timestep_number % 1 == 0)
          output_results();
      }
  }
} // namespace Step58



// @sect4{The main() function}
//
// The rest is again boiler plate and exactly as in almost all of the previous
// tutorial programs:
int main()
{
  try
    {
      using namespace Step58;

      NonlinearSchroedingerEquation<2> nse;
      nse.run();
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
