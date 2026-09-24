/* ------------------------------------------------------------------------
 *
 * SPDX-License-Identifier: LGPL-2.1-or-later
 * Copyright (C) 2010 - 2024 by the deal.II authors
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
 * Author: Arthur Bawin, Polytechnique Montreal, 2026.
 */

#include <deal.II/lac/generic_linear_algebra.h>

// This program only uses PETSc parallel matrix and vector wrappers, although
// it is quite straightforward to modify it to use Trilinos wrappers.
namespace LA
{
#if defined(DEAL_II_WITH_PETSC) && !defined(DEAL_II_PETSC_WITH_COMPLEX)
  using namespace dealii::LinearAlgebraPETSc;
#else
#  error DEAL_II_WITH_PETSC required
#endif
} // namespace LA

// The include files for this program should all be familiar.
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/timer.h>
#include <deal.II/base/utilities.h>
#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/index_set.h>
#include <deal.II/base/work_stream.h>
#include <deal.II/base/discrete_time.h>
#include <deal.II/base/table_handler.h>

#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/sparsity_tools.h>
#include <deal.II/lac/petsc_matrix_base.h>
#include <deal.II/lac/petsc_vector_base.h>
#include <deal.II/lac/petsc_solver.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/grid/grid_tools.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/vector_tools_boundary.h>
#include <deal.II/numerics/data_out.h>

#include <fstream>
#include <iostream>

namespace Step106
{
  using namespace dealii;

  // @sect3{Inlet velocity}
  //
  // This small class template defines the velocity field at the inlet.
  // To be a bit modular, this function takes as arguments the total number of
  // components in the FESystem, and the starting component for the velocity
  // field, @u_lower. It only assigns value to the velocity components.
  template <int dim>
  class Inlet : public Function<dim>
  {
  public:
    Inlet(const unsigned int u_lower, const unsigned int n_components)
      : Function<dim>(n_components)
      , u_lower(u_lower)
    {}

    virtual void vector_value(const Point<dim> & /*p*/,
                              Vector<double> &values) const override
    {
      for (unsigned int d = 0; d < dim; ++d)
        // Small perturbation to trigger the Von Karman streets earlier
        values[u_lower + d] = 0;
      values[u_lower] = 1.;
    }

  public:
    // Lower bound (first component) of the velocity variable
    const unsigned int u_lower;
  };

  // @sect3{Velocity initial condition}
  //
  // This class is not strictly required, as we could simply reuse the inlet
  // boundary condition to set the initial velocity, but this allows us to
  // cheat a little beat and give an initial vertical velocity. This creates
  // a slight asymmetry in the flow, and helps trigger the Von Karman streets
  // earlier.
  template <int dim>
  class InitialVelocity : public Function<dim>
  {
  public:
    InitialVelocity(const unsigned int u_lower, const unsigned int n_components)
      : Function<dim>(n_components)
      , u_lower(u_lower)
    {}

    virtual void vector_value(const Point<dim> & /*p*/,
                              Vector<double> &values) const override
    {
      for (unsigned int d = 0; d < dim; ++d)
        // Perturbation
        values[u_lower + d] = 0.01;
      values[u_lower] = 1.;
    }

  public:
    const unsigned int u_lower;
  };

  // @sect3{Time handler for BDF1 time stepping method}
  //
  // We'll be using a first order backward differentiation formula method
  // (BDF1), also simply known as the implicit Euler method. The class
  // below extends the <code>DiscreteTime<code> class, which handles the
  // beginning, ending, and step size information, by also storing BDF
  // methods-related data, i.e., the number of previous solutions to store, and
  // the coefficients of the BDF expansion. We will use this class to compute
  // the approximation of the time derivative of the velocity field.
  //
  // Currently, only the first order method is implemented: it is easy to extend
  // this class to handle order $n$ BDF methods by adding the associated
  // constant or variable time step expansion coefficients. However, since these
  // methods require the previous $n$ solutions to approximate the time
  // derivative, one then has to deal with the problem of properly initializing
  // the previous solutions to start the simulation.
  class BDF1DiscreteTime : public DiscreteTime
  {
  public:
    BDF1DiscreteTime(const unsigned int bdf_order,
                     const double       start_time,
                     const double       end_time,
                     const double       desired_start_step_size = 0.)
      : DiscreteTime(start_time, end_time, desired_start_step_size)
    {
      const double dt = this->get_next_step_size();

      if (bdf_order == 1)
        {
          n_previous_sol   = 1;
          bdf_coefficients = {1. / dt, -1. / dt};
        }
      else
        DEAL_II_NOT_IMPLEMENTED();
    }

    // Return the coefficients of the BDF expansion.
    const std::vector<double> &get_bdf_coefficients() const
    {
      return bdf_coefficients;
    }

    // Rotates the current and previous solutions, that is, assign $u^{n}$ to
    // $u^{n-1}$, and so on.
    void
    rotate_solutions(const LA::MPI::Vector        &present_solution,
                     std::vector<LA::MPI::Vector> &previous_solutions) const
    {
      for (unsigned int j = previous_solutions.size() - 1; j >= 1; --j)
        previous_solutions[j] = previous_solutions[j - 1];
      previous_solutions[0] = present_solution;
    }

    // Compute the time derivative of a quantity of type ValueType at the
    // <code>quadrature_node_index<code>-th quadrature node. In this program,
    // this is only called to evaluate dudt, in which case ValueType is a
    // Tensor<1, dim>.
    template <typename ValueType>
    ValueType compute_time_derivative(
      const unsigned int                         quadrature_node_index,
      const std::vector<ValueType>              &present_solution,
      const std::vector<std::vector<ValueType>> &previous_solutions) const
    {
      ValueType value_dot =
        bdf_coefficients[0] * present_solution[quadrature_node_index];
      for (unsigned int i = 1; i < bdf_coefficients.size(); ++i)
        value_dot += bdf_coefficients[i] *
                     previous_solutions[i - 1][quadrature_node_index];
      return value_dot;
    }

  public:
    unsigned int        n_previous_sol;
    std::vector<double> bdf_coefficients;
  };

  // @sect3{Assembly scratch data}
  //
  // This program, as many others, uses the WorkStream facilities to compute and
  // assemble the local matrix and right-hand side contributions to the global
  // linear system. See for example step-9 for an introduction to these
  // concepts.
  //
  // Although this machinery allows for multithreaded assembly, it is currently
  // not used in this program, as PETSc vectors are unfortunately not
  // thread-safe
  template <int dim>
  class ScratchData
  {
  public:
    ScratchData(const FESystem<dim>       &fe,
                const Mapping<dim>        &mapping,
                const Quadrature<dim>     &cell_quadrature,
                const Quadrature<dim - 1> &face_quadrature,
                const BDF1DiscreteTime    &time_handler)
      : fe_values(mapping,
                  fe,
                  cell_quadrature,
                  update_values | update_gradients | update_JxW_values)
      , fe_face_values(mapping,
                       fe,
                       face_quadrature,
                       update_values | update_gradients | update_JxW_values)
      , n_q_points(cell_quadrature.size())
      , n_faces(fe.reference_cell().n_faces())
      , n_faces_q_points(face_quadrature.size())
      , dofs_per_cell(fe.dofs_per_cell)
      , time_handler(time_handler)
    {
      allocate();
    }

    // Copy constructor, needed by the WorkStream::run() routine.
    ScratchData(const ScratchData &other)
      : fe_values(other.fe_values.get_mapping(),
                  other.fe_values.get_fe(),
                  other.fe_values.get_quadrature(),
                  other.fe_values.get_update_flags())
      , fe_face_values(other.fe_face_values.get_mapping(),
                       other.fe_face_values.get_fe(),
                       other.fe_face_values.get_quadrature(),
                       other.fe_face_values.get_update_flags())
      , n_q_points(other.n_q_points)
      , n_faces(other.n_faces)
      , n_faces_q_points(other.n_faces_q_points)
      , dofs_per_cell(other.dofs_per_cell)
      , time_handler(other.time_handler)
    {
      allocate();
    }

  private:
    // Allocate (resize and zero) the various vectors.
    void allocate()
    {
      velocity.first_vector_component = 0;
      pressure.component              = dim;
      lambda.first_vector_component   = dim + 1;

      JxW.resize(n_q_points);
      face_JxW.resize(n_faces, std::vector<double>(n_faces_q_points));

      // Navier-Stokes-related fields at time n
      present_velocity_values.resize(n_q_points);
      present_velocity_gradients.resize(n_q_points);
      present_velocity_sym_gradients.resize(n_q_points);
      present_velocity_divergence.resize(n_q_points);
      present_pressure_values.resize(n_q_points);
      previous_velocity_values.resize(time_handler.n_previous_sol,
                                      std::vector<Tensor<1, dim>>(n_q_points));
      present_velocity_time_derivatives.resize(n_q_points);
      present_face_velocity_values.resize(
        n_faces, std::vector<Tensor<1, dim>>(n_faces_q_points));

      // Navier-Stokes-related shape functions
      phi_u.resize(n_q_points, std::vector<Tensor<1, dim>>(dofs_per_cell));
      grad_phi_u.resize(n_q_points, std::vector<Tensor<2, dim>>(dofs_per_cell));
      sym_grad_phi_u.resize(
        n_q_points, std::vector<SymmetricTensor<2, dim>>(dofs_per_cell));
      div_phi_u.resize(n_q_points, std::vector<double>(dofs_per_cell));
      phi_p.resize(n_q_points, std::vector<double>(dofs_per_cell));
      phi_u_face.resize(n_faces,
                        std::vector<std::vector<Tensor<1, dim>>>(
                          n_faces_q_points,
                          std::vector<Tensor<1, dim>>(dofs_per_cell)));

      // Lagrange multiplier and its shape functions on faces
      present_face_lambda_values.resize(
        n_faces, std::vector<Tensor<1, dim>>(n_faces_q_points));
      phi_l_face.resize(n_faces,
                        std::vector<std::vector<Tensor<1, dim>>>(
                          n_faces_q_points,
                          std::vector<Tensor<1, dim>>(dofs_per_cell)));
    }

    // Fill the Navier-Stokes field vectors with volume data
    template <typename VectorType>
    void reinit_cell(const VectorType              &current_solution,
                     const std::vector<VectorType> &previous_solutions)
    {
      fe_values[velocity].get_function_values(current_solution,
                                              present_velocity_values);
      fe_values[velocity].get_function_gradients(current_solution,
                                                 present_velocity_gradients);
      fe_values[velocity].get_function_symmetric_gradients(
        current_solution, present_velocity_sym_gradients);
      fe_values[velocity].get_function_divergences(current_solution,
                                                   present_velocity_divergence);
      fe_values[pressure].get_function_values(current_solution,
                                              present_pressure_values);

      // Previous solutions
      for (unsigned int i = 0; i < previous_solutions.size(); ++i)
        fe_values[velocity].get_function_values(previous_solutions[i],
                                                previous_velocity_values[i]);

      for (unsigned int q = 0; q < n_q_points; ++q)
        {
          JxW[q] = fe_values.JxW(q);

          // Time derivatives
          present_velocity_time_derivatives[q] =
            time_handler.compute_time_derivative(q,
                                                 present_velocity_values,
                                                 previous_velocity_values);

          for (unsigned int k = 0; k < dofs_per_cell; ++k)
            {
              phi_u[q][k]          = fe_values[velocity].value(k, q);
              grad_phi_u[q][k]     = fe_values[velocity].gradient(k, q);
              sym_grad_phi_u[q][k] = symmetrize(grad_phi_u[q][k]);
              div_phi_u[q][k]      = fe_values[velocity].divergence(k, q);
              phi_p[q][k]          = fe_values[pressure].value(k, q);
            }
        }
    }

    // Fill the Navier-Stokes and Lagrange multiplier vectors with face data
    template <typename VectorType>
    void reinit_face(const unsigned int i_face,
                     const VectorType  &current_solution)
    {
      fe_face_values[velocity].get_function_values(
        current_solution, present_face_velocity_values[i_face]);
      fe_face_values[lambda].get_function_values(
        current_solution, present_face_lambda_values[i_face]);

      for (unsigned int q = 0; q < n_faces_q_points; ++q)
        {
          face_JxW[i_face][q] = fe_face_values.JxW(q);
          for (unsigned int k = 0; k < dofs_per_cell; ++k)
            {
              phi_u_face[i_face][q][k] = fe_face_values[velocity].value(k, q);
              fe_face_values[velocity].divergence(k, q);
              phi_l_face[i_face][q][k] = fe_face_values[lambda].value(k, q);
            }
        }
    }

  public:
    // Public reinit function to call when setting up the ScratchData on a cell.
    template <typename VectorType>
    void reinit(const typename DoFHandler<dim>::active_cell_iterator &cell,
                const VectorType              &current_solution,
                const std::vector<VectorType> &previous_solutions)
    {
      bdf_c0 = time_handler.get_bdf_coefficients()[0];

      // Volume contributions
      fe_values.reinit(cell);
      reinit_cell(current_solution, previous_solutions);

      // Face contributions
      if (cell->at_boundary())
        for (const auto i_face : cell->face_indices())
          {
            const auto &face = cell->face(i_face);
            if (face->at_boundary())
              {
                fe_face_values.reinit(cell, i_face);
                reinit_face(i_face, current_solution);
              }
          }
    }

  public:
    FEValues<dim>     fe_values;
    FEFaceValues<dim> fe_face_values;

    const unsigned int n_q_points;
    const unsigned int n_faces;
    const unsigned int n_faces_q_points;
    const unsigned int dofs_per_cell;

    const BDF1DiscreteTime &time_handler;

    std::vector<double>              JxW;
    std::vector<std::vector<double>> face_JxW;

    // First of the BDF coefficients
    double bdf_c0;

    FEValuesExtractors::Vector velocity;
    FEValuesExtractors::Scalar pressure;
    FEValuesExtractors::Vector lambda;

    // Navier-Stokes fields and shape functions.
    // Volume data, stored at each quadrature node
    std::vector<Tensor<1, dim>>              present_velocity_values;
    std::vector<Tensor<2, dim>>              present_velocity_gradients;
    std::vector<SymmetricTensor<2, dim>>     present_velocity_sym_gradients;
    std::vector<double>                      present_velocity_divergence;
    std::vector<Tensor<1, dim>>              present_velocity_time_derivatives;
    std::vector<double>                      present_pressure_values;
    std::vector<std::vector<Tensor<1, dim>>> previous_velocity_values;

    // Face data, stored at each face and each quadrature node
    std::vector<std::vector<Tensor<1, dim>>> present_face_velocity_values;

    // Volume shape functions, stored at each quadrature node and each dof
    std::vector<std::vector<Tensor<1, dim>>>          phi_u;
    std::vector<std::vector<Tensor<2, dim>>>          grad_phi_u;
    std::vector<std::vector<SymmetricTensor<2, dim>>> sym_grad_phi_u;
    std::vector<std::vector<double>>                  div_phi_u;
    std::vector<std::vector<double>>                  phi_p;
    std::vector<std::vector<Tensor<1, dim>>>          grad_phi_p;

    // Face shape functions, stored at each face, quadrature node, and dof
    std::vector<std::vector<std::vector<Tensor<1, dim>>>> phi_u_face;

    // Lagrange multiplier data, stored similarly
    std::vector<std::vector<Tensor<1, dim>>> present_face_lambda_values;
    std::vector<std::vector<std::vector<Tensor<1, dim>>>> phi_l_face;
  };

  class CopyData
  {
  public:
    CopyData(const unsigned int n_dofs_per_cell)
      : local_matrix(n_dofs_per_cell, n_dofs_per_cell)
      , local_rhs(n_dofs_per_cell)
      , local_dof_indices(n_dofs_per_cell)
    {}

  public:
    FullMatrix<double>                   local_matrix;
    Vector<double>                       local_rhs;
    std::vector<types::global_dof_index> local_dof_indices;
    bool                                 cell_is_locally_owned;
  };

  // @sect3{Parameters}
  //
  // A small struct to hold the simulation parameters.
  struct Parameters
  {
    const unsigned int velocity_degree = 2;
    const unsigned int pressure_degree = 1;
    const unsigned int lambda_degree   = 2;

    const unsigned int mapping_degree = 2;

    const double density             = 1.;
    const double kinematic_viscosity = 0.005;

    const double start_time = 0.;
    const double end_time   = 1.;
    const double time_step  = 0.05;

    const unsigned int max_iterations       = 20;
    const double       tolerance            = 1e-10;
    const double       divergence_tolerance = 1e4;
  };

  // @sect3{The <code>NavierStokesWithWeakNoSlip<code> class template}
  //
  // The main class is similar to the main class in step-57, which solves the
  // steady Navier-Stokes equations using Newton's method. The
  // principal differences are:
  //
  // - Instead of block preconditioning, we simply use a direct solver (Mumps).
  // Since we additionally solve the problem in parallel, we use PETSc
  // <code>LA::MPI::Vector</code>s instead of deal.II's serial
  // <code>BlockVector</code>s.
  // - Since we solve the <i>unsteady</i> Navier-Stokes equations, we start from
  // an initial velocity condition rather than from an educated initial guess,
  // as one would when solving the steady-state solution (e.g., from the
  // solution of the Stokes equations for the same boundary conditions).
  // - We need to store additional constraints for the Lagrange multiplier, to
  // specify that the lambda dofs not on the cylinder will be set to zero.
  // - For clarity, the set up functions for numbering the degrees of freedom,
  // creating the boundary conditions and Lagrange multiplier constraints, and
  // creating the sparsity pattern are split into dedicated functions.
  template <int dim>
  class NavierStokesWithWeakNoSlip
  {
  public:
    NavierStokesWithWeakNoSlip();
    void run();

  private:
    void create_grid();
    void setup_dofs();
    void create_lagrange_multiplier_constraints();
    void create_constraints(const bool                 homogeneous,
                            AffineConstraints<double> &constraints);
    void create_sparsity_pattern();
    void set_initial_conditions();

    void assemble_rhs();
    void assemble_local_rhs(
      const typename DoFHandler<dim>::active_cell_iterator &cell,
      ScratchData<dim>                                     &scratch_data,
      CopyData                                             &copy_data);
    void copy_local_to_global_rhs(const CopyData &copy_data);

    void assemble_matrix();
    void assemble_local_matrix(
      const typename DoFHandler<dim>::active_cell_iterator &cell,
      ScratchData<dim>                                     &scratch_data,
      CopyData                                             &copy_data);
    void copy_local_to_global_matrix(const CopyData &copy_data);

    void solve_linear_system();
    void solve_nonlinear_problem();

    void output_results();
    void check_no_slip_constraint();
    void compute_forces_on_cylinder();

  public:
    // These are constants and callbacks defined for convenience, which store
    // the ordering of the unknown variables in the finite element system
    static constexpr unsigned int n_components = 2 * dim + 1;
    static constexpr unsigned int u_lower      = 0;
    static constexpr unsigned int p_lower      = dim;
    static constexpr unsigned int l_lower      = dim + 1;

    inline bool is_velocity(const unsigned int component) const
    {
      return u_lower <= component && component < u_lower + dim;
    }
    inline bool is_pressure(const unsigned int component) const
    {
      return p_lower == component;
    }
    inline bool is_lambda(const unsigned int component) const
    {
      return l_lower <= component && component < l_lower + dim;
    }

  public:
    MPI_Comm mpi_communicator;

    const Parameters param;

    parallel::distributed::Triangulation<dim> triangulation;
    DoFHandler<dim>                           dof_handler;
    FESystem<dim>                             fe;
    MappingQ<dim>                             mapping;

    QGauss<dim>     cell_quadrature;
    QGauss<dim - 1> face_quadrature;

    IndexSet locally_owned_dofs;
    IndexSet locally_relevant_dofs;

    LA::MPI::SparseMatrix system_matrix;

    // Parallel vectors with ghosts
    LA::MPI::Vector              present_solution;
    std::vector<LA::MPI::Vector> previous_solutions;
    LA::MPI::Vector              evaluation_point;

    // Parallel, fully distributed vectors (without ghosts)
    LA::MPI::Vector system_rhs;
    LA::MPI::Vector local_evaluation_point;
    LA::MPI::Vector newton_update;

    BDF1DiscreteTime time_handler;

    ConditionalOStream pcout;
    TimerOutput        computing_timer;

    std::vector<std::pair<double, std::string>> visualization_times_and_names;

    AffineConstraints<double> lambda_constraints;
    AffineConstraints<double> zero_constraints;
    AffineConstraints<double> nonzero_constraints;

    // The ID of the boundary on which the weakly enforced no-slip condition is
    // applied
    types::boundary_id weak_no_slip_boundary_id = numbers::invalid_unsigned_int;

    // Table to store the total force on the cylinder at each time step
    TableHandler forces_table;
  };

  // @sect3{The <code>NavierStokesWithWeakNoSlip</code> class implementation}

  // @sect4{NavierStokesWithWeakNoSlip::NavierStokesWithWeakNoSlip}

  // In addition to the standard MPI, <code>pcout</code>, and timer
  // initializations, the constructor initializes the <code>FESystem</code> with
  // the velocity, pressure, and Lagrange multiplier finite element spaces, the
  // (possibly higher-order) mapping, the quadrature rules, and the time
  // handler.
  //
  // Note that although we defined convenience constants for the ordering of the
  // variable, this is not enforced when creating the <code>FESystem</code>, so
  // this is something to be wary of if you experiment and change the ordering.
  template <int dim>
  NavierStokesWithWeakNoSlip<dim>::NavierStokesWithWeakNoSlip()
    : mpi_communicator(MPI_COMM_WORLD)
    , param(Parameters())
    , triangulation(mpi_communicator)
    , dof_handler(triangulation)
    , fe(FE_Q<dim>(param.velocity_degree) ^ dim,
         FE_Q<dim>(param.pressure_degree),
         FE_Q<dim>(param.lambda_degree) ^ dim)
    , mapping(param.mapping_degree)
    , cell_quadrature(param.velocity_degree + 1)
    , face_quadrature(param.velocity_degree + 1)
    , time_handler(/* bdf_order = */ 1,
                   param.start_time,
                   param.end_time,
                   param.time_step)
    , pcout(std::cout,
            (Utilities::MPI::this_mpi_process(mpi_communicator) == 0))
    , computing_timer(mpi_communicator,
                      pcout,
                      TimerOutput::summary,
                      TimerOutput::wall_times)
  {}

  // @sect4{NavierStokesWithWeakNoSlip<dim>::create_grid}

  // This function creates the mesh using the ready-made
  // <code>uniform_channel_with_cylinder</code> function from the GridGenerator
  // namespace. The obtained triangulation is further refined around the
  // cylinder, to better capture the geometry, and in the wake, for which a
  // basic rectangular box is used to assign refinement flags.
  template <int dim>
  void NavierStokesWithWeakNoSlip<dim>::create_grid()
  {
    TimerOutput::Scope t(computing_timer, "Create mesh");

    const std::vector<unsigned int> lengths_and_heights = {8, 24, 8, 8};

    // Use the default settings of the grid genarator, except that we would like
    // to color the boundary entities, to assign the boundary conditions.
    // This assigns the boundary id 2 to the cylinder.
    GridGenerator::uniform_channel_with_cylinder(
      triangulation, lengths_and_heights, 2, 3, 0.75, 2, 2, false, true);

    weak_no_slip_boundary_id = 2;

    if constexpr (dim == 3)
      {
        std::vector<GridTools::PeriodicFacePair<
          typename parallel::distributed::Triangulation<dim>::cell_iterator>>
          periodicity_vector;

        GridTools::collect_periodic_faces(triangulation,
                                          5,
                                          6,
                                          /* direction = */ 2,
                                          periodicity_vector);

        triangulation.add_periodicity(periodicity_vector);
      }

    // Refine the starting mesh around the cylinder, and in a rectangular wake
    // behind it.
    const Point<dim> center;
    const double     inner_radius = 1.5;
    for (unsigned int step = 0; step < ((dim == 2) ? 2 : 1); ++step)
      {
        for (const auto &cell : triangulation.active_cell_iterators())
          for (const auto v : cell->vertex_indices())
            {
              const Point<dim> p = cell->vertex(v);

              // Refine around the cylinder and in a rectangular wake
              if (center.distance(p) <= inner_radius ||
                  (p[0] > 0. && -1.5 <= p[1] && p[1] <= 1.5))
                {
                  cell->set_refine_flag();
                  break;
                }
            }

        triangulation.execute_coarsening_and_refinement();
      }
  }

  // @sect4{NavierStokesWithWeakNoSlip<dim>::setup_dofs}

  // This function numbers the degrees of freedom (dofs) on the mesh, and
  // extracts the list of locally owned and relevant dofs, which are then used
  // to initialize the PETSc parallel vectors.
  //
  // Because the sparsity pattern is slightly more involved than the standard
  // <code>DynamicSparsityPattern</code>, it is computed in its dedicated
  // function, where the system matrix is also initialized.
  template <int dim>
  void NavierStokesWithWeakNoSlip<dim>::setup_dofs()
  {
    TimerOutput::Scope t(computing_timer, "Setup dofs");

    auto &comm = mpi_communicator;

    dof_handler.distribute_dofs(fe);

    pcout << "Number of degrees of freedom: " << dof_handler.n_dofs()
          << std::endl;

    locally_owned_dofs = dof_handler.locally_owned_dofs();
    locally_relevant_dofs =
      DoFTools::extract_locally_relevant_dofs(dof_handler);

    // Initialize parallel vectors with ghosts
    present_solution.reinit(locally_owned_dofs, locally_relevant_dofs, comm);
    evaluation_point.reinit(locally_owned_dofs, locally_relevant_dofs, comm);

    // Initialize parallel vectors without ghosts
    local_evaluation_point.reinit(locally_owned_dofs, comm);
    newton_update.reinit(locally_owned_dofs, comm);
    system_rhs.reinit(locally_owned_dofs, comm);

    // Allocate for previous BDF solutions (which include ghosts)
    previous_solutions.resize(time_handler.n_previous_sol);
    for (auto &previous_sol : previous_solutions)
      previous_sol.reinit(locally_owned_dofs, locally_relevant_dofs, comm);
  }

  // @sect4{NavierStokesWithWeakNoSlip<dim>::create_lagrange_multiplier_constraints}

  // This functions creates the linear constraints for the Lagrange multiplier
  // (lambda). In the current case, we want to mimic a field defined only on a
  // codimension 1 entity (the cylinder), although the field is actually defined
  // in the whole volume. To this end, we identify the lambda dofs located on
  // the cylinder using <code>DoFTools::extract_boundary_dofs</code>, and
  // constrain all the other lambda dofs to zero.
  //
  // We would like the resulting constraints to be consistent in parallel,
  // however. For this exact problem, it is actually relatively easy, because
  // the constrained Lagrange multiplier dofs do not actually intervene in any
  // computation. Indeed, only the lambda dofs on the faces on the cylinder are
  // needed in the assembly routines, and the constrained volume dofs do not
  // participate. This is because we are using Lagrange finite elements, for
  // which the shape functions associated with the constrained volume dofs are
  // identically zero on the faces lying on the cylinder. Thus, there is no risk
  // of computing erroneous integrals because the constrained values of owned
  // and ghosted dofs are inconsistent. In general, however, we should make sure
  // that the constraints are consistent between owned and ghosted dofs. The
  // strategy used here is to constrain only the *owned* lambda dofs not on the
  // cylinder, and then let the <code>make_consistent_in_parallel</code>
  // function apply these constraints to the ghost dofs accordingly.
  template <int dim>
  void NavierStokesWithWeakNoSlip<dim>::create_lagrange_multiplier_constraints()
  {
    lambda_constraints.reinit(locally_owned_dofs, locally_relevant_dofs);

    const FEValuesExtractors::Vector lambda(l_lower);
    const ComponentMask              lambda_mask = fe.component_mask(lambda);

    // Get the relevant Lagrange multiplier dofs on the cylinder
    // (DoFTools::extract_boundary_dofs returns a list of *relevant* dofs)
    IndexSet relevant_boundary_dofs =
      DoFTools::extract_boundary_dofs(dof_handler,
                                      lambda_mask,
                                      {weak_no_slip_boundary_id});

    // Get all owned Lagrange multiplier dofs
    IndexSet owned_lambda_dofs =
      DoFTools::extract_dofs(dof_handler, lambda_mask);

    // Constrain the owned lambda dofs not on the boundary
    for (const auto dof : owned_lambda_dofs)
      if (!relevant_boundary_dofs.is_element(dof))
        lambda_constraints.constrain_dof_to_zero(dof);

    lambda_constraints.close();

    // Make the constraints consistent in parallel
    lambda_constraints.make_consistent_in_parallel(
      locally_owned_dofs,
      DoFTools::extract_locally_active_dofs(dof_handler),
      mpi_communicator);
  }

  // @sect4{NavierStokesWithWeakNoSlip<dim>::create_constraints}
  //
  // This function creates the whole set of constraints, which includes the
  // hanging node constraints, setting the Lagrange multiplier to zero outside
  // the cylinder, and enforcing the (non)homogeneous boundary conditions.
  //
  // The argument <code>homogeneous</code> determines whether the passed
  // <code>constraints</code> are the zero or nonzero constraints.
  template <int dim>
  void NavierStokesWithWeakNoSlip<dim>::create_constraints(
    const bool                 homogeneous,
    AffineConstraints<double> &constraints)
  {
    constraints.clear();
    constraints.reinit(locally_owned_dofs, locally_relevant_dofs);

    DoFTools::make_hanging_node_constraints(dof_handler, constraints);

    Functions::ZeroFunction<dim> zero_fun(n_components);

    // Velocity boundary conditions
    {
      const FEValuesExtractors::Vector velocity(u_lower);
      const ComponentMask velocity_mask = fe.component_mask(velocity);

      Inlet<dim> inlet_fun(u_lower, n_components);

      // There is nothing to do here for boundary 1 (outflow, natural
      // condition), and for boundary 2 (cylinder, no-slip is weakly enforced).
      types::boundary_id           inlet_boundary  = 0;
      std::set<types::boundary_id> slip_boundaries = {3, 4};

      // Inlet velocity
      {
        Function<dim> *fun_ptr;
        if (homogeneous)
          fun_ptr = &zero_fun;
        else
          fun_ptr = &inlet_fun;
        VectorTools::interpolate_boundary_values(mapping,
                                                 dof_handler,
                                                 inlet_boundary,
                                                 *fun_ptr,
                                                 constraints,
                                                 velocity_mask);
      }

      // Slip boundaries (zero normal flux)
      VectorTools::compute_no_normal_flux_constraints(
        dof_handler,
        u_lower,
        slip_boundaries,
        constraints,
        mapping,
        /*use_manifold_for_normal=*/false);

      if constexpr (dim == 3)
        {
          const FEValuesExtractors::Vector lambda(l_lower);
          const ComponentMask lambda_mask = fe.component_mask(lambda);

          std::vector<GridTools::PeriodicFacePair<
            typename DoFHandler<dim>::cell_iterator>>
            periodicity_vector;
          GridTools::collect_periodic_faces(dof_handler,
                                            5,
                                            6,
                                            /* direction = */ 2,
                                            periodicity_vector);

          // Set all fields as periodic (by not providing any component mask)
          DoFTools::make_periodicity_constraints<dim, dim>(periodicity_vector,
                                                           constraints);
        }
    }

    constraints.close();

    // Merge the zero lambda constraints.
    // Some lambda dofs will already be constrained by the hanging node
    // constraints, in which case there is a clash between both sets of
    // constraints. But since there is no hanging node on the cylinder, these
    // dofs will be constrained to zero anyway, so we can simply overwrite the
    // constraints by the lambda zero constraints.
    constraints.merge(
      lambda_constraints,
      AffineConstraints<double>::MergeConflictBehavior::right_object_wins);
  }

  // @sect4{NavierStokesWithWeakNoSlip<dim>::create_sparsity_pattern}

  // We use a coupling table to create the sparsity pattern. By setting
  // <code>keep_constrained_dofs</code> to false, we also eliminate the
  // couplings from all the constrained lambda dofs in the volume.
  template <int dim>
  void NavierStokesWithWeakNoSlip<dim>::create_sparsity_pattern()
  {
    DynamicSparsityPattern dsp(locally_relevant_dofs);

    Table<2, DoFTools::Coupling> coupling_table(n_components, n_components);
    for (unsigned int c = 0; c < n_components; ++c)
      for (unsigned int d = 0; d < n_components; ++d)
        {
          coupling_table[c][d] = DoFTools::none;

          // u couples to all variables
          if (is_velocity(c))
            coupling_table[c][d] = DoFTools::always;
          // p couples to u only
          else if (is_pressure(c) && is_velocity(d))
            coupling_table[c][d] = DoFTools::always;
          // lambda couples to u only
          else if (is_lambda(c) && is_velocity(d))
            coupling_table[c][d] = DoFTools::always;
        }

    DoFTools::make_sparsity_pattern(dof_handler,
                                    coupling_table,
                                    dsp,
                                    nonzero_constraints,
                                    /* keep_constrained_dofs = */ false);
    SparsityTools::distribute_sparsity_pattern(dsp,
                                               locally_owned_dofs,
                                               mpi_communicator,
                                               locally_relevant_dofs);
    system_matrix.reinit(locally_owned_dofs,
                         locally_owned_dofs,
                         dsp,
                         mpi_communicator);
  }

  // @sect4{NavierStokesWithWeakNoSlip<dim>::set_initial_conditions}

  // Next we handle the initial condition. Only the velocity has initial
  // conditions, as the pressure and Lagrange multiplier adjust to the
  // instantaneous velocity. There is not much to say here: the initial
  // condition is interpolated into the local_evaluation_point, then we apply
  // the nonhomogeneous constraints, and copy the local solution into the
  // ghosted vectors. The solution is copied into present_solution to allow
  // outputting the initial condition for visualization, and to
  // evaluation_point, as this is the vector used to compute the assembly during
  // the Newton iterations.
  template <int dim>
  void NavierStokesWithWeakNoSlip<dim>::set_initial_conditions()
  {
    const FEValuesExtractors::Vector velocity(u_lower);
    const ComponentMask velocity_mask = fe.component_mask(velocity);

    InitialVelocity<dim> initial_velocity(u_lower, n_components);
    VectorTools::interpolate(mapping,
                             dof_handler,
                             initial_velocity,
                             local_evaluation_point,
                             velocity_mask);

    // Apply non-homogeneous Dirichlet BC and set as current solution
    nonzero_constraints.distribute(local_evaluation_point);
    present_solution = local_evaluation_point;
    evaluation_point = local_evaluation_point;
  }

  // @sect4{NavierStokesWithWeakNoSlip<dim>::assemble_rhs}

  // This function assembles the right-hand side of a single Newton iteration.
  // Aside from resetting the global RHS, scratch and copy data, the function
  // itself does not do much, and simply forwards the work to the cell workers
  // <code>assemble_local_rhs</code> and <code>copy_local_to_global_rhs</code>.
  template <int dim>
  void NavierStokesWithWeakNoSlip<dim>::assemble_rhs()
  {
    TimerOutput::Scope t(computing_timer, "Assemble RHS");

    system_rhs = 0;

    ScratchData<dim> scratch_data(
      fe, mapping, cell_quadrature, face_quadrature, time_handler);
    CopyData copy_data(fe.n_dofs_per_cell());

    // Assemble RHS
    WorkStream::run(dof_handler.begin_active(),
                    dof_handler.end(),
                    *this,
                    &NavierStokesWithWeakNoSlip::assemble_local_rhs,
                    &NavierStokesWithWeakNoSlip::copy_local_to_global_rhs,
                    scratch_data,
                    copy_data);
    system_rhs.compress(VectorOperation::add);
  }

  // @sect4{NavierStokesWithWeakNoSlip<dim>::assemble_local_rhs}

  // This worker function computes the local Newton RHS on an owned cell, and
  // fills the local_rhs in the passed <code>CopyData</code>. This is done by
  // first assembling the volume integrals, then the face integrals on the faces
  // lying on the cylinder.
  template <int dim>
  void NavierStokesWithWeakNoSlip<dim>::assemble_local_rhs(
    const typename DoFHandler<dim>::active_cell_iterator &cell,
    ScratchData<dim>                                     &scratch_data,
    CopyData                                             &copy_data)
  {
    copy_data.cell_is_locally_owned = cell->is_locally_owned();
    if (!cell->is_locally_owned())
      return;

    cell->get_dof_indices(copy_data.local_dof_indices);

    auto &sd = scratch_data;
    sd.reinit(cell, evaluation_point, previous_solutions);

    auto &local_rhs = copy_data.local_rhs;
    local_rhs       = 0;

    // Cell assembly
    for (unsigned int q = 0; q < sd.n_q_points; ++q)
      {
        // Navier-Stokes data
        const double JxW        = sd.JxW[q];
        const auto  &dudt       = sd.present_velocity_time_derivatives[q];
        const auto  &u          = sd.present_velocity_values[q];
        const auto  &grad_u     = sd.present_velocity_gradients[q];
        const auto  &sym_grad_u = sd.present_velocity_sym_gradients[q];
        const double div_u      = sd.present_velocity_divergence[q];
        const double p          = sd.present_pressure_values[q];

        const auto to_multiply_by_phi_u_i = dudt + grad_u * u;

        const auto &phi_u          = sd.phi_u[q];
        const auto &sym_grad_phi_u = sd.sym_grad_phi_u[q];
        const auto &div_phi_u      = sd.div_phi_u[q];
        const auto &phi_p          = sd.phi_p[q];

        for (unsigned int i = 0; i < sd.dofs_per_cell; ++i)
          {
            // Time derivative, convective acceleration, pressure gradient,
            // and diffusion
            local_rhs(i) -=
              (phi_u[i] * to_multiply_by_phi_u_i - div_phi_u[i] * p +
               2. * param.kinematic_viscosity *
                 scalar_product(sym_grad_u, sym_grad_phi_u[i])) *
              JxW;

            // Mass equation
            local_rhs(i) -= phi_p[i] * (-div_u) * JxW;
          }
      }

    // Face assembly
    if (cell->at_boundary())
      for (unsigned int i_face = 0; i_face < sd.n_faces; ++i_face)
        {
          const auto &face = cell->face(i_face);
          if (face->at_boundary() &&
              face->boundary_id() == weak_no_slip_boundary_id)
            for (unsigned int q = 0; q < sd.n_faces_q_points; ++q)
              {
                const double face_JxW = sd.face_JxW[i_face][q];
                const auto  &phi_u    = sd.phi_u_face[i_face][q];
                const auto  &phi_l    = sd.phi_l_face[i_face][q];
                const auto  &lambda = sd.present_face_lambda_values[i_face][q];

                // Compute the no-slip velocity constraint.
                // This is typically u - g = 0, with g the prescribed velocity
                // field, which is 0 for this example.
                const auto &fluid_velocity =
                  sd.present_face_velocity_values[i_face][q];
                const auto velocity_constraint = fluid_velocity;

                for (unsigned int i = 0; i < sd.dofs_per_cell; ++i)
                  local_rhs(i) -=
                    (-phi_u[i] * lambda + -velocity_constraint * phi_l[i]) *
                    face_JxW;
              }
        }
  }

  // @sect4{NavierStokesWithWeakNoSlip<dim>::copy_local_to_global_rhs}

  // After having filled the local rhs on an owned cell, this function is
  // called by WorkStream::run to copy its content into the global residual.
  // Since Newton increments must have homogeneous boundary conditions,
  // the set of constraints used to distribute the local rhs is the *zero*
  // constraints.
  template <int dim>
  void NavierStokesWithWeakNoSlip<dim>::copy_local_to_global_rhs(
    const CopyData &copy_data)
  {
    if (copy_data.cell_is_locally_owned)
      zero_constraints.distribute_local_to_global(copy_data.local_rhs,
                                                  copy_data.local_dof_indices,
                                                  system_rhs);
  }

  // @sect4{NavierStokesWithWeakNoSlip<dim>::assemble_matrix}

  // This function is basically identical to <code>assemble_rhs</code>.
  template <int dim>
  void NavierStokesWithWeakNoSlip<dim>::assemble_matrix()
  {
    TimerOutput::Scope t(computing_timer, "Assemble matrix");

    system_matrix = 0;

    ScratchData<dim> scratch_data(
      fe, mapping, cell_quadrature, face_quadrature, time_handler);
    CopyData copy_data(fe.n_dofs_per_cell());

    // Assemble matrix
    WorkStream::run(dof_handler.begin_active(),
                    dof_handler.end(),
                    *this,
                    &NavierStokesWithWeakNoSlip::assemble_local_matrix,
                    &NavierStokesWithWeakNoSlip::copy_local_to_global_matrix,
                    scratch_data,
                    copy_data);
    system_matrix.compress(VectorOperation::add);
  }

  // @sect4{NavierStokesWithWeakNoSlip<dim>::assemble_local_matrix}

  // This function assembles the local matrix on an owned cell, and follows the
  // same structure as its rhs counterpart, assembling first the volume
  // integrals, then the face integrals on the cylinder.
  template <int dim>
  void NavierStokesWithWeakNoSlip<dim>::assemble_local_matrix(
    const typename DoFHandler<dim>::active_cell_iterator &cell,
    ScratchData<dim>                                     &scratch_data,
    CopyData                                             &copy_data)
  {
    copy_data.cell_is_locally_owned = cell->is_locally_owned();
    if (!cell->is_locally_owned())
      return;

    cell->get_dof_indices(copy_data.local_dof_indices);

    auto &sd = scratch_data;
    sd.reinit(cell, evaluation_point, previous_solutions);

    auto &local_matrix = copy_data.local_matrix;
    local_matrix       = 0;

    std::vector<Tensor<1, dim>> to_multiply_by_phi_u_i_momentum(
      sd.dofs_per_cell);

    // Cell assembly
    for (unsigned int q = 0; q < sd.n_q_points; ++q)
      {
        // Navier-Stokes data
        const double JxW    = sd.JxW[q];
        const auto  &u      = sd.present_velocity_values[q];
        const auto  &grad_u = sd.present_velocity_gradients[q];

        const auto &phi_u          = sd.phi_u[q];
        const auto &grad_phi_u     = sd.grad_phi_u[q];
        const auto &sym_grad_phi_u = sd.sym_grad_phi_u[q];
        const auto &div_phi_u      = sd.div_phi_u[q];
        const auto &phi_p          = sd.phi_p[q];

        // Precompute quantities depending only on j
        for (unsigned int j = 0; j < sd.dofs_per_cell; ++j)
          to_multiply_by_phi_u_i_momentum[j] =
            sd.bdf_c0 * phi_u[j] + grad_phi_u[j] * u + grad_u * phi_u[j];

        for (unsigned int i = 0; i < sd.dofs_per_cell; ++i)
          for (unsigned int j = 0; j < sd.dofs_per_cell; ++j)
            local_matrix(i, j) +=
              (-div_phi_u[i] * phi_p[j] +
               phi_u[i] * to_multiply_by_phi_u_i_momentum[j] +
               2. * param.kinematic_viscosity *
                 scalar_product(sym_grad_phi_u[j], sym_grad_phi_u[i]) -
               phi_p[i] * div_phi_u[j]) *
              JxW;
      }

    // Face assembly
    if (cell->at_boundary())
      for (unsigned int i_face = 0; i_face < sd.n_faces; ++i_face)
        {
          const auto &face = cell->face(i_face);
          if (face->at_boundary() &&
              face->boundary_id() == weak_no_slip_boundary_id)
            for (unsigned int q = 0; q < sd.n_faces_q_points; ++q)
              {
                const double face_JxW = sd.face_JxW[i_face][q];
                const auto  &phi_u    = sd.phi_u_face[i_face][q];
                const auto  &phi_l    = sd.phi_l_face[i_face][q];

                for (unsigned int i = 0; i < sd.dofs_per_cell; ++i)
                  for (unsigned int j = 0; j < sd.dofs_per_cell; ++j)
                    local_matrix(i, j) +=
                      (-phi_l[j] * phi_u[i] - phi_u[j] * phi_l[i]) * face_JxW;
              }
        }
  }

  // @sect4{NavierStokesWithWeakNoSlip<dim>::copy_local_to_global_matrix}

  // This function is, once again, functionally identical to its right-hand side
  // counterpart.
  template <int dim>
  void NavierStokesWithWeakNoSlip<dim>::copy_local_to_global_matrix(
    const CopyData &copy_data)
  {
    if (copy_data.cell_is_locally_owned)
      zero_constraints.distribute_local_to_global(copy_data.local_matrix,
                                                  copy_data.local_dof_indices,
                                                  system_matrix);
  }

  // @sect4{NavierStokesWithWeakNoSlip<dim>::solve_linear_system}
  //
  // The function to solve the linear system at each Newton iteration is very
  // simple, as we use Mumps as a silver bullet solver, which avoids the problem
  // of finding a good preconditioner for the monolithic matrix.
  // This is, of course, viable for relatively small problems, typically in two
  // dimensions, but quickly becomes limited in three dimension and for larger
  // dof counts.
  template <int dim>
  void NavierStokesWithWeakNoSlip<dim>::solve_linear_system()
  {
    TimerOutput::Scope t(computing_timer, "Solve direct");

    LA::MPI::Vector completely_distributed_solution(locally_owned_dofs,
                                                    mpi_communicator);

    // Solve with MUMPS
    SolverControl                    solver_control;
    PETScWrappers::SparseDirectMUMPS linear_solver(solver_control);

    linear_solver.solve(system_matrix,
                        completely_distributed_solution,
                        system_rhs);

    newton_update = completely_distributed_solution;
    zero_constraints.distribute(newton_update);
  }

  // @sect4{NavierStokesWithWeakNoSlip<dim>::solve_nonlinear_problem}

  // This function solves the coupled nonlinear problem at a given time step
  // $t^n$, and iterates to find the triplet $(u^n, p^n, \lambda^n)$. The
  // implementation is inspired from <code>newton_iteration</code> in the
  // step-57 tutorial, except that the present function is a "raw" Newton's
  // method, without line search. In the present implementation, the matrix is
  // assembled at each iteration, but in practice, heuristics can be used to
  // limit the number of re-assemblies at the price of a few additional, cheaper
  // iterations.
  template <int dim>
  void NavierStokesWithWeakNoSlip<dim>::solve_nonlinear_problem()
  {
    bool         stop          = false;
    unsigned int iter          = 0;
    double       norm_residual = 0;

    evaluation_point = present_solution;

    while (!stop)
      {
        // Assemble residual and check if tolerance is reached
        assemble_rhs();
        norm_residual = system_rhs.l2_norm();

        pcout << "Newton iter. " << std::setw(2) << iter << ": "
              << std::scientific << std::setprecision(8)
              << " nonlinear residual = " << norm_residual << std::endl;

        // Abort if residual is too high
        const bool diverged =
          iter > 0 && norm_residual > param.divergence_tolerance;
        AssertThrow(!diverged, ExcMessage("Nonlinear solver diverged!"));

        if (norm_residual <= param.tolerance)
          {
            pcout << "Stopping because residual is below prescribed tolerance ("
                  << std::setprecision(2) << param.tolerance << ")"
                  << std::endl;
            break;
          }

        // Assemble matrix and solve
        assemble_matrix();
        solve_linear_system();

        // Increment solution, apply nonhomogeneous constraints, and go to next
        // iteration
        local_evaluation_point = evaluation_point;
        local_evaluation_point.add(1., newton_update);
        nonzero_constraints.distribute(local_evaluation_point);
        evaluation_point = local_evaluation_point;

        if (++iter > param.max_iterations)
          stop = true;
      }

    const bool solution_found = norm_residual <= param.tolerance;
    AssertThrow(solution_found,
                ExcMessage("Nonlinear solver did not converge"));

    // Update present solution
    present_solution = evaluation_point;
  }

  // @sect4{NavierStokesWithWeakNoSlip<dim>::output_results}

  // The function to output the results is relatively standard. Its
  // specificities are:
  //
  // - As for most MPI examples, the partition ID (subdomain) of the cells is
  // exported in addition to the solved fields.
  // - If the mapping is high-order, high-order VTK cells are written instead of
  // subdividing the mesh for visualization.
  // - A single VTU file is written, even with multiple MPI ranks
  // - A PVD file is created at each time step: this allows only opening this
  // file to view the unsteady solution, without pesky mismatched VTU files from
  // previous runs.
  template <int dim>
  void NavierStokesWithWeakNoSlip<dim>::output_results()
  {
    TimerOutput::Scope t(computing_timer, "Write outputs");

    std::string output_directory = "./results/";

    // ID of the partition
    Vector<float> subdomain(triangulation.n_active_cells());
    for (unsigned int i = 0; i < subdomain.size(); ++i)
      subdomain(i) = triangulation.locally_owned_subdomain();

    std::vector<std::string> solution_names(dim, "velocity");
    solution_names.push_back("pressure");
    for (unsigned int d = 0; d < dim; ++d)
      solution_names.push_back("lagrange_multiplier");

    std::vector<DataComponentInterpretation::DataComponentInterpretation>
      data_component_interpretation(
        dim, DataComponentInterpretation::component_is_part_of_vector);
    data_component_interpretation.push_back(
      DataComponentInterpretation::component_is_scalar);
    for (unsigned int d = 0; d < dim; ++d)
      data_component_interpretation.push_back(
        DataComponentInterpretation::component_is_part_of_vector);

    DataOut<dim> data_out;

    // Write high-order elements if needed
    if (param.mapping_degree > 1)
      {
        DataOutBase::VtkFlags flags;
        flags.write_higher_order_cells = true;
        data_out.set_flags(flags);
      }

    data_out.attach_dof_handler(dof_handler);
    data_out.add_data_vector(present_solution,
                             solution_names,
                             DataOut<dim>::type_dof_data,
                             data_component_interpretation);
    data_out.add_data_vector(subdomain, "subdomain");
    data_out.build_patches(mapping, 2, DataOut<dim>::curved_inner_cells);

    const std::string pvtu_file =
      data_out.write_vtu_with_pvtu_record(output_directory,
                                          "solution",
                                          time_handler.get_step_number(),
                                          mpi_communicator,
                                          /* n_digits_for_counter = */ 4,
                                          /* n_groups = */ 1);
    visualization_times_and_names.emplace_back(time_handler.get_current_time(),
                                               pvtu_file);
    std::ofstream pvd_output(output_directory + "solution.pvd");
    DataOutBase::write_pvd_record(pvd_output, visualization_times_and_names);
  }

  // @sect4{NavierStokesWithWeakNoSlip<dim>::check_no_slip_constraint}

  // Since we are enforcing the no-slip constraint weakly, it is a good idea
  // to quantiy how strictly this constraint is effectively enforced. The
  // function below loops over the cylinder faces and evaluates the constraint,
  // that is, evaluates the velocity norm on these faces, and measures how it
  // deviates from zero.
  //
  // When the polynomial degree of the Lagrange multiplier is lower than the
  // velocity, the accuracy of the no-slip constraint is mesh-dependent. When
  // their degrees are equal, there is a one-to-one match between their dofs,
  // and one shows that the constraint is enforced exactly (that is, up to
  // the precision of the linear solver), however it comes with inf-sup
  // considerations.
  template <int dim>
  void NavierStokesWithWeakNoSlip<dim>::check_no_slip_constraint()
  {
    double l2_local = 0, li_local = 0;

    FEFaceValues<dim> fe_face_values(mapping,
                                     fe,
                                     face_quadrature,
                                     update_values | update_JxW_values);

    const unsigned int               n_faces_q_points = face_quadrature.size();
    std::vector<Tensor<1, dim>>      velocity_values(n_faces_q_points);
    const FEValuesExtractors::Vector velocity(u_lower);

    for (auto cell : dof_handler.active_cell_iterators())
      if (cell->is_locally_owned())
        for (const auto i_face : cell->face_indices())
          {
            const auto &face = cell->face(i_face);
            if (face->at_boundary() &&
                face->boundary_id() == weak_no_slip_boundary_id)
              {
                fe_face_values.reinit(cell, i_face);
                fe_face_values[velocity].get_function_values(present_solution,
                                                             velocity_values);
                for (unsigned int q = 0; q < n_faces_q_points; ++q)
                  {
                    Tensor<1, dim> constraint = velocity_values[q];

                    // Measure constraint enforcement (target is zero)
                    l2_local += constraint * constraint * fe_face_values.JxW(q);
                    li_local = std::max(li_local, constraint.norm());
                  }
              }
          }

    const double l2_error =
      std::sqrt(Utilities::MPI::sum(l2_local, mpi_communicator));
    const double li_error = Utilities::MPI::max(li_local, mpi_communicator);

    if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
      {
        std::cout << "Checking no-slip enforcement on cylinder:" << std::endl;
        std::cout << "||uh||_L2   = " << l2_error << std::endl;
        std::cout << "||uh||_Linf = " << li_error << std::endl;
      }
  }

  // @sect4{NavierStokesWithWeakNoSlip<dim>::compute_forces_on_cylinder}

  // This last postprocessing function evaluates the total fluid forces on the
  // cylinder, which is simply the negative of the integral of the Lagrange
  // multiplier. Each force component is written to a table, which is written to
  // a file at each time step.
  template <int dim>
  void NavierStokesWithWeakNoSlip<dim>::compute_forces_on_cylinder()
  {
    FEFaceValues<dim> fe_face_values(mapping,
                                     fe,
                                     face_quadrature,
                                     update_values | update_JxW_values);

    const unsigned int               n_faces_q_points = face_quadrature.size();
    std::vector<Tensor<1, dim>>      lambda_values(n_faces_q_points);
    const FEValuesExtractors::Vector lambda(l_lower);

    Tensor<1, dim> lambda_integral_local;
    for (auto cell : dof_handler.active_cell_iterators())
      if (cell->is_locally_owned())
        for (unsigned int i_face = 0; i_face < cell->n_faces(); ++i_face)
          {
            const auto &face = cell->face(i_face);
            if (face->at_boundary() &&
                face->boundary_id() == weak_no_slip_boundary_id)
              {
                fe_face_values.reinit(cell, i_face);
                fe_face_values[lambda].get_function_values(present_solution,
                                                           lambda_values);
                for (unsigned int q = 0; q < n_faces_q_points; ++q)
                  lambda_integral_local +=
                    lambda_values[q] * fe_face_values.JxW(q);
              }
          }
    // Reduce each component, and take the negative to get the force
    Tensor<1, dim> forces;
    for (unsigned int d = 0; d < dim; ++d)
      forces[d] = -param.density * Utilities::MPI::sum(lambda_integral_local[d],
                                                       mpi_communicator);

    // Write forces to table
    std::vector<std::string> dim_str = {"x", "y", "z"};
    forces_table.add_value("time", time_handler.get_current_time());
    for (unsigned int d = 0; d < dim; ++d)
      {
        forces_table.add_value("F" + dim_str[d], forces[d]);
        forces_table.set_precision("F" + dim_str[d], 3);
        forces_table.set_scientific("F" + dim_str[d], true);
      }
    if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
      {
        std::ofstream out("./results/forces.txt");
        out << std::scientific << std::setprecision(3);
        forces_table.write_text(out);
      }
  }

  // @sect4{NavierStokesWithWeakNoSlip<dim>::run}

  // The run function sets up the problem and controls the time integration
  // loop. After applying the initial condition, the solutions are rotated, so
  // that the previous solution is set to the initial condition as we enter the
  // time loop.
  //
  // The initial condition is postprocessed as well to store information at t =
  // 0, although at that point the Lagrange multipliers have not yet been
  // computed, thus the no-slip is not satisfied (it is the norm of the initial
  // velocity on the cylinder), and the fluid forces are zero.
  //
  // Note that the inhomogeneous boundary conditions (the nonzero constraints)
  // are created only once. This is valid only if these conditions are time
  // independent, otherwise we should re-create those in the time integration
  // loop, after advancing the time.
  template <int dim>
  void NavierStokesWithWeakNoSlip<dim>::run()
  {
    pcout << "Running on " << Utilities::MPI::n_mpi_processes(mpi_communicator)
          << " MPI rank(s)..." << std::endl;

    create_grid();
    setup_dofs();
    create_lagrange_multiplier_constraints();
    create_constraints(/* homogeneous = */ true, zero_constraints);
    create_constraints(/* homogeneous = */ false, nonzero_constraints);
    create_sparsity_pattern();
    set_initial_conditions();
    time_handler.rotate_solutions(present_solution, previous_solutions);

    // Postprocess initial condition
    output_results();
    check_no_slip_constraint();
    compute_forces_on_cylinder();

    while (!time_handler.is_at_end())
      {
        time_handler.advance_time();

        pcout << std::endl;
        pcout << "Time step " << time_handler.get_step_number()
              << " : time = " << time_handler.get_current_time() << std::endl;

        solve_nonlinear_problem();
        output_results();
        check_no_slip_constraint();
        compute_forces_on_cylinder();

        time_handler.rotate_solutions(present_solution, previous_solutions);
      }
  }
} // namespace Step106

int main(int argc, char *argv[])
{
  try
    {
      using namespace dealii;
      using namespace Step106;

      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

      NavierStokesWithWeakNoSlip<2> flow;
      flow.run();
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
