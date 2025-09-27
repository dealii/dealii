// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2022 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

//
// Description:
//
// A performance benchmark based on step 67 and step 76, extended to the
// compressible Navier-Stokes equations. We measure timings for grid creation,
// setup of unknowns, explicit Runge-Kutta time stepping with face-centric
// loop through MatrixFree::loop() (step-67 style) as well as cell-centric
// loop through MatrixFree::loop_cell_centric() (step-76 style). We use a
// problem with periodic boundary conditions to avoid defining complicated
// definitions of boundary data and inflow/outflow treatment.
//
// Status: stable
//
// Note: this test is marked "stable" and used for performance
// instrumentation in our testsuite, https://dealii.org/performance_tests
//

#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/time_stepping.h>
#include <deal.II/base/timer.h>
#include <deal.II/base/utilities.h>
#include <deal.II/base/vectorization.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_system.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/la_parallel_vector.h>

#include <deal.II/matrix_free/fe_evaluation.h>
#include <deal.II/matrix_free/matrix_free.h>
#include <deal.II/matrix_free/operators.h>

#include <deal.II/numerics/data_out.h>

#include <fstream>
#include <iomanip>
#include <iostream>


#define ENABLE_MPI

#include "performance_test_driver.h"

namespace NavierStokes_DG
{
  using namespace dealii;

  constexpr unsigned int dimension     = 3;
  constexpr unsigned int fe_degree     = 4;
  constexpr unsigned int n_q_points_1d = fe_degree + 2;

  constexpr unsigned int group_size = numbers::invalid_unsigned_int;

  using Number = double;

  using VectorizedArrayType = VectorizedArray<Number>;

  constexpr double gamma     = 1.4;
  constexpr double R         = 287.;
  constexpr double c_v       = R / (gamma - 1.);
  constexpr double c_p       = gamma / c_v;
  constexpr double viscosity = 1. / 1600;
  constexpr double lambda    = viscosity * c_p / 0.71;
  constexpr double Ma        = 0.1;

  const double courant_number = 0.07 / std::pow(fe_degree, 1.5);



  class LowStorageRungeKuttaIntegrator
  {
  public:
    LowStorageRungeKuttaIntegrator()
    {
      TimeStepping::LowStorageRungeKutta<
        LinearAlgebra::distributed::Vector<Number>>
        rk_integrator(TimeStepping::LOW_STORAGE_RK_STAGE3_ORDER3);
      std::vector<double> ci; // not used
      rk_integrator.get_coefficients(ai, bi, ci);
    }

    unsigned int
    n_stages() const
    {
      return bi.size();
    }

    template <typename VectorType, typename Operator>
    void
    perform_time_step(const Operator &pde_operator,
                      const double    current_time,
                      const double    time_step,
                      VectorType     &solution,
                      VectorType     &vec_ri,
                      VectorType     &vec_ki) const
    {
      AssertDimension(ai.size() + 1, bi.size());

      vec_ki.swap(solution);

      double sum_previous_bi = 0;
      for (unsigned int stage = 0; stage < bi.size(); ++stage)
        {
          const double c_i = stage == 0 ? 0 : sum_previous_bi + ai[stage - 1];

          pde_operator.perform_stage(stage,
                                     current_time + c_i * time_step,
                                     bi[stage] * time_step,
                                     (stage == bi.size() - 1 ?
                                        0 :
                                        ai[stage] * time_step),
                                     (stage % 2 == 0 ? vec_ki : vec_ri),
                                     (stage % 2 == 0 ? vec_ri : vec_ki),
                                     solution);

          if (stage > 0)
            sum_previous_bi += bi[stage - 1];
        }
    }

  private:
    std::vector<double> bi;
    std::vector<double> ai;
  };



  template <int dim>
  class ExactSolution : public Function<dim>
  {
  public:
    ExactSolution(const double time)
      : Function<dim>(dim + 2, time)
    {}

    virtual double
    value(const Point<dim> &p, const unsigned int component = 0) const override;
  };



  template <int dim>
  double
  ExactSolution<dim>::value(const Point<dim>  &x,
                            const unsigned int component) const
  {
    const double c0 = 1. / Ma;
    const double T0 = c0 * c0 / (gamma * R);
    if (component == 0)
      return 1 + 1. / (R * T0) * 1. / 16. *
                   (std::cos(2 * x[0]) + std::cos(2 * x[1])) *
                   (std::cos(2 * x[2]) + 2.);
    else if (component == 1)
      return std::sin(x[0]) * std::cos(x[1]) * std::cos(x[2]);
    else if (component == 2)
      return -std::cos(x[0]) * std::sin(x[1]) * std::cos(x[2]);
    else if (component == dim + 1)
      return c_v * T0 +
             0.5 * (Utilities::fixed_power<2>(std::sin(x[0]) * std::cos(x[1]) *
                                              std::cos(x[2])) +
                    Utilities::fixed_power<2>(std::cos(x[0]) * std::sin(x[1]) *
                                              std::cos(x[2])));
    else
      return 0.;
  }



  template <int dim, typename Number>
  inline DEAL_II_ALWAYS_INLINE //
    Tensor<1, dim, Number>
    fluid_velocity(const Tensor<1, dim + 2, Number> &conserved_variables)
  {
    const Number inverse_density = Number(1.) / conserved_variables[0];

    Tensor<1, dim, Number> velocity;
    for (unsigned int d = 0; d < dim; ++d)
      velocity[d] = conserved_variables[1 + d] * inverse_density;

    return velocity;
  }

  template <int dim, typename Number>
  inline DEAL_II_ALWAYS_INLINE //
    Number
    fluid_pressure(const Tensor<1, dim + 2, Number> &conserved_variables)
  {
    const Tensor<1, dim, Number> velocity =
      fluid_velocity<dim>(conserved_variables);

    Number rho_u_dot_u = conserved_variables[1] * velocity[0];
    for (unsigned int d = 1; d < dim; ++d)
      rho_u_dot_u += conserved_variables[1 + d] * velocity[d];

    return (gamma - 1.) * (conserved_variables[dim + 1] - 0.5 * rho_u_dot_u);
  }

  template <int dim, typename Number>
  inline DEAL_II_ALWAYS_INLINE //
    Tensor<1, dim + 2, Tensor<1, dim, Number>>
    euler_flux(const Tensor<1, dim + 2, Number> &conserved_variables)
  {
    const Tensor<1, dim, Number> velocity =
      fluid_velocity<dim>(conserved_variables);
    const Number pressure = fluid_pressure<dim>(conserved_variables);

    Tensor<1, dim + 2, Tensor<1, dim, Number>> flux;
    for (unsigned int d = 0; d < dim; ++d)
      {
        flux[0][d] = conserved_variables[1 + d];
        for (unsigned int e = 0; e < dim; ++e)
          flux[e + 1][d] = conserved_variables[e + 1] * velocity[d];
        flux[d + 1][d] += pressure;
        flux[dim + 1][d] =
          velocity[d] * (conserved_variables[dim + 1] + pressure);
      }

    return flux;
  }

  template <int n_components, int dim, typename Number>
  inline DEAL_II_ALWAYS_INLINE //
    Tensor<1, n_components, Number>
    operator*(const Tensor<1, n_components, Tensor<1, dim, Number>> &matrix,
              const Tensor<1, dim, Number>                          &vector)
  {
    Tensor<1, n_components, Number> result;
    for (unsigned int d = 0; d < n_components; ++d)
      result[d] = matrix[d] * vector;
    return result;
  }

  template <int dim, typename Number>
  inline DEAL_II_ALWAYS_INLINE //
    Tensor<1, dim + 2, Number>
    euler_numerical_flux(const Tensor<1, dim + 2, Number> &u_m,
                         const Tensor<1, dim + 2, Number> &u_p,
                         const Tensor<1, dim, Number>     &normal)
  {
    const auto velocity_m = fluid_velocity<dim>(u_m);
    const auto velocity_p = fluid_velocity<dim>(u_p);

    const auto pressure_m = fluid_pressure<dim>(u_m);
    const auto pressure_p = fluid_pressure<dim>(u_p);

    const auto flux_m = euler_flux<dim>(u_m);
    const auto flux_p = euler_flux<dim>(u_p);

    const auto avg_velocity_normal = 0.5 * ((velocity_m + velocity_p) * normal);
    const auto avg_c               = std::sqrt(std::abs(
      0.5 * gamma * (pressure_p * (1. / u_p[0]) + pressure_m * (1. / u_m[0]))));
    const Number s_pos     = std::max(Number(), avg_velocity_normal + avg_c);
    const Number s_neg     = std::min(Number(), avg_velocity_normal - avg_c);
    const Number inverse_s = Number(1.) / (s_pos - s_neg);

    return inverse_s *
           ((s_pos * (flux_m * normal) - s_neg * (flux_p * normal)) -
            s_pos * s_neg * (u_m - u_p));
  }

  template <int dim, typename Number>
  inline DEAL_II_ALWAYS_INLINE //
    Tensor<2, dim, Number>
    fluid_velocity_gradient(
      const Tensor<1, dim + 2, Number>                 &conserved_variables,
      const Tensor<1, dim + 2, Tensor<1, dim, Number>> &gradients)
  {
    const Number inverse_density = Number(1.) / conserved_variables[0];
    const Tensor<1, dim, Number> velocity =
      fluid_velocity<dim>(conserved_variables);

    Tensor<2, dim, Number> gradient;
    for (unsigned int d = 0; d < dim; ++d)
      for (unsigned int e = 0; e < dim; ++e)
        gradient[d][e] = inverse_density *
                         (gradients[d + 1][e] - velocity[d] * gradients[0][e]);

    return gradient;
  }

  template <int dim, typename Number>
  inline DEAL_II_ALWAYS_INLINE //
    Number
    fluid_temperature(const Tensor<1, dim + 2, Number> &conserved_variables)
  {
    const Number inverse_density = Number(1.) / conserved_variables[0];
    const Number inverse_R       = 1. / R;
    return fluid_pressure(conserved_variables) * inverse_density * inverse_R;
  }

  template <int dim, typename Number>
  inline DEAL_II_ALWAYS_INLINE //
    Tensor<1, dim, Number>
    fluid_temperature_gradient(
      const Tensor<1, dim + 2, Number>                 &conserved_variables,
      const Tensor<1, dim + 2, Tensor<1, dim, Number>> &gradients)
  {
    const Number inverse_R = 1. / R;
    return (gamma - 1.) * inverse_R *
           (gradients[dim + 1] -
            fluid_velocity<dim>(conserved_variables) *
              fluid_velocity_gradient(conserved_variables, gradients));
  }

  template <int dim, typename Number>
  inline DEAL_II_ALWAYS_INLINE //
    Tensor<1, dim + 2, Tensor<1, dim, Number>>
    viscous_flux(const Tensor<1, dim + 2, Number> &conserved_variables,
                 const Tensor<1, dim + 2, Tensor<1, dim, Number>> &gradients)
  {
    const Tensor<1, dim, Number> velocity =
      fluid_velocity<dim>(conserved_variables);
    const Tensor<2, dim, Number> grad_u =
      fluid_velocity_gradient(conserved_variables, gradients);
    const Number scaled_div_u = viscosity * (2. / 3.) * trace(grad_u);

    Tensor<1, dim + 2, Tensor<1, dim, Number>> result;
    for (unsigned int d = 0; d < dim; ++d)
      {
        for (unsigned int e = d; e < dim; ++e)
          {
            result[d + 1][e] = viscosity * (grad_u[d][e] + grad_u[e][d]);
            result[e + 1][d] = result[d + 1][e];
          }
        result[d + 1][d] -= scaled_div_u;
      }

    result[dim + 1] =
      lambda * fluid_temperature_gradient(conserved_variables, gradients);
    for (unsigned int d = 0; d < dim; ++d)
      result[dim + 1][d] += result[d + 1] * velocity;

    return result;
  }



  template <int dim, typename VectorizedArrayType>
  VectorizedArrayType
  evaluate_function(const Function<dim>                   &function,
                    const Point<dim, VectorizedArrayType> &p_vectorized,
                    const unsigned int                     component)
  {
    VectorizedArrayType result;
    for (unsigned int v = 0; v < VectorizedArrayType::size(); ++v)
      {
        Point<dim> p;
        for (unsigned int d = 0; d < dim; ++d)
          p[d] = p_vectorized[d][v];
        result[v] = function.value(p, component);
      }
    return result;
  }



  template <int dim, typename VectorizedArrayType, int n_components = dim + 2>
  Tensor<1, n_components, VectorizedArrayType>
  evaluate_function(const Function<dim>                   &function,
                    const Point<dim, VectorizedArrayType> &p_vectorized)
  {
    AssertDimension(function.n_components, n_components);
    Tensor<1, n_components, VectorizedArrayType> result;
    for (unsigned int v = 0; v < VectorizedArrayType::size(); ++v)
      {
        Point<dim> p;
        for (unsigned int d = 0; d < dim; ++d)
          p[d] = p_vectorized[d][v];
        for (unsigned int d = 0; d < n_components; ++d)
          result[d][v] = function.value(p, d);
      }
    return result;
  }



  template <int dim, int degree, int n_points_1d>
  class NavierStokesOperator
  {
  public:
    static constexpr unsigned int n_quadrature_points_1d = n_points_1d;

    NavierStokesOperator();

    ~NavierStokesOperator();

    void
    reinit(const Mapping<dim> &mapping, const DoFHandler<dim> &dof_handler);

    void
    set_inflow_boundary(const types::boundary_id       boundary_id,
                        std::unique_ptr<Function<dim>> inflow_function);

    void
    set_subsonic_outflow_boundary(
      const types::boundary_id       boundary_id,
      std::unique_ptr<Function<dim>> outflow_energy);

    void
    set_wall_boundary(const types::boundary_id boundary_id);

    void
    set_body_force(std::unique_ptr<Function<dim>> body_force);

    void
    perform_stage(const unsigned int                                stage,
                  const Number                                      cur_time,
                  const Number                                      bi,
                  const Number                                      ai,
                  const LinearAlgebra::distributed::Vector<Number> &current_ri,
                  LinearAlgebra::distributed::Vector<Number>       &vec_ki,
                  LinearAlgebra::distributed::Vector<Number> &solution) const;

    void
    perform_stage_face(
      const unsigned int                                stage,
      const Number                                      cur_time,
      const Number                                      bi,
      const Number                                      ai,
      const LinearAlgebra::distributed::Vector<Number> &current_ri,
      LinearAlgebra::distributed::Vector<Number>       &vec_ki,
      LinearAlgebra::distributed::Vector<Number>       &solution) const;

    void
    project(const Function<dim>                        &function,
            LinearAlgebra::distributed::Vector<Number> &solution) const;

    std::array<double, 3>
    compute_errors(
      const Function<dim>                              &function,
      const LinearAlgebra::distributed::Vector<Number> &solution) const;

    std::array<double, 2>
    compute_kinetic_energy(
      const LinearAlgebra::distributed::Vector<Number> &solution) const;

    double
    compute_cell_transport_speed(
      const LinearAlgebra::distributed::Vector<Number> &solution) const;

    void
    initialize_vector(LinearAlgebra::distributed::Vector<Number> &vector) const;

    mutable double time_loop;
    mutable double time_rk_update;

  private:
    MPI_Comm subcommunicator;

    MatrixFree<dim, Number, VectorizedArrayType> data;

    std::map<types::boundary_id, std::unique_ptr<Function<dim>>>
      inflow_boundaries;
    std::map<types::boundary_id, std::unique_ptr<Function<dim>>>
                                   subsonic_outflow_boundaries;
    std::set<types::boundary_id>   wall_boundaries;
    std::unique_ptr<Function<dim>> body_force;

    void
    operation_on_cell(const MatrixFree<dim, Number, VectorizedArrayType> &mf,
                      LinearAlgebra::distributed::Vector<Number>         &dst,
                      const LinearAlgebra::distributed::Vector<Number>   &src,
                      const std::pair<unsigned int, unsigned int> &range) const;

    void
    operation_cell(const MatrixFree<dim, Number, VectorizedArrayType> &mf,
                   LinearAlgebra::distributed::Vector<Number>         &dst,
                   const LinearAlgebra::distributed::Vector<Number>   &src,
                   const std::pair<unsigned int, unsigned int> &range) const;

    void
    operation_face(const MatrixFree<dim, Number, VectorizedArrayType> &mf,
                   LinearAlgebra::distributed::Vector<Number>         &dst,
                   const LinearAlgebra::distributed::Vector<Number>   &src,
                   const std::pair<unsigned int, unsigned int> &range) const;

    void
    operation_boundary(
      const MatrixFree<dim, Number, VectorizedArrayType> &mf,
      LinearAlgebra::distributed::Vector<Number>         &dst,
      const LinearAlgebra::distributed::Vector<Number>   &src,
      const std::pair<unsigned int, unsigned int>        &range) const;

    void
    local_apply_inverse_mass_matrix(
      const MatrixFree<dim, Number> &,
      LinearAlgebra::distributed::Vector<Number>       &dst,
      const LinearAlgebra::distributed::Vector<Number> &src,
      const std::pair<unsigned int, unsigned int>      &cell_range) const;

    mutable double                                      ai;
    mutable double                                      bi;
    mutable LinearAlgebra::distributed::Vector<Number> *solution;
    mutable unsigned int                                stage;
  };



  template <int dim, int degree, int n_points_1d>
  NavierStokesOperator<dim, degree, n_points_1d>::NavierStokesOperator()
  {
#ifdef DEAL_II_WITH_MPI
    if (group_size == 1)
      {
        this->subcommunicator = MPI_COMM_SELF;
      }
    else if (group_size == numbers::invalid_unsigned_int)
      {
        const auto rank = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);

        MPI_Comm_split_type(MPI_COMM_WORLD,
                            MPI_COMM_TYPE_SHARED,
                            rank,
                            MPI_INFO_NULL,
                            &subcommunicator);
      }
    else
      {
        DEAL_II_NOT_IMPLEMENTED();
      }
#else
    (void)subcommunicator;
    (void)group_size;
    this->subcommunicator = MPI_COMM_SELF;
#endif

    time_loop      = 0.;
    time_rk_update = 0.;
  }


  template <int dim, int degree, int n_points_1d>
  NavierStokesOperator<dim, degree, n_points_1d>::~NavierStokesOperator()
  {
#ifdef DEAL_II_WITH_MPI
    if (this->subcommunicator != MPI_COMM_SELF)
      MPI_Comm_free(&subcommunicator);
#endif
  }


  template <int dim, int degree, int n_points_1d>
  void
  NavierStokesOperator<dim, degree, n_points_1d>::reinit(
    const Mapping<dim>    &mapping,
    const DoFHandler<dim> &dof_handler)
  {
    const std::vector<const DoFHandler<dim> *> dof_handlers = {&dof_handler};
    const AffineConstraints<double>            dummy;
    const std::vector<const AffineConstraints<double> *> constraints = {&dummy};
    const std::vector<Quadrature<1>> quadratures = {QGauss<1>(n_q_points_1d),
                                                    QGauss<1>(fe_degree + 1)};

    typename MatrixFree<dim, Number, VectorizedArrayType>::AdditionalData
      additional_data;
    additional_data.mapping_update_flags =
      (update_gradients | update_JxW_values | update_quadrature_points |
       update_values);
    additional_data.mapping_update_flags_inner_faces =
      (update_JxW_values | update_quadrature_points | update_normal_vectors |
       update_values);
    additional_data.mapping_update_flags_boundary_faces =
      (update_JxW_values | update_quadrature_points | update_normal_vectors |
       update_values);
    additional_data.tasks_parallel_scheme =
      MatrixFree<dim, Number, VectorizedArrayType>::AdditionalData::none;

    MatrixFreeTools::categorize_by_boundary_ids(dof_handler.get_triangulation(),
                                                additional_data);

    additional_data.communicator_sm = subcommunicator;

    data.reinit(
      mapping, dof_handlers, constraints, quadratures, additional_data);
  }



  template <int dim, int degree, int n_points_1d>
  void
  NavierStokesOperator<dim, degree, n_points_1d>::perform_stage(
    const unsigned int                                stage,
    const Number                                      current_time,
    const Number                                      bi,
    const Number                                      ai,
    const LinearAlgebra::distributed::Vector<Number> &current_ri,
    LinearAlgebra::distributed::Vector<Number>       &vec_ki,
    LinearAlgebra::distributed::Vector<Number>       &solution) const
  {
    for (auto &i : inflow_boundaries)
      i.second->set_time(current_time);
    for (auto &i : subsonic_outflow_boundaries)
      i.second->set_time(current_time);

    this->ai       = ai;
    this->bi       = bi;
    this->solution = &solution;
    this->stage    = stage;

    data.loop_cell_centric(
      &NavierStokesOperator::operation_on_cell,
      this,
      vec_ki,
      current_ri,
      true,
      MatrixFree<dim, Number, VectorizedArrayType>::DataAccessOnFaces::values);
  }



  template <int dim, int degree, int n_points_1d>
  void
  NavierStokesOperator<dim, degree, n_points_1d>::operation_on_cell(
    const MatrixFree<dim, Number, VectorizedArrayType> &data,
    LinearAlgebra::distributed::Vector<Number>         &vec_ki,
    const LinearAlgebra::distributed::Vector<Number>   &current_ri,
    const std::pair<unsigned int, unsigned int>        &cell_range) const
  {
    using FECellIntegral = FEEvaluation<dim,
                                        degree,
                                        n_points_1d,
                                        dim + 2,
                                        Number,
                                        VectorizedArrayType>;
    using FEFaceIntegral = FEFaceEvaluation<dim,
                                            degree,
                                            n_points_1d,
                                            dim + 2,
                                            Number,
                                            VectorizedArrayType>;

    FECellIntegral phi(data);
    FECellIntegral phi_temp(data);
    FEFaceIntegral phi_m(data, true);
    FEFaceIntegral phi_p(data, false);

    Tensor<1, dim, VectorizedArrayType>     constant_body_force;
    const Functions::ConstantFunction<dim> *constant_function =
      dynamic_cast<Functions::ConstantFunction<dim> *>(body_force.get());

    if (constant_function)
      constant_body_force = evaluate_function<dim, VectorizedArrayType, dim>(
        *constant_function, Point<dim, VectorizedArrayType>());

    const dealii::internal::EvaluatorTensorProduct<
      dealii::internal::EvaluatorVariant::evaluate_evenodd,
      dim,
      n_points_1d,
      n_points_1d,
      VectorizedArrayType,
      Number>
      eval({},
           data.get_shape_info().data[0].shape_gradients_collocation_eo,
           {});

    internal::EvaluatorTensorProduct<
      internal::EvaluatorVariant::evaluate_evenodd,
      dim - 1,
      n_q_points_1d,
      n_q_points_1d,
      VectorizedArrayType,
      Number>
      eval_face({},
                data.get_shape_info().data[0].shape_gradients_collocation_eo,
                {});

    AlignedVector<VectorizedArrayType> buffer;
    buffer.resize_fast(phi.static_n_q_points * phi.n_components);
    AlignedVector<VectorizedArrayType> buffer_face;
    buffer_face.resize_fast(phi_m.static_n_q_points * 2);

    for (unsigned int cell = cell_range.first; cell < cell_range.second; ++cell)
      {
        phi.reinit(cell);

        if (ai != Number())
          phi_temp.reinit(cell);

        if (ai != Number() && stage == 0)
          {
            phi.read_dof_values(current_ri);

            for (unsigned int i = 0;
                 i < phi.static_dofs_per_component * (dim + 2);
                 ++i)
              phi_temp.begin_dof_values()[i] = phi.begin_dof_values()[i];

            phi.evaluate(EvaluationFlags::values | EvaluationFlags::gradients);
          }
        else
          {
            phi.gather_evaluate(current_ri,
                                EvaluationFlags::values |
                                  EvaluationFlags::gradients);
          }

        for (unsigned int i = 0; i < phi.static_n_q_points * (dim + 2); ++i)
          buffer[i] = phi.begin_values()[i];

        for (const unsigned int q : phi.quadrature_point_indices())
          {
            const auto w_q      = phi.get_value(q);
            const auto grad_w_q = phi.get_gradient(q);
            auto       flux     = euler_flux<dim>(w_q);
            const auto viscous  = viscous_flux(w_q, grad_w_q);
            for (unsigned int d = 0; d < dim + 2; ++d)
              flux[d] = flux[d] - viscous[d];
            phi.submit_gradient(flux, q);
            if (body_force.get() != nullptr)
              {
                const Tensor<1, dim, VectorizedArrayType> force =
                  constant_function ?
                    constant_body_force :
                    evaluate_function<dim, VectorizedArrayType, dim>(
                      *body_force, phi.quadrature_point(q));

                Tensor<1, dim + 2, VectorizedArrayType> forcing;
                for (unsigned int d = 0; d < dim; ++d)
                  forcing[d + 1] = w_q[0] * force[d];
                for (unsigned int d = 0; d < dim; ++d)
                  forcing[dim + 1] += force[d] * w_q[d + 1];

                phi.submit_value(forcing, q);
              }
          }

        {
          auto *values_ptr   = phi.begin_values();
          auto *gradient_ptr = phi.begin_gradients();

          for (unsigned int c = 0; c < dim + 2; ++c)
            {
              if (dim >= 1 && body_force.get() == nullptr)
                eval.template gradients<0, false, false, dim>(gradient_ptr + 0,
                                                              values_ptr);
              else if (dim >= 1)
                eval.template gradients<0, false, true, dim>(gradient_ptr + 0,
                                                             values_ptr);
              if (dim >= 2)
                eval.template gradients<1, false, true, dim>(gradient_ptr + 1,
                                                             values_ptr);
              if (dim >= 3)
                eval.template gradients<2, false, true, dim>(gradient_ptr + 2,
                                                             values_ptr);

              values_ptr += phi.static_n_q_points;
              gradient_ptr += phi.static_n_q_points * dim;
            }
        }

        for (unsigned int face = 0; face < GeometryInfo<dim>::faces_per_cell;
             ++face)
          {
            const auto boundary_ids =
              data.get_faces_by_cells_boundary_id(cell, face);

            Assert(std::equal(boundary_ids.begin(),
                              boundary_ids.begin() +
                                data.n_active_entries_per_cell_batch(cell),
                              boundary_ids.begin()),
                   ExcMessage("Boundary IDs of lanes differ."));

            const auto boundary_id = boundary_ids[0];

            phi_m.reinit(cell, face);

            const AlignedVector<Number> &shape_data =
              data.get_shape_info().data.front().quadrature_data_on_face[face %
                                                                         2];
            const std::array<int, 2> n_blocks{
              {(dim > 1 ? n_q_points_1d : 1), (dim > 2 ? n_q_points_1d : 1)}};

            std::array<int, 2> steps;
            if (face / 2 == 0)
              steps = {{n_q_points_1d, 0}};
            else if (dim == 2 && face / 2 == 1)
              steps = {{1, 0}};
            else if (face / 2 == 1)
              // in 3d, the coordinate system is zx, not xz -> switch indices
              steps = {{n_q_points_1d * n_q_points_1d,
                        -static_cast<int>(n_q_points_1d * n_q_points_1d *
                                          n_q_points_1d) +
                          1}};
            else
              steps = {{1, 0}};

            for (unsigned int d = 0; d < dim + 2; ++d)
              {
                const unsigned int n_q_points_face = phi_m.static_n_q_points;
                if (face / 2 == 0)
                  internal::
                    interpolate_to_face<n_q_points_1d, 1, true, false, 1>(
                      shape_data.begin(),
                      n_blocks,
                      steps,
                      buffer.data() + d * phi.static_n_q_points,
                      buffer_face.data());
                else if (face / 2 == 1)
                  internal::interpolate_to_face<n_q_points_1d,
                                                n_q_points_1d,
                                                true,
                                                false,
                                                1>(shape_data.begin(),
                                                   n_blocks,
                                                   steps,
                                                   buffer.data() +
                                                     d * phi.static_n_q_points,
                                                   buffer_face.data());
                else if (face / 2 == 2)
                  internal::interpolate_to_face<n_q_points_1d,
                                                n_q_points_1d * n_q_points_1d,
                                                true,
                                                false,
                                                1>(shape_data.begin(),
                                                   n_blocks,
                                                   steps,
                                                   buffer.data() +
                                                     d * phi.static_n_q_points,
                                                   buffer_face.data());

                if (dim > 1)
                  eval_face.template gradients<0, true, false, dim>(
                    buffer_face.data(),
                    phi_m.begin_gradients() + (d * dim) * n_q_points_face);
                if (dim > 2)
                  eval_face.template gradients<1, true, false, dim>(
                    buffer_face.data(),
                    phi_m.begin_gradients() + (d * dim) * n_q_points_face + 1);

                for (unsigned int i = 0; i < n_q_points_face; ++i)
                  {
                    phi_m.begin_values()[d * n_q_points_face + i] =
                      buffer_face[i];

                    phi_m.begin_gradients()[(d * dim) * n_q_points_face +
                                            i * dim + dim - 1] =
                      buffer_face[n_q_points_face + i];
                  }
              }

            if (boundary_id == numbers::internal_face_boundary_id)
              {
                phi_p.reinit(cell, face);
                phi_p.gather_evaluate(current_ri,
                                      EvaluationFlags::values |
                                        EvaluationFlags::gradients);

                const auto tau_ip =
                  (std::abs((phi_m.normal_vector(0) *
                             phi_m.inverse_jacobian(0))[dim - 1]) +
                   std::abs((phi_p.normal_vector(0) *
                             phi_p.inverse_jacobian(0))[dim - 1])) *
                  Number(viscosity * (degree + 1) * (degree + 1));

                for (const unsigned int q : phi_m.quadrature_point_indices())
                  {
                    const auto w_m    = phi_m.get_value(q);
                    const auto w_p    = phi_p.get_value(q);
                    const auto normal = phi_m.normal_vector(q);
                    auto       numerical_flux =
                      -euler_numerical_flux<dim>(w_m, w_p, normal);
                    const auto grad_w_m = phi_m.get_gradient(q);
                    const auto grad_w_p = phi_p.get_gradient(q);

                    const auto flux_q1 = viscous_flux(w_m, grad_w_m);
                    for (unsigned int d = 0; d < dim + 2; ++d)
                      numerical_flux[d] += 0.5 * (flux_q1[d] * normal);
                    const auto flux_q2 = viscous_flux(w_p, grad_w_p);
                    for (unsigned int d = 0; d < dim + 2; ++d)
                      numerical_flux[d] += 0.5 * (flux_q2[d] * normal);
                    numerical_flux -= tau_ip * (w_m - w_p);
                    phi_m.submit_value(numerical_flux, q);

                    Tensor<1, dim + 2, Tensor<1, dim, VectorizedArrayType>>
                      w_jump;
                    for (unsigned int d = 0; d < dim + 2; ++d)
                      for (unsigned int e = 0; e < dim; ++e)
                        w_jump[d][e] =
                          (w_m[d] - w_p[d]) * (Number(0.5) * normal[e]);
                    phi_m.submit_gradient(viscous_flux(w_m, w_jump), q);
                  }
              }
            else
              {
                const auto tau_ip =
                  std::abs((phi_m.normal_vector(0) *
                            phi_m.inverse_jacobian(0))[dim - 1]) *
                  Number(2. * viscosity * (degree + 1) * (degree + 1));

                for (const unsigned int q : phi_m.quadrature_point_indices())
                  {
                    const auto w_m      = phi_m.get_value(q);
                    const auto normal   = phi_m.normal_vector(q);
                    const auto grad_w_m = phi_m.get_gradient(q);
                    const auto grad_w_p = grad_w_m;

                    auto rho_u_dot_n = w_m[1] * normal[0];
                    for (unsigned int d = 1; d < dim; ++d)
                      rho_u_dot_n += w_m[1 + d] * normal[d];

                    bool at_outflow = false;

                    Tensor<1, dim + 2, VectorizedArrayType> w_p;

                    if (wall_boundaries.find(boundary_id) !=
                        wall_boundaries.end())
                      {
                        w_p[0] = w_m[0];
                        for (unsigned int d = 0; d < dim; ++d)
                          w_p[d + 1] =
                            w_m[d + 1] - 2. * rho_u_dot_n * normal[d];
                        w_p[dim + 1] = w_m[dim + 1];
                      }
                    else if (inflow_boundaries.find(boundary_id) !=
                             inflow_boundaries.end())
                      w_p = evaluate_function(
                        *inflow_boundaries.find(boundary_id)->second,
                        phi_m.quadrature_point(q));
                    else if (subsonic_outflow_boundaries.find(boundary_id) !=
                             subsonic_outflow_boundaries.end())
                      {
                        w_p = w_m;
                        w_p[dim + 1] =
                          evaluate_function(*subsonic_outflow_boundaries
                                               .find(boundary_id)
                                               ->second,
                                            phi_m.quadrature_point(q),
                                            dim + 1);
                        at_outflow = true;
                      }
                    else
                      AssertThrow(false,
                                  ExcMessage(
                                    "Unknown boundary id, did "
                                    "you set a boundary condition for "
                                    "this part of the domain boundary?"));

                    auto flux = -euler_numerical_flux<dim>(w_m, w_p, normal);

                    if (at_outflow)
                      for (unsigned int v = 0; v < VectorizedArrayType::size();
                           ++v)
                        {
                          if (rho_u_dot_n[v] < -1e-12)
                            for (unsigned int d = 0; d < dim; ++d)
                              flux[d + 1][v] = 0.;
                        }

                    const auto flux_q1 = viscous_flux(w_m, grad_w_m);
                    for (unsigned int d = 0; d < dim + 2; ++d)
                      flux[d] += 0.5 * (flux_q1[d] * normal);
                    const auto flux_q2 = viscous_flux(w_p, grad_w_p);
                    for (unsigned int d = 0; d < dim + 2; ++d)
                      flux[d] += 0.5 * (flux_q2[d] * normal);
                    flux -= tau_ip * (w_m - w_p);
                    phi_m.submit_value(flux, q);

                    Tensor<1, dim + 2, Tensor<1, dim, VectorizedArrayType>>
                      w_jump;
                    for (unsigned int d = 0; d < dim + 2; ++d)
                      for (unsigned int e = 0; e < dim; ++e)
                        w_jump[d][e] =
                          (w_m[d] - w_p[d]) * (Number(0.5) * normal[e]);
                    phi_m.submit_gradient(viscous_flux(w_m, w_jump), q);
                  }
              }

            for (unsigned int d = 0; d < dim + 2; ++d)
              {
                const unsigned int n_q_points_face = phi_m.static_n_q_points;
                for (unsigned int i = 0; i < n_q_points_face; ++i)
                  {
                    buffer_face[i] =
                      phi_m.begin_values()[d * n_q_points_face + i];
                    buffer_face[n_q_points_face + i] =
                      phi_m.begin_gradients()[d * dim * n_q_points_face +
                                              i * dim + dim - 1];
                  }

                if (dim > 2)
                  eval_face.template gradients<1, false, true, dim>(
                    phi_m.begin_gradients() + d * dim * n_q_points_face + 1,
                    buffer_face.data());
                if (dim > 1)
                  eval_face.template gradients<0, false, true, dim>(
                    phi_m.begin_gradients() + d * dim * n_q_points_face,
                    buffer_face.data());

                if (face / 2 == 0)
                  internal::
                    interpolate_to_face<n_q_points_1d, 1, false, true, 1>(
                      shape_data.begin(),
                      n_blocks,
                      steps,
                      buffer_face.data(),
                      phi.begin_values() + d * phi.static_n_q_points);
                else if (face / 2 == 1)
                  internal::interpolate_to_face<n_q_points_1d,
                                                n_q_points_1d,
                                                false,
                                                true,
                                                1>(shape_data.begin(),
                                                   n_blocks,
                                                   steps,
                                                   buffer_face.data(),
                                                   phi.begin_values() +
                                                     d * phi.static_n_q_points);
                else if (face / 2 == 2)
                  internal::interpolate_to_face<n_q_points_1d,
                                                n_q_points_1d * n_q_points_1d,
                                                false,
                                                true,
                                                1>(shape_data.begin(),
                                                   n_blocks,
                                                   steps,
                                                   buffer_face.data(),
                                                   phi.begin_values() +
                                                     d * phi.static_n_q_points);
              }
          }

        for (unsigned int q = 0; q < phi.static_n_q_points; ++q)
          {
            const auto factor = VectorizedArrayType(1.0) / phi.JxW(q);
            for (unsigned int c = 0; c < dim + 2; ++c)
              phi.begin_values()[c * phi.static_n_q_points + q] =
                phi.begin_values()[c * phi.static_n_q_points + q] * factor;
          }

        internal::FEEvaluationImplBasisChange<
          dealii::internal::EvaluatorVariant::evaluate_evenodd,
          internal::EvaluatorQuantity::hessian,
          dim,
          degree + 1,
          n_points_1d>::do_backward(dim + 2,
                                    data.get_shape_info()
                                      .data[0]
                                      .inverse_shape_values_eo,
                                    false,
                                    phi.begin_values(),
                                    phi.begin_dof_values());

        if (ai == Number())
          {
            for (unsigned int q = 0; q < phi.static_dofs_per_cell; ++q)
              phi.begin_dof_values()[q] = bi * phi.begin_dof_values()[q];
            phi.distribute_local_to_global(*solution);
          }
        else
          {
            if (stage != 0)
              phi_temp.read_dof_values(*solution);

            for (unsigned int q = 0; q < phi.static_dofs_per_cell; ++q)
              {
                const auto K_i = phi.begin_dof_values()[q];

                phi.begin_dof_values()[q] =
                  phi_temp.begin_dof_values()[q] + (ai * K_i);

                phi_temp.begin_dof_values()[q] += bi * K_i;
              }
            phi.set_dof_values(vec_ki);
            phi_temp.set_dof_values(*solution);
          }
      }
  }



  template <int dim, int degree, int n_points_1d>
  void
  NavierStokesOperator<dim, degree, n_points_1d>::perform_stage_face(
    const unsigned int                                stage,
    const Number                                      current_time,
    const Number                                      bi,
    const Number                                      ai,
    const LinearAlgebra::distributed::Vector<Number> &current_ri,
    LinearAlgebra::distributed::Vector<Number>       &vec_ki,
    LinearAlgebra::distributed::Vector<Number>       &solution) const
  {
    for (auto &i : inflow_boundaries)
      i.second->set_time(current_time);
    for (auto &i : subsonic_outflow_boundaries)
      i.second->set_time(current_time);

    {
      Timer timer;
      data.loop(&NavierStokesOperator::operation_cell,
                &NavierStokesOperator::operation_face,
                &NavierStokesOperator::operation_boundary,
                this,
                vec_ki,
                current_ri,
                true,
                MatrixFree<dim, Number, VectorizedArrayType>::
                  DataAccessOnFaces::gradients,
                MatrixFree<dim, Number, VectorizedArrayType>::
                  DataAccessOnFaces::gradients);
      time_loop += timer.wall_time();
    }

    {
      Timer timer;
      data.cell_loop(
        &NavierStokesOperator::local_apply_inverse_mass_matrix,
        this,
        vec_ki,
        vec_ki,
        std::function<void(const unsigned int, const unsigned int)>(),
        [&](const unsigned int start_range, const unsigned int end_range) {
          if (ai == Number())
            {
              /* DEAL_II_OPENMP_SIMD_PRAGMA */
              for (unsigned int i = start_range; i < end_range; ++i)
                {
                  const Number k_i          = vec_ki.local_element(i);
                  const Number sol_i        = solution.local_element(i);
                  solution.local_element(i) = sol_i + bi * k_i;
                }
            }
          else
            {
              /* DEAL_II_OPENMP_SIMD_PRAGMA */
              if (stage == 0)
                for (unsigned int i = start_range; i < end_range; ++i)
                  {
                    const Number k_i          = vec_ki.local_element(i);
                    const Number sol_i        = current_ri.local_element(i);
                    solution.local_element(i) = sol_i + bi * k_i;
                    vec_ki.local_element(i)   = sol_i + ai * k_i;
                  }
              else
                for (unsigned int i = start_range; i < end_range; ++i)
                  {
                    const Number k_i          = vec_ki.local_element(i);
                    const Number sol_i        = solution.local_element(i);
                    solution.local_element(i) = sol_i + bi * k_i;
                    vec_ki.local_element(i)   = sol_i + ai * k_i;
                  }
            }
        });
      time_rk_update += timer.wall_time();
    }
  }



  template <int dim, int degree, int n_points_1d>
  void
  NavierStokesOperator<dim, degree, n_points_1d>::operation_cell(
    const MatrixFree<dim, Number, VectorizedArrayType> &data,
    LinearAlgebra::distributed::Vector<Number>         &dst,
    const LinearAlgebra::distributed::Vector<Number>   &src,
    const std::pair<unsigned int, unsigned int>        &cell_range) const
  {
    using FECellIntegral = FEEvaluation<dim,
                                        degree,
                                        n_points_1d,
                                        dim + 2,
                                        Number,
                                        VectorizedArrayType>;

    FECellIntegral phi(data);

    Tensor<1, dim, VectorizedArrayType>     constant_body_force;
    const Functions::ConstantFunction<dim> *constant_function =
      dynamic_cast<Functions::ConstantFunction<dim> *>(body_force.get());

    if (constant_function)
      constant_body_force = evaluate_function<dim, VectorizedArrayType, dim>(
        *constant_function, Point<dim, VectorizedArrayType>());

    for (unsigned int cell = cell_range.first; cell < cell_range.second; ++cell)
      {
        phi.reinit(cell);
        phi.gather_evaluate(src,
                            EvaluationFlags::values |
                              EvaluationFlags::gradients);

        for (const unsigned int q : phi.quadrature_point_indices())
          {
            const auto w_q      = phi.get_value(q);
            const auto grad_w_q = phi.get_gradient(q);
            auto       flux     = euler_flux<dim>(w_q);
            const auto viscous  = viscous_flux(w_q, grad_w_q);
            for (unsigned int d = 0; d < dim + 2; ++d)
              flux[d] = flux[d] - viscous[d];
            phi.submit_gradient(flux, q);
            if (body_force.get() != nullptr)
              {
                const Tensor<1, dim, VectorizedArrayType> force =
                  constant_function ?
                    constant_body_force :
                    evaluate_function<dim, VectorizedArrayType, dim>(
                      *body_force, phi.quadrature_point(q));

                Tensor<1, dim + 2, VectorizedArrayType> forcing;
                for (unsigned int d = 0; d < dim; ++d)
                  forcing[d + 1] = w_q[0] * force[d];
                for (unsigned int d = 0; d < dim; ++d)
                  forcing[dim + 1] += force[d] * w_q[d + 1];

                phi.submit_value(forcing, q);
              }
          }

        phi.integrate_scatter(((body_force.get() != nullptr) ?
                                 EvaluationFlags::values :
                                 EvaluationFlags::nothing) |
                                EvaluationFlags::gradients,
                              dst);
      }
  }



  template <int dim, int degree, int n_points_1d>
  void
  NavierStokesOperator<dim, degree, n_points_1d>::operation_face(
    const MatrixFree<dim, Number, VectorizedArrayType> &data,
    LinearAlgebra::distributed::Vector<Number>         &dst,
    const LinearAlgebra::distributed::Vector<Number>   &src,
    const std::pair<unsigned int, unsigned int>        &face_range) const
  {
    using FEFaceIntegral = FEFaceEvaluation<dim,
                                            degree,
                                            n_points_1d,
                                            dim + 2,
                                            Number,
                                            VectorizedArrayType>;
    FEFaceIntegral phi_m(data, true);
    FEFaceIntegral phi_p(data, false);
    for (unsigned int face = face_range.first; face < face_range.second; ++face)
      {
        phi_p.reinit(face);
        phi_p.gather_evaluate(src,
                              EvaluationFlags::values |
                                EvaluationFlags::gradients);

        phi_m.reinit(face);
        phi_m.gather_evaluate(src,
                              EvaluationFlags::values |
                                EvaluationFlags::gradients);

        const auto tau_ip = (std::abs((phi_m.normal_vector(0) *
                                       phi_m.inverse_jacobian(0))[dim - 1]) +
                             std::abs((phi_p.normal_vector(0) *
                                       phi_p.inverse_jacobian(0))[dim - 1])) *
                            Number(viscosity * (degree + 1) * (degree + 1));

        for (const unsigned int q : phi_m.quadrature_point_indices())
          {
            const auto w_m      = phi_m.get_value(q);
            const auto w_p      = phi_p.get_value(q);
            const auto normal   = phi_m.normal_vector(q);
            auto numerical_flux = -euler_numerical_flux<dim>(w_m, w_p, normal);
            const auto grad_w_m = phi_m.get_gradient(q);
            const auto grad_w_p = phi_p.get_gradient(q);

            const auto flux_q1 = viscous_flux(w_m, grad_w_m);
            for (unsigned int d = 0; d < dim + 2; ++d)
              numerical_flux[d] += 0.5 * (flux_q1[d] * normal);
            const auto flux_q2 = viscous_flux(w_p, grad_w_p);
            for (unsigned int d = 0; d < dim + 2; ++d)
              numerical_flux[d] += 0.5 * (flux_q2[d] * normal);
            numerical_flux -= tau_ip * (w_m - w_p);
            phi_m.submit_value(numerical_flux, q);
            phi_p.submit_value(-numerical_flux, q);

            Tensor<1, dim + 2, Tensor<1, dim, VectorizedArrayType>> w_jump;
            for (unsigned int d = 0; d < dim + 2; ++d)
              for (unsigned int e = 0; e < dim; ++e)
                w_jump[d][e] = (w_m[d] - w_p[d]) * (Number(0.5) * normal[e]);
            phi_m.submit_gradient(viscous_flux(w_m, w_jump), q);
            phi_p.submit_gradient(viscous_flux(w_p, w_jump), q);
          }

        phi_m.integrate_scatter(EvaluationFlags::values |
                                  EvaluationFlags::gradients,
                                dst);
        phi_p.integrate_scatter(EvaluationFlags::values |
                                  EvaluationFlags::gradients,
                                dst);
      }
  }



  template <int dim, int degree, int n_points_1d>
  void
  NavierStokesOperator<dim, degree, n_points_1d>::operation_boundary(
    const MatrixFree<dim, Number, VectorizedArrayType> &data,
    LinearAlgebra::distributed::Vector<Number>         &dst,
    const LinearAlgebra::distributed::Vector<Number>   &src,
    const std::pair<unsigned int, unsigned int>        &face_range) const
  {
    AssertThrow(false, ExcNotImplemented());
    FEFaceEvaluation<dim,
                     degree,
                     n_points_1d,
                     dim + 2,
                     Number,
                     VectorizedArrayType>
      phi_m(data, true);
    for (unsigned int face = face_range.first; face < face_range.second; ++face)
      {
        phi_m.reinit(face);
        phi_m.gather_evaluate(src,
                              EvaluationFlags::values |
                                EvaluationFlags::gradients);

        const auto tau_ip =
          std::abs(
            (phi_m.normal_vector(0) * phi_m.inverse_jacobian(0))[dim - 1]) *
          Number(2. * viscosity * (degree + 1) * (degree + 1));

        const auto boundary_id = data.get_boundary_id(face);

        for (const unsigned int q : phi_m.quadrature_point_indices())
          {
            const auto w_m      = phi_m.get_value(q);
            const auto normal   = phi_m.normal_vector(q);
            const auto grad_w_m = phi_m.get_gradient(q);
            const auto grad_w_p = grad_w_m;

            auto rho_u_dot_n = w_m[1] * normal[0];
            for (unsigned int d = 1; d < dim; ++d)
              rho_u_dot_n += w_m[1 + d] * normal[d];

            bool at_outflow = false;

            Tensor<1, dim + 2, VectorizedArrayType> w_p;

            if (wall_boundaries.find(boundary_id) != wall_boundaries.end())
              {
                w_p[0] = w_m[0];
                for (unsigned int d = 0; d < dim; ++d)
                  w_p[d + 1] = w_m[d + 1] - 2. * rho_u_dot_n * normal[d];
                w_p[dim + 1] = w_m[dim + 1];
              }
            else if (inflow_boundaries.find(boundary_id) !=
                     inflow_boundaries.end())
              w_p =
                evaluate_function(*inflow_boundaries.find(boundary_id)->second,
                                  phi_m.quadrature_point(q));
            else if (subsonic_outflow_boundaries.find(boundary_id) !=
                     subsonic_outflow_boundaries.end())
              {
                w_p          = w_m;
                w_p[dim + 1] = evaluate_function(
                  *subsonic_outflow_boundaries.find(boundary_id)->second,
                  phi_m.quadrature_point(q),
                  dim + 1);
                at_outflow = true;
              }
            else
              AssertThrow(false,
                          ExcMessage("Unknown boundary id, did "
                                     "you set a boundary condition for "
                                     "this part of the domain boundary?"));

            auto flux = -euler_numerical_flux<dim>(w_m, w_p, normal);

            if (at_outflow)
              for (unsigned int v = 0; v < VectorizedArrayType::size(); ++v)
                {
                  if (rho_u_dot_n[v] < -1e-12)
                    for (unsigned int d = 0; d < dim; ++d)
                      flux[d + 1][v] = 0.;
                }

            const auto flux_q1 = viscous_flux(w_m, grad_w_m);
            for (unsigned int d = 0; d < dim + 2; ++d)
              flux[d] += 0.5 * (flux_q1[d] * normal);
            const auto flux_q2 = viscous_flux(w_p, grad_w_p);
            for (unsigned int d = 0; d < dim + 2; ++d)
              flux[d] += 0.5 * (flux_q2[d] * normal);
            flux -= tau_ip * (w_m - w_p);
            phi_m.submit_value(flux, q);

            Tensor<1, dim + 2, Tensor<1, dim, VectorizedArrayType>> w_jump;
            for (unsigned int d = 0; d < dim + 2; ++d)
              for (unsigned int e = 0; e < dim; ++e)
                w_jump[d][e] = (w_m[d] - w_p[d]) * (Number(0.5) * normal[e]);

            phi_m.submit_gradient(viscous_flux(w_m, w_jump), q);
          }

        phi_m.integrate_scatter(EvaluationFlags::values |
                                  EvaluationFlags::gradients,
                                dst);
      }
  }



  template <int dim, int degree, int n_points_1d>
  void
  NavierStokesOperator<dim, degree, n_points_1d>::
    local_apply_inverse_mass_matrix(
      const MatrixFree<dim, Number> &,
      LinearAlgebra::distributed::Vector<Number>       &dst,
      const LinearAlgebra::distributed::Vector<Number> &src,
      const std::pair<unsigned int, unsigned int>      &cell_range) const
  {
    FEEvaluation<dim, degree, degree + 1, dim + 2, Number> phi(data, 0, 1);
    MatrixFreeOperators::CellwiseInverseMassMatrix<dim, degree, dim + 2, Number>
      inverse(phi);

    for (unsigned int cell = cell_range.first; cell < cell_range.second; ++cell)
      {
        phi.reinit(cell);
        phi.read_dof_values(src);

        inverse.apply(phi.begin_dof_values(), phi.begin_dof_values());

        phi.set_dof_values(dst);
      }
  }



  template <int dim, int degree, int n_points_1d>
  void
  NavierStokesOperator<dim, degree, n_points_1d>::initialize_vector(
    LinearAlgebra::distributed::Vector<Number> &vector) const
  {
    data.initialize_dof_vector(vector);
  }



  template <int dim, int degree, int n_points_1d>
  void
  NavierStokesOperator<dim, degree, n_points_1d>::set_inflow_boundary(
    const types::boundary_id       boundary_id,
    std::unique_ptr<Function<dim>> inflow_function)
  {
    AssertThrow(subsonic_outflow_boundaries.find(boundary_id) ==
                    subsonic_outflow_boundaries.end() &&
                  wall_boundaries.find(boundary_id) == wall_boundaries.end(),
                ExcMessage("You already set the boundary with id " +
                           std::to_string(static_cast<int>(boundary_id)) +
                           " to another type of boundary before now setting " +
                           "it as inflow"));
    AssertThrow(inflow_function->n_components == dim + 2,
                ExcMessage("Expected function with dim+2 components"));

    inflow_boundaries[boundary_id] = std::move(inflow_function);
  }



  template <int dim, int degree, int n_points_1d>
  void
  NavierStokesOperator<dim, degree, n_points_1d>::set_subsonic_outflow_boundary(
    const types::boundary_id       boundary_id,
    std::unique_ptr<Function<dim>> outflow_function)
  {
    AssertThrow(inflow_boundaries.find(boundary_id) ==
                    inflow_boundaries.end() &&
                  wall_boundaries.find(boundary_id) == wall_boundaries.end(),
                ExcMessage("You already set the boundary with id " +
                           std::to_string(static_cast<int>(boundary_id)) +
                           " to another type of boundary before now setting " +
                           "it as subsonic outflow"));
    AssertThrow(outflow_function->n_components == dim + 2,
                ExcMessage("Expected function with dim+2 components"));

    subsonic_outflow_boundaries[boundary_id] = std::move(outflow_function);
  }



  template <int dim, int degree, int n_points_1d>
  void
  NavierStokesOperator<dim, degree, n_points_1d>::set_wall_boundary(
    const types::boundary_id boundary_id)
  {
    AssertThrow(inflow_boundaries.find(boundary_id) ==
                    inflow_boundaries.end() &&
                  subsonic_outflow_boundaries.find(boundary_id) ==
                    subsonic_outflow_boundaries.end(),
                ExcMessage("You already set the boundary with id " +
                           std::to_string(static_cast<int>(boundary_id)) +
                           " to another type of boundary before now setting " +
                           "it as wall boundary"));

    wall_boundaries.insert(boundary_id);
  }



  template <int dim, int degree, int n_points_1d>
  void
  NavierStokesOperator<dim, degree, n_points_1d>::set_body_force(
    std::unique_ptr<Function<dim>> body_force)
  {
    AssertDimension(body_force->n_components, dim);

    this->body_force = std::move(body_force);
  }



  template <int dim, int degree, int n_points_1d>
  void
  NavierStokesOperator<dim, degree, n_points_1d>::project(
    const Function<dim>                        &function,
    LinearAlgebra::distributed::Vector<Number> &solution) const
  {
    FEEvaluation<dim, degree, degree + 1, dim + 2, Number, VectorizedArrayType>
      phi(data, 0, 1);
    MatrixFreeOperators::CellwiseInverseMassMatrix<dim,
                                                   degree,
                                                   dim + 2,
                                                   Number,
                                                   VectorizedArrayType>
      inverse(phi);
    solution.zero_out_ghost_values();
    for (unsigned int cell = 0; cell < data.n_cell_batches(); ++cell)
      {
        phi.reinit(cell);
        for (const unsigned int q : phi.quadrature_point_indices())
          phi.submit_dof_value(evaluate_function(function,
                                                 phi.quadrature_point(q)),
                               q);
        inverse.transform_from_q_points_to_basis(dim + 2,
                                                 phi.begin_dof_values(),
                                                 phi.begin_dof_values());
        phi.set_dof_values(solution);
      }
  }



  template <int dim, int degree, int n_points_1d>
  std::array<double, 3>
  NavierStokesOperator<dim, degree, n_points_1d>::compute_errors(
    const Function<dim>                              &function,
    const LinearAlgebra::distributed::Vector<Number> &solution) const
  {
    double errors_squared[3] = {};
    FEEvaluation<dim, degree, n_points_1d, dim + 2, Number, VectorizedArrayType>
      phi(data, 0, 0);

    for (unsigned int cell = 0; cell < data.n_cell_batches(); ++cell)
      {
        phi.reinit(cell);
        phi.gather_evaluate(solution, EvaluationFlags::values);
        VectorizedArrayType local_errors_squared[3] = {};
        for (const unsigned int q : phi.quadrature_point_indices())
          {
            const auto error =
              evaluate_function(function, phi.quadrature_point(q)) -
              phi.get_value(q);
            const auto JxW = phi.JxW(q);

            local_errors_squared[0] += error[0] * error[0] * JxW;
            for (unsigned int d = 0; d < dim; ++d)
              local_errors_squared[1] += (error[d + 1] * error[d + 1]) * JxW;
            local_errors_squared[2] += (error[dim + 1] * error[dim + 1]) * JxW;
          }
        for (unsigned int v = 0; v < data.n_active_entries_per_cell_batch(cell);
             ++v)
          for (unsigned int d = 0; d < 3; ++d)
            errors_squared[d] += local_errors_squared[d][v];
      }

    Utilities::MPI::sum(errors_squared, MPI_COMM_WORLD, errors_squared);

    std::array<double, 3> errors;
    for (unsigned int d = 0; d < 3; ++d)
      errors[d] = std::sqrt(errors_squared[d]);

    return errors;
  }



  template <int dim, int degree, int n_points_1d>
  std::array<double, 2>
  NavierStokesOperator<dim, degree, n_points_1d>::compute_kinetic_energy(
    const LinearAlgebra::distributed::Vector<Number> &solution) const
  {
    double squared[2] = {};
    FEEvaluation<dim, degree, n_points_1d, dim + 2, Number, VectorizedArrayType>
      phi(data, 0, 0);

    for (unsigned int cell = 0; cell < data.n_cell_batches(); ++cell)
      {
        phi.reinit(cell);
        phi.gather_evaluate(solution,
                            EvaluationFlags::values |
                              EvaluationFlags::gradients);
        VectorizedArrayType local_squared[2] = {};
        for (const unsigned int q : phi.quadrature_point_indices())
          {
            const auto JxW      = phi.JxW(q);
            const auto w_q      = phi.get_value(q);
            const auto velocity = fluid_velocity<dim>(w_q);
            const auto velocity_grad =
              fluid_velocity_gradient(w_q, phi.get_gradient(q));
            local_squared[0] += velocity.norm_square() * JxW;
            local_squared[1] +=
              scalar_product(velocity_grad, velocity_grad) * JxW;
          }
        for (unsigned int v = 0; v < data.n_active_entries_per_cell_batch(cell);
             ++v)
          for (unsigned int d = 0; d < 2; ++d)
            squared[d] += local_squared[d][v];
      }

    Utilities::MPI::sum(squared, MPI_COMM_WORLD, squared);

    std::array<double, 2> result{
      {0.5 * squared[0] / Utilities::fixed_power<dim>(2. * numbers::PI),
       viscosity * squared[1] / Utilities::fixed_power<dim>(2. * numbers::PI)}};

    return result;
  }



  template <int dim, int degree, int n_points_1d>
  double
  NavierStokesOperator<dim, degree, n_points_1d>::compute_cell_transport_speed(
    const LinearAlgebra::distributed::Vector<Number> &solution) const
  {
    Number max_transport = 0;
    FEEvaluation<dim, degree, degree + 1, dim + 2, Number, VectorizedArrayType>
      phi(data, 0, 1);

    for (unsigned int cell = 0; cell < data.n_cell_batches(); ++cell)
      {
        phi.reinit(cell);
        phi.gather_evaluate(solution, EvaluationFlags::values);
        VectorizedArrayType local_max = 0.;
        for (const unsigned int q : phi.quadrature_point_indices())
          {
            const auto solution = phi.get_value(q);
            const auto velocity = fluid_velocity<dim>(solution);
            const auto pressure = fluid_pressure<dim>(solution);

            const auto          inverse_jacobian = phi.inverse_jacobian(q);
            const auto          convective_speed = inverse_jacobian * velocity;
            VectorizedArrayType convective_limit = 0.;
            for (unsigned int d = 0; d < dim; ++d)
              convective_limit =
                std::max(convective_limit, std::abs(convective_speed[d]));

            const auto speed_of_sound =
              std::sqrt(gamma * pressure * (1. / solution[0]));

            Tensor<1, dim, VectorizedArrayType> eigenvector;
            for (unsigned int d = 0; d < dim; ++d)
              eigenvector[d] = 1.;
            for (unsigned int i = 0; i < 5; ++i)
              {
                eigenvector = transpose(inverse_jacobian) *
                              (inverse_jacobian * eigenvector);
                VectorizedArrayType eigenvector_norm = 0.;
                for (unsigned int d = 0; d < dim; ++d)
                  eigenvector_norm =
                    std::max(eigenvector_norm, std::abs(eigenvector[d]));
                eigenvector /= eigenvector_norm;
              }
            const auto jac_times_ev   = inverse_jacobian * eigenvector;
            const auto max_eigenvalue = std::sqrt(
              (jac_times_ev * jac_times_ev) / (eigenvector * eigenvector));
            local_max =
              std::max(local_max,
                       max_eigenvalue * speed_of_sound + convective_limit);
          }

        for (unsigned int v = 0; v < data.n_active_entries_per_cell_batch(cell);
             ++v)
          for (unsigned int d = 0; d < 3; ++d)
            max_transport = std::max(max_transport, local_max[v]);
      }

    max_transport = Utilities::MPI::max(max_transport, MPI_COMM_WORLD);

    return max_transport;
  }



  template <int dim, int degree, int n_points_1d>
  class NavierStokesOperatorFaceCentric
  {
  public:
    NavierStokesOperatorFaceCentric(
      const NavierStokesOperator<dim, degree, n_points_1d> &ns_operator)
      : ns_operator(ns_operator)
    {}

    void
    perform_stage(const unsigned int stage,
                  const Number       current_time,
                  const Number       bi,
                  const Number       ai,
                  const LinearAlgebra::distributed::Vector<Number> &current_ri,
                  LinearAlgebra::distributed::Vector<Number>       &vec_ki,
                  LinearAlgebra::distributed::Vector<Number> &solution) const
    {
      ns_operator.perform_stage_face(
        stage, current_time, bi, ai, current_ri, vec_ki, solution);
    }

  private:
    const NavierStokesOperator<dim, degree, n_points_1d> &ns_operator;
  };



  template <int dim>
  class FlowProblem
  {
  public:
    FlowProblem();

    Measurement
    run();

  private:
    void
    make_grid();

    void
    output_results(const unsigned int result_number);

    LinearAlgebra::distributed::Vector<Number> solution;

#ifdef DEAL_II_WITH_P4EST
    parallel::distributed::Triangulation<dim> triangulation;
#else
    Triangulation<dim> triangulation;
#endif

    FESystem<dim>   fe;
    MappingQ<dim>   mapping;
    DoFHandler<dim> dof_handler;

    NavierStokesOperator<dim, fe_degree, n_q_points_1d> flow_operator;

    double time, time_step;

    class Postprocessor : public DataPostprocessor<dim>
    {
    public:
      Postprocessor();

      virtual void
      evaluate_vector_field(
        const DataPostprocessorInputs::Vector<dim> &inputs,
        std::vector<Vector<double>> &computed_quantities) const override;

      virtual std::vector<std::string>
      get_names() const override;

      virtual std::vector<
        DataComponentInterpretation::DataComponentInterpretation>
      get_data_component_interpretation() const override;

      virtual UpdateFlags
      get_needed_update_flags() const override;

    private:
      const bool do_schlieren_plot;
    };
  };



  template <int dim>
  FlowProblem<dim>::Postprocessor::Postprocessor()
    : do_schlieren_plot(dim == 2)
  {}



  template <int dim>
  void
  FlowProblem<dim>::Postprocessor::evaluate_vector_field(
    const DataPostprocessorInputs::Vector<dim> &inputs,
    std::vector<Vector<double>>                &computed_quantities) const
  {
    const unsigned int n_evaluation_points = inputs.solution_values.size();

    if (do_schlieren_plot == true)
      Assert(inputs.solution_gradients.size() == n_evaluation_points,
             ExcInternalError());

    Assert(computed_quantities.size() == n_evaluation_points,
           ExcInternalError());
    Assert(inputs.solution_values[0].size() == dim + 2, ExcInternalError());
    Assert(computed_quantities[0].size() ==
             dim + 2 + (do_schlieren_plot == true ? 1 : 0),
           ExcInternalError());

    for (unsigned int q = 0; q < n_evaluation_points; ++q)
      {
        Tensor<1, dim + 2> solution;
        for (unsigned int d = 0; d < dim + 2; ++d)
          solution[d] = inputs.solution_values[q](d);

        const double         density  = solution[0];
        const Tensor<1, dim> velocity = fluid_velocity<dim>(solution);
        const double         pressure = fluid_pressure<dim>(solution);

        for (unsigned int d = 0; d < dim; ++d)
          computed_quantities[q](d) = velocity[d];
        computed_quantities[q](dim)     = pressure;
        computed_quantities[q](dim + 1) = std::sqrt(gamma * pressure / density);

        if (do_schlieren_plot == true)
          computed_quantities[q](dim + 2) =
            inputs.solution_gradients[q][0] * inputs.solution_gradients[q][0];
      }
  }



  template <int dim>
  std::vector<std::string>
  FlowProblem<dim>::Postprocessor::get_names() const
  {
    std::vector<std::string> names;
    for (unsigned int d = 0; d < dim; ++d)
      names.emplace_back("velocity");
    names.emplace_back("pressure");
    names.emplace_back("speed_of_sound");

    if (do_schlieren_plot == true)
      names.emplace_back("schlieren_plot");

    return names;
  }



  template <int dim>
  std::vector<DataComponentInterpretation::DataComponentInterpretation>
  FlowProblem<dim>::Postprocessor::get_data_component_interpretation() const
  {
    std::vector<DataComponentInterpretation::DataComponentInterpretation>
      interpretation;
    for (unsigned int d = 0; d < dim; ++d)
      interpretation.push_back(
        DataComponentInterpretation::component_is_part_of_vector);
    interpretation.push_back(DataComponentInterpretation::component_is_scalar);
    interpretation.push_back(DataComponentInterpretation::component_is_scalar);

    if (do_schlieren_plot == true)
      interpretation.push_back(
        DataComponentInterpretation::component_is_scalar);

    return interpretation;
  }



  template <int dim>
  UpdateFlags
  FlowProblem<dim>::Postprocessor::get_needed_update_flags() const
  {
    if (do_schlieren_plot == true)
      return update_values | update_gradients;
    else
      return update_values;
  }



  template <int dim>
  FlowProblem<dim>::FlowProblem()
    :
#ifdef DEAL_II_WITH_P4EST
    triangulation(MPI_COMM_WORLD)
#endif
    , fe(FE_DGQHermite<dim>(fe_degree), dim + 2)
    , mapping(fe_degree)
    , dof_handler(triangulation)
    , flow_operator()
    , time(0)
    , time_step(0)
  {}



  template <int dim>
  void
  FlowProblem<dim>::make_grid()
  {
    Point<dim> lower_left, upper_right;
    for (unsigned int d = 0; d < dim; ++d)
      lower_left[d] = -numbers::PI;

    for (unsigned int d = 0; d < dim; ++d)
      upper_right[d] = numbers::PI;

    GridGenerator::hyper_rectangle(triangulation, lower_left, upper_right);
    for (const auto &cell : triangulation.cell_iterators())
      for (const unsigned int face : cell->face_indices())
        if (cell->at_boundary(face))
          cell->face(face)->set_boundary_id(face);
    std::vector<
      GridTools::PeriodicFacePair<typename Triangulation<dim>::cell_iterator>>
      periodic_faces;
    for (unsigned int d = 0; d < dim; ++d)
      GridTools::collect_periodic_faces(
        triangulation, 2 * d, 2 * d + 1, d, periodic_faces);
    triangulation.add_periodicity(periodic_faces);

    triangulation.refine_global(2);

    switch (get_testing_environment())
      {
        case TestingEnvironment::light:
          triangulation.refine_global(1);
          break;
        case TestingEnvironment::medium:
          triangulation.refine_global(2);
          break;
        case TestingEnvironment::heavy:
          triangulation.refine_global(3);
          break;
      }
  }



  template <int dim>
  void
  FlowProblem<dim>::output_results(const unsigned int)
  {
    Postprocessor postprocessor;
    DataOut<dim>  data_out;

    DataOutBase::VtkFlags flags;
    flags.write_higher_order_cells = true;
    data_out.set_flags(flags);

    data_out.attach_dof_handler(dof_handler);
    {
      std::vector<std::string> names;
      names.emplace_back("density");
      for (unsigned int d = 0; d < dim; ++d)
        names.emplace_back("momentum");
      names.emplace_back("energy");

      std::vector<DataComponentInterpretation::DataComponentInterpretation>
        interpretation;
      interpretation.push_back(
        DataComponentInterpretation::component_is_scalar);
      for (unsigned int d = 0; d < dim; ++d)
        interpretation.push_back(
          DataComponentInterpretation::component_is_part_of_vector);
      interpretation.push_back(
        DataComponentInterpretation::component_is_scalar);

      data_out.add_data_vector(dof_handler, solution, names, interpretation);
    }
    data_out.add_data_vector(solution, postprocessor);

    Vector<double> mpi_owner(triangulation.n_active_cells());
    mpi_owner = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
    data_out.add_data_vector(mpi_owner, "owner");

    data_out.build_patches(mapping,
                           fe.degree,
                           DataOut<dim>::curved_inner_cells);

    // Do not write a file here to be independent of file system
  }



  template <int dim>
  Measurement
  FlowProblem<dim>::run()
  {
    std::map<std::string, dealii::Timer> timer;

    timer["setup_grid"].start();
    make_grid();
    timer["setup_grid"].stop();

    timer["setup_dofs"].start();
    dof_handler.distribute_dofs(fe);
    timer["setup_dofs"].stop();

    timer["setup_matrix_free"].start();
    flow_operator.reinit(mapping, dof_handler);
    flow_operator.initialize_vector(solution);
    LinearAlgebra::distributed::Vector<Number> rk_register_1;
    LinearAlgebra::distributed::Vector<Number> rk_register_2;
    rk_register_1.reinit(solution);
    rk_register_2.reinit(solution);
    timer["setup_matrix_free"].stop();

    const LowStorageRungeKuttaIntegrator integrator;

    timer["project_initial"].start();
    flow_operator.project(ExactSolution<dim>(time), solution);

    double min_vertex_distance = std::numeric_limits<double>::max();
    for (const auto &cell : triangulation.active_cell_iterators())
      if (cell->is_locally_owned())
        min_vertex_distance =
          std::min(min_vertex_distance, cell->minimum_vertex_distance());
    min_vertex_distance =
      Utilities::MPI::min(min_vertex_distance, MPI_COMM_WORLD);

    time_step = courant_number * integrator.n_stages() /
                flow_operator.compute_cell_transport_speed(solution);
    timer["project_initial"].stop();

    time = 0;

    timer["write_output"].start();
    output_results(0);
    timer["write_output"].stop();

    unsigned int timestep_number = 0;
    while (timestep_number < 20)
      {
        timer["rk_timestep_cellbased"].start();
        integrator.perform_time_step(flow_operator,
                                     time,
                                     time_step,
                                     solution,
                                     rk_register_1,
                                     rk_register_2);
        timer["rk_timestep_cellbased"].stop();

        timer["analyze_solution"].start();
        const std::array<double, 2> energy =
          flow_operator.compute_kinetic_energy(solution);
        AssertThrow(energy[0] > 0 && energy[1] > 0, ExcInternalError());
        timer["analyze_solution"].stop();

        time += time_step;
        ++timestep_number;
      }

    while (timestep_number < 40)
      {
        timer["rk_timestep_facebased"].start();
        NavierStokesOperatorFaceCentric<dim, fe_degree, n_q_points_1d>
          flow_operator_face(flow_operator);
        integrator.perform_time_step(flow_operator_face,
                                     time,
                                     time_step,
                                     solution,
                                     rk_register_1,
                                     rk_register_2);
        timer["rk_timestep_facebased"].stop();

        timer["analyze_solution"].start();
        const std::array<double, 2> energy =
          flow_operator.compute_kinetic_energy(solution);
        AssertThrow(energy[0] > 0 && energy[1] > 0, ExcInternalError());
        timer["analyze_solution"].stop();

        time += time_step;
        ++timestep_number;
      }

    return {timer["setup_grid"].wall_time(),
            timer["setup_dofs"].wall_time(),
            timer["setup_matrix_free"].wall_time(),
            timer["project_initial"].wall_time(),
            timer["write_output"].wall_time(),
            timer["analyze_solution"].wall_time(),
            timer["rk_timestep_cellbased"].wall_time(),
            timer["rk_timestep_facebased"].wall_time(),
            flow_operator.time_loop,
            flow_operator.time_rk_update};
  }

} // namespace NavierStokes_DG



std::tuple<Metric, unsigned int, std::vector<std::string>>
describe_measurements()
{
  return {Metric::timing,
          4,
          {"setup_grid",
           "setup_dofs",
           "setup_matrix_free",
           "project_initial",
           "write_output",
           "analyze_solution",
           "rk_timestep_cellbased",
           "rk_timestep_facebased",
           "rk_timestep_facebased_loop",
           "rk_timestep_facebased_update"}};
}



Measurement
perform_single_measurement()
{
  return NavierStokes_DG::FlowProblem<NavierStokes_DG::dimension>().run();
}
