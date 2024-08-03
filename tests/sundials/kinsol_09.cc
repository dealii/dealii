// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2021 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#include <deal.II/base/parameter_handler.h>

#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/vector.h>

#include <deal.II/sundials/kinsol.h>

#include "../tests.h"

// A test for a 1d nonlinear problem taken from work of Bob Myhill in
// ASPECT. It shows that we are calling the residual() function more
// than perhaps necessary.


const double pressure    = 59448242.437;
const double temperature = 327.2685405;
const double log_edot_ii = -40.053535387;
const double grain_size  = 0.001;

const double gas_constant = 8.31446;

using namespace dealii;

namespace DiffusionCreepParameters
{
  const double prefactor           = 1.5e-16;
  const double grain_size_exponent = 4;
  const double stress_exponent     = 1;
  const double activation_energy   = 375000;
  const double activation_volume   = 6e-06;
}; // namespace DiffusionCreepParameters

namespace DislocationCreepParameters
{
  const double prefactor         = 2e-15;
  const double stress_exponent   = 3.5;
  const double activation_energy = 480000;
  const double activation_volume = 8e-06;
}; // namespace DislocationCreepParameters


std::pair<double, double>
compute_diffusion_log_strain_rate_and_derivative(const double log_stress,
                                                 const double pressure,
                                                 const double temperature)
{
  const double log_strain_rate_diffusion =
    std::log(DiffusionCreepParameters::prefactor) + log_stress -
    DiffusionCreepParameters::grain_size_exponent * std::log(grain_size) -
    (DiffusionCreepParameters::activation_energy +
     pressure * DiffusionCreepParameters::activation_volume) /
      (gas_constant * temperature);

  const double dlog_strain_rate_dlog_stress_diffusion = 1.0;

  return std::make_pair(log_strain_rate_diffusion,
                        dlog_strain_rate_dlog_stress_diffusion);
}



std::pair<double, double>
compute_dislocation_log_strain_rate_and_derivative(const double log_stress,
                                                   const double pressure,
                                                   const double temperature)
{
  const double log_strain_rate_dislocation =
    std::log(DislocationCreepParameters::prefactor) +
    DislocationCreepParameters::stress_exponent * log_stress -
    (DislocationCreepParameters::activation_energy +
     pressure * DislocationCreepParameters::activation_volume) /
      (gas_constant * temperature);

  const double dlog_strain_rate_dlog_stress_dislocation =
    DislocationCreepParameters::stress_exponent;

  return std::make_pair(log_strain_rate_dislocation,
                        dlog_strain_rate_dlog_stress_dislocation);
}



std::pair<double, double>
compute_log_strain_rate_residual_and_derivative(
  const double current_log_stress_ii,
  const double pressure,
  const double temperature,
  const double log_edot_ii)
{
  const std::pair<double, double> log_diff_edot_and_deriv =
    compute_diffusion_log_strain_rate_and_derivative(current_log_stress_ii,
                                                     pressure,
                                                     temperature);
  const std::pair<double, double> log_disl_edot_and_deriv =
    compute_dislocation_log_strain_rate_and_derivative(current_log_stress_ii,
                                                       pressure,
                                                       temperature);

  const double strain_rate_diffusion = std::exp(log_diff_edot_and_deriv.first);
  const double strain_rate_dislocation =
    std::exp(log_disl_edot_and_deriv.first);
  double log_strain_rate_deriv =
    (strain_rate_diffusion * log_diff_edot_and_deriv.second +
     strain_rate_dislocation * log_disl_edot_and_deriv.second) /
    (strain_rate_diffusion + strain_rate_dislocation);
  const double log_strain_rate_iterate =
    std::log(strain_rate_diffusion + strain_rate_dislocation);
  return std::make_pair(log_strain_rate_iterate - log_edot_ii,
                        log_strain_rate_deriv);
}


int
main()
{
  initlog();

  const double maximum_viscosity = 1.e30;
  double       log_strain_rate_deriv;

  // For diffusion creep, viscosity is grain size dependent
  const double prefactor_stress_diffusion =
    DiffusionCreepParameters::prefactor *
    std::pow(grain_size, -DiffusionCreepParameters::grain_size_exponent) *
    std::exp(
      -(std::max(DiffusionCreepParameters::activation_energy +
                   pressure * DiffusionCreepParameters::activation_volume,
                 0.0)) /
      (gas_constant * temperature));

  SUNDIALS::KINSOL<Vector<double>>::AdditionalData additional_data;
  additional_data.strategy = dealii::SUNDIALS::KINSOL<>::AdditionalData::newton;
  additional_data.function_tolerance            = 1e-10;
  additional_data.maximum_non_linear_iterations = 200;
  additional_data.maximum_setup_calls           = 10;

  int                              n_residual_evaluations = 0;
  SUNDIALS::KINSOL<Vector<double>> nonlinear_solver(additional_data);


  nonlinear_solver.reinit_vector = [&](Vector<double> &x) { x.reinit(1); };


  nonlinear_solver.residual = [&](const Vector<double> &current_log_stress_ii,
                                  Vector<double>       &residual) {
    std::tie(residual(0), log_strain_rate_deriv) =
      compute_log_strain_rate_residual_and_derivative(current_log_stress_ii[0],
                                                      pressure,
                                                      temperature,
                                                      log_edot_ii);

    deallog << std::setprecision(11)
            << "     Computing residual at x=" << current_log_stress_ii[0]
            << ", f(x)=" << residual(0) << std::endl;
    n_residual_evaluations += 1;
  };


  nonlinear_solver.setup_jacobian =
    [&](const Vector<double> &current_log_stress_ii,
        const Vector<double> & /*current_f*/) {
      // Do nothing here, because we calculate the Jacobian in the residual
      // function
      deallog << "     Recomputing J at x=" << current_log_stress_ii[0]
              << std::endl;
    };


  nonlinear_solver.solve_with_jacobian = [&](const Vector<double> &residual,
                                             Vector<double>       &solution,
                                             const double /*tolerance*/) {
    deallog << "     Solving for dx with residual=" << residual[0] << std::endl;

    solution(0) = residual(0) / log_strain_rate_deriv;
  };



  // Start with the assumption that all strain is accommodated by diffusion
  // creep: If the diffusion creep prefactor is very small, that means that the
  // diffusion viscosity is very large. In this case, use the maximum viscosity
  // instead to compute the starting guess.
  const double stress_ii =
    (prefactor_stress_diffusion > (0.5 / maximum_viscosity) ?
       std::exp(log_edot_ii) / prefactor_stress_diffusion :
       0.5 / maximum_viscosity);

  // KINSOL works on vectors and so the scalar (log) stress is inserted into
  // a vector of length 1
  Vector<double> log_stress_ii(1);
  log_stress_ii[0] = std::log(stress_ii);

  nonlinear_solver.solve(log_stress_ii);
  deallog << n_residual_evaluations << " residual evaluations" << std::endl;
}
