// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2014 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// test Runge-Kutta methods in TimeStepping with a) a polynomial with expected
// error 0 and b) convergence order for y=0.1*exp(t^2)
#include <deal.II/base/time_stepping.h>

#include <deal.II/lac/vector.h>

#include "../tests.h"

Vector<double>
f1(const double t, const Vector<double> &y)
{
  Vector<double> values(y);
  for (unsigned int i = 0; i < values.size(); ++i)
    values[i] = 1.0;

  return values;
}

Vector<double>
f2(const double t, const Vector<double> &y)
{
  Vector<double> values(y);
  for (unsigned int i = 0; i < values.size(); ++i)
    values[i] = 2.0 * t;

  return values;
}

Vector<double>
f3(const double t, const Vector<double> &y)
{
  Vector<double> values(y);
  for (unsigned int i = 0; i < values.size(); ++i)
    values[i] = 3.0 * t * t;

  return values;
}

Vector<double>
f4(const double t, const Vector<double> &y)
{
  Vector<double> values(y);
  for (unsigned int i = 0; i < values.size(); ++i)
    values[i] = 4.0 * t * t * t;

  return values;
}

Vector<double>
f5(const double t, const Vector<double> &y)
{
  Vector<double> values(y);
  for (unsigned int i = 0; i < values.size(); ++i)
    values[i] = 5.0 * t * t * t * t;

  return values;
}

Vector<double>
f6(const double t, const Vector<double> &y)
{
  Vector<double> values(y);
  for (unsigned int i = 0; i < values.size(); ++i)
    values[i] = 6.0 * t * t * t * t * t;

  return values;
}


Vector<double>
my_rhs_function(const double t, const Vector<double> &y)
{
  Vector<double> values(y);
  for (unsigned int i = 0; i < values.size(); ++i)
    values[i] = y[i] * 2.0 * t;

  return values;
}

Vector<double>
id_minus_tau_J_inv1(const double t, const double tau, const Vector<double> &y)
{
  return y;
}

Vector<double>
id_minus_tau_J_inv2(const double t, const double tau, const Vector<double> &y)
{
  return y;
}

Vector<double>
id_minus_tau_J_inv3(const double t, const double tau, const Vector<double> &y)
{
  return y;
}

Vector<double>
id_minus_tau_J_inv4(const double t, const double tau, const Vector<double> &y)
{
  return y;
}

Vector<double>
id_minus_tau_J_inv5(const double t, const double tau, const Vector<double> &y)
{
  return y;
}

Vector<double>
id_minus_tau_J_inv6(const double t, const double tau, const Vector<double> &y)
{
  return y;
}

double
my1(const double t)
{
  return t;
}

double
my2(const double t)
{
  return t * t;
}

double
my3(const double t)
{
  return t * t * t;
}

double
my4(const double t)
{
  return t * t * t * t;
}

double
my5(const double t)
{
  return t * t * t * t * t;
}

double
my6(const double t)
{
  return t * t * t * t * t * t;
}

double
my_exact_solution(const double t)
{
  return 0.1 * std::exp(t * t);
}

void
test(TimeStepping::RungeKutta<Vector<double>>                           &solver,
     std::function<Vector<double>(const double, const Vector<double> &)> f,
     std::function<Vector<double>(const double,
                                  const double,
                                  const Vector<double> &)> id_minus_tau_J_inv,
     std::function<double(const double)>                   my)
{
  unsigned int n_time_steps = 1;
  unsigned int size         = 1;
  double       initial_time = 0.0, final_time = 1.0;
  double       time_step =
    (final_time - initial_time) / static_cast<double>(n_time_steps);
  double         time = initial_time;
  Vector<double> solution(size);
  Vector<double> exact_solution(size);
  for (unsigned int i = 0; i < size; ++i)
    {
      solution[i]       = my(initial_time);
      exact_solution[i] = my(final_time);
    }

  for (unsigned int i = 0; i < n_time_steps; ++i)
    time = solver.evolve_one_time_step(
      f, id_minus_tau_J_inv, time, time_step, solution);

  Vector<double> error(exact_solution);
  error.sadd(1.0, -1.0, solution);
  double error_norm = error.l2_norm();
  deallog << error_norm << std::endl;
}

void
test2(TimeStepping::EmbeddedExplicitRungeKutta<Vector<double>> &solver,
      std::function<Vector<double>(const double, const Vector<double> &)> f,
      std::function<Vector<double>(const double,
                                   const double,
                                   const Vector<double> &)> id_minus_tau_J_inv,
      std::function<double(const double)>                   my)
{
  double         initial_time = 0.0, final_time = 1.0;
  double         time_step    = 1.0;
  unsigned int   size         = 1;
  unsigned int   n_time_steps = 0;
  double         time         = initial_time;
  Vector<double> solution(size);
  Vector<double> exact_solution(size);
  for (unsigned int i = 0; i < size; ++i)
    {
      solution[i]       = my(initial_time);
      exact_solution[i] = my(final_time);
    }


  while (time < final_time)
    {
      if (time + time_step > final_time)
        time_step = final_time - time;
      time = solver.evolve_one_time_step(
        f, id_minus_tau_J_inv, time, time_step, solution);
      time_step = solver.get_status().delta_t_guess;
    }

  Vector<double> error(exact_solution);
  error.sadd(1.0, -1.0, solution);
  double error_norm = error.l2_norm();
  deallog << error_norm << std::endl;
}


void
test_convergence(
  TimeStepping::RungeKutta<Vector<double>>                           &solver,
  std::function<Vector<double>(const double, const Vector<double> &)> f,
  std::function<Vector<double>(const double,
                               const double,
                               const Vector<double> &)> id_minus_tau_J_inv,
  std::function<double(const double)>                   my)
{
  std::vector<double> errors;
  double              initial_time = 0.0, final_time = 1.0;
  unsigned int        size = 1;
  Vector<double>      solution(size);
  Vector<double>      exact_solution(size);
  for (unsigned int i = 0; i < size; ++i)
    {
      exact_solution[i] = my(final_time);
    }

  deallog << "convergence rate" << std::endl;
  for (unsigned int cycle = 0; cycle < 8; ++cycle)
    {
      unsigned int n_time_steps = std::pow(2., static_cast<double>(cycle));
      double       time_step =
        (final_time - initial_time) / static_cast<double>(n_time_steps);
      double time = initial_time;
      for (unsigned int i = 0; i < size; ++i)
        solution[i] = my(initial_time);

      for (unsigned int i = 0; i < n_time_steps; ++i)
        time = solver.evolve_one_time_step(
          f, id_minus_tau_J_inv, time, time_step, solution);

      Vector<double> error(exact_solution);
      error.sadd(1.0, -1.0, solution);
      double error_norm = error.l2_norm();
      errors.push_back(error_norm);
      if (cycle > 1)
        deallog << std::log(std::fabs(errors[cycle - 1] / errors[cycle])) /
                     std::log(2.)
                << std::endl;
    }
}

int
main()
{
  initlog();
  // deallog.precision(4);
  {
    deallog << "Forward Euler" << std::endl;
    TimeStepping::ExplicitRungeKutta<Vector<double>> fe(
      TimeStepping::FORWARD_EULER);
    test(fe, f1, id_minus_tau_J_inv1, my1);

    deallog << "Runge-Kutta third order" << std::endl;
    TimeStepping::ExplicitRungeKutta<Vector<double>> rk3(
      TimeStepping::RK_THIRD_ORDER);
    test(rk3, f3, id_minus_tau_J_inv3, my3);

    deallog << "Strong Stability Preserving Runge-Kutta third order"
            << std::endl;
    TimeStepping::ExplicitRungeKutta<Vector<double>> ssp_rk3(
      TimeStepping::SSP_THIRD_ORDER);
    test(ssp_rk3, f3, id_minus_tau_J_inv3, my3);

    deallog << "Runge-Kutta fourth order" << std::endl;
    TimeStepping::ExplicitRungeKutta<Vector<double>> rk4(
      TimeStepping::RK_CLASSIC_FOURTH_ORDER);
    test(rk4, f4, id_minus_tau_J_inv4, my4);

    deallog << "Runge-Kutta fifth order" << std::endl;
    TimeStepping::ExplicitRungeKutta<Vector<double>> rk5(
      TimeStepping::RK_FIFTH_ORDER);
    test(rk5, f5, id_minus_tau_J_inv5, my5);

    deallog << "Runge-Kutta sixth order" << std::endl;
    TimeStepping::ExplicitRungeKutta<Vector<double>> rk6(
      TimeStepping::RK_SIXTH_ORDER);
    test(rk5, f6, id_minus_tau_J_inv6, my6);

    deallog << "Low-storage Runge-Kutta stage 3 order 3" << std::endl;
    TimeStepping::LowStorageRungeKutta<Vector<double>> lsrk33(
      TimeStepping::LOW_STORAGE_RK_STAGE3_ORDER3);
    test(lsrk33, f3, id_minus_tau_J_inv3, my3);

    deallog << "Low-storage Runge-Kutta stage 5 order 4" << std::endl;
    TimeStepping::LowStorageRungeKutta<Vector<double>> lsrk54(
      TimeStepping::LOW_STORAGE_RK_STAGE5_ORDER4);
    test(lsrk54, f4, id_minus_tau_J_inv4, my4);

    deallog << "Low-storage Runge-Kutta stage 7 order 4" << std::endl;
    TimeStepping::LowStorageRungeKutta<Vector<double>> lsrk74(
      TimeStepping::LOW_STORAGE_RK_STAGE7_ORDER4);
    test(lsrk74, f4, id_minus_tau_J_inv4, my4);

    deallog << "Low-storage Runge-Kutta stage 9 order 5" << std::endl;
    TimeStepping::LowStorageRungeKutta<Vector<double>> lsrk95(
      TimeStepping::LOW_STORAGE_RK_STAGE9_ORDER5);
    test(lsrk95, f5, id_minus_tau_J_inv5, my5);

    deallog << "Backward Euler" << std::endl;
    TimeStepping::ImplicitRungeKutta<Vector<double>> be(
      TimeStepping::BACKWARD_EULER);
    test(be, f1, id_minus_tau_J_inv1, my1);

    deallog << "Implicit midpoint" << std::endl;
    TimeStepping::ImplicitRungeKutta<Vector<double>> im(
      TimeStepping::IMPLICIT_MIDPOINT);
    test(im, f2, id_minus_tau_J_inv2, my2);

    deallog << "Crank-Nicolson" << std::endl;
    TimeStepping::ImplicitRungeKutta<Vector<double>> cn(
      TimeStepping::CRANK_NICOLSON);
    test(cn, f2, id_minus_tau_J_inv2, my2);

    deallog << "SDIRK" << std::endl;
    TimeStepping::ImplicitRungeKutta<Vector<double>> sdirk(
      TimeStepping::SDIRK_TWO_STAGES);
    test(sdirk, f2, id_minus_tau_J_inv2, my2);

    deallog << "Heun-Euler" << std::endl;
    TimeStepping::EmbeddedExplicitRungeKutta<Vector<double>> he(
      TimeStepping::HEUN_EULER);
    test2(he, f2, id_minus_tau_J_inv2, my2);

    deallog << "Bogacki-Shampine" << std::endl;
    TimeStepping::EmbeddedExplicitRungeKutta<Vector<double>> bs(
      TimeStepping::BOGACKI_SHAMPINE);
    test2(bs, f3, id_minus_tau_J_inv3, my3);
    bs.free_memory();

    deallog << "DOPRI" << std::endl;
    TimeStepping::EmbeddedExplicitRungeKutta<Vector<double>> dopri(
      TimeStepping::DOPRI);
    test2(dopri, f5, id_minus_tau_J_inv5, my5);
    dopri.free_memory();

    deallog << "Fehlberg" << std::endl;
    TimeStepping::EmbeddedExplicitRungeKutta<Vector<double>> fehlberg(
      TimeStepping::FEHLBERG);
    test2(fehlberg, f5, id_minus_tau_J_inv5, my5);

    deallog << "Cash-Karp" << std::endl;
    TimeStepping::EmbeddedExplicitRungeKutta<Vector<double>> ck(
      TimeStepping::CASH_KARP);
    test2(ck, f5, id_minus_tau_J_inv5, my5);
  }

  {
    deallog << "Forward Euler first order" << std::endl;
    TimeStepping::ExplicitRungeKutta<Vector<double>> rk1(
      TimeStepping::FORWARD_EULER);
    test_convergence(rk1,
                     my_rhs_function,
                     id_minus_tau_J_inv1,
                     my_exact_solution);

    deallog << "Runge-Kutta third order" << std::endl;
    TimeStepping::ExplicitRungeKutta<Vector<double>> rk3(
      TimeStepping::RK_THIRD_ORDER);
    test_convergence(rk3,
                     my_rhs_function,
                     id_minus_tau_J_inv3,
                     my_exact_solution);

    deallog << "Strong Stability Preserving Runge-Kutta third order"
            << std::endl;
    TimeStepping::ExplicitRungeKutta<Vector<double>> ssp_rk3(
      TimeStepping::SSP_THIRD_ORDER);
    test_convergence(ssp_rk3,
                     my_rhs_function,
                     id_minus_tau_J_inv3,
                     my_exact_solution);

    deallog << "Runge-Kutta fourth order" << std::endl;
    TimeStepping::ExplicitRungeKutta<Vector<double>> rk4(
      TimeStepping::RK_CLASSIC_FOURTH_ORDER);
    test_convergence(rk4,
                     my_rhs_function,
                     id_minus_tau_J_inv4,
                     my_exact_solution);

    deallog << "Runge-Kutta fifth order" << std::endl;
    TimeStepping::ExplicitRungeKutta<Vector<double>> rk5(
      TimeStepping::RK_FIFTH_ORDER);
    test_convergence(rk5,
                     my_rhs_function,
                     id_minus_tau_J_inv5,
                     my_exact_solution);

    deallog << "Runge-Kutta sixth order" << std::endl;
    TimeStepping::ExplicitRungeKutta<Vector<double>> rk6(
      TimeStepping::RK_SIXTH_ORDER);
    test_convergence(rk6,
                     my_rhs_function,
                     id_minus_tau_J_inv6,
                     my_exact_solution);

    deallog << "Low-storage Runge-Kutta stage 3 order 3" << std::endl;
    TimeStepping::LowStorageRungeKutta<Vector<double>> lsrk33(
      TimeStepping::LOW_STORAGE_RK_STAGE3_ORDER3);
    test_convergence(lsrk33,
                     my_rhs_function,
                     id_minus_tau_J_inv3,
                     my_exact_solution);

    deallog << "Low-storage Runge-Kutta stage 5 order 4" << std::endl;
    TimeStepping::LowStorageRungeKutta<Vector<double>> lsrk54(
      TimeStepping::LOW_STORAGE_RK_STAGE5_ORDER4);
    test_convergence(lsrk54,
                     my_rhs_function,
                     id_minus_tau_J_inv4,
                     my_exact_solution);

    deallog << "Low-storage Runge-Kutta stage 7 order 4" << std::endl;
    TimeStepping::LowStorageRungeKutta<Vector<double>> lsrk74(
      TimeStepping::LOW_STORAGE_RK_STAGE7_ORDER4);
    test_convergence(lsrk74,
                     my_rhs_function,
                     id_minus_tau_J_inv4,
                     my_exact_solution);

    deallog << "Low-storage Runge-Kutta stage 9 order 5" << std::endl;
    TimeStepping::LowStorageRungeKutta<Vector<double>> lsrk95(
      TimeStepping::LOW_STORAGE_RK_STAGE9_ORDER5);
    test_convergence(lsrk95,
                     my_rhs_function,
                     id_minus_tau_J_inv5,
                     my_exact_solution);

    deallog << "Backward Euler first order" << std::endl;
    TimeStepping::ImplicitRungeKutta<Vector<double>> be(
      TimeStepping::BACKWARD_EULER);
    test_convergence(be,
                     my_rhs_function,
                     id_minus_tau_J_inv1,
                     my_exact_solution);

    deallog << "Implicit midpoint second order" << std::endl;
    TimeStepping::ImplicitRungeKutta<Vector<double>> im(
      TimeStepping::IMPLICIT_MIDPOINT);
    test_convergence(im,
                     my_rhs_function,
                     id_minus_tau_J_inv2,
                     my_exact_solution);

    deallog << "Crank-Nicolson second order" << std::endl;
    TimeStepping::ImplicitRungeKutta<Vector<double>> cn(
      TimeStepping::CRANK_NICOLSON);
    test_convergence(cn,
                     my_rhs_function,
                     id_minus_tau_J_inv2,
                     my_exact_solution);

    deallog << "SDIRK second order" << std::endl;
    TimeStepping::ImplicitRungeKutta<Vector<double>> sdirk(
      TimeStepping::SDIRK_TWO_STAGES);
    test_convergence(sdirk,
                     my_rhs_function,
                     id_minus_tau_J_inv2,
                     my_exact_solution);
  }

  return 0;
}
