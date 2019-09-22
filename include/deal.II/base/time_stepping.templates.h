// ---------------------------------------------------------------------
//
// Copyright (C) 2014 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------

#ifndef dealii_time_stepping_templates_h
#define dealii_time_stepping_templates_h

#include <deal.II/base/exceptions.h>
#include <deal.II/base/time_stepping.h>

#include <functional>

DEAL_II_NAMESPACE_OPEN

namespace TimeStepping
{
  // ----------------------------------------------------------------------
  // RungeKutta
  // ----------------------------------------------------------------------

  template <typename VectorType>
  double
  RungeKutta<VectorType>::evolve_one_time_step(
    std::vector<std::function<VectorType(const double, const VectorType &)>> &F,
    std::vector<
      std::function<VectorType(const double, const double, const VectorType &)>>
      &J_inverse,

    double      t,
    double      delta_t,
    VectorType &y)
  {
    AssertThrow(
      F.size() == 0,
      ExcMessage(
        "RungeKutta methods cannot handle more that one function to integate."));
    AssertThrow(
      J_inverse.size() == 0,
      ExcMessage(
        "RungeKutta methods cannot handle more that one function to integate."));

    return evolve_one_time_step(F[0], J_inverse[0], t, delta_t, y);
  }



  // ----------------------------------------------------------------------
  // ExplicitRungeKutta
  // ----------------------------------------------------------------------

  template <typename VectorType>
  ExplicitRungeKutta<VectorType>::ExplicitRungeKutta(
    const runge_kutta_method method)
  {
    // virtual functions called in constructors and destructors never use the
    // override in a derived class
    // for clarity be explicit on which function is called
    ExplicitRungeKutta<VectorType>::initialize(method);
  }



  template <typename VectorType>
  void
  ExplicitRungeKutta<VectorType>::initialize(const runge_kutta_method method)
  {
    status.method = method;

    switch (method)
      {
        case (FORWARD_EULER):
          {
            this->n_stages = 1;
            this->b.push_back(1.0);
            this->c.push_back(0.0);

            break;
          }
        case (RK_THIRD_ORDER):
          {
            this->n_stages = 3;
            this->b.reserve(this->n_stages);
            this->c.reserve(this->n_stages);
            this->b.push_back(1.0 / 6.0);
            this->b.push_back(2.0 / 3.0);
            this->b.push_back(1.0 / 6.0);
            this->c.push_back(0.0);
            this->c.push_back(0.5);
            this->c.push_back(1.0);
            std::vector<double> tmp;
            this->a.push_back(tmp);
            tmp.resize(1);
            tmp[0] = 0.5;
            this->a.push_back(tmp);
            tmp.resize(2);
            tmp[0] = -1.0;
            tmp[1] = 2.0;
            this->a.push_back(tmp);

            break;
          }
        case (RK_CLASSIC_FOURTH_ORDER):
          {
            this->n_stages = 4;
            this->b.reserve(this->n_stages);
            this->c.reserve(this->n_stages);
            std::vector<double> tmp;
            this->a.push_back(tmp);
            tmp.resize(1);
            tmp[0] = 0.5;
            this->a.push_back(tmp);
            tmp.resize(2);
            tmp[0] = 0.0;
            tmp[1] = 0.5;
            this->a.push_back(tmp);
            tmp.resize(3);
            tmp[1] = 0.0;
            tmp[2] = 1.0;
            this->a.push_back(tmp);
            this->b.push_back(1.0 / 6.0);
            this->b.push_back(1.0 / 3.0);
            this->b.push_back(1.0 / 3.0);
            this->b.push_back(1.0 / 6.0);
            this->c.push_back(0.0);
            this->c.push_back(0.5);
            this->c.push_back(0.5);
            this->c.push_back(1.0);

            break;
          }
        default:
          {
            AssertThrow(
              false, ExcMessage("Unimplemented explicit Runge-Kutta method."));
          }
      }
  }



  template <typename VectorType>
  double
  ExplicitRungeKutta<VectorType>::evolve_one_time_step(
    const std::function<VectorType(const double, const VectorType &)> &f,
    const std::function<
      VectorType(const double, const double, const VectorType &)>
      & /*id_minus_tau_J_inverse*/,
    double      t,
    double      delta_t,
    VectorType &y)
  {
    return evolve_one_time_step(f, t, delta_t, y);
  }



  template <typename VectorType>
  double
  ExplicitRungeKutta<VectorType>::evolve_one_time_step(
    const std::function<VectorType(const double, const VectorType &)> &f,
    double                                                             t,
    double                                                             delta_t,
    VectorType &                                                       y)
  {
    std::vector<VectorType> f_stages(this->n_stages, y);
    // Compute the different stages needed.
    compute_stages(f, t, delta_t, y, f_stages);

    // Linear combinations of the stages.
    for (unsigned int i = 0; i < this->n_stages; ++i)
      y.sadd(1., delta_t * this->b[i], f_stages[i]);

    return (t + delta_t);
  }



  template <typename VectorType>
  const typename ExplicitRungeKutta<VectorType>::Status &
  ExplicitRungeKutta<VectorType>::get_status() const
  {
    return status;
  }



  template <typename VectorType>
  void
  ExplicitRungeKutta<VectorType>::compute_stages(
    const std::function<VectorType(const double, const VectorType &)> &f,
    const double                                                       t,
    const double                                                       delta_t,
    const VectorType &                                                 y,
    std::vector<VectorType> &f_stages) const
  {
    for (unsigned int i = 0; i < this->n_stages; ++i)
      {
        VectorType Y(y);
        for (unsigned int j = 0; j < i; ++j)
          Y.sadd(1., delta_t * this->a[i][j], f_stages[j]);
        // Evaluate the function f at the point (t+c[i]*delta_t,Y).
        f_stages[i] = f(t + this->c[i] * delta_t, Y);
      }
  }



  // ----------------------------------------------------------------------
  // ImplicitRungeKutta
  // ----------------------------------------------------------------------

  template <typename VectorType>
  ImplicitRungeKutta<VectorType>::ImplicitRungeKutta(
    const runge_kutta_method method,
    const unsigned int       max_it,
    const double             tolerance)
    : RungeKutta<VectorType>()
    , skip_linear_combi(false)
    , max_it(max_it)
    , tolerance(tolerance)
  {
    // virtual functions called in constructors and destructors never use the
    // override in a derived class
    // for clarity be explicit on which function is called
    ImplicitRungeKutta<VectorType>::initialize(method);
  }



  template <typename VectorType>
  void
  ImplicitRungeKutta<VectorType>::initialize(const runge_kutta_method method)
  {
    status.method = method;

    switch (method)
      {
        case (BACKWARD_EULER):
          {
            this->n_stages = 1;
            this->a.push_back(std::vector<double>(1, 1.0));
            this->b.push_back(1.0);
            this->c.push_back(1.0);

            break;
          }
        case (IMPLICIT_MIDPOINT):
          {
            this->a.push_back(std::vector<double>(1, 0.5));
            this->b.push_back(1.0);
            this->c.push_back(0.5);
            this->n_stages = 1;

            break;
          }
        case (CRANK_NICOLSON):
          {
            this->n_stages = 2;
            this->b.reserve(this->n_stages);
            this->c.reserve(this->n_stages);
            this->a.push_back(std::vector<double>(1, 0.0));
            this->a.push_back(std::vector<double>(2, 0.5));
            this->b.push_back(0.5);
            this->b.push_back(0.5);
            this->c.push_back(0.0);
            this->c.push_back(1.0);

            break;
          }
        case (SDIRK_TWO_STAGES):
          {
            this->n_stages = 2;
            this->b.reserve(this->n_stages);
            this->c.reserve(this->n_stages);
            double const gamma = 1.0 - 1.0 / std::sqrt(2.0);
            this->b.push_back(1.0 - gamma);
            this->b.push_back(gamma);
            this->a.push_back(std::vector<double>(1, gamma));
            this->a.push_back(this->b);
            this->c.push_back(gamma);
            this->c.push_back(1.0);

            break;
          }
        default:
          {
            AssertThrow(
              false, ExcMessage("Unimplemented implicit Runge-Kutta method."));
          }
      }
  }



  template <typename VectorType>
  double
  ImplicitRungeKutta<VectorType>::evolve_one_time_step(
    const std::function<VectorType(const double, const VectorType &)> &f,
    const std::function<VectorType(const double,
                                   const double,
                                   const VectorType &)> &id_minus_tau_J_inverse,
    double                                               t,
    double                                               delta_t,
    VectorType &                                         y)
  {
    VectorType              old_y(y);
    std::vector<VectorType> f_stages(this->n_stages, y);
    // Compute the different stages needed.
    compute_stages(f, id_minus_tau_J_inverse, t, delta_t, y, f_stages);

    // If necessary, compute the linear combinations of the stages.
    if (skip_linear_combi == false)
      {
        y = old_y;
        for (unsigned int i = 0; i < this->n_stages; ++i)
          y.sadd(1., delta_t * this->b[i], f_stages[i]);
      }

    return (t + delta_t);
  }



  template <typename VectorType>
  void
  ImplicitRungeKutta<VectorType>::set_newton_solver_parameters(
    unsigned int max_it_,
    double       tolerance_)
  {
    max_it    = max_it_;
    tolerance = tolerance_;
  }



  template <typename VectorType>
  const typename ImplicitRungeKutta<VectorType>::Status &
  ImplicitRungeKutta<VectorType>::get_status() const
  {
    return status;
  }



  template <typename VectorType>
  void
  ImplicitRungeKutta<VectorType>::compute_stages(
    const std::function<VectorType(const double, const VectorType &)> &f,
    const std::function<VectorType(const double,
                                   const double,
                                   const VectorType &)> &id_minus_tau_J_inverse,
    double                                               t,
    double                                               delta_t,
    VectorType &                                         y,
    std::vector<VectorType> &                            f_stages)
  {
    VectorType z(y);
    for (unsigned int i = 0; i < this->n_stages; ++i)
      {
        VectorType old_y(z);
        for (unsigned int j = 0; j < i; ++j)
          old_y.sadd(1., delta_t * this->a[i][j], f_stages[j]);

        // Solve the nonlinear system using Newton's method
        const double new_t       = t + this->c[i] * delta_t;
        const double new_delta_t = this->a[i][i] * delta_t;
        VectorType & f_stage     = f_stages[i];
        newton_solve(
          [this, &f, new_t, new_delta_t, &old_y, &f_stage](
            const VectorType &y, VectorType &residual) {
            this->compute_residual(
              f, new_t, new_delta_t, old_y, y, f_stage, residual);
          },
          [&id_minus_tau_J_inverse, new_t, new_delta_t](const VectorType &y) {
            return id_minus_tau_J_inverse(new_t, new_delta_t, y);
          },
          y);
      }
  }



  template <typename VectorType>
  void
  ImplicitRungeKutta<VectorType>::newton_solve(
    const std::function<void(const VectorType &, VectorType &)> &get_residual,
    const std::function<VectorType(const VectorType &)> &id_minus_tau_J_inverse,
    VectorType &                                         y)
  {
    VectorType residual(y);
    get_residual(y, residual);
    unsigned int i                     = 0;
    const double initial_residual_norm = residual.l2_norm();
    double       norm_residual         = initial_residual_norm;
    while (i < max_it)
      {
        y.sadd(1.0, -1.0, id_minus_tau_J_inverse(residual));
        get_residual(y, residual);
        norm_residual = residual.l2_norm();
        if (norm_residual < tolerance)
          break;
        ++i;
      }
    status.n_iterations  = i + 1;
    status.norm_residual = norm_residual;
  }



  template <typename VectorType>
  void
  ImplicitRungeKutta<VectorType>::compute_residual(
    const std::function<VectorType(const double, const VectorType &)> &f,
    double                                                             t,
    double                                                             delta_t,
    const VectorType &                                                 old_y,
    const VectorType &                                                 y,
    VectorType &                                                       tendency,
    VectorType &residual) const
  {
    // The tendency is stored to save one evaluation of f.
    tendency = f(t, y);
    residual = tendency;
    residual.sadd(delta_t, 1.0, old_y);
    residual.sadd(-1.0, 1., y);
  }



  // ----------------------------------------------------------------------
  // EmbeddedExplicitRungeKutta
  // ----------------------------------------------------------------------

  template <typename VectorType>
  EmbeddedExplicitRungeKutta<VectorType>::EmbeddedExplicitRungeKutta(
    const runge_kutta_method method,
    const double             coarsen_param,
    const double             refine_param,
    const double             min_delta,
    const double             max_delta,
    const double             refine_tol,
    const double             coarsen_tol)
    : coarsen_param(coarsen_param)
    , refine_param(refine_param)
    , min_delta_t(min_delta)
    , max_delta_t(max_delta)
    , refine_tol(refine_tol)
    , coarsen_tol(coarsen_tol)
    , last_same_as_first(false)
    , last_stage(nullptr)
    , status{}
  {
    // virtual functions called in constructors and destructors never use the
    // override in a derived class
    // for clarity be explicit on which function is called
    EmbeddedExplicitRungeKutta<VectorType>::initialize(method);
  }



  template <typename VectorType>
  void
  EmbeddedExplicitRungeKutta<VectorType>::initialize(
    const runge_kutta_method method)
  {
    status.method = method;

    switch (method)
      {
        case (HEUN_EULER):
          {
            this->n_stages = 2;
            this->a.push_back(std::vector<double>());
            this->a.push_back(std::vector<double>(1, 1.0));
            this->c.push_back(0.0);
            this->c.push_back(1.0);
            b1.push_back(0.5);
            b1.push_back(0.5);
            b2.push_back(1.0);
            b2.push_back(0.0);

            break;
          }
        case (BOGACKI_SHAMPINE):
          {
            last_same_as_first = true;
            this->n_stages     = 4;
            this->c.reserve(this->n_stages);
            this->b1.reserve(this->n_stages);
            this->b2.reserve(this->n_stages);
            std::vector<double> tmp;
            this->a.push_back(tmp);
            tmp.resize(1);
            tmp[0] = 0.5;
            this->a.push_back(tmp);
            tmp.resize(2);
            tmp[0] = 0.0;
            tmp[1] = 0.75;
            this->a.push_back(tmp);
            tmp.resize(3);
            tmp[0] = 2.0 / 9.0;
            tmp[1] = 1.0 / 3.0;
            tmp[2] = 4.0 / 9.0;
            this->a.push_back(tmp);
            this->c.push_back(0.0);
            this->c.push_back(0.5);
            this->c.push_back(0.75);
            this->c.push_back(1.0);
            this->b1.push_back(2.0 / 9.0);
            this->b1.push_back(1.0 / 3.0);
            this->b1.push_back(4.0 / 9.0);
            this->b1.push_back(0.0);
            this->b2.push_back(7.0 / 24.0);
            this->b2.push_back(0.25);
            this->b2.push_back(1.0 / 3.0);
            this->b2.push_back(0.125);

            break;
          }
        case (DOPRI):
          {
            last_same_as_first = true;
            this->n_stages     = 7;
            this->c.reserve(this->n_stages);
            this->b1.reserve(this->n_stages);
            this->b2.reserve(this->n_stages);
            std::vector<double> tmp;
            this->a.push_back(tmp);
            tmp.resize(1);
            tmp[0] = 1. / 5.;
            this->a.push_back(tmp);
            tmp.resize(2);
            tmp[0] = 3. / 40.;
            tmp[1] = 9. / 40.;
            this->a.push_back(tmp);
            tmp.resize(3);
            tmp[0] = 44. / 45.;
            tmp[1] = -56. / 15.;
            tmp[2] = 32. / 9.;
            this->a.push_back(tmp);
            tmp.resize(4);
            tmp[0] = 19372. / 6561.;
            tmp[1] = -25360. / 2187.;
            tmp[2] = 64448. / 6561.;
            tmp[3] = -212. / 729.;
            this->a.push_back(tmp);
            tmp.resize(5);
            tmp[0] = 9017. / 3168.;
            tmp[1] = -355. / 33.;
            tmp[2] = 46732. / 5247.;
            tmp[3] = 49. / 176.;
            tmp[4] = -5103. / 18656;
            this->a.push_back(tmp);
            tmp.resize(6);
            tmp[0] = 35. / 384.;
            tmp[1] = 0.;
            tmp[2] = 500. / 1113.;
            tmp[3] = 125. / 192.;
            tmp[4] = -2187. / 6784.;
            tmp[5] = 11. / 84.;
            this->a.push_back(tmp);
            this->c.push_back(0.);
            this->c.push_back(1. / 5.);
            this->c.push_back(3. / 10.);
            this->c.push_back(4. / 5.);
            this->c.push_back(8. / 9.);
            this->c.push_back(1.);
            this->c.push_back(1.);
            this->b1.push_back(35. / 384.);
            this->b1.push_back(0.);
            this->b1.push_back(500. / 1113.);
            this->b1.push_back(125. / 192.);
            this->b1.push_back(-2187. / 6784.);
            this->b1.push_back(11. / 84.);
            this->b1.push_back(0.);
            this->b2.push_back(5179. / 57600.);
            this->b2.push_back(0.);
            this->b2.push_back(7571. / 16695.);
            this->b2.push_back(393. / 640.);
            this->b2.push_back(-92097. / 339200.);
            this->b2.push_back(187. / 2100.);
            this->b2.push_back(1. / 40.);

            break;
          }
        case (FEHLBERG):
          {
            this->n_stages = 6;
            this->c.reserve(this->n_stages);
            this->b1.reserve(this->n_stages);
            this->b2.reserve(this->n_stages);
            std::vector<double> tmp;
            this->a.push_back(tmp);
            tmp.resize(1);
            tmp[0] = 0.25;
            this->a.push_back(tmp);
            tmp.resize(2);
            tmp[0] = 0.09375;
            tmp[1] = 0.28125;
            this->a.push_back(tmp);
            tmp.resize(3);
            tmp[0] = 1932.0 / 2197.0;
            tmp[1] = -7200.0 / 2197.0;
            tmp[2] = 7296.0 / 2197.0;
            this->a.push_back(tmp);
            tmp.resize(4);
            tmp[0] = 439.0 / 216.0;
            tmp[1] = -8.0;
            tmp[2] = 3680.0 / 513.0;
            tmp[3] = -845.0 / 4104.0;
            this->a.push_back(tmp);
            tmp.resize(5);
            tmp[0] = -8.0 / 27.0;
            tmp[1] = 2.0;
            tmp[2] = -3544.0 / 2565.0;
            tmp[3] = 1859.0 / 4104.0;
            tmp[4] = -0.275;
            this->a.push_back(tmp);
            this->c.push_back(0.0);
            this->c.push_back(0.25);
            this->c.push_back(0.375);
            this->c.push_back(12.0 / 13.0);
            this->c.push_back(1.0);
            this->c.push_back(0.5);
            this->b1.push_back(16.0 / 135.0);
            this->b1.push_back(0.0);
            this->b1.push_back(6656.0 / 12825.0);
            this->b1.push_back(28561.0 / 56430.0);
            this->b1.push_back(-0.18);
            this->b1.push_back(2.0 / 55.0);
            this->b2.push_back(25.0 / 216.0);
            this->b2.push_back(0.0);
            this->b2.push_back(1408.0 / 2565.0);
            this->b2.push_back(2197.0 / 4104.0);
            this->b2.push_back(-0.2);
            this->b2.push_back(0.0);

            break;
          }
        case (CASH_KARP):
          {
            this->n_stages = 6;
            this->c.reserve(this->n_stages);
            this->b1.reserve(this->n_stages);
            this->b2.reserve(this->n_stages);
            std::vector<double> tmp;
            this->a.push_back(tmp);
            tmp.resize(1);
            tmp[0] = 0.2;
            this->a.push_back(tmp);
            tmp.resize(2);
            tmp[0] = 0.075;
            tmp[1] = 0.225;
            this->a.push_back(tmp);
            tmp.resize(3);
            tmp[0] = 0.3;
            tmp[1] = -0.9;
            tmp[2] = 1.2;
            this->a.push_back(tmp);
            tmp.resize(4);
            tmp[0] = -11.0 / 54.0;
            tmp[1] = 2.5;
            tmp[2] = -70.0 / 27.0;
            tmp[3] = 35.0 / 27.0;
            this->a.push_back(tmp);
            tmp.resize(5);
            tmp[0] = 1631.0 / 55296.0;
            tmp[1] = 175.0 / 512.0;
            tmp[2] = 575.0 / 13824.0;
            tmp[3] = 44275.0 / 110592.0;
            tmp[4] = 253.0 / 4096.0;
            this->a.push_back(tmp);
            this->c.push_back(0.0);
            this->c.push_back(0.2);
            this->c.push_back(0.3);
            this->c.push_back(0.6);
            this->c.push_back(1.0);
            this->c.push_back(0.875);
            this->b1.push_back(37.0 / 378.0);
            this->b1.push_back(0.0);
            this->b1.push_back(250.0 / 621.0);
            this->b1.push_back(125.0 / 594.0);
            this->b1.push_back(0.0);
            this->b1.push_back(512.0 / 1771.0);
            this->b2.push_back(2825.0 / 27648.0);
            this->b2.push_back(0.0);
            this->b2.push_back(18575.0 / 48384.0);
            this->b2.push_back(13525.0 / 55296.0);
            this->b2.push_back(277.0 / 14336.0);
            this->b2.push_back(0.25);

            break;
          }
        default:
          {
            AssertThrow(
              false, ExcMessage("Unimplemented Embedded Runge-Kutta method."));
          }
      }
  }



  template <typename VectorType>
  void
  EmbeddedExplicitRungeKutta<VectorType>::free_memory()
  {
    if (last_stage != nullptr)
      delete last_stage;

    last_stage = nullptr;
  }



  template <typename VectorType>
  double
  EmbeddedExplicitRungeKutta<VectorType>::evolve_one_time_step(
    const std::function<VectorType(const double, const VectorType &)> &f,
    const std::function<
      VectorType(const double, const double, const VectorType &)>
      & /*id_minus_tau_J_inverse*/,
    double      t,
    double      delta_t,
    VectorType &y)
  {
    return evolve_one_time_step(f, t, delta_t, y);
  }



  template <typename VectorType>
  double
  EmbeddedExplicitRungeKutta<VectorType>::evolve_one_time_step(
    const std::function<VectorType(const double, const VectorType &)> &f,
    double                                                             t,
    double                                                             delta_t,
    VectorType &                                                       y)
  {
    bool                    done       = false;
    unsigned int            count      = 0;
    double                  error_norm = 0.;
    VectorType              old_y(y);
    VectorType              error(y);
    std::vector<VectorType> f_stages(this->n_stages, y);

    while (!done)
      {
        error = 0.;
        y     = old_y;
        // Compute the different stages needed.
        compute_stages(f, t, delta_t, y, f_stages);

        for (unsigned int i = 0; i < this->n_stages; ++i)
          {
            y.sadd(1., delta_t * this->b1[i], f_stages[i]);
            error.sadd(1., delta_t * (b2[i] - b1[i]), f_stages[i]);
          }

        error_norm = error.l2_norm();
        // Check if the norm of error is less than the coarsening tolerance
        if (error_norm < coarsen_tol)
          {
            done = true;
            // Increase the guessed time step
            double new_delta_t = delta_t * coarsen_param;
            // Check that the guessed time step is smaller than the maximum time
            // step allowed.
            if (new_delta_t > max_delta_t)
              {
                status.exit_delta_t  = MAX_DELTA_T;
                status.delta_t_guess = max_delta_t;
              }
            else
              {
                status.exit_delta_t  = DELTA_T;
                status.delta_t_guess = new_delta_t;
              }
          }
        // Check if the norm of error is less than the refining tolerance
        else if (error_norm < refine_tol)
          {
            done                 = true;
            status.exit_delta_t  = DELTA_T;
            status.delta_t_guess = delta_t;
          }
        else
          {
            // If the time step is already the smallest acceptable, exit.
            if (delta_t == min_delta_t)
              {
                done                 = true;
                status.exit_delta_t  = MIN_DELTA_T;
                status.delta_t_guess = delta_t;
              }
            // Reduce the time step.
            else
              {
                delta_t *= refine_param;
                if (delta_t < min_delta_t)
                  delta_t = min_delta_t;
              }
          }

        ++count;
      }

    // Save the last stage if necessary
    if (last_same_as_first == true)
      {
        if (last_stage == nullptr)
          last_stage = new VectorType(f_stages.back());
        else
          *last_stage = f_stages.back();
      }

    status.n_iterations = count;
    status.error_norm   = error_norm;

    return (t + delta_t);
  }



  template <typename VectorType>
  void
  EmbeddedExplicitRungeKutta<VectorType>::set_time_adaptation_parameters(
    const double coarsen_param_,
    const double refine_param_,
    const double min_delta_,
    const double max_delta_,
    const double refine_tol_,
    const double coarsen_tol_)
  {
    coarsen_param = coarsen_param_;
    refine_param  = refine_param_;
    min_delta_t   = min_delta_;
    max_delta_t   = max_delta_;
    refine_tol    = refine_tol_;
    coarsen_tol   = coarsen_tol_;
  }



  template <typename VectorType>
  const typename EmbeddedExplicitRungeKutta<VectorType>::Status &
  EmbeddedExplicitRungeKutta<VectorType>::get_status() const
  {
    return status;
  }


  template <typename VectorType>
  void
  EmbeddedExplicitRungeKutta<VectorType>::compute_stages(
    const std::function<VectorType(const double, const VectorType &)> &f,
    const double                                                       t,
    const double                                                       delta_t,
    const VectorType &                                                 y,
    std::vector<VectorType> &                                          f_stages)
  {
    VectorType   Y(y);
    unsigned int i = 0;

    // If the last stage is the same as the first, we can skip the evaluation
    // of the first stage.
    if (last_same_as_first == true)
      {
        if (last_stage != nullptr)
          {
            f_stages[0] = *last_stage;
            i           = 1;
          }
      }

    for (; i < this->n_stages; ++i)
      {
        Y = y;
        for (unsigned int j = 0; j < i; ++j)
          Y.sadd(1.0, delta_t * this->a[i][j], f_stages[j]);
        f_stages[i] = f(t + this->c[i] * delta_t, Y);
      }
  }
} // namespace TimeStepping

DEAL_II_NAMESPACE_CLOSE

#endif
