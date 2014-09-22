// ---------------------------------------------------------------------
//
// Copyright (C) 2014 by the deal.II authors
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

#ifndef __deal2__time_stepping_h
#define __deal2__time_stepping_h


#include <deal.II/base/config.h>
#include <deal.II/base/std_cxx11/function.h>

#include <vector>

DEAL_II_NAMESPACE_OPEN

/**
 * Namespace containing the time stepping methods.
 *
 * @author Bruno Turcksin
 * @date 2014
 */

namespace TimeStepping
{
  /**
   * Runge-Kutta methods available:
   *   - Explicit methods:
   *     - FORWARD_EULER: first order
   *     - RK_THIRD_ORDER: third order Runge-Kutta
   *     - RK_CLASSIC_FOURTH_ORDER: classical fourth order Runge-Kutta
   *   - Implicit methods:
   *     - BACKWARD_EULER: first order
   *     - IMPLICIT_MIDPOINT: second order
   *     - CRANK_NICOLSON: second order
   *     - SDIRK_TWO_STAGES: second order
   *   - Embedded explicit methods:
   *     - HEUN_EULER: second order
   *     - BOGACKI_SHAMPINE: third order
   *     - DOPRI: Dormand-Prince fifth order (method used by ode45 in MATLAB)
   *     - FEHLBERG: fifth order
   *     - CASH_KARP: firth order
   */
  enum runge_kutta_method { FORWARD_EULER, RK_THIRD_ORDER, RK_CLASSIC_FOURTH_ORDER,
                            BACKWARD_EULER, IMPLICIT_MIDPOINT, CRANK_NICOLSON,
                            SDIRK_TWO_STAGES, HEUN_EULER, BOGACKI_SHAMPINE, DOPRI,
                            FEHLBERG, CASH_KARP
                          };



  /**
   * Reason for exiting evolve_one_time_step when using an embedded
   * method: DELTA_T (the time step is in the valid range), MIN_DELTA_T (the
   * time step was increased to the minimum acceptable time step), MAX_DELTA_T
   * (the time step was reduced to the maximum acceptable time step).
   */
  enum embedded_runge_kutta_time_step { DELTA_T, MIN_DELTA_T, MAX_DELTA_T };



  /**
   * Abstract class for time stepping methods. These methods assume that the
   * equation has the form: $ \frac{\partial y}{\partial t} = f(t,y) $.
   */
  template <typename VECTOR>
  class TimeStepping
  {
  public:
    /**
     * Virtual destructor.
     */
    virtual ~TimeStepping() {};

    /**
     * Purely virtual function. This function is used to advance from time @p
     * t to t+ @p delta_t. @p F is a vector of functions $ f(t,y) $ that should be
     * integrated, the input parameters are the time t and the vector y and the
     * output is value of f at this point. @p J_inverse is a vector
     * functions that compute the inverse of the Jacobians associated to the
     * implicit problems. The input parameters are the
     * time, $ \tau $, and a vector. The output is the value of function
     * at this point. This function returns the time at the end of the
     * time step.
     */
    virtual double evolve_one_time_step(
      std::vector<std_cxx11::function<VECTOR (const double, const VECTOR &)> > &F,
      std::vector<std_cxx11::function<VECTOR (const double, const double, const VECTOR &)> > &J_inverse,
      double t,
      double delta_t,
      VECTOR &y) = 0;

    /**
     * Empty structure used to store informations.
     */
    struct Status {};

    /**
     * Purely virtual function that return Status.
     */
    virtual const Status &get_status() const = 0;
  };



  /**
   * Base class for the Runge-Kutta method
   *
   * @author Damien Lebrun-Grandie, Bruno Turcksin
   * @date 2014
   */
  template <typename VECTOR>
  class RungeKutta : public TimeStepping<VECTOR>
  {
  public:
    /**
     * Virtual destructor.
     */
    virtual ~RungeKutta() {};

    /**
     * Purely virtual method used to initialize the Runge-Kutta method.
     */
    virtual void initialize(runge_kutta_method method) = 0;
    /**
     * This function is used to advance from time @p
     * t to t+ @p delta_t. @p F is a vector of functions $ f(t,y) $ that should be
     * integrated, the input parameters are the time t and the vector y and the
     * output is value of f at this point. @p J_inverse is a vector
     * functions that compute the inverse of the Jacobians associated to the
     * implicit problems. The input parameters are the
     * time, $ \tau $, and a vector. The output is the value of function
     * at this point. This function returns the time at the end of the
     * time step. When using Runge-Kutta methods, @p F and @ J_inverse can
     * only contain one element.
     */
    double evolve_one_time_step(
      std::vector<std_cxx11::function<VECTOR (const double, const VECTOR &)> > &F,
      std::vector<std_cxx11::function<VECTOR (const double, const double, const VECTOR &)> > &J_inverse,
      double t,
      double delta_t,
      VECTOR &y);

    /**
     * Purely virtual function. This function is used to advance from time @p t
     * to t+ @p delta_t. @p f  is the function $ f(t,y) $ that should be
     * integrated, the input parameters are the time t and the vector y and the
     * output is value of f at this point. @p id_minus_tau_J_inverse is a function
     * that computes $ inv(I-\tau J)$ where $ I $ is the identity matrix,
     * $ \tau $ is given, and $ J $ is the Jacobian $ \frac{\partial
     * J}{\partial y} $. The input parameters are the time, $ \tau $, and
     * a vector. The output is the value of function at this point.
     * evolve_one_time_step returns the time at the end of the time step.
     */
    virtual double evolve_one_time_step(
      std_cxx11::function<VECTOR (const double, const VECTOR &)> f,
      std_cxx11::function<VECTOR (const double, const double, const VECTOR &)> id_minus_tau_J_inverse,
      double t,
      double delta_t,
      VECTOR &y) = 0;

  protected:
    /**
     * Number of stages of the Runge-Kutta method.
     */
    unsigned int n_stages;

    /**
     * Butcher tableau coefficients.
     */
    std::vector<double> b;

    /**
     * Butcher tableau coefficients.
     */
    std::vector<double> c;

    /**
     * Butcher tableau coefficients.
     */
    std::vector<std::vector<double> > a;
  };



  /**
   * ExplicitRungeKutta is derived from RungeKutta and implement the explicit methods.
   */
  template <typename VECTOR>
  class ExplicitRungeKutta : public RungeKutta<VECTOR>
  {
  public:
    using RungeKutta<VECTOR>::evolve_one_time_step;

    /**
     * Default constructor. initialize(runge_kutta_method) needs to be called
     * before the object can be used.
     */
    ExplicitRungeKutta() {};

    /**
     * Constructor. This function calls initialize(runge_kutta_method).
     */
    ExplicitRungeKutta(runge_kutta_method method);

    /**
     * Initialize the explicit Runge-Kutta method.
     */
    void initialize(runge_kutta_method method);

    /**
     * This function is used to advance from time @p t to t+ @p delta_t. @p f
     * is the function $ f(t,y) $ that should be integrated, the input
     * parameters are the time t and the vector y and the output is value of
     * f at this point. @p id_minus_tau_J_inverse is a function that computes
     * $ inv(I-\tau J)$ where $ I $ is the identity matrix, $ \tau
     * $ is given, and $ J $ is the Jacobian $ \frac{\partial
     * J}{\partial y} $. The input parameter are the time, $ \tau $, and
     * a vector. The output is the value of function at this point.
     * evolve_one_time_step returns the time at the end of the time step.
     */
    double evolve_one_time_step(
      std_cxx11::function<VECTOR (const double, const VECTOR &)> f,
      std_cxx11::function<VECTOR (const double, const double, const VECTOR &)> id_minus_tau_J_inverse,
      double t,
      double delta_t,
      VECTOR &y);

    /**
     * This function is used to advance from time @p t to t+ @p delta_t.
     * This function is similar to the one derived from RungeKutta, but
     * does not required id_minus_tau_J_inverse because it is not used for
     * explicit methods. evolve_one_time_step returns the time at the end of the
     * time step.
     */
    double evolve_one_time_step(std_cxx11::function<VECTOR (const double, const VECTOR &)> f,
                                double t,
                                double delta_t,
                                VECTOR &y);

    /**
     * This structure stores the name of the method used.
     */
    struct Status : public TimeStepping<VECTOR>::Status
    {
      runge_kutta_method method;
    };

    /**
     * Return the status of the current object.
     */
    const Status &get_status() const;

  private:
    /**
     * Compute the different stages needed.
     */
    void compute_stages(std_cxx11::function<VECTOR (const double, const VECTOR &)> f,
                        const double t,
                        const double delta_t,
                        const VECTOR &y,
                        std::vector<VECTOR> &f_stages) const;

    /**
     * Status structure of the object.
     */
    Status status;
  };



  /**
   * This class is derived from RungeKutta and implement the implicit
   * methods. This class works only for Diagonal Implicit Runge-Kutta
   * (DIRK) methods.
   */
  template <typename VECTOR>
  class ImplicitRungeKutta : public RungeKutta<VECTOR>
  {
  public:
    using RungeKutta<VECTOR>::evolve_one_time_step;

    /**
     * Default constructor. initialize(runge_kutta_method) and
     * set_newton_solver_parameters(unsigned int,double)
     * need to be called before the object can be used.
     */
    ImplicitRungeKutta() {};

    /**
     * Constructor. This function calls initialize(runge_kutta_method) and
     * initialize the maximum number of iterations and the tolerance of the
     * Newton solver.
     */
    ImplicitRungeKutta(runge_kutta_method method, unsigned int max_it=100, double tolerance=1e-6);

    /**
     * Initialize the implicit Runge-Kutta method.
     */
    void initialize(runge_kutta_method method);

    /**
     * This function is used to advance from time @p t to t+ @p delta_t. @p f
     * is the function $ f(t,y) $ that should be integrated, the input
     * parameters are the time t and the vector y and the output is value of
     * f at this point. @p id_minus_tau_J_inverse is a function that computes
     * $ inv(I-\tau J)$ where $ I $ is the identity matrix, $ \tau
     * $ is given, and $ J $ is the Jacobian $ \frac{\partial
     * J}{\partial y} $. The input parameters are the time, $ \tau $, and
     * a vector. The output is the value of function at this point.
     * evolve_one_time_step returns the time at the end of the time step.
     */
    double evolve_one_time_step(
      std_cxx11::function<VECTOR (const double, const VECTOR &)> f,
      std_cxx11::function<VECTOR (const double, const double, const VECTOR &)> id_minus_tau_J_inverse,
      double t,
      double delta_t,
      VECTOR &y);

    /**
     * Set the maximum number of iterations and the tolerance used by the
     * Newton solver.
     */
    void set_newton_solver_parameters(unsigned int max_it, double tolerance);

    /**
     * Structure that stores the name of the method, the number of
     * Newton iterations and the norm of the residual when exiting the
     * Newton solver.
     */
    struct Status : public TimeStepping<VECTOR>::Status
    {
      runge_kutta_method method;
      unsigned int       n_iterations;
      double             norm_residual;
    };

    /**
     * Return the status of the current object.
     */
    const Status &get_status() const;

  private:
    /**
     * Compute the different stages needed.
     */
    void compute_stages(
      std_cxx11::function<VECTOR (const double, const VECTOR &)> f,
      std_cxx11::function<VECTOR (const double, const double, const VECTOR &)> id_minus_tau_J_inverse,
      double t,
      double delta_t,
      VECTOR &y,
      std::vector<VECTOR> &f_stages);

    /**
     * Newton solver used for the implicit stages.
     */
    void newton_solve(std_cxx11::function<void (const VECTOR &,VECTOR &)> get_residual,
                      std_cxx11::function<VECTOR (const VECTOR &)> id_minus_tau_J_inverse,
                      VECTOR &y);

    /**
     * Compute the residual needed by the Newton solver.
     */
    void compute_residual(std_cxx11::function<VECTOR (const double, const VECTOR &)> f,
                          double t,
                          double delta_t,
                          const VECTOR &old_y,
                          const VECTOR &y,
                          VECTOR &tendency,
                          VECTOR &residual) const;

    /**
     * When using SDIRK, there is no need to compute the linear combination
     * of the stages. Thus, when this flag is true, the linear combination
     * is skipped.
     */
    bool skip_linear_combi;

    /**
     * Maximum number of iterations of the Newton solver.
     */
    unsigned int max_it;

    /**
     * Tolerance of the Newton solver.
     */
    double tolerance;

    /**
     * Status structure of the object.
     */
    Status status;
  };



  /**
   * This is class is derived from RungeKutta and implement embedded explicit
   * methods.
   */
  template <typename VECTOR>
  class EmbeddedExplicitRungeKutta : public RungeKutta<VECTOR>
  {
  public:
    using RungeKutta<VECTOR>::evolve_one_time_step;

    /**
     * Default constructor. initialize(runge_kutta_method) and
     * set_time_adaptation_parameters(double, double, double, double, double, double)
     * need to be called before the object can be used.
     */
    EmbeddedExplicitRungeKutta() {};

    /**
     * Constructor. This function calls initialize(runge_kutta_method) and
     * initialize the parameters needed for time adaptation.
     */
    EmbeddedExplicitRungeKutta(runge_kutta_method method,
                               double coarsen_param = 1.2,
                               double refine_param = 0.8,
                               double min_delta = 1e-14,
                               double max_delta = 1e100,
                               double refine_tol = 1e-8,
                               double coarsen_tol = 1e-12);

    /**
     * Destructor.
     */
    ~EmbeddedExplicitRungeKutta()
    {
      free_memory();
    }

    /**
     * If necessary, deallocate memory allocated by the object.
     */
    void free_memory();

    /**
     * Initialize the embedded explicit Runge-Kutta method.
     */
    void initialize(runge_kutta_method method);

    /**
     * This function is used to advance from time @p t to t+ @p delta_t. @p f
     * is the function $ f(t,y) $ that should be integrated, the input
     * parameters are the time t and the vector y and the output is value of
     * f at this point. @p id_minus_tau_J_inverse is a function that computes
     * $ inv(I-\tau J)$ where $ I $ is the identity matrix, $ \tau
     * $ is given, and $ J $ is the Jacobian $ \frac{\partial
     * J}{\partial y} $. The input parameters are the time, $ \tau $, and
     * a vector. The output is the value of function at this point.
     * evolve_one_time_step returns the time at the end of the time step.
     */
    double evolve_one_time_step(
      std_cxx11::function<VECTOR (const double, const VECTOR &)> f,
      std_cxx11::function<VECTOR (const double, const double, const VECTOR &)> id_minus_tau_J_inverse,
      double t,
      double delta_t,
      VECTOR &y);

    /**
     * This function is used to advance from time @p t to t+ @p delta_t.
     * This function is similar to the one derived from TimeStepping, but
     * does not required id_minus_tau_J_inverse because it is not used for
     * explicit methods. evolve_one_time_step returns the time at the end of the
     * time step.
     */
    double evolve_one_time_step(std_cxx11::function<VECTOR (const double, const VECTOR &)> f,
                                double t,
                                double delta_t,
                                VECTOR &y);

    /**
     * Set the parameters necessary for the time adaptation.
     */
    void set_time_adaptation_parameters(double coarsen_param,
                                        double refine_param,
                                        double min_delta,
                                        double max_delta,
                                        double refine_tol,
                                        double coarsen_tol);

    /**
     * Structure that stores the name of the method, the reason to exit
     * evolve_one_time_step, the number of iteration inside n_iterations, a guess
     * of what the next time step should be, and an estimate of the norm of
     * the error.
     */
    struct Status : public TimeStepping<VECTOR>::Status
    {
      runge_kutta_method method;
      embedded_runge_kutta_time_step exit_delta_t;
      unsigned int n_iterations;
      double delta_t_guess;
      double error_norm;
    };

    /**
     * Return the status of the current object.
     */
    const Status &get_status() const;

  private:
    /**
     * Compute the different stages needed.
     */
    void compute_stages(std_cxx11::function<VECTOR (const double, const VECTOR &)> f,
                        const double t,
                        const double delta_t,
                        const VECTOR &y,
                        std::vector<VECTOR> &f_stages);

    /**
     * This parameter is the factor (>1) by which the time step is
     * multiplied when the time stepping can be coarsen.
     */
    double coarsen_param;

    /**
     * This parameter is the factor (<1) by which the time step is
     * multiplied when the time stepping must be refined.
     */
    double refine_param;

    /**
     * Smallest time step allowed.
     */
    double min_delta_t;

    /**
     * Largest time step allowed.
     */
    double max_delta_t;

    /**
     * Refinement tolerance: if the error estimate is larger than
     * refine_tol, the time step is refined.
     */
    double refine_tol;

    /**
     * Coarsening tolerance: if the error estimate is smaller than
     * coarse_tol, the time step is coarsen.
     */
    double coarsen_tol;

    /**
     * If the flag is true, the last stage is the same as the first stage
     * and one evaluation of f can be saved.
     */
    bool last_same_as_first;

    /**
     * Butcher tableau coefficients.
     */
    std::vector<double> b1;

    /**
     * Butcher tableau coefficients.
     */
    std::vector<double> b2;

    /**
     * If the last_same_as_first flag is set to true, the last stage is
     * saved and reused as the first stage of the next time step.
     */
    VECTOR *last_stage;

    /**
     * Status structure of the object.
     */
    Status status;
  };
}

DEAL_II_NAMESPACE_CLOSE

#endif
