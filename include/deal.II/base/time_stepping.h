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

#ifndef dealii_time_stepping_h
#define dealii_time_stepping_h


#include <deal.II/base/config.h>

#include <deal.II/base/signaling_nan.h>

#include <functional>
#include <vector>

DEAL_II_NAMESPACE_OPEN

/**
 * Namespace containing the time stepping methods.
 */

namespace TimeStepping
{
  /**
   * The following Runge-Kutta methods are available:
   * - Explicit methods (see ExplicitRungeKutta::initialize):
   *   - FORWARD_EULER (first order)
   *   - RK_THIRD_ORDER (third order Runge-Kutta)
   *   - SSP_THIRD_ORDER (third order SSP Runge-Kutta)
   *   - RK_CLASSIC_FOURTH_ORDER (classical fourth order Runge-Kutta)
   * - Low-storage (explicit) Runge-Kutta methods
   *   - LOW_STORAGE_RK_STAGE3_ORDER3 (Three stages and third order)
   *   - LOW_STORAGE_RK_STAGE5_ORDER4 (Five stages and fourth order)
   *   - LOW_STORAGE_RK_STAGE7_ORDER4 (Seven stages and fourth order)
   *   - LOW_STORAGE_RK_STAGE9_ORDER5 (Nine stages and fifth order)
   * - Implicit methods (see ImplicitRungeKutta::initialize):
   *   - BACKWARD_EULER (first order)
   *   - IMPLICIT_MIDPOINT (second order)
   *   - CRANK_NICOLSON (second order)
   *   - SDIRK_TWO_STAGES (second order)
   * - Embedded explicit methods (see EmbeddedExplicitRungeKutta::initialize):
   *   - HEUN_EULER (second order)
   *   - BOGACKI_SHAMPINE (third order)
   *   - DOPRI (Dormand-Prince method, fifth order; this is the method used by
   * ode45 in MATLAB)
   *   - FEHLBERG (fifth order)
   *   - CASH_KARP (fifth order)
   */
  enum runge_kutta_method
  {
    /**
     * Forward Euler method, first order.
     */
    FORWARD_EULER,
    /**
     * Third order Runge-Kutta method.
     */
    RK_THIRD_ORDER,
    /**
     * Third order Strong Stability Preserving (SSP) Runge-Kutta method
     * (SSP time discretizations are also called Total Variation Diminishing
     * (TVD) methods in the literature, see @cite gottlieb2001strong).
     */
    SSP_THIRD_ORDER,
    /**
     * Classical fourth order Runge-Kutta method.
     */
    RK_CLASSIC_FOURTH_ORDER,
    /**
     * Fifth order Runge-Kutta method.
     */
    RK_FIFTH_ORDER,
    /**
     * Sixth order Runge-Kutta method.
     */
    RK_SIXTH_ORDER,
    /**
     * Three-stage scheme of order three by Kennedy et al.
     * @cite KennedyCarpenterLewis2000. Its stability region is
     * significantly smaller than the higher order schemes, but due to three
     * stages only, it is very competitive in terms of the work per stage.
     */
    LOW_STORAGE_RK_STAGE3_ORDER3,
    /**
     * Five-stage scheme of order four,
     * defined in the paper by Kennedy et al. @cite KennedyCarpenterLewis2000.
     */
    LOW_STORAGE_RK_STAGE5_ORDER4,
    /**
     * Seven-stage scheme of order four defined in the paper by Tselios and
     * Simos @cite TseliosSimos2007.
     */
    LOW_STORAGE_RK_STAGE7_ORDER4,
    /**
     * Nine-stage scheme of order five
     * defined in the paper by Kennedy et al. @cite KennedyCarpenterLewis2000.
     */
    LOW_STORAGE_RK_STAGE9_ORDER5,
    /**
     * Backward Euler method, first order.
     */
    BACKWARD_EULER,
    /**
     * Implicit midpoint method, second order.
     */
    IMPLICIT_MIDPOINT,
    /**
     * Crank-Nicolson method, second order.
     */
    CRANK_NICOLSON,
    /**
     * Two stage SDIRK method (short for "singly diagonally implicit
     * Runge-Kutta"), second order.
     */
    SDIRK_TWO_STAGES,
    /**
     * Heun's method (improved Euler's method), second order.
     */
    HEUN_EULER,
    /**
     * Bogacki–Shampine method, third-order.
     */
    BOGACKI_SHAMPINE,
    /**
     * Dormand-Prince method, fifth order; this is the method used by
     * ode45 in MATLAB.
     */
    DOPRI,
    /**
     * Fehlberg method, fifth order.
     */
    FEHLBERG,
    /**
     * Cash–Karp method, fifth order.
     */
    CASH_KARP,
    /**
     * Invalid.
     */
    invalid
  };



  /**
   * Reason for exiting evolve_one_time_step when using an embedded method:
   * DELTA_T, MIN_DELTA_T, MAX_DELTA_T.
   */
  enum embedded_runge_kutta_time_step
  {
    /**
     * The time step is in the valid range.
     */
    DELTA_T,
    /**
     * The time step was increased to the minimum acceptable time step.
     */
    MIN_DELTA_T,
    /**
     * The time step was reduced to the maximum acceptable time step.
     */
    MAX_DELTA_T
  };



  /**
   * Abstract class for time stepping methods. These methods assume that the
   * equation has the form: $ \frac{\partial y}{\partial t} = f(t,y) $.
   */
  template <typename VectorType>
  class TimeStepping
  {
  public:
    /**
     * Virtual destructor.
     */
    virtual ~TimeStepping() = default;

    /**
     * Purely virtual function. This function is used to advance from time @p
     * t to t+ @p delta_t. @p F is a vector of functions $ f(t,y) $ that
     * should be integrated, the input parameters are the time t and the
     * vector y and the output is value of f at this point. @p J_inverse is a
     * vector functions that compute the inverse of the Jacobians associated
     * to the implicit problems. The input parameters are the time, $ \tau $,
     * and a vector. The output is the value of function at this point. This
     * function returns the time at the end of the time step.
     */
    virtual double
    evolve_one_time_step(
      std::vector<std::function<VectorType(const double, const VectorType &)>>
                                                                     &F,
      std::vector<std::function<
        VectorType(const double, const double, const VectorType &)>> &J_inverse,
      double                                                          t,
      double                                                          delta_t,
      VectorType                                                     &y) = 0;

    /**
     * Empty structure used to store information.
     */
    struct Status
    {};

    /**
     * Purely virtual function that return Status.
     */
    virtual const Status &
    get_status() const = 0;
  };



  /**
   * Base class for the Runge-Kutta method
   */
  template <typename VectorType>
  class RungeKutta : public TimeStepping<VectorType>
  {
  public:
    /**
     * Virtual destructor.
     */
    virtual ~RungeKutta() override = default;

    /**
     * Purely virtual method used to initialize the Runge-Kutta method.
     */
    virtual void
    initialize(const runge_kutta_method method) = 0;

    /**
     * This function is used to advance from time @p t to t+ @p delta_t. @p F
     * is a vector of functions $ f(t,y) $ that should be integrated, the
     * input parameters are the time t and the vector y and the output is
     * value of f at this point. @p J_inverse is a vector functions that
     * compute the inverse of the Jacobians associated to the implicit
     * problems. The input parameters are the time, $ \tau $, and a vector.
     * The output is the value of function at this point. This function
     * returns the time at the end of the time step. When using Runge-Kutta
     * methods, @p F and @p J_inverse can only contain one element.
     */
    double
    evolve_one_time_step(
      std::vector<std::function<VectorType(const double, const VectorType &)>>
                                                                     &F,
      std::vector<std::function<
        VectorType(const double, const double, const VectorType &)>> &J_inverse,
      double                                                          t,
      double                                                          delta_t,
      VectorType &y) override;

    /**
     * Purely virtual function. This function is used to advance from time @p
     * t to t+ @p delta_t. @p f  is the function $ f(t,y) $ that should be
     * integrated, the input parameters are the time t and the vector y and
     * the output is value of f at this point. @p id_minus_tau_J_inverse is a
     * function that computes $ inv(I-\tau J)$ where $ I $ is the identity
     * matrix, $ \tau $ is given, and $ J $ is the Jacobian $ \frac{\partial
     * f}{\partial y} $. The input parameters are the time, $ \tau $, and a
     * vector. The output is the value of function at this point.
     * evolve_one_time_step returns the time at the end of the time step.
     */
    virtual double
    evolve_one_time_step(
      const std::function<VectorType(const double, const VectorType &)> &f,
      const std::function<
        VectorType(const double, const double, const VectorType &)>
                 &id_minus_tau_J_inverse,
      double      t,
      double      delta_t,
      VectorType &y) = 0;

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
    std::vector<std::vector<double>> a;
  };



  /**
   * ExplicitRungeKutta is derived from RungeKutta and implement the explicit
   * methods.
   */
  template <typename VectorType>
  class ExplicitRungeKutta : public RungeKutta<VectorType>
  {
  public:
    using RungeKutta<VectorType>::evolve_one_time_step;

    /**
     * Default constructor. This constructor creates an object for which
     * you will want to call <code>initialize(runge_kutta_method)</code>
     * before it can be used.
     */
    ExplicitRungeKutta() = default;

    /**
     * Constructor. This function calls initialize(runge_kutta_method).
     */
    ExplicitRungeKutta(const runge_kutta_method method);

    /**
     * Initialize the explicit Runge-Kutta method.
     */
    void
    initialize(const runge_kutta_method method) override;

    /**
     * This function is used to advance from time @p t to t+ @p delta_t. @p f
     * is the function $ f(t,y) $ that should be integrated, the input
     * parameters are the time t and the vector y and the output is value of f
     * at this point. @p id_minus_tau_J_inverse is a function that computes $
     * inv(I-\tau J)$ where $ I $ is the identity matrix, $ \tau $ is given,
     * and $ J $ is the Jacobian $ \frac{\partial f}{\partial y} $. The input
     * parameter are the time, $ \tau $, and a vector. The output is the value
     * of function at this point. evolve_one_time_step returns the time at the
     * end of the time step.
     *
     * @note @p id_minus_tau_J_inverse is ignored since the method is explicit.
     */
    double
    evolve_one_time_step(
      const std::function<VectorType(const double, const VectorType &)> &f,
      const std::function<
        VectorType(const double, const double, const VectorType &)>
                 &id_minus_tau_J_inverse,
      double      t,
      double      delta_t,
      VectorType &y) override;

    /**
     * This function is used to advance from time @p t to t+ @p delta_t. This
     * function is similar to the one derived from RungeKutta, but does not
     * required id_minus_tau_J_inverse because it is not used for explicit
     * methods. evolve_one_time_step returns the time at the end of the time
     * step.
     */
    double
    evolve_one_time_step(
      const std::function<VectorType(const double, const VectorType &)> &f,
      double                                                             t,
      double      delta_t,
      VectorType &y);

    /**
     * This structure stores the name of the method used.
     */
    struct Status : public TimeStepping<VectorType>::Status
    {
      Status()
        : method(invalid)
      {}

      runge_kutta_method method;
    };

    /**
     * Return the status of the current object.
     */
    const Status &
    get_status() const override;

  private:
    /**
     * Compute the different stages needed.
     */
    void
    compute_stages(
      const std::function<VectorType(const double, const VectorType &)> &f,
      const double                                                       t,
      const double             delta_t,
      const VectorType        &y,
      std::vector<VectorType> &f_stages) const;

    /**
     * Status structure of the object.
     */
    Status status;
  };



  /**
   * The LowStorageRungeKutta class is derived from RungeKutta and implements a
   * specific class of explicit methods. The main advantages of low-storage
   * methods are the reduced memory consumption and the reduced memory access.
   */
  template <typename VectorType>
  class LowStorageRungeKutta : public RungeKutta<VectorType>
  {
  public:
    using RungeKutta<VectorType>::evolve_one_time_step;

    /**
     * Default constructor. This constructor creates an object for which
     * you will want to call <code>initialize(runge_kutta_method)</code>
     * before it can be used.
     */
    LowStorageRungeKutta() = default;

    /**
     * Constructor. This function calls initialize(runge_kutta_method).
     */
    LowStorageRungeKutta(const runge_kutta_method method);

    /**
     * Initialize the explicit Runge-Kutta method.
     */
    void
    initialize(const runge_kutta_method method) override;

    /**
     * This function is used to advance from time @p t to t+ @p delta_t. @p f
     * is the function $ f(t,y) $ that should be integrated, the input
     * parameters are the time t and the vector y and the output is value of f
     * at this point. @p id_minus_tau_J_inverse is a function that computes $
     * inv(I-\tau J)$ where $ I $ is the identity matrix, $ \tau $ is given,
     * and $ J $ is the Jacobian $ \frac{\partial f}{\partial y} $. The input
     * parameters are the time, $ \tau $, and a vector. The output is the value
     * of function at this point. evolve_one_time_step returns the time at the
     * end of the time step.
     */
    double
    evolve_one_time_step(
      const std::function<VectorType(const double, const VectorType &)> &f,
      const std::function<
        VectorType(const double, const double, const VectorType &)>
                 &id_minus_tau_J_inverse,
      double      t,
      double      delta_t,
      VectorType &y) override;

    /**
     * This function is used to advance from time @p t to t+ @p delta_t. This
     * function is similar to the one derived from RungeKutta, but does not
     * required id_minus_tau_J_inverse because it is not used for explicit
     * methods. evolve_one_time_step returns the time at the end of the time
     * step. Note that vec_ki holds the evaluation of the differential operator,
     * and vec_ri holds the right-hand side for the differential operator
     * application.
     */
    double
    evolve_one_time_step(
      const std::function<VectorType(const double, const VectorType &)> &f,
      double                                                             t,
      double      delta_t,
      VectorType &solution,
      VectorType &vec_ri,
      VectorType &vec_ki);

    /**
     * Get the coefficients of the scheme.
     * Note that here vector @p a is not the conventional definition in terms of a
     * Butcher tableau but merely one of the sub-diagonals. More details can be
     * found in step-67 and the references therein.
     */
    void
    get_coefficients(std::vector<double> &a,
                     std::vector<double> &b,
                     std::vector<double> &c) const;

    /**
     * This structure stores the name of the method used.
     */
    struct Status : public TimeStepping<VectorType>::Status
    {
      Status()
        : method(invalid)
      {}

      runge_kutta_method method;
    };

    /**
     * Return the status of the current object.
     */
    const Status &
    get_status() const override;

  private:
    /**
     * Compute  one stage of low storage rk.
     */
    void
    compute_one_stage(
      const std::function<VectorType(const double, const VectorType &)> &f,
      const double                                                       t,
      const double      factor_solution,
      const double      factor_ai,
      const VectorType &current_ri,
      VectorType       &vec_ki,
      VectorType       &solution,
      VectorType       &next_ri) const;

    /**
     * Status structure of the object.
     */
    Status status;
  };



  /**
   * This class is derived from RungeKutta and implement the implicit methods.
   * This class works only for Diagonal Implicit Runge-Kutta (DIRK) methods.
   */
  template <typename VectorType>
  class ImplicitRungeKutta : public RungeKutta<VectorType>
  {
  public:
    using RungeKutta<VectorType>::evolve_one_time_step;

    /**
     * Default constructor. initialize(runge_kutta_method) and
     * set_newton_solver_parameters(unsigned int,double) need to be called
     * before the object can be used.
     */
    ImplicitRungeKutta() = default;

    /**
     * Constructor. This function calls initialize(runge_kutta_method) and
     * initialize the maximum number of iterations and the tolerance of the
     * Newton solver.
     */
    ImplicitRungeKutta(const runge_kutta_method method,
                       const unsigned int       max_it    = 100,
                       const double             tolerance = 1e-6);

    /**
     * Initialize the implicit Runge-Kutta method.
     */
    void
    initialize(const runge_kutta_method method) override;

    /**
     * This function is used to advance from time @p t to t+ @p delta_t. @p f
     * is the function $ f(t,y) $ that should be integrated, the input
     * parameters are the time t and the vector y and the output is value of f
     * at this point. @p id_minus_tau_J_inverse is a function that computes $
     * (I-\tau J)^{-1}$ where $ I $ is the identity matrix, $ \tau $ is given,
     * and $ J $ is the Jacobian $ \frac{\partial f}{\partial y} $. The input
     * parameters this function receives are the time, $ \tau $, and a vector.
     * The output is the value of function at this point. evolve_one_time_step
     * returns the time at the end of the time step.
     */
    double
    evolve_one_time_step(
      const std::function<VectorType(const double, const VectorType &)> &f,
      const std::function<
        VectorType(const double, const double, const VectorType &)>
                 &id_minus_tau_J_inverse,
      double      t,
      double      delta_t,
      VectorType &y) override;

    /**
     * Set the maximum number of iterations and the tolerance used by the
     * Newton solver.
     */
    void
    set_newton_solver_parameters(const unsigned int max_it,
                                 const double       tolerance);

    /**
     * Structure that stores the name of the method, the number of Newton
     * iterations and the norm of the residual when exiting the Newton solver.
     */
    struct Status : public TimeStepping<VectorType>::Status
    {
      Status()
        : method(invalid)
        , n_iterations(numbers::invalid_unsigned_int)
        , norm_residual(numbers::signaling_nan<double>())
      {}

      runge_kutta_method method;
      unsigned int       n_iterations;
      double             norm_residual;
    };

    /**
     * Return the status of the current object.
     */
    const Status &
    get_status() const override;

  private:
    /**
     * Compute the different stages needed.
     */
    void
    compute_stages(
      const std::function<VectorType(const double, const VectorType &)> &f,
      const std::function<
        VectorType(const double, const double, const VectorType &)>
                              &id_minus_tau_J_inverse,
      double                   t,
      double                   delta_t,
      VectorType              &y,
      std::vector<VectorType> &f_stages);

    /**
     * Newton solver used for the implicit stages.
     */
    void
    newton_solve(
      const std::function<void(const VectorType &, VectorType &)> &get_residual,
      const std::function<VectorType(const VectorType &)>
                 &id_minus_tau_J_inverse,
      VectorType &y);

    /**
     * Compute the residual needed by the Newton solver.
     */
    void
    compute_residual(
      const std::function<VectorType(const double, const VectorType &)> &f,
      double                                                             t,
      double            delta_t,
      const VectorType &new_y,
      const VectorType &y,
      VectorType       &tendency,
      VectorType       &residual) const;

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
   * This class is derived from RungeKutta and implements embedded explicit
   * methods.
   */
  template <typename VectorType>
  class EmbeddedExplicitRungeKutta : public RungeKutta<VectorType>
  {
  public:
    using RungeKutta<VectorType>::evolve_one_time_step;

    /**
     * Default constructor. initialize(runge_kutta_method) and
     * set_time_adaptation_parameters(double, double, double, double, double,
     * double) need to be called before the object can be used.
     */
    EmbeddedExplicitRungeKutta() = default;

    /**
     * Constructor. This function calls initialize(runge_kutta_method) and
     * initialize the parameters needed for time adaptation.
     */
    EmbeddedExplicitRungeKutta(const runge_kutta_method method,
                               const double             coarsen_param = 1.2,
                               const double             refine_param  = 0.8,
                               const double             min_delta     = 1e-14,
                               const double             max_delta     = 1e100,
                               const double             refine_tol    = 1e-8,
                               const double             coarsen_tol   = 1e-12);

    /**
     * Destructor.
     */
    ~EmbeddedExplicitRungeKutta() override
    {
      free_memory();
    }

    /**
     * If necessary, deallocate memory allocated by the object.
     */
    void
    free_memory();

    /**
     * Initialize the embedded explicit Runge-Kutta method.
     */
    void
    initialize(const runge_kutta_method method) override;

    /**
     * This function is used to advance from time @p t to t+ @p delta_t. @p f
     * is the function $ f(t,y) $ that should be integrated, the input
     * parameters are the time t and the vector y and the output is value of f
     * at this point. @p id_minus_tau_J_inverse is a function that computes $
     * inv(I-\tau J)$ where $ I $ is the identity matrix, $ \tau $ is given,
     * and $ J $ is the Jacobian $ \frac{\partial f}{\partial y} $. The input
     * parameters are the time, $ \tau $, and a vector. The output is the
     * value of function at this point. evolve_one_time_step returns the time
     * at the end of the time step.
     *
     * @note @p id_minus_tau_J_inverse is ignored since the method is explicit.
     */
    double
    evolve_one_time_step(
      const std::function<VectorType(const double, const VectorType &)> &f,
      const std::function<
        VectorType(const double, const double, const VectorType &)>
                 &id_minus_tau_J_inverse,
      double      t,
      double      delta_t,
      VectorType &y) override;

    /**
     * This function is used to advance from time @p t to t+ @p delta_t. This
     * function is similar to the one derived from TimeStepping, but does not
     * required id_minus_tau_J_inverse because it is not used for explicit
     * methods. evolve_one_time_step returns the time at the end of the time
     * step.
     */
    double
    evolve_one_time_step(
      const std::function<VectorType(const double, const VectorType &)> &f,
      double                                                             t,
      double      delta_t,
      VectorType &y);

    /**
     * Set the parameters necessary for the time adaptation.
     */
    void
    set_time_adaptation_parameters(const double coarsen_param,
                                   const double refine_param,
                                   const double min_delta,
                                   const double max_delta,
                                   const double refine_tol,
                                   const double coarsen_tol);

    /**
     * Structure that stores the name of the method, the reason to exit
     * evolve_one_time_step, the number of iteration inside n_iterations, a
     * guess of what the next time step should be, and an estimate of the norm
     * of the error.
     */
    struct Status : public TimeStepping<VectorType>::Status
    {
      Status()
        : method(invalid)
      {}

      runge_kutta_method             method;
      embedded_runge_kutta_time_step exit_delta_t;
      unsigned int                   n_iterations;
      double                         delta_t_guess;
      double                         error_norm;
    };

    /**
     * Return the status of the current object.
     */
    const Status &
    get_status() const override;

  private:
    /**
     * Compute the different stages needed.
     */
    void
    compute_stages(
      const std::function<VectorType(const double, const VectorType &)> &f,
      const double                                                       t,
      const double             delta_t,
      const VectorType        &y,
      std::vector<VectorType> &f_stages);

    /**
     * This parameter is the factor (>1) by which the time step is multiplied
     * when the time stepping can be coarsen.
     */
    double coarsen_param;

    /**
     * This parameter is the factor (<1) by which the time step is multiplied
     * when the time stepping must be refined.
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
     * Refinement tolerance: if the error estimate is larger than refine_tol,
     * the time step is refined.
     */
    double refine_tol;

    /**
     * Coarsening tolerance: if the error estimate is smaller than coarse_tol,
     * the time step is coarsen.
     */
    double coarsen_tol;

    /**
     * If the flag is true, the last stage is the same as the first stage and
     * one evaluation of f can be saved.
     */
    bool last_same_as_first = false;

    /**
     * Butcher tableau coefficients.
     */
    std::vector<double> b1;

    /**
     * Butcher tableau coefficients.
     */
    std::vector<double> b2;

    /**
     * If the last_same_as_first flag is set to true, the last stage is saved
     * and reused as the first stage of the next time step.
     */
    VectorType *last_stage = nullptr;

    /**
     * Status structure of the object.
     */
    Status status;
  };
} // namespace TimeStepping

DEAL_II_NAMESPACE_CLOSE

#endif
