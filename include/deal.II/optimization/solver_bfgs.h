//-----------------------------------------------------------
//
//    Copyright (C) 2018 - 2020 by the deal.II authors
//
//    This file is part of the deal.II library.
//
//    The deal.II library is free software; you can use it, redistribute
//    it, and/or modify it under the terms of the GNU Lesser General
//    Public License as published by the Free Software Foundation; either
//    version 2.1 of the License, or (at your option) any later version.
//    The full text of the license can be found in the file LICENSE.md at
//    the top level directory of deal.II.
//
//---------------------------------------------------------------

#ifndef dealii_solver_bfgs_h
#define dealii_solver_bfgs_h

#include <deal.II/base/config.h>

#include <deal.II/lac/solver.h>

#include <deal.II/numerics/history.h>

#include <deal.II/optimization/line_minimization.h>

DEAL_II_NAMESPACE_OPEN

/**
 * Implement the limited memory BFGS minimization method.
 *
 * This class implements a method to minimize a given function for which only
 * the values of the function and its derivatives, but not its second
 * derivatives are available. The BFGS method is a variation of the Newton
 * method for function minimization in which the Hessian matrix is only
 * approximated. In particular, the Hessian is updated using the formula of
 * Broyden, Fletcher, Goldfarb, and Shanno (BFGS):
 * @f{align*}{
 * H^{(k+1)} &= \left[
 * I-\rho_{(k)} s^{(k)} \otimes y^{(k)}
 * \right]
 * H^{(k)}
 * \left[
 * I -\rho^{(k)} y^{(k)} \otimes s^{(k)}
 * \right]
 * +
 * \rho^{(k)} s^{(k)} \otimes s^{(k)}  \\
 * y^{(k)} &\dealcoloneq g^{(k+1)} - g^{(k)} \\
 * s^{(k)} &\dealcoloneq x^{(k+1)} - x^{(k)} \\
 * \rho^{(k)} &\dealcoloneq \frac{1}{y^{(k)} \cdot s^{(k)}}
 * @f}
 * for a symmetric positive definite $H$. Limited memory variant is
 * implemented via the two-loop recursion.
 *
 * @author Denis Davydov, 2018
 */
template <typename VectorType>
class SolverBFGS : public SolverBase<VectorType>
{
public:
  /**
   * Number type.
   */
  using Number = typename VectorType::value_type;


  /**
   * Standardized data struct to pipe additional data to the solver.
   */
  struct AdditionalData
  {
    /**
     * Constructor.
     */
    explicit AdditionalData(const unsigned int max_history_size = 5,
                            const bool         debug_output     = false);

    /**
     * Maximum history size.
     */
    unsigned int max_history_size;

    /**
     * Print extra debug output to deallog.
     */
    bool debug_output;
  };


  /**
   * Constructor.
   */
  explicit SolverBFGS(SolverControl &       residual_control,
                      const AdditionalData &data = AdditionalData());

  /**
   * Solve the unconstrained minimization problem
   * \f[
   * \min_{\mathbf x} f(\mathbf x)
   * \f]
   * starting from initial state @p x.
   *
   * The function @p compute takes two arguments indicating the values of $x$
   * and of the gradient $g=\nabla f(\mathbf x)=\frac{\partial f}{\partial
   * \mathbf x}$. When called, it needs to update the gradient $g$ at the given
   * location $x$ and return the value of the function being minimized, i.e.,
   * $f(\mathbf x)$.
   */
  void
  solve(
    const std::function<Number(const VectorType &x, VectorType &g)> &compute,
    VectorType &                                                     x);

  /**
   * Connect a slot to perform a custom line-search.
   *
   * Given the value of function @p f, the current value of unknown @p x,
   * the gradient @p g and the search direction @p p,
   * return the size $\alpha$ of the step $x \leftarrow x + \alpha p$,
   * and update @p x, @p g and @p f accordingly.
   */
  boost::signals2::connection
  connect_line_search_slot(
    const std::function<
      Number(Number &f, VectorType &x, VectorType &g, const VectorType &p)>
      &slot);

  /**
   * Connect a slot to perform a custom preconditioning.
   *
   * The preconditioner is applied inside the two loop recursion to
   * vector `g` using the history of position increments `s` and
   * gradient increments `y`.
   *
   * One possibility is to use the oldest `s,y` pair:
   * @code
   *  const auto preconditioner = [](VectorType &                         g,
   *                                 const FiniteSizeHistory<VectorType> &s,
   *                                 const FiniteSizeHistory<VectorType> &y) {
   *    if (s.size() > 0)
   *      {
   *        const unsigned int i  = s.size() - 1;
   *        const auto         yy = y[i] * y[i];
   *        const auto         sy = s[i] * y[i];
   *        Assert(yy > 0 && sy > 0, ExcInternalError());
   *        g *= sy / yy;
   *      }
   *  };
   * @endcode
   *
   * No preconditioning is performed if the code using this class has not
   * attached anything to the signal.
   */
  boost::signals2::connection
  connect_preconditioner_slot(
    const std::function<void(VectorType &                         g,
                             const FiniteSizeHistory<VectorType> &s,
                             const FiniteSizeHistory<VectorType> &y)> &slot);


protected:
  /**
   * Additional data to the solver.
   */
  const AdditionalData additional_data;

  /**
   * Signal used to perform line search.
   */
  boost::signals2::signal<
    Number(Number &f, VectorType &x, VectorType &g, const VectorType &p)>
    line_search_signal;

  /**
   * Signal used to perform preconditioning.
   */
  boost::signals2::signal<void(VectorType &                         g,
                               const FiniteSizeHistory<VectorType> &s,
                               const FiniteSizeHistory<VectorType> &y)>
    preconditioner_signal;
};


// -------------------  inline and template functions ----------------
#ifndef DOXYGEN

template <typename VectorType>
SolverBFGS<VectorType>::AdditionalData::AdditionalData(
  const unsigned int max_history_size_,
  const bool         debug_output_)
  : max_history_size(max_history_size_)
  , debug_output(debug_output_)
{}



template <typename VectorType>
SolverBFGS<VectorType>::SolverBFGS(SolverControl &       solver_control,
                                   const AdditionalData &data)
  : SolverBase<VectorType>(solver_control)
  , additional_data(data)
{}



template <class VectorType>
boost::signals2::connection
SolverBFGS<VectorType>::connect_line_search_slot(
  const std::function<
    Number(Number &f, VectorType &x, VectorType &g, const VectorType &p)> &slot)
{
  Assert(line_search_signal.empty(),
         ExcMessage("One should not attach more than one line search signal."));
  return line_search_signal.connect(slot);
}



template <class VectorType>
boost::signals2::connection
SolverBFGS<VectorType>::connect_preconditioner_slot(
  const std::function<void(VectorType &                         g,
                           const FiniteSizeHistory<VectorType> &s,
                           const FiniteSizeHistory<VectorType> &y)> &slot)
{
  Assert(preconditioner_signal.empty(),
         ExcMessage(
           "One should not attach more than one preconditioner signal."));
  return preconditioner_signal.connect(slot);
}



template <typename VectorType>
void
SolverBFGS<VectorType>::solve(
  const std::function<typename VectorType::value_type(const VectorType &x,
                                                      VectorType &f)> &compute,
  VectorType &                                                         x)
{
  // Also see scipy Fortran implementation
  // https://github.com/scipy/scipy/blob/master/scipy/optimize/lbfgsb_src/lbfgsb.f
  // and Octave-optim implementation:
  // https://sourceforge.net/p/octave/optim/ci/default/tree/src/__bfgsmin.cc
  LogStream::Prefix prefix("BFGS");

  // default line search:
  bool   first_step = true;
  Number f_prev     = 0.;
  // provide default line search if no signal was attached
  VectorType x0;
  if (line_search_signal.empty())
    {
      x0.reinit(x);
      const auto default_line_min =
        [&](Number &f, VectorType &x, VectorType &g, const VectorType &p) {
          const Number f0 = f;
          const Number g0 = g * p;
          Assert(g0 < 0,
                 ExcMessage(
                   "Function does not decrease along the current direction"));

          // save current solution value (to be used in line_search):
          x0 = x;

          // see scipy implementation
          // https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.line_search.html#scipy.optimize.line_search
          // and Eq. 2.6.8 in Fletcher 2013, Practical methods of optimization
          Number df = f_prev - f;
          Assert(first_step || df >= 0.,
                 ExcMessage("Function value is not decreasing"));
          df = std::max(df, 100. * std::numeric_limits<Number>::epsilon());
          // guess a reasonable first step:
          const Number a1 =
            (first_step ? 1. : std::min(1., -1.01 * 2. * df / g0));
          Assert(a1 > 0., ExcInternalError());
          f_prev = f;

          // 1D line-search function
          const auto line_func =
            [&](const Number &x_line) -> std::pair<Number, Number> {
            x = x0;
            x.add(x_line, p);
            f                   = compute(x, g);
            const Number g_line = g * p;
            return std::make_pair(f, g_line);
          };

          // loose line search:
          const auto res = LineMinimization::line_search<Number>(
            line_func,
            f0,
            g0,
            LineMinimization::poly_fit<Number>,
            a1,
            0.9,
            0.001);

          if (first_step)
            first_step = false;

          return res.first;
        };
      this->connect_line_search_slot(default_line_min);
    }

  // FIXME: Octave has convergence in terms of:
  // function change tolerance, default 1e-12
  // parameter change tolerance, default 1e-6
  // gradient tolerance, default 1e-5
  // SolverBase and/or SolverControl need extension

  VectorType g(x), p(x), y_k(x), s_k(x);

  std::vector<Number> c1;
  c1.reserve(additional_data.max_history_size);

  // limited history
  FiniteSizeHistory<VectorType> y(additional_data.max_history_size);
  FiniteSizeHistory<VectorType> s(additional_data.max_history_size);
  FiniteSizeHistory<Number>     rho(additional_data.max_history_size);

  unsigned int m = 0;
  Number       f;

  SolverControl::State conv = SolverControl::iterate;
  unsigned int         k    = 0;

  f = compute(x, g);

  conv = this->iteration_status(k, g.l2_norm(), x);
  if (conv != SolverControl::iterate)
    return;

  while (conv == SolverControl::iterate)
    {
      if (additional_data.debug_output)
        deallog << "Iteration " << k << " history " << m << std::endl
                << "f=" << f << std::endl;

      // 1. Two loop recursion to calculate p = - H*g
      c1.resize(m);
      p = g;
      // first loop:
      for (unsigned int i = 0; i < m; ++i)
        {
          c1[i] = rho[i] * (s[i] * p);
          p.add(-c1[i], y[i]);
        }
      // H0
      if (!preconditioner_signal.empty())
        preconditioner_signal(p, s, y);

      // second loop:
      for (int i = m - 1; i >= 0; --i)
        {
          Assert(i >= 0, ExcInternalError());
          const Number c2 = rho[i] * (y[i] * p);
          p.add(c1[i] - c2, s[i]);
        }
      p *= -1.;

      // 2. Line search
      s_k                = x;
      y_k                = g;
      const Number alpha = line_search_signal(f, x, g, p)
                             .get(); // <-- signals return boost::optional
      s_k.sadd(-1, 1, x);
      y_k.sadd(-1, 1, g);

      if (additional_data.debug_output)
        deallog << "Line search a=" << alpha << " f=" << f << std::endl;

      // 3. Check convergence
      k++;
      const Number g_l2 = g.l2_norm();
      conv              = this->iteration_status(k, g_l2, x);
      if (conv != SolverControl::iterate)
        break;

      // 4. Store s, y, rho
      const Number curvature = s_k * y_k;
      if (additional_data.debug_output)
        deallog << "Curvature " << curvature << std::endl;

      if (curvature > 0. && additional_data.max_history_size > 0)
        {
          s.add(s_k);
          y.add(y_k);
          rho.add(1. / curvature);
          m = s.size();

          Assert(y.size() == m, ExcInternalError());
          Assert(rho.size() == m, ExcInternalError());
        }

      Assert(m <= additional_data.max_history_size, ExcInternalError());
    }

  // In the case of failure: throw exception.
  AssertThrow(conv == SolverControl::success,
              SolverControl::NoConvergence(k, g.l2_norm()));
}

#endif

DEAL_II_NAMESPACE_CLOSE

#endif
