//-----------------------------------------------------------
//
//    Copyright (C) 2018 by the deal.II authors
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

#ifndef dealii_line_minimization_h
#define dealii_line_minimization_h

#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/numbers.h>
#include <deal.II/base/utilities.h>

#include <deal.II/numerics/history.h>

#include <boost/optional.hpp>

#include <errno.h>
#include <sys/stat.h>

#include <fstream>
#include <string>

DEAL_II_NAMESPACE_OPEN

using namespace dealii;

/**
 * A namespace for various algorithms related to minimization over line.
 */
namespace LineMinimization
{
  /**
   * Given $x\_low$ and $x\_hi$ together with values of function
   * $f(x\_low)$ and $f(x\_hi)$ and the gradient $g(x\_low)$, return the local
   * minimizer of the quadratic interpolation function.
   *
   * The return type is optional to fit with similar function which may
   * not have solution for given parameters.
   */
  template <typename NumberType>
  boost::optional<NumberType>
  quadratic_fit(const NumberType x_low,
                const NumberType f_low,
                const NumberType g_low,
                const NumberType x_hi,
                const NumberType f_hi);

  /**
   * Given $x\_low$ and $x\_hi$ together with values of function
   * $f(x\_low)$ and $f(x\_hi)$) and its gradients ($g(x\_low)*g(x\_hi) < 0$) at
   * those points, return the local minimizer of the cubic interpolation
   * function. That is the location where the cubic interpolation function
   * attains its minimum.
   *
   * The return type is optional as the real-valued solution might not exist.
   */
  template <typename NumberType>
  boost::optional<NumberType>
  cubic_fit(const NumberType x_low,
            const NumberType f_low,
            const NumberType g_low,
            const NumberType x_hi,
            const NumberType f_hi,
            const NumberType g_hi);

  /**
   * Find the minimizer of a cubic polynomial that goes through the
   * points $f\_low=f(x\_low)$, $f\_hi=f(x\_hi)$ and $f\_rec(x\_rec)$
   * and has derivatve $g\_low$ at $x\_low$.
   *
   * The return type is optional as the real-valued solution might not exist.
   */
  template <typename NumberType>
  boost::optional<NumberType>
  cubic_fit_three_points(const NumberType x_low,
                         const NumberType f_low,
                         const NumberType g_low,
                         const NumberType x_hi,
                         const NumberType f_hi,
                         const NumberType x_rec,
                         const NumberType f_rec);

  /**
   * Return the minimizer of a polynomial using function values @p f_low @p f_hi @p f_rec[0] at three points
   * @p x_low @p x_hi @p x_rec[0] and derivatives at two points @p g_low and @p g_hi. The returned point
   * should be within the bounds @p bounds .
   *
   * This function will first try the cubic fit (see cubic_fit ). If it's
   * unsuccessfull or not
   * within the provided @p bounds the quadratic fit will be performed (see quadratic_fit ). The function will
   * fallback to bisection if quadratic fit fails.
   */
  template <typename NumberType>
  NumberType
  poly_fit(const NumberType                        x_low,
           const NumberType                        f_low,
           const NumberType                        g_low,
           const NumberType                        x_hi,
           const NumberType                        f_hi,
           const NumberType                        g_hi,
           const FiniteSizeHistory<NumberType> &   x_rec,
           const FiniteSizeHistory<NumberType> &   f_rec,
           const FiniteSizeHistory<NumberType> &   g_rec,
           const std::pair<NumberType, NumberType> bounds);

  /**
   * Same as above but doing cubic fit with three points (see
   * cubic_fit_three_points ).
   */
  template <typename NumberType>
  NumberType
  poly_fit_three_points(const NumberType                        x_low,
                        const NumberType                        f_low,
                        const NumberType                        g_low,
                        const NumberType                        x_hi,
                        const NumberType                        f_hi,
                        const NumberType                        g_hi,
                        const FiniteSizeHistory<NumberType> &   x_rec,
                        const FiniteSizeHistory<NumberType> &   f_rec,
                        const FiniteSizeHistory<NumberType> &   g_rec,
                        const std::pair<NumberType, NumberType> bounds);


  /**
   * Perform a line search in $(0,max]$ with strong Wolfe conditions
   * \f[
   * f(\alpha) \le f(0) + \alpha \mu f'(0) \\
   * |f'(\alpha)| \le \eta |f'(0)|
   * \f]
   * using one dimensional
   * functions @p func and a function @p interpolate to choose a new point
   * from the interval based on the function values and derivatives at its ends.
   * @p a1 is a trial estimate of the first step.
   * Interpolation can be done using poly_fit or poly_fit_three_points .
   *
   * The function implements Algorithms
   * 2.6.2 and 2.6.4 on pages 34-35 in Fletcher, 2013, Practical methods of
   * optimization. these are minor variations of the Algorithm 3.5 and 3.6 on
   * pages 60-61 in Nocedal and Wright, Numerical optimization.
   * It consists of a bracketing phase and a zoom phase, where @p interpolate is used.
   *
   * The function returns the step size and the number of times function @p func was called.
   *
   * @param func a one dimensional function which returns value and derivative
   * at the given point.
   * @param f0 function value the origin
   * @param g0 function derivative the origin
   * @param interpolate a function which determines how interpolation is done
   * during the zoom phase. It takes values and derivatives at the current
   * interval/bracket ($f\_low$, $f\_hi$) as well as up to 5 values and
   * derivatives at previous steps. The returned value is to be provided within
   * the given bounds.
   * @param a1 initial trial step for bracketing phase
   * @param eta a parameter in the second Wolfe condition (curvature condition)
   * @param mu a parameter in the first Wolfe condition (sufficient decrease)
   * @param a_max maximum allowed step size
   * @param max_evaluations maximum allowed number of function evaluations
   * @param debug_output a flag to do extra debug output into deallog static
   * object
   */
  template <typename NumberType>
  std::pair<NumberType, unsigned int>
  line_search(
    const std::function<std::pair<NumberType, NumberType>(const NumberType x)>
      &              func,
    const NumberType f0,
    const NumberType g0,
    const std::function<
      NumberType(const NumberType                        x_low,
                 const NumberType                        f_low,
                 const NumberType                        g_low,
                 const NumberType                        x_hi,
                 const NumberType                        f_hi,
                 const NumberType                        g_hi,
                 const FiniteSizeHistory<NumberType> &   x_rec,
                 const FiniteSizeHistory<NumberType> &   f_rec,
                 const FiniteSizeHistory<NumberType> &   g_rec,
                 const std::pair<NumberType, NumberType> bounds)> &interpolate,
    const NumberType                                               a1,
    const NumberType                                               eta = 0.9,
    const NumberType                                               mu  = 0.01,
    const NumberType   a_max           = std::numeric_limits<NumberType>::max(),
    const unsigned int max_evaluations = 20,
    const bool         debug_output    = false);

  // -------------------  inline and template functions ----------------
#ifndef DOXYGEN

  template <typename NumberType>
  boost::optional<NumberType>
  quadratic_fit(const NumberType x1,
                const NumberType f1,
                const NumberType g1,
                const NumberType x2,
                const NumberType f2)
  {
    Assert(x1 != x2, ExcMessage("Point are the same"));
    const NumberType denom = (2. * g1 * x2 - 2. * g1 * x1 - 2. * f2 + 2. * f1);
    if (denom == 0)
      return boost::none;
    else
      return (g1 * (x2 * x2 - x1 * x1) + 2. * (f1 - f2) * x1) / denom;
  }

  template <typename NumberType>
  boost::optional<NumberType>
  cubic_fit(const NumberType x1,
            const NumberType f1,
            const NumberType g1,
            const NumberType x2,
            const NumberType f2,
            const NumberType g2)
  {
    Assert(x1 != x2, ExcMessage("Points are the same"));
    const NumberType beta1 = g1 + g2 - 3. * (f1 - f2) / (x1 - x2);
    const NumberType s     = beta1 * beta1 - g1 * g2;
    if (s < 0)
      return boost::none;

    const NumberType beta2 = std::sqrt(s);
    const NumberType denom =
      x1 < x2 ? g2 - g1 + 2. * beta2 : g1 - g2 + 2. * beta2;
    if (denom == 0.)
      return boost::none;

    return x1 < x2 ? x2 - (x2 - x1) * (g2 + beta2 - beta1) / denom :
                     x1 - (x1 - x2) * (g1 + beta2 - beta1) / denom;
  }



  template <typename NumberType>
  boost::optional<NumberType>
  cubic_fit_three_points(const NumberType x1,
                         const NumberType f1,
                         const NumberType g1,
                         const NumberType x2,
                         const NumberType f2,
                         const NumberType x3,
                         const NumberType f3)
  {
    Assert(x1 != x2, ExcMessage("Points are the same"));
    Assert(x1 != x3, ExcMessage("Points are the same"));
    // f(x) = A *(x-x1)^3 + B*(x-x1)^2 + C*(x-x1) + D
    // =>
    // D = f1
    // C = g1

    // the rest is a system of 2 equations:

    const NumberType x2_shift = x2 - x1;
    const NumberType x3_shift = x3 - x1;
    const NumberType r1       = f2 - f1 - g1 * x2_shift;
    const NumberType r2       = f3 - f1 - g1 * x3_shift;
    const NumberType denom =
      std::pow(x2_shift * x3_shift, 2) * (x2_shift - x3_shift);
    if (denom == 0.)
      return boost::none;

    const NumberType A =
      (r1 * std::pow(x3_shift, 2) - r2 * std::pow(x2_shift, 2)) / denom;
    const NumberType B =
      (r2 * std::pow(x2_shift, 3) - r1 * std::pow(x3_shift, 3)) / denom;
    const NumberType &C = g1;

    // now get the minimizer:
    const NumberType radical = B * B - A * C * 3;
    if (radical < 0)
      return boost::none;

    return x1 + (-B + std::sqrt(radical)) / (A * 3);
  }


  template <typename NumberType>
  NumberType
  poly_fit(const NumberType x1,
           const NumberType f1,
           const NumberType g1,
           const NumberType x2,
           const NumberType f2,
           const NumberType g2,
           const FiniteSizeHistory<NumberType> &,
           const FiniteSizeHistory<NumberType> &,
           const FiniteSizeHistory<NumberType> &,
           const std::pair<NumberType, NumberType> bounds)
  {
    Assert(bounds.first < bounds.second, ExcMessage("Incorrect bounds"));

    // Similar to scipy implementation but we fit based on two points
    // with their gradients and do bisection on bounds.
    // https://github.com/scipy/scipy/blob/v1.0.0/scipy/optimize/linesearch.py#L555-L563

    // First try cubic interpolation
    boost::optional<NumberType> res = cubic_fit(x1, f1, g1, x2, f2, g2);
    if (res && *res >= bounds.first && *res <= bounds.second)
      return *res;

    // cubic either fails or outside of safe region, do quadratic:
    res = quadratic_fit(x1, f1, g1, x2, f2);
    if (res && *res >= bounds.first && *res <= bounds.second)
      return *res;

    // quadratic either failed or outside of safe region. Do bisection
    // on safe region
    return (bounds.first + bounds.second) * 0.5;
  }



  template <typename NumberType>
  NumberType
  poly_fit_three_points(const NumberType                     x1,
                        const NumberType                     f1,
                        const NumberType                     g1,
                        const NumberType                     x2,
                        const NumberType                     f2,
                        const NumberType                     g2,
                        const FiniteSizeHistory<NumberType> &x_rec,
                        const FiniteSizeHistory<NumberType> &f_rec,
                        const FiniteSizeHistory<NumberType> & /*g_rec*/,
                        const std::pair<NumberType, NumberType> bounds)
  {
    Assert(bounds.first < bounds.second, ExcMessage("Incorrect bounds"));
    AssertDimension(x_rec.size(), f_rec.size());

    // Same as scipy implementation where cubic fit is using 3 points
    // https://github.com/scipy/scipy/blob/v1.0.0/scipy/optimize/linesearch.py#L555-L563

    // First try cubic interpolation after first iteration
    boost::optional<NumberType> res =
      x_rec.size() > 0 ?
        cubic_fit_three_points(x1, f1, g1, x2, f2, x_rec[0], f_rec[0]) :
        boost::none;
    if (res && *res >= bounds.first && *res <= bounds.second)
      return *res;

    // cubic either fails or outside of safe region, do quadratic:
    res = quadratic_fit(x1, f1, g1, x2, f2);
    if (res && *res >= bounds.first && *res <= bounds.second)
      return *res;

    // quadratic either failed or outside of safe region. Do bisection
    // on safe region
    return (bounds.first + bounds.second) * 0.5;
  }


  template <typename NumberType>
  std::pair<NumberType, unsigned int>
  line_search(
    const std::function<std::pair<NumberType, NumberType>(const NumberType x)>
      &              func,
    const NumberType f0,
    const NumberType g0,
    const std::function<
      NumberType(const NumberType                        x_low,
                 const NumberType                        f_low,
                 const NumberType                        g_low,
                 const NumberType                        x_hi,
                 const NumberType                        f_hi,
                 const NumberType                        g_hi,
                 const FiniteSizeHistory<NumberType> &   x_rec,
                 const FiniteSizeHistory<NumberType> &   f_rec,
                 const FiniteSizeHistory<NumberType> &   g_rec,
                 const std::pair<NumberType, NumberType> bounds)> &choose,
    const NumberType                                               a1,
    const NumberType                                               eta,
    const NumberType                                               mu,
    const NumberType                                               a_max,
    const unsigned int max_evaluations,
    const bool         debug_output)
  {
    // Note that scipy use dcsrch() from Minpack2 Fortran lib for line search
    Assert(mu < 0.5 && mu > 0, ExcMessage("mu is not in (0,1/2)."));
    Assert(eta < 1. && eta > mu, ExcMessage("eta is not in (mu,1)."));
    Assert(a_max > 0, ExcMessage("max is not positive."));
    Assert(a1 > 0 && a1 <= a_max, ExcMessage("a1 is not in (0,max]."));
    Assert(g0 < 0, ExcMessage("Initial slope is not negative"));

    // Growth parameter for bracketing phase:
    // 1 < tau1
    const NumberType tau1 = 9.;
    // shrink parameters for sectioning phase to prevent ai from being
    // arbitrary close to the extremes of the interval.
    // 0 < tau2 < tau3 <= 1/2
    // tau2 <= eta is advisable
    const NumberType tau2 = 0.1; // bound for closeness to a_lo
    const NumberType tau3 = 0.5; // bound for closeness to a_hi

    const NumberType g0_abs = std::abs(g0);
    const NumberType f_min  = f0 + a_max * mu * g0;

    // return True if the first Wolfe condition (sufficient decrease) is
    // satisfied
    const auto w1 = [&](const NumberType a, const NumberType f) {
      return f <= f0 + a * mu * g0;
    };

    // return True if the second Wolfe condition (curvature condition) is
    // satisfied
    const auto w2 = [&](const NumberType g) {
      return std::abs(g) <= eta * g0_abs;
    };

    // Bracketing phase (Algorithm 2.6.2): look for a non-trivial interval
    // which is known to contain an interval of acceptable points.
    // We adopt notation of Noceal.
    const NumberType x    = std::numeric_limits<NumberType>::signaling_NaN();
    NumberType       a_lo = x, f_lo = x, g_lo = x;
    NumberType       a_hi = x, f_hi = x, g_hi = x;
    NumberType       ai = x, fi = x, gi = x;

    // count function calls in i:
    unsigned int i = 0;
    {
      NumberType f_prev, g_prev, a_prev;
      ai     = a1;
      f_prev = f0;
      g_prev = g0;
      a_prev = 0;

      while (i < max_evaluations)
        {
          const auto fgi = func(ai);
          fi             = fgi.first;
          gi             = fgi.second;
          i++;

          if (debug_output)
            deallog << "Bracketing phase: " << i << std::endl
                    << ai << " " << fi << " " << gi << " " << w1(ai, fi) << " "
                    << w2(gi) << " " << f_min << std::endl;

          // first check if we can stop bracketing or the whole line search:
          if (fi <= f_min || ai == a_max)
            {
              if (debug_output)
                deallog << "Reached the maximum step size." << std::endl;
              return std::make_pair(ai, i);
            }

          if (!w1(ai, fi) ||
              (fi >= f_prev && i > 1)) // violate first Wolfe or not descending
            {
              a_lo = a_prev;
              f_lo = f_prev;
              g_lo = g_prev;

              a_hi = ai;
              f_hi = fi;
              g_hi = gi;
              break; // end bracketing
            }

          if (w2(gi)) // satisfies both Wolfe conditions
            {
              if (debug_output)
                deallog << "Satisfied both Wolfe conditions during Bracketing."
                        << std::endl;

              Assert(w1(ai, fi), ExcInternalError());
              return std::make_pair(ai, i);
            }

          if (gi >= 0) // not descending
            {
              a_lo = ai;
              f_lo = fi;
              g_lo = gi;

              a_hi = a_prev;
              f_hi = f_prev;
              g_hi = g_prev;
              break; // end bracketing
            }

          // extrapolation step with the bounds
          const auto bounds =
            std::make_pair(2. * ai - a_prev,
                           std::min(a_max, ai + tau1 * (ai - a_prev)));

          a_prev = ai;
          f_prev = fi;
          g_prev = gi;

          // NOTE: Fletcher's 2.6.2 includes optional extrapolation, we
          // simply take the upper bound
          // Scipy increases by factor of two:
          // https://github.com/scipy/scipy/blob/v1.0.0/scipy/optimize/linesearch.py#L447
          ai = bounds.second;
        }
    }

    AssertThrow(
      i < max_evaluations,
      ExcMessage(
        "Could not find the initial bracket within the given number of iterations."));

    // Check properties of the bracket (Theorem 3.2 in More and Thuente, 94
    // and Eq. 2.6.3 in Fletcher 2013

    // FIXME: these conditions are actually violated for Fig3 and a1=10^3 in
    // More and Thorenton, 94.

    /*
    Assert((f_lo < f_hi) && w1(a_lo, f_lo), ExcInternalError());
    Assert(((a_hi - a_lo) * g_lo < 0) && !w2(g_lo), ExcInternalError());
    Assert((w1(a_hi, f_hi) || f_hi >= f_lo), ExcInternalError());
    */

    // keep short history of last points to improve interpolation
    FiniteSizeHistory<NumberType> a_rec(5), f_rec(5), g_rec(5);
    // if neither a_lo nor a_hi are zero:
    if (std::abs(a_lo) > std::numeric_limits<NumberType>::epsilon() &&
        std::abs(a_hi) > std::numeric_limits<NumberType>::epsilon())
      {
        a_rec.add(0);
        f_rec.add(f0);
        g_rec.add(g0);
      }

    // Now sectioning phase: we allow both [a_lo, a_hi] and [a_hi, a_lo]
    while (i < max_evaluations)
      {
        const NumberType a_lo_safe = a_lo + tau2 * (a_hi - a_lo);
        const NumberType a_hi_safe = a_hi - tau3 * (a_hi - a_lo);
        const auto       bounds    = std::minmax(a_lo_safe, a_hi_safe);

        ai = choose(
          a_lo, f_lo, g_lo, a_hi, f_hi, g_hi, a_rec, f_rec, g_rec, bounds);

        const std::pair<NumberType, NumberType> fgi = func(ai);
        fi                                          = fgi.first;
        gi                                          = fgi.second;
        i++;

        if (debug_output)
          deallog << "Sectioning phase: " << i << std::endl
                  << a_lo << " " << f_lo << " " << g_lo << " " << w1(a_lo, f_lo)
                  << " " << w2(g_lo) << std::endl
                  << a_hi << " " << f_hi << " " << g_hi << " " << w1(a_hi, f_hi)
                  << " " << w2(g_hi) << std::endl
                  << ai << " " << fi << " " << gi << " " << w1(ai, fi) << " "
                  << w2(gi) << std::endl;

        if (!w1(ai, fi) || fi >= f_lo)
          // take [a_lo, ai]
          {
            a_rec.add(a_hi);
            f_rec.add(f_hi);
            g_rec.add(g_hi);

            a_hi = ai;
            f_hi = fi;
            g_hi = gi;
          }
        else
          {
            if (w2(gi)) // satisfies both wolf
              {
                if (debug_output)
                  deallog << "Satisfied both Wolfe conditions." << std::endl;
                Assert(w1(ai, fi), ExcInternalError());
                return std::make_pair(ai, i);
              }

            if (gi * (a_hi - a_lo) >= 0)
              // take [ai, a_lo]
              {
                a_rec.add(a_hi);
                f_rec.add(f_hi);
                g_rec.add(g_hi);

                a_hi = a_lo;
                f_hi = f_lo;
                g_hi = g_lo;
              }
            else
              // take [ai, a_hi]
              {
                a_rec.add(a_lo);
                f_rec.add(f_lo);
                g_rec.add(g_lo);
              }

            a_lo = ai;
            f_lo = fi;
            g_lo = gi;
          }
      }

    // if we got here, we could not find the solution
    AssertThrow(
      false,
      ExcMessage(
        "Could not could complete the sectioning phase within the given number of iterations."));
    return std::make_pair(std::numeric_limits<NumberType>::signaling_NaN(), i);
  }

#endif

} // namespace LineMinimization

DEAL_II_NAMESPACE_CLOSE

#endif // dealii_line_minimization_h
