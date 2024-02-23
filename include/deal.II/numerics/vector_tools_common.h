// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2020 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_vector_tools_common_h
#define dealii_vector_tools_common_h


#include <deal.II/base/config.h>

#include <deal.II/base/patterns.h>

DEAL_II_NAMESPACE_OPEN

namespace VectorTools
{
  /**
   * Denote which norm/integral is to be computed by the
   * integrate_difference() function on each cell and compute_global_error()
   * for the whole domain.
   * Let $f:\Omega \rightarrow \mathbb{R}^c$ be a finite element function
   * with $c$ components where component $c$ is denoted by $f_c$ and $\hat{f}$
   * be the reference function (the @p fe_function and @p exact_solution
   * arguments to integrate_difference()). Let $e_c = \hat{f}_c - f_c$
   * be the difference or error between the two. Further,
   * let  $w:\Omega \rightarrow \mathbb{R}^c$ be the @p weight function of integrate_difference(), which is
   * assumed to be equal to one if not supplied. Finally, let $p$ be the
   * @p exponent argument (for $L_p$-norms).
   *
   * In the following,we denote by $E_K$ the local error computed by
   * integrate_difference() on cell $K$, whereas $E$ is the global error
   * computed by compute_global_error(). Note that integrals are
   * approximated by quadrature in the usual way:
   * @f[
   * \int_A f(x) dx \approx \sum_q f(x_q) \omega_q.
   * @f]
   * Similarly for suprema over a cell $T$:
   * @f[
   * \sup_{x\in T} |f(x)| dx \approx \max_q |f(x_q)|.
   * @f]
   */
  enum NormType
  {
    /**
     * The function or difference of functions is integrated on each cell $K$:
     * @f[
     *   E_K
     * = \int_K \sum_c (\hat{f}_c - f_c) \, w_c
     * = \int_K \sum_c e_c \, w_c
     * @f]
     * and summed up to get
     * @f[
     *   E = \sum_K E_K
     *     = \int_\Omega \sum_c (\hat{f}_c - f_c) \, w_c
     * @f]
     * or, for $w \equiv 1$:
     * @f[
     *   E = \int_\Omega (\hat{f} - f)
     *     = \int_\Omega e.
     * @f]
     *
     * Note: This differs from what is typically known as
     * the mean of a function by a factor of $\frac{1}{|\Omega|}$. To
     * compute the mean you can also use compute_mean_value(). Finally,
     * pay attention to the sign: if $\hat{f}=0$, this will compute the
     * negative of the mean of $f$.
     */
    mean,

    /**
     * The absolute value of the function is integrated:
     * @f[
     *   E_K = \int_K \sum_c |e_c| \, w_c
     * @f]
     * and
     * @f[
     *   E = \sum_K E_K = \int_\Omega \sum_c |e_c| w_c,
     * @f]
     * or, for $w \equiv 1$:
     * @f[
     *   E  = \| e \|_{L^1}.
     * @f]
     */
    L1_norm,

    /**
     * The square of the function is integrated and the square root of the
     * result is computed on each cell:
     * @f[
     *   E_K = \sqrt{ \int_K \sum_c e_c^2 \, w_c }
     * @f]
     * and
     * @f[
     *   E = \sqrt{\sum_K E_K^2} = \sqrt{ \int_\Omega  \sum_c e_c^2 \, w_c }
     * @f]
     * or, for $w \equiv 1$:
     * @f[
     *   E = \sqrt{ \int_\Omega e^2 }
     *     = \| e \|_{L^2}
     * @f]
     */
    L2_norm,

    /**
     * The absolute value to the $p$-th power is integrated and the $p$-th
     * root is computed on each cell. The exponent $p$ is the @p
     * exponent argument of integrate_difference() and compute_mean_value():
     * @f[
     *   E_K = \left( \int_K \sum_c |e_c|^p \, w_c \right)^{1/p}
     * @f]
     * and
     * @f[
     *   E = \left( \sum_K E_K^p \right)^{1/p}
     * @f]
     * or, for $w \equiv 1$:
     * @f[
     *   E = \| e \|_{L^p}.
     * @f]
     */
    Lp_norm,

    /**
     * The maximum absolute value of the function:
     * @f[
     *   E_K = \sup_K \max_c |e_c| \, w_c
     * @f]
     * and
     * @f[
     *   E = \max_K E_K
     * = \sup_\Omega \max_c |e_c| \, w_c
     * @f]
     * or, for $w \equiv 1$:
     * @f[
     *   E  = \sup_\Omega \|e\|_\infty = \| e \|_{L^\infty}.
     * @f]
     */
    Linfty_norm,

    /**
     * #L2_norm of the gradient:
     * @f[
     *   E_K = \sqrt{ \int_K \sum_c (\nabla e_c)^2 \, w_c }
     * @f]
     * and
     * @f[
     *   E = \sqrt{\sum_K E_K^2} = \sqrt{ \int_\Omega \sum_c (\nabla e_c)^2 \,
     * w_c }
     * @f]
     * or, for $w \equiv 1$:
     * @f[
     *   E = \| \nabla e \|_{L^2}.
     * @f]
     */
    H1_seminorm,

    /**
     * #L2_norm of the divergence of a vector field. The function $f$ is
     * expected to have $c \geq \text{dim}$ components and the first @p dim
     * will be used to compute the divergence:
     * @f[
     *   E_K = \sqrt{ \int_K \left( \sum_c \frac{\partial e_c}{\partial x_c} \,
     * \sqrt{w_c} \right)^2 }
     * @f]
     * and
     * @f[
     *   E = \sqrt{\sum_K E_K^2}
     *     = \sqrt{ \int_\Omega \left( \sum_c \frac{\partial e_c}{\partial x_c}
     * \, \sqrt{w_c} \right)^2  }
     * @f]
     * or, for $w \equiv 1$:
     * @f[
     *   E = \| \nabla \cdot e \|_{L^2}.
     * @f]
     */
    Hdiv_seminorm,

    /**
     * The square of this norm is the square of the #L2_norm plus the square
     * of the #H1_seminorm:
     * @f[
     *   E_K = \sqrt{ \int_K \sum_c (e_c^2 + (\nabla e_c)^2) \, w_c }
     * @f]
     * and
     * @f[
     *   E = \sqrt{\sum_K E_K^2} = \sqrt{ \int_\Omega \sum_c (e_c^2 + (\nabla
     * e_c)^2) \, w_c }
     * @f]
     * or, for $w \equiv 1$:
     * @f[
     *   E = \left( \| e \|_{L^2}^2 + \| \nabla e \|_{L^2}^2 \right)^{1/2}.
     * @f]
     */
    H1_norm,

    /**
     * #Lp_norm of the gradient:
     * @f[
     *   E_K = \left( \int_K \sum_c |\nabla e_c|^p \, w_c \right)^{1/p}
     * @f]
     * and
     * @f[
     *   E = \left( \sum_K E_K^p \right)^{1/p}
     *     = \left( \int_\Omega \sum_c |\nabla e_c|^p \, w_c \right)^{1/p}
     * @f]
     * or, for $w \equiv 1$:
     * @f[
     *   E = \| \nabla e \|_{L^p}.
     * @f]
     */
    W1p_seminorm,

    /**
     * The same as the #H1_norm but using <i>L<sup>p</sup></i>:
     * @f[
     *   E_K = \left( \int_K \sum_c (|e_c|^p + |\nabla e_c|^p) \, w_c
     * \right)^{1/p}
     * @f]
     * and
     * @f[
     *   E = \left( \sum_K E_K^p \right)^{1/p}
     *     = \left( \int_\Omega \sum_c (|e_c|^p + |\nabla e_c|^p) \, w_c
     * \right)^{1/p}
     * @f]
     * or, for $w \equiv 1$:
     * @f[
     *   E = \left( \| e \|_{L^p}^p + \| \nabla e \|_{L^p}^p \right)^{1/p}.
     * @f]
     */
    W1p_norm,

    /**
     * #Linfty_norm of the gradient:
     * @f[
     *   E_K = \sup_K \max_c |\nabla e_c| \, w_c
     * @f]
     * and
     * @f[
     *   E = \max_K E_K
     *     = \sup_\Omega \max_c |\nabla e_c| \, w_c
     * @f]
     * or, for $w \equiv 1$:
     * @f[
     *   E = \| \nabla e \|_{L^\infty}.
     * @f]
     */
    W1infty_seminorm,

    /**
     * The sum of #Linfty_norm and #W1infty_seminorm:
     * @f[
     *   E_K = \sup_K \max_c |e_c| \, w_c + \sup_K \max_c |\nabla e_c| \, w_c.
     * @f]
     * The global norm is not implemented in compute_global_error(),
     * because it is impossible to compute the sum of the global
     * norms from the values $E_K$. As a work-around, you can compute the
     * global #Linfty_norm and #W1infty_seminorm separately and then add them
     * to get (with $w \equiv 1$):
     * @f[
     *   E = \| e \|_{L^\infty} + \| \nabla e \|_{L^\infty}.
     * @f]
     */
    W1infty_norm
  };

  /**
   * Exception
   */
  DeclExceptionMsg(ExcPointNotAvailableHere,
                   "The given point is inside a cell of a "
                   "parallel::distributed::Triangulation that is not "
                   "locally owned by this processor.");
} // namespace VectorTools

// Make sure we can use NormType with Patterns.
namespace Patterns
{
  namespace Tools
  {
    template <>
    struct Convert<VectorTools::NormType, void>
    {
      /**
       * Return the Correct pattern for NormType.
       */
      static std::unique_ptr<Patterns::PatternBase>
      to_pattern()
      {
        return std::make_unique<Patterns::Selection>(
          "mean|L1_norm|L2_norm|Lp_norm|"
          "Linfty_norm|H1_seminorm|Hdiv_seminorm|"
          "H1_norm|W1p_seminorm|W1p_norm|"
          "W1infty_seminorm|W1infty_norm");
      }



      /**
       * Convert a NormType to a string.
       */
      static std::string
      to_string(const VectorTools::NormType &s,
                const Patterns::PatternBase &p =
                  *Convert<VectorTools::NormType>::to_pattern())
      {
        std::string str;
        if (s == VectorTools::mean)
          str = "mean";
        else if (s == VectorTools::L1_norm)
          str = "L1_norm";
        else if (s == VectorTools::L2_norm)
          str = "L2_norm";
        else if (s == VectorTools::Lp_norm)
          str = "Lp_norm";
        else if (s == VectorTools::Linfty_norm)
          str = "Linfty_norm";
        else if (s == VectorTools::H1_seminorm)
          str = "H1_seminorm";
        else if (s == VectorTools::Hdiv_seminorm)
          str = "Hdiv_seminorm";
        else if (s == VectorTools::H1_norm)
          str = "H1_norm";
        else if (s == VectorTools::W1p_seminorm)
          str = "W1p_seminorm";
        else if (s == VectorTools::W1infty_seminorm)
          str = "W1infty_seminorm";
        else if (s == VectorTools::W1infty_norm)
          str = "W1infty_norm";
        else if (s == VectorTools::W1p_norm)
          str = "W1p_norm";
        else
          {
            AssertThrow(false, ExcMessage("Didn't recognize a norm type."));
          }
        AssertThrow(p.match(str), ExcInternalError());
        return str;
      }


      /**
       * Convert a string to a NormType.
       */
      static VectorTools::NormType
      to_value(const std::string           &str,
               const Patterns::PatternBase &p =
                 *Convert<VectorTools::NormType>::to_pattern())
      {
        VectorTools::NormType norm = VectorTools::mean;
        AssertThrow(p.match(str),
                    ExcMessage(
                      "String " + str +
                      " cannot be converted to VectorTools::NormType"));

        if (str == "mean")
          norm = VectorTools::mean;
        else if (str == "L1_norm")
          norm = VectorTools::L1_norm;
        else if (str == "L2_norm")
          norm = VectorTools::L2_norm;
        else if (str == "Lp_norm")
          norm = VectorTools::Lp_norm;
        else if (str == "Linfty_norm")
          norm = VectorTools::Linfty_norm;
        else if (str == "H1_seminorm")
          norm = VectorTools::H1_seminorm;
        else if (str == "Hdiv_seminorm")
          norm = VectorTools::Hdiv_seminorm;
        else if (str == "H1_norm")
          norm = VectorTools::H1_norm;
        else if (str == "W1p_seminorm")
          norm = VectorTools::W1p_seminorm;
        else if (str == "W1infty_seminorm")
          norm = VectorTools::W1infty_seminorm;
        else if (str == "W1infty_norm")
          norm = VectorTools::W1infty_norm;
        else if (str == "W1p_norm")
          norm = VectorTools::W1p_norm;
        else
          {
            AssertThrow(false, ExcMessage("Didn't recognize a norm type."));
          }
        return norm;
      }
    };
  } // namespace Tools
} // namespace Patterns

DEAL_II_NAMESPACE_CLOSE

#endif // dealii_vector_tools_common_h
