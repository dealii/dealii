// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2018 by the deal.II authors
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

#ifndef dealii_smoothness_estimator_h
#define dealii_smoothness_estimator_h


#include <deal.II/base/config.h>

#include <deal.II/lac/vector.h>

#include <deal.II/numerics/vector_tools.h>

#include <vector>


DEAL_II_NAMESPACE_OPEN


/**
 * Estimate the smoothness of a solution based on the decay of coefficients from
 * a series expansion.
 *
 * From the definition, we can write our series expansion $\hat U_{\bf k}$ as a
 * matrix product
 * @f[
 *    \hat U_{\bf k}
 *    = {\cal F}_{{\bf k},j} u_j,
 * @f]
 * with $u_j$ the coefficients and ${\cal F}_{{\bf k},j}$ the transformation
 * matrix from the expansion. We use the classes FESeries::Fourier and
 * FESeries::Legendre to determine all coefficients $u_j$.
 *
 * The next step is that we have to estimate how fast these coefficients
 * decay with $|{\bf k}|$. Thus, we perform a least-squares fit
 * @f[
 *    \min_{\alpha,\mu}
 *    \frac 12 \sum_{{\bf k}, |{\bf k}|\le N}
 *    \left( |\hat U_{\bf k}| - \alpha |{\bf k}|^{-\mu}\right)^2
 * @f]
 * with linear regressions coefficients $\alpha$ and $\mu$. For simplification,
 * we apply a logarithm on our minimization problem
 * @f[
 *    \min_{\beta,\mu}
 *    Q(\beta,\mu) =
 *    \frac 12 \sum_{{\bf k}, |{\bf k}|\le N}
 *    \left( \ln |\hat U_{\bf k}| - \beta + \mu \ln |{\bf k}|\right)^2,
 * @f]
 * where $\beta=\ln \alpha$. This is now a problem for which the
 * optimality conditions $\frac{\partial Q}{\partial\beta}=0,
 * \frac{\partial Q}{\partial\mu}=0$, are linear in $\beta,\mu$. We can
 * write these conditions as follows:
 * @f[
 *    \left(\begin{array}{cc}
 *    \sum_{{\bf k}, |{\bf k}|\le N} 1 &
 *    \sum_{{\bf k}, |{\bf k}|\le N} \ln |{\bf k}|
 *    \\
 *    \sum_{{\bf k}, |{\bf k}|\le N} \ln |{\bf k}| &
 *    \sum_{{\bf k}, |{\bf k}|\le N} (\ln |{\bf k}|)^2
 *    \end{array}\right)
 *    \left(\begin{array}{c}
 *    \beta \\ -\mu
 *    \end{array}\right)
 *    =
 *    \left(\begin{array}{c}
 *    \sum_{{\bf k}, |{\bf k}|\le N} \ln |\hat U_{{\bf k}}|
 *    \\
 *    \sum_{{\bf k}, |{\bf k}|\le N} \ln |\hat U_{{\bf k}}| \ln |{\bf k}|
 *    \end{array}\right)
 * @f]
 * Solving for $\beta$ and $\mu$ is nothing else but a linear regression fit and
 * to do that we will use FESeries::linear_regression().
 *
 * While we are not particularly interested in the actual value of
 * $\beta$, the formula above gives us a mean to calculate the value of
 * the exponent $\mu$ that we can then use to determine that
 * $\hat u(\hat{\bf x})$ is in $H^s(\hat K)$ with $s=\mu-\frac d2$. These
 * Sobolev indices $s$ will suffice as our smoothness estimators and will be
 * calculated on each cell for any provided solution.
 *
 * @note An extensive demonstration of the use of these functions is provided in step-27.
 *
 * @ingroup numerics
 * @author Denis Davydov, 2016, Marc Fehling, 2018
 */
namespace SmoothnessEstimator
{
  /**
   * Estimates the smoothness of the provided solutions using the method of
   * decaying coefficents as outlined above.
   *
   * The @p regression_strategy parameter determines which norm will be used on the subset of
   * coeffiecients $\mathbf{k}$ with the same absolute value $|\mathbf{k}|$.
   * Default is VectorTools::Linfty_norm for a maximum approximation.
   *
   * Smoothness indicators will be calculated for each solution in @p all_solutions
   * and stored in @p all_smoothness_indicators in the same order.
   *
   * An individual @p fe_series object can be supplied, which has to be constructed with the
   * same FECollection object as the @p dof_handler.
   */
  template <typename FESeriesType, typename DoFHandlerType, typename VectorType>
  void
  estimate_by_coeff_decay(
    FESeriesType &                         fe_series,
    const DoFHandlerType &                 dof_handler,
    const std::vector<const VectorType *> &all_solutions,
    const std::vector<Vector<float> *> &   all_smoothness_indicators,
    const VectorTools::NormType regression_strategy = VectorTools::Linfty_norm);

  /**
   * Same as the function above, only for one @p solution vector.
   */
  template <typename FESeriesType, typename DoFHandlerType, typename VectorType>
  void
  estimate_by_coeff_decay(
    FESeriesType &              fe_series,
    const DoFHandlerType &      dof_handler,
    const VectorType &          solution,
    Vector<float> &             smoothness_indicators,
    const VectorTools::NormType regression_strategy = VectorTools::Linfty_norm);

  /**
   * Same as the function above, but with a default configuration for the chosen
   * series expansion, using 2-point Gaussian quadrature for each finite
   * element.
   *
   * Provide the desired series expansion as a template argument, i.e.
   * @code
   *    SmoothnessEstimator::estimate_by_coeff_decay<FESeries::Fourier<dim>>(
   *      dof_handler, all_solutions, all_smoothness_indicators);
   * @endcode
   */
  template <typename FESeriesType, typename DoFHandlerType, typename VectorType>
  void
  estimate_by_coeff_decay(
    const DoFHandlerType &                 dof_handler,
    const std::vector<const VectorType *> &all_solutions,
    const std::vector<Vector<float> *> &   all_smoothness_indicators,
    const VectorTools::NormType regression_strategy = VectorTools::Linfty_norm);

  /**
   * Same as the function above, only for one @p solution vector.
   */
  template <typename FESeriesType, typename DoFHandlerType, typename VectorType>
  void
  estimate_by_coeff_decay(
    const DoFHandlerType &      dof_handler,
    const VectorType &          solution,
    Vector<float> &             smoothness_indicators,
    const VectorTools::NormType regression_strategy = VectorTools::Linfty_norm);
} // namespace SmoothnessEstimator


DEAL_II_NAMESPACE_CLOSE

#endif
