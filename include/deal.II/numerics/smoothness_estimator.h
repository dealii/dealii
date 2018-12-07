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

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_series.h>

#include <deal.II/hp/dof_handler.h>

#include <deal.II/lac/vector.h>

#include <deal.II/numerics/vector_tools.h>

#include <vector>


DEAL_II_NAMESPACE_OPEN


/**
 * A namespace for various smoothness estimation strategies for hp-adaptive FEM.
 */
namespace SmoothnessEstimator
{
  /**
   * Estimate smoothness from decay of Legendre absolute values of coefficients
   * on the reference cell.
   *
   * In one dimension, the finite element solution on the reference element with
   * polynomial degree $p$ can be written as
   * @f[
   *    u_h(\hat x) = \sum_{j=0}^{p} a_j P_j(\hat x)
   * @f]
   * where $\{P_j(x)\}$ are Legendre polynomials. The decay of the coefficients
   * is estimated by performing the linear regression fit of
   * @f[
   *   \ln |a_j| \sim C - \sigma j
   * @f]
   * for $j=0,..,p$. The rate of the decay $\sigma$ can be used to estimate the
   * smoothness. For example, one strategy to implement hp-refinement
   * criteria is to perform p-refinement if $\sigma>1$.
   *
   * Extension to higher dimension is done by performing the fit in each
   * coordinate direction separately and then taking the lowest value of
   * $\sigma$.
   *
   * For each input vector of degrees of freedom defined on a DoFHandler,
   * this function returns a vector with as many elements as there are cells
   * where each element contains $\exp(-\sigma)$, which is a so-called
   * analyticity (see references below).
   *
   * @param [in] fe_series FESeries::Legendre object to calculate coefficients.
   * This object needs to be initialized to have at least $p+1$ coefficients in
   * each direction, where $p$ is the maximum polynomial degree to be used.
   * @param [in] dof_hander An hp::DoFHandler
   * @param [in] all_solutions A vector of pointers to the solution vectors
   * @param [out] all_smoothness_indicators A vector of pointers to the smoothness indicators for each @p all_solutions.
   * @param [in] coefficients_predicate A predicate to select Legendre
   * coefficients $a_j \;\; j=0\dots p$ for linear regression in each coordinate
   * direction. The user is responsible for updating the vector of `flags`
   * provided to this function. Note that its size is $p+1$, where $p$ is the
   * polynomial degree of the FE basis on a given element. Default
   * implementation will use all Legendre coefficients in each coordinate
   * direction, i.e. set all elements of the vector to `true`.
   * @param [in] smallest_abs_coefficient The smallest absolute value of the
   * coefficient to be used in linear regression in each coordinate direction.
   * Note that Legendre coefficients of some functions may have a repeating
   * pattern of zero coefficients (i.e. for functions that are locally symmetric
   * or antisymmetric about the midpoint of the element in any coordinate
   * direction). Thus this parameters allows to ingore small (in absolute value)
   * coefficients within the linear regression fit. In case there are less than
   * two non-zero coefficients for a coordinate direction, this direction will
   * be skipped. If all coefficients are zero, the returned value for this cell
   * will be zero (i.e. corresponding to the $\sigma=\infty$).
   *
   * For more theoretical details see
   * @code{.bib}
   * @Article{Mavriplis1994,
   *  author    = {Mavriplis, Catherine},
   *  title     = {Adaptive mesh strategies for the spectral element method},
   *  journal   = {{Computer Methods in Applied Mechanics and Engineering}},
   *  year      = {1994},
   *  volume    = {116},
   *  number    = {1},
   *  pages     = {77--86},
   *  publisher = {Elsevier},
   * }
   * @article{Houston2005,
   *  author    = {Houston, Paul and S{\"u}li, Endre},
   *  title     = {A note on the design of hp-adaptive finite element
   *               methods for elliptic partial differential equations},
   *  journal   = {{Computer Methods in Applied Mechanics and Engineering}},
   *  number    = {2},
   *  pages     = {229--243},
   *  publisher = {Elsevier},
   *  volume    = {194},
   *  year      = {2005}
   * }
   * @article{Eibner2007,
   *  author    = {Eibner, Tino and Melenk, Jens Markus},
   *  title     = {An adaptive strategy for hp-FEM based on testing for
   *               analyticity},
   *  journal   = {{Computational Mechanics}},
   *  year      = {2007},
   *  volume    = {39},
   *  number    = {5},
   *  pages     = {575--595},
   *  publisher = {Springer},
   * }
   * @endcode
   * and for the application within the deal.II:
   * @code{.bib}
   * @article{Davydov2017,
   *   author  = {Denis Davydov and Tymofiy Gerasimov and Jean-Paul Pelteret and
   *              Paul Steinmann},
   *   title   = {Convergence study of the h-adaptive PUM and the hp-adaptive
   *              FEM applied to eigenvalue problems in quantum mechanics},
   *   journal = {{Advanced Modeling and Simulation in Engineering Sciences}},
   *   year    = {2017},
   *   volume  = {4},
   *   number  = {1},
   *   pages   = {7},
   *   issn    = {2213-7467},
   *   doi     = {10.1186/s40323-017-0093-0},
   * }
   * @endcode
   *
   * @ingroup numerics
   * @author Denis Davydov, 2018
   */
  template <int dim, int spacedim, typename VectorType>
  void
  legendre_coefficient_decay(
    FESeries::Legendre<dim, spacedim> &    fe_series,
    const hp::DoFHandler<dim, spacedim> &  dof_handler,
    const std::vector<const VectorType *> &all_solutions,
    const std::vector<Vector<float> *> &   all_smoothness_indicators,
    const std::function<void(std::vector<bool> &flags)> coefficients_predicate =
      [](std::vector<bool> &flags) -> void {
      std::fill(flags.begin(), flags.end(), true);
    },
    const double smallest_abs_coefficient = 1e-10);

  /**
   * Same as above, but for a single solution vector.
   */
  template <int dim, int spacedim, typename VectorType>
  void
  legendre_coefficient_decay(
    FESeries::Legendre<dim, spacedim> &                 fe_series,
    const hp::DoFHandler<dim, spacedim> &               dof_handler,
    const VectorType &                                  solution,
    Vector<float> &                                     smoothness_indicators,
    const std::function<void(std::vector<bool> &flags)> coefficients_predicate =
      [](std::vector<bool> &flags) -> void {
      std::fill(flags.begin(), flags.end(), true);
    },
    const double smallest_abs_coefficient = 1e-10);

  /**
   * Same as above, but for a single solution vector and with the default
   * FESeries::Legendre.
   */
  template <int dim, int spacedim, typename VectorType>
  void
  legendre_coefficient_decay(
    const hp::DoFHandler<dim, spacedim> &               dof_handler,
    const VectorType &                                  solution,
    Vector<float> &                                     smoothness_indicators,
    const std::function<void(std::vector<bool> &flags)> coefficients_predicate =
      [](std::vector<bool> &flags) -> void {
      std::fill(flags.begin(), flags.end(), true);
    },
    const double smallest_abs_coefficient = 1e-10);

  /**
   * Estimate the smoothness of a solution based on the decay of coefficients
   * from a series expansion.
   *
   * From the definition, we can write our series expansion $\hat U_{\bf k}$ as
   * a matrix product
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
   * with linear regressions coefficients $\alpha$ and $\mu$. For
   * simplification, we apply a logarithm on our minimization problem
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
   * Solving for $\beta$ and $\mu$ is nothing else but a linear regression fit
   * and to do that we will use FESeries::linear_regression().
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
   * The @p regression_strategy parameter determines which norm will be used on the subset of
   * coeffiecients $\mathbf{k}$ with the same absolute value $|\mathbf{k}|$.
   * Default is VectorTools::Linfty_norm for a maximum approximation.
   *
   * Smoothness indicators will be calculated for each solution in @p all_solutions
   * and stored in @p all_smoothness_indicators in the same order.
   *
   * An individual @p fe_series object can be supplied, which has to be constructed with the
   * same FECollection object as the @p dof_handler.
   *
   * @ingroup numerics
   * @author Denis Davydov, 2016, Marc Fehling, 2018
   *
   */
  template <int dim, int spacedim, typename VectorType>
  void
  fourier_coefficient_decay(
    FESeries::Fourier<dim, spacedim> &     fe_series,
    const hp::DoFHandler<dim, spacedim> &  dof_handler,
    const std::vector<const VectorType *> &all_solutions,
    const std::vector<Vector<float> *> &   all_smoothness_indicators,
    const VectorTools::NormType regression_strategy = VectorTools::Linfty_norm);

  /**
   * Same as the function above, only for one @p solution vector.
   */
  template <int dim, int spacedim, typename VectorType>
  void
  fourier_coefficient_decay(
    FESeries::Fourier<dim, spacedim> &   fe_series,
    const hp::DoFHandler<dim, spacedim> &dof_handler,
    const VectorType &                   solution,
    Vector<float> &                      smoothness_indicators,
    const VectorTools::NormType regression_strategy = VectorTools::Linfty_norm);

  /**
   * Same as the function above, but with a default configuration for the chosen
   * series expansion, using 2-point Gaussian quadrature for each finite
   * element.
   *
   * Provide the desired series expansion as a template argument, i.e.
   * @code
   *    SmoothnessEstimator::estimate_by_coefficient_decay<FESeries::Fourier<dim>>(
   *      dof_handler, all_solutions, all_smoothness_indicators);
   * @endcode
   */
  template <int dim, int spacedim, typename VectorType>
  void
  fourier_coefficient_decay(
    const hp::DoFHandler<dim, spacedim> &  dof_handler,
    const std::vector<const VectorType *> &all_solutions,
    const std::vector<Vector<float> *> &   all_smoothness_indicators,
    const VectorTools::NormType regression_strategy = VectorTools::Linfty_norm);

  /**
   * Same as the function above, only for one @p solution vector.
   */
  template <int dim, int spacedim, typename VectorType>
  void
  fourier_coefficient_decay(
    const hp::DoFHandler<dim, spacedim> &dof_handler,
    const VectorType &                   solution,
    Vector<float> &                      smoothness_indicators,
    const VectorTools::NormType regression_strategy = VectorTools::Linfty_norm);
} // namespace SmoothnessEstimator


DEAL_II_NAMESPACE_CLOSE

#endif
