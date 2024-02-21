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

#ifndef dealii_smoothness_estimator_h
#define dealii_smoothness_estimator_h


#include <deal.II/base/config.h>

#include <deal.II/fe/component_mask.h>

#include <deal.II/numerics/vector_tools_common.h>

#include <functional>
#include <vector>


DEAL_II_NAMESPACE_OPEN


// forward declarations
#ifndef DOXYGEN
template <typename Number>
class Vector;

template <int dim, int spacedim>
DEAL_II_CXX20_REQUIRES((concepts::is_valid_dim_spacedim<dim, spacedim>))
class DoFHandler;

namespace FESeries
{
  template <int dim, int spacedim>
  class Fourier;
  template <int dim, int spacedim>
  class Legendre;
} // namespace FESeries

namespace hp
{
  template <int dim, int spacedim>
  class FECollection;
} // namespace hp
#endif


/**
 * A namespace for various smoothness estimation strategies for hp-adaptive FEM.
 *
 * Smoothness estimation is one strategy to decide whether a cell with a large
 * error estimate should undergo h- or p-refinement. Typical strategies decide
 * to increase the polynomial degree on a cell if the solution is particularly
 * smooth, whereas one would refine the mesh if the solution on the cell is
 * singular, has kinks in some derivative, or is otherwise not particularly
 * smooth. All of these strategies rely on a way to identify how "smooth" a
 * function is on a given cell.
 */
namespace SmoothnessEstimator
{
  /**
   * Smoothness estimation strategy based on the decay of Legendre expansion
   * coefficients.
   *
   * In one dimension, the finite element solution on cell $K$ with polynomial
   * degree $p$ can be written as
   * @f{eqnarray*}{
   *    u_h(x) &=& \sum_j u_j \varphi_j (x) \\
   *    u_{h, k}(x) &=& \sum_{k=0}^{p} a_k \widetilde P_k (x),
   *    \quad a_k = \sum_j {\cal L}_{k,j} u_j
   * @f}
   * where $u_j$ are degrees of freedom and $\varphi_j$ are the corresponding
   * shape functions. $\{\widetilde P_k(x)\}$ are Legendre polynomials on cell
   * $K$. $a_k$ and ${\cal L}_{k,j}$ are coefficients and transformation
   * matrices from the Legendre expansion of each shape function. For practical
   * reasons, we will perform the calculation of these matrices and coefficients
   * only on the reference cell $\hat K$. We only have to calculate the
   * transformation matrices once this way. However, results are only applicable
   * if the mapping from the reference cell to the actual cell is affine. We use
   * the class FESeries::Legendre to determine all coefficients $a_k$.
   *
   * A function is analytic, i.e., representable by a power series, if and only
   * if their Legendre expansion coefficients decay as (see @cite eibner2007hp)
   * @f[
   *   |a_k| \sim c \, \exp(-\sigma k)
   * @f]
   * We determine their decay rate $\sigma$ by performing the linear regression
   * fit of
   * @f[
   *   \ln |a_k| \sim C - \sigma k
   * @f]
   * for $k=0,\ldots,p$, with $p$ the polynomial degree of the finite element.
   * The rate of the decay $\sigma$ can be used to estimate the smoothness. For
   * example, one strategy to implement hp-refinement criteria is to perform
   * p-refinement if $\sigma>1$ (see @cite mavriplis1994hp).
   */
  namespace Legendre
  {
    /**
     * In this variant of the estimation strategy for higher dimensions, we will
     * consider all mode vectors $\bf k$ describing Legendre polynomials
     * $\widetilde P_{\bf k}$ and perform one least-squares fit over all
     * coefficients at once. If there are multiple coefficients corresponding to
     * the same absolute value of modes $\|{\bf k}\|_1$, we take the maximum
     * among those. Thus, the least-squares fit is performed on
     * @f{eqnarray*}{
     *   \widetilde P_{\bf k}({\bf x}) &=&
     *     \widetilde P_{k_1} (x_1) \ldots \widetilde P_{k_d} (x_d) \\
     *   \ln \left( \max\limits_{\|{\bf k}\|_1} |a_{\bf k}| \right) &\sim&
     *     C - \sigma \|{\bf k}\|_1
     * @f}
     * for ${\bf k}=(k_1,\ldots,k_d)$ and $k_i=0,\ldots,p$, with $p$ the
     * polynomial degree of the finite element.
     *
     * For a finite element approximation @p solution, this function writes the
     * decay rate for every cell into the output vector @p smoothness_indicators.
     *
     * @param [in] fe_legendre FESeries::Legendre object to calculate coefficients.
     * This object needs to be initialized to have at least $p+1$ coefficients
     * in each direction for every finite element in the collection, where $p$
     * is its polynomial degree.
     * @param [in] dof_handler A DoFHandler.
     * @param [in] solution A solution vector.
     * @param [out] smoothness_indicators A vector for smoothness indicators.
     * @param [in] regression_strategy Determines which norm will be used on the
     * subset of coefficients $\mathbf{k}$ with the same absolute value
     * $\|{\bf k}\|_1$. Default is VectorTools::Linfty_norm for a maximum
     * approximation.
     * @param [in] smallest_abs_coefficient The smallest absolute value of the
     * coefficient to be used in linear regression. Note that Legendre
     * coefficients of some functions may have a repeating pattern of zero
     * coefficients (i.e. for functions that are locally symmetric or
     * antisymmetric about the midpoint of the element in any coordinate
     * direction). Thus this parameters allows to ignore small (in absolute
     * value) coefficients within the linear regression fit. In case there are
     * less than two nonzero coefficients, the returned value for this cell will
     * be $\sigma=\infty$.
     * @param [in] only_flagged_cells Smoothness indicators are usually used to
     * decide whether to perform h- or p-adaptation. So in most cases, we only
     * need to calculate those indicators on cells flagged for refinement or
     * coarsening. This parameter controls whether this particular subset or all
     * cells will be considered. By default, all cells will be considered. When
     * only flagged cells are supposed to be considered, smoothness indicators
     * will only be set on those vector entries of flagged cells; the others
     * will be set to a signaling NaN.
     *
     * For more theoretical details see @cite mavriplis1994hp
     * @cite houston2005hp @cite eibner2007hp.
     */
    template <int dim, int spacedim, typename VectorType>
    void
    coefficient_decay(FESeries::Legendre<dim, spacedim> &fe_legendre,
                      const DoFHandler<dim, spacedim>   &dof_handler,
                      const VectorType                  &solution,
                      Vector<float>                     &smoothness_indicators,
                      const VectorTools::NormType        regression_strategy =
                        VectorTools::Linfty_norm,
                      const double smallest_abs_coefficient = 1e-10,
                      const bool   only_flagged_cells       = false);

    /**
     * In this variant of the estimation strategy for higher dimensions, we only
     * consider modes in each coordinate direction, i.e., only mode vectors
     * $\bf k$ with one nonzero entry. We perform the least-squares fit in
     * each coordinate direction separately and take the lowest decay rate
     * $\sigma$ among them.
     *
     * For a finite element approximation @p solution, this function writes the
     * decay rate for every cell into the output vector @p smoothness_indicators.
     *
     * @param [in] fe_legendre FESeries::Legendre object to calculate coefficients.
     * This object needs to be initialized to have at least $p+1$ coefficients
     * in each direction, where $p$ is the maximum polynomial degree to be used.
     * @param [in] dof_handler A DoFHandler
     * @param [in] solution A solution vector
     * @param [out] smoothness_indicators A vector for smoothness indicators
     * @param [in] coefficients_predicate A predicate to select Legendre
     * coefficients $a_j$, $j=0,\ldots,p$ for linear regression in each
     * coordinate direction. The user is responsible for updating the vector of
     * `flags` provided to this function. Note that its size is $p+1$, where $p$
     * is the polynomial degree of the FE basis on a given element. The default
     * implementation will use all Legendre coefficients in each coordinate
     * direction, i.e. set all elements of the vector to `true`.
     * @param [in] smallest_abs_coefficient The smallest absolute value of the
     * coefficient to be used in linear regression in each coordinate direction.
     * Note that Legendre coefficients of some functions may have a repeating
     * pattern of zero coefficients (i.e. for functions that are locally
     * symmetric or antisymmetric about the midpoint of the element in any
     * coordinate direction). Thus this parameters allows to ignore small (in
     * absolute value) coefficients within the linear regression fit. In case
     * there are less than two nonzero coefficients for a coordinate direction,
     * this direction will be skipped. If all coefficients are zero, the
     * returned value for this cell will be $\sigma=\infty$.
     * @param [in] only_flagged_cells Smoothness indicators are usually used to
     * decide whether to perform h- or p-adaptation. So in most cases, we only
     * need to calculate those indicators on cells flagged for refinement or
     * coarsening. This parameter controls whether this particular subset or all
     * cells will be considered. By default, all cells will be considered. When
     * only flagged cells are supposed to be considered, smoothness indicators
     * will only be set on those vector entries of flagged cells; the others
     * will be set to NaN.
     *
     * For more theoretical details and the application within the deal.II
     * library see @cite davydov2017hp.
     */
    template <int dim, int spacedim, typename VectorType>
    void
    coefficient_decay_per_direction(
      FESeries::Legendre<dim, spacedim> &fe_legendre,
      const DoFHandler<dim, spacedim>   &dof_handler,
      const VectorType                  &solution,
      Vector<float>                     &smoothness_indicators,
      const ComponentMask               &coefficients_predicate   = {},
      const double                       smallest_abs_coefficient = 1e-10,
      const bool                         only_flagged_cells       = false);

    /**
     * Returns a FESeries::Legendre object for Legendre series expansions with
     * the default configuration for smoothness estimation purposes.
     *
     * For each finite element of the provided @p fe_collection, we use as many
     * modes as its polynomial degree plus two. This includes the first Legendre
     * polynomial which is just a constant. Further for each element, we use a
     * Gaussian quadrature designed to yield exact results for the highest order
     * Legendre polynomial used.
     *
     * As the Legendre expansion can only be performed on scalar fields, this
     * class does not operate on vector-valued finite elements and will
     * therefore throw an assertion. However, each component of a finite element
     * field can be treated as a scalar field, respectively, on which Legendre
     * expansions are again possible. For this purpose, the optional parameter
     * @p component defines which component of each FiniteElement will be used.
     * The default value of @p component only applies to scalar FEs, in which
     * case it indicates that the sole component is to be decomposed. For
     * vector-valued FEs, a non-default value must be explicitly provided.
     */
    template <int dim, int spacedim>
    FESeries::Legendre<dim, spacedim>
    default_fe_series(
      const hp::FECollection<dim, spacedim> &fe_collection,
      const unsigned int component = numbers::invalid_unsigned_int);
  } // namespace Legendre



  /**
   * Smoothness estimation strategy based on the decay of Fourier expansion
   * coefficients.
   *
   * From the definition, we can write our Fourier series expansion
   * $a_{\bf k}$ of the finite element solution on cell $K$ with polynomial
   * degree $p$ as a matrix product
   * @f{eqnarray*}{
   *    u_h({\bf x}) &=& \sum_j u_j \varphi_j ({\bf x}) \\
   *    u_{h, {\bf k}}({\bf x}) &=&
   *      \sum_{{\bf k}, \|{\bf k}\|\le p} a_{\bf k} \phi_{\bf k}({\bf x}),
   *      \quad a_{\bf k} = \sum_j {\cal F}_{{\bf k},j} u_j
   * @f}
   * with $u_j$ the degrees of freedom and $\varphi_j$ the corresponding shape
   * functions. $\{\phi_{\bf k}({\bf x}) = \exp(i \, 2 \pi \, {\bf k} \cdot
   * {\bf x}) \}$ are exponential functions on cell $K$. $a_{\bf k}$ and ${\cal
   * F}_{{\bf k},j}$ are coefficients and transformation matrices from the
   * Fourier expansion of each shape function. For practical reasons, we will
   * perform the calculation of these matrices and coefficients only on the
   * reference cell $\hat K$. We only have to calculate the transformation
   * matrices once this way. However, results are only applicable if mapping
   * from the reference cell to the actual cell is linear. We use the class
   * FESeries::Fourier to determine all coefficients $a_{\bf k}$.
   *
   * If the finite element approximation on cell $K$ is part of the Hilbert
   * space $H^s(K)$, then the following integral must exist for both the finite
   * element and spectral representation of our solution
   * @f{eqnarray*}{
   *   \| \nabla^s u_h({\bf x}) \|_{L^2(K)}^2 &=&
   *     \int\limits_K \left| \nabla^s u_h({\bf x}) \right|^2 d{\bf x} <
   *     \infty \\
   *   \| \nabla^s u_{h, {\bf k}}({\bf x}) \|_{L^2(K)}^2 &=&
   *     \int\limits_K \left| \sum\limits_{\bf k} (-i \, 2 \pi {\bf k})^s \,
   *     a_{\bf k} \, \phi_{\bf k}({\bf x}) \right|^2 d{\bf x} =
   *     (2 \pi)^{2s} \sum\limits_{\bf k} \left| a_{\bf k} \right|^2
   *     \|{\bf k}\|_2^{2s} < \infty
   * @f}
   * The sum is finite only if the summands decay at least with order
   * @f[
   *   |a_{\bf k}|^2 \|{\bf k}\|_2^{2s} \|{\bf k}\|_2^{d - 1} =
   *     {\cal O}\left( \|{\bf k}\|_2^{-1-\epsilon} \right)
   * @f]
   * for all $\epsilon > 0$. The additional factor stems from the fact that,
   * since we sum over all multi-indices ${\bf k}$ that are located on a
   * dim-dimensional sphere, we actually have, up to a constant,
   * $\|{\bf k}\|_2^{d-1}$ modes located in each increment $\|{\bf k}\|_2 +
   * d\|{\bf k}\|_2$ that need to be taken into account. By a comparison of
   * exponents, we can rewrite this condition as
   * @f[
   *   |a_{\bf k}| = {\cal O}\left(\|{\bf k}\|_2^
   *     {-\left(s + \frac d2 + \epsilon \right)} \right)
   * @f]
   *
   * The next step is to estimate how fast these coefficients
   * decay with $\|{\bf k}\|_2$. Thus, we perform a least-squares fit
   * @f[
   *    \min_{\alpha,\sigma}
   *    \frac 12 \sum_{{\bf k}, \|{\bf k}\|_2 \le p}
   *    \left( |a_{\bf k}| - \alpha \|{\bf k}\|_2^{-\sigma}\right)^2
   * @f]
   * with regression coefficients $\alpha$ and $\sigma$. For simplification, we
   * apply a logarithm on our minimization problem
   * @f[
   *    \min_{\beta,\sigma}
   *    Q(\beta,\sigma) =
   *    \frac 12 \sum_{{\bf k}, \|{\bf k}\|_2 \le p}
   *    \left( \ln |a_{\bf k}| - \beta + \sigma \ln \|{\bf k}\|_2
   * \right)^2,
   * @f]
   * where $\beta=\ln \alpha$. This is now a problem for which the optimality
   * conditions $\frac{\partial Q}{\partial\beta}=0,
   * \frac{\partial Q}{\partial\sigma}=0$, are linear in $\beta,\sigma$. We can
   * write these conditions as follows:
   * @f[
   *    \left(\begin{array}{cc}
   *    \sum_{{\bf k}, \|{\bf k}\|_2 \le p} 1 &
   *    \sum_{{\bf k}, \|{\bf k}\|_2 \le p} \ln \|{\bf k}\|_2
   *    \\
   *    \sum_{{\bf k}, \|{\bf k}\|_2 \le p} \ln \|{\bf k}\|_2 &
   *    \sum_{{\bf k}, \|{\bf k}\|_2 \le p} (\ln \|{\bf k}\|_2)^2
   *    \end{array}\right)
   *    \left(\begin{array}{c}
   *    \beta \\ -\sigma
   *    \end{array}\right)
   *    =
   *    \left(\begin{array}{c}
   *    \sum_{{\bf k}, \|{\bf k}\|_2\le p} \ln |a_{{\bf k}}|
   *    \\
   *    \sum_{{\bf k}, \|{\bf k}\|_2\le p} \ln |a_{{\bf k}}| \ln \|{\bf
   * k}\|_2 \end{array}\right)
   * @f]
   * Solving for $\beta$ and $\sigma$ is just a linear regression fit and to do
   * that we will use FESeries::linear_regression().
   *
   * While we are not particularly interested in the actual value of
   * $\beta$, the formula above gives us a means to calculate the value of
   * the exponent $\sigma$ that we can then use to determine that
   * $u(\hat{\bf x})$ is in $H^s(K)$ with $s=\sigma-\frac d2$. The
   * decay rates $\sigma$ will suffice as our smoothness indicators and
   * will be calculated on each cell for any provided solution.
   *
   * @note An extensive demonstration of the use of these functions is
   * provided in step-27.
   */
  namespace Fourier
  {
    /**
     * In this variant of the estimation strategy for higher dimensions, we will
     * consider all mode vectors $\bf k$ describing Fourier polynomials
     * $P_{\bf k}$ and perform one least-squares fit over all coefficients
     * at once. If there are multiple coefficients corresponding to the same
     * absolute value of modes $\|\bf k\|_2$, we take the maximum among those.
     * Thus, the least-squares fit is performed on
     * @f[
     *   \ln \left( \max\limits_{\|{\bf k}\|_2} |a_{\bf k}| \right) \sim
     *     C - \sigma \ln \|{\bf k}\|_2
     * @f]
     * for ${\bf k}=(k_1,\ldots,k_d)$ and $k_i=0,\ldots,p$, with $p$ the
     * polynomial degree of the finite element. We exclude the $\|{\bf k}\|_2=0$
     * modes to avoid the singularity of the logarithm.
     *
     * The @p regression_strategy parameter determines which norm will be used
     * on the subset of coefficients $\bf k$ with the same absolute value
     * $\|{\bf k}\|_2$. Default is VectorTools::Linfty_norm for a maximum
     * approximation.
     *
     * For a provided solution vector @p solution defined on a DoFHandler
     * @p dof_handler, this function returns a vector @p smoothness_indicators
     * with as many elements as there are cells where each element contains the
     * estimated regularity $\sigma$.
     *
     * A series expansion object @p fe_fourier has to be supplied, which needs
     * to be constructed with the same FECollection object as the @p dof_handler.
     *
     * The parameter @p smallest_abs_coefficient allows to ignore small (in
     * absolute value) coefficients within the linear regression fit. In case
     * there are less than two nonzero coefficients for a coordinate direction,
     * this direction will be skipped. If all coefficients are zero, the
     * returned value for this cell will be $\sigma=\infty$.
     *
     * Smoothness indicators are usually used to decide whether to perform h- or
     * p-adaptation. So in most cases, we only need to calculate those
     * indicators on cells flagged for refinement or coarsening. The parameter
     * @p only_flagged_cells controls whether this particular subset or all
     * cells will be considered. By default, all cells will be considered.
     * When only flagged cells are supposed to be considered, smoothness
     * indicators will only be set on those vector entries of flagged cells;
     * the others will be set to a signaling NaN.
     */
    template <int dim, int spacedim, typename VectorType>
    void
    coefficient_decay(FESeries::Fourier<dim, spacedim> &fe_fourier,
                      const DoFHandler<dim, spacedim>  &dof_handler,
                      const VectorType                 &solution,
                      Vector<float>                    &smoothness_indicators,
                      const VectorTools::NormType       regression_strategy =
                        VectorTools::Linfty_norm,
                      const double smallest_abs_coefficient = 1e-10,
                      const bool   only_flagged_cells       = false);

    /**
     * In this variant of the estimation strategy for higher dimensions, we only
     * consider modes in each coordinate direction, i.e., only mode vectors
     * $\bf k$ with one nonzero entry. We perform the least-squares fit in
     * each coordinate direction separately and take the lowest decay rate
     * $\sigma$ among them.
     *
     * The @p coefficients_predicate parameter selects Fourier coefficients
     * $a_j$, $j=0,\ldots,p$ for linear regression in each coordinate
     * direction. The user is responsible for updating the vector of `flags`
     * provided to this function. Note that its size is $p+1$, where $p$ is the
     * polynomial degree of the FE basis on a given element. The default
     * implementation will use all Fourier coefficients in each coordinate
     * direction, i.e., set all the elements of the vector to `true`.
     *
     * For a provided solution vector @p solution defined on a DoFHandler
     * @p dof_handler, this function returns a vector @p smoothness_indicators
     * with as many elements as there are cells where each element contains the
     * estimated regularity $\sigma$.
     *
     * A series expansion object @p fe_fourier has to be supplied, which needs
     * to be constructed with the same FECollection object as the @p dof_handler.
     *
     * The parameter @p smallest_abs_coefficient allows to ignore small (in
     * absolute value) coefficients within the linear regression fit. In case
     * there are fewer than two nonzero coefficients for a coordinate direction,
     * this direction will be skipped. If all coefficients are zero, the
     * returned value for this cell will be $\sigma=\infty$.
     *
     * Smoothness indicators are usually used to decide whether to perform h- or
     * p-adaptation. So in most cases, we only need to calculate those
     * indicators on cells flagged for refinement or coarsening. The parameter
     * @p only_flagged_cells controls whether this particular subset or all
     * cells will be considered. By default, all cells will be considered.
     * When only flagged cells are supposed to be considered, smoothness
     * indicators will only be set on those vector entries of flagged cells;
     * the others will be set to a signaling NaN.
     */
    template <int dim, int spacedim, typename VectorType>
    void
    coefficient_decay_per_direction(
      FESeries::Fourier<dim, spacedim> &fe_fourier,
      const DoFHandler<dim, spacedim>  &dof_handler,
      const VectorType                 &solution,
      Vector<float>                    &smoothness_indicators,
      const ComponentMask              &coefficients_predicate   = {},
      const double                      smallest_abs_coefficient = 1e-10,
      const bool                        only_flagged_cells       = false);

    /**
     * Returns a FESeries::Fourier object for Fourier series expansions with
     * the default configuration for smoothness estimation purposes.
     *
     * For each finite element of the provided @p fe_collection, we use as many
     * modes as its polynomial degree plus two. Further for each element, we use
     * a 5-point Gaussian quarature iterated in each dimension by the maximal
     * wave number, which is the number of modes decreased by one since we start
     * with $k = 0$.
     *
     * As the Fourier expansion can only be performed on scalar fields, this
     * class does not operate on vector-valued finite elements and will
     * therefore throw an assertion. However, each component of a finite element
     * field can be treated as a scalar field, respectively, on which Fourier
     * expansions are again possible. For this purpose, the optional parameter
     * @p component defines which component of each FiniteElement will be used.
     * The default value of @p component only applies to scalar FEs, in which
     * case it indicates that the sole component is to be decomposed. For
     * vector-valued FEs, a non-default value must be explicitly provided.
     */
    template <int dim, int spacedim>
    FESeries::Fourier<dim, spacedim>
    default_fe_series(
      const hp::FECollection<dim, spacedim> &fe_collection,
      const unsigned int component = numbers::invalid_unsigned_int);
  } // namespace Fourier
} // namespace SmoothnessEstimator


DEAL_II_NAMESPACE_CLOSE

#endif
