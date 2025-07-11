// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2016 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_fe_series_h
#define dealii_fe_series_h



#include <deal.II/base/config.h>

#include <deal.II/base/enable_observer_pointer.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/table.h>
#include <deal.II/base/table_indices.h>
#include <deal.II/base/tensor.h>

#include <deal.II/hp/fe_collection.h>
#include <deal.II/hp/q_collection.h>

#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/vector.h>

#include <deal.II/numerics/vector_tools_common.h>

#include <memory>
#include <string>
#include <vector>


DEAL_II_NAMESPACE_OPEN


/**
 * @addtogroup fe
 * @{
 */


/**
 * This namespace offers functions to calculate expansion series of the
 * solution on the reference element. Coefficients of expansion are often used
 * to estimate local smoothness of the underlying FiniteElement field to decide
 * on h- or p-adaptive refinement strategy.
 */
namespace FESeries
{
  /**
   * A class to calculate expansion of a scalar FE (or a single component
   * of vector-valued FE) field into Fourier series on a reference element.
   * The exponential form of the Fourier series is  based on completeness
   * and Hermitian orthogonality of the set of exponential
   * functions $ \phi_{\bf k}({\bf x}) = \exp(2 \pi i\, {\bf k} \cdot {\bf x})$.
   * For example in 1d the L2-orthogonality condition reads
   * @f[
   *   \int_0^1 \phi_k(x) \phi_l^\ast(x) dx=\delta_{kl}.
   * @f]
   * Note that $ \phi_{\bf k} = \phi_{-\bf k}^\ast $.
   *
   * The arbitrary scalar FE field on the reference element can be expanded in
   * the complete orthogonal exponential basis as
   * @f[
   *    u({\bf x})
   *    = \sum_{\bf k} c_{\bf k} \phi_{\bf k}({\bf x}).
   * @f]
   * From the orthogonality property of the basis, it follows that
   * @f[
   *    c_{\bf k} =
   *    \int_{[0,1]^d} u({\bf x}) \phi_{\bf k}^\ast ({\bf x}) d{\bf x}\,.
   * @f]
   * It is this complex-valued expansion coefficients, that are calculated by
   * this class. Note that $ u({\bf x}) = \sum_i u_i N_i({\bf x})$,
   * where $ N_i({\bf x}) $ are real-valued FiniteElement shape functions.
   * Consequently $ c_{\bf k} \equiv c_{-\bf k}^\ast $ and
   * we only need to compute $ c_{\bf k} $ for positive indices
   * $ \bf k $ .
   */
  template <int dim, int spacedim = dim>
  class Fourier : public EnableObserverPointer
  {
  public:
    using CoefficientType = typename std::complex<double>;

    /**
     * Constructor that initializes all required data structures.
     *
     * The @p n_coefficients_per_direction defines the number of coefficients in
     * each direction, @p fe_collection is the hp::FECollection for which
     * expansion will be used and @p q_collection is the hp::QCollection used to
     * integrate the expansion for each FiniteElement in @p fe_collection.
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
    Fourier(const std::vector<unsigned int>       &n_coefficients_per_direction,
            const hp::FECollection<dim, spacedim> &fe_collection,
            const hp::QCollection<dim>            &q_collection,
            const unsigned int component = numbers::invalid_unsigned_int);

    /**
     * Calculate @p fourier_coefficients of the cell vector field given by
     * @p local_dof_values corresponding to FiniteElement with
     * @p cell_active_fe_index .
     */
    template <typename Number>
    void
    calculate(const dealii::Vector<Number> &local_dof_values,
              const unsigned int            cell_active_fe_index,
              Table<dim, CoefficientType>  &fourier_coefficients);

    /**
     * Return the number of coefficients in each coordinate direction for the
     * finite element associated with @p index in the provided hp::FECollection.
     */
    unsigned int
    get_n_coefficients_per_direction(const unsigned int index) const;

    /**
     * Calculate all transformation matrices to transfer the finite element
     * solution to the series expansion representation.
     *
     * These matrices will be generated on demand by calling calculate() and
     * stored for recurring purposes. Usually, this operation consumes a lot of
     * workload. With this function, all matrices will be calculated in advance.
     * This way, we can separate their costly generation from the actual
     * application.
     */
    void
    precalculate_all_transformation_matrices();

    /**
     * Write all transformation matrices of this object to a stream for the
     * purpose of serialization.
     *
     * Since any of its transformation matrices has to be generated only once
     * for a given scenario, it is common practice to determine them in advance
     * calling precalculate_all_transformation_matrices() and keep them via
     * serialization.
     */
    template <class Archive>
    void
    save_transformation_matrices(Archive &ar, const unsigned int version);

    /**
     * Read all transformation matrices from a stream and recover them for this
     * object.
     */
    template <class Archive>
    void
    load_transformation_matrices(Archive &ar, const unsigned int version);

    /**
     * Test for equality of two series expansion objects.
     */
    bool
    operator==(const Fourier<dim, spacedim> &fourier) const;

  private:
    /**
     * Number of coefficients in each direction for each finite element in the
     * registered hp::FECollection.
     */
    const std::vector<unsigned int> n_coefficients_per_direction;

    /**
     * hp::FECollection for which transformation matrices will be calculated.
     */
    ObserverPointer<const hp::FECollection<dim, spacedim>> fe_collection;

    /**
     * hp::QCollection used in calculation of transformation matrices.
     */
    const hp::QCollection<dim> q_collection;

    /**
     * Angular frequencies $ 2 \pi {\bf k} $ .
     */
    Table<dim, Tensor<1, dim>> k_vectors;

    /**
     * Transformation matrices for each FiniteElement.
     */
    std::vector<FullMatrix<CoefficientType>> fourier_transform_matrices;

    /**
     * Auxiliary vector to store unrolled coefficients.
     */
    std::vector<CoefficientType> unrolled_coefficients;

    /**
     * Which component of FiniteElement should be used to calculate the
     * expansion.
     */
    const unsigned int component;
  };



  /**
   * A class to calculate expansion of a scalar FE (or a single component
   * of vector-valued FE) field into series of Legendre functions on a
   * reference element.
   *
   * Legendre functions are solutions to Legendre's differential equation
   * @f[
   *    \frac{d}{dx}\left([1-x^2] \frac{d}{dx} P_n(x)\right) +
   *    n[n+1] P_n(x) = 0
   * @f]
   * and can be expressed using Rodrigues' formula
   * @f[
   *    P_n(x) = \frac{1}{2^n n!} \frac{d^n}{dx^n}[x^2-1]^n.
   * @f]
   * These polynomials are orthogonal with respect to the $ L^2 $ inner
   * product on the interval $ [-1;1] $
   * @f[
   *    \int_{-1}^1 P_m(x) P_n(x) = \frac{2}{2n + 1} \delta_{mn}
   * @f]
   * and are complete.
   * A family of $ L^2 $-orthogonal polynomials on $ [0;1] $ can be
   * constructed via
   * @f[
   *    \widetilde P_m = \sqrt{2} P_m(2x-1).
   * @f]
   *
   *
   * An arbitrary scalar FE field on the reference element $ [0;1] $ can be
   * expanded in the complete orthogonal basis as
   * @f[
   *    u(x)
   *    = \sum_{m} c_m \widetilde P_{m}(x).
   * @f]
   * From the orthogonality property of the basis, it follows that
   * @f[
   *    c_m = \frac{2m+1}{2}
   *    \int_0^1 u(x) \widetilde P_m(x) dx .
   * @f]
   * This class calculates coefficients $ c_{\bf k} $ using
   * $ dim $-dimensional Legendre polynomials constructed from
   * $ \widetilde P_m(x) $ using tensor product rule.
   */
  template <int dim, int spacedim = dim>
  class Legendre : public EnableObserverPointer
  {
  public:
    using CoefficientType = double;

    /**
     * Constructor that initializes all required data structures.
     *
     * The @p n_coefficients_per_direction defines the number of coefficients in
     * each direction, @p fe_collection is the hp::FECollection for which
     * expansion will be used and @p q_collection is the hp::QCollection used to
     * integrate the expansion for each FiniteElement in @p fe_collection.
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
    Legendre(const std::vector<unsigned int> &n_coefficients_per_direction,
             const hp::FECollection<dim, spacedim> &fe_collection,
             const hp::QCollection<dim>            &q_collection,
             const unsigned int component = numbers::invalid_unsigned_int);

    /**
     * Calculate @p legendre_coefficients of the cell vector field given by
     * @p local_dof_values corresponding to FiniteElement with
     * @p cell_active_fe_index .
     */
    template <typename Number>
    void
    calculate(const dealii::Vector<Number> &local_dof_values,
              const unsigned int            cell_active_fe_index,
              Table<dim, CoefficientType>  &legendre_coefficients);

    /**
     * Return the number of coefficients in each coordinate direction for the
     * finite element associated with @p index in the provided hp::FECollection.
     */
    unsigned int
    get_n_coefficients_per_direction(const unsigned int index) const;

    /**
     * Calculate all transformation matrices to transfer the finite element
     * solution to the series expansion representation.
     *
     * These matrices will be generated on demand by calling calculate() and
     * stored for recurring purposes. Usually, this operation consumes a lot of
     * workload. With this function, all matrices will be calculated in advance.
     * This way, we can separate their costly generation from the actual
     * application.
     */
    void
    precalculate_all_transformation_matrices();

    /**
     * Write all transformation matrices of this object to a stream for the
     * purpose of serialization.
     *
     * Since any of its transformation matrices has to be generated only once
     * for a given scenario, it is common practice to determine them in advance
     * calling precalculate_all_transformation_matrices() and keep them via
     * serialization.
     */
    template <class Archive>
    void
    save_transformation_matrices(Archive &ar, const unsigned int version);

    /**
     * Read all transformation matrices from a stream and recover them for this
     * object.
     */
    template <class Archive>
    void
    load_transformation_matrices(Archive &ar, const unsigned int version);

    /**
     * Test for equality of two series expansion objects.
     */
    bool
    operator==(const Legendre<dim, spacedim> &legendre) const;

  private:
    /**
     * Number of coefficients in each direction for each finite element in the
     * registered hp::FECollection.
     */
    const std::vector<unsigned int> n_coefficients_per_direction;

    /**
     * hp::FECollection for which transformation matrices will be calculated.
     */
    ObserverPointer<const hp::FECollection<dim, spacedim>> fe_collection;

    /**
     * hp::QCollection used in calculation of transformation matrices.
     */
    const hp::QCollection<dim> q_collection;

    /**
     * Transformation matrices for each FiniteElement.
     */
    std::vector<FullMatrix<CoefficientType>> legendre_transform_matrices;

    /**
     * Auxiliary vector to store unrolled coefficients.
     */
    std::vector<CoefficientType> unrolled_coefficients;

    /**
     * Which component of FiniteElement should be used to calculate the
     * expansion.
     */
    const unsigned int component;
  };



  /**
   * Calculate the @p norm of subsets of @p coefficients defined by
   * @p predicate being constant. Return the pair of vectors of predicate values
   * and the vector of calculated subset norms.
   *
   * @p predicate should return a pair of <code>bool</code> and <code>unsigned
   * int</code>. The former is a flag whether a given TableIndices should be
   * used in calculation, whereas the latter is the unrolled value of indices
   * according to which the subsets of coefficients will be formed.
   *
   * Only those coefficients will be considered which are larger than
   * @p smallest_abs_coefficient.
   *
   * @note Only the following values of @p norm_type are implemented and make
   * sense in this case: mean, L1_norm, L2_norm, Linfty_norm. The mean norm ca
   * only be applied to real valued coefficients.
   */
  template <int dim, typename CoefficientType>
  std::pair<std::vector<unsigned int>, std::vector<double>>
  process_coefficients(const Table<dim, CoefficientType> &coefficients,
                       const std::function<std::pair<bool, unsigned int>(
                         const TableIndices<dim> &)>     &predicate,
                       const VectorTools::NormType        norm_type,
                       const double smallest_abs_coefficient = 1e-10);

  /**
   * Linear regression least-square fit of $y = k \, x + b$.
   * The size of the input vectors should be equal and more than 1.
   * The returned pair will contain $k$ (first) and $b$ (second).
   */
  std::pair<double, double>
  linear_regression(const std::vector<double> &x, const std::vector<double> &y);

} // namespace FESeries

/** @} */



#ifndef DOXYGEN

// -------------------  inline and template functions ----------------

namespace internal
{
  namespace FESeriesImplementation
  {
    template <int dim, typename CoefficientType>
    void
    fill_map_index(
      const Table<dim, CoefficientType> &coefficients,
      const TableIndices<dim>           &ind,
      const std::function<
        std::pair<bool, unsigned int>(const TableIndices<dim> &)> &predicate,
      std::map<unsigned int, std::vector<CoefficientType>> &pred_to_values)
    {
      const std::pair<bool, unsigned int> pred_pair = predicate(ind);
      // don't add a value if predicate is false
      if (pred_pair.first == false)
        return;

      const unsigned int     pred_value  = pred_pair.second;
      const CoefficientType &coeff_value = coefficients(ind);
      // If pred_value is not in the pred_to_values map, the element will be
      // created. Otherwise a reference to the existing element is returned.
      pred_to_values[pred_value].push_back(coeff_value);
    }



    template <typename CoefficientType>
    void
    fill_map(
      const Table<1, CoefficientType> &coefficients,
      const std::function<
        std::pair<bool, unsigned int>(const TableIndices<1> &)> &predicate,
      std::map<unsigned int, std::vector<CoefficientType>>      &pred_to_values)
    {
      for (unsigned int i = 0; i < coefficients.size(0); ++i)
        {
          const TableIndices<1> ind(i);
          fill_map_index(coefficients, ind, predicate, pred_to_values);
        }
    }



    template <typename CoefficientType>
    void
    fill_map(
      const Table<2, CoefficientType> &coefficients,
      const std::function<
        std::pair<bool, unsigned int>(const TableIndices<2> &)> &predicate,
      std::map<unsigned int, std::vector<CoefficientType>>      &pred_to_values)
    {
      for (unsigned int i = 0; i < coefficients.size(0); ++i)
        for (unsigned int j = 0; j < coefficients.size(1); ++j)
          {
            const TableIndices<2> ind(i, j);
            fill_map_index(coefficients, ind, predicate, pred_to_values);
          }
    }



    template <typename CoefficientType>
    void
    fill_map(
      const Table<3, CoefficientType> &coefficients,
      const std::function<
        std::pair<bool, unsigned int>(const TableIndices<3> &)> &predicate,
      std::map<unsigned int, std::vector<CoefficientType>>      &pred_to_values)
    {
      for (unsigned int i = 0; i < coefficients.size(0); ++i)
        for (unsigned int j = 0; j < coefficients.size(1); ++j)
          for (unsigned int k = 0; k < coefficients.size(2); ++k)
            {
              const TableIndices<3> ind(i, j, k);
              fill_map_index(coefficients, ind, predicate, pred_to_values);
            }
    }



    template <typename Number>
    double
    complex_mean_value(const Number &value)
    {
      return value;
    }



    template <typename Number>
    double
    complex_mean_value(const std::complex<Number> &value)
    {
      AssertThrow(false,
                  ExcMessage(
                    "FESeries::process_coefficients() can not be used with "
                    "complex-valued coefficients and VectorTools::mean norm."));
      return std::abs(value);
    }
  } // namespace FESeriesImplementation
} // namespace internal



template <int dim, typename CoefficientType>
std::pair<std::vector<unsigned int>, std::vector<double>>
FESeries::process_coefficients(
  const Table<dim, CoefficientType> &coefficients,
  const std::function<std::pair<bool, unsigned int>(const TableIndices<dim> &)>
                             &predicate,
  const VectorTools::NormType norm_type,
  const double                smallest_abs_coefficient)
{
  Assert(smallest_abs_coefficient >= 0.,
         ExcMessage("smallest_abs_coefficient should be non-negative."));

  std::vector<unsigned int> predicate_values;
  std::vector<double>       norm_values;

  // first, parse all table elements into a map of predicate values and
  // coefficients. We could have stored (predicate values ->TableIndices) map,
  // but its processing would have been much harder later on.
  std::map<unsigned int, std::vector<CoefficientType>> pred_to_values;
  internal::FESeriesImplementation::fill_map(coefficients,
                                             predicate,
                                             pred_to_values);

  // now go through the map and populate the @p norm_values based on @p norm:
  for (const auto &pred_to_value : pred_to_values)
    {
      Vector<CoefficientType> values(pred_to_value.second.cbegin(),
                                     pred_to_value.second.cend());

      double norm_value = 0;
      switch (norm_type)
        {
          case VectorTools::L2_norm:
            {
              norm_value = values.l2_norm();
              break;
            }
          case VectorTools::L1_norm:
            {
              norm_value = values.l1_norm();
              break;
            }
          case VectorTools::Linfty_norm:
            {
              norm_value = values.linfty_norm();
              break;
            }
          case VectorTools::mean:
            {
              norm_value = internal::FESeriesImplementation::complex_mean_value(
                values.mean_value());
              break;
            }
          default:
            AssertThrow(false, ExcNotImplemented());
            break;
        }

      // will use all non-zero coefficients
      if (std::abs(norm_value) > smallest_abs_coefficient)
        {
          predicate_values.push_back(pred_to_value.first);
          norm_values.push_back(norm_value);
        }
    }

  return std::make_pair(predicate_values, norm_values);
}



template <int dim, int spacedim>
template <class Archive>
inline void
FESeries::Fourier<dim, spacedim>::save_transformation_matrices(
  Archive &ar,
  const unsigned int /*version*/)
{
  // Store information about those resources which have been used to generate
  // the transformation matrices.
  // mode vector
  ar &n_coefficients_per_direction;

  // finite element collection
  unsigned int size = fe_collection->size();
  ar          &size;
  for (unsigned int i = 0; i < size; ++i)
    ar &(*fe_collection)[i].get_name();

  // quadrature collection
  size = q_collection.size();
  ar &size;
  for (unsigned int i = 0; i < size; ++i)
    ar &q_collection[i];

  // Store the actual transform matrices.
  ar &fourier_transform_matrices;
}



template <int dim, int spacedim>
template <class Archive>
inline void
FESeries::Fourier<dim, spacedim>::load_transformation_matrices(
  Archive &ar,
  const unsigned int /*version*/)
{
  // Check whether the currently registered resources are compatible with
  // the transformation matrices to load.
  // mode vector
  std::vector<unsigned int> compare_coefficients;
  ar                       &compare_coefficients;
  Assert(compare_coefficients == n_coefficients_per_direction,
         ExcMessage("A different number of coefficients vector has been used "
                    "to generate the transformation matrices you are about "
                    "to load!"));

  // finite element collection
  unsigned int size;
  ar          &size;
  AssertDimension(size, fe_collection->size());
  std::string name;
  for (unsigned int i = 0; i < size; ++i)
    {
      ar &name;
      Assert(name.compare((*fe_collection)[i].get_name()) == 0,
             ExcMessage("A different FECollection has been used to generate "
                        "the transformation matrices you are about to load!"));
    }

  // quadrature collection
  ar &size;
  AssertDimension(size, q_collection.size());
  Quadrature<dim> quadrature;
  for (unsigned int i = 0; i < size; ++i)
    {
      ar &quadrature;
      Assert(quadrature == q_collection[i],
             ExcMessage("A different QCollection has been used to generate "
                        "the transformation matrices you are about to load!"));
    }

  // Restore the transform matrices since all prerequisites are fulfilled.
  ar &fourier_transform_matrices;
}



template <int dim, int spacedim>
template <class Archive>
inline void
FESeries::Legendre<dim, spacedim>::save_transformation_matrices(
  Archive &ar,
  const unsigned int /*version*/)
{
  // Store information about those resources which have been used to generate
  // the transformation matrices.
  // mode vector
  ar &n_coefficients_per_direction;

  // finite element collection
  unsigned int size = fe_collection->size();
  ar          &size;
  for (unsigned int i = 0; i < size; ++i)
    ar &(*fe_collection)[i].get_name();

  // quadrature collection
  size = q_collection.size();
  ar &size;
  for (unsigned int i = 0; i < size; ++i)
    ar &q_collection[i];

  // Store the actual transform matrices.
  ar &legendre_transform_matrices;
}



template <int dim, int spacedim>
template <class Archive>
inline void
FESeries::Legendre<dim, spacedim>::load_transformation_matrices(
  Archive &ar,
  const unsigned int /*version*/)
{
  // Check whether the currently registered resources are compatible with
  // the transformation matrices to load.
  // mode vector
  std::vector<unsigned int> compare_coefficients;
  ar                       &compare_coefficients;
  Assert(compare_coefficients == n_coefficients_per_direction,
         ExcMessage("A different number of coefficients vector has been used "
                    "to generate the transformation matrices you are about "
                    "to load!"));

  // finite element collection
  unsigned int size;
  ar          &size;
  AssertDimension(size, fe_collection->size());
  std::string name;
  for (unsigned int i = 0; i < size; ++i)
    {
      ar &name;
      Assert(name.compare((*fe_collection)[i].get_name()) == 0,
             ExcMessage("A different FECollection has been used to generate "
                        "the transformation matrices you are about to load!"));
    }

  // quadrature collection
  ar &size;
  AssertDimension(size, q_collection.size());
  Quadrature<dim> quadrature;
  for (unsigned int i = 0; i < size; ++i)
    {
      ar &quadrature;
      Assert(quadrature == q_collection[i],
             ExcMessage("A different QCollection has been used to generate "
                        "the transformation matrices you are about to load!"));
    }

  // Restore the transform matrices since all prerequisites are fulfilled.
  ar &legendre_transform_matrices;
}


#endif // DOXYGEN

DEAL_II_NAMESPACE_CLOSE

#endif // dealii_fe_series_h
