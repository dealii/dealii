// ---------------------------------------------------------------------
//
// Copyright (C) 2016 by the deal.II authors
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

#ifndef dealii__fe_series_H
#define dealii__fe_series_H



#include <deal.II/base/std_cxx11/shared_ptr.h>
#include <deal.II/base/config.h>
#include <deal.II/base/subscriptor.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/table.h>
#include <deal.II/base/table_indices.h>
#include <deal.II/base/tensor.h>
#include <deal.II/hp/fe_collection.h>
#include <deal.II/hp/q_collection.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/vector.h>
#include <deal.II/numerics/vector_tools.h>

#include <vector>
#include <string>


DEAL_II_NAMESPACE_OPEN


/*!@addtogroup feall */
/*@{*/


/**
 * This namespace offers functions to calculate expansion series of the
 * solution on the reference element. Coefficients of expansion are often used
 * to estimate local smoothness of the underlying FiniteElement field to decide
 * on h- or p-adaptive refinement strategy.
 *
 * @author Denis Davydov, 2016;
 */
namespace FESeries
{
  /**
   * A class to calculate expansion of a scalar FE field into Fourier series
   * on a reference element. The exponential form of the Fourier series is
   * based on completeness and Hermitian orthogonality of the set of exponential
   * functions \f$ \phi_{\bf k}({\bf x}) = \exp(2 \pi i\, {\bf k} \cdot {\bf x})\f$.
   * For example in 1D the L2-orthogonality condition reads
   * @f[
   *   \int_0^1 \phi_k(x) \phi_l^\ast(x) dx=\delta_{kl}.
   * @f]
   * Note that \f$ \phi_{\bf k} = \phi_{-\bf k}^\ast \f$.
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
   * this class. Note that \f$ u({\bf x}) = \sum_i u_i N_i({\bf x})\f$,
   * where \f$ N_i({\bf x}) \f$ are real-valued FiniteElement shape functions.
   * Consequently \f$ c_{\bf k} \equiv c_{-\bf k}^\ast \f$ and
   * we only need to compute \f$ c_{\bf k} \f$ for positive indices
   * \f$ \bf k \f$ .
   *
   * @author Denis Davydov, 2016.
   */
  template <int dim>
  class Fourier : public Subscriptor
  {
  public:
    /**
     * A non-default constructor. The @p size_in_each_direction defines the number
     * of modes in each direction, @p fe_collection is the hp::FECollection
     * for which expansion will be used and @p q_collection is the hp::QCollection
     * used to integrate the expansion for each FiniteElement
     * in @p fe_collection.
     */
    Fourier(const unsigned int size_in_each_direction,
            const hp::FECollection<dim> &fe_collection,
            const hp::QCollection<dim> &q_collection);

    /**
     * Calculate @p fourier_coefficients of the cell vector field given by
     * @p local_dof_values corresponding to FiniteElement with
     * @p cell_active_fe_index .
     */
    void calculate(const dealii::Vector<double>     &local_dof_values,
                   const unsigned int                cell_active_fe_index,
                   Table<dim,std::complex<double> > &fourier_coefficients);

  private:
    /**
     * hp::FECollection for which transformation matrices will be calculated.
     */
    SmartPointer<const hp::FECollection<dim> > fe_collection;

    /**
     * hp::QCollection used in calculation of transformation matrices.
     */
    SmartPointer<const hp::QCollection<dim> >  q_collection;

    /**
     * Ensure that the transformation matrix for FiniteElement index
     * @p fe_index is calculated. If not, calculate it.
     */
    void ensure_existence(const unsigned int fe_index);

    /**
     * Angular frequencies \f$ 2 \pi {\bf k} \f$ .
     */
    Table<dim, Tensor<1,dim> > k_vectors;

    /**
     * Transformation matrices for each FiniteElement.
     */
    std::vector<FullMatrix<std::complex<double> > > fourier_transform_matrices;

    /**
     * Auxiliary vector to store unrolled coefficients.
     */
    std::vector<std::complex<double> > unrolled_coefficients;

  };

  /**
   * A class to calculate expansion of a scalar FE field into series of Legendre
   * functions on a reference element.
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
   * These polynomials are orthogonal with respect to the \f$ L^2 \f$ inner
   * product on the interval \f$ [-1;1] \f$
   * @f[
   *    \int_{-1}^1 P_m(x) P_n(x) = \frac{2}{2n + 1} \delta_{mn}
   * @f]
   * and are complete.
   * A family of \f$ L^2 \f$-orthogonal polynomials on \f$ [0;1] \f$ can be
   * constructed via
   * @f[
   *    \widetilde P_m = \sqrt{2} P_m(2x-1).
   * @f]
   *
   *
   * An arbitrary scalar FE field on the reference element \f$ [0;1] \f$ can be
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
   * This class calculates coefficients \f$ c_{\bf k} \f$ using
   * \f$ dim \f$-dimensional Legendre polynomials constructed from
   * \f$ \widetilde P_m(x) \f$ using tensor product rule.
   *
   * @author Denis Davydov, 2016.
   */
  template <int dim>
  class Legendre : public Subscriptor
  {
  public:
    /**
     * A non-default constructor. The @p size_in_each_direction defines the number
     * of coefficients in each direction, @p fe_collection is the hp::FECollection
     * for which expansion will be used and @p q_collection is the hp::QCollection
     * used to integrate the expansion for each FiniteElement
     * in @p fe_collection.
     */
    Legendre(const unsigned int size_in_each_direction,
             const hp::FECollection<dim> &fe_collection,
             const hp::QCollection<dim> &q_collection);

    /**
     * Calculate @p legendre_coefficients of the cell vector field given by
     * @p local_dof_values corresponding to FiniteElement with
     * @p cell_active_fe_index .
     */
    void calculate(const dealii::Vector<double> &local_dof_values,
                   const unsigned int            cell_active_fe_index,
                   Table<dim,double>            &legendre_coefficients);

  private:
    /**
     * Number of coefficients in each direction
     */
    const unsigned int N;

    /**
     * hp::FECollection for which transformation matrices will be calculated.
     */
    SmartPointer<const hp::FECollection<dim> > fe_collection;

    /**
     * hp::QCollection used in calculation of transformation matrices.
     */
    SmartPointer<const hp::QCollection<dim> > q_collection;

    /**
     * Ensure that the transformation matrix for FiniteElement index
     * @p fe_index is calculated. If not, calculate it.
     */
    void ensure_existence(const unsigned int fe_index);

    /**
     * Transformation matrices for each FiniteElement.
     */
    std::vector<FullMatrix<double> > legendre_transform_matrices;

    /**
     * Auxiliary vector to store unrolled coefficients.
     */
    std::vector<double> unrolled_coefficients;

  };


  /**
   * Calculate the @p norm of subsets of @p coefficients defined by
   * @p predicate being constant. Return the pair of vectors of predicate values
   * and the vector of calculated subset norms.
   *
   * @p predicate should return a pair of <code>bool</code> and <code>unsigned int</code>.
   * The former is a flag whether a given TableIndices should be used in
   * calculation, whereas the latter is the unrolled value of indices according
   * to which the subsets of coefficients will be formed.
   *
   * @note Only the following values of @p norm are implemented and make sense
   * in this case: mean, L1_norm, L2_norm, Linfty_norm. The mean norm can only
   * be applied to real valued coefficients.
   */
  template <int dim, typename T>
  std::pair<std::vector<unsigned int>,std::vector<double> >
  process_coefficients(const Table<dim,T> &coefficients,
                       const std_cxx11::function<std::pair<bool,unsigned int>(const TableIndices<dim> &)> &predicate,
                       const VectorTools::NormType norm);



  /**
   * Linear regression least-square fit of $y = k \, x + b$.
   * The size of the input vectors should be equal and more than 1.
   * The returned pair will contain $k$ (first) and $b$ (second).
   */
  std::pair<double,double> linear_regression(const std::vector<double> &x,
                                             const std::vector<double> &y);

}

/*@}*/

#ifndef DOXYGEN

// -------------------  inline and template functions ----------------

namespace
{
  template <int dim,typename T>
  void fill_map_index(const Table<dim,T> &coefficients,
                      const TableIndices<dim> &ind,
                      const std_cxx11::function<std::pair<bool,unsigned int>(const TableIndices<dim> &)> &predicate,
                      std::map<unsigned int, std::vector<T> > &pred_to_values)
  {
    const std::pair<bool,unsigned int> pred_pair = predicate(ind);
    // don't add a value if predicate is false
    if (pred_pair.first == false)
      return;

    const unsigned int &pred_value = pred_pair.second;
    const T &coeff_value = coefficients(ind);
    // If pred_value is not in the pred_to_values map, the element will be created.
    // Otherwise a reference to the existing element is returned.
    pred_to_values[pred_value].push_back(coeff_value);
  }

  template <typename T>
  void fill_map(const Table<1,T> &coefficients,
                const std_cxx11::function<std::pair<bool,unsigned int>(const TableIndices<1> &)> &predicate,
                std::map<unsigned int, std::vector<T> > &pred_to_values)
  {
    for (unsigned int i = 0; i < coefficients.size(0); i++)
      {
        const TableIndices<1> ind(i);
        fill_map_index(coefficients,ind,predicate,pred_to_values);
      }

  }

  template <typename T>
  void fill_map(const Table<2,T> &coefficients,
                const std_cxx11::function<std::pair<bool,unsigned int>(const TableIndices<2> &)> &predicate,
                std::map<unsigned int, std::vector<T> > &pred_to_values)
  {
    for (unsigned int i = 0; i < coefficients.size(0); i++)
      for (unsigned int j = 0; j < coefficients.size(1); j++)
        {
          const TableIndices<2> ind(i,j);
          fill_map_index(coefficients,ind,predicate,pred_to_values);
        }

  }

  template <typename T>
  void fill_map(const Table<3,T> &coefficients,
                const std_cxx11::function<std::pair<bool,unsigned int>(const TableIndices<3> &)> &predicate,
                std::map<unsigned int, std::vector<T> > &pred_to_values)
  {
    for (unsigned int i = 0; i < coefficients.size(0); i++)
      for (unsigned int j = 0; j < coefficients.size(1); j++)
        for (unsigned int k = 0; k < coefficients.size(2); k++)
          {
            const TableIndices<3> ind(i,j,k);
            fill_map_index(coefficients,ind,predicate,pred_to_values);
          }
  }


  template <typename T>
  double complex_mean_value(const T &value)
  {
    return value;
  }

  template <typename T>
  double complex_mean_value(const std::complex<T> &value)
  {
    AssertThrow(false, ExcMessage("FESeries::process_coefficients() can not be used with"
                                  "complex-valued coefficients and VectorTools::mean norm."));
    return std::abs(value);
  }

}


template <int dim, typename T>
std::pair<std::vector<unsigned int>,std::vector<double> >
FESeries::process_coefficients(const Table<dim,T> &coefficients,
                               const std_cxx11::function<std::pair<bool,unsigned int>(const TableIndices<dim> &)> &predicate,
                               const VectorTools::NormType norm)
{
  std::vector<unsigned int> predicate_values;
  std::vector<double> norm_values;

  // first, parse all table elements into a map of predicate values and coefficients.
  // We could have stored (predicate values ->TableIndicies) map, but its
  // processing would have been much harder later on.
  std::map<unsigned int, std::vector<T> > pred_to_values;
  fill_map(coefficients,predicate,pred_to_values);

  // now go through the map and populate the @p norm_values based on @p norm:
  for (typename std::map<unsigned int, std::vector<T> >::const_iterator it = pred_to_values.begin();
       it != pred_to_values.end(); it++)
    {
      predicate_values.push_back(it->first);
      Vector<T> values(it->second.begin(),
                       it->second.end());

      switch (norm)
        {
        case VectorTools::L2_norm:
        {
          norm_values.push_back(values.l2_norm());
          break;
        }
        case VectorTools::L1_norm:
        {
          norm_values.push_back(values.l1_norm());
          break;
        }
        case VectorTools::Linfty_norm:
        {
          norm_values.push_back(values.linfty_norm());
          break;
        }
        case VectorTools::mean:
        {
          norm_values.push_back(complex_mean_value(values.mean_value()));
          break;
        }
        default:
          AssertThrow(false, ExcNotImplemented());
          break;
        }
    }

  return std::make_pair(predicate_values,norm_values);
}


#endif // DOXYGEN

DEAL_II_NAMESPACE_CLOSE

#endif // dealii__fe_series_H
