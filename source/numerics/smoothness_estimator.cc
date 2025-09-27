// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2018 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/signaling_nan.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_series.h>

#include <deal.II/grid/filtered_iterator.h>

#include <deal.II/hp/q_collection.h>

#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/la_parallel_block_vector.h>
#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/lac/petsc_block_vector.h>
#include <deal.II/lac/petsc_vector.h>
#include <deal.II/lac/trilinos_epetra_vector.h>
#include <deal.II/lac/trilinos_parallel_block_vector.h>
#include <deal.II/lac/trilinos_vector.h>
#include <deal.II/lac/vector.h>

#include <deal.II/numerics/smoothness_estimator.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <utility>


DEAL_II_NAMESPACE_OPEN


namespace SmoothnessEstimator
{
  namespace
  {
    /**
     * Resizes @p coeff to @p N in each dimension.
     */
    template <int dim, typename CoefficientType>
    void
    resize(Table<dim, CoefficientType> &coeff, const unsigned int N)
    {
      TableIndices<dim> size;
      for (unsigned int d = 0; d < dim; ++d)
        size[d] = N;
      coeff.reinit(size);
    }
  } // namespace



  namespace Legendre
  {
    namespace
    {
      /**
       * We will need to take the maximum absolute value of Legendre
       * coefficients which correspond to $k$-vector $|{\bf k}|= const$. To
       * filter the coefficients Table we will use the
       * FESeries::process_coefficients() which requires a predicate to be
       * specified. The predicate should operate on TableIndices and return a
       * pair of <code>bool</code> and <code>unsigned int</code>. The latter is
       * the value of the map from TableIndices to unsigned int.  It is used to
       * define subsets of coefficients from which we search for the one with
       * highest absolute value, i.e. $l^\infty$-norm. The <code>bool</code>
       * parameter defines which indices should be used in processing. In the
       * current case we are interested in coefficients which correspond to $0
       * <= i+j < N$ and $0 <= i+j+k < N$ in 2d and 3d, respectively.
       */
      template <int dim>
      std::pair<bool, unsigned int>
      index_sum_less_than_N(const TableIndices<dim> &ind, const unsigned int N)
      {
        unsigned int v = 0;
        for (unsigned int i = 0; i < dim; ++i)
          v += ind[i];

        return std::make_pair((v < N), v);
      }
    } // namespace



    template <int dim, int spacedim, typename VectorType>
    void
    coefficient_decay(FESeries::Legendre<dim, spacedim> &fe_legendre,
                      const DoFHandler<dim, spacedim>   &dof_handler,
                      const VectorType                  &solution,
                      Vector<float>                     &smoothness_indicators,
                      const VectorTools::NormType        regression_strategy,
                      const double smallest_abs_coefficient,
                      const bool   only_flagged_cells)
    {
      using number = typename VectorType::value_type;
      using number_coeff =
        typename FESeries::Legendre<dim, spacedim>::CoefficientType;

      smoothness_indicators.reinit(
        dof_handler.get_triangulation().n_active_cells());

      unsigned int             n_modes;
      Table<dim, number_coeff> expansion_coefficients;

      Vector<number>      local_dof_values;
      std::vector<double> converted_indices;
      std::pair<std::vector<unsigned int>, std::vector<double>> res;
      for (const auto &cell : dof_handler.active_cell_iterators() |
                                IteratorFilters::LocallyOwnedCell())
        {
          if (!only_flagged_cells || cell->refine_flag_set() ||
              cell->coarsen_flag_set())
            {
              n_modes = fe_legendre.get_n_coefficients_per_direction(
                cell->active_fe_index());
              resize(expansion_coefficients, n_modes);

              local_dof_values.reinit(cell->get_fe().n_dofs_per_cell());
              cell->get_dof_values(solution, local_dof_values);

              fe_legendre.calculate(local_dof_values,
                                    cell->active_fe_index(),
                                    expansion_coefficients);

              // We fit our exponential decay of expansion coefficients to the
              // provided regression_strategy on each possible value of |k|.
              // To this end, we use FESeries::process_coefficients() to
              // rework coefficients into the desired format.
              res = FESeries::process_coefficients<dim>(
                expansion_coefficients,
                [n_modes](const TableIndices<dim> &indices) {
                  return index_sum_less_than_N(indices, n_modes);
                },
                regression_strategy,
                smallest_abs_coefficient);

              Assert(res.first.size() == res.second.size(), ExcInternalError());

              // Last, do the linear regression.
              float regularity = std::numeric_limits<float>::infinity();
              if (res.first.size() > 1)
                {
                  // Prepare linear equation for the logarithmic least squares
                  // fit.
                  converted_indices.assign(res.first.begin(), res.first.end());

                  for (auto &residual_element : res.second)
                    residual_element = std::log(residual_element);

                  const std::pair<double, double> fit =
                    FESeries::linear_regression(converted_indices, res.second);
                  regularity = static_cast<float>(-fit.first);
                }

              smoothness_indicators(cell->active_cell_index()) = regularity;
            }
          else
            smoothness_indicators(cell->active_cell_index()) =
              numbers::signaling_nan<float>();
        }
    }



    template <int dim, int spacedim, typename VectorType>
    void
    coefficient_decay_per_direction(
      FESeries::Legendre<dim, spacedim> &fe_legendre,
      const DoFHandler<dim, spacedim>   &dof_handler,
      const VectorType                  &solution,
      Vector<float>                     &smoothness_indicators,
      const ComponentMask               &coefficients_predicate,
      const double                       smallest_abs_coefficient,
      const bool                         only_flagged_cells)
    {
      Assert(smallest_abs_coefficient >= 0.,
             ExcMessage("smallest_abs_coefficient should be non-negative."));

      using number = typename VectorType::value_type;
      using number_coeff =
        typename FESeries::Legendre<dim, spacedim>::CoefficientType;

      smoothness_indicators.reinit(
        dof_handler.get_triangulation().n_active_cells());

      unsigned int             n_modes;
      Table<dim, number_coeff> expansion_coefficients;
      Vector<number>           local_dof_values;

      // auxiliary vector to do linear regression
      const unsigned int max_degree =
        dof_handler.get_fe_collection().max_degree();

      std::vector<double> x, y;
      x.reserve(max_degree);
      y.reserve(max_degree);

      for (const auto &cell : dof_handler.active_cell_iterators() |
                                IteratorFilters::LocallyOwnedCell())
        {
          if (!only_flagged_cells || cell->refine_flag_set() ||
              cell->coarsen_flag_set())
            {
              n_modes = fe_legendre.get_n_coefficients_per_direction(
                cell->active_fe_index());
              resize(expansion_coefficients, n_modes);

              const unsigned int pe = cell->get_fe().degree;
              Assert(pe > 0, ExcInternalError());

              // since we use coefficients with indices [1,pe] in each
              // direction, the number of coefficients we need to calculate is
              // at least N=pe+1
              AssertIndexRange(pe, n_modes);

              local_dof_values.reinit(cell->get_fe().n_dofs_per_cell());
              cell->get_dof_values(solution, local_dof_values);

              fe_legendre.calculate(local_dof_values,
                                    cell->active_fe_index(),
                                    expansion_coefficients);

              // choose the smallest decay of coefficients in each direction,
              // i.e. the maximum decay slope k_v as in exp(-k_v)
              double k_v = std::numeric_limits<double>::max();
              for (unsigned int d = 0; d < dim; ++d)
                {
                  x.resize(0);
                  y.resize(0);

                  // will use all non-zero coefficients allowed by the
                  // predicate function
                  for (unsigned int i = 0; i <= pe; ++i)
                    if (coefficients_predicate[i])
                      {
                        TableIndices<dim> ind;
                        ind[d] = i;
                        const double coeff_abs =
                          std::abs(expansion_coefficients(ind));

                        if (coeff_abs > smallest_abs_coefficient)
                          {
                            x.push_back(i);
                            y.push_back(std::log(coeff_abs));
                          }
                      }

                  // in case we don't have enough non-zero coefficient to fit,
                  // skip this direction
                  if (x.size() < 2)
                    continue;

                  const std::pair<double, double> fit =
                    FESeries::linear_regression(x, y);

                  // decay corresponds to negative slope
                  // take the lesser negative slope along each direction
                  k_v = std::min(k_v, -fit.first);
                }

              smoothness_indicators(cell->active_cell_index()) =
                static_cast<float>(k_v);
            }
          else
            smoothness_indicators(cell->active_cell_index()) =
              numbers::signaling_nan<float>();
        }
    }



    template <int dim, int spacedim>
    FESeries::Legendre<dim, spacedim>
    default_fe_series(const hp::FECollection<dim, spacedim> &fe_collection,
                      const unsigned int                     component)
    {
      // Default number of coefficients per direction.
      //
      // With a number of modes equal to the polynomial degree plus two for each
      // finite element, the smoothness estimation algorithm tends to produce
      // stable results.
      std::vector<unsigned int> n_coefficients_per_direction;
      n_coefficients_per_direction.reserve(fe_collection.size());
      for (unsigned int i = 0; i < fe_collection.size(); ++i)
        n_coefficients_per_direction.push_back(fe_collection[i].degree + 2);

      // Default quadrature collection.
      //
      // We initialize a FESeries::Legendre expansion object object which will
      // be used to calculate the expansion coefficients. In addition to the
      // hp::FECollection, we need to provide quadrature rules hp::QCollection
      // for integration on the reference cell.
      // We will need to assemble the expansion matrices for each of the finite
      // elements we deal with, i.e. the matrices F_k,j. We have to do that for
      // each of the finite elements in use. To that end we need a quadrature
      // rule. As a default, we use the same quadrature formula for each finite
      // element, namely a Gauss formula that yields exact results for the
      // highest order Legendre polynomial used.
      //
      // We start with the zeroth Legendre polynomial which is just a constant,
      // so the highest Legendre polynomial will be of order (n_modes - 1).
      hp::QCollection<dim> q_collection;
      for (unsigned int i = 0; i < fe_collection.size(); ++i)
        {
          const QGauss<dim>  quadrature(n_coefficients_per_direction[i]);
          const QSorted<dim> quadrature_sorted(quadrature);
          q_collection.push_back(quadrature_sorted);
        }

      return FESeries::Legendre<dim, spacedim>(n_coefficients_per_direction,
                                               fe_collection,
                                               q_collection,
                                               component);
    }
  } // namespace Legendre



  namespace Fourier
  {
    namespace
    {
      /**
       * We will need to take the maximum absolute value of Fourier coefficients
       * which correspond to $k$-vector $|{\bf k}|= const$. To filter the
       * coefficients Table we will use the FESeries::process_coefficients()
       * which requires a predicate to be specified. The predicate should
       * operate on TableIndices and return a pair of <code>bool</code> and
       * <code>unsigned int</code>. The latter is the value of the map from
       * TableIndices to unsigned int.  It is used to define subsets of
       * coefficients from which we search for the one with highest absolute
       * value, i.e. $l^\infty$-norm. The <code>bool</code> parameter defines
       * which indices should be used in processing. In the current case we are
       * interested in coefficients which correspond to $0 < i^2+j^2 < N^2$ and
       * $0 < i^2+j^2+k^2 < N^2$ in 2d and 3d, respectively.
       */
      template <int dim>
      std::pair<bool, unsigned int>
      index_norm_greater_than_zero_and_less_than_N_squared(
        const TableIndices<dim> &ind,
        const unsigned int       N)
      {
        unsigned int v = 0;
        for (unsigned int i = 0; i < dim; ++i)
          v += ind[i] * ind[i];

        return std::make_pair((v > 0 && v < N * N), v);
      }
    } // namespace



    template <int dim, int spacedim, typename VectorType>
    void
    coefficient_decay(FESeries::Fourier<dim, spacedim> &fe_fourier,
                      const DoFHandler<dim, spacedim>  &dof_handler,
                      const VectorType                 &solution,
                      Vector<float>                    &smoothness_indicators,
                      const VectorTools::NormType       regression_strategy,
                      const double smallest_abs_coefficient,
                      const bool   only_flagged_cells)
    {
      using number = typename VectorType::value_type;
      using number_coeff =
        typename FESeries::Fourier<dim, spacedim>::CoefficientType;

      smoothness_indicators.reinit(
        dof_handler.get_triangulation().n_active_cells());

      unsigned int             n_modes;
      Table<dim, number_coeff> expansion_coefficients;

      Vector<number>      local_dof_values;
      std::vector<double> ln_k;
      std::pair<std::vector<unsigned int>, std::vector<double>> res;
      for (const auto &cell : dof_handler.active_cell_iterators() |
                                IteratorFilters::LocallyOwnedCell())
        {
          if (!only_flagged_cells || cell->refine_flag_set() ||
              cell->coarsen_flag_set())
            {
              n_modes = fe_fourier.get_n_coefficients_per_direction(
                cell->active_fe_index());
              resize(expansion_coefficients, n_modes);

              // Inside the loop, we first need to get the values of the local
              // degrees of freedom and then need to compute the series
              // expansion by multiplying this vector with the matrix ${\cal
              // F}$ corresponding to this finite element.
              local_dof_values.reinit(cell->get_fe().n_dofs_per_cell());
              cell->get_dof_values(solution, local_dof_values);

              fe_fourier.calculate(local_dof_values,
                                   cell->active_fe_index(),
                                   expansion_coefficients);

              // We fit our exponential decay of expansion coefficients to the
              // provided regression_strategy on each possible value of |k|.
              // To this end, we use FESeries::process_coefficients() to
              // rework coefficients into the desired format.
              res = FESeries::process_coefficients<dim>(
                expansion_coefficients,
                [n_modes](const TableIndices<dim> &indices) {
                  return index_norm_greater_than_zero_and_less_than_N_squared(
                    indices, n_modes);
                },
                regression_strategy,
                smallest_abs_coefficient);

              Assert(res.first.size() == res.second.size(), ExcInternalError());

              // Last, do the linear regression.
              float regularity = std::numeric_limits<float>::infinity();
              if (res.first.size() > 1)
                {
                  // Prepare linear equation for the logarithmic least squares
                  // fit.
                  //
                  // First, calculate ln(|k|).
                  //
                  // For Fourier expansion, this translates to
                  // ln(2*pi*sqrt(predicate)) = ln(2*pi) + 0.5*ln(predicate).
                  // Since we are just interested in the slope of a linear
                  // regression later, we omit the ln(2*pi) factor.
                  ln_k.resize(res.first.size());
                  for (unsigned int f = 0; f < res.first.size(); ++f)
                    ln_k[f] = 0.5 * std::log(static_cast<double>(res.first[f]));

                  // Second, calculate ln(U_k).
                  for (auto &residual_element : res.second)
                    residual_element = std::log(residual_element);

                  const std::pair<double, double> fit =
                    FESeries::linear_regression(ln_k, res.second);
                  // Compute regularity s = mu - dim/2
                  regularity = static_cast<float>(-fit.first) -
                               ((dim > 1) ? (.5 * dim) : 0);
                }

              // Store result in the vector of estimated values for each cell.
              smoothness_indicators(cell->active_cell_index()) = regularity;
            }
          else
            smoothness_indicators(cell->active_cell_index()) =
              numbers::signaling_nan<float>();
        }
    }



    template <int dim, int spacedim, typename VectorType>
    void
    coefficient_decay_per_direction(
      FESeries::Fourier<dim, spacedim> &fe_fourier,
      const DoFHandler<dim, spacedim>  &dof_handler,
      const VectorType                 &solution,
      Vector<float>                    &smoothness_indicators,
      const ComponentMask              &coefficients_predicate,
      const double                      smallest_abs_coefficient,
      const bool                        only_flagged_cells)
    {
      Assert(smallest_abs_coefficient >= 0.,
             ExcMessage("smallest_abs_coefficient should be non-negative."));

      using number = typename VectorType::value_type;
      using number_coeff =
        typename FESeries::Fourier<dim, spacedim>::CoefficientType;

      smoothness_indicators.reinit(
        dof_handler.get_triangulation().n_active_cells());

      unsigned int             n_modes;
      Table<dim, number_coeff> expansion_coefficients;
      Vector<number>           local_dof_values;

      // auxiliary vector to do linear regression
      const unsigned int max_degree =
        dof_handler.get_fe_collection().max_degree();

      std::vector<double> x, y;
      x.reserve(max_degree);
      y.reserve(max_degree);

      for (const auto &cell : dof_handler.active_cell_iterators() |
                                IteratorFilters::LocallyOwnedCell())
        {
          if (!only_flagged_cells || cell->refine_flag_set() ||
              cell->coarsen_flag_set())
            {
              n_modes = fe_fourier.get_n_coefficients_per_direction(
                cell->active_fe_index());
              resize(expansion_coefficients, n_modes);

              const unsigned int pe = cell->get_fe().degree;
              Assert(pe > 0, ExcInternalError());

              // since we use coefficients with indices [1,pe] in each
              // direction, the number of coefficients we need to calculate is
              // at least N=pe+1
              AssertIndexRange(pe, n_modes);

              local_dof_values.reinit(cell->get_fe().n_dofs_per_cell());
              cell->get_dof_values(solution, local_dof_values);

              fe_fourier.calculate(local_dof_values,
                                   cell->active_fe_index(),
                                   expansion_coefficients);

              // choose the smallest decay of coefficients in each direction,
              // i.e. the maximum decay slope k_v as in exp(-k_v)
              double k_v = std::numeric_limits<double>::max();
              for (unsigned int d = 0; d < dim; ++d)
                {
                  x.resize(0);
                  y.resize(0);

                  // will use all non-zero coefficients allowed by the
                  // predicate function
                  //
                  // skip i=0 because of logarithm
                  for (unsigned int i = 1; i <= pe; ++i)
                    if (coefficients_predicate[i])
                      {
                        TableIndices<dim> ind;
                        ind[d] = i;
                        const double coeff_abs =
                          std::abs(expansion_coefficients(ind));

                        if (coeff_abs > smallest_abs_coefficient)
                          {
                            x.push_back(std::log(i));
                            y.push_back(std::log(coeff_abs));
                          }
                      }

                  // in case we don't have enough non-zero coefficient to fit,
                  // skip this direction
                  if (x.size() < 2)
                    continue;

                  const std::pair<double, double> fit =
                    FESeries::linear_regression(x, y);

                  // decay corresponds to negative slope
                  // take the lesser negative slope along each direction
                  k_v = std::min(k_v, -fit.first);
                }

              smoothness_indicators(cell->active_cell_index()) =
                static_cast<float>(k_v);
            }
          else
            smoothness_indicators(cell->active_cell_index()) =
              numbers::signaling_nan<float>();
        }
    }



    template <int dim, int spacedim>
    FESeries::Fourier<dim, spacedim>
    default_fe_series(const hp::FECollection<dim, spacedim> &fe_collection,
                      const unsigned int                     component)
    {
      // Default number of coefficients per direction.
      //
      // Since we omit the zero-th mode in the Fourier decay strategy, make sure
      // that we have at least two modes to work with per finite element. With a
      // number of modes equal to the polynomial degree plus two for each finite
      // element, the smoothness estimation algorithm tends to produce stable
      // results.
      std::vector<unsigned int> n_coefficients_per_direction;
      n_coefficients_per_direction.reserve(fe_collection.size());
      for (unsigned int i = 0; i < fe_collection.size(); ++i)
        n_coefficients_per_direction.push_back(fe_collection[i].degree + 2);

      // Default quadrature collection.
      //
      // We initialize a series expansion object object which will be used to
      // calculate the expansion coefficients. In addition to the
      // hp::FECollection, we need to provide quadrature rules hp::QCollection
      // for integration on the reference cell.
      // We will need to assemble the expansion matrices for each of the finite
      // elements we deal with, i.e. the matrices F_k,j. We have to do that for
      // each of the finite elements in use. To that end we need a quadrature
      // rule. As a default, we use the same quadrature formula for each finite
      // element, namely one that is obtained by iterating a 5-point Gauss
      // formula as many times as the maximal exponent we use for the term
      // exp(ikx). Since the first mode corresponds to k = 0, the maximal wave
      // number is k = n_modes - 1.
      const QGauss<1>      base_quadrature(5);
      hp::QCollection<dim> q_collection;
      for (unsigned int i = 0; i < fe_collection.size(); ++i)
        {
          const QIterated<dim> quadrature(base_quadrature,
                                          n_coefficients_per_direction[i] - 1);
          const QSorted<dim>   quadrature_sorted(quadrature);
          q_collection.push_back(quadrature_sorted);
        }

      return FESeries::Fourier<dim, spacedim>(n_coefficients_per_direction,
                                              fe_collection,
                                              q_collection,
                                              component);
    }
  } // namespace Fourier
} // namespace SmoothnessEstimator


// explicit instantiations
#include "numerics/smoothness_estimator.inst"

DEAL_II_NAMESPACE_CLOSE
