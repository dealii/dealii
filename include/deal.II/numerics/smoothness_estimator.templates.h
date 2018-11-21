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

#ifndef dealii_smoothness_estimator_templates_h
#define dealii_smoothness_estimator_templates_h

#include <deal.II/base/quadrature.h>

#include <deal.II/fe/fe_series.h>

#include <deal.II/hp/q_collection.h>

#include <deal.II/numerics/smoothness_estimator.h>

#include <algorithm>
#include <cmath>
#include <functional>
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
      for (unsigned int d = 0; d < dim; d++)
        size[d] = N;
      coeff.reinit(size);
    }



    /**
     * Calculates predicates of @p ind in the form
     * \f$
     *    v = \sum\limits_{d=0}^{dim} ind[d]^2
     * \f$.
     *
     * We flag the predicate whether it fulfills the criterion
     * \f$
     *    0 < v < max_degree^2
     * \f$
     * using @p max_degree.
     */
    template <int dim>
    std::pair<bool, unsigned int>
    predicate(const TableIndices<dim> &ind, const unsigned int max_degree)
    {
      unsigned int v = 0;
      for (unsigned int i = 0; i < dim; i++)
        v += ind[i] * ind[i];
      if (v > 0 && v < max_degree * max_degree)
        return std::make_pair(true, v);
      else
        return std::make_pair(false, v);
    }
  } // namespace



  template <typename FESeriesType, typename DoFHandlerType, typename VectorType>
  void
  estimate_by_coeff_decay(
    FESeriesType &                         fe_series,
    const DoFHandlerType &                 dof_handler,
    const std::vector<const VectorType *> &all_solutions,
    const std::vector<Vector<float> *> &   all_smoothness_indicators,
    const VectorTools::NormType            regression_strategy)
  {
    AssertDimension(all_solutions.size(), all_smoothness_indicators.size());

    for (auto &smoothness_indicator : all_smoothness_indicators)
      smoothness_indicator->reinit(
        dof_handler.get_triangulation().n_active_cells());

    const unsigned int dim = DoFHandlerType::dimension;
    const unsigned int max_degree =
      dof_handler.get_fe_collection().max_degree();

    Table<dim, typename FESeriesType::CoefficientType> expansion_coefficients;
    resize(expansion_coefficients, max_degree);

    Vector<typename VectorType::value_type>                   local_dof_values;
    std::vector<double>                                       ln_k;
    std::pair<std::vector<unsigned int>, std::vector<double>> res;
    for (auto &cell : dof_handler.active_cell_iterators())
      if (cell->is_locally_owned())
        {
          local_dof_values.reinit(cell->get_fe().dofs_per_cell);

          auto solution_it              = all_solutions.cbegin();
          auto smoothness_indicators_it = all_smoothness_indicators.begin();
          for (; solution_it != all_solutions.cend();
               ++solution_it, ++smoothness_indicators_it)
            {
              // Inside the loop, we first need to get the values of the local
              // degrees of freedom and then need to compute the series
              // expansion by multiplying this vector with the matrix ${\cal F}$
              // corresponding to this finite element.
              cell->get_dof_values(*(*solution_it), local_dof_values);

              fe_series.calculate(local_dof_values,
                                  cell->active_fe_index(),
                                  expansion_coefficients);

              // We fit our exponential decay of expansion coefficients to the
              // provided regression_strategy on each possible value of |k|. To
              // this end, we use FESeries::process_coefficients() to rework
              // coefficients into the desired format.
              res = FESeries::process_coefficients<dim>(
                expansion_coefficients,
                std::bind(&predicate<dim>, std::placeholders::_1, max_degree),
                regression_strategy);

              Assert(res.first.size() == res.second.size(), ExcInternalError());

              // Prepare linear equation for the logarithmic least squares fit.
              //
              // First, calculate ln(|k|). This vector will be the same for all
              // the cells so we can calculate ln(|k|) only once.
              //
              // For Fourier expansion, this translates to
              // ln(2*pi*sqrt(predicate)) = ln(2*pi) + 0.5*ln(predicate). Since
              // we are just interested in a linear regression later, we omit
              // the ln(2*pi) factor.
              // For Legendre expansion, this translates to
              // 0.5*ln(predicate) as well, without the pi factor.
              if (ln_k.empty())
                {
                  ln_k.resize(res.first.size());
                  for (unsigned int f = 0; f < res.first.size(); ++f)
                    ln_k[f] = 0.5 * std::log((double)res.first[f]);
                }

              // Second, calculate ln(U_k).
              for (double &residual_element : res.second)
                residual_element = std::log(residual_element);

              // Last, do the linear regression.
              std::pair<double, double> fit =
                FESeries::linear_regression(ln_k, res.second);

              // Compute the Sobolev index s=mu-dim/2 and store it in the vector
              // of estimated values for each cell.
              (*(*smoothness_indicators_it))(cell->active_cell_index()) =
                (float)(-fit.first - .5 * dim);
            }
        }
  }



  template <typename FESeriesType, typename DoFHandlerType, typename VectorType>
  void
  estimate_by_coeff_decay(FESeriesType &              fe_series,
                          const DoFHandlerType &      dof_handler,
                          const VectorType &          solution,
                          Vector<float> &             smoothness_indicators,
                          const VectorTools::NormType regression_strategy)
  {
    const std::vector<const VectorType *> all_solutions(1, &solution);
    const std::vector<Vector<float> *>    all_smoothness_indicators(
      1, &smoothness_indicators);

    estimate_by_coeff_decay(fe_series,
                            dof_handler,
                            all_solutions,
                            all_smoothness_indicators,
                            regression_strategy);
  }



  template <typename FESeriesType, typename DoFHandlerType, typename VectorType>
  void
  estimate_by_coeff_decay(
    const DoFHandlerType &                 dof_handler,
    const std::vector<const VectorType *> &all_solutions,
    const std::vector<Vector<float> *> &   all_smoothness_indicators,
    const VectorTools::NormType            regression_strategy)
  {
    const unsigned int dim = DoFHandlerType::dimension;
    const unsigned int max_degree =
      dof_handler.get_fe_collection().max_degree();

    // We initialize a series expansion object object which will be used to
    // calculate the expansion coefficients. In addition to the
    // hp::FECollection, we need to provide quadrature rules hp::QCollection for
    // integration on the reference cell.
    // We will need to assemble the expansion matrices for each of the finite
    // elements we deal with, i.e. the matrices F_k,j. We have to do that for
    // each of the finite elements in use. To that end we need a quadrature
    // rule. As a default, we use the same quadrature formula for each finite
    // element, namely one that is obtained by iterating a 2-point Gauss formula
    // as many times as the maximal exponent we use for the term exp(ikx).
    QGauss<1>      base_quadrature(2);
    QIterated<dim> quadrature(base_quadrature, max_degree);

    hp::QCollection<dim> expansion_q_collection;
    for (unsigned int i = 0; i < dof_handler.get_fe_collection().size(); ++i)
      expansion_q_collection.push_back(quadrature);

    // The FESeries::Fourier class' constructor first parameter $N$ defines the
    // number of coefficients in 1D with the total number of coefficients being
    // $N^{dim}$.
    FESeriesType fe_series(max_degree,
                           dof_handler.get_fe_collection(),
                           expansion_q_collection);

    estimate_by_coeff_decay(fe_series,
                            dof_handler,
                            all_solutions,
                            all_smoothness_indicators,
                            regression_strategy);
  }



  template <typename FESeriesType, typename DoFHandlerType, typename VectorType>
  void
  estimate_by_coeff_decay(const DoFHandlerType &      dof_handler,
                          const VectorType &          solution,
                          Vector<float> &             smoothness_indicators,
                          const VectorTools::NormType regression_strategy)
  {
    const std::vector<const VectorType *> all_solutions(1, &solution);
    const std::vector<Vector<float> *>    all_smoothness_indicators(
      1, &smoothness_indicators);

    estimate_by_coeff_decay<FESeriesType>(dof_handler,
                                          all_solutions,
                                          all_smoothness_indicators,
                                          regression_strategy);
  }
} // namespace SmoothnessEstimator


DEAL_II_NAMESPACE_CLOSE

#endif
