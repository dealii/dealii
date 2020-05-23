// ---------------------------------------------------------------------
//
// Copyright (C) 2016 - 2020 by the deal.II authors
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



#include <deal.II/base/thread_management.h>

#include <deal.II/fe/fe_series.h>

#include <iostream>
#ifdef DEAL_II_WITH_GSL
#  include <gsl/gsl_sf_legendre.h>
#endif


DEAL_II_NAMESPACE_OPEN

namespace
{
  DeclException2(ExcLegendre,
                 int,
                 double,
                 << "x[" << arg1 << "] = " << arg2 << " is not in [0,1]");

  /*
   * dim dimensional Legendre function with indices @p indices
   * evaluated at @p x_q in [0,1]^dim.
   */
  template <int dim>
  double
  Lh(const Point<dim> &x_q, const TableIndices<dim> &indices)
  {
#ifdef DEAL_II_WITH_GSL
    double res = 1.0;
    for (unsigned int d = 0; d < dim; d++)
      {
        const double x = 2.0 * (x_q[d] - 0.5);
        Assert((x_q[d] <= 1.0) && (x_q[d] >= 0.), ExcLegendre(d, x_q[d]));
        const int ind = indices[d];
        res *= std::sqrt(2.0) * gsl_sf_legendre_Pl(ind, x);
      }
    return res;

#else

    (void)x_q;
    (void)indices;
    AssertThrow(false,
                ExcMessage("deal.II has to be configured with GSL "
                           "in order to use Legendre transformation."));
    return 0;
#endif
  }



  /*
   * Multiplier in Legendre coefficients
   */
  template <int dim>
  double
  multiplier(const TableIndices<dim> &indices)
  {
    double res = 1.0;
    for (unsigned int d = 0; d < dim; d++)
      res *= (0.5 + indices[d]);

    return res;
  }



  template <int dim, int spacedim>
  double
  integrate(const FiniteElement<dim, spacedim> &fe,
            const Quadrature<dim> &             quadrature,
            const TableIndices<dim> &           indices,
            const unsigned int                  dof)
  {
    double sum = 0;
    for (unsigned int q = 0; q < quadrature.size(); ++q)
      {
        const Point<dim> &x_q = quadrature.point(q);
        sum +=
          Lh(x_q, indices) * fe.shape_value(dof, x_q) * quadrature.weight(q);
      }
    return sum * multiplier(indices);
  }



  /**
   * Ensure that the transformation matrix for FiniteElement index
   * @p fe_index is calculated. If not, calculate it.
   */
  template <int spacedim>
  void
  ensure_existence(
    const std::vector<unsigned int> &    n_coefficients_per_direction,
    const hp::FECollection<1, spacedim> &fe_collection,
    const hp::QCollection<1> &           q_collection,
    const unsigned int                   fe,
    std::vector<FullMatrix<double>> &    legendre_transform_matrices)
  {
    AssertIndexRange(fe, fe_collection.size());

    if (legendre_transform_matrices[fe].m() == 0)
      {
        legendre_transform_matrices[fe].reinit(n_coefficients_per_direction[fe],
                                               fe_collection[fe].dofs_per_cell);

        for (unsigned int k = 0; k < n_coefficients_per_direction[fe]; ++k)
          for (unsigned int j = 0; j < fe_collection[fe].dofs_per_cell; ++j)
            legendre_transform_matrices[fe](k, j) = integrate(
              fe_collection[fe], q_collection[fe], TableIndices<1>(k), j);
      }
  }

  template <int spacedim>
  void
  ensure_existence(
    const std::vector<unsigned int> &    n_coefficients_per_direction,
    const hp::FECollection<2, spacedim> &fe_collection,
    const hp::QCollection<2> &           q_collection,
    const unsigned int                   fe,
    std::vector<FullMatrix<double>> &    legendre_transform_matrices)
  {
    AssertIndexRange(fe, fe_collection.size());

    if (legendre_transform_matrices[fe].m() == 0)
      {
        legendre_transform_matrices[fe].reinit(
          Utilities::fixed_power<2>(n_coefficients_per_direction[fe]),
          fe_collection[fe].dofs_per_cell);

        unsigned int k = 0;
        for (unsigned int k1 = 0; k1 < n_coefficients_per_direction[fe]; ++k1)
          for (unsigned int k2 = 0; k2 < n_coefficients_per_direction[fe];
               ++k2, k++)
            for (unsigned int j = 0; j < fe_collection[fe].dofs_per_cell; ++j)
              legendre_transform_matrices[fe](k, j) =
                integrate(fe_collection[fe],
                          q_collection[fe],
                          TableIndices<2>(k1, k2),
                          j);
      }
  }

  template <int spacedim>
  void
  ensure_existence(
    const std::vector<unsigned int> &    n_coefficients_per_direction,
    const hp::FECollection<3, spacedim> &fe_collection,
    const hp::QCollection<3> &           q_collection,
    const unsigned int                   fe,
    std::vector<FullMatrix<double>> &    legendre_transform_matrices)
  {
    AssertIndexRange(fe, fe_collection.size());

    if (legendre_transform_matrices[fe].m() == 0)
      {
        legendre_transform_matrices[fe].reinit(
          Utilities::fixed_power<3>(n_coefficients_per_direction[fe]),
          fe_collection[fe].dofs_per_cell);

        unsigned int k = 0;
        for (unsigned int k1 = 0; k1 < n_coefficients_per_direction[fe]; ++k1)
          for (unsigned int k2 = 0; k2 < n_coefficients_per_direction[fe]; ++k2)
            for (unsigned int k3 = 0; k3 < n_coefficients_per_direction[fe];
                 ++k3, k++)
              for (unsigned int j = 0; j < fe_collection[fe].dofs_per_cell; ++j)
                legendre_transform_matrices[fe](k, j) =
                  integrate(fe_collection[fe],
                            q_collection[fe],
                            TableIndices<3>(k1, k2, k3),
                            j);
      }
  }
} // namespace



namespace FESeries
{
  template <int dim, int spacedim>
  Legendre<dim, spacedim>::Legendre(
    const std::vector<unsigned int> &      n_coefficients_per_direction,
    const hp::FECollection<dim, spacedim> &fe_collection,
    const hp::QCollection<dim> &           q_collection)
    : n_coefficients_per_direction(n_coefficients_per_direction)
    , fe_collection(&fe_collection)
    , q_collection(q_collection)
    , legendre_transform_matrices(fe_collection.size())
  {
    Assert(n_coefficients_per_direction.size() == fe_collection.size() &&
             n_coefficients_per_direction.size() == q_collection.size(),
           ExcMessage("All parameters are supposed to have the same size."));

    // reserve sufficient memory
    const unsigned int max_n_coefficients_per_direction =
      *std::max_element(n_coefficients_per_direction.cbegin(),
                        n_coefficients_per_direction.cend());
    unrolled_coefficients.reserve(
      Utilities::fixed_power<dim>(max_n_coefficients_per_direction));
  }



  template <int dim, int spacedim>
  Legendre<dim, spacedim>::Legendre(
    const unsigned int                     n_coefficients_per_direction,
    const hp::FECollection<dim, spacedim> &fe_collection,
    const hp::QCollection<dim> &           q_collection)
    : Legendre<dim, spacedim>(
        std::vector<unsigned int>(fe_collection.size(),
                                  n_coefficients_per_direction),
        fe_collection,
        q_collection)
  {}



  template <int dim, int spacedim>
  inline bool
  Legendre<dim, spacedim>::
  operator==(const Legendre<dim, spacedim> &legendre) const
  {
    return (
      (n_coefficients_per_direction == legendre.n_coefficients_per_direction) &&
      (*fe_collection == *(legendre.fe_collection)) &&
      (q_collection == legendre.q_collection) &&
      (legendre_transform_matrices == legendre.legendre_transform_matrices));
  }



  template <int dim, int spacedim>
  void
  Legendre<dim, spacedim>::precalculate_all_transformation_matrices()
  {
    Threads::TaskGroup<> task_group;
    for (unsigned int fe = 0; fe < fe_collection->size(); ++fe)
      task_group += Threads::new_task([&, fe]() {
        ensure_existence(n_coefficients_per_direction,
                         *fe_collection,
                         q_collection,
                         fe,
                         legendre_transform_matrices);
      });

    task_group.join_all();
  }



  template <int dim, int spacedim>
  unsigned int
  Legendre<dim, spacedim>::get_n_coefficients_per_direction(
    const unsigned int index) const
  {
    return n_coefficients_per_direction[index];
  }



  template <int dim, int spacedim>
  template <typename Number>
  void
  Legendre<dim, spacedim>::calculate(
    const dealii::Vector<Number> &local_dof_values,
    const unsigned int            cell_active_fe_index,
    Table<dim, CoefficientType> & legendre_coefficients)
  {
    for (unsigned int d = 0; d < dim; ++d)
      AssertDimension(legendre_coefficients.size(d),
                      n_coefficients_per_direction[cell_active_fe_index]);

    ensure_existence(n_coefficients_per_direction,
                     *fe_collection,
                     q_collection,
                     cell_active_fe_index,
                     legendre_transform_matrices);

    const FullMatrix<CoefficientType> &matrix =
      legendre_transform_matrices[cell_active_fe_index];

    unrolled_coefficients.resize(Utilities::fixed_power<dim>(
      n_coefficients_per_direction[cell_active_fe_index]));
    std::fill(unrolled_coefficients.begin(),
              unrolled_coefficients.end(),
              CoefficientType(0.));

    Assert(unrolled_coefficients.size() == matrix.m(), ExcInternalError());

    Assert(local_dof_values.size() == matrix.n(),
           ExcDimensionMismatch(local_dof_values.size(), matrix.n()));

    for (unsigned int i = 0; i < unrolled_coefficients.size(); i++)
      for (unsigned int j = 0; j < local_dof_values.size(); j++)
        unrolled_coefficients[i] += matrix[i][j] * local_dof_values[j];

    legendre_coefficients.fill(unrolled_coefficients.begin());
  }
} // namespace FESeries


// explicit instantiations
#include "fe_series_legendre.inst"

DEAL_II_NAMESPACE_CLOSE
