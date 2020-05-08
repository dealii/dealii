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



#include <deal.II/base/numbers.h>
#include <deal.II/base/thread_management.h>

#include <deal.II/fe/fe_series.h>

#include <iostream>


DEAL_II_NAMESPACE_OPEN

namespace
{
  void set_k_vectors(Table<1, Tensor<1, 1>> &k_vectors, const unsigned int N)
  {
    k_vectors.reinit(TableIndices<1>(N));
    for (unsigned int i = 0; i < N; ++i)
      k_vectors(i)[0] = 2. * numbers::PI * i;
  }

  void set_k_vectors(Table<2, Tensor<1, 2>> &k_vectors, const unsigned int N)
  {
    k_vectors.reinit(TableIndices<2>(N, N));
    for (unsigned int i = 0; i < N; ++i)
      for (unsigned int j = 0; j < N; ++j)
        {
          k_vectors(i, j)[0] = 2. * numbers::PI * i;
          k_vectors(i, j)[1] = 2. * numbers::PI * j;
        }
  }

  void set_k_vectors(Table<3, Tensor<1, 3>> &k_vectors, const unsigned int N)
  {
    k_vectors.reinit(TableIndices<3>(N, N, N));
    for (unsigned int i = 0; i < N; ++i)
      for (unsigned int j = 0; j < N; ++j)
        for (unsigned int k = 0; k < N; ++k)
          {
            k_vectors(i, j, k)[0] = 2. * numbers::PI * i;
            k_vectors(i, j, k)[1] = 2. * numbers::PI * j;
            k_vectors(i, j, k)[2] = 2. * numbers::PI * k;
          }
  }



  template <int dim, int spacedim>
  std::complex<double>
  integrate(const FiniteElement<dim, spacedim> &fe,
            const Quadrature<dim> &             quadrature,
            const Tensor<1, dim> &              k_vector,
            const unsigned int                  j)
  {
    std::complex<double> sum = 0;
    for (unsigned int q = 0; q < quadrature.size(); ++q)
      {
        const Point<dim> &x_q = quadrature.point(q);
        sum += std::exp(std::complex<double>(0, 1) * (k_vector * x_q)) *
               fe.shape_value(j, x_q) * quadrature.weight(q);
      }
    return sum;
  }



  /*
   * Ensure that the transformation matrix for FiniteElement index
   * @p fe_index is calculated. If not, calculate it.
   */
  template <int spacedim>
  void
  ensure_existence(
    const hp::FECollection<1, spacedim> &          fe_collection,
    const hp::QCollection<1> &                     q_collection,
    const Table<1, Tensor<1, 1>> &                 k_vectors,
    const unsigned int                             fe,
    std::vector<FullMatrix<std::complex<double>>> &fourier_transform_matrices)
  {
    AssertIndexRange(fe, fe_collection.size());

    if (fourier_transform_matrices[fe].m() == 0)
      {
        fourier_transform_matrices[fe].reinit(k_vectors.n_elements(),
                                              fe_collection[fe].dofs_per_cell);

        for (unsigned int k = 0; k < k_vectors.size(0); ++k)
          for (unsigned int j = 0; j < fe_collection[fe].dofs_per_cell; ++j)
            fourier_transform_matrices[fe](k, j) =
              integrate(fe_collection[fe], q_collection[fe], k_vectors(k), j);
      }
  }

  template <int spacedim>
  void
  ensure_existence(
    const hp::FECollection<2, spacedim> &          fe_collection,
    const hp::QCollection<2> &                     q_collection,
    const Table<2, Tensor<1, 2>> &                 k_vectors,
    const unsigned int                             fe,
    std::vector<FullMatrix<std::complex<double>>> &fourier_transform_matrices)
  {
    AssertIndexRange(fe, fe_collection.size());

    if (fourier_transform_matrices[fe].m() == 0)
      {
        fourier_transform_matrices[fe].reinit(k_vectors.n_elements(),
                                              fe_collection[fe].dofs_per_cell);

        unsigned int k = 0;
        for (unsigned int k1 = 0; k1 < k_vectors.size(0); ++k1)
          for (unsigned int k2 = 0; k2 < k_vectors.size(1); ++k2, k++)
            for (unsigned int j = 0; j < fe_collection[fe].dofs_per_cell; ++j)
              fourier_transform_matrices[fe](k, j) = integrate(
                fe_collection[fe], q_collection[fe], k_vectors(k1, k2), j);
      }
  }

  template <int spacedim>
  void
  ensure_existence(
    const hp::FECollection<3, spacedim> &          fe_collection,
    const hp::QCollection<3> &                     q_collection,
    const Table<3, Tensor<1, 3>> &                 k_vectors,
    const unsigned int                             fe,
    std::vector<FullMatrix<std::complex<double>>> &fourier_transform_matrices)
  {
    AssertIndexRange(fe, fe_collection.size());

    if (fourier_transform_matrices[fe].m() == 0)
      {
        fourier_transform_matrices[fe].reinit(k_vectors.n_elements(),
                                              fe_collection[fe].dofs_per_cell);

        unsigned int k = 0;
        for (unsigned int k1 = 0; k1 < k_vectors.size(0); ++k1)
          for (unsigned int k2 = 0; k2 < k_vectors.size(1); ++k2)
            for (unsigned int k3 = 0; k3 < k_vectors.size(2); ++k3, k++)
              for (unsigned int j = 0; j < fe_collection[fe].dofs_per_cell; ++j)
                fourier_transform_matrices[fe](k, j) =
                  integrate(fe_collection[fe],
                            q_collection[fe],
                            k_vectors(k1, k2, k3),
                            j);
      }
  }
} // namespace



namespace FESeries
{
  template <int dim, int spacedim>
  Fourier<dim, spacedim>::Fourier(
    const unsigned int                     N,
    const hp::FECollection<dim, spacedim> &fe_collection,
    const hp::QCollection<dim> &           q_collection)
    : fe_collection(&fe_collection)
    , q_collection(&q_collection)
    , fourier_transform_matrices(fe_collection.size())
  {
    set_k_vectors(k_vectors, N);
    unrolled_coefficients.resize(k_vectors.n_elements());
  }



  template <int dim, int spacedim>
  inline bool
  Fourier<dim, spacedim>::
  operator==(const Fourier<dim, spacedim> &fourier) const
  {
    return ((*fe_collection == *(fourier.fe_collection)) &&
            (*q_collection == *(fourier.q_collection)) &&
            (k_vectors == fourier.k_vectors) &&
            (fourier_transform_matrices == fourier.fourier_transform_matrices));
  }



  template <int dim, int spacedim>
  void
  Fourier<dim, spacedim>::precalculate_all_transformation_matrices()
  {
    Threads::TaskGroup<> task_group;
    for (unsigned int fe = 0; fe < fe_collection->size(); ++fe)
      task_group += Threads::new_task([&, fe]() {
        ensure_existence(*fe_collection,
                         *q_collection,
                         k_vectors,
                         fe,
                         fourier_transform_matrices);
      });

    task_group.join_all();
  }



  template <int dim, int spacedim>
  template <typename Number>
  void
  Fourier<dim, spacedim>::calculate(
    const Vector<Number> &       local_dof_values,
    const unsigned int           cell_active_fe_index,
    Table<dim, CoefficientType> &fourier_coefficients)
  {
    ensure_existence(*fe_collection,
                     *q_collection,
                     k_vectors,
                     cell_active_fe_index,
                     fourier_transform_matrices);

    const FullMatrix<CoefficientType> &matrix =
      fourier_transform_matrices[cell_active_fe_index];

    std::fill(unrolled_coefficients.begin(),
              unrolled_coefficients.end(),
              CoefficientType(0.));

    Assert(unrolled_coefficients.size() == matrix.m(), ExcInternalError());

    Assert(local_dof_values.size() == matrix.n(),
           ExcDimensionMismatch(local_dof_values.size(), matrix.n()));

    for (unsigned int i = 0; i < unrolled_coefficients.size(); i++)
      for (unsigned int j = 0; j < local_dof_values.size(); j++)
        unrolled_coefficients[i] += matrix[i][j] * local_dof_values[j];

    fourier_coefficients.fill(unrolled_coefficients.begin());
  }
} // namespace FESeries


// explicit instantiations
#include "fe_series_fourier.inst"

DEAL_II_NAMESPACE_CLOSE
