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



#include <deal.II/fe/fe_series.h>
#include <deal.II/base/numbers.h>

#include <cctype>
#include <iostream>

#include <deal.II/base/config.h>
#ifdef DEAL_II_WITH_GSL
#include <gsl/gsl_sf_legendre.h>
#endif


DEAL_II_NAMESPACE_OPEN

namespace FESeries
{

  /*-------------- Fourier -------------------------------*/

  void set_k_vectors(Table<1, Tensor<1,1> > &k_vectors,
                     const unsigned int N)
  {
    k_vectors.reinit(TableIndices<1>(N));
    for (unsigned int i=0; i<N; ++i)
      k_vectors(i)[0] = 2. * numbers::PI * i;

  }

  void set_k_vectors(Table<2, Tensor<1,2> > &k_vectors,
                     const unsigned int N)
  {
    k_vectors.reinit(TableIndices<2>(N,N));
    for (unsigned int i=0; i<N; ++i)
      for (unsigned int j=0; j<N; ++j)
        {
          k_vectors(i,j)[0] = 2. * numbers::PI * i;
          k_vectors(i,j)[1] = 2. * numbers::PI * j;
        }
  }

  void set_k_vectors(Table<3, Tensor<1,3> > &k_vectors,
                     const unsigned int N)
  {
    k_vectors.reinit(TableIndices<3>(N,N,N));
    for (unsigned int i=0; i<N; ++i)
      for (unsigned int j=0; j<N; ++j)
        for (unsigned int k=0; k<N; ++k)
          {
            k_vectors(i,j,k)[0] = 2. * numbers::PI * i;
            k_vectors(i,j,k)[0] = 2. * numbers::PI * j;
            k_vectors(i,j,k)[0] = 2. * numbers::PI * k;
          }
  }


  template <int dim>
  Fourier<dim>::Fourier(const unsigned int N,
                        const hp::FECollection<dim> &fe_collection,
                        const hp::QCollection<dim> &q_collection)
    :
    fe_collection(&fe_collection),
    q_collection(&q_collection),
    fourier_transform_matrices(fe_collection.size())
  {
    set_k_vectors(k_vectors,N);
    unrolled_coefficients.resize(k_vectors.n_elements());
  }

  template <int dim>
  void Fourier<dim>::calculate(const Vector<double>             &local_dof_values,
                               const unsigned int                cell_active_fe_index,
                               Table<dim,std::complex<double> > &fourier_coefficients)
  {
    ensure_existence(cell_active_fe_index);
    const FullMatrix<std::complex<double> > &matrix = fourier_transform_matrices[cell_active_fe_index];

    std::fill(unrolled_coefficients.begin(),
              unrolled_coefficients.end(),
              std::complex<double>(0.));

    Assert (unrolled_coefficients.size() == matrix.m(),
            ExcInternalError());

    Assert (local_dof_values.size() == matrix.n(),
            ExcDimensionMismatch(local_dof_values.size(),matrix.n()));

    for (unsigned int i = 0; i < unrolled_coefficients.size(); i++)
      for (unsigned int j = 0; j < local_dof_values.size(); j++)
        unrolled_coefficients[i] += matrix[i][j] *
                                    local_dof_values[j];

    fourier_coefficients.fill(unrolled_coefficients.begin());
  }

  template <int dim>
  std::complex<double> integrate(const FiniteElement<dim> &fe,
                                 const Quadrature<dim> &quadrature,
                                 const Tensor<1,dim> &k_vector,
                                 const unsigned int j)
  {
    std::complex<double> sum = 0;
    for (unsigned int q=0; q<quadrature.size(); ++q)
      {
        const Point<dim> &x_q = quadrature.point(q);
        sum += std::exp(std::complex<double>(0,1) *
                        (k_vector * x_q)) *
               fe.shape_value(j,x_q) *
               quadrature.weight(q);
      }
    return sum;
  }


  template <>
  void Fourier<2>::ensure_existence(const unsigned int fe)
  {
    Assert (fe < fe_collection->size(),
            ExcIndexRange(fe,0,fe_collection->size()))

    if (fourier_transform_matrices[fe].m() == 0)
      {
        fourier_transform_matrices[fe].reinit(k_vectors.n_elements(),
                                              (*fe_collection)[fe].dofs_per_cell);

        unsigned int k = 0;
        for (unsigned int k1=0; k1<k_vectors.size(0); ++k1)
          for (unsigned int k2=0; k2<k_vectors.size(1); ++k2,k++)
            for (unsigned int j=0; j<(*fe_collection)[fe].dofs_per_cell; ++j)
              fourier_transform_matrices[fe](k,j) = integrate((*fe_collection)[fe],
                                                              (*q_collection)[fe],
                                                              k_vectors(k1,k2),
                                                              j);
      }
  }

  template <>
  void Fourier<3>::ensure_existence(const unsigned int fe)
  {
    Assert (fe < fe_collection->size(),
            ExcIndexRange(fe,0,fe_collection->size()))

    if (fourier_transform_matrices[fe].m() == 0)
      {
        fourier_transform_matrices[fe].reinit(k_vectors.n_elements(),
                                              (*fe_collection)[fe].dofs_per_cell);

        unsigned int k = 0;
        for (unsigned int k1=0; k1<k_vectors.size(0); ++k1)
          for (unsigned int k2=0; k2<k_vectors.size(1); ++k2)
            for (unsigned int k3=0; k3<k_vectors.size(2); ++k3, k++)
              for (unsigned int j=0; j<(*fe_collection)[fe].dofs_per_cell; ++j)
                fourier_transform_matrices[fe](k,j) = integrate((*fe_collection)[fe],
                                                                (*q_collection)[fe],
                                                                k_vectors(k1,k2,k3),
                                                                j);
      }
  }

  template <>
  void Fourier<1>::ensure_existence(const unsigned int fe)
  {
    Assert (fe < fe_collection->size(),
            ExcIndexRange(fe,0,fe_collection->size()))

    if (fourier_transform_matrices[fe].m() == 0)
      {
        fourier_transform_matrices[fe].reinit(k_vectors.n_elements(),
                                              (*fe_collection)[fe].dofs_per_cell);

        for (unsigned int k=0; k<k_vectors.size(0); ++k)
          for (unsigned int j=0; j<(*fe_collection)[fe].dofs_per_cell; ++j)
            fourier_transform_matrices[fe](k,j) = integrate((*fe_collection)[fe],
                                                            (*q_collection)[fe],
                                                            k_vectors(k),
                                                            j);
      }
  }


  /*-------------- Legendre -------------------------------*/
  DeclException2 (ExcLegendre, int, double,
                  << "x["<<arg1 << "] = "<<arg2 << " is not in [0,1]");

  /* dim dimensional Legendre function with indices @p indices
   * evaluated at @p x_q in [0,1]^dim.
   */
  template <int dim>
  double Lh(const Point<dim>  &x_q,
            const TableIndices<dim> &indices)
  {
#ifdef DEAL_II_WITH_GSL
    double res = 1.0;
    for (unsigned int d = 0; d < dim; d++)
      {
        const double x = 2.0*(x_q[d]-0.5);
        Assert ( (x_q[d] <= 1.0) && (x_q[d] >= 0.),
                 ExcLegendre(d,x_q[d]));
        const int ind = indices[d];
        res *= sqrt(2.0) * gsl_sf_legendre_Pl (ind, x);
      }
    return res;

#else

    (void)x_q;
    (void)indices;
    AssertThrow(false, ExcMessage("deal.II has to be configured with GSL"
                                  "in order to use Legendre transformation."));
    return 0;
#endif
  }

  /*
   * Multiplier in Legendre coefficients
   */
  template <int dim>
  double multiplier(const TableIndices<dim> &indices)
  {
    double res = 1.0;
    for (unsigned int d = 0; d < dim; d++)
      res *= (0.5+indices[d]);

    return res;
  }


  template <int dim>
  Legendre<dim>::Legendre(const unsigned int size_in_each_direction,
                          const hp::FECollection<dim> &fe_collection,
                          const hp::QCollection<dim> &q_collection)
    :
    N(size_in_each_direction),
    fe_collection(&fe_collection),
    q_collection(&q_collection),
    legendre_transform_matrices(fe_collection.size()),
    unrolled_coefficients(Utilities::fixed_power<dim>(N),
                          0.)
  {
  }

  template <int dim>
  void Legendre<dim>::calculate(const dealii::Vector<double> &local_dof_values,
                                const unsigned int            cell_active_fe_index,
                                Table<dim,double>            &legendre_coefficients)
  {
    ensure_existence(cell_active_fe_index);
    const FullMatrix<double> &matrix = legendre_transform_matrices[cell_active_fe_index];

    std::fill(unrolled_coefficients.begin(),
              unrolled_coefficients.end(),
              0.);

    Assert (unrolled_coefficients.size() == matrix.m(),
            ExcInternalError());

    Assert (local_dof_values.size() == matrix.n(),
            ExcDimensionMismatch(local_dof_values.size(),matrix.n()));

    for (unsigned int i = 0; i < unrolled_coefficients.size(); i++)
      for (unsigned int j = 0; j < local_dof_values.size(); j++)
        unrolled_coefficients[i] += matrix[i][j] *
                                    local_dof_values[j];

    legendre_coefficients.fill(unrolled_coefficients.begin());
  }


  template <int dim>
  double integrate_Legendre(const FiniteElement<dim> &fe,
                            const Quadrature<dim> &quadrature,
                            const TableIndices<dim> &indices,
                            const unsigned int dof)
  {
    double sum = 0;
    for (unsigned int q=0; q<quadrature.size(); ++q)
      {
        const Point<dim> &x_q = quadrature.point(q);
        sum += Lh(x_q, indices) *
               fe.shape_value(dof,x_q) *
               quadrature.weight(q);
      }
    return sum * multiplier(indices);
  }

  template <>
  void Legendre<1>::ensure_existence(const unsigned int fe)
  {
    Assert (fe < fe_collection->size(),
            ExcIndexRange(fe,0,fe_collection->size()))
    if (legendre_transform_matrices[fe].m() == 0)
      {
        legendre_transform_matrices[fe].reinit(N,
                                               (*fe_collection)[fe].dofs_per_cell);

        for (unsigned int k=0; k<N; ++k)
          for (unsigned int j=0; j<(*fe_collection)[fe].dofs_per_cell; ++j)
            legendre_transform_matrices[fe](k,j) =
              integrate_Legendre((*fe_collection)[fe],
                                 (*q_collection)[fe],
                                 TableIndices<1>(k),
                                 j);
      }
  }



  template <>
  void Legendre<2>::ensure_existence(const unsigned int fe)
  {
    Assert (fe < fe_collection->size(),
            ExcIndexRange(fe,0,fe_collection->size()))

    if (legendre_transform_matrices[fe].m() == 0)
      {
        legendre_transform_matrices[fe].reinit(N*N,
                                               (*fe_collection)[fe].dofs_per_cell);

        unsigned int k = 0;
        for (unsigned int k1=0; k1<N; ++k1)
          for (unsigned int k2=0; k2<N; ++k2,k++)
            for (unsigned int j=0; j<(*fe_collection)[fe].dofs_per_cell; ++j)
              legendre_transform_matrices[fe](k,j) =
                integrate_Legendre((*fe_collection)[fe],
                                   (*q_collection)[fe],
                                   TableIndices<2>(k1,k2),
                                   j);
      }
  }

  template <>
  void Legendre<3>::ensure_existence(const unsigned int fe)
  {
    Assert (fe < fe_collection->size(),
            ExcIndexRange(fe,0,fe_collection->size()))

    if (legendre_transform_matrices[fe].m() == 0)
      {
        legendre_transform_matrices[fe].reinit(N*N*N,
                                               (*fe_collection)[fe].dofs_per_cell);

        unsigned int k = 0;
        for (unsigned int k1=0; k1<N; ++k1)
          for (unsigned int k2=0; k2<N; ++k2)
            for (unsigned int k3=0; k3<N; ++k3, k++)
              for (unsigned int j=0; j<(*fe_collection)[fe].dofs_per_cell; ++j)
                legendre_transform_matrices[fe](k,j) =
                  integrate_Legendre((*fe_collection)[fe],
                                     (*q_collection)[fe],
                                     TableIndices<3>(k1,k2,k3),
                                     j);
      }
  }

  /*-------------- linear_regression -------------------------------*/
  std::pair<double,double> linear_regression(const std::vector<double> &x,
                                             const std::vector<double> &y)
  {
    FullMatrix<double> K(2,2), invK(2,2);
    Vector<double>     X(2), B(2);

    Assert (x.size() == y.size(),
            ExcMessage("x and y are expected to have the same size"));

    Assert (x.size() >= 2,
            dealii::ExcMessage("at least two points are required for linear regression fit"));

    double sum_1 = 0.0,
           sum_x = 0.0,
           sum_x2= 0.0,
           sum_y = 0.0,
           sum_xy= 0.0;

    for (unsigned int i = 0; i < x.size(); i++)
      {
        sum_1  += 1.0;
        sum_x  += x[i];
        sum_x2 += x[i]*x[i];
        sum_y  += y[i];
        sum_xy += x[i]*y[i];
      }

    K(0,0) = sum_1;
    K(0,1) = sum_x;
    K(1,0) = sum_x;
    K(1,1) = sum_x2;

    B(0)   = sum_y;
    B(1)   = sum_xy;

    invK.invert(K);
    invK.vmult(X,B,false);

    return std::make_pair(X(1),X(0));
  }


} // end of namespace FESeries



/*-------------- Explicit Instantiations -------------------------------*/
template class FESeries::Fourier<1>;
template class FESeries::Fourier<2>;
template class FESeries::Fourier<3>;
template class FESeries::Legendre<1>;
template class FESeries::Legendre<2>;
template class FESeries::Legendre<3>;



DEAL_II_NAMESPACE_CLOSE
