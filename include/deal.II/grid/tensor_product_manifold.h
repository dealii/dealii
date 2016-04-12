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

#ifndef dealii__tensor_product_manifold_h
#define dealii__tensor_product_manifold_h

#include <deal.II/base/config.h>
#include <deal.II/base/subscriptor.h>
#include <deal.II/base/point.h>
#include <deal.II/grid/manifold.h>

DEAL_II_NAMESPACE_OPEN



/**
  * @brief Tensor product manifold of two ChartManifolds.
  *
  * This manifold will combine the ChartManifolds @p A and @p B given in the
  * constructor to form a new ChartManifold by building the tensor product
  * $A\cross B$. The first @p spacedim_A
  * dimensions in the real space and the first @p chartdim_A dimensions
  * of the chart will be given by manifold @p A, while the remaining
  * coordinates are given by @p B. The manifold is to be used by a
  * <tt>Triangulation@<dim, space_dim_A+space_dim_B@></tt>.
  *
  * An example usage would be the combination of a SphericalManifold with
  * space dimension 2 and a FlatManifold with space dimension 1 to form
  * a cylindrical manifold.
  *
  * pull_back(), push_forward(), and push_forward_gradient() are implemented
  * by splitting the input argument into inputs for @p A and @p B according
  * to the given dimensions and applying the corresponding operations before
  * concatenating the result.
  *
  * @note The dimension arguments @p dim_A and @p dim_B are not used.
  *
  * @tparam dim Dimension of cells (needs to match first template argument of
  * the Triangulation to be attached to.
  * @tparam dim_A Dimension of ChartManifold A.
  * @tparam spacedim_A Spacial dimension of ChartManifold A.
  * @tparam chartdim_A Chart dimension of ChartManifold A.
  * @tparam dim_B Dimension of ChartManifold B.
  * @tparam spacedim_B Spacial dimension of ChartManifold B.
  * @tparam chartdim_B Chart dimension of ChartManifold B.
  *
  * @author Luca Heltai, Timo Heister, 2016
  */
template <int dim,
          int dim_A, int spacedim_A, int chartdim_A,
          int dim_B, int spacedim_B, int chartdim_B>
class TensorProductManifold :
  public ChartManifold<dim,spacedim_A+spacedim_B,chartdim_A+chartdim_B>
{
public:
  /**
   * The chart dimension is the sum of the chart dimensions of the manifolds
   * @p A and @p B.
   */
  static const unsigned int chartdim = chartdim_A+chartdim_B;
  /**
   * The space dimension is the sum of the space dimensions of the manifolds
   * @p A and @p B.
   */
  static const unsigned int spacedim = spacedim_A+spacedim_B;

  /**
   * Constructor.
   */
  TensorProductManifold (
    const ChartManifold<dim_A, spacedim_A, chartdim_A> &manifold_A,
    const ChartManifold<dim_B, spacedim_B, chartdim_B> &manifold_B);

  /**
   * Pull back operation.
   */
  virtual
  Point<chartdim>
  pull_back(const Point<spacedim> &space_point) const;

  /**
   * Push forward operation.
   */
  virtual
  Point<spacedim>
  push_forward(const Point<chartdim> &chart_point) const;

  /**
   * Gradient.
   */
  virtual
  DerivativeForm<1,chartdim,spacedim>
  push_forward_gradient(const Point<chartdim> &chart_point) const;

private:
  SmartPointer<const ChartManifold<dim_A, spacedim_A, chartdim_A>,
               TensorProductManifold<dim, dim_A, spacedim_A, chartdim_A, dim_B, spacedim_B, chartdim_B> > manifold_A;

  SmartPointer<const ChartManifold<dim_B, spacedim_B, chartdim_B>,
               TensorProductManifold<dim, dim_A, spacedim_A, chartdim_A, dim_B, spacedim_B, chartdim_B> > manifold_B;
};



/*------------------Template Implementations------------------------*/



namespace internal
{
  namespace TensorProductManifold
  {
    template <int dim1, int dim2>
    Tensor<1,dim1+dim2> concat(const Tensor<1,dim1> &p1, const Tensor<1,dim2> &p2)
    {
      Tensor<1,dim1+dim2> r;
      for (unsigned int d=0; d<dim1; ++d)
        r[d] = p1[d];
      for (unsigned int d=0; d<dim2; ++d)
        r[dim1+d] = p2[d];
      return r;
    }

    template <int dim1, int dim2>
    Point<dim1+dim2> concat(const Point<dim1> &p1, const Point<dim2> &p2)
    {
      Point<dim1+dim2> r;
      for (unsigned int d=0; d<dim1; ++d)
        r[d] = p1[d];
      for (unsigned int d=0; d<dim2; ++d)
        r[dim1+d] = p2[d];
      return r;
    }

    template <int dim1, int dim2>
    void split_point(const Point<dim1+dim2> &source, Point<dim1> &p1, Point<dim2> &p2)
    {
      for (unsigned int d=0; d<dim1; ++d)
        p1[d] = source[d];
      for (unsigned int d=0; d<dim2; ++d)
        p2[d] = source[dim1+d];
    }

  }
}

template <int dim,
          int dim_A, int spacedim_A, int chartdim_A,
          int dim_B, int spacedim_B, int chartdim_B>
TensorProductManifold<dim, dim_A, spacedim_A, chartdim_A, dim_B, spacedim_B, chartdim_B>
::TensorProductManifold(
  const ChartManifold<dim_A, spacedim_A, chartdim_A> &manifold_A,
  const ChartManifold<dim_B, spacedim_B, chartdim_B> &manifold_B)
  : ChartManifold<dim,spacedim_A+spacedim_B,chartdim_A+chartdim_B> (
    internal::TensorProductManifold::concat(
      manifold_A.get_periodicity(),
      manifold_B.get_periodicity())),
  manifold_A (&manifold_A),
  manifold_B (&manifold_B)
{}

template <int dim,
          int dim_A, int spacedim_A, int chartdim_A,
          int dim_B, int spacedim_B, int chartdim_B>
Point<TensorProductManifold<dim, dim_A, spacedim_A, chartdim_A, dim_B, spacedim_B, chartdim_B>::chartdim>
TensorProductManifold<dim, dim_A, spacedim_A, chartdim_A, dim_B, spacedim_B, chartdim_B>
::pull_back(const Point<TensorProductManifold<dim, dim_A, spacedim_A, chartdim_A, dim_B, spacedim_B, chartdim_B>::spacedim> &space_point) const
{
  Point<spacedim_A> space_point_A;
  Point<spacedim_B> space_point_B;
  internal::TensorProductManifold::split_point(space_point, space_point_A, space_point_B);

  Point<chartdim_A> result_A = manifold_A->pull_back(space_point_A);
  Point<chartdim_B> result_B = manifold_B->pull_back(space_point_B);

  return internal::TensorProductManifold::concat(result_A, result_B);
}

template <int dim,
          int dim_A, int spacedim_A, int chartdim_A,
          int dim_B, int spacedim_B, int chartdim_B>
Point<TensorProductManifold<dim, dim_A, spacedim_A, chartdim_A, dim_B, spacedim_B, chartdim_B>::spacedim>
TensorProductManifold<dim, dim_A, spacedim_A, chartdim_A, dim_B, spacedim_B, chartdim_B>
::push_forward(const Point<TensorProductManifold<dim, dim_A, spacedim_A, chartdim_A, dim_B, spacedim_B, chartdim_B>::chartdim> &chart_point) const
{
  Point<chartdim_A> chart_point_A;
  Point<chartdim_B> chart_point_B;
  internal::TensorProductManifold::split_point(chart_point, chart_point_A, chart_point_B);

  Point<spacedim_A> result_A = manifold_A->push_forward(chart_point_A);
  Point<spacedim_B> result_B = manifold_B->push_forward(chart_point_B);

  return internal::TensorProductManifold::concat(result_A, result_B);
}

template <int dim,
          int dim_A, int spacedim_A, int chartdim_A,
          int dim_B, int spacedim_B, int chartdim_B>
DerivativeForm<1,
               TensorProductManifold<dim, dim_A, spacedim_A, chartdim_A, dim_B, spacedim_B, chartdim_B>::chartdim,
               TensorProductManifold<dim, dim_A, spacedim_A, chartdim_A, dim_B, spacedim_B, chartdim_B>::spacedim>

               TensorProductManifold<dim, dim_A, spacedim_A, chartdim_A, dim_B, spacedim_B, chartdim_B>
               ::push_forward_gradient(const Point<TensorProductManifold<dim, dim_A, spacedim_A, chartdim_A, dim_B, spacedim_B, chartdim_B>::chartdim> &chart_point) const
{
  Point<chartdim_A> chart_point_A;
  Point<chartdim_B> chart_point_B;
  internal::TensorProductManifold::split_point(chart_point, chart_point_A, chart_point_B);

  DerivativeForm<1,chartdim_A,spacedim_A> result_A
    = manifold_A->push_forward_gradient(chart_point_A);
  DerivativeForm<1,chartdim_B,spacedim_B> result_B
    = manifold_B->push_forward_gradient(chart_point_B);


  DerivativeForm<1,chartdim,spacedim> result;
  for (unsigned int i = 0; i<chartdim_A; ++i)
    for (unsigned int j = 0; j<spacedim_A; ++j)
      result[j][i] = result_A[j][i];
  for (unsigned int i = 0; i<chartdim_B; ++i)
    for (unsigned int j = 0; j<spacedim_B; ++j)
      result[j+spacedim_A][i+chartdim_A] = result_B[j][i];

  return result;
}




DEAL_II_NAMESPACE_CLOSE

#endif
