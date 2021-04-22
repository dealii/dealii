// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2020 by the deal.II authors
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

#include <deal.II/base/memory_consumption.h>
#include <deal.II/base/quadrature.h>
#include <deal.II/base/utilities.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <memory>
#include <vector>

DEAL_II_NAMESPACE_OPEN


template <>
Quadrature<0>::Quadrature(const unsigned int n_q)
  : quadrature_points(n_q)
  , weights(n_q, 0)
  , is_tensor_product_flag(false)
{}



template <int dim>
Quadrature<dim>::Quadrature(const unsigned int n_q)
  : quadrature_points(n_q, Point<dim>())
  , weights(n_q, 0)
  , is_tensor_product_flag(dim == 1)
{}



template <int dim>
void
Quadrature<dim>::initialize(const std::vector<Point<dim>> &p,
                            const std::vector<double> &    w)
{
  AssertDimension(w.size(), p.size());
  quadrature_points = p;
  weights           = w;
}



template <int dim>
Quadrature<dim>::Quadrature(const std::vector<Point<dim>> &points,
                            const std::vector<double> &    weights)
  : quadrature_points(points)
  , weights(weights)
  , is_tensor_product_flag(dim == 1)
{
  Assert(weights.size() == points.size(),
         ExcDimensionMismatch(weights.size(), points.size()));
}



template <int dim>
Quadrature<dim>::Quadrature(const std::vector<Point<dim>> &points)
  : quadrature_points(points)
  , weights(points.size(), std::numeric_limits<double>::infinity())
  , is_tensor_product_flag(dim == 1)
{
  Assert(weights.size() == points.size(),
         ExcDimensionMismatch(weights.size(), points.size()));
}



template <int dim>
Quadrature<dim>::Quadrature(const Point<dim> &point)
  : quadrature_points(std::vector<Point<dim>>(1, point))
  , weights(std::vector<double>(1, 1.))
  , is_tensor_product_flag(true)
  , tensor_basis(new std::array<Quadrature<1>, dim>())
{
  for (unsigned int i = 0; i < dim; ++i)
    {
      const std::vector<Point<1>> quad_vec_1d(1, Point<1>(point[i]));
      (*tensor_basis)[i] = Quadrature<1>(quad_vec_1d, weights);
    }
}



#ifndef DOXYGEN
template <>
Quadrature<1>::Quadrature(const Point<1> &point)
  : quadrature_points(std::vector<Point<1>>(1, point))
  , weights(std::vector<double>(1, 1.))
  , is_tensor_product_flag(true)
{}



template <>
Quadrature<0>::Quadrature(const Point<0> &)
  : is_tensor_product_flag(false)
{
  Assert(false, ExcImpossibleInDim(0));
}



template <>
Quadrature<0>::Quadrature(const SubQuadrature &, const Quadrature<1> &)
{
  Assert(false, ExcImpossibleInDim(0));
}
#endif // DOXYGEN



template <int dim>
Quadrature<dim>::Quadrature(const SubQuadrature &q1, const Quadrature<1> &q2)
  : quadrature_points(q1.size() * q2.size())
  , weights(q1.size() * q2.size())
  , is_tensor_product_flag(q1.is_tensor_product())
{
  unsigned int present_index = 0;
  for (unsigned int i2 = 0; i2 < q2.size(); ++i2)
    for (unsigned int i1 = 0; i1 < q1.size(); ++i1)
      {
        // compose coordinates of new quadrature point by tensor product in the
        // last component
        for (unsigned int d = 0; d < dim - 1; ++d)
          quadrature_points[present_index](d) = q1.point(i1)(d);
        quadrature_points[present_index](dim - 1) = q2.point(i2)(0);

        weights[present_index] = q1.weight(i1) * q2.weight(i2);

        ++present_index;
      }

#ifdef DEBUG
  if (size() > 0)
    {
      double sum = 0;
      for (unsigned int i = 0; i < size(); ++i)
        sum += weights[i];
      // we cannot guarantee the sum of weights to be exactly one, but it should
      // be near that.
      Assert((sum > 0.999999) && (sum < 1.000001), ExcInternalError());
    }
#endif

  if (is_tensor_product_flag)
    {
      tensor_basis = std::make_unique<std::array<Quadrature<1>, dim>>();
      for (unsigned int i = 0; i < dim - 1; ++i)
        (*tensor_basis)[i] = q1.get_tensor_basis()[i];
      (*tensor_basis)[dim - 1] = q2;
    }
}


#ifndef DOXYGEN
template <>
Quadrature<1>::Quadrature(const SubQuadrature &, const Quadrature<1> &q2)
  : quadrature_points(q2.size())
  , weights(q2.size())
  , is_tensor_product_flag(true)
{
  unsigned int present_index = 0;
  for (unsigned int i2 = 0; i2 < q2.size(); ++i2)
    {
      // compose coordinates of new quadrature point by tensor product in the
      // last component
      quadrature_points[present_index](0) = q2.point(i2)(0);

      weights[present_index] = q2.weight(i2);

      ++present_index;
    }

#  ifdef DEBUG
  if (size() > 0)
    {
      double sum = 0;
      for (unsigned int i = 0; i < size(); ++i)
        sum += weights[i];
      // we cannot guarantee the sum of weights to be exactly one, but it should
      // be near that.
      Assert((sum > 0.999999) && (sum < 1.000001), ExcInternalError());
    }
#  endif
}



template <>
Quadrature<0>::Quadrature(const Quadrature<1> &)
  : Subscriptor()
  ,
  // quadrature_points(1),
  weights(1, 1.)
  , is_tensor_product_flag(false)
{}


template <>
Quadrature<1>::Quadrature(const Quadrature<0> &)
  : Subscriptor()
{
  // this function should never be called -- this should be the copy constructor
  // in 1d...
  Assert(false, ExcImpossibleInDim(1));
}
#endif // DOXYGEN



template <int dim>
Quadrature<dim>::Quadrature(const Quadrature<dim != 1 ? 1 : 0> &q)
  : Subscriptor()
  , quadrature_points(Utilities::fixed_power<dim>(q.size()))
  , weights(Utilities::fixed_power<dim>(q.size()))
  , is_tensor_product_flag(true)
{
  Assert(dim <= 3, ExcNotImplemented());

  const unsigned int n0 = q.size();
  const unsigned int n1 = (dim > 1) ? n0 : 1;
  const unsigned int n2 = (dim > 2) ? n0 : 1;

  unsigned int k = 0;
  for (unsigned int i2 = 0; i2 < n2; ++i2)
    for (unsigned int i1 = 0; i1 < n1; ++i1)
      for (unsigned int i0 = 0; i0 < n0; ++i0)
        {
          quadrature_points[k](0) = q.point(i0)(0);
          if (dim > 1)
            quadrature_points[k](1) = q.point(i1)(0);
          if (dim > 2)
            quadrature_points[k](2) = q.point(i2)(0);
          weights[k] = q.weight(i0);
          if (dim > 1)
            weights[k] *= q.weight(i1);
          if (dim > 2)
            weights[k] *= q.weight(i2);
          ++k;
        }

  tensor_basis = std::make_unique<std::array<Quadrature<1>, dim>>();
  for (unsigned int i = 0; i < dim; ++i)
    (*tensor_basis)[i] = q;
}



template <int dim>
Quadrature<dim>::Quadrature(const Quadrature<dim> &q)
  : Subscriptor()
  , quadrature_points(q.quadrature_points)
  , weights(q.weights)
  , is_tensor_product_flag(q.is_tensor_product_flag)
{
  if (dim > 1 && is_tensor_product_flag)
    tensor_basis =
      std::make_unique<std::array<Quadrature<1>, dim>>(*q.tensor_basis);
}



template <int dim>
Quadrature<dim> &
Quadrature<dim>::operator=(const Quadrature<dim> &q)
{
  weights                = q.weights;
  quadrature_points      = q.quadrature_points;
  is_tensor_product_flag = q.is_tensor_product_flag;
  if (dim > 1 && is_tensor_product_flag)
    {
      if (tensor_basis == nullptr)
        tensor_basis =
          std::make_unique<std::array<Quadrature<1>, dim>>(*q.tensor_basis);
      else
        *tensor_basis = *q.tensor_basis;
    }
  return *this;
}



template <int dim>
bool
Quadrature<dim>::operator==(const Quadrature<dim> &q) const
{
  return ((quadrature_points == q.quadrature_points) && (weights == q.weights));
}



template <int dim>
std::size_t
Quadrature<dim>::memory_consumption() const
{
  return (MemoryConsumption::memory_consumption(quadrature_points) +
          MemoryConsumption::memory_consumption(weights));
}



template <int dim>
typename std::conditional<dim == 1,
                          std::array<Quadrature<1>, dim>,
                          const std::array<Quadrature<1>, dim> &>::type
Quadrature<dim>::get_tensor_basis() const
{
  Assert(this->is_tensor_product_flag == true,
         ExcMessage("This function only makes sense if "
                    "this object represents a tensor product!"));
  Assert(tensor_basis != nullptr, ExcInternalError());

  return *tensor_basis;
}


#ifndef DOXYGEN
template <>
std::array<Quadrature<1>, 1>
Quadrature<1>::get_tensor_basis() const
{
  Assert(this->is_tensor_product_flag == true,
         ExcMessage("This function only makes sense if "
                    "this object represents a tensor product!"));

  return std::array<Quadrature<1>, 1>{{*this}};
}
#endif



//---------------------------------------------------------------------------
template <int dim>
QAnisotropic<dim>::QAnisotropic(const Quadrature<1> &qx)
  : Quadrature<dim>(qx.size())
{
  Assert(dim == 1, ExcImpossibleInDim(dim));
  unsigned int k = 0;
  for (unsigned int k1 = 0; k1 < qx.size(); ++k1)
    {
      this->quadrature_points[k](0) = qx.point(k1)(0);
      this->weights[k++]            = qx.weight(k1);
    }
  Assert(k == this->size(), ExcInternalError());
  this->is_tensor_product_flag = true;
}



template <int dim>
QAnisotropic<dim>::QAnisotropic(const Quadrature<1> &qx,
                                const Quadrature<1> &qy)
  : Quadrature<dim>(qx.size() * qy.size())
{
  Assert(dim == 2, ExcImpossibleInDim(dim));
}



template <>
QAnisotropic<2>::QAnisotropic(const Quadrature<1> &qx, const Quadrature<1> &qy)
  : Quadrature<2>(qx.size() * qy.size())
{
  unsigned int k = 0;
  for (unsigned int k2 = 0; k2 < qy.size(); ++k2)
    for (unsigned int k1 = 0; k1 < qx.size(); ++k1)
      {
        this->quadrature_points[k](0) = qx.point(k1)(0);
        this->quadrature_points[k](1) = qy.point(k2)(0);
        this->weights[k++]            = qx.weight(k1) * qy.weight(k2);
      }
  Assert(k == this->size(), ExcInternalError());
  this->is_tensor_product_flag = true;
  const std::array<Quadrature<1>, 2> q_array{{qx, qy}};
  this->tensor_basis = std::make_unique<std::array<Quadrature<1>, 2>>(q_array);
}



template <int dim>
QAnisotropic<dim>::QAnisotropic(const Quadrature<1> &qx,
                                const Quadrature<1> &qy,
                                const Quadrature<1> &qz)
  : Quadrature<dim>(qx.size() * qy.size() * qz.size())
{
  Assert(dim == 3, ExcImpossibleInDim(dim));
}



template <>
QAnisotropic<3>::QAnisotropic(const Quadrature<1> &qx,
                              const Quadrature<1> &qy,
                              const Quadrature<1> &qz)
  : Quadrature<3>(qx.size() * qy.size() * qz.size())
{
  unsigned int k = 0;
  for (unsigned int k3 = 0; k3 < qz.size(); ++k3)
    for (unsigned int k2 = 0; k2 < qy.size(); ++k2)
      for (unsigned int k1 = 0; k1 < qx.size(); ++k1)
        {
          this->quadrature_points[k](0) = qx.point(k1)(0);
          this->quadrature_points[k](1) = qy.point(k2)(0);
          this->quadrature_points[k](2) = qz.point(k3)(0);
          this->weights[k++] = qx.weight(k1) * qy.weight(k2) * qz.weight(k3);
        }
  Assert(k == this->size(), ExcInternalError());
  this->is_tensor_product_flag = true;
  const std::array<Quadrature<1>, 3> q_array{{qx, qy, qz}};
  this->tensor_basis = std::make_unique<std::array<Quadrature<1>, 3>>(q_array);
}



// ------------------------------------------------------------ //

namespace internal
{
  namespace QIteratedImplementation
  {
    namespace
    {
      bool
      uses_both_endpoints(const Quadrature<1> &base_quadrature)
      {
        const bool at_left =
          std::any_of(base_quadrature.get_points().cbegin(),
                      base_quadrature.get_points().cend(),
                      [](const Point<1> &p) { return p == Point<1>{0.}; });
        const bool at_right =
          std::any_of(base_quadrature.get_points().cbegin(),
                      base_quadrature.get_points().cend(),
                      [](const Point<1> &p) { return p == Point<1>{1.}; });
        return (at_left && at_right);
      }
    } // namespace
  }   // namespace QIteratedImplementation
} // namespace internal



template <>
QIterated<0>::QIterated(const Quadrature<1> &, const unsigned int)
  : Quadrature<0>()
{
  Assert(false, ExcNotImplemented());
}



template <>
QIterated<1>::QIterated(const Quadrature<1> &base_quadrature,
                        const unsigned int   n_copies)
  : Quadrature<1>(
      internal::QIteratedImplementation::uses_both_endpoints(base_quadrature) ?
        (base_quadrature.size() - 1) * n_copies + 1 :
        base_quadrature.size() * n_copies)
{
  Assert(base_quadrature.size() > 0, ExcNotInitialized());
  Assert(n_copies > 0, ExcZero());

  if (!internal::QIteratedImplementation::uses_both_endpoints(base_quadrature))
    // we don't have to skip some points in order to get a reasonable quadrature
    // formula
    {
      unsigned int next_point = 0;
      for (unsigned int copy = 0; copy < n_copies; ++copy)
        for (unsigned int q_point = 0; q_point < base_quadrature.size();
             ++q_point)
          {
            this->quadrature_points[next_point] =
              Point<1>(base_quadrature.point(q_point)(0) / n_copies +
                       (1.0 * copy) / n_copies);
            this->weights[next_point] =
              base_quadrature.weight(q_point) / n_copies;

            ++next_point;
          }
    }
  else
    // skip doubly available points
    {
      unsigned int next_point = 0;

      // first find out the weights of the left and the right boundary points.
      // note that these usually are but need not necessarily be the same
      double       double_point_weight = 0;
      unsigned int n_end_points        = 0;
      for (unsigned int i = 0; i < base_quadrature.size(); ++i)
        // add up the weight if this is an endpoint
        if ((base_quadrature.point(i) == Point<1>(0.0)) ||
            (base_quadrature.point(i) == Point<1>(1.0)))
          {
            double_point_weight += base_quadrature.weight(i);
            ++n_end_points;
          }
      // scale the weight correctly
      double_point_weight /= n_copies;

      // make sure the base quadrature formula has only one quadrature point per
      // end point
      Assert(n_end_points == 2, ExcInvalidQuadratureFormula());


      for (unsigned int copy = 0; copy < n_copies; ++copy)
        for (unsigned int q_point = 0; q_point < base_quadrature.size();
             ++q_point)
          {
            // skip the left point of this copy since we have already entered it
            // the last time
            if ((copy > 0) && (base_quadrature.point(q_point) == Point<1>(0.0)))
              continue;

            this->quadrature_points[next_point] =
              Point<1>(base_quadrature.point(q_point)(0) / n_copies +
                       (1.0 * copy) / n_copies);

            // if this is the rightmost point of one of the non-last copies:
            // give it the double weight
            if ((copy != n_copies - 1) &&
                (base_quadrature.point(q_point) == Point<1>(1.0)))
              this->weights[next_point] = double_point_weight;
            else
              this->weights[next_point] =
                base_quadrature.weight(q_point) / n_copies;

            ++next_point;
          }
    }

#if DEBUG
  double sum_of_weights = 0;
  for (unsigned int i = 0; i < this->size(); ++i)
    sum_of_weights += this->weight(i);
  Assert(std::fabs(sum_of_weights - 1) < 1e-13, ExcInternalError());
#endif
}



// construct higher dimensional quadrature formula by tensor product
// of lower dimensional iterated quadrature formulae
template <int dim>
QIterated<dim>::QIterated(const Quadrature<1> &base_quadrature,
                          const unsigned int   N)
  : Quadrature<dim>(QIterated<dim - 1>(base_quadrature, N),
                    QIterated<1>(base_quadrature, N))
{}



// explicit instantiations; note: we need them all for all dimensions
template class Quadrature<0>;
template class Quadrature<1>;
template class Quadrature<2>;
template class Quadrature<3>;
template class QAnisotropic<1>;
template class QAnisotropic<2>;
template class QAnisotropic<3>;
template class QIterated<1>;
template class QIterated<2>;
template class QIterated<3>;

DEAL_II_NAMESPACE_CLOSE
