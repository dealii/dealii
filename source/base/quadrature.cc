// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 1998 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

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


#ifndef DOXYGEN
template <>
Quadrature<0>::Quadrature()
  : is_tensor_product_flag(false)
{}

template <>
Quadrature<0>::Quadrature(const unsigned int n_q)
  : quadrature_points(n_q)
  , weights(n_q, 0)
  , is_tensor_product_flag(false)
{}
#endif



template <int dim>
Quadrature<dim>::Quadrature()
  : is_tensor_product_flag(dim == 1)
{}



template <int dim>
Quadrature<dim>::Quadrature(const unsigned int n_q)
  : quadrature_points(n_q, Point<dim>())
  , weights(n_q, 0)
  , is_tensor_product_flag(dim == 1)
{}



template <int dim>
void
Quadrature<dim>::initialize(const ArrayView<const Point<dim>> &points,
                            const ArrayView<const double>     &weights)
{
  this->weights.clear();
  if (weights.size() > 0)
    {
      AssertDimension(weights.size(), points.size());
      this->weights.insert(this->weights.end(), weights.begin(), weights.end());
    }
  else
    this->weights.resize(points.size(),
                         std::numeric_limits<double>::infinity());

  quadrature_points.clear();
  quadrature_points.insert(quadrature_points.end(),
                           points.begin(),
                           points.end());

  is_tensor_product_flag = dim == 1;
}



template <int dim>
Quadrature<dim>::Quadrature(const std::vector<Point<dim>> &points,
                            const std::vector<double>     &weights)
  : quadrature_points(points)
  , weights(weights)
  , is_tensor_product_flag(dim == 1)
{
  Assert(weights.size() == points.size(),
         ExcDimensionMismatch(weights.size(), points.size()));
}



template <int dim>
Quadrature<dim>::Quadrature(std::vector<Point<dim>> &&points,
                            std::vector<double>     &&weights)
  : quadrature_points(std::move(points))
  , weights(std::move(weights))
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
          quadrature_points[present_index][d] = q1.point(i1)[d];
        quadrature_points[present_index][dim - 1] = q2.point(i2)[0];

        weights[present_index] = q1.weight(i1) * q2.weight(i2);

        ++present_index;
      }

  if constexpr (running_in_debug_mode())
    {
      if (size() > 0)
        {
          double sum = 0;
          for (unsigned int i = 0; i < size(); ++i)
            sum += weights[i];
          // we cannot guarantee the sum of weights to be exactly one, but it
          // should be near that.
          Assert((sum > 0.999999) && (sum < 1.000001), ExcInternalError());
        }
    }

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
      quadrature_points[present_index][0] = q2.point(i2)[0];

      weights[present_index] = q2.weight(i2);

      ++present_index;
    }

  if constexpr (running_in_debug_mode())
    {
      if (size() > 0)
        {
          double sum = 0;
          for (unsigned int i = 0; i < size(); ++i)
            sum += weights[i];
          // we cannot guarantee the sum of weights to be exactly one, but it
          // should be near that.
          Assert((sum > 0.999999) && (sum < 1.000001), ExcInternalError());
        }
    }
}



template <>
Quadrature<0>::Quadrature(const Quadrature<1> &)
  : EnableObserverPointer()
  , quadrature_points(1)
  , weights(1, 1.)
  , is_tensor_product_flag(false)
{}


template <>
Quadrature<1>::Quadrature(const Quadrature<0> &)
  : EnableObserverPointer()
{
  // this function should never be called -- this should be the copy constructor
  // in 1d...
  Assert(false, ExcImpossibleInDim(1));
}
#endif // DOXYGEN



template <int dim>
Quadrature<dim>::Quadrature(const Quadrature<dim != 1 ? 1 : 0> &q)
  : EnableObserverPointer()
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
          quadrature_points[k][0] = q.point(i0)[0];
          if (dim > 1)
            quadrature_points[k][1] = q.point(i1)[0];
          if (dim > 2)
            quadrature_points[k][2] = q.point(i2)[0];
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
  : EnableObserverPointer()
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
typename std::conditional_t<dim == 1,
                            std::array<Quadrature<1>, dim>,
                            const std::array<Quadrature<1>, dim> &>
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
      this->quadrature_points[k][0] = qx.point(k1)[0];
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

  // placate compiler in the dim == 1 case
  constexpr int dim_1 = dim == 2 ? 1 : 0;

  unsigned int k = 0;
  for (unsigned int k2 = 0; k2 < qy.size(); ++k2)
    for (unsigned int k1 = 0; k1 < qx.size(); ++k1)
      {
        this->quadrature_points[k][0]     = qx.point(k1)[0];
        this->quadrature_points[k][dim_1] = qy.point(k2)[0];
        this->weights[k++]                = qx.weight(k1) * qy.weight(k2);
      }
  Assert(k == this->size(), ExcInternalError());

  this->is_tensor_product_flag = true;
  this->tensor_basis       = std::make_unique<std::array<Quadrature<1>, dim>>();
  (*this->tensor_basis)[0] = qx;
  (*this->tensor_basis)[dim_1] = qy;
}



template <int dim>
QAnisotropic<dim>::QAnisotropic(const Quadrature<1> &qx,
                                const Quadrature<1> &qy,
                                const Quadrature<1> &qz)
  : Quadrature<dim>(qx.size() * qy.size() * qz.size())
{
  Assert(dim == 3, ExcImpossibleInDim(dim));

  // placate compiler in lower dimensions
  constexpr int dim_1 = dim == 3 ? 1 : 0;
  constexpr int dim_2 = dim == 3 ? 2 : 0;

  unsigned int k = 0;
  for (unsigned int k3 = 0; k3 < qz.size(); ++k3)
    for (unsigned int k2 = 0; k2 < qy.size(); ++k2)
      for (unsigned int k1 = 0; k1 < qx.size(); ++k1)
        {
          this->quadrature_points[k][0]     = qx.point(k1)[0];
          this->quadrature_points[k][dim_1] = qy.point(k2)[0];
          this->quadrature_points[k][dim_2] = qz.point(k3)[0];
          this->weights[k++] = qx.weight(k1) * qy.weight(k2) * qz.weight(k3);
        }
  Assert(k == this->size(), ExcInternalError());

  this->is_tensor_product_flag = true;
  this->tensor_basis       = std::make_unique<std::array<Quadrature<1>, dim>>();
  (*this->tensor_basis)[0] = qx;
  (*this->tensor_basis)[dim_1] = qy;
  (*this->tensor_basis)[dim_2] = qz;
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

      std::vector<Point<1>>
      create_equidistant_interval_points(const unsigned int n_copies)
      {
        std::vector<Point<1>> support_points(n_copies + 1);

        for (unsigned int copy = 0; copy < n_copies; ++copy)
          support_points[copy][0] =
            static_cast<double>(copy) / static_cast<double>(n_copies);

        support_points[n_copies][0] = 1.0;

        return support_points;
      }
    } // namespace
  }   // namespace QIteratedImplementation
} // namespace internal



template <>
QIterated<0>::QIterated(const Quadrature<1> &, const std::vector<Point<1>> &)
  : Quadrature<0>()
{
  DEAL_II_NOT_IMPLEMENTED();
}



template <>
QIterated<0>::QIterated(const Quadrature<1> &, const unsigned int)
  : Quadrature<0>()
{
  DEAL_II_NOT_IMPLEMENTED();
}



template <>
QIterated<1>::QIterated(const Quadrature<1>         &base_quadrature,
                        const std::vector<Point<1>> &intervals)
  : Quadrature<1>(
      internal::QIteratedImplementation::uses_both_endpoints(base_quadrature) ?
        (base_quadrature.size() - 1) * (intervals.size() - 1) + 1 :
        base_quadrature.size() * (intervals.size() - 1))
{
  Assert(base_quadrature.size() > 0, ExcNotInitialized());
  Assert(intervals.size() > 1, ExcZero());

  const unsigned int n_copies = intervals.size() - 1;

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
              Point<1>(base_quadrature.point(q_point)[0] *
                         (intervals[copy + 1][0] - intervals[copy][0]) +
                       intervals[copy][0]);
            this->weights[next_point] =
              base_quadrature.weight(q_point) *
              (intervals[copy + 1][0] - intervals[copy][0]);

            ++next_point;
          }
    }
  else
    // skip doubly available points
    {
      const unsigned int left_index =
        std::distance(base_quadrature.get_points().begin(),
                      std::find_if(base_quadrature.get_points().cbegin(),
                                   base_quadrature.get_points().cend(),
                                   [](const Point<1> &p) {
                                     return p == Point<1>{0.};
                                   }));

      const unsigned int right_index =
        std::distance(base_quadrature.get_points().begin(),
                      std::find_if(base_quadrature.get_points().cbegin(),
                                   base_quadrature.get_points().cend(),
                                   [](const Point<1> &p) {
                                     return p == Point<1>{1.};
                                   }));

      const unsigned double_point_offset =
        left_index + (base_quadrature.size() - right_index);

      for (unsigned int copy = 0, next_point = 0; copy < n_copies; ++copy)
        for (unsigned int q_point = 0; q_point < base_quadrature.size();
             ++q_point)
          {
            // skip the left point of this copy since we have already entered it
            // the last time
            if ((copy > 0) && (base_quadrature.point(q_point) == Point<1>(0.0)))
              {
                Assert(this->quadrature_points[next_point - double_point_offset]
                           .distance(Point<1>(
                             base_quadrature.point(q_point)[0] *
                               (intervals[copy + 1][0] - intervals[copy][0]) +
                             intervals[copy][0])) < 1e-10 /*tolerance*/,
                       ExcInternalError());

                this->weights[next_point - double_point_offset] +=
                  base_quadrature.weight(q_point) *
                  (intervals[copy + 1][0] - intervals[copy][0]);

                continue;
              }

            this->quadrature_points[next_point] =
              Point<1>(base_quadrature.point(q_point)[0] *
                         (intervals[copy + 1][0] - intervals[copy][0]) +
                       intervals[copy][0]);

            // if this is the rightmost point of one of the non-last copies:
            // give it the double weight
            this->weights[next_point] =
              base_quadrature.weight(q_point) *
              (intervals[copy + 1][0] - intervals[copy][0]);

            ++next_point;
          }
    }

  // make sure that there is no rounding error for 0.0 and 1.0, since there
  // are multiple asserts in the library checking for equality without
  // tolerances
  for (auto &i : this->quadrature_points)
    if (std::abs(i[0] - 0.0) < 1e-12)
      i[0] = 0.0;
    else if (std::abs(i[0] - 1.0) < 1e-12)
      i[0] = 1.0;

  if constexpr (running_in_debug_mode())
    {
      double sum_of_weights = 0;
      for (unsigned int i = 0; i < this->size(); ++i)
        sum_of_weights += this->weight(i);
      Assert(std::fabs(sum_of_weights - 1) < 1e-13, ExcInternalError());
    }
}



template <>
QIterated<1>::QIterated(const Quadrature<1> &base_quadrature,
                        const unsigned int   n_copies)
  : QIterated<1>(
      base_quadrature,
      internal::QIteratedImplementation::create_equidistant_interval_points(
        n_copies))
{
  Assert(base_quadrature.size() > 0, ExcNotInitialized());
  Assert(n_copies > 0, ExcZero());
}



// construct higher dimensional quadrature formula by tensor product
// of lower dimensional iterated quadrature formulae
template <int dim>
QIterated<dim>::QIterated(const Quadrature<1>         &base_quadrature,
                          const std::vector<Point<1>> &intervals)
  : Quadrature<dim>(QIterated<dim - 1>(base_quadrature, intervals),
                    QIterated<1>(base_quadrature, intervals))
{}



template <int dim>
QIterated<dim>::QIterated(const Quadrature<1> &base_quadrature,
                          const unsigned int   n_copies)
  : Quadrature<dim>(QIterated<dim - 1>(base_quadrature, n_copies),
                    QIterated<1>(base_quadrature, n_copies))
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
