// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2001 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#include <deal.II/base/function_lib.h>
#include <deal.II/base/point.h>
#include <deal.II/base/tensor.h>

#include <deal.II/lac/vector.h>

#include <cmath>

DEAL_II_NAMESPACE_OPEN


namespace Functions
{
  template <int dim>
  CutOffFunctionBase<dim>::CutOffFunctionBase(
    const double       r,
    const Point<dim>   p,
    const unsigned int n_components,
    const unsigned int select,
    const bool         integrate_to_one,
    const double       unitary_integral_value)
    : Function<dim>(n_components)
    , center(p)
    , radius(r)
    , selected(select)
    , integrate_to_one(integrate_to_one)
    , unitary_integral_value(unitary_integral_value)
    , rescaling(integrate_to_one ? 1. / (unitary_integral_value *
                                         Utilities::fixed_power<dim>(radius)) :
                                   1.0)
  {
    Assert(r > 0, ExcMessage("You must specify a radius > 0."));
  }



  template <int dim>
  void
  CutOffFunctionBase<dim>::set_center(const Point<dim> &p)
  {
    center = p;
  }



  template <int dim>
  const Point<dim> &
  CutOffFunctionBase<dim>::get_center() const
  {
    return center;
  }



  template <int dim>
  void
  CutOffFunctionBase<dim>::set_radius(const double r)
  {
    radius = r;
    Assert(r > 0, ExcMessage("You must specify a radius > 0."));
    if (integrate_to_one)
      rescaling =
        1. / (unitary_integral_value * Utilities::fixed_power<dim>(radius));
    else
      rescaling = 1.0;
  }



  template <int dim>
  double
  CutOffFunctionBase<dim>::get_radius() const
  {
    return radius;
  }



  template <int dim>
  bool
  CutOffFunctionBase<dim>::integrates_to_one() const
  {
    return integrate_to_one;
  }



  template <int dim>
  CutOffFunctionTensorProduct<dim>::CutOffFunctionTensorProduct(
    double             radius,
    const Point<dim>  &center,
    const unsigned int n_components,
    const unsigned int select,
    const bool         integrate_to_one)
    : CutOffFunctionBase<dim>(radius,
                              center,
                              n_components,
                              select,
                              integrate_to_one)
    , initialized(false)
  {}



  template <int dim>
  void
  CutOffFunctionTensorProduct<dim>::set_center(const Point<dim> &p)
  {
    Assert(initialized, ExcNotInitialized());
    for (unsigned int i = 0; i < dim; ++i)
      base[i]->set_center(Point<1>(p[i]));
    CutOffFunctionBase<dim>::set_center(p);
  }



  template <int dim>
  void
  CutOffFunctionTensorProduct<dim>::set_radius(const double r)
  {
    Assert(initialized, ExcNotInitialized());
    for (unsigned int i = 0; i < dim; ++i)
      base[i]->set_radius(r);
    CutOffFunctionBase<dim>::set_radius(r);
  }



  template <int dim>
  double
  CutOffFunctionTensorProduct<dim>::value(const Point<dim>  &p,
                                          const unsigned int component) const
  {
    Assert(initialized, ExcNotInitialized());
    double ret = 1.0;
    for (unsigned int i = 0; i < dim; ++i)
      ret *= base[i]->value(Point<1>(p[i]), component);
    return ret;
  }



  template <int dim>
  Tensor<1, dim>
  CutOffFunctionTensorProduct<dim>::gradient(const Point<dim>  &p,
                                             const unsigned int component) const
  {
    Assert(initialized, ExcNotInitialized());
    Tensor<1, dim> ret;
    for (unsigned int d = 0; d < dim; ++d)
      {
        ret[d] = base[d]->gradient(Point<1>(p[d]), component)[0];
        for (unsigned int i = 0; i < dim; ++i)
          if (i != d)
            ret[d] *= base[i]->value(Point<1>(p[i]), component);
      }
    return ret;
  }



  //--------------------------------------------------------------------
  namespace
  {
    // Integral of CutOffFunctionLinfty in dimension 1, 2, and 3 when the radius
    // is one
    const double integral_Linfty[] = {2.0,
                                      3.14159265358979323846264338328,
                                      4.18879020478639098461685784437};

    // Integral of CutOffFunctionW1 in dimension 1, 2, and 3 when the radius
    // is one
    const double integral_W1[] = {1.0,
                                  1.04719755119659774615421446109,
                                  1.04719755119659774615421446109};

    // Integral of CutOffFunctionCinfty in dimension 1, 2, and 3 when the radius
    // is one
    const double integral_Cinfty[] = {1.20690032243787617533623799633,
                                      1.26811216112759608094632335664,
                                      1.1990039070192139033798473858};

    // Integral of CutOffFunctionC1 in dimension 1, 2, and 3 when the radius
    // is one
    const double integral_C1[] = {1.0,
                                  0.93417655442731527615578663815,
                                  0.821155557658032806157358815206};
  } // namespace


  template <int dim>
  CutOffFunctionLinfty<dim>::CutOffFunctionLinfty(
    const double       r,
    const Point<dim>   p,
    const unsigned int n_components,
    const unsigned int select,
    const bool         integrate_to_one)
    : CutOffFunctionBase<dim>(r,
                              p,
                              n_components,
                              select,
                              integrate_to_one,
                              integral_Linfty[dim - 1])
  {}


  template <int dim>
  double
  CutOffFunctionLinfty<dim>::value(const Point<dim>  &p,
                                   const unsigned int component) const
  {
    if (this->selected == CutOffFunctionBase<dim>::no_component ||
        component == this->selected)
      return ((this->center.distance(p) < this->radius) ? this->rescaling : 0.);
    return 0.;
  }


  template <int dim>
  void
  CutOffFunctionLinfty<dim>::value_list(const std::vector<Point<dim>> &points,
                                        std::vector<double>           &values,
                                        const unsigned int component) const
  {
    Assert(values.size() == points.size(),
           ExcDimensionMismatch(values.size(), points.size()));
    AssertIndexRange(component, this->n_components);


    if (this->selected == CutOffFunctionBase<dim>::no_component ||
        component == this->selected)
      for (unsigned int k = 0; k < values.size(); ++k)
        values[k] = (this->center.distance(points[k]) < this->radius) ?
                      this->rescaling :
                      0.;
    else
      std::fill(values.begin(), values.end(), 0.);
  }


  template <int dim>
  void
  CutOffFunctionLinfty<dim>::vector_value_list(
    const std::vector<Point<dim>> &points,
    std::vector<Vector<double>>   &values) const
  {
    Assert(values.size() == points.size(),
           ExcDimensionMismatch(values.size(), points.size()));

    for (unsigned int k = 0; k < values.size(); ++k)
      {
        const double val = (this->center.distance(points[k]) < this->radius) ?
                             this->rescaling :
                             0.;
        if (this->selected == CutOffFunctionBase<dim>::no_component)
          values[k] = val;
        else
          {
            values[k]                 = 0;
            values[k](this->selected) = val;
          }
      }
  }

  template <int dim>
  CutOffFunctionW1<dim>::CutOffFunctionW1(const double       r,
                                          const Point<dim>   p,
                                          const unsigned int n_components,
                                          const unsigned int select,
                                          const bool         integrate_to_one)
    : CutOffFunctionBase<dim>(r,
                              p,
                              n_components,
                              select,
                              integrate_to_one,
                              integral_W1[dim - 1])
  {}


  template <int dim>
  double
  CutOffFunctionW1<dim>::value(const Point<dim>  &p,
                               const unsigned int component) const
  {
    if (this->selected == CutOffFunctionBase<dim>::no_component ||
        component == this->selected)
      {
        const double d = this->center.distance(p);
        return ((d < this->radius) ?
                  (this->radius - d) / this->radius * this->rescaling :
                  0.);
      }
    return 0.;
  }


  template <int dim>
  void
  CutOffFunctionW1<dim>::value_list(const std::vector<Point<dim>> &points,
                                    std::vector<double>           &values,
                                    const unsigned int component) const
  {
    Assert(values.size() == points.size(),
           ExcDimensionMismatch(values.size(), points.size()));

    if (this->selected == CutOffFunctionBase<dim>::no_component ||
        component == this->selected)
      for (unsigned int i = 0; i < values.size(); ++i)
        {
          const double d = this->center.distance(points[i]);
          values[i]      = ((d < this->radius) ?
                              (this->radius - d) / this->radius * this->rescaling :
                              0.);
        }
    else
      std::fill(values.begin(), values.end(), 0.);
  }



  template <int dim>
  void
  CutOffFunctionW1<dim>::vector_value_list(
    const std::vector<Point<dim>> &points,
    std::vector<Vector<double>>   &values) const
  {
    Assert(values.size() == points.size(),
           ExcDimensionMismatch(values.size(), points.size()));

    for (unsigned int k = 0; k < values.size(); ++k)
      {
        const double d = this->center.distance(points[k]);
        const double val =
          (d < this->radius) ?
            (this->radius - d) / this->radius * this->rescaling :
            0.;
        if (this->selected == CutOffFunctionBase<dim>::no_component)
          values[k] = val;
        else
          {
            values[k]                 = 0;
            values[k](this->selected) = val;
          }
      }
  }


  template <int dim>
  CutOffFunctionCinfty<dim>::CutOffFunctionCinfty(
    const double       r,
    const Point<dim>   p,
    const unsigned int n_components,
    const unsigned int select,
    bool               integrate_to_one)
    : CutOffFunctionBase<dim>(r,
                              p,
                              n_components,
                              select,
                              integrate_to_one,
                              integral_Cinfty[dim - 1])
  {}


  template <int dim>
  double
  CutOffFunctionCinfty<dim>::value(const Point<dim>  &p,
                                   const unsigned int component) const
  {
    if (this->selected == CutOffFunctionBase<dim>::no_component ||
        component == this->selected)
      {
        const double d = this->center.distance(p);
        const double r = this->radius;
        if (d >= r)
          return 0.;
        const double e = -r * r / (r * r - d * d);
        return ((e < -50) ? 0. : numbers::E * std::exp(e) * this->rescaling);
      }
    return 0.;
  }


  template <int dim>
  void
  CutOffFunctionCinfty<dim>::value_list(const std::vector<Point<dim>> &points,
                                        std::vector<double>           &values,
                                        const unsigned int component) const
  {
    Assert(values.size() == points.size(),
           ExcDimensionMismatch(values.size(), points.size()));

    const double r = this->radius;

    if (this->selected == CutOffFunctionBase<dim>::no_component ||
        component == this->selected)
      for (unsigned int i = 0; i < values.size(); ++i)
        {
          const double d = this->center.distance(points[i]);
          if (d >= r)
            {
              values[i] = 0.;
            }
          else
            {
              const double e = -r * r / (r * r - d * d);
              values[i] =
                (e < -50) ? 0. : numbers::E * std::exp(e) * this->rescaling;
            }
        }
    else
      std::fill(values.begin(), values.end(), 0.);
  }


  template <int dim>
  void
  CutOffFunctionCinfty<dim>::vector_value_list(
    const std::vector<Point<dim>> &points,
    std::vector<Vector<double>>   &values) const
  {
    Assert(values.size() == points.size(),
           ExcDimensionMismatch(values.size(), points.size()));

    for (unsigned int k = 0; k < values.size(); ++k)
      {
        const double d   = this->center.distance(points[k]);
        const double r   = this->radius;
        double       val = 0.;
        if (d < this->radius)
          {
            const double e = -r * r / (r * r - d * d);
            if (e > -50)
              val = numbers::E * std::exp(e) * this->rescaling;
          }

        if (this->selected == CutOffFunctionBase<dim>::no_component)
          values[k] = val;
        else
          {
            values[k]                 = 0;
            values[k](this->selected) = val;
          }
      }
  }



  template <int dim>
  Tensor<1, dim>
  CutOffFunctionCinfty<dim>::gradient(const Point<dim> &p,
                                      const unsigned int) const
  {
    const double d = this->center.distance(p);
    const double r = this->radius;
    if (d >= r)
      return Tensor<1, dim>();
    const double e = -d * d / (r - d) / (r + d);
    return ((e < -50) ?
              Point<dim>() :
              (p - this->center) / d *
                (-2.0 * r * r / Utilities::fixed_power<2>(-r * r + d * d) * d *
                 std::exp(e)) *
                this->rescaling);
  }



  template <int dim>
  CutOffFunctionC1<dim>::CutOffFunctionC1(const double       r,
                                          const Point<dim>   p,
                                          const unsigned int n_components,
                                          const unsigned int select,
                                          bool               integrate_to_one)
    : CutOffFunctionBase<dim>(r,
                              p,
                              n_components,
                              select,
                              integrate_to_one,
                              integral_C1[dim - 1])
  {}


  template <int dim>
  double
  CutOffFunctionC1<dim>::value(const Point<dim>  &p,
                               const unsigned int component) const
  {
    if (this->selected == CutOffFunctionBase<dim>::no_component ||
        component == this->selected)
      {
        const double d = this->center.distance(p);
        const double r = this->radius;
        if (d >= r)
          return 0.;
        return .5 * (std::cos(numbers::PI * d / r) + 1) * this->rescaling;
      }
    return 0.;
  }


  template <int dim>
  void
  CutOffFunctionC1<dim>::value_list(const std::vector<Point<dim>> &points,
                                    std::vector<double>           &values,
                                    const unsigned int component) const
  {
    Assert(values.size() == points.size(),
           ExcDimensionMismatch(values.size(), points.size()));

    const double r = this->radius;

    if (this->selected == CutOffFunctionBase<dim>::no_component ||
        component == this->selected)
      for (unsigned int i = 0; i < values.size(); ++i)
        {
          const double d = this->center.distance(points[i]);
          if (d >= r)
            {
              values[i] = 0.;
            }
          else
            {
              values[i] =
                .5 * (std::cos(numbers::PI * d / r) + 1) * this->rescaling;
            }
        }
    else
      std::fill(values.begin(), values.end(), 0.);
  }


  template <int dim>
  void
  CutOffFunctionC1<dim>::vector_value_list(
    const std::vector<Point<dim>> &points,
    std::vector<Vector<double>>   &values) const
  {
    Assert(values.size() == points.size(),
           ExcDimensionMismatch(values.size(), points.size()));

    for (unsigned int k = 0; k < values.size(); ++k)
      {
        const double d   = this->center.distance(points[k]);
        const double r   = this->radius;
        double       val = 0.;
        if (d < this->radius)
          {
            val = .5 * (std::cos(numbers::PI * d / r) + 1) * this->rescaling;
          }

        if (this->selected == CutOffFunctionBase<dim>::no_component)
          values[k] = val;
        else
          {
            values[k]                 = 0;
            values[k](this->selected) = val;
          }
      }
  }



  template <int dim>
  Tensor<1, dim>
  CutOffFunctionC1<dim>::gradient(const Point<dim> &p, const unsigned int) const
  {
    const double d = this->center.distance(p);
    const double r = this->radius;
    if (d >= r)
      return Tensor<1, dim>();
    return (-0.5 * numbers::PI * std::sin(numbers::PI * d / r) / r) *
           (p - this->center) / d * this->rescaling;
  }


  // explicit instantiations
  template class CutOffFunctionBase<1>;
  template class CutOffFunctionBase<2>;
  template class CutOffFunctionBase<3>;

  template class CutOffFunctionLinfty<1>;
  template class CutOffFunctionLinfty<2>;
  template class CutOffFunctionLinfty<3>;

  template class CutOffFunctionW1<1>;
  template class CutOffFunctionW1<2>;
  template class CutOffFunctionW1<3>;

  template class CutOffFunctionCinfty<1>;
  template class CutOffFunctionCinfty<2>;
  template class CutOffFunctionCinfty<3>;

  template class CutOffFunctionC1<1>;
  template class CutOffFunctionC1<2>;
  template class CutOffFunctionC1<3>;

  template class CutOffFunctionTensorProduct<1>;
  template class CutOffFunctionTensorProduct<2>;
  template class CutOffFunctionTensorProduct<3>;
} // namespace Functions

DEAL_II_NAMESPACE_CLOSE
