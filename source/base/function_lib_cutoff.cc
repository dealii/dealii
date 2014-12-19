// ---------------------------------------------------------------------
//
// Copyright (C) 2001 - 2013 by the deal.II authors
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

#include <deal.II/base/tensor.h>
#include <deal.II/base/point.h>
#include <deal.II/base/function_lib.h>
#include <deal.II/lac/vector.h>

#include <cmath>

DEAL_II_NAMESPACE_OPEN


namespace Functions
{
  template<int dim>
  CutOffFunctionBase<dim>::CutOffFunctionBase (const double r,
                                               const Point<dim> p,
                                               const unsigned int n_components,
                                               const unsigned int select)
    :
    Function<dim> (n_components),
    center(p),
    radius(r),
    selected(select)
  {}


  template<int dim>
  void
  CutOffFunctionBase<dim>::new_center (const Point<dim> &p)
  {
    center = p;
  }


  template<int dim>
  void
  CutOffFunctionBase<dim>::new_radius (const double r)
  {
    radius = r;
  }

//////////////////////////////////////////////////////////////////////

  template<int dim>
  CutOffFunctionLinfty<dim>::CutOffFunctionLinfty (const double r,
                                                   const Point<dim> p,
                                                   const unsigned int n_components,
                                                   const unsigned int select)
    :
    CutOffFunctionBase<dim> (r, p, n_components, select)
  {}


  template<int dim>
  double
  CutOffFunctionLinfty<dim>::value (const Point<dim>   &p,
                                    const unsigned int component) const
  {
    if (this->selected==CutOffFunctionBase<dim>::no_component
        ||
        component==this->selected)
      return ((this->center.distance(p)<this->radius) ? 1. : 0.);
    return 0.;
  }


  template<int dim>
  void
  CutOffFunctionLinfty<dim>::value_list (const std::vector<Point<dim> > &points,
                                         std::vector<double>            &values,
                                         const unsigned int component) const
  {
    Assert (values.size() == points.size(),
            ExcDimensionMismatch(values.size(), points.size()));
    Assert (component < this->n_components,
            ExcIndexRange(component,0,this->n_components));


    if (this->selected==CutOffFunctionBase<dim>::no_component
        ||
        component==this->selected)
      for (unsigned int k=0; k<values.size(); ++k)
        values[k] = (this->center.distance(points[k])<this->radius) ? 1. : 0.;
    else
      std::fill (values.begin(), values.end(), 0.);
  }


  template<int dim>
  void
  CutOffFunctionLinfty<dim>::vector_value_list (
    const std::vector<Point<dim> > &points,
    std::vector<Vector<double> >           &values) const
  {
    Assert (values.size() == points.size(),
            ExcDimensionMismatch(values.size(), points.size()));

    for (unsigned int k=0; k<values.size(); ++k)
      {
        const double
        val = (this->center.distance(points[k])<this->radius) ? 1. : 0.;
        if (this->selected==CutOffFunctionBase<dim>::no_component)
          values[k] = val;
        else
          {
            values[k] = 0;
            values[k](this->selected) = val;
          }
      }
  }


  template<int dim>
  CutOffFunctionW1<dim>::CutOffFunctionW1 (const double     r,
                                           const Point<dim> p,
                                           const unsigned int n_components,
                                           const unsigned int select)
    :
    CutOffFunctionBase<dim> (r, p, n_components, select)
  {}


  template<int dim>
  double
  CutOffFunctionW1<dim>::value (const Point<dim>   &p,
                                const unsigned int component) const
  {
    if (this->selected==CutOffFunctionBase<dim>::no_component
        ||
        component==this->selected)
      {
        const double d = this->center.distance(p);
        return ((d<this->radius) ? (this->radius-d) : 0.);
      }
    return 0.;
  }


  template<int dim>
  void
  CutOffFunctionW1<dim>::value_list (const std::vector<Point<dim> > &points,
                                     std::vector<double>            &values,
                                     const unsigned int component) const
  {
    Assert (values.size() == points.size(),
            ExcDimensionMismatch(values.size(), points.size()));

    if (this->selected==CutOffFunctionBase<dim>::no_component
        ||
        component==this->selected)
      for (unsigned int i=0; i<values.size(); ++i)
        {
          const double d = this->center.distance(points[i]);
          values[i] = ((d<this->radius) ? (this->radius-d) : 0.);
        }
    else
      std::fill (values.begin(), values.end(), 0.);
  }



  template<int dim>
  void
  CutOffFunctionW1<dim>::vector_value_list (
    const std::vector<Point<dim> > &points,
    std::vector<Vector<double> >           &values) const
  {
    Assert (values.size() == points.size(),
            ExcDimensionMismatch(values.size(), points.size()));

    for (unsigned int k=0; k<values.size(); ++k)
      {
        const double d = this->center.distance(points[k]);
        const double
        val = (d<this->radius) ? (this->radius-d) : 0.;
        if (this->selected==CutOffFunctionBase<dim>::no_component)
          values[k] = val;
        else
          {
            values[k] = 0;
            values[k](this->selected) = val;
          }
      }
  }


  template<int dim>
  CutOffFunctionCinfty<dim>::CutOffFunctionCinfty (const double     r,
                                                   const Point<dim> p,
                                                   const unsigned int n_components,
                                                   const unsigned int select)
    :
    CutOffFunctionBase<dim> (r, p, n_components, select)
  {}


  template<int dim>
  double
  CutOffFunctionCinfty<dim>::value (const Point<dim>   &p,
                                    const unsigned int component) const
  {
    if (this->selected==CutOffFunctionBase<dim>::no_component
        ||
        component==this->selected)
      {
        const double d = this->center.distance(p);
        const double r = this->radius;
        if (d>=r)
          return 0.;
        const double e = -r*r/(r*r-d*d);
        return  ((e<-50) ? 0. : numbers::E * exp(e));
      }
    return 0.;
  }


  template<int dim>
  void
  CutOffFunctionCinfty<dim>::value_list (const std::vector<Point<dim> > &points,
                                         std::vector<double>            &values,
                                         const unsigned int component) const
  {
    Assert (values.size() == points.size(),
            ExcDimensionMismatch(values.size(), points.size()));

    const double r = this->radius;

    if (this->selected==CutOffFunctionBase<dim>::no_component
        ||
        component==this->selected)
      for (unsigned int i=0; i<values.size(); ++i)
        {
          const double d = this->center.distance(points[i]);
          if (d>=r)
            {
              values[i] = 0.;
            }
          else
            {
              const double e = -r*r/(r*r-d*d);
              values[i] = (e<-50) ? 0. : numbers::E * exp(e);
            }
        }
    else
      std::fill (values.begin(), values.end(), 0.);
  }


  template<int dim>
  void
  CutOffFunctionCinfty<dim>::vector_value_list (
    const std::vector<Point<dim> > &points,
    std::vector<Vector<double> >           &values) const
  {
    Assert (values.size() == points.size(),
            ExcDimensionMismatch(values.size(), points.size()));

    for (unsigned int k=0; k<values.size(); ++k)
      {
        const double d = this->center.distance(points[k]);
        const double r = this->radius;
        double e = 0.;
        double val = 0.;
        if (d<this->radius)
          {
            e = -r*r/(r*r-d*d);
            if (e>-50)
              val = numbers::E * exp(e);
          }

        if (this->selected==CutOffFunctionBase<dim>::no_component)
          values[k] = val;
        else
          {
            values[k] = 0;
            values[k](this->selected) = val;
          }
      }
  }



  template<int dim>
  Tensor<1,dim>
  CutOffFunctionCinfty<dim>::gradient (const Point<dim>   &p,
                                       const unsigned int) const
  {
    const double d = this->center.distance(p);
    const double r = this->radius;
    if (d>=r)
      return Tensor<1,dim>();
    const double e = -d*d/(r-d)/(r+d);
    return  ((e<-50) ?
             Point<dim>() :
             (p-this->center)/d*(-2.0*r*r/pow(-r*r+d*d,2.0)*d*exp(e)));
  }


// explicit instantiations
  template class CutOffFunctionLinfty <1>;
  template class CutOffFunctionLinfty <2>;
  template class CutOffFunctionLinfty <3>;

  template class CutOffFunctionW1 <1>;
  template class CutOffFunctionW1 <2>;
  template class CutOffFunctionW1 <3>;

  template class CutOffFunctionCinfty <1>;
  template class CutOffFunctionCinfty <2>;
  template class CutOffFunctionCinfty <3>;
}

DEAL_II_NAMESPACE_CLOSE
