// ---------------------------------------------------------------------
//
// Copyright (C) 2015 by the deal.II authors
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


#include <deal.II/fe/fe_rannacher_turek.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/lac/vector.h>
#include <algorithm>

#include <sstream>


DEAL_II_NAMESPACE_OPEN


template <int dim>
FE_RannacherTurek<dim>::FE_RannacherTurek(const unsigned int degree,
                                          const unsigned int n_face_support_points) :
  FE_Poly<PolynomialsRannacherTurek<dim>, dim>(
    PolynomialsRannacherTurek<dim>(),
    FiniteElementData<dim>(this->get_dpo_vector(),
                           1,
                           2,
                           FiniteElementData<dim>::L2),
    std::vector<bool>(4, false), // restriction not implemented
    std::vector<ComponentMask>(4, std::vector<bool>(1, true))),
  degree(degree),
  n_face_support_points(n_face_support_points)
{
  Assert(dim == 2, ExcNotImplemented());
  Assert(degree == 0, ExcNotImplemented());
  this->initialize_support_points();
}



template <int dim>
std::vector<unsigned int> FE_RannacherTurek<dim>::get_dpo_vector()
{
  std::vector<unsigned int> dpo(dim + 1, 0);
  dpo[dim - 1] = 1;

  return dpo;
}



template <int dim>
std::string FE_RannacherTurek<dim>::get_name() const
{
  std::ostringstream namebuf;
  namebuf << "FE_RannacherTurek"
          << "<" << dim << ">"
          << "(" << this->degree << ", " << this->n_face_support_points << ")";
  return namebuf.str();
}



template <int dim>
FiniteElement<dim> *FE_RannacherTurek<dim>::clone() const
{
  return new FE_RannacherTurek<dim>(this->degree, this->n_face_support_points);
}



template <int dim>
void FE_RannacherTurek<dim>::initialize_support_points()
{
  Assert(dim == 2, ExcNotImplemented());
  dealii::QGauss<dim-1> face_quadrature(this->n_face_support_points);
  this->weights = face_quadrature.get_weights();
  this->generalized_support_points.resize(4*face_quadrature.size());
  for (unsigned int q = 0;
       q < face_quadrature.size();
       ++q)
    {
      this->generalized_support_points[0*face_quadrature.size() + q] =
        dealii::Point<dim>(0, 1 - face_quadrature.point(q)(0));
      this->generalized_support_points[1*face_quadrature.size() + q] =
        dealii::Point<dim>(1, 1 - face_quadrature.point(q)(0));
      this->generalized_support_points[2*face_quadrature.size() + q] =
        dealii::Point<dim>(face_quadrature.point(q)(0), 0);
      this->generalized_support_points[3*face_quadrature.size() + q] =
        dealii::Point<dim>(face_quadrature.point(q)(0), 1);
    }
}



template <int dim>
void FE_RannacherTurek<dim>::interpolate(
  std::vector<double> &local_dofs,
  const std::vector<double> &values) const
{
  AssertDimension(values.size(), this->generalized_support_points.size());
  AssertDimension(local_dofs.size(), this->dofs_per_cell);

  const unsigned int q_points_per_face = this->weights.size();
  std::fill(local_dofs.begin(), local_dofs.end(), 0.0);

  std::vector<double>::const_iterator value = values.begin();
  for (unsigned int face = 0;
       face < dealii::GeometryInfo<dim>::faces_per_cell;
       ++face)
    {
      for (unsigned int q = 0;
           q < q_points_per_face;
           ++q)
        {
          local_dofs[face] += (*value) * this->weights[q];
          ++value;
        }
    }
}



template <int dim>
void FE_RannacherTurek<dim>::interpolate(
  std::vector<double> &local_dofs,
  const std::vector<Vector<double> > &values,
  unsigned int offset) const
{
  AssertDimension(values.size(), this->generalized_support_points.size());
  AssertDimension(local_dofs.size(), this->dofs_per_cell);

  // extract component at offset and call scalar version of this function
  std::vector<double> scalar_values(values.size());
  for (unsigned int q = 0; q < values.size(); ++q)
    {
      scalar_values[q] = values[q][offset];
    }
  this->interpolate(local_dofs, scalar_values);
}



template <int dim>
void FE_RannacherTurek<dim>::interpolate(
  std::vector<double> &local_dofs,
  const VectorSlice<const std::vector<std::vector<double> > > &values) const
{
  AssertDimension(values.size(), 1);
  AssertDimension(values[0].size(), this->generalized_support_points.size());
  AssertDimension(local_dofs.size(), this->dofs_per_cell);

  // convert data structure to use scalar version of this function
  std::vector<double> scalar_values(values[0].size());
  for (unsigned int q = 0; q < values[0].size(); ++q)
    {
      scalar_values[q] = values[0][q];
    }
  this->interpolate(local_dofs, scalar_values);
}



// explicit instantiations
#include "fe_rannacher_turek.inst"

DEAL_II_NAMESPACE_CLOSE
