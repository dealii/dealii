// ---------------------------------------------------------------------
//
// Copyright (C) 2015 - 2020 by the deal.II authors
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


#include <deal.II/base/quadrature_lib.h>

#include <deal.II/fe/fe_rannacher_turek.h>

#include <deal.II/lac/vector.h>

#include <algorithm>
#include <memory>
#include <sstream>


DEAL_II_NAMESPACE_OPEN


template <int dim>
FE_RannacherTurek<dim>::FE_RannacherTurek(
  const unsigned int order,
  const unsigned int n_face_support_points)
  : FE_Poly<dim>(PolynomialsRannacherTurek<dim>(),
                 FiniteElementData<dim>(this->get_dpo_vector(),
                                        1,
                                        2,
                                        FiniteElementData<dim>::L2),
                 std::vector<bool>(4, false), // restriction not implemented
                 std::vector<ComponentMask>(4, std::vector<bool>(1, true)))
  , order(order)
  , n_face_support_points(n_face_support_points)
{
  Assert(dim == 2, ExcNotImplemented());
  Assert(order == 0, ExcNotImplemented());
  this->initialize_support_points();
}



template <int dim>
std::vector<unsigned int>
FE_RannacherTurek<dim>::get_dpo_vector()
{
  std::vector<unsigned int> dpo(dim + 1, 0);
  dpo[dim - 1] = 1;

  return dpo;
}



template <int dim>
std::string
FE_RannacherTurek<dim>::get_name() const
{
  std::ostringstream namebuf;
  namebuf << "FE_RannacherTurek"
          << "<" << dim << ">"
          << "(" << this->order << ", " << this->n_face_support_points << ")";
  return namebuf.str();
}



template <int dim>
std::unique_ptr<FiniteElement<dim, dim>>
FE_RannacherTurek<dim>::clone() const
{
  return std::make_unique<FE_RannacherTurek<dim>>(this->order,
                                                  this->n_face_support_points);
}



template <int dim>
void
FE_RannacherTurek<dim>::initialize_support_points()
{
  Assert(dim == 2, ExcNotImplemented());
  dealii::QGauss<dim - 1> face_quadrature(this->n_face_support_points);
  this->weights = face_quadrature.get_weights();
  this->generalized_support_points.resize(4 * face_quadrature.size());
  for (unsigned int q = 0; q < face_quadrature.size(); ++q)
    {
      this->generalized_support_points[0 * face_quadrature.size() + q] =
        dealii::Point<dim>(0, 1 - face_quadrature.point(q)(0));
      this->generalized_support_points[1 * face_quadrature.size() + q] =
        dealii::Point<dim>(1, 1 - face_quadrature.point(q)(0));
      this->generalized_support_points[2 * face_quadrature.size() + q] =
        dealii::Point<dim>(face_quadrature.point(q)(0), 0);
      this->generalized_support_points[3 * face_quadrature.size() + q] =
        dealii::Point<dim>(face_quadrature.point(q)(0), 1);
    }
}



template <int dim>
void
FE_RannacherTurek<dim>::convert_generalized_support_point_values_to_dof_values(
  const std::vector<Vector<double>> &support_point_values,
  std::vector<double> &              nodal_values) const
{
  AssertDimension(support_point_values.size(),
                  this->generalized_support_points.size());
  AssertDimension(nodal_values.size(), this->dofs_per_cell);

  const unsigned int q_points_per_face = this->weights.size();
  std::fill(nodal_values.begin(), nodal_values.end(), 0.0);

  std::vector<Vector<double>>::const_iterator value =
    support_point_values.begin();
  for (const unsigned int face : dealii::GeometryInfo<dim>::face_indices())
    {
      for (unsigned int q = 0; q < q_points_per_face; ++q)
        {
          nodal_values[face] += (*value)[0] * this->weights[q];
          ++value;
        }
    }
}



// explicit instantiations
#include "fe_rannacher_turek.inst"

DEAL_II_NAMESPACE_CLOSE
