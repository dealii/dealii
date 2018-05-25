// ---------------------------------------------------------------------
//
// Copyright (C) 2000 - 2017 by the deal.II authors
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


#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/std_cxx14/memory.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/lac/vector.h>

#include <sstream>
#include <vector>

DEAL_II_NAMESPACE_OPEN


namespace internal
{
  namespace FE_Q
  {
    namespace
    {
      std::vector<Point<1>>
      get_QGaussLobatto_points(const unsigned int degree)
      {
        if (degree > 0)
          return QGaussLobatto<1>(degree + 1).get_points();
        else
          {
            typedef dealii::FE_Q_Base<TensorProductPolynomials<1>, 1, 1> FEQ;
            AssertThrow(false, FEQ::ExcFEQCannotHaveDegree0());
          }
        return std::vector<Point<1>>();
      }
    } // namespace
  }   // namespace FE_Q
} // namespace internal



template <int dim, int spacedim>
FE_Q<dim, spacedim>::FE_Q(const unsigned int degree) :
  FE_Q_Base<TensorProductPolynomials<dim>, dim, spacedim>(
    TensorProductPolynomials<dim>(Polynomials::generate_complete_Lagrange_basis(
      internal::FE_Q::get_QGaussLobatto_points(degree))),
    FiniteElementData<dim>(this->get_dpo_vector(degree),
                           1,
                           degree,
                           FiniteElementData<dim>::H1),
    std::vector<bool>(1, false))
{
  this->initialize(internal::FE_Q::get_QGaussLobatto_points(degree));
}



template <int dim, int spacedim>
FE_Q<dim, spacedim>::FE_Q(const Quadrature<1> &points) :
  FE_Q_Base<TensorProductPolynomials<dim>, dim, spacedim>(
    TensorProductPolynomials<dim>(
      Polynomials::generate_complete_Lagrange_basis(points.get_points())),
    FiniteElementData<dim>(this->get_dpo_vector(points.size() - 1),
                           1,
                           points.size() - 1,
                           FiniteElementData<dim>::H1),
    std::vector<bool>(1, false))
{
  this->initialize(points.get_points());
}



template <int dim, int spacedim>
std::string
FE_Q<dim, spacedim>::get_name() const
{
  // note that the FETools::get_fe_by_name function depends on the
  // particular format of the string this function returns, so they have to be
  // kept in synch

  std::ostringstream  namebuf;
  bool                equidistant = true;
  std::vector<double> points(this->degree + 1);

  // Decode the support points in one coordinate direction.
  std::vector<unsigned int> lexicographic =
    this->poly_space.get_numbering_inverse();
  for (unsigned int j = 0; j <= this->degree; j++)
    points[j] = this->unit_support_points[lexicographic[j]][0];

  // Check whether the support points are equidistant.
  for (unsigned int j = 0; j <= this->degree; j++)
    if (std::fabs(points[j] - (double)j / this->degree) > 1e-15)
      {
        equidistant = false;
        break;
      }

  if (equidistant == true)
    {
      if (this->degree > 2)
        namebuf << "FE_Q<" << Utilities::dim_string(dim, spacedim)
                << ">(QIterated(QTrapez()," << this->degree << "))";
      else
        namebuf << "FE_Q<" << Utilities::dim_string(dim, spacedim) << ">("
                << this->degree << ")";
    }
  else
    {
      // Check whether the support points come from QGaussLobatto.
      const QGaussLobatto<1> points_gl(this->degree + 1);
      bool                   gauss_lobatto = true;
      for (unsigned int j = 0; j <= this->degree; j++)
        if (points[j] != points_gl.point(j)(0))
          {
            gauss_lobatto = false;
            break;
          }
      if (gauss_lobatto == true)
        namebuf << "FE_Q<" << Utilities::dim_string(dim, spacedim) << ">("
                << this->degree << ")";
      else
        namebuf << "FE_Q<" << Utilities::dim_string(dim, spacedim)
                << ">(QUnknownNodes(" << this->degree << "))";
    }
  return namebuf.str();
}



template <int dim, int spacedim>
void
FE_Q<dim, spacedim>::convert_generalized_support_point_values_to_dof_values(
  const std::vector<Vector<double>> &support_point_values,
  std::vector<double> &              nodal_values) const
{
  AssertDimension(support_point_values.size(),
                  this->get_unit_support_points().size());
  AssertDimension(support_point_values.size(), nodal_values.size());
  AssertDimension(this->dofs_per_cell, nodal_values.size());

  for (unsigned int i = 0; i < this->dofs_per_cell; ++i)
    {
      AssertDimension(support_point_values[i].size(), 1);

      nodal_values[i] = support_point_values[i](0);
    }
}



template <int dim, int spacedim>
std::unique_ptr<FiniteElement<dim, spacedim>>
FE_Q<dim, spacedim>::clone() const
{
  return std_cxx14::make_unique<FE_Q<dim, spacedim>>(*this);
}


// explicit instantiations
#include "fe_q.inst"

DEAL_II_NAMESPACE_CLOSE
