// ---------------------------------------------------------------------
//
// Copyright (C) 2000 - 2013 by the deal.II authors
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
#include <deal.II/fe/fe_q.h>

#include <vector>
#include <sstream>

DEAL_II_NAMESPACE_OPEN



template <int dim, int spacedim>
FE_Q<dim,spacedim>::FE_Q (const unsigned int degree)
  :
  FE_Q_Base<TensorProductPolynomials<dim>, dim, spacedim> (
    TensorProductPolynomials<dim>(Polynomials::LagrangeEquidistant::generate_complete_basis(degree)),
    FiniteElementData<dim>(this->get_dpo_vector(degree),
                           1, degree,
                           FiniteElementData<dim>::H1),
    std::vector<bool> (1, false))
{
  Assert (degree > 0,
          ExcMessage ("This element can only be used for polynomial degrees "
                      "greater than zero"));
  std::vector<Point<1> > support_points_1d(degree+1);
  for (unsigned int i=0; i<=degree; ++i)
    support_points_1d[i][0] = static_cast<double>(i)/degree;

  this->initialize(support_points_1d);
}



template <int dim, int spacedim>
FE_Q<dim,spacedim>::FE_Q (const Quadrature<1> &points)
  :
  FE_Q_Base<TensorProductPolynomials<dim>, dim, spacedim> (
    TensorProductPolynomials<dim>(Polynomials::generate_complete_Lagrange_basis(points.get_points())),
    FiniteElementData<dim>(this->get_dpo_vector(points.size()-1),
                           1, points.size()-1,
                           FiniteElementData<dim>::H1),
    std::vector<bool> (1, false))
{
  const unsigned int degree = points.size()-1;
  Assert (degree > 0,
          ExcMessage ("This element can only be used for polynomial degrees "
                      "at least zero"));

  this->initialize(points.get_points());
}



template <int dim, int spacedim>
std::string
FE_Q<dim,spacedim>::get_name () const
{
  // note that the FETools::get_fe_from_name function depends on the
  // particular format of the string this function returns, so they have to be
  // kept in synch

  std::ostringstream namebuf;
  bool equidistant = true;
  std::vector<double> points(this->degree+1);

  // Decode the support points in one coordinate direction.
  std::vector<unsigned int> lexicographic = this->poly_space.get_numbering_inverse();
  for (unsigned int j=0; j<=this->degree; j++)
    points[j] = this->unit_support_points[lexicographic[j]][0];

  // Check whether the support points are equidistant.
  for (unsigned int j=0; j<=this->degree; j++)
    if (std::fabs(points[j] - (double)j/this->degree) > 1e-15)
      {
        equidistant = false;
        break;
      }

  if (equidistant == true)
    namebuf << "FE_Q<"
            << Utilities::dim_string(dim,spacedim)
            << ">(" << this->degree << ")";
  else
    {
      // Check whether the support points come from QGaussLobatto.
      const QGaussLobatto<1> points_gl(this->degree+1);
      bool gauss_lobatto = true;
      for (unsigned int j=0; j<=this->degree; j++)
        if (points[j] != points_gl.point(j)(0))
          {
            gauss_lobatto = false;
            break;
          }
      if (gauss_lobatto == true)
        namebuf << "FE_Q<"
                << Utilities::dim_string(dim,spacedim)
                << ">(QGaussLobatto(" << this->degree+1 << "))";
      else
        namebuf << "FE_Q<"
                << Utilities::dim_string(dim,spacedim)
                << ">(QUnknownNodes(" << this->degree << "))";
    }
  return namebuf.str();
}



template <int dim, int spacedim>
FiniteElement<dim,spacedim> *
FE_Q<dim,spacedim>::clone() const
{
  return new FE_Q<dim,spacedim>(*this);
}


// explicit instantiations
#include "fe_q.inst"

DEAL_II_NAMESPACE_CLOSE
