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
#include <deal.II/fe/fe_q_iso_q1.h>
#include <deal.II/fe/fe_nothing.h>

#include <vector>
#include <sstream>

DEAL_II_NAMESPACE_OPEN





template <int dim, int spacedim>
FE_Q_iso_Q1<dim,spacedim>::FE_Q_iso_Q1 (const unsigned int subdivisions)
  :
  FE_Q_Base<TensorProductPolynomials<dim,Polynomials::PiecewisePolynomial<double> >, dim, spacedim> (
    TensorProductPolynomials<dim,Polynomials::PiecewisePolynomial<double> >
    (Polynomials::generate_complete_Lagrange_basis_on_subdivisions(subdivisions, 1)),
    FiniteElementData<dim>(this->get_dpo_vector(subdivisions),
                           1, subdivisions,
                           FiniteElementData<dim>::H1),
    std::vector<bool> (1, false))
{
  Assert (subdivisions > 0,
          ExcMessage ("This element can only be used with a positive number of "
                      "subelements"));

  QTrapez<1> trapez;
  QIterated<1> points (trapez, subdivisions);

  this->initialize(points.get_points());
}



template <int dim, int spacedim>
std::string
FE_Q_iso_Q1<dim,spacedim>::get_name () const
{
  // note that the FETools::get_fe_from_name function depends on the
  // particular format of the string this function returns, so they have to be
  // kept in synch

  std::ostringstream namebuf;
  namebuf << "FE_Q_iso_Q1<"
          << Utilities::dim_string(dim,spacedim)
          << ">(" << this->degree << ")";
  return namebuf.str();
}



template <int dim, int spacedim>
FiniteElement<dim,spacedim> *
FE_Q_iso_Q1<dim,spacedim>::clone() const
{
  return new FE_Q_iso_Q1<dim,spacedim>(*this);
}



template <int dim, int spacedim>
FiniteElementDomination::Domination
FE_Q_iso_Q1<dim,spacedim>::
compare_for_face_domination (const FiniteElement<dim,spacedim> &fe_other) const
{
  if (const FE_Q_iso_Q1<dim,spacedim> *fe_q_iso_q1_other
      = dynamic_cast<const FE_Q_iso_Q1<dim,spacedim>*>(&fe_other))
    {
      // different behavior as in FE_Q: as FE_Q_iso_Q1(2) is not a subspace of
      // FE_Q_iso_Q1(3), need that the element degrees are multiples of each
      // other
      if (this->degree < fe_q_iso_q1_other->degree &&
          fe_q_iso_q1_other->degree % this->degree == 0)
        return FiniteElementDomination::this_element_dominates;
      else if (this->degree == fe_q_iso_q1_other->degree)
        return FiniteElementDomination::either_element_can_dominate;
      else if (this->degree > fe_q_iso_q1_other->degree &&
               this->degree % fe_q_iso_q1_other->degree == 0)
        return FiniteElementDomination::other_element_dominates;
      else
        return FiniteElementDomination::neither_element_dominates;
    }
  else if (dynamic_cast<const FE_Nothing<dim>*>(&fe_other) != 0)
    {
      // the FE_Nothing has no degrees of freedom and it is typically used in
      // a context where we don't require any continuity along the interface
      return FiniteElementDomination::no_requirements;
    }

  Assert (false, ExcNotImplemented());
  return FiniteElementDomination::neither_element_dominates;
}


// explicit instantiations
#include "fe_q_iso_q1.inst"

DEAL_II_NAMESPACE_CLOSE
