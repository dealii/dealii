// ---------------------------------------------------------------------
//
// Copyright (C) 2000 - 2018 by the deal.II authors
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

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_nothing.h>
#include <deal.II/fe/fe_q_iso_q1.h>

#include <deal.II/lac/vector.h>

#include <memory>
#include <sstream>
#include <vector>

DEAL_II_NAMESPACE_OPEN



template <int dim, int spacedim>
FE_Q_iso_Q1<dim, spacedim>::FE_Q_iso_Q1(const unsigned int subdivisions)
  : FE_Q_Base<dim, spacedim>(
      TensorProductPolynomials<dim, Polynomials::PiecewisePolynomial<double>>(
        Polynomials::generate_complete_Lagrange_basis_on_subdivisions(
          subdivisions,
          1)),
      FiniteElementData<dim>(this->get_dpo_vector(subdivisions),
                             1,
                             subdivisions,
                             FiniteElementData<dim>::H1),
      std::vector<bool>(1, false))
{
  Assert(subdivisions > 0,
         ExcMessage("This element can only be used with a positive number of "
                    "subelements"));

  QTrapezoid<1> trapez;
  QIterated<1>  points(trapez, subdivisions);

  this->initialize(points.get_points());
}



template <int dim, int spacedim>
std::string
FE_Q_iso_Q1<dim, spacedim>::get_name() const
{
  // note that the FETools::get_fe_by_name function depends on the
  // particular format of the string this function returns, so they have to be
  // kept in sync

  std::ostringstream namebuf;
  namebuf << "FE_Q_iso_Q1<" << Utilities::dim_string(dim, spacedim) << ">("
          << this->degree << ")";
  return namebuf.str();
}



template <int dim, int spacedim>
void
FE_Q_iso_Q1<dim, spacedim>::
  convert_generalized_support_point_values_to_dof_values(
    const std::vector<Vector<double>> &support_point_values,
    std::vector<double> &              nodal_values) const
{
  AssertDimension(support_point_values.size(),
                  this->get_unit_support_points().size());
  AssertDimension(support_point_values.size(), nodal_values.size());
  AssertDimension(this->n_dofs_per_cell(), nodal_values.size());

  for (unsigned int i = 0; i < this->n_dofs_per_cell(); ++i)
    {
      AssertDimension(support_point_values[i].size(), 1);

      nodal_values[i] = support_point_values[i](0);
    }
}



template <int dim, int spacedim>
std::unique_ptr<FiniteElement<dim, spacedim>>
FE_Q_iso_Q1<dim, spacedim>::clone() const
{
  return std::make_unique<FE_Q_iso_Q1<dim, spacedim>>(*this);
}



template <int dim, int spacedim>
FiniteElementDomination::Domination
FE_Q_iso_Q1<dim, spacedim>::compare_for_domination(
  const FiniteElement<dim, spacedim> &fe_other,
  const unsigned int                  codim) const
{
  Assert(codim <= dim, ExcImpossibleInDim(dim));
  (void)codim;

  // vertex/line/face domination
  // (if fe_other is derived from FE_DGQ)
  // ------------------------------------
  if (codim > 0)
    if (dynamic_cast<const FE_DGQ<dim, spacedim> *>(&fe_other) != nullptr)
      // there are no requirements between continuous and discontinuous elements
      return FiniteElementDomination::no_requirements;

  // vertex/line/face domination
  // (if fe_other is not derived from FE_DGQ)
  // & cell domination
  // ----------------------------------------
  if (const FE_Q_iso_Q1<dim, spacedim> *fe_q_iso_q1_other =
        dynamic_cast<const FE_Q_iso_Q1<dim, spacedim> *>(&fe_other))
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
  else if (const FE_Nothing<dim> *fe_nothing =
             dynamic_cast<const FE_Nothing<dim> *>(&fe_other))
    {
      if (fe_nothing->is_dominating())
        return FiniteElementDomination::other_element_dominates;
      else
        // the FE_Nothing has no degrees of freedom and it is typically used
        // in a context where we don't require any continuity along the
        // interface
        return FiniteElementDomination::no_requirements;
    }

  Assert(false, ExcNotImplemented());
  return FiniteElementDomination::neither_element_dominates;
}


// explicit instantiations
#include "fe_q_iso_q1.inst"

DEAL_II_NAMESPACE_CLOSE
