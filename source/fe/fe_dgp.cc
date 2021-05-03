// ---------------------------------------------------------------------
//
// Copyright (C) 2002 - 2018 by the deal.II authors
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


#include <deal.II/fe/fe_dgp.h>
#include <deal.II/fe/fe_nothing.h>
#include <deal.II/fe/fe_tools.h>

#include <memory>
#include <sstream>


DEAL_II_NAMESPACE_OPEN

template <int dim, int spacedim>
FE_DGP<dim, spacedim>::FE_DGP(const unsigned int degree)
  : FE_Poly<dim, spacedim>(
      PolynomialSpace<dim>(
        Polynomials::Legendre::generate_complete_basis(degree)),
      FiniteElementData<dim>(get_dpo_vector(degree),
                             1,
                             degree,
                             FiniteElementData<dim>::L2),
      std::vector<bool>(
        FiniteElementData<dim>(get_dpo_vector(degree), 1, degree)
          .n_dofs_per_cell(),
        true),
      std::vector<ComponentMask>(
        FiniteElementData<dim>(get_dpo_vector(degree), 1, degree)
          .n_dofs_per_cell(),
        std::vector<bool>(1, true)))
{
  // Reinit the vectors of restriction and prolongation matrices to the right
  // sizes
  this->reinit_restriction_and_prolongation_matrices();
  // Fill prolongation matrices with embedding operators
  if (dim == spacedim)
    {
      FETools::compute_embedding_matrices(*this, this->prolongation);
      // Fill restriction matrices with L2-projection
      FETools::compute_projection_matrices(*this, this->restriction);
    }
}


template <int dim, int spacedim>
std::string
FE_DGP<dim, spacedim>::get_name() const
{
  // note that the FETools::get_fe_by_name function depends on the
  // particular format of the string this function returns, so they have to be
  // kept in sync

  std::ostringstream namebuf;
  namebuf << "FE_DGP<" << Utilities::dim_string(dim, spacedim) << ">("
          << this->degree << ")";

  return namebuf.str();
}



template <int dim, int spacedim>
std::unique_ptr<FiniteElement<dim, spacedim>>
FE_DGP<dim, spacedim>::clone() const
{
  return std::make_unique<FE_DGP<dim, spacedim>>(*this);
}



//---------------------------------------------------------------------------
// Auxiliary functions
//---------------------------------------------------------------------------


template <int dim, int spacedim>
std::vector<unsigned int>
FE_DGP<dim, spacedim>::get_dpo_vector(const unsigned int deg)
{
  std::vector<unsigned int> dpo(dim + 1, 0U);
  dpo[dim] = deg + 1;
  for (unsigned int i = 1; i < dim; ++i)
    {
      dpo[dim] *= deg + 1 + i;
      dpo[dim] /= i + 1;
    }
  return dpo;
}



template <int dim, int spacedim>
void
FE_DGP<dim, spacedim>::get_face_interpolation_matrix(
  const FiniteElement<dim, spacedim> &x_source_fe,
  FullMatrix<double> &                interpolation_matrix,
  const unsigned int) const
{
  // this is only implemented, if the source FE is also a DGP element. in that
  // case, both elements have no dofs on their faces and the face
  // interpolation matrix is necessarily empty -- i.e. there isn't much we
  // need to do here.
  (void)interpolation_matrix;
  using FE    = FiniteElement<dim, spacedim>;
  using FEDGP = FE_DGP<dim, spacedim>;
  AssertThrow((x_source_fe.get_name().find("FE_DGP<") == 0) ||
                (dynamic_cast<const FEDGP *>(&x_source_fe) != nullptr),
              typename FE::ExcInterpolationNotImplemented());

  Assert(interpolation_matrix.m() == 0,
         ExcDimensionMismatch(interpolation_matrix.m(), 0));
  Assert(interpolation_matrix.n() == 0,
         ExcDimensionMismatch(interpolation_matrix.n(), 0));
}



template <int dim, int spacedim>
void
FE_DGP<dim, spacedim>::get_subface_interpolation_matrix(
  const FiniteElement<dim, spacedim> &x_source_fe,
  const unsigned int,
  FullMatrix<double> &interpolation_matrix,
  const unsigned int) const
{
  // this is only implemented, if the source FE is also a DGP element. in that
  // case, both elements have no dofs on their faces and the face
  // interpolation matrix is necessarily empty -- i.e. there isn't much we
  // need to do here.
  (void)interpolation_matrix;
  using FE    = FiniteElement<dim, spacedim>;
  using FEDGP = FE_DGP<dim, spacedim>;
  AssertThrow((x_source_fe.get_name().find("FE_DGP<") == 0) ||
                (dynamic_cast<const FEDGP *>(&x_source_fe) != nullptr),
              typename FE::ExcInterpolationNotImplemented());

  Assert(interpolation_matrix.m() == 0,
         ExcDimensionMismatch(interpolation_matrix.m(), 0));
  Assert(interpolation_matrix.n() == 0,
         ExcDimensionMismatch(interpolation_matrix.n(), 0));
}



template <int dim, int spacedim>
bool
FE_DGP<dim, spacedim>::hp_constraints_are_implemented() const
{
  return true;
}



template <int dim, int spacedim>
std::vector<std::pair<unsigned int, unsigned int>>
FE_DGP<dim, spacedim>::hp_vertex_dof_identities(
  const FiniteElement<dim, spacedim> &fe_other) const
{
  // there are no such constraints for DGP elements at all
  if (dynamic_cast<const FE_DGP<dim, spacedim> *>(&fe_other) != nullptr)
    return std::vector<std::pair<unsigned int, unsigned int>>();
  else
    {
      Assert(false, ExcNotImplemented());
      return std::vector<std::pair<unsigned int, unsigned int>>();
    }
}



template <int dim, int spacedim>
std::vector<std::pair<unsigned int, unsigned int>>
FE_DGP<dim, spacedim>::hp_line_dof_identities(
  const FiniteElement<dim, spacedim> &fe_other) const
{
  // there are no such constraints for DGP elements at all
  if (dynamic_cast<const FE_DGP<dim, spacedim> *>(&fe_other) != nullptr)
    return std::vector<std::pair<unsigned int, unsigned int>>();
  else
    {
      Assert(false, ExcNotImplemented());
      return std::vector<std::pair<unsigned int, unsigned int>>();
    }
}



template <int dim, int spacedim>
std::vector<std::pair<unsigned int, unsigned int>>
FE_DGP<dim, spacedim>::hp_quad_dof_identities(
  const FiniteElement<dim, spacedim> &fe_other,
  const unsigned int) const
{
  // there are no such constraints for DGP elements at all
  if (dynamic_cast<const FE_DGP<dim, spacedim> *>(&fe_other) != nullptr)
    return std::vector<std::pair<unsigned int, unsigned int>>();
  else
    {
      Assert(false, ExcNotImplemented());
      return std::vector<std::pair<unsigned int, unsigned int>>();
    }
}



template <int dim, int spacedim>
FiniteElementDomination::Domination
FE_DGP<dim, spacedim>::compare_for_domination(
  const FiniteElement<dim, spacedim> &fe_other,
  const unsigned int                  codim) const
{
  Assert(codim <= dim, ExcImpossibleInDim(dim));

  // vertex/line/face domination
  // ---------------------------
  if (codim > 0)
    // this is a discontinuous element, so by definition there will
    // be no constraints wherever this element comes together with
    // any other kind of element
    return FiniteElementDomination::no_requirements;

  // cell domination
  // ---------------
  if (const FE_DGP<dim, spacedim> *fe_dgp_other =
        dynamic_cast<const FE_DGP<dim, spacedim> *>(&fe_other))
    {
      if (this->degree < fe_dgp_other->degree)
        return FiniteElementDomination::this_element_dominates;
      else if (this->degree == fe_dgp_other->degree)
        return FiniteElementDomination::either_element_can_dominate;
      else
        return FiniteElementDomination::other_element_dominates;
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



template <int dim, int spacedim>
bool
FE_DGP<dim, spacedim>::has_support_on_face(const unsigned int,
                                           const unsigned int) const
{
  // all shape functions have support on all faces
  return true;
}



template <int dim, int spacedim>
std::pair<Table<2, bool>, std::vector<unsigned int>>
FE_DGP<dim, spacedim>::get_constant_modes() const
{
  Table<2, bool> constant_modes(1, this->n_dofs_per_cell());
  constant_modes(0, 0) = true;
  return std::pair<Table<2, bool>, std::vector<unsigned int>>(
    constant_modes, std::vector<unsigned int>(1, 0));
}



template <int dim, int spacedim>
std::size_t
FE_DGP<dim, spacedim>::memory_consumption() const
{
  Assert(false, ExcNotImplemented());
  return 0;
}



// explicit instantiations
#include "fe_dgp.inst"


DEAL_II_NAMESPACE_CLOSE
