// ---------------------------------------------------------------------
//
// Copyright (C) 2002 - 2013 by the deal.II authors
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


#include <deal.II/fe/fe_dgp.h>
#include <deal.II/fe/fe_tools.h>

#include <sstream>

DEAL_II_NAMESPACE_OPEN

template <int dim, int spacedim>
FE_DGP<dim,spacedim>::FE_DGP (const unsigned int degree)
  :
  FE_Poly<PolynomialSpace<dim>, dim, spacedim> (
    PolynomialSpace<dim>(Polynomials::Legendre::generate_complete_basis(degree)),
    FiniteElementData<dim>(get_dpo_vector(degree), 1, degree, FiniteElementData<dim>::L2),
    std::vector<bool>(FiniteElementData<dim>(get_dpo_vector(degree), 1, degree).dofs_per_cell,true),
    std::vector<ComponentMask>(FiniteElementData<dim>(
                                 get_dpo_vector(degree), 1, degree).dofs_per_cell, std::vector<bool>(1,true)))
{
  // Reinit the vectors of restriction and prolongation matrices to the right
  // sizes
  this->reinit_restriction_and_prolongation_matrices();
  // Fill prolongation matrices with embedding operators
  if (dim == spacedim)
    {
      FETools::compute_embedding_matrices (*this, this->prolongation);
      // Fill restriction matrices with L2-projection
      FETools::compute_projection_matrices (*this, this->restriction);
    }
}


template <int dim, int spacedim>
std::string
FE_DGP<dim,spacedim>::get_name () const
{
  // note that the FETools::get_fe_from_name function depends on the
  // particular format of the string this function returns, so they have to be
  // kept in synch

  std::ostringstream namebuf;
  namebuf << "FE_DGP<"
          << Utilities::dim_string(dim,spacedim)
          << ">(" << this->degree << ")";

  return namebuf.str();
}



template <int dim, int spacedim>
FiniteElement<dim,spacedim> *
FE_DGP<dim,spacedim>::clone() const
{
  return new FE_DGP<dim,spacedim>(*this);
}



//---------------------------------------------------------------------------
// Auxiliary functions
//---------------------------------------------------------------------------


template <int dim, int spacedim>
std::vector<unsigned int>
FE_DGP<dim,spacedim>::get_dpo_vector (const unsigned int deg)
{
  std::vector<unsigned int> dpo(dim+1, 0U);
  dpo[dim] = deg+1;
  for (unsigned int i=1; i<dim; ++i)
    {
      dpo[dim] *= deg+1+i;
      dpo[dim] /= i+1;
    }
  return dpo;
}



template <int dim, int spacedim>
void
FE_DGP<dim,spacedim>::
get_face_interpolation_matrix (const FiniteElement<dim,spacedim> &x_source_fe,
                               FullMatrix<double>       &interpolation_matrix) const
{
  // this is only implemented, if the source FE is also a DGP element. in that
  // case, both elements have no dofs on their faces and the face
  // interpolation matrix is necessarily empty -- i.e. there isn't much we
  // need to do here.
  typedef FiniteElement<dim,spacedim> FE;
  typedef FE_DGP<dim,spacedim> FEDGP;
  AssertThrow ((x_source_fe.get_name().find ("FE_DGP<") == 0)
               ||
               (dynamic_cast<const FEDGP *>(&x_source_fe) != 0),
               typename FE::
               ExcInterpolationNotImplemented());

  Assert (interpolation_matrix.m() == 0,
          ExcDimensionMismatch (interpolation_matrix.m(),
                                0));
  Assert (interpolation_matrix.n() == 0,
          ExcDimensionMismatch (interpolation_matrix.n(),
                                0));
}



template <int dim, int spacedim>
void
FE_DGP<dim,spacedim>::
get_subface_interpolation_matrix (const FiniteElement<dim,spacedim> &x_source_fe,
                                  const unsigned int ,
                                  FullMatrix<double>           &interpolation_matrix) const
{
  // this is only implemented, if the source FE is also a DGP element. in that
  // case, both elements have no dofs on their faces and the face
  // interpolation matrix is necessarily empty -- i.e. there isn't much we
  // need to do here.
  typedef FiniteElement<dim,spacedim> FE;
  typedef FE_DGP<dim,spacedim> FEDGP;
  AssertThrow ((x_source_fe.get_name().find ("FE_DGP<") == 0)
               ||
               (dynamic_cast<const FEDGP *>(&x_source_fe) != 0),
               typename FE::
               ExcInterpolationNotImplemented());

  Assert (interpolation_matrix.m() == 0,
          ExcDimensionMismatch (interpolation_matrix.m(),
                                0));
  Assert (interpolation_matrix.n() == 0,
          ExcDimensionMismatch (interpolation_matrix.n(),
                                0));
}



template <int dim, int spacedim>
bool
FE_DGP<dim,spacedim>::hp_constraints_are_implemented () const
{
  return true;
}



template <int dim, int spacedim>
std::vector<std::pair<unsigned int, unsigned int> >
FE_DGP<dim,spacedim>::
hp_vertex_dof_identities (const FiniteElement<dim,spacedim> &fe_other) const
{
  // there are no such constraints for DGP elements at all
  if (dynamic_cast<const FE_DGP<dim,spacedim>*>(&fe_other) != 0)
    return
      std::vector<std::pair<unsigned int, unsigned int> > ();
  else
    {
      Assert (false, ExcNotImplemented());
      return std::vector<std::pair<unsigned int, unsigned int> > ();
    }
}



template <int dim, int spacedim>
std::vector<std::pair<unsigned int, unsigned int> >
FE_DGP<dim,spacedim>::
hp_line_dof_identities (const FiniteElement<dim,spacedim> &fe_other) const
{
  // there are no such constraints for DGP elements at all
  if (dynamic_cast<const FE_DGP<dim,spacedim>*>(&fe_other) != 0)
    return
      std::vector<std::pair<unsigned int, unsigned int> > ();
  else
    {
      Assert (false, ExcNotImplemented());
      return std::vector<std::pair<unsigned int, unsigned int> > ();
    }
}



template <int dim, int spacedim>
std::vector<std::pair<unsigned int, unsigned int> >
FE_DGP<dim,spacedim>::
hp_quad_dof_identities (const FiniteElement<dim,spacedim>        &fe_other) const
{
  // there are no such constraints for DGP elements at all
  if (dynamic_cast<const FE_DGP<dim,spacedim>*>(&fe_other) != 0)
    return
      std::vector<std::pair<unsigned int, unsigned int> > ();
  else
    {
      Assert (false, ExcNotImplemented());
      return std::vector<std::pair<unsigned int, unsigned int> > ();
    }
}



template <int dim, int spacedim>
FiniteElementDomination::Domination
FE_DGP<dim,spacedim>::compare_for_face_domination (const FiniteElement<dim,spacedim> &fe_other) const
{
  // check whether both are discontinuous elements, see the description of
  // FiniteElementDomination::Domination
  if (dynamic_cast<const FE_DGP<dim,spacedim>*>(&fe_other) != 0)
    return FiniteElementDomination::no_requirements;

  Assert (false, ExcNotImplemented());
  return FiniteElementDomination::neither_element_dominates;
}



template <int dim, int spacedim>
bool
FE_DGP<dim,spacedim>::has_support_on_face (const unsigned int,
                                           const unsigned int) const
{
  // all shape functions have support on all faces
  return true;
}



template <int dim, int spacedim>
std::pair<Table<2,bool>, std::vector<unsigned int> >
FE_DGP<dim,spacedim>::get_constant_modes () const
{
  Table<2,bool> constant_modes(1, this->dofs_per_cell);
  constant_modes(0,0) = true;
  return std::pair<Table<2,bool>, std::vector<unsigned int> >
         (constant_modes, std::vector<unsigned int>(1, 0));
}



template <int dim, int spacedim>
std::size_t
FE_DGP<dim,spacedim>::memory_consumption () const
{
  Assert (false, ExcNotImplemented ());
  return 0;
}



// explicit instantiations
#include "fe_dgp.inst"


DEAL_II_NAMESPACE_CLOSE
