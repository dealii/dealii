//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2002, 2003, 2004, 2005, 2006, 2007, 2008 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------

#include <fe/fe_dgp.h>
#include <fe/fe_tools.h>

#include <sstream>

DEAL_II_NAMESPACE_OPEN

template <int dim>
FE_DGP<dim>::FE_DGP (const unsigned int degree)
		:
		FE_Poly<PolynomialSpace<dim>, dim> (
		  PolynomialSpace<dim>(Polynomials::Legendre::generate_complete_basis(degree)),
		  FiniteElementData<dim>(get_dpo_vector(degree), 1, degree, FiniteElementData<dim>::L2),
		  std::vector<bool>(FiniteElementData<dim>(get_dpo_vector(degree), 1, degree).dofs_per_cell,true),
		  std::vector<std::vector<bool> >(FiniteElementData<dim>(
		    get_dpo_vector(degree), 1, degree).dofs_per_cell, std::vector<bool>(1,true)))
{
				   // Reinit the vectors of
				   // restriction and prolongation
				   // matrices to the right sizes
  this->reinit_restriction_and_prolongation_matrices();
				   // Fill prolongation matrices with embedding operators
  FETools::compute_embedding_matrices (*this, this->prolongation);
				   // Fill restriction matrices with L2-projection
  FETools::compute_projection_matrices (*this, this->restriction);
}



template <int dim>
std::string
FE_DGP<dim>::get_name () const
{
				   // note that the
				   // FETools::get_fe_from_name
				   // function depends on the
				   // particular format of the string
				   // this function returns, so they
				   // have to be kept in synch

  std::ostringstream namebuf;  
  namebuf << "FE_DGP<" << dim << ">(" << this->degree << ")";

  return namebuf.str();
}



template <int dim>
FiniteElement<dim> *
FE_DGP<dim>::clone() const
{
  return new FE_DGP<dim>(*this);
}



//---------------------------------------------------------------------------
// Auxiliary functions
//---------------------------------------------------------------------------


template <int dim>
std::vector<unsigned int>
FE_DGP<dim>::get_dpo_vector (const unsigned int deg)
{
  std::vector<unsigned int> dpo(dim+1, 0U);
  dpo[dim] = deg+1;
  for (unsigned int i=1;i<dim;++i)
    {
      dpo[dim] *= deg+1+i;
      dpo[dim] /= i+1;
    }
  return dpo;
}



template <int dim>
void
FE_DGP<dim>::
get_face_interpolation_matrix (const FiniteElement<dim> &x_source_fe,
			       FullMatrix<double>       &interpolation_matrix) const
{
				   // this is only implemented, if the source
				   // FE is also a DGP element. in that case,
				   // both elements have no dofs on their
				   // faces and the face interpolation matrix
				   // is necessarily empty -- i.e. there isn't
				   // much we need to do here.
  AssertThrow ((x_source_fe.get_name().find ("FE_DGP<") == 0)
               ||
               (dynamic_cast<const FE_DGP<dim>*>(&x_source_fe) != 0),
               typename FiniteElement<dim>::
               ExcInterpolationNotImplemented());
  
  Assert (interpolation_matrix.m() == 0,
	  ExcDimensionMismatch (interpolation_matrix.m(),
				0));
  Assert (interpolation_matrix.n() == 0,
	  ExcDimensionMismatch (interpolation_matrix.n(),
				0));
}



template <int dim>
void
FE_DGP<dim>::
get_subface_interpolation_matrix (const FiniteElement<dim> &x_source_fe,
				  const unsigned int ,
				  FullMatrix<double>           &interpolation_matrix) const
{
				   // this is only implemented, if the source
				   // FE is also a DGP element. in that case,
				   // both elements have no dofs on their
				   // faces and the face interpolation matrix
				   // is necessarily empty -- i.e. there isn't
				   // much we need to do here.
  AssertThrow ((x_source_fe.get_name().find ("FE_DGP<") == 0)
               ||
               (dynamic_cast<const FE_DGP<dim>*>(&x_source_fe) != 0),
               typename FiniteElement<dim>::
               ExcInterpolationNotImplemented());
  
  Assert (interpolation_matrix.m() == 0,
	  ExcDimensionMismatch (interpolation_matrix.m(),
				0));
  Assert (interpolation_matrix.n() == 0,
	  ExcDimensionMismatch (interpolation_matrix.n(),
				0));
}



template <int dim>
bool
FE_DGP<dim>::hp_constraints_are_implemented () const
{
  return true;
}



template <int dim>
std::vector<std::pair<unsigned int, unsigned int> >
FE_DGP<dim>::
hp_vertex_dof_identities (const FiniteElement<dim> &fe_other) const
{
				   // there are no such constraints for DGP
				   // elements at all
  if (dynamic_cast<const FE_DGP<dim>*>(&fe_other) != 0)
    return
      std::vector<std::pair<unsigned int, unsigned int> > ();
  else
    {
      Assert (false, ExcNotImplemented());
      return std::vector<std::pair<unsigned int, unsigned int> > ();
    }
}



template <int dim>
std::vector<std::pair<unsigned int, unsigned int> >
FE_DGP<dim>::
hp_line_dof_identities (const FiniteElement<dim> &fe_other) const
{
				   // there are no such constraints for DGP
				   // elements at all
  if (dynamic_cast<const FE_DGP<dim>*>(&fe_other) != 0)
    return
      std::vector<std::pair<unsigned int, unsigned int> > ();
  else
    {
      Assert (false, ExcNotImplemented());
      return std::vector<std::pair<unsigned int, unsigned int> > ();
    }
}



template <int dim>
std::vector<std::pair<unsigned int, unsigned int> >
FE_DGP<dim>::
hp_quad_dof_identities (const FiniteElement<dim>        &fe_other) const
{
				   // there are no such constraints for DGP
				   // elements at all
  if (dynamic_cast<const FE_DGP<dim>*>(&fe_other) != 0)
    return
      std::vector<std::pair<unsigned int, unsigned int> > ();
  else
    {
      Assert (false, ExcNotImplemented());
      return std::vector<std::pair<unsigned int, unsigned int> > ();
    }
}



template <int dim>
FiniteElementDomination::Domination
FE_DGP<dim>::compare_for_face_domination (const FiniteElement<dim> &fe_other) const
{
				   // check whether both are discontinuous
				   // elements and both could dominate, see
				   // the description of
				   // FiniteElementDomination::Domination
  if (dynamic_cast<const FE_DGP<dim>*>(&fe_other) != 0)
    return FiniteElementDomination::either_element_can_dominate;

  Assert (false, ExcNotImplemented());
  return FiniteElementDomination::neither_element_dominates;
}



template <int dim>
bool
FE_DGP<dim>::has_support_on_face (const unsigned int,
				  const unsigned int) const
{
                                   // all shape functions have support on all
                                   // faces
  return true;
}



template <int dim>
unsigned int
FE_DGP<dim>::memory_consumption () const
{
  Assert (false, ExcNotImplemented ());
  return 0;
}



template class FE_DGP<deal_II_dimension>;

DEAL_II_NAMESPACE_CLOSE
