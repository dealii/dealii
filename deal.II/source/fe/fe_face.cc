//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2009, 2011 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------

#include <deal.II/fe/fe_face.h>
#include <deal.II/fe/fe_poly_face.templates.h>
#include <sstream>

DEAL_II_NAMESPACE_OPEN

template <int dim, int spacedim>
FE_FaceQ<dim,spacedim>::FE_FaceQ (const unsigned int degree)
		:
                FE_PolyFace<TensorProductPolynomials<dim-1>, dim, spacedim> (
		  TensorProductPolynomials<dim-1>(Polynomials::LagrangeEquidistant::generate_complete_basis(degree)),
		  FiniteElementData<dim>(get_dpo_vector(degree), 1, degree, FiniteElementData<dim>::L2),
		  std::vector<bool>(1,true))
{}


template <int dim, int spacedim>
FiniteElement<dim,spacedim>*
FE_FaceQ<dim,spacedim>::clone() const
{
  return new FE_FaceQ<dim,spacedim>(this->degree);
}


template <int dim, int spacedim>
std::string
FE_FaceQ<dim,spacedim>::get_name () const
{
				   // note that the
				   // FETools::get_fe_from_name
				   // function depends on the
				   // particular format of the string
				   // this function returns, so they
				   // have to be kept in synch

  std::ostringstream namebuf;
  namebuf << "FE_FaceQ<" << dim << ">(" << this->degree << ")";

  return namebuf.str();
}



template <int dim, int spacedim>
bool
FE_FaceQ<dim,spacedim>::has_support_on_face (
  const unsigned int shape_index,
  const unsigned int face_index) const
{
  return (face_index == (shape_index/this->dofs_per_face));
}


template <int dim, int spacedim>
std::vector<unsigned int>
FE_FaceQ<dim,spacedim>::get_dpo_vector (const unsigned int deg)
{
  std::vector<unsigned int> dpo(dim+1, 0U);
  dpo[dim-1] = deg+1;
  for (unsigned int i=1;i<dim-1;++i)
    dpo[dim-1] *= deg+1;
  return dpo;
}



// explicit instantiations
#include "fe_face.inst"


DEAL_II_NAMESPACE_CLOSE
