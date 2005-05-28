//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2003, 2004, 2005 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------

#include <base/quadrature_lib.h>
#include <base/table.h>
#include <grid/tria.h>
#include <grid/tria_iterator.h>
#include <dofs/dof_accessor.h>
#include <fe/fe.h>
#include <fe/mapping.h>
#include <fe/fe_raviart_thomas.h>
#include <fe/fe_values.h>
#include <fe/fe_tools.h>

#ifdef HAVE_STD_STRINGSTREAM
#  include <sstream>
#else
#  include <strstream>
#endif

#include <iostream>

template <int dim>
FE_RaviartThomasNodal<dim>::FE_RaviartThomasNodal (const unsigned int deg)
		:
		FE_PolyTensor<PolynomialsRaviartThomas<dim>, dim> (
		  deg,
		  FiniteElementData<dim>(get_dpo_vector(deg),
					 dim, deg+1),
		  get_ria_vector (deg),
		  std::vector<std::vector<bool> >(
		    FiniteElementData<dim>(get_dpo_vector(deg),
					   dim,deg+1).dofs_per_cell,
		    std::vector<bool>(dim,true)))
{
  Assert (dim >= 2, ExcImpossibleInDim(dim));
  
  this->mapping_type = this->independent_on_cartesian;
				   // These must be done first, since
				   // they change the evaluation of
				   // basis functions
  initialize_unit_support_points (deg);
  initialize_node_matrix();

  for (unsigned int i=0; i<GeometryInfo<dim>::children_per_cell; ++i)
    this->prolongation[i].reinit (this->dofs_per_cell,
				  this->dofs_per_cell);
  FETools::compute_embedding_matrices (*this, &this->prolongation[0]);
}



template <int dim>
std::string
FE_RaviartThomasNodal<dim>::get_name () const
{
				   // note that the
				   // FETools::get_fe_from_name
				   // function depends on the
				   // particular format of the string
				   // this function returns, so they
				   // have to be kept in synch

#ifdef HAVE_STD_STRINGSTREAM
  std::ostringstream namebuf;
#else
  std::ostrstream namebuf;
#endif
  
  namebuf << "FE_RaviartThomasNodal<" << dim << ">(" << this->degree-1 << ")";

#ifndef HAVE_STD_STRINGSTREAM
  namebuf << std::ends;
#endif
  return namebuf.str();
}



template <int dim>
bool
FE_RaviartThomasNodal<dim>::has_support_on_face (unsigned int, unsigned int) const
{
  return true;
}



template <int dim>
FiniteElement<dim> *
FE_RaviartThomasNodal<dim>::clone() const
{
  return new FE_RaviartThomasNodal<dim>(this->degree-1);
}


#if deal_II_dimension == 1

template <>
std::vector<unsigned int>
FE_RaviartThomasNodal<1>::get_dpo_vector (const unsigned int deg)
{
  std::vector<unsigned int> dpo(2);
  dpo[0] = 1;
  dpo[1] = deg;
  return dpo;
}

#endif


template <int dim>
std::vector<unsigned int>
FE_RaviartThomasNodal<dim>::get_dpo_vector (const unsigned int deg)
{
                                   // the element is face-based and we have
                                   // (deg+1)^(dim-1) DoFs per face
  unsigned int dofs_per_face = 1;
  for (unsigned int d=1; d<dim; ++d)
    dofs_per_face *= deg+1;

                                   // and then there are interior dofs
  const unsigned int
    interior_dofs = dim*deg*dofs_per_face;
  
  std::vector<unsigned int> dpo(dim+1);
  dpo[dim-1] = dofs_per_face;
  dpo[dim]   = interior_dofs;
  
  return dpo;
}



#if deal_II_dimension == 1

template <>
std::vector<bool>
FE_RaviartThomasNodal<1>::get_ria_vector (const unsigned int)
{
  Assert (false, ExcImpossibleInDim(1));
  return std::vector<bool>();
}

#endif


template <int dim>
std::vector<bool>
FE_RaviartThomasNodal<dim>::get_ria_vector (const unsigned int deg)
{
  unsigned int dofs_per_cell, dofs_per_face;
  switch (dim)
    {
      case 2:
	    dofs_per_face = deg+1;
	    dofs_per_cell = 2*(deg+1)*(deg+2);
	    break;
      case 3:
	    dofs_per_face = (deg+1)*(deg+1);
	    dofs_per_cell = 3*(deg+1)*(deg+1)*(deg+2);
	    break;
      default:
	    Assert (false, ExcNotImplemented());
    }
  Assert (FiniteElementData<dim>(get_dpo_vector(deg),dim).dofs_per_cell ==
	  dofs_per_cell,
	  ExcInternalError());
  Assert (FiniteElementData<dim>(get_dpo_vector(deg),dim).dofs_per_face ==
	  dofs_per_face,
	  ExcInternalError());
  
				   // all face dofs need to be
				   // non-additive, since they have
				   // continuity requirements.
				   // however, the interior dofs are
				   // made additive
  std::vector<bool> ret_val(dofs_per_cell,false);
  for (unsigned int i=GeometryInfo<dim>::faces_per_cell*dofs_per_face;
       i < dofs_per_cell; ++i)
    ret_val[i] = true;

  return ret_val;
}


template <int dim>
void
FE_RaviartThomasNodal<dim>::initialize_unit_support_points (const unsigned int deg)
{
  this->unit_support_points.resize (this->dofs_per_cell);
  this->unit_face_support_points.resize (this->dofs_per_face);
  
  unsigned int current = 0;

				   // On the faces, we choose as many
				   // Gauss points as necessary to
				   // determine the normal component
				   // uniquely. This is the deg of
				   // the Raviart-Thomas element plus
				   // one.
  if (dim>1)
    {
      QGauss<dim-1> face_points (deg+1);
      Assert (face_points.n_quadrature_points == this->dofs_per_face,
	      ExcInternalError());
      for (unsigned int k=0;k<this->dofs_per_face;++k)
	this->unit_face_support_points[k] = face_points.point(k);
      Quadrature<dim> faces = QProjector<dim>::project_to_all_faces(face_points);
      for (unsigned int k=0;k<//faces.n_quadrature_points
			  this->dofs_per_face*GeometryInfo<dim>::faces_per_cell;++k)
	this->unit_support_points[k] = faces.point(k);

      current = this->dofs_per_face*GeometryInfo<dim>::faces_per_cell;
    }
  
  if (deg==0) return;
				   // In the interior, we need
				   // anisotropic Gauss quadratures,
				   // different for each direction.
  QGauss<1> high(deg+1);
  QGauss<1> low(deg);

  for (unsigned int d=0;d<dim;++d)
    {
      QAnisotropic<dim>* quadrature;
      if (dim == 1) quadrature = new QAnisotropic<dim>(high);
      if (dim == 2) quadrature = new QAnisotropic<dim>(((d==0) ? low : high),
						       ((d==1) ? low : high));
      if (dim == 3) quadrature = new QAnisotropic<dim>(((d==0) ? low : high),
						       ((d==1) ? low : high),
						       ((d==2) ? low : high));
      for (unsigned int k=0;k<quadrature->n_quadrature_points;++k)
	this->unit_support_points[current++] = quadrature->point(k);
      delete quadrature;
    }
  Assert (current == this->dofs_per_cell, ExcInternalError());
}


template <int dim>
void
FE_RaviartThomasNodal<dim>::initialize_node_matrix ()
{
  const unsigned int n_dofs = this->dofs_per_cell;

				   // We use an auxiliary matrix in
				   // this function. Therefore,
				   // inverse_node_matrix is still
				   // empty and shape_value_component
				   // returns the 'raw' shape values.
  FullMatrix<double> N(n_dofs, n_dofs);

				   // The curent node functional index
  unsigned int current = 0;

				   // For each face and all quadrature
				   // points on it, the node value is
				   // the normal component of the
				   // shape function, possibly
				   // pointing in negative direction.
  for (unsigned int face = 0; face<GeometryInfo<dim>::faces_per_cell; ++face)
    for (unsigned int k=0;k<this->dofs_per_face;++k)
      {
	for (unsigned int i=0;i<n_dofs;++i)
	  N(current,i) = this->shape_value_component(
	    i, this->unit_support_points[current],
	    GeometryInfo< dim >::unit_normal_direction[face]);
//			 * GeometryInfo< dim >::unit_normal_orientation[face];
	++current;
      }
				   // Interior degrees of freedom in each direction
  const unsigned int n_cell = (n_dofs - current) / dim;
  
  for (unsigned int d=0;d<dim;++d)
    for (unsigned int k=0;k<n_cell;++k)
      {
	for (unsigned int i=0;i<n_dofs;++i)
	  N(current,i) = this->shape_value_component(i, this->unit_support_points[current], d);
	++current;
      }
  Assert (current == n_dofs, ExcInternalError());
  this->inverse_node_matrix.reinit(n_dofs, n_dofs);
  this->inverse_node_matrix.invert(N);
}



template class FE_RaviartThomasNodal<deal_II_dimension>;
