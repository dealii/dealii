//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2003, 2004, 2005, 2006, 2007, 2008 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------

#include <base/quadrature_lib.h>
#include <base/qprojector.h>
#include <base/table.h>
#include <grid/tria.h>
#include <grid/tria_iterator.h>
#include <dofs/dof_accessor.h>
#include <fe/fe.h>
#include <fe/mapping.h>
#include <fe/fe_raviart_thomas.h>
#include <fe/fe_values.h>
#include <fe/fe_tools.h>

#include <sstream>

//TODO: implement the adjust_quad_dof_index_for_face_orientation_table and
//adjust_line_dof_index_for_line_orientation_table fields, and write tests
//similar to bits/face_orientation_and_fe_q_*


DEAL_II_NAMESPACE_OPEN


template <int dim>
FE_RaviartThomasNodal<dim>::FE_RaviartThomasNodal (const unsigned int deg)
		:
		FE_PolyTensor<PolynomialsRaviartThomas<dim>, dim> (
		  deg,
		  FiniteElementData<dim>(get_dpo_vector(deg),
					 dim, deg+1, FiniteElementData<dim>::Hdiv, 1),
		  get_ria_vector (deg),
		  std::vector<std::vector<bool> >(PolynomialsRaviartThomas<dim>::compute_n_pols(deg),
						  std::vector<bool>(dim,true)))
{
  Assert (dim >= 2, ExcImpossibleInDim(dim));
  const unsigned int n_dofs = this->dofs_per_cell;
  
  this->mapping_type = this->independent_on_cartesian;
				   // These must be done first, since
				   // they change the evaluation of
				   // basis functions

				   // Set up the generalized support
				   // points
  initialize_support_points (deg);
				   //Now compute the inverse node
				   //matrix, generating the correct
				   //basis functions from the raw
				   //ones.
  
				   // We use an auxiliary matrix in
				   // this function. Therefore,
				   // inverse_node_matrix is still
				   // empty and shape_value_component
				   // returns the 'raw' shape values.
  FullMatrix<double> M(n_dofs, n_dofs);
  FETools::compute_node_matrix(M, *this);
  this->inverse_node_matrix.reinit(n_dofs, n_dofs);
  this->inverse_node_matrix.invert(M);
				   // From now on, the shape functions
				   // will be the correct ones, not
				   // the raw shape functions anymore.
  
  for (unsigned int i=0; i<GeometryInfo<dim>::children_per_cell; ++i)
    this->prolongation[i].reinit (this->dofs_per_cell,
				  this->dofs_per_cell);
  FETools::compute_embedding_matrices (*this, this->prolongation);
  
  FullMatrix<double> face_embeddings[GeometryInfo<dim>::subfaces_per_face];
  for (unsigned int i=0; i<GeometryInfo<dim>::subfaces_per_face; ++i)
    face_embeddings[i].reinit (this->dofs_per_face, this->dofs_per_face);
  FETools::compute_face_embedding_matrices(*this, face_embeddings, 0, 0);
  this->interface_constraints.reinit((1<<(dim-1)) * this->dofs_per_face,
				     this->dofs_per_face);
  unsigned int target_row=0;
  for (unsigned int d=0;d<GeometryInfo<dim>::subfaces_per_face;++d)
    for (unsigned int i=0;i<face_embeddings[d].m();++i)
      {
	for (unsigned int j=0;j<face_embeddings[d].n();++j)
	  this->interface_constraints(target_row,j) = face_embeddings[d](i,j);
	++target_row;
      }
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

  std::ostringstream namebuf;  
  namebuf << "FE_RaviartThomasNodal<" << dim << ">(" << this->degree-1 << ")";

  return namebuf.str();
}


template <int dim>
FiniteElement<dim> *
FE_RaviartThomasNodal<dim>::clone() const
{
  return new FE_RaviartThomasNodal<dim>(*this);
}



template <int dim>
void
FE_RaviartThomasNodal<dim>::interpolate(
  std::vector<double>&,
  const std::vector<double>&) const
{
  Assert(false, ExcNotImplemented());
}


template <int dim>
void
FE_RaviartThomasNodal<dim>::interpolate(
  std::vector<double>&    local_dofs,
  const std::vector<Vector<double> >& values,
  unsigned int offset) const
{
  Assert (values.size() == this->generalized_support_points.size(),
	  ExcDimensionMismatch(values.size(), this->generalized_support_points.size()));
  Assert (local_dofs.size() == this->dofs_per_cell,
	  ExcDimensionMismatch(local_dofs.size(),this->dofs_per_cell));
  Assert (values[0].size() >= offset+this->n_components(),
	  ExcDimensionMismatch(values[0].size(),offset+this->n_components()));

				   // First do interpolation on
				   // faces. There, the component
				   // evaluated depends on the face
				   // direction and orientation.
  unsigned int fbase = 0;
  unsigned int f=0;
  for (;f<GeometryInfo<dim>::faces_per_cell;
       ++f, fbase+=this->dofs_per_face)
    {
      for (unsigned int i=0;i<this->dofs_per_face;++i)
	{
	  local_dofs[fbase+i] = values[fbase+i](offset+GeometryInfo<dim>::unit_normal_direction[f]);
	}
    }

				   // The remaining points form dim
				   // chunks, one for each component.
  const unsigned int istep = (this->dofs_per_cell - fbase) / dim;
  Assert ((this->dofs_per_cell - fbase) % dim == 0, ExcInternalError());

  f = 0;
  while (fbase < this->dofs_per_cell)
    {
      for (unsigned int i=0;i<istep;++i)
	{
	  local_dofs[fbase+i] = values[fbase+i](offset+f);
	}
      fbase+=istep;
      ++f;
    }
  Assert (fbase == this->dofs_per_cell, ExcInternalError());
}



template <int dim>
void
FE_RaviartThomasNodal<dim>::interpolate(
  std::vector<double>& local_dofs,
  const VectorSlice<const std::vector<std::vector<double> > >& values) const
{
  Assert (values.size() == this->n_components(),
	  ExcDimensionMismatch(values.size(), this->n_components()));
  Assert (values[0].size() == this->generalized_support_points.size(),
	  ExcDimensionMismatch(values.size(), this->generalized_support_points.size()));
  Assert (local_dofs.size() == this->dofs_per_cell,
	  ExcDimensionMismatch(local_dofs.size(),this->dofs_per_cell));
				   // First do interpolation on
				   // faces. There, the component
				   // evaluated depends on the face
				   // direction and orientation.
  unsigned int fbase = 0;
  unsigned int f=0;
  for (;f<GeometryInfo<dim>::faces_per_cell;
       ++f, fbase+=this->dofs_per_face)
    {
      for (unsigned int i=0;i<this->dofs_per_face;++i)
	{
	  local_dofs[fbase+i] = values[GeometryInfo<dim>::unit_normal_direction[f]][fbase+i];
	}
    }
				   // The remaining points form dim
				   // chunks, one for each component.
  const unsigned int istep = (this->dofs_per_cell - fbase) / dim;
  Assert ((this->dofs_per_cell - fbase) % dim == 0, ExcInternalError());

  f = 0;
  while (fbase < this->dofs_per_cell)
    {
      for (unsigned int i=0;i<istep;++i)
	{
	  local_dofs[fbase+i] = values[f][fbase+i];
	}
      fbase+=istep;
      ++f;
    }
  Assert (fbase == this->dofs_per_cell, ExcInternalError());
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
  const unsigned int dofs_per_cell = PolynomialsRaviartThomas<dim>::compute_n_pols(deg);
  unsigned int dofs_per_face = deg+1;
  for (unsigned int d=2;d<dim;++d)
    dofs_per_face *= deg+1;
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
FE_RaviartThomasNodal<dim>::initialize_support_points (const unsigned int deg)
{
  this->generalized_support_points.resize (this->dofs_per_cell);
  this->generalized_face_support_points.resize (this->dofs_per_face);

				   // Number of the point being entered
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
      Assert (face_points.size() == this->dofs_per_face,
	      ExcInternalError());
      for (unsigned int k=0;k<this->dofs_per_face;++k)
	this->generalized_face_support_points[k] = face_points.point(k);
      Quadrature<dim> faces = QProjector<dim>::project_to_all_faces(face_points);
      for (unsigned int k=0;
	   k<this->dofs_per_face*GeometryInfo<dim>::faces_per_cell;
	   ++k)
	this->generalized_support_points[k] = faces.point(k+QProjector<dim>
																	  ::DataSetDescriptor::face(0,
																										 true,
																										 false,
																										 false,
																										 this->dofs_per_face));

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
      Assert(dim<=3, ExcNotImplemented());
      
      for (unsigned int k=0;k<quadrature->size();++k)
	this->generalized_support_points[current++] = quadrature->point(k);
      delete quadrature;
    }
  Assert (current == this->dofs_per_cell, ExcInternalError());
}


template class FE_RaviartThomasNodal<deal_II_dimension>;

DEAL_II_NAMESPACE_CLOSE
