//----------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------------------------------------------

#include <base/quadrature.h>
#include <base/polynomial.h>
#include <base/tensor_product_polynomials.h>
#include <grid/tria.h>
#include <grid/tria_iterator.h>
#include <dofs/dof_accessor.h>
#include <fe/fe.h>
#include <fe/mapping.h>
#include <fe/mapping_q1.h>
#include <fe/fe_dgq.h>
#include <fe/fe_values.h>

static std::vector<bool> dummy(1,true);

// Embedding matrices (produced by tests/fe/embedding with postprocessing)

#if (deal_II_dimension == 1)
#include "mat_dgq.1"
#endif

#if (deal_II_dimension == 2)
#include "mat_dgq.2"
#endif

#if (deal_II_dimension == 3)
#include "mat_dgq.3"
#endif


template <int dim>
FE_DGQ<dim>::InternalData::~InternalData ()
{
  for (unsigned int i=0;i<differences.size ();++i)
    if (differences[i] != 0)
      delete differences[i];
}




template <int dim>
FE_DGQ<dim>::FE_DGQ (unsigned int degree)
		:
		FiniteElement<dim> (FiniteElementData<dim>(get_dpo_vector(degree),1), dummy),
				      degree(degree),
				      polynomials(degree+1),
				      poly(0)
{
  std::vector<SmartPointer<Polynomial<double> > > v(degree+1);
  if (degree==0)
    {
      std::vector<double> coeff(1);
      coeff[0] = 1.;
      polynomials[0] = Polynomial<double> (coeff);
      v[0] = &(polynomials[0]);
    } else {
      for (unsigned int i=0;i<=degree;++i)
	{
	  LagrangeEquidistant p(degree, i);
	  polynomials[i] = p;
	  v[i] = &(polynomials[i]);
	}
    }
  
  poly = new TensorProductPolynomials<dim> (v);

  Assert (degree <= 10, ExcNotImplemented());
  
  std::vector<unsigned int> right;
  std::vector<unsigned int> top;
  rotate_indices (right, 'Z');
  if (dim>2)
    rotate_indices (top, 'X');

  if (embedding[degree] != 0)
    {
      prolongation[0].fill (embedding[degree]);
      switch (dim)
	{
	  case 1:
		prolongation[1].fill_permutation (prolongation[0],
						  right, right);
		break;
	  case 2:
		prolongation[1].fill_permutation (prolongation[0],
						  right, right);
		prolongation[2].fill_permutation (prolongation[1],
						  right, right);
		prolongation[3].fill_permutation (prolongation[2],
						  right, right);
		break;
	  case 3:
		prolongation[1].fill_permutation (prolongation[0],
						  right, right);
		prolongation[5].fill_permutation (prolongation[1],
						  right, right);
		prolongation[4].fill_permutation (prolongation[5],
						  right, right);
		prolongation[7].fill_permutation (prolongation[4],
						  top, top);
		prolongation[3].fill_permutation (prolongation[7],
						  top, top);
		prolongation[6].fill_permutation (prolongation[5],
						  top, top);
		prolongation[2].fill_permutation (prolongation[6],
						  top, top);
		break;
	  default:
		for (unsigned int i=0;i<GeometryInfo<dim>::children_per_cell;++i)
		  prolongation[i].reinit(0);
	}
    } else {
      for (unsigned int i=0;i<GeometryInfo<dim>::children_per_cell;++i)
	prolongation[i].reinit(0);
    }
  
  if (projection_matrices[degree] != 0)
    {
      restriction[0].fill (projection_matrices[degree]);
      switch (dim)
	{
	  case 1:
		restriction[1].fill_permutation (restriction[0],
						 right, right);
		break;
	  case 2:
		restriction[1].fill_permutation (restriction[0],
						 right, right);
		restriction[2].fill_permutation (restriction[1],
						 right, right);
		restriction[3].fill_permutation (restriction[2],
						 right, right);
		break;
	  case 3:
		restriction[1].fill_permutation (restriction[0],
						 right, right);
		restriction[5].fill_permutation (restriction[1],
						 right, right);
		restriction[4].fill_permutation (restriction[5],
						 right, right);
		restriction[7].fill_permutation (restriction[4],
						 top, top);
		restriction[3].fill_permutation (restriction[7],
						 top, top);
		restriction[6].fill_permutation (restriction[5],
						 top, top);
		restriction[2].fill_permutation (restriction[6],
						 top, top);
		break;
	  default:
		for (unsigned int i=0;i<GeometryInfo<dim>::children_per_cell;++i)
		  restriction[i].reinit(0);
	}
    }
}


template <int dim>
FE_DGQ<dim>::~FE_DGQ ()
{
  delete poly;
}


template <int dim>
FiniteElement<dim> *
FE_DGQ<dim>::clone() const
{
  return new FE_DGQ<dim>(degree);
}


template <int dim>
void
FE_DGQ<dim>::get_unit_support_points (std::vector<Point<dim> > &points) const
{
  if (degree>0)
    compute_support_points (points, degree);
  else
    {
      points.resize(1);
      for (unsigned int i=0; i<dim; ++i)
	points[0](i)=0.5;
    }
}

    
template <int dim>
void
FE_DGQ<dim>::get_unit_face_support_points (std::vector<Point<dim-1> > &points) const
{
  points.resize(0);
}

//----------------------------------------------------------------------
// Auxilliary functions
//----------------------------------------------------------------------


template <int dim>
std::vector<unsigned int>
FE_DGQ<dim>::get_dpo_vector(unsigned int deg)
{
  std::vector<unsigned int> dpo(dim+1, 0);
  dpo[dim] = ++deg;
  for (unsigned int i=1;i<dim;++i)
    dpo[dim] *= deg;
  return dpo;
}


template <int dim>
UpdateFlags
FE_DGQ<dim>::update_once (UpdateFlags flags) const
{
  UpdateFlags out = update_default;

//  cerr << "FE_DGQ:Once: " << hex << flags << ' ';

  if (flags & update_values)
    out |= update_values;

//  cerr << out << dec << endl;

  return out;
}


template <int dim>
UpdateFlags
FE_DGQ<dim>::update_each (UpdateFlags flags) const
{
  UpdateFlags out = update_default;

//  cerr << "FE_DGQ:Each: " << hex << flags << ' ';

  if (flags & update_gradients)
    out |= update_gradients | update_covariant_transformation;
  if (flags & update_second_derivatives)
    out |= update_second_derivatives;

//  cerr << out << dec << endl;

  return out;
}


template <int dim>
void
FE_DGQ<dim>::compute_support_points (std::vector<Point<dim> >& support_points,
				     unsigned int degree)
{
  Assert(degree>0, ExcInternalError());
				   // number of points: (degree+1)^dim
  unsigned int n= degree+1;
  for (unsigned int i=1;i<dim;++i)
    n *= degree+1;

  support_points.resize(n);
  
  double step = 1./degree;
  Point<dim> p;
  
  unsigned int k=0;
  for (unsigned int iz=0;iz <= ((dim>2) ? degree : 0) ; ++iz)
    for (unsigned int iy=0;iy <= ((dim>1) ? degree : 0) ; ++iy)
      for (unsigned int ix=0;ix<=degree;++ix)
        {
	  p(0) = ix * step;
	  if (dim>1)
	    p(1) = iy * step;
	  if (dim>2)
	    p(2) = iz * step;
	  
	  support_points[k++] = p;
	}
}


template <int dim>
void
FE_DGQ<dim>::rotate_indices (std::vector<unsigned int>& numbers,
			     const char direction) const
{
  const unsigned int n = degree+1;
  unsigned int s = n;
  for (unsigned int i=1;i<dim;++i)
    s *= n;
  numbers.resize (s);
  
  unsigned int l = 0;

  if (dim==1)
    {
				       // Mirror around midpoint
      for (unsigned int i=n;i>0;)
	numbers[l++]=--i;
    } else {
      switch (direction)
	{
					   // Rotate xy-plane
					   // counter-clockwise
	  case 'z':
		for (unsigned int iz=0;iz<((dim>2) ? n:1);++iz)
		  for (unsigned int j=0;j<n;++j)
		    for (unsigned int i=0;i<n;++i)
		      {
			unsigned int k = n*i-j+n-1 + n*n*iz;
			numbers[l++] = k;
		      }
		break;
						 // Rotate xy-plane
						 // clockwise
	  case 'Z':
		for (unsigned int iz=0;iz<((dim>2) ? n:1);++iz)
		  for (unsigned int iy=0;iy<n;++iy)
		    for (unsigned int ix=0;ix<n;++ix)
		      {
			unsigned int k = n*ix-iy+n-1 + n*n*iz;
			numbers[k] = l++;
		      }
		break;
						 // Rotate yz-plane
						 // counter-clockwise
	  case 'x':
		Assert (dim>2, ExcDimensionMismatch (dim,3));
		for (unsigned int iz=0;iz<n;++iz)
		  for (unsigned int iy=0;iy<n;++iy)
		    for (unsigned int ix=0;ix<n;++ix)
		      {
			unsigned int k = n*(n*iy-iz+n-1) + ix;
			numbers[l++] = k;
		      }
		break;
						 // Rotate yz-plane
						 // clockwise
	  case 'X':
		Assert (dim>2, ExcDimensionMismatch (dim,3));
		for (unsigned int iz=0;iz<n;++iz)
		  for (unsigned int iy=0;iy<n;++iy)
		    for (unsigned int ix=0;ix<n;++ix)
		      {
			unsigned int k = n*(n*iy-iz+n-1) + ix;
			numbers[k] = l++;
		      }
		break;
	  default:
		Assert (false, ExcNotImplemented ());
	}
    }
}




//----------------------------------------------------------------------
// Data field initialization
//----------------------------------------------------------------------
template <int dim>
Mapping<dim>::InternalDataBase*
FE_DGQ<dim>::get_data (const UpdateFlags update_flags,
		       const Mapping<dim>& mapping,
		       const Quadrature<dim> &quadrature) const
{
  InternalData* data = new InternalData;
  std::vector<double> values(0);
  std::vector<Tensor<1,dim> > grads(0);
  std::vector<Tensor<2,dim> > grad_grads(0);
  
  data->update_once = update_once(update_flags);
  data->update_each = update_each(update_flags);
  data->update_flags = data->update_once | data->update_each;

  const UpdateFlags flags(data->update_flags);
  
  if (flags & update_values)
    {
      values.resize (dofs_per_cell);
      data->shape_values.resize(dofs_per_cell,
				std::vector<double>(quadrature.n_quadrature_points));
    }

  if (flags & update_gradients)
    {
      grads.resize (dofs_per_cell);
      data->shape_gradients.resize(dofs_per_cell,
				   std::vector<Tensor<1,dim> >(quadrature.n_quadrature_points));
    }

  if (flags & update_second_derivatives)
    {
      data->initialize (this, mapping, quadrature);
    }
  
  
  if (flags & (update_values | update_gradients))
    for (unsigned int i=0;i<quadrature.n_quadrature_points;++i)
      {
	poly->compute(quadrature.point(i), values, grads, grad_grads);
	for (unsigned int k=0;k<dofs_per_cell; ++k)
	  {
	    if (flags & update_values)
	      data->shape_values[k][i] = values[k];
	    if (flags & update_gradients)
	      data->shape_gradients[k][i] = grads[k];
	  }
      }
  return data;
}



//----------------------------------------------------------------------
// Fill data of FEValues
//----------------------------------------------------------------------
template <int dim>
void
FE_DGQ<dim>::fill_fe_values (const Mapping<dim>                   &mapping,
			     const DoFHandler<dim>::cell_iterator &cell,
			     const Quadrature<dim>                &quadrature,
			     Mapping<dim>::InternalDataBase      &mapping_data,
			     Mapping<dim>::InternalDataBase      &fedata,
			     FEValuesData<dim>                    &data) const
{
  InternalData &fe_data = dynamic_cast<InternalData &> (fedata);
  
  const UpdateFlags flags(fe_data.current_update_flags());

  for (unsigned int k=0;k<dofs_per_cell;++k)
    {
      for (unsigned int i=0;i<quadrature.n_quadrature_points;++i)
	{
	  if (flags & update_values)
	    data.shape_values(k,i)
	      = fe_data.shape_values[k][i];
	}
      
      if (flags & update_gradients)
	mapping.transform_covariant(data.shape_gradients[k],
				    fe_data.shape_gradients[k],
				    mapping_data, 0);
    }
  if (flags & update_second_derivatives)
    compute_2nd (mapping, cell, 0, mapping_data, fe_data, data);
  
  fe_data.first_cell = false;
}



template <int dim>
void
FE_DGQ<dim>::fill_fe_face_values (const Mapping<dim> &mapping,
				  const DoFHandler<dim>::cell_iterator &cell,
				  const unsigned int face,
				  const Quadrature<dim-1>& quadrature,
				  Mapping<dim>::InternalDataBase& mapping_data,
				  Mapping<dim>::InternalDataBase& fedata,
				  FEValuesData<dim>& data) const
{
  InternalData &fe_data = dynamic_cast<InternalData &> (fedata);
  unsigned int offset = face * quadrature.n_quadrature_points;
  
  const UpdateFlags flags(fe_data.update_once | fe_data.update_each);

  for (unsigned int k=0;k<dofs_per_cell;++k)
    {
      for (unsigned int i=0;i<quadrature.n_quadrature_points;++i)
	{
	  if (flags & update_values)
	    data.shape_values(k,i)
	      = fe_data.shape_values[k][i+offset];
	}
      
      if (flags & update_gradients)
	mapping.transform_covariant(data.shape_gradients[k],
				    fe_data.shape_gradients[k],
				    mapping_data, offset);
      if (flags & update_second_derivatives)
	compute_2nd (mapping, cell, offset, mapping_data, fe_data, data);
    }
  
  fe_data.first_cell = false;
}



template <int dim>
void
FE_DGQ<dim>::fill_fe_subface_values (const Mapping<dim> &mapping,
				     const DoFHandler<dim>::cell_iterator &cell,
				     const unsigned int face,
				     const unsigned int subface,
				     const Quadrature<dim-1>& quadrature,
				     Mapping<dim>::InternalDataBase& mapping_data,
				     Mapping<dim>::InternalDataBase& fedata,
				     FEValuesData<dim>& data) const
{
  InternalData &fe_data = dynamic_cast<InternalData &> (fedata);
  unsigned int offset = (face * GeometryInfo<dim>::subfaces_per_face
			 + subface) * quadrature.n_quadrature_points;

  const UpdateFlags flags(fe_data.update_once | fe_data.update_each);

  for (unsigned int k=0;k<dofs_per_cell;++k)
    {
      for (unsigned int i=0;i<quadrature.n_quadrature_points;++i)
	{
	  if (flags & update_values)
	    data.shape_values(k,i)
	      = fe_data.shape_values[k][i+offset];
	}
      
      if (flags & update_gradients)
	mapping.transform_covariant(data.shape_gradients[k],
				    fe_data.shape_gradients[k],
				    mapping_data, offset);
      if (flags & update_second_derivatives)
	compute_2nd (mapping, cell, offset, mapping_data, fe_data, data);
    }
  
  fe_data.first_cell = false;
}


template <int dim>
unsigned int
FE_DGQ<dim>::memory_consumption () const
{
  Assert (false, ExcNotImplemented ());
  return 0;
}



template FE_DGQ<deal_II_dimension>;
