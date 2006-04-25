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

#include <base/quadrature.h>
#include <base/quadrature_lib.h>
#include <base/qprojector.h>
#include <base/table.h>
#include <grid/tria.h>
#include <grid/tria_iterator.h>
#include <dofs/dof_accessor.h>
#include <fe/fe.h>
#include <fe/mapping.h>
#include <fe/fe_abf.h>
#include <fe/fe_values.h>
#include <fe/fe_tools.h>

#ifdef HAVE_STD_STRINGSTREAM
#  include <sstream>
#else
#  include <strstream>
#endif

#include <iostream>
using namespace std;


template <int dim>
FE_ABF<dim>::FE_ABF (const unsigned int deg)
		:
		FE_PolyTensor<PolynomialsABF<dim>, dim> (
		  deg,
		  FiniteElementData<dim>(get_dpo_vector(deg),
					 dim, deg+1, FiniteElementData<dim>::Hdiv, 1),
		  std::vector<bool>(PolynomialsABF<dim>::compute_n_pols(deg), true),
		  std::vector<std::vector<bool> >(PolynomialsABF<dim>::compute_n_pols(deg),
						  std::vector<bool>(dim,true))),
		rt_order(deg)
{
  Assert (dim >= 2, ExcImpossibleInDim(dim));
  const unsigned int n_dofs = this->dofs_per_cell;
  
  this->mapping_type = this->contravariant;
				   // First, initialize the
				   // generalized support points and
				   // quadrature weights, since they
				   // are required for interpolation.
  initialize_support_points(deg);
				   // Now compute the inverse node
				   //matrix, generating the correct
				   //basis functions from the raw
				   //ones.
  FullMatrix<double> M(n_dofs, n_dofs);
  FETools::compute_node_matrix(M, *this);

  M.print (std::cout);

  this->inverse_node_matrix.reinit(n_dofs, n_dofs);
  this->inverse_node_matrix.invert(M);
				   // From now on, the shape functions
				   // will be the correct ones, not
				   // the raw shape functions anymore.
  

				   // initialize the various matrices
  for (unsigned int i=0; i<GeometryInfo<dim>::children_per_cell; ++i)
    {
      this->prolongation[i].reinit (n_dofs, n_dofs);
      this->restriction[i].reinit (n_dofs, n_dofs);
    }

  FETools::compute_embedding_matrices (*this, &this->prolongation[0]);
  //  initialize_restriction ();

  // TODO
  std::vector<FullMatrix<double> >
    face_embeddings(1<<(dim-1), FullMatrix<double>(this->dofs_per_face,
						   this->dofs_per_face));
  //FETools::compute_face_embedding_matrices(*this, &face_embeddings[0], 0, 0);
  this->interface_constraints.reinit((1<<(dim-1)) * this->dofs_per_face,
				     this->dofs_per_face);
  unsigned int target_row=0;
  for (unsigned int d=0;d<face_embeddings.size();++d)
    for (unsigned int i=0;i<face_embeddings[d].m();++i)
      {
	for (unsigned int j=0;j<face_embeddings[d].n();++j)
	  this->interface_constraints(target_row,j) = face_embeddings[d](i,j);
	++target_row;
      }
}



template <int dim>
std::string
FE_ABF<dim>::get_name () const
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
  
  namebuf << "FE_ABF<" << dim << ">(" << rt_order << ")";

#ifndef HAVE_STD_STRINGSTREAM
  namebuf << std::ends;
#endif
  return namebuf.str();
}



template <int dim>
FiniteElement<dim> *
FE_ABF<dim>::clone() const
{
  return new FE_ABF<dim>(rt_order);
}


//---------------------------------------------------------------------------
// Auxiliary and internal functions
//---------------------------------------------------------------------------


#if deal_II_dimension == 1

template <int dim>
void
FE_ABF<dim>::initialize_support_points (const unsigned int deg)
{
  return;
  
  Assert (false, ExcNotImplemented());
  
  QGauss<dim> cell_quadrature(deg+1);
  const unsigned int n_interior_points
    = (deg>0) ? cell_quadrature.n_quadrature_points : 0;
  
  this->generalized_support_points.resize (2 + n_interior_points);
  
				   // Number of the point being entered
  unsigned int current = 0;

  
  if (deg==0) return;

  interior_weights.reinit(TableIndices<3>(2+n_interior_points, 0, dim));
  
  for (unsigned int k=0;k<cell_quadrature.n_quadrature_points;++k)
    this->generalized_support_points[current++] = cell_quadrature.point(k);
  
  Assert (current == this->generalized_support_points.size(),
	  ExcInternalError());
}

#else

// Version for 2d and higher. See above for 1d version
template <int dim>
void
FE_ABF<dim>::initialize_support_points (const unsigned int deg)
{
  QGauss<dim> cell_quadrature(deg+2);
  const unsigned int n_interior_points = cell_quadrature.n_quadrature_points;

  unsigned int n_face_points = (dim>1) ? 1 : 0;
				   // compute (deg+1)^(dim-1)
  for (unsigned int d=1;d<dim;++d)
    n_face_points *= deg+1;

  this->generalized_support_points.resize (GeometryInfo<dim>::faces_per_cell*n_face_points
					   + n_interior_points);
  this->generalized_face_support_points.resize (n_face_points);


  // These might be required when the faces contribution is computed
  // Therefore they will be initialised at this point.
  std::vector<AnisotropicPolynomials<dim>* > polynomials_abf(dim);

  // Generate x_1^{i} x_2^{r+1} ...
  for (unsigned int dd=0; dd<dim; ++dd)
    {
      std::vector<std::vector<Polynomials::Polynomial<double> > > poly(dim);
      for (unsigned int d=0;d<dim;++d)
	poly[d].push_back (Polynomials::Monomial<double> (deg+1));
      poly[dd] = Polynomials::Monomial<double>::generate_complete_basis(deg);

      polynomials_abf[dd] = new AnisotropicPolynomials<dim>(poly);
    }

				   // Number of the point being entered
  unsigned int current = 0;

  if (dim>1)
    {
      QGauss<dim-1> face_points (deg+1);
      TensorProductPolynomials<dim-1> legendre
	= Polynomials::Legendre::generate_complete_basis(deg);

      boundary_weights.reinit(n_face_points, legendre.n());
      
//       Assert (face_points.n_quadrature_points == this->dofs_per_face,
// 	      ExcInternalError());
      
      for (unsigned int k=0;k<n_face_points;++k)
	{
	  this->generalized_face_support_points[k] = face_points.point(k);
					   // Compute its quadrature
					   // contribution for each
					   // moment.
	  for (unsigned int i=0;i<legendre.n();++i)
	    {
	      boundary_weights(k, i)
		= face_points.weight(k)
		* legendre.compute_value(i, face_points.point(k));
	    }
	}

      Quadrature<dim> faces = QProjector<dim>::project_to_all_faces(face_points);
      for (;current<GeometryInfo<dim>::faces_per_cell*n_face_points;
	   ++current)
	{
					   // Enter the support point
					   // into the vector
	  this->generalized_support_points[current] = faces.point(current);
	}


      // Now initialise edge interior weights for the ABF elements.
      // These are completely independent from the usual edge moments. They
      // stem from applying the Gauss theorem to the nodal values, which
      // was necessary to cast the ABF elements into the deal.II framework
      // for vector valued elements.
      boundary_weights_abf.reinit(faces.n_quadrature_points, polynomials_abf[0]->n() * dim);
      for (unsigned int k=0;k < faces.n_quadrature_points;++k)
	{
	  for (unsigned int i=0;i<polynomials_abf[0]->n() * dim;++i)
	    {
	      boundary_weights_abf(k,i) = polynomials_abf[i%dim]->
		compute_value(i / dim, faces.point(k)) * faces.weight(k);
	    }
	}
    }
  
  // Create Legendre basis for the
  // space D_xi Q_k
  if (deg>0)
    {
      std::vector<AnisotropicPolynomials<dim>* > polynomials(dim);
      
      for (unsigned int dd=0;dd<dim;++dd)
	{
	  std::vector<std::vector<Polynomials::Polynomial<double> > > poly(dim);
	  for (unsigned int d=0;d<dim;++d)
	    poly[d] = Polynomials::Legendre::generate_complete_basis(deg);
	  poly[dd] = Polynomials::Legendre::generate_complete_basis(deg-1);

	  polynomials[dd] = new AnisotropicPolynomials<dim>(poly);
	}
      
      interior_weights.reinit(TableIndices<3>(n_interior_points, polynomials[0]->n(), dim));
      
      for (unsigned int k=0;k<cell_quadrature.n_quadrature_points;++k)
	{
	  for (unsigned int i=0;i<polynomials[0]->n();++i)
	    for (unsigned int d=0;d<dim;++d)
	      interior_weights(k,i,d) = cell_quadrature.weight(k)
		* polynomials[d]->compute_value(i,cell_quadrature.point(k));
	}
      
      for (unsigned int d=0;d<dim;++d)
	delete polynomials[d];
    }
  

  // Decouple the creation of the generalized support points 
  // from computation of interior weights.
  for (unsigned int k=0;k<cell_quadrature.n_quadrature_points;++k)
    this->generalized_support_points[current++] = cell_quadrature.point(k);

  // Additional functionality for the ABF elements
  // TODO: Here the canonical extension of the principle
  // behind the ABF elements is implemented. It is unclear,
  // if this really leads to the ABF spaces in 3D!
  interior_weights_abf.reinit(TableIndices<3>(cell_quadrature.n_quadrature_points, 
					      polynomials_abf[0]->n() * dim, dim));
  Tensor<1, dim> poly_grad;

  for (unsigned int k=0;k<cell_quadrature.n_quadrature_points;++k)
    {
      for (unsigned int i=0;i<polynomials_abf[0]->n() * dim;++i)
	{
	  poly_grad = polynomials_abf[i%dim]->compute_grad(i / dim,cell_quadrature.point(k))
	    * cell_quadrature.weight(k);
	  // The minus sign comes from the use of the Gauss theorem to replace the divergence.
	  for (unsigned int d=0;d<dim;++d)
	    interior_weights_abf(k,i,d) = -poly_grad[d];
	}
    }

  for (unsigned int d=0;d<dim;++d)
    delete polynomials_abf[d];

  Assert (current == this->generalized_support_points.size(),
	  ExcInternalError());
}

#endif


#if deal_II_dimension == 1

template <>
std::vector<unsigned int>
FE_ABF<1>::get_dpo_vector (const unsigned int)
{
  Assert (false, ExcImpossibleInDim(1));
  return std::vector<unsigned int>();
}

#endif


template <int dim>
std::vector<unsigned int>
FE_ABF<dim>::get_dpo_vector (const unsigned int rt_order)
{
                                   // the element is face-based (not
                                   // to be confused with George
                                   // W. Bush's Faith Based
                                   // Initiative...), and we have
                                   // (rt_order+1)^(dim-1) DoFs per face
  unsigned int dofs_per_face = 1;
  for (unsigned int d=0; d<dim-1; ++d)
    dofs_per_face *= rt_order+1;

                                   // and then there are interior dofs
  const unsigned int
    interior_dofs = dim*(rt_order+1)*dofs_per_face;
  
  std::vector<unsigned int> dpo(dim+1);
  dpo[dim-1] = dofs_per_face;
  dpo[dim]   = interior_dofs;
  
  return dpo;
}



template <int dim>
UpdateFlags
FE_ABF<dim>::update_once (const UpdateFlags) const
{
				   // even the values have to be
				   // computed on the real cell, so
				   // nothing can be done in advance
  return update_default;
}



template <int dim>
UpdateFlags
FE_ABF<dim>::update_each (const UpdateFlags flags) const
{
  UpdateFlags out = update_default;


  if (flags & update_values)
    out |= update_values             | update_covariant_transformation
                                     | update_contravariant_transformation 
                                     | update_JxW_values;
  if (flags & update_gradients)
    out |= update_gradients          | update_covariant_transformation 
                                     | update_contravariant_transformation
                                     | update_JxW_values;
  //TODO: Set update flags appropriately and figure out, how the second
  // derivatives for the RT elements can be computed correctly.
  if (flags & update_second_derivatives)
    out |= update_second_derivatives | update_contravariant_transformation;

  return out;
}

//---------------------------------------------------------------------------
// Data field initialization
//---------------------------------------------------------------------------




template <int dim>
unsigned int
FE_ABF<dim>::n_base_elements () const
{
  return 1;
}



template <int dim>
const FiniteElement<dim> &
FE_ABF<dim>::base_element (const unsigned int index) const
{
  Assert (index==0, ExcIndexRange(index, 0, 1));
  return *this;
}



template <int dim>
unsigned int
FE_ABF<dim>::element_multiplicity (const unsigned int index) const
{
  Assert (index==0, ExcIndexRange(index, 0, 1));
  return 1;
}



template <int dim>
bool
FE_ABF<dim>::has_support_on_face (const unsigned int shape_index,
                                            const unsigned int face_index) const
{
  Assert (shape_index < this->dofs_per_cell,
	  ExcIndexRange (shape_index, 0, this->dofs_per_cell));
  Assert (face_index < GeometryInfo<dim>::faces_per_cell,
	  ExcIndexRange (face_index, 0, GeometryInfo<dim>::faces_per_cell));

				   // Return computed values if we
				   // know them easily. Otherwise, it
				   // is always safe to return true.
  switch (rt_order)
    {
      case 0:
      {
        switch (dim)
          {
            case 2:
            {
                                               // only on the one
                                               // non-adjacent face
                                               // are the values
                                               // actually zero. list
                                               // these in a table
              return (face_index != GeometryInfo<dim>::opposite_face[shape_index]);
            }
            
            default:
	      return true;
          };
      };
      
      default:  // other rt_order
	return true;
    };
  
  return true;
}



template <int dim>
void
FE_ABF<dim>::interpolate(
  std::vector<double>&,
  const std::vector<double>&) const
{
  Assert(false, ExcNotImplemented());
}



template <int dim>
void
FE_ABF<dim>::interpolate(
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

  std::fill(local_dofs.begin(), local_dofs.end(), 0.);

  const unsigned int n_face_points = boundary_weights.size(0);
  for (unsigned int face=0;face<GeometryInfo<dim>::faces_per_cell;++face)
    for (unsigned int k=0;k<n_face_points;++k)
      for (unsigned int i=0;i<boundary_weights.size(1);++i)
      {
	local_dofs[i+face*this->dofs_per_face] += boundary_weights(k,i)
			 * values[face*n_face_points+k](GeometryInfo<dim>::unit_normal_direction[face]+offset);
      }
  
  const unsigned start_cell_dofs = GeometryInfo<dim>::faces_per_cell*this->dofs_per_face;
  const unsigned start_cell_points = GeometryInfo<dim>::faces_per_cell*n_face_points;
  
  for (unsigned int k=0;k<interior_weights.size(0);++k)
    for (unsigned int i=0;i<interior_weights.size(1);++i)
      for (unsigned int d=0;d<dim;++d)
	local_dofs[start_cell_dofs+i*dim+d] += interior_weights(k,i,d) * values[k+start_cell_points](d+offset);

  //TODO: Insert missing code for ABF elements. (cf. other interpolate method)
}


template <int dim>
void
FE_ABF<dim>::interpolate(
  std::vector<double>& local_dofs,
  const VectorSlice<const std::vector<std::vector<double> > >& values) const
{
  Assert (values.size() == this->n_components(),
	  ExcDimensionMismatch(values.size(), this->n_components()));
  Assert (values[0].size() == this->generalized_support_points.size(),
	  ExcDimensionMismatch(values[0].size(), this->generalized_support_points.size()));
  Assert (local_dofs.size() == this->dofs_per_cell,
	  ExcDimensionMismatch(local_dofs.size(),this->dofs_per_cell));

  std::fill(local_dofs.begin(), local_dofs.end(), 0.);

  const unsigned int n_face_points = boundary_weights.size(0);
  for (unsigned int face=0;face<GeometryInfo<dim>::faces_per_cell;++face)
    for (unsigned int k=0;k<n_face_points;++k)
      for (unsigned int i=0;i<boundary_weights.size(1);++i)
      {
	local_dofs[i+face*this->dofs_per_face] += boundary_weights(k,i)
			 * values[GeometryInfo<dim>::unit_normal_direction[face]][face*n_face_points+k];
      }
  
  const unsigned start_cell_dofs = GeometryInfo<dim>::faces_per_cell*this->dofs_per_face;
  const unsigned start_cell_points = GeometryInfo<dim>::faces_per_cell*n_face_points;

  for (unsigned int k=0;k<interior_weights.size(0);++k)
    for (unsigned int i=0;i<interior_weights.size(1);++i)
      for (unsigned int d=0;d<dim;++d)
	local_dofs[start_cell_dofs+i*dim+d] += interior_weights(k,i,d) * values[d][k+start_cell_points];

  const unsigned start_abf_dofs = start_cell_dofs + interior_weights.size(1) * dim;

  // Cell integral of ABF terms
  for (unsigned int k=0;k<interior_weights_abf.size(0);++k)
    for (unsigned int i=0;i<interior_weights_abf.size(1);++i)
      for (unsigned int d=0;d<dim;++d)
	local_dofs[start_abf_dofs+i] += interior_weights_abf(k,i,d) * values[d][k+start_cell_points];

  // Face integral of ABF terms
  for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
    {
      double n_orient = (double) GeometryInfo<dim>::unit_normal_orientation[face];
      for (unsigned int fp=0; fp < n_face_points; ++fp)
	{
	  // TODO: Check what the face_orientation has to be in 3D
	  unsigned int k = QProjector<dim>::DataSetDescriptor::face (face, false, n_face_points);
	  for (unsigned int i=0; i<boundary_weights_abf.size(1); ++i)
	    local_dofs[start_abf_dofs+i] += n_orient * boundary_weights_abf(k + fp, i) 
	      * values[GeometryInfo<dim>::unit_normal_direction[face]][k + fp];
	}
    }

  // TODO: Check if this "correction" can be removed.
  for (unsigned int i=0; i<boundary_weights_abf.size(1); ++i)
    if (fabs (local_dofs[start_abf_dofs+i]) < 1.0e-16)
      local_dofs[start_abf_dofs+i] = 0.0;
}


template <int dim>
unsigned int
FE_ABF<dim>::memory_consumption () const
{
  Assert (false, ExcNotImplemented ());
  return 0;
}


template class FE_ABF<deal_II_dimension>;
