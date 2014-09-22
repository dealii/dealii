// ---------------------------------------------------------------------
//
// Copyright (C) 2003 - 2013 by the deal.II authors
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


#include <deal.II/base/utilities.h>
#include <deal.II/base/quadrature.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/qprojector.h>
#include <deal.II/base/table.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/fe/fe.h>
#include <deal.II/fe/mapping.h>
#include <deal.II/fe/fe_abf.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_tools.h>

#include <sstream>
#include <iostream>

//TODO: implement the adjust_quad_dof_index_for_face_orientation_table and
//adjust_line_dof_index_for_line_orientation_table fields, and write tests
//similar to bits/face_orientation_and_fe_q_*


DEAL_II_NAMESPACE_OPEN


template <int dim>
FE_ABF<dim>::FE_ABF (const unsigned int deg)
  :
  FE_PolyTensor<PolynomialsABF<dim>, dim> (
    deg,
    FiniteElementData<dim>(get_dpo_vector(deg),
                           dim, deg+1, FiniteElementData<dim>::Hdiv, 1),
    std::vector<bool>(PolynomialsABF<dim>::compute_n_pols(deg), true),
    std::vector<ComponentMask>(PolynomialsABF<dim>::compute_n_pols(deg),
                               std::vector<bool>(dim,true))),
  rt_order(deg)
{
  Assert (dim >= 2, ExcImpossibleInDim(dim));
  const unsigned int n_dofs = this->dofs_per_cell;

  this->mapping_type = mapping_raviart_thomas;
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

  this->inverse_node_matrix.reinit(n_dofs, n_dofs);
  this->inverse_node_matrix.invert(M);
  // From now on, the shape functions
  // will be the correct ones, not
  // the raw shape functions anymore.

  // Reinit the vectors of
  // restriction and prolongation
  // matrices to the right sizes.
  // Restriction only for isotropic
  // refinement
  this->reinit_restriction_and_prolongation_matrices(true);
  // Fill prolongation matrices with embedding operators
  FETools::compute_embedding_matrices (*this, this->prolongation);

  initialize_restriction ();

  // TODO[TL]: for anisotropic refinement we will probably need a table of submatrices with an array for each refine case
  std::vector<FullMatrix<double> >
  face_embeddings(1<<(dim-1), FullMatrix<double>(this->dofs_per_face,
                                                 this->dofs_per_face));
  // TODO: Something goes wrong there. The error of the least squares fit
  // is to large ...
  // FETools::compute_face_embedding_matrices(*this, &face_embeddings[0], 0, 0);
  this->interface_constraints.reinit((1<<(dim-1)) * this->dofs_per_face,
                                     this->dofs_per_face);
  unsigned int target_row=0;
  for (unsigned int d=0; d<face_embeddings.size(); ++d)
    for (unsigned int i=0; i<face_embeddings[d].m(); ++i)
      {
        for (unsigned int j=0; j<face_embeddings[d].n(); ++j)
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

  std::ostringstream namebuf;

  namebuf << "FE_ABF<" << dim << ">(" << rt_order << ")";

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



// Version for 2d and higher. See above for 1d version
template <int dim>
void
FE_ABF<dim>::initialize_support_points (const unsigned int deg)
{
  QGauss<dim> cell_quadrature(deg+2);
  const unsigned int n_interior_points = cell_quadrature.size();

  unsigned int n_face_points = (dim>1) ? 1 : 0;
  // compute (deg+1)^(dim-1)
  for (unsigned int d=1; d<dim; ++d)
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
      for (unsigned int d=0; d<dim; ++d)
        poly[d].push_back (Polynomials::Monomial<double> (deg+1));
      poly[dd] = Polynomials::Monomial<double>::generate_complete_basis(deg);

      polynomials_abf[dd] = new AnisotropicPolynomials<dim>(poly);
    }

  // Number of the point being entered
  unsigned int current = 0;

  if (dim>1)
    {
      QGauss<dim-1> face_points (deg+1);
      TensorProductPolynomials<dim-1> legendre =
        Polynomials::Legendre::generate_complete_basis(deg);

      boundary_weights.reinit(n_face_points, legendre.n());

//       Assert (face_points.size() == this->dofs_per_face,
//            ExcInternalError());

      for (unsigned int k=0; k<n_face_points; ++k)
        {
          this->generalized_face_support_points[k] = face_points.point(k);
          // Compute its quadrature
          // contribution for each
          // moment.
          for (unsigned int i=0; i<legendre.n(); ++i)
            {
              boundary_weights(k, i)
                = face_points.weight(k)
                  * legendre.compute_value(i, face_points.point(k));
            }
        }

      Quadrature<dim> faces = QProjector<dim>::project_to_all_faces(face_points);
      for (; current<GeometryInfo<dim>::faces_per_cell*n_face_points;
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
      boundary_weights_abf.reinit(faces.size(), polynomials_abf[0]->n() * dim);
      for (unsigned int k=0; k < faces.size(); ++k)
        {
          for (unsigned int i=0; i<polynomials_abf[0]->n() * dim; ++i)
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

      for (unsigned int dd=0; dd<dim; ++dd)
        {
          std::vector<std::vector<Polynomials::Polynomial<double> > > poly(dim);
          for (unsigned int d=0; d<dim; ++d)
            poly[d] = Polynomials::Legendre::generate_complete_basis(deg);
          poly[dd] = Polynomials::Legendre::generate_complete_basis(deg-1);

          polynomials[dd] = new AnisotropicPolynomials<dim>(poly);
        }

      interior_weights.reinit(TableIndices<3>(n_interior_points, polynomials[0]->n(), dim));

      for (unsigned int k=0; k<cell_quadrature.size(); ++k)
        {
          for (unsigned int i=0; i<polynomials[0]->n(); ++i)
            for (unsigned int d=0; d<dim; ++d)
              interior_weights(k,i,d) = cell_quadrature.weight(k)
                                        * polynomials[d]->compute_value(i,cell_quadrature.point(k));
        }

      for (unsigned int d=0; d<dim; ++d)
        delete polynomials[d];
    }


  // Decouple the creation of the generalized support points
  // from computation of interior weights.
  for (unsigned int k=0; k<cell_quadrature.size(); ++k)
    this->generalized_support_points[current++] = cell_quadrature.point(k);

  // Additional functionality for the ABF elements
  // TODO: Here the canonical extension of the principle
  // behind the ABF elements is implemented. It is unclear,
  // if this really leads to the ABF spaces in 3D!
  interior_weights_abf.reinit(TableIndices<3>(cell_quadrature.size(),
                                              polynomials_abf[0]->n() * dim, dim));
  Tensor<1, dim> poly_grad;

  for (unsigned int k=0; k<cell_quadrature.size(); ++k)
    {
      for (unsigned int i=0; i<polynomials_abf[0]->n() * dim; ++i)
        {
          poly_grad = polynomials_abf[i%dim]->compute_grad(i / dim,cell_quadrature.point(k))
                      * cell_quadrature.weight(k);
          // The minus sign comes from the use of the Gauss theorem to replace the divergence.
          for (unsigned int d=0; d<dim; ++d)
            interior_weights_abf(k,i,d) = -poly_grad[d];
        }
    }

  for (unsigned int d=0; d<dim; ++d)
    delete polynomials_abf[d];

  Assert (current == this->generalized_support_points.size(),
          ExcInternalError());
}



// This function is the same Raviart-Thomas interpolation performed by
// interpolate. Still, we cannot use interpolate, since it was written
// for smooth functions. The functions interpolated here are not
// smooth, maybe even not continuous. Therefore, we must double the
// number of quadrature points in each direction in order to integrate
// only smooth functions.

// Then again, the interpolated function is chosen such that the
// moments coincide with the function to be interpolated.

template <int dim>
void
FE_ABF<dim>::initialize_restriction()
{
  if (dim==1)
    {
      unsigned int iso=RefinementCase<dim>::isotropic_refinement-1;
      for (unsigned int i=0; i<GeometryInfo<dim>::max_children_per_cell; ++i)
        this->restriction[iso][i].reinit(0,0);
      return;
    }
  unsigned int iso=RefinementCase<dim>::isotropic_refinement-1;
  QGauss<dim-1> q_base (rt_order+1);
  const unsigned int n_face_points = q_base.size();
  // First, compute interpolation on
  // subfaces
  for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
    {
      // The shape functions of the
      // child cell are evaluated
      // in the quadrature points
      // of a full face.
      Quadrature<dim> q_face
        = QProjector<dim>::project_to_face(q_base, face);
      // Store shape values, since the
      // evaluation suffers if not
      // ordered by point
      Table<2,double> cached_values(this->dofs_per_cell, q_face.size());
      for (unsigned int k=0; k<q_face.size(); ++k)
        for (unsigned int i = 0; i < this->dofs_per_cell; ++i)
          cached_values(i,k)
            = this->shape_value_component(i, q_face.point(k),
                                          GeometryInfo<dim>::unit_normal_direction[face]);

      for (unsigned int sub=0; sub<GeometryInfo<dim>::max_children_per_face; ++sub)
        {
          // The weight fuctions for
          // the coarse face are
          // evaluated on the subface
          // only.
          Quadrature<dim> q_sub
            = QProjector<dim>::project_to_subface(q_base, face, sub);
          const unsigned int child
            = GeometryInfo<dim>::child_cell_on_face(
                RefinementCase<dim>::isotropic_refinement, face, sub);

          // On a certain face, we must
          // compute the moments of ALL
          // fine level functions with
          // the coarse level weight
          // functions belonging to
          // that face. Due to the
          // orthogonalization process
          // when building the shape
          // functions, these weights
          // are equal to the
          // corresponding shpe
          // functions.
          for (unsigned int k=0; k<n_face_points; ++k)
            for (unsigned int i_child = 0; i_child < this->dofs_per_cell; ++i_child)
              for (unsigned int i_face = 0; i_face < this->dofs_per_face; ++i_face)
                {
                  // The quadrature
                  // weights on the
                  // subcell are NOT
                  // transformed, so we
                  // have to do it here.
                  this->restriction[iso][child](face*this->dofs_per_face+i_face,
                                                i_child)
                  += Utilities::fixed_power<dim-1>(.5) * q_sub.weight(k)
                     * cached_values(i_child, k)
                     * this->shape_value_component(face*this->dofs_per_face+i_face,
                                                   q_sub.point(k),
                                                   GeometryInfo<dim>::unit_normal_direction[face]);
                }
        }
    }

  if (rt_order==0) return;

  // Create Legendre basis for the
  // space D_xi Q_k. Here, we cannot
  // use the shape functions
  std::vector<AnisotropicPolynomials<dim>* > polynomials(dim);
  for (unsigned int dd=0; dd<dim; ++dd)
    {
      std::vector<std::vector<Polynomials::Polynomial<double> > > poly(dim);
      for (unsigned int d=0; d<dim; ++d)
        poly[d] = Polynomials::Legendre::generate_complete_basis(rt_order);
      poly[dd] = Polynomials::Legendre::generate_complete_basis(rt_order-1);

      polynomials[dd] = new AnisotropicPolynomials<dim>(poly);
    }

  QGauss<dim> q_cell(rt_order+1);
  const unsigned int start_cell_dofs
    = GeometryInfo<dim>::faces_per_cell*this->dofs_per_face;

  // Store shape values, since the
  // evaluation suffers if not
  // ordered by point
  Table<3,double> cached_values(this->dofs_per_cell, q_cell.size(), dim);
  for (unsigned int k=0; k<q_cell.size(); ++k)
    for (unsigned int i = 0; i < this->dofs_per_cell; ++i)
      for (unsigned int d=0; d<dim; ++d)
        cached_values(i,k,d) = this->shape_value_component(i, q_cell.point(k), d);

  for (unsigned int child=0; child<GeometryInfo<dim>::max_children_per_cell; ++child)
    {
      Quadrature<dim> q_sub = QProjector<dim>::project_to_child(q_cell, child);

      for (unsigned int k=0; k<q_sub.size(); ++k)
        for (unsigned int i_child = 0; i_child < this->dofs_per_cell; ++i_child)
          for (unsigned int d=0; d<dim; ++d)
            for (unsigned int i_weight=0; i_weight<polynomials[d]->n(); ++i_weight)
              {
                this->restriction[iso][child](start_cell_dofs+i_weight*dim+d,
                                              i_child)
                += q_sub.weight(k)
                   * cached_values(i_child, k, d)
                   * polynomials[d]->compute_value(i_weight, q_sub.point(k));
              }
    }

  for (unsigned int d=0; d<dim; ++d)
    delete polynomials[d];
}



template <int dim>
std::vector<unsigned int>
FE_ABF<dim>::get_dpo_vector (const unsigned int rt_order)
{
  if (dim == 1)
    {
      Assert (false, ExcImpossibleInDim(1));
      return std::vector<unsigned int>();
    }

  // the element is face-based (not
  // to be confused with George
  // W. Bush's Faith Based
  // Initiative...), and we have
  // (rt_order+1)^(dim-1) DoFs per face
  unsigned int dofs_per_face = 1;
  for (int d=0; d<dim-1; ++d)
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
           | update_piola
           | update_JxW_values;
  if (flags & update_gradients)
    out |= update_gradients          | update_covariant_transformation
           | update_piola
           | update_JxW_values;
  //TODO: Set update flags appropriately and figure out, how the second
  // derivatives for the RT elements can be computed correctly.
  if (flags & update_hessians)
    out |= update_hessians | update_contravariant_transformation;

  return out;
}

//---------------------------------------------------------------------------
// Data field initialization
//---------------------------------------------------------------------------

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
  std::vector<double> &,
  const std::vector<double> &) const
{
  Assert(false, ExcNotImplemented());
}



template <int dim>
void
FE_ABF<dim>::interpolate(
  std::vector<double>    &local_dofs,
  const std::vector<Vector<double> > &values,
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
  for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
    for (unsigned int k=0; k<n_face_points; ++k)
      for (unsigned int i=0; i<boundary_weights.size(1); ++i)
        {
          local_dofs[i+face*this->dofs_per_face] += boundary_weights(k,i)
                                                    * values[face*n_face_points+k](GeometryInfo<dim>::unit_normal_direction[face]+offset);
        }

  const unsigned int start_cell_dofs = GeometryInfo<dim>::faces_per_cell*this->dofs_per_face;
  const unsigned int start_cell_points = GeometryInfo<dim>::faces_per_cell*n_face_points;

  for (unsigned int k=0; k<interior_weights.size(0); ++k)
    for (unsigned int i=0; i<interior_weights.size(1); ++i)
      for (unsigned int d=0; d<dim; ++d)
        local_dofs[start_cell_dofs+i*dim+d] += interior_weights(k,i,d) * values[k+start_cell_points](d+offset);

  const unsigned int start_abf_dofs = start_cell_dofs + interior_weights.size(1) * dim;

  // Cell integral of ABF terms
  for (unsigned int k=0; k<interior_weights_abf.size(0); ++k)
    for (unsigned int i=0; i<interior_weights_abf.size(1); ++i)
      for (unsigned int d=0; d<dim; ++d)
        local_dofs[start_abf_dofs+i] += interior_weights_abf(k,i,d) * values[k+start_cell_points](d+offset);

  // Face integral of ABF terms
  for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
    {
      double n_orient = (double) GeometryInfo<dim>::unit_normal_orientation[face];
      for (unsigned int fp=0; fp < n_face_points; ++fp)
        {
          // TODO: Check what the face_orientation, face_flip and face_rotation  have to be in 3D
          unsigned int k = QProjector<dim>::DataSetDescriptor::face (face, false, false, false, n_face_points);
          for (unsigned int i=0; i<boundary_weights_abf.size(1); ++i)
            local_dofs[start_abf_dofs+i] += n_orient * boundary_weights_abf(k + fp, i)
                                            * values[k + fp](GeometryInfo<dim>::unit_normal_direction[face]+offset);
        }
    }

  // TODO: Check if this "correction" can be removed.
  for (unsigned int i=0; i<boundary_weights_abf.size(1); ++i)
    if (std::fabs (local_dofs[start_abf_dofs+i]) < 1.0e-16)
      local_dofs[start_abf_dofs+i] = 0.0;
}

template <int dim>
void
FE_ABF<dim>::interpolate(
  std::vector<double> &local_dofs,
  const VectorSlice<const std::vector<std::vector<double> > > &values) const
{
  Assert (values.size() == this->n_components(),
          ExcDimensionMismatch(values.size(), this->n_components()));
  Assert (values[0].size() == this->generalized_support_points.size(),
          ExcDimensionMismatch(values[0].size(), this->generalized_support_points.size()));
  Assert (local_dofs.size() == this->dofs_per_cell,
          ExcDimensionMismatch(local_dofs.size(),this->dofs_per_cell));

  std::fill(local_dofs.begin(), local_dofs.end(), 0.);

  const unsigned int n_face_points = boundary_weights.size(0);
  for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
    for (unsigned int k=0; k<n_face_points; ++k)
      for (unsigned int i=0; i<boundary_weights.size(1); ++i)
        {
          local_dofs[i+face*this->dofs_per_face] += boundary_weights(k,i)
                                                    * values[GeometryInfo<dim>::unit_normal_direction[face]][face*n_face_points+k];
        }

  const unsigned int start_cell_dofs = GeometryInfo<dim>::faces_per_cell*this->dofs_per_face;
  const unsigned int start_cell_points = GeometryInfo<dim>::faces_per_cell*n_face_points;

  for (unsigned int k=0; k<interior_weights.size(0); ++k)
    for (unsigned int i=0; i<interior_weights.size(1); ++i)
      for (unsigned int d=0; d<dim; ++d)
        local_dofs[start_cell_dofs+i*dim+d] += interior_weights(k,i,d) * values[d][k+start_cell_points];

  const unsigned int start_abf_dofs = start_cell_dofs + interior_weights.size(1) * dim;

  // Cell integral of ABF terms
  for (unsigned int k=0; k<interior_weights_abf.size(0); ++k)
    for (unsigned int i=0; i<interior_weights_abf.size(1); ++i)
      for (unsigned int d=0; d<dim; ++d)
        local_dofs[start_abf_dofs+i] += interior_weights_abf(k,i,d) * values[d][k+start_cell_points];

  // Face integral of ABF terms
  for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
    {
      double n_orient = (double) GeometryInfo<dim>::unit_normal_orientation[face];
      for (unsigned int fp=0; fp < n_face_points; ++fp)
        {
          // TODO: Check what the face_orientation, face_flip and face_rotation have to be in 3D
          unsigned int k = QProjector<dim>::DataSetDescriptor::face (face, false, false, false, n_face_points);
          for (unsigned int i=0; i<boundary_weights_abf.size(1); ++i)
            local_dofs[start_abf_dofs+i] += n_orient * boundary_weights_abf(k + fp, i)
                                            * values[GeometryInfo<dim>::unit_normal_direction[face]][k + fp];
        }
    }

  // TODO: Check if this "correction" can be removed.
  for (unsigned int i=0; i<boundary_weights_abf.size(1); ++i)
    if (std::fabs (local_dofs[start_abf_dofs+i]) < 1.0e-16)
      local_dofs[start_abf_dofs+i] = 0.0;
}


template <int dim>
std::size_t
FE_ABF<dim>::memory_consumption () const
{
  Assert (false, ExcNotImplemented ());
  return 0;
}



/*-------------- Explicit Instantiations -------------------------------*/
#include "fe_abf.inst"

DEAL_II_NAMESPACE_CLOSE
