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


#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/qprojector.h>
#include <deal.II/base/polynomials_p.h>
#include <deal.II/base/table.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/fe/fe.h>
#include <deal.II/fe/mapping.h>
#include <deal.II/fe/fe_bdm.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_tools.h>

#include <iostream>
#include <sstream>


DEAL_II_NAMESPACE_OPEN

template <int dim>
FE_BDM<dim>::FE_BDM (const unsigned int deg)
  :
  FE_PolyTensor<PolynomialsBDM<dim>, dim> (
    deg,
    FiniteElementData<dim>(get_dpo_vector(deg),
                           dim, deg+1, FiniteElementData<dim>::Hdiv, 1),
    get_ria_vector (deg),
    std::vector<ComponentMask>(PolynomialsBDM<dim>::compute_n_pols(deg),
                               std::vector<bool>(dim,true)))
{
  Assert (dim >= 2, ExcImpossibleInDim(dim));
  Assert (dim<3, ExcNotImplemented());
  Assert (deg > 0, ExcMessage("Lowest order BDM element are degree 1, but you asked for degree 0"));

  const unsigned int n_dofs = this->dofs_per_cell;

  this->mapping_type = mapping_bdm;
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

//   std::cout << std::endl;
//   M.print_formatted(std::cout, 2, true);

  this->inverse_node_matrix.reinit(n_dofs, n_dofs);
  this->inverse_node_matrix.invert(M);
  // From now on, the shape functions
  // will be the correct ones, not
  // the raw shape functions anymore.

  this->reinit_restriction_and_prolongation_matrices(true, true);
  FETools::compute_embedding_matrices (*this, this->prolongation, true);

  FullMatrix<double> face_embeddings[GeometryInfo<dim>::max_children_per_face];
  for (unsigned int i=0; i<GeometryInfo<dim>::max_children_per_face; ++i)
    face_embeddings[i].reinit (this->dofs_per_face, this->dofs_per_face);
  FETools::compute_face_embedding_matrices(*this, face_embeddings, 0, 0);
  this->interface_constraints.reinit((1<<(dim-1)) * this->dofs_per_face,
                                     this->dofs_per_face);
  unsigned int target_row=0;
  for (unsigned int d=0; d<GeometryInfo<dim>::max_children_per_face; ++d)
    for (unsigned int i=0; i<face_embeddings[d].m(); ++i)
      {
        for (unsigned int j=0; j<face_embeddings[d].n(); ++j)
          this->interface_constraints(target_row,j) = face_embeddings[d](i,j);
        ++target_row;
      }
}



template <int dim>
std::string
FE_BDM<dim>::get_name () const
{
  // note that the
  // FETools::get_fe_from_name
  // function depends on the
  // particular format of the string
  // this function returns, so they
  // have to be kept in synch

  // note that this->degree is the maximal
  // polynomial degree and is thus one higher
  // than the argument given to the
  // constructor
  std::ostringstream namebuf;
  namebuf << "FE_BDM<" << dim << ">(" << this->degree-1 << ")";

  return namebuf.str();
}


template <int dim>
FiniteElement<dim> *
FE_BDM<dim>::clone() const
{
  return new FE_BDM<dim>(*this);
}



template <int dim>
void
FE_BDM<dim>::interpolate(
  std::vector<double> &,
  const std::vector<double> &) const
{
  Assert(false, ExcNotImplemented());
}


template <int dim>
void
FE_BDM<dim>::interpolate(
  std::vector<double> &,
  const std::vector<Vector<double> > &,
  unsigned int) const
{
  Assert(false, ExcNotImplemented());
}



template <int dim>
void
FE_BDM<dim>::interpolate(
  std::vector<double> &local_dofs,
  const VectorSlice<const std::vector<std::vector<double> > > &values) const
{
  AssertDimension (values.size(), dim);
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
  for (; f<GeometryInfo<dim>::faces_per_cell;
       ++f, fbase+=this->dofs_per_face)
    {
      for (unsigned int i=0; i<this->dofs_per_face; ++i)
        {
          local_dofs[fbase+i] = values[GeometryInfo<dim>::unit_normal_direction[f]][fbase+i];
        }
    }

  // Done for BDM1
  if (fbase == this->dofs_per_cell) return;

  // What's missing are the interior
  // degrees of freedom. In each
  // point, we take all components of
  // the solution.
  Assert ((this->dofs_per_cell - fbase) % dim == 0, ExcInternalError());

  // Here, the number of the point
  // and of the shape function
  // coincides. This will change
  // below, since we have more
  // support points than test
  // functions in the interior.
  const unsigned int pbase = fbase;
  for (unsigned int d=0; d<dim; ++d, fbase += test_values[0].size())
    {
      for (unsigned int i=0; i<test_values[0].size(); ++i)
        {
          local_dofs[fbase+i] = 0.;
          for (unsigned int k=0; k<test_values.size(); ++k)
            local_dofs[fbase+i] += values[d][pbase+k] * test_values[k][i];
        }
    }

  Assert (fbase == this->dofs_per_cell, ExcInternalError());
}




template <int dim>
std::vector<unsigned int>
FE_BDM<dim>::get_dpo_vector (const unsigned int deg)
{
  // the element is face-based and we have
  // (deg+1)^(dim-1) DoFs per face
  unsigned int dofs_per_face = 1;
  for (unsigned int d=1; d<dim; ++d)
    dofs_per_face *= deg+1;

  // and then there are interior dofs
  unsigned int
  interior_dofs = dim==1?deg:dim*deg*(deg-1)/2;
  if (dim>2)
    {
      interior_dofs *= deg-2;
      interior_dofs /= 3;
    }

  std::vector<unsigned int> dpo(dim+1);
  dpo[dim-1] = dofs_per_face;
  dpo[dim]   = interior_dofs;

  return dpo;
}



template <int dim>
std::vector<bool>
FE_BDM<dim>::get_ria_vector (const unsigned int deg)
{
  if (dim==1)
    {
      Assert (false, ExcImpossibleInDim(1));
      return std::vector<bool>();
    }

  Assert(dim==2, ExcNotImplemented());
  const unsigned int dofs_per_cell = PolynomialsBDM<dim>::compute_n_pols(deg);
  const unsigned int dofs_per_face = deg+1;
  // all dofs need to be
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
FE_BDM<dim>::initialize_support_points (const unsigned int deg)
{
  // interior point in 1d
  unsigned int npoints = deg;
  // interior point in 2d
  if (dim >= 2)
    {
      npoints *= deg;
//      npoints /= 2;
    }
  // interior point in 2d
  if (dim >= 3)
    {
      npoints *= deg;
//      npoints /= 3;
    }
  npoints += GeometryInfo<dim>::faces_per_cell * this->dofs_per_face;

  this->generalized_support_points.resize (npoints);
  this->generalized_face_support_points.resize (this->dofs_per_face);

  // Number of the point being entered
  unsigned int current = 0;

  // On the faces, we choose as many
  // Gauss points as necessary to
  // determine the normal component
  // uniquely. This is the deg of
  // the BDM element plus
  // one.
  if (dim>1)
    {
      QGauss<dim-1> face_points (deg+1);
      Assert (face_points.size() == this->dofs_per_face,
              ExcInternalError());
      for (unsigned int k=0; k<this->dofs_per_face; ++k)
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

  if (deg<=1) return;
  // Although the polynomial space is
  // only P_{k-2}, we use the tensor
  // product points for Q_{k-2}
  QGauss<dim> quadrature(deg);

  // Remember where interior points start
  const unsigned int ibase=current;
//  for (unsigned int k=0;k<deg-1;++k)
  for (unsigned int j=0; j<deg; ++j)
    for (unsigned int i=0; i<deg; ++i)
      {
        this->generalized_support_points[current] = quadrature.point(current-ibase);
        ++current;
      }
  Assert(current == npoints, ExcInternalError());

  // Finaly, compute the values of
  // the test functios in the
  // interior quadrature points
  PolynomialsP<dim> poly(deg-2);

  test_values.resize(quadrature.size());
  std::vector<Tensor<1,dim> > dummy1;
  std::vector<Tensor<2,dim> > dummy2;

  for (unsigned int k=0; k<quadrature.size(); ++k)
    {
      test_values[k].resize(poly.n());
      poly.compute(quadrature.point(k), test_values[k], dummy1, dummy2);
      for (unsigned int i=0; i < poly.n(); ++i)
        {
          test_values[k][i] *= quadrature.weight(k);
        }
    }
}



/*-------------- Explicit Instantiations -------------------------------*/
#include "fe_bdm.inst"

DEAL_II_NAMESPACE_CLOSE
