// ---------------------------------------------------------------------
//
// Copyright (C) 2003 - 2016 by the deal.II authors
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
                           dim, deg+1, FiniteElementData<dim>::Hdiv),
    get_ria_vector (deg),
    std::vector<ComponentMask>(PolynomialsBDM<dim>::compute_n_pols(deg),
                               std::vector<bool>(dim,true)))
{
  Assert (dim >= 2, ExcImpossibleInDim(dim));
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

  // Embedding errors become pretty large, so we just replace the
  // regular threshold in both "computing_..." functions by 1.
  this->reinit_restriction_and_prolongation_matrices(true, true);
  FETools::compute_embedding_matrices (*this, this->prolongation, true, 1.);

  FullMatrix<double> face_embeddings[GeometryInfo<dim>::max_children_per_face];
  for (unsigned int i=0; i<GeometryInfo<dim>::max_children_per_face; ++i)
    face_embeddings[i].reinit (this->dofs_per_face, this->dofs_per_face);
  FETools::compute_face_embedding_matrices(*this, face_embeddings, 0, 0, 1.);
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
  // FETools::get_fe_by_name
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

  // First do interpolation on faces. There, the component evaluated
  // depends on the face direction and orientation.

  // The index of the first dof on this face or the cell
  unsigned int dbase = 0;
  // The index of the first generalized support point on this face or the cell
  unsigned int pbase = 0;
  for (unsigned int f = 0; f<GeometryInfo<dim>::faces_per_cell; ++f)
    {
      // Old version with no moments in 2D. See comment below in
      // initialize_support_points()
      if (test_values_face.size() == 0)
        {
          for (unsigned int i=0; i<this->dofs_per_face; ++i)
            local_dofs[dbase+i] = values[GeometryInfo<dim>::unit_normal_direction[f]][pbase+i];
          pbase += this->dofs_per_face;
        }
      else
        {
          for (unsigned int i=0; i<this->dofs_per_face; ++i)
            {
              double s = 0.;
              for (unsigned int k=0; k<test_values_face.size(); ++k)
                s += values[GeometryInfo<dim>::unit_normal_direction[f]][pbase+k] * test_values_face[k][i];
              local_dofs[dbase+i] = s;
            }
          pbase += test_values_face.size();
        }
      dbase += this->dofs_per_face;
    }

  AssertDimension (dbase, this->dofs_per_face * GeometryInfo<dim>::faces_per_cell);
  AssertDimension (pbase, this->generalized_support_points.size() - test_values_cell.size());

  // Done for BDM1
  if (dbase == this->dofs_per_cell) return;

  // What's missing are the interior
  // degrees of freedom. In each
  // point, we take all components of
  // the solution.
  Assert ((this->dofs_per_cell - dbase) % dim == 0, ExcInternalError());

  for (unsigned int d=0; d<dim; ++d, dbase += test_values_cell[0].size())
    {
      for (unsigned int i=0; i<test_values_cell[0].size(); ++i)
        {
          double s = 0.;
          for (unsigned int k=0; k<test_values_cell.size(); ++k)
            s += values[d][pbase+k] * test_values_cell[k][i];
          local_dofs[dbase+i] = s;
        }
    }

  Assert (dbase == this->dofs_per_cell, ExcInternalError());
}




template <int dim>
std::vector<unsigned int>
FE_BDM<dim>::get_dpo_vector (const unsigned int deg)
{
  // the element is face-based and we have as many degrees of freedom
  // on the faces as there are polynomials of degree up to
  // deg. Observe the odd convention of
  // PolynomialSpace::compute_n_pols()!
  unsigned int dofs_per_face = PolynomialSpace<dim-1>::compute_n_pols(deg+1);

  // and then there are interior dofs, namely the number of
  // polynomials up to degree deg-2 in dim dimensions.
  unsigned int interior_dofs = 0;
  if (deg>1)
    interior_dofs = dim * PolynomialSpace<dim>::compute_n_pols(deg-1);

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

  const unsigned int dofs_per_cell = PolynomialsBDM<dim>::compute_n_pols(deg);
  const unsigned int dofs_per_face = PolynomialSpace<dim-1>::compute_n_pols(deg);

  Assert(GeometryInfo<dim>::faces_per_cell*dofs_per_face < dofs_per_cell,
         ExcInternalError());

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


namespace
{
  // This function sets up the values of the polynomials we want to
  // take moments with in the quadrature points. In fact, we multiply
  // thos by the weights, such that the sum of function values and
  // test_values over quadrature points yields the interpolated degree
  // of freedom.
  template <int dim>
  void
  initialize_test_values (std::vector<std::vector<double> > &test_values,
                          const Quadrature<dim> &quadrature,
                          const unsigned int deg)
  {
    PolynomialsP<dim> poly(deg);
    std::vector<Tensor<1,dim> > dummy1;
    std::vector<Tensor<2,dim> > dummy2;
    std::vector<Tensor<3,dim> > dummy3;
    std::vector<Tensor<4,dim> > dummy4;

    test_values.resize(quadrature.size());

    for (unsigned int k=0; k<quadrature.size(); ++k)
      {
        test_values[k].resize(poly.n());
        poly.compute(quadrature.point(k), test_values[k], dummy1, dummy2,
                     dummy3, dummy4);
        for (unsigned int i=0; i < poly.n(); ++i)
          {
            test_values[k][i] *= quadrature.weight(k);
          }
      }
  }

  // This specialization only serves to avoid error messages. Nothing
  // useful can be computed in dimension zero and thus the vector
  // length stays zero.
  template <>
  void
  initialize_test_values (std::vector<std::vector<double> > &,
                          const Quadrature<0> &,
                          const unsigned int)
  {}
}


template <int dim>
void
FE_BDM<dim>::initialize_support_points (const unsigned int deg)
{
  // Our support points are quadrature points on faces and inside the
  // cell. First on the faces, we have to test polynomials of degree
  // up to deg, which means we need dg+1 points in each direction. The
  // fact that we do not have tensor product polynomials will be
  // considered later. In 2D, we can use point values.
  QGauss<dim-1> face_points (deg+1);

  // Copy the quadrature formula to the face points.
  this->generalized_face_support_points.resize (face_points.size());
  for (unsigned int k=0; k<face_points.size(); ++k)
    this->generalized_face_support_points[k] = face_points.point(k);

  // In the interior, we only test with polynomials of degree up to
  // deg-2, thus we use deg points. Note that deg>=1 and the lowest
  // order element has no points in the cell, such that we have to
  // distinguish this case.
  QGauss<dim> cell_points(deg==1 ? 0 : deg);

  // Compute the size of the whole support point set
  const unsigned int npoints
    = cell_points.size() + GeometryInfo<dim>::faces_per_cell * face_points.size();

  this->generalized_support_points.resize (npoints);

  Quadrature<dim> faces = QProjector<dim>::project_to_all_faces(face_points);
  for (unsigned int k=0; k < face_points.size()*GeometryInfo<dim>::faces_per_cell; ++k)
    this->generalized_support_points[k]
      = faces.point(k+QProjector<dim>
                    ::DataSetDescriptor::face(0, true, false, false,
                                              this->dofs_per_face));

  // Currently, for backward compatibility, we do not use moments, but
  // point values on faces in 2D. In 3D, this is impossible, since the
  // moments are only taken with respect to PolynomialsP.
  if (dim>2)
    initialize_test_values(test_values_face, face_points, deg);

  if (deg<=1) return;

  // Remember where interior points start
  const unsigned int ibase = face_points.size()*GeometryInfo<dim>::faces_per_cell;
  for (unsigned int k=0; k<cell_points.size(); ++k)
    {
      this->generalized_support_points[ibase+k] = cell_points.point(k);
    }
  // Finally, compute the values of
  // the test functions in the
  // interior quadrature points

  initialize_test_values(test_values_cell, cell_points, deg-2);
}



/*-------------- Explicit Instantiations -------------------------------*/
#include "fe_bdm.inst"

DEAL_II_NAMESPACE_CLOSE

