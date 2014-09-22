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
#include <deal.II/base/table.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/fe/fe.h>
#include <deal.II/fe/mapping.h>
#include <deal.II/fe/fe_raviart_thomas.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_tools.h>

#include <sstream>


DEAL_II_NAMESPACE_OPEN


template <int dim>
FE_RaviartThomasNodal<dim>::FE_RaviartThomasNodal (const unsigned int deg)
  :
  FE_PolyTensor<PolynomialsRaviartThomas<dim>, dim> (
    deg,
    FiniteElementData<dim>(get_dpo_vector(deg),
                           dim, deg+1, FiniteElementData<dim>::Hdiv, 1),
    get_ria_vector (deg),
    std::vector<ComponentMask>(PolynomialsRaviartThomas<dim>::compute_n_pols(deg),
                               std::vector<bool>(dim,true)))
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

  // Reinit the vectors of
  // prolongation matrices to the
  // right sizes. There are no
  // restriction matrices implemented
  for (unsigned int ref_case=RefinementCase<dim>::cut_x;
       ref_case<RefinementCase<dim>::isotropic_refinement+1; ++ref_case)
    {
      const unsigned int nc = GeometryInfo<dim>::n_children(RefinementCase<dim>(ref_case));

      for (unsigned int i=0; i<nc; ++i)
        this->prolongation[ref_case-1][i].reinit (n_dofs, n_dofs);
    }
  // Fill prolongation matrices with embedding operators
  FETools::compute_embedding_matrices (*this, this->prolongation);
  // TODO[TL]: for anisotropic refinement we will probably need a table of submatrices with an array for each refine case
  FullMatrix<double> face_embeddings[GeometryInfo<dim>::max_children_per_face];
  for (unsigned int i=0; i<GeometryInfo<dim>::max_children_per_face; ++i)
    face_embeddings[i].reinit (this->dofs_per_face, this->dofs_per_face);
  FETools::compute_face_embedding_matrices<dim,double>(*this, face_embeddings, 0, 0);
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
FE_RaviartThomasNodal<dim>::get_name () const
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
  namebuf << "FE_RaviartThomasNodal<" << dim << ">(" << this->degree-1 << ")";

  return namebuf.str();
}


template <int dim>
FiniteElement<dim> *
FE_RaviartThomasNodal<dim>::clone() const
{
  return new FE_RaviartThomasNodal<dim>(*this);
}


//---------------------------------------------------------------------------
// Auxiliary and internal functions
//---------------------------------------------------------------------------



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

  if (deg==0) return;
  // In the interior, we need
  // anisotropic Gauss quadratures,
  // different for each direction.
  QGauss<1> high(deg+1);
  QGauss<1> low(deg);

  for (unsigned int d=0; d<dim; ++d)
    {
      QAnisotropic<dim> *quadrature;
      if (dim == 1) quadrature = new QAnisotropic<dim>(high);
      if (dim == 2) quadrature = new QAnisotropic<dim>(((d==0) ? low : high),
                                                         ((d==1) ? low : high));
      if (dim == 3) quadrature = new QAnisotropic<dim>(((d==0) ? low : high),
                                                         ((d==1) ? low : high),
                                                         ((d==2) ? low : high));
      Assert(dim<=3, ExcNotImplemented());

      for (unsigned int k=0; k<quadrature->size(); ++k)
        this->generalized_support_points[current++] = quadrature->point(k);
      delete quadrature;
    }
  Assert (current == this->dofs_per_cell, ExcInternalError());
}



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



template <>
std::vector<bool>
FE_RaviartThomasNodal<1>::get_ria_vector (const unsigned int)
{
  Assert (false, ExcImpossibleInDim(1));
  return std::vector<bool>();
}



template <int dim>
std::vector<bool>
FE_RaviartThomasNodal<dim>::get_ria_vector (const unsigned int deg)
{
  const unsigned int dofs_per_cell = PolynomialsRaviartThomas<dim>::compute_n_pols(deg);
  unsigned int dofs_per_face = deg+1;
  for (unsigned int d=2; d<dim; ++d)
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
bool
FE_RaviartThomasNodal<dim>::has_support_on_face (
  const unsigned int shape_index,
  const unsigned int face_index) const
{
  Assert (shape_index < this->dofs_per_cell,
          ExcIndexRange (shape_index, 0, this->dofs_per_cell));
  Assert (face_index < GeometryInfo<dim>::faces_per_cell,
          ExcIndexRange (face_index, 0, GeometryInfo<dim>::faces_per_cell));

  // The first degrees of freedom are
  // on the faces and each face has
  // degree degrees.
  const unsigned int support_face = shape_index / this->degree;

  // The only thing we know for sure
  // is that shape functions with
  // support on one face are zero on
  // the opposite face.
  if (support_face < GeometryInfo<dim>::faces_per_cell)
    return (face_index != GeometryInfo<dim>::opposite_face[support_face]);

  // In all other cases, return true,
  // which is safe
  return true;
}


template <int dim>
void
FE_RaviartThomasNodal<dim>::interpolate(
  std::vector<double> &,
  const std::vector<double> &) const
{
  Assert(false, ExcNotImplemented());
}


template <int dim>
void
FE_RaviartThomasNodal<dim>::interpolate(
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
      for (unsigned int i=0; i<istep; ++i)
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
  std::vector<double> &local_dofs,
  const VectorSlice<const std::vector<std::vector<double> > > &values) const
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
  for (; f<GeometryInfo<dim>::faces_per_cell;
       ++f, fbase+=this->dofs_per_face)
    {
      for (unsigned int i=0; i<this->dofs_per_face; ++i)
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
      for (unsigned int i=0; i<istep; ++i)
        {
          local_dofs[fbase+i] = values[f][fbase+i];
        }
      fbase+=istep;
      ++f;
    }
  Assert (fbase == this->dofs_per_cell, ExcInternalError());
}


//TODO: There are tests that check that the following few functions don't produce assertion failures, but none that actually check whether they do the right thing. one example for such a test would be to project a function onto an hp space and make sure that the convergence order is correct with regard to the lowest used polynomial degree

template <int dim>
bool
FE_RaviartThomasNodal<dim>::hp_constraints_are_implemented () const
{
  return true;
}


template <int dim>
std::vector<std::pair<unsigned int, unsigned int> >
FE_RaviartThomasNodal<dim>::hp_vertex_dof_identities (
  const FiniteElement<dim> &fe_other) const
{
  // we can presently only compute these
  // identities if both FEs are
  // FE_RaviartThomasNodals. in that case, no
  // dofs are assigned on the vertex, so we
  // shouldn't be getting here at all.
  if (dynamic_cast<const FE_RaviartThomasNodal<dim>*>(&fe_other)!=0)
    return std::vector<std::pair<unsigned int, unsigned int> > ();
  else
    {
      Assert(false, ExcNotImplemented());
      return std::vector<std::pair<unsigned int, unsigned int> > ();
    }
}



template <int dim>
std::vector<std::pair<unsigned int, unsigned int> >
FE_RaviartThomasNodal<dim>::
hp_line_dof_identities (const FiniteElement<dim> &fe_other) const
{
  // we can presently only compute
  // these identities if both FEs are
  // FE_RaviartThomasNodals
  if (const FE_RaviartThomasNodal<dim> *fe_q_other
      = dynamic_cast<const FE_RaviartThomasNodal<dim>*>(&fe_other))
    {
      // dofs are located on faces; these are
      // only lines in 2d
      if (dim != 2)
        return std::vector<std::pair<unsigned int, unsigned int> >();

      // dofs are located along lines, so two
      // dofs are identical only if in the
      // following two cases (remember that
      // the face support points are Gauss
      // points):
      //1. this->degree = fe_q_other->degree,
      //   in the case, all the dofs on
      //   the line are identical
      //2. this->degree-1 and fe_q_other->degree-1
      //   are both even, i.e. this->dof_per_line
      //   and fe_q_other->dof_per_line are both odd,
      //   there exists only one point (the middle one)
      //   such that dofs are identical on this point
      //
      // to understand this, note that
      // this->degree is the *maximal*
      // polynomial degree, and is thus one
      // higher than the argument given to
      // the constructor
      const unsigned int p = this->degree-1;
      const unsigned int q = fe_q_other->degree-1;

      std::vector<std::pair<unsigned int, unsigned int> > identities;

      if (p==q)
        for (unsigned int i=0; i<p+1; ++i)
          identities.push_back (std::make_pair(i,i));

      else if (p%2==0 && q%2==0)
        identities.push_back(std::make_pair(p/2,q/2));

      return identities;
    }
  else
    {
      Assert (false, ExcNotImplemented());
      return std::vector<std::pair<unsigned int, unsigned int> > ();
    }
}


template <int dim>
std::vector<std::pair<unsigned int, unsigned int> >
FE_RaviartThomasNodal<dim>::hp_quad_dof_identities (
  const FiniteElement<dim> &fe_other) const
{
  // we can presently only compute
  // these identities if both FEs are
  // FE_RaviartThomasNodals
  if (const FE_RaviartThomasNodal<dim> *fe_q_other
      = dynamic_cast<const FE_RaviartThomasNodal<dim>*>(&fe_other))
    {
      // dofs are located on faces; these are
      // only quads in 3d
      if (dim != 3)
        return std::vector<std::pair<unsigned int, unsigned int> >();

      // this works exactly like the line
      // case above
      const unsigned int p = this->dofs_per_quad;
      const unsigned int q = fe_q_other->dofs_per_quad;

      std::vector<std::pair<unsigned int, unsigned int> > identities;

      if (p==q)
        for (unsigned int i=0; i<p; ++i)
          identities.push_back (std::make_pair(i,i));

      else if (p%2!=0 && q%2!=0)
        identities.push_back(std::make_pair(p/2, q/2));

      return identities;
    }
  else
    {
      Assert (false, ExcNotImplemented());
      return std::vector<std::pair<unsigned int, unsigned int> > ();
    }
}


template <int dim>
FiniteElementDomination::Domination
FE_RaviartThomasNodal<dim>::compare_for_face_domination (
  const FiniteElement<dim> &fe_other) const
{
  if (const FE_RaviartThomasNodal<dim> *fe_q_other
      = dynamic_cast<const FE_RaviartThomasNodal<dim>*>(&fe_other))
    {
      if (this->degree < fe_q_other->degree)
        return FiniteElementDomination::this_element_dominates;
      else if (this->degree == fe_q_other->degree)
        return FiniteElementDomination::either_element_can_dominate;
      else
        return FiniteElementDomination::other_element_dominates;
    }

  Assert (false, ExcNotImplemented());
  return FiniteElementDomination::neither_element_dominates;
}



template <>
void
FE_RaviartThomasNodal<1>::get_face_interpolation_matrix (
  const FiniteElement<1,1> &/*x_source_fe*/,
  FullMatrix<double>     &/*interpolation_matrix*/) const
{
  Assert (false, ExcImpossibleInDim(1));
}


template <>
void
FE_RaviartThomasNodal<1>::get_subface_interpolation_matrix (
  const FiniteElement<1,1> &/*x_source_fe*/,
  const unsigned int      /*subface*/,
  FullMatrix<double>     &/*interpolation_matrix*/) const
{
  Assert (false, ExcImpossibleInDim(1));
}



template <int dim>
void
FE_RaviartThomasNodal<dim>::get_face_interpolation_matrix (
  const FiniteElement<dim> &x_source_fe,
  FullMatrix<double>       &interpolation_matrix) const
{
  // this is only implemented, if the
  // source FE is also a
  // RaviartThomasNodal element
  AssertThrow ((x_source_fe.get_name().find ("FE_RaviartThomasNodal<") == 0)
               ||
               (dynamic_cast<const FE_RaviartThomasNodal<dim>*>(&x_source_fe) != 0),
               typename FiniteElement<dim>::
               ExcInterpolationNotImplemented());

  Assert (interpolation_matrix.n() == this->dofs_per_face,
          ExcDimensionMismatch (interpolation_matrix.n(),
                                this->dofs_per_face));
  Assert (interpolation_matrix.m() == x_source_fe.dofs_per_face,
          ExcDimensionMismatch (interpolation_matrix.m(),
                                x_source_fe.dofs_per_face));

  // ok, source is a RaviartThomasNodal element, so
  // we will be able to do the work
  const FE_RaviartThomasNodal<dim> &source_fe
    = dynamic_cast<const FE_RaviartThomasNodal<dim>&>(x_source_fe);

  // Make sure, that the element,
  // for which the DoFs should be
  // constrained is the one with
  // the higher polynomial degree.
  // Actually the procedure will work
  // also if this assertion is not
  // satisfied. But the matrices
  // produced in that case might
  // lead to problems in the
  // hp procedures, which use this
  // method.
  Assert (this->dofs_per_face <= source_fe.dofs_per_face,
          typename FiniteElement<dim>::
          ExcInterpolationNotImplemented ());

  // generate a quadrature
  // with the generalized support points.
  // This is later based as a
  // basis for the QProjector,
  // which returns the support
  // points on the face.
  Quadrature<dim-1> quad_face_support (source_fe.get_generalized_face_support_points ());

  // Rule of thumb for FP accuracy,
  // that can be expected for a
  // given polynomial degree.
  // This value is used to cut
  // off values close to zero.
  double eps = 2e-13*this->degree*(dim-1);

  // compute the interpolation
  // matrix by simply taking the
  // value at the support points.
  const Quadrature<dim> face_projection
    = QProjector<dim>::project_to_face (quad_face_support, 0);

  for (unsigned int i=0; i<source_fe.dofs_per_face; ++i)
    {
      const Point<dim> &p = face_projection.point (i);

      for (unsigned int j=0; j<this->dofs_per_face; ++j)
        {
          double matrix_entry
            = this->shape_value_component (this->face_to_cell_index(j, 0),
                                           p, 0);

          // Correct the interpolated
          // value. I.e. if it is close
          // to 1 or 0, make it exactly
          // 1 or 0. Unfortunately, this
          // is required to avoid problems
          // with higher order elements.
          if ( std::fabs(matrix_entry - 1.0) < eps )
            matrix_entry = 1.0;
          if ( std::fabs(matrix_entry) < eps )
            matrix_entry = 0.0;

          interpolation_matrix(i,j) = matrix_entry;
        }
    }

  // make sure that the row sum of
  // each of the matrices is 1 at
  // this point. this must be so
  // since the shape functions sum up
  // to 1
  for (unsigned int j=0; j<source_fe.dofs_per_face; ++j)
    {
      double sum = 0.;

      for (unsigned int i=0; i<this->dofs_per_face; ++i)
        sum += interpolation_matrix(j,i);

      Assert (std::fabs(sum-1) < 2e-13*this->degree*(dim-1),
              ExcInternalError());
    }
}


template <int dim>
void
FE_RaviartThomasNodal<dim>::get_subface_interpolation_matrix (
  const FiniteElement<dim> &x_source_fe,
  const unsigned int subface,
  FullMatrix<double>       &interpolation_matrix) const
{
  // this is only implemented, if the
  // source FE is also a
  // RaviartThomasNodal element
  AssertThrow ((x_source_fe.get_name().find ("FE_RaviartThomasNodal<") == 0)
               ||
               (dynamic_cast<const FE_RaviartThomasNodal<dim>*>(&x_source_fe) != 0),
               typename FiniteElement<dim>::
               ExcInterpolationNotImplemented());

  Assert (interpolation_matrix.n() == this->dofs_per_face,
          ExcDimensionMismatch (interpolation_matrix.n(),
                                this->dofs_per_face));
  Assert (interpolation_matrix.m() == x_source_fe.dofs_per_face,
          ExcDimensionMismatch (interpolation_matrix.m(),
                                x_source_fe.dofs_per_face));

  // ok, source is a RaviartThomasNodal element, so
  // we will be able to do the work
  const FE_RaviartThomasNodal<dim> &source_fe
    = dynamic_cast<const FE_RaviartThomasNodal<dim>&>(x_source_fe);

  // Make sure, that the element,
  // for which the DoFs should be
  // constrained is the one with
  // the higher polynomial degree.
  // Actually the procedure will work
  // also if this assertion is not
  // satisfied. But the matrices
  // produced in that case might
  // lead to problems in the
  // hp procedures, which use this
  // method.
  Assert (this->dofs_per_face <= source_fe.dofs_per_face,
          typename FiniteElement<dim>::
          ExcInterpolationNotImplemented ());

  // generate a quadrature
  // with the generalized support points.
  // This is later based as a
  // basis for the QProjector,
  // which returns the support
  // points on the face.
  Quadrature<dim-1> quad_face_support (source_fe.get_generalized_face_support_points ());

  // Rule of thumb for FP accuracy,
  // that can be expected for a
  // given polynomial degree.
  // This value is used to cut
  // off values close to zero.
  double eps = 2e-13*this->degree*(dim-1);

  // compute the interpolation
  // matrix by simply taking the
  // value at the support points.

  const Quadrature<dim> subface_projection
    = QProjector<dim>::project_to_subface (quad_face_support, 0, subface);

  for (unsigned int i=0; i<source_fe.dofs_per_face; ++i)
    {
      const Point<dim> &p = subface_projection.point (i);

      for (unsigned int j=0; j<this->dofs_per_face; ++j)
        {
          double matrix_entry
            = this->shape_value_component (this->face_to_cell_index(j, 0), p, 0);

          // Correct the interpolated
          // value. I.e. if it is close
          // to 1 or 0, make it exactly
          // 1 or 0. Unfortunately, this
          // is required to avoid problems
          // with higher order elements.
          if ( std::fabs(matrix_entry - 1.0) < eps )
            matrix_entry = 1.0;
          if ( std::fabs(matrix_entry) < eps )
            matrix_entry = 0.0;

          interpolation_matrix(i,j) = matrix_entry;
        }
    }

  // make sure that the row sum of
  // each of the matrices is 1 at
  // this point. this must be so
  // since the shape functions sum up
  // to 1
  for (unsigned int j=0; j<source_fe.dofs_per_face; ++j)
    {
      double sum = 0.;

      for (unsigned int i=0; i<this->dofs_per_face; ++i)
        sum += interpolation_matrix(j,i);

      Assert (std::fabs(sum-1) < 2e-13*this->degree*(dim-1),
              ExcInternalError());
    }
}



// explicit instantiations
#include "fe_raviart_thomas_nodal.inst"


DEAL_II_NAMESPACE_CLOSE
