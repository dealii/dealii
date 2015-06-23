// ---------------------------------------------------------------------
//
// Copyright (C) 2001 - 2014 by the deal.II authors
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


#include <deal.II/base/quadrature.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/template_constraints.h>
#include <deal.II/fe/fe_dgabf_scalar.h>
#include <deal.II/fe/fe_tools.h>


#include <iostream>
#include <sstream>

DEAL_II_NAMESPACE_OPEN


// namespace for some functions that are used in this file. they are
// specific to numbering conventions used for the FE_DGABF_Scalar element, and
// are thus not very interesting to the outside world
namespace
{
  // given an integer N, compute its
  // integer square root (if it
  // exists, otherwise give up)
  inline unsigned int int_sqrt (const unsigned int N)
  {
    for (unsigned int i=0; i<=N; ++i)
      if (i*i == N)
        return i;
    Assert (false, ExcInternalError());
    return numbers::invalid_unsigned_int;
  }


  // given an integer N, compute its
  // integer cube root (if it
  // exists, otherwise give up)
  inline unsigned int int_cuberoot (const unsigned int N)
  {
    for (unsigned int i=0; i<=N; ++i)
      if (i*i*i == N)
        return i;
    Assert (false, ExcInternalError());
    return numbers::invalid_unsigned_int;
  }


  // given N, generate i=0...N-1
  // equidistant points in the
  // interior of the interval [0,1]
  inline Point<1>
  generate_unit_point (const unsigned int i,
                       const unsigned int N,
                       const dealii::internal::int2type<1>  )
  {
    Assert (i<N, ExcInternalError());
    if (N==1)
      return Point<1> (.5);
    else
      {
        const double h = 1./(N-1);
        return Point<1>(i*h);
      }
  }


  // given N, generate i=0...N-1
  // equidistant points in the domain
  // [0,1]^2
  inline Point<2>
  generate_unit_point (const unsigned int i,
                       const unsigned int N,
                       const dealii::internal::int2type<2>  )
  {
    Assert (i<N, ExcInternalError());

    if (N==1)
      return Point<2> (.5, .5);
    else
      {
        Assert (N>=4, ExcInternalError());
        const unsigned int N1d = int_sqrt(N);
        const double h = 1./(N1d-1);

        return Point<2> (i%N1d * h,
                         i/N1d * h);
      }
  }




  // given N, generate i=0...N-1
  // equidistant points in the domain
  // [0,1]^3
  inline Point<3>
  generate_unit_point (const unsigned int i,
                       const unsigned int N,
                       const dealii::internal::int2type<3>  )
  {
    Assert (i<N, ExcInternalError());
    if (N==1)
      return Point<3> (.5, .5, .5);
    else
      {
        Assert (N>=8, ExcInternalError());

        const unsigned int N1d = int_cuberoot(N);
        const double h = 1./(N1d-1);

        return Point<3> (i%N1d * h,
                         (i/N1d)%N1d * h,
                         i/(N1d*N1d) * h);
      }
  }
}




template <int dim, int spacedim>
FE_DGABF_Scalar<dim, spacedim>::FE_DGABF_Scalar (const unsigned int degree)
  :
  FE_Poly<ABFScalarPolynomials<dim>, dim, spacedim> (
    ABFScalarPolynomials<dim>(Polynomials::Legendre::generate_complete_basis(degree+1)),
    FiniteElementData<dim>(get_dpo_vector(degree), 1, degree, FiniteElementData<dim>::L2),
    std::vector<bool>(FiniteElementData<dim>(get_dpo_vector(degree),1, degree).dofs_per_cell, true),
    std::vector<ComponentMask>(FiniteElementData<dim>(
                                 get_dpo_vector(degree),1, degree).dofs_per_cell, std::vector<bool>(1,true)))
{
  
  // do not initialize embedding and restriction here. these matrices are
  // initialized on demand in get_restriction_matrix and
  // get_prolongation_matrix

  // note: no face support points for DG elements
}


template <int dim, int spacedim>
std::string
FE_DGABF_Scalar<dim, spacedim>::get_name () const
{
  // note that the
  // FETools::get_fe_from_name
  // function depends on the
  // particular format of the string
  // this function returns, so they
  // have to be kept in synch

  std::ostringstream namebuf;
  namebuf << "FE_DGABF_Scalar<"
          << Utilities::dim_string(dim,spacedim)
          << ">(" << this->degree << ")";
  return namebuf.str();
}



template <int dim, int spacedim>
FiniteElement<dim,spacedim> *
FE_DGABF_Scalar<dim, spacedim>::clone() const
{
  return new FE_DGABF_Scalar<dim, spacedim>(*this);
}


//---------------------------------------------------------------------------
// Auxiliary functions
//---------------------------------------------------------------------------


template <int dim, int spacedim>
std::vector<unsigned int>
FE_DGABF_Scalar<dim, spacedim>::get_dpo_vector (const unsigned int deg)
{
  std::vector<unsigned int> dpo(dim+1, 0U);
  dpo[dim] = deg+1;
  for (unsigned int i=1; i<dim; ++i)
    dpo[dim] *= deg+1;

  unsigned int x = 1;
  for(unsigned int i=0; i<dim-1; ++i)
    x *= deg+1;

  dpo[dim] += dim*x;

  return dpo;
}



template <int dim, int spacedim>
void
FE_DGABF_Scalar<dim, spacedim>::rotate_indices (std::vector<unsigned int> &numbers,
                                       const char                 direction) const
{
  Assert (false, ExcNotImplemented ());
}



template <int dim, int spacedim>
void
FE_DGABF_Scalar<dim, spacedim>::
get_interpolation_matrix (const FiniteElement<dim, spacedim> &x_source_fe,
                          FullMatrix<double>           &interpolation_matrix) const
{
  // this is only implemented, if the
  // source FE is also a
  // DGQ element
  typedef FiniteElement<dim, spacedim> FE;
  AssertThrow ((dynamic_cast<const FE_DGABF_Scalar<dim, spacedim>*>(&x_source_fe) != 0),
               typename FE::ExcInterpolationNotImplemented() );

  // ok, source is a Q element, so
  // we will be able to do the work
  const FE_DGABF_Scalar<dim, spacedim> &source_fe
    = dynamic_cast<const FE_DGABF_Scalar<dim, spacedim>&>(x_source_fe);

  Assert (interpolation_matrix.m() == this->dofs_per_cell,
          ExcDimensionMismatch (interpolation_matrix.m(),
                                this->dofs_per_cell));
  Assert (interpolation_matrix.n() == source_fe.dofs_per_cell,
          ExcDimensionMismatch (interpolation_matrix.n(),
                                source_fe.dofs_per_cell));


  // compute the interpolation
  // matrices in much the same way as
  // we do for the embedding matrices
  // from mother to child.
  FullMatrix<double> cell_interpolation (this->dofs_per_cell,
                                         this->dofs_per_cell);
  FullMatrix<double> source_interpolation (this->dofs_per_cell,
                                           source_fe.dofs_per_cell);
  FullMatrix<double> tmp (this->dofs_per_cell,
                          source_fe.dofs_per_cell);
  for (unsigned int j=0; j<this->dofs_per_cell; ++j)
    {
      // generate a point on this
      // cell and evaluate the
      // shape functions there
      const Point<dim> p = generate_unit_point (j, this->dofs_per_cell,
                                                dealii::internal::int2type<dim>());
      for (unsigned int i=0; i<this->dofs_per_cell; ++i)
        cell_interpolation(j,i)
          = this->poly_space.compute_value (i, p);

      for (unsigned int i=0; i<source_fe.dofs_per_cell; ++i)
        source_interpolation(j,i)
          = source_fe.poly_space.compute_value (i, p);
    }

  // then compute the
  // interpolation matrix matrix
  // for this coordinate
  // direction
  cell_interpolation.gauss_jordan ();
  cell_interpolation.mmult (interpolation_matrix,
                            source_interpolation);

  // cut off very small values
  for (unsigned int i=0; i<this->dofs_per_cell; ++i)
    for (unsigned int j=0; j<source_fe.dofs_per_cell; ++j)
      if (std::fabs(interpolation_matrix(i,j)) < 1e-15)
        interpolation_matrix(i,j) = 0.;

  // make sure that the row sum of
  // each of the matrices is 1 at
  // this point. this must be so
  // since the shape functions sum up
  // to 1
  for (unsigned int i=0; i<this->dofs_per_cell; ++i)
    {
      double sum = 0.;
      for (unsigned int j=0; j<source_fe.dofs_per_cell; ++j)
        sum += interpolation_matrix(i,j);

      Assert (std::fabs(sum-1) < 5e-14*std::max(this->degree,1U)*dim,
              ExcInternalError());
    }
}



template <int dim, int spacedim>
void
FE_DGABF_Scalar<dim, spacedim>::
get_face_interpolation_matrix (const FiniteElement<dim, spacedim> &x_source_fe,
                               FullMatrix<double>       &interpolation_matrix) const
{
  // this is only implemented, if the source
  // FE is also a DGQ element. in that case,
  // both elements have no dofs on their
  // faces and the face interpolation matrix
  // is necessarily empty -- i.e. there isn't
  // much we need to do here.
  (void)interpolation_matrix;
  typedef FiniteElement<dim,spacedim> FE;
  AssertThrow ((dynamic_cast<const FE_DGABF_Scalar<dim, spacedim>*>(&x_source_fe) != 0),
               typename FE::ExcInterpolationNotImplemented());

  Assert (interpolation_matrix.m() == 0,
          ExcDimensionMismatch (interpolation_matrix.m(),
                                0));
  Assert (interpolation_matrix.n() == 0,
          ExcDimensionMismatch (interpolation_matrix.m(),
                                0));
}



template <int dim, int spacedim>
void
FE_DGABF_Scalar<dim, spacedim>::
get_subface_interpolation_matrix (const FiniteElement<dim, spacedim> &x_source_fe,
                                  const unsigned int ,
                                  FullMatrix<double>           &interpolation_matrix) const
{
  // this is only implemented, if the source
  // FE is also a DGQ element. in that case,
  // both elements have no dofs on their
  // faces and the face interpolation matrix
  // is necessarily empty -- i.e. there isn't
  // much we need to do here.
  (void)interpolation_matrix;
  typedef FiniteElement<dim, spacedim> FE;
  AssertThrow ((dynamic_cast<const FE_DGABF_Scalar<dim, spacedim>*>(&x_source_fe) != 0),
               typename FE::ExcInterpolationNotImplemented());

  Assert (interpolation_matrix.m() == 0,
          ExcDimensionMismatch (interpolation_matrix.m(),
                                0));
  Assert (interpolation_matrix.n() == 0,
          ExcDimensionMismatch (interpolation_matrix.m(),
                                0));
}



template <int dim, int spacedim>
const FullMatrix<double> &
FE_DGABF_Scalar<dim,spacedim>
::get_prolongation_matrix (const unsigned int child,
                           const RefinementCase<dim> &refinement_case) const
{
  Assert (refinement_case<RefinementCase<dim>::isotropic_refinement+1,
          ExcIndexRange(refinement_case,0,RefinementCase<dim>::isotropic_refinement+1));
  Assert (refinement_case!=RefinementCase<dim>::no_refinement,
          ExcMessage("Prolongation matrices are only available for refined cells!"));
  Assert (child<GeometryInfo<dim>::n_children(refinement_case),
          ExcIndexRange(child,0,GeometryInfo<dim>::n_children(refinement_case)));

  // initialization upon first request
  if (this->prolongation[refinement_case-1][child].n() == 0)
    {
      Threads::Mutex::ScopedLock lock(this->mutex);

      // if matrix got updated while waiting for the lock
      if (this->prolongation[refinement_case-1][child].n() ==
          this->dofs_per_cell)
        return this->prolongation[refinement_case-1][child];

      // now do the work. need to get a non-const version of data in order to
      // be able to modify them inside a const function
      FE_DGABF_Scalar<dim,spacedim> &this_nonconst = const_cast<FE_DGABF_Scalar<dim,spacedim>& >(*this);
      if (refinement_case == RefinementCase<dim>::isotropic_refinement)
        {
          std::vector<std::vector<FullMatrix<double> > >
          isotropic_matrices(RefinementCase<dim>::isotropic_refinement);
          isotropic_matrices.back().
          resize(GeometryInfo<dim>::n_children(RefinementCase<dim>(refinement_case)),
                 FullMatrix<double>(this->dofs_per_cell, this->dofs_per_cell));
          if (dim == spacedim)
            FETools::compute_embedding_matrices (*this, isotropic_matrices, true);
          else
            FETools::compute_embedding_matrices (FE_DGABF_Scalar<dim>(this->degree),
                                                 isotropic_matrices, true);
          this_nonconst.prolongation[refinement_case-1].swap(isotropic_matrices.back());
        }
      else
        {
          // must compute both restriction and prolongation matrices because
          // we only check for their size and the reinit call initializes them
          // all
          this_nonconst.reinit_restriction_and_prolongation_matrices();
          if (dim == spacedim)
            {
              FETools::compute_embedding_matrices (*this, this_nonconst.prolongation);
              FETools::compute_projection_matrices (*this, this_nonconst.restriction);
            }
          else
            {
              FE_DGABF_Scalar<dim> tmp(this->degree);
              FETools::compute_embedding_matrices (tmp, this_nonconst.prolongation);
              FETools::compute_projection_matrices (tmp, this_nonconst.restriction);
            }
        }
    }

  // finally return the matrix
  return this->prolongation[refinement_case-1][child];
}



template <int dim, int spacedim>
const FullMatrix<double> &
FE_DGABF_Scalar<dim,spacedim>
::get_restriction_matrix (const unsigned int child,
                          const RefinementCase<dim> &refinement_case) const
{
  Assert (refinement_case<RefinementCase<dim>::isotropic_refinement+1,
          ExcIndexRange(refinement_case,0,RefinementCase<dim>::isotropic_refinement+1));
  Assert (refinement_case!=RefinementCase<dim>::no_refinement,
          ExcMessage("Restriction matrices are only available for refined cells!"));
  Assert (child<GeometryInfo<dim>::n_children(refinement_case),
          ExcIndexRange(child,0,GeometryInfo<dim>::n_children(refinement_case)));

  // initialization upon first request
  if (this->restriction[refinement_case-1][child].n() == 0)
    {
      Threads::Mutex::ScopedLock lock(this->mutex);

      // if matrix got updated while waiting for the lock...
      if (this->restriction[refinement_case-1][child].n() ==
          this->dofs_per_cell)
        return this->restriction[refinement_case-1][child];

      // now do the work. need to get a non-const version of data in order to
      // be able to modify them inside a const function
      FE_DGABF_Scalar<dim,spacedim> &this_nonconst = const_cast<FE_DGABF_Scalar<dim,spacedim>& >(*this);
      if (refinement_case == RefinementCase<dim>::isotropic_refinement)
        {
          std::vector<std::vector<FullMatrix<double> > >
          isotropic_matrices(RefinementCase<dim>::isotropic_refinement);
          isotropic_matrices.back().
          resize(GeometryInfo<dim>::n_children(RefinementCase<dim>(refinement_case)),
                 FullMatrix<double>(this->dofs_per_cell, this->dofs_per_cell));
          if (dim == spacedim)
            FETools::compute_projection_matrices (*this, isotropic_matrices, true);
          else
            FETools::compute_projection_matrices (FE_DGABF_Scalar<dim>(this->degree),
                                                  isotropic_matrices, true);
          this_nonconst.restriction[refinement_case-1].swap(isotropic_matrices.back());
        }
      else
        {
          // must compute both restriction and prolongation matrices because
          // we only check for their size and the reinit call initializes them
          // all
          this_nonconst.reinit_restriction_and_prolongation_matrices();
          if (dim == spacedim)
            {
              FETools::compute_embedding_matrices (*this, this_nonconst.prolongation);
              FETools::compute_projection_matrices (*this, this_nonconst.restriction);
            }
          else
            {
              FE_DGABF_Scalar<dim> tmp(this->degree);
              FETools::compute_embedding_matrices (tmp, this_nonconst.prolongation);
              FETools::compute_projection_matrices (tmp, this_nonconst.restriction);
            }
        }
    }

  // finally return the matrix
  return this->restriction[refinement_case-1][child];
}



template <int dim, int spacedim>
bool
FE_DGABF_Scalar<dim, spacedim>::hp_constraints_are_implemented () const
{
  return true;
}



template <int dim, int spacedim>
std::vector<std::pair<unsigned int, unsigned int> >
FE_DGABF_Scalar<dim, spacedim>::
hp_vertex_dof_identities (const FiniteElement<dim, spacedim> &/*fe_other*/) const
{
  // this element is discontinuous, so by definition there can
  // be no identities between its dofs and those of any neighbor
  // (of whichever type the neighbor may be -- after all, we have
  // no face dofs on this side to begin with)
  return
    std::vector<std::pair<unsigned int, unsigned int> > ();
}



template <int dim, int spacedim>
std::vector<std::pair<unsigned int, unsigned int> >
FE_DGABF_Scalar<dim, spacedim>::
hp_line_dof_identities (const FiniteElement<dim, spacedim> &/*fe_other*/) const
{
  // this element is discontinuous, so by definition there can
  // be no identities between its dofs and those of any neighbor
  // (of whichever type the neighbor may be -- after all, we have
  // no face dofs on this side to begin with)
  return
    std::vector<std::pair<unsigned int, unsigned int> > ();
}



template <int dim, int spacedim>
std::vector<std::pair<unsigned int, unsigned int> >
FE_DGABF_Scalar<dim, spacedim>::
hp_quad_dof_identities (const FiniteElement<dim, spacedim> &/*fe_other*/) const
{
  // this element is discontinuous, so by definition there can
  // be no identities between its dofs and those of any neighbor
  // (of whichever type the neighbor may be -- after all, we have
  // no face dofs on this side to begin with)
  return
    std::vector<std::pair<unsigned int, unsigned int> > ();
}



template <int dim, int spacedim>
FiniteElementDomination::Domination
FE_DGABF_Scalar<dim, spacedim>::compare_for_face_domination (const FiniteElement<dim, spacedim> &/*fe_other*/) const
{
  // this is a discontinuous element, so by definition there will
  // be no constraints wherever this element comes together with
  // any other kind of element
  return FiniteElementDomination::no_requirements;
}



template <int dim, int spacedim>
bool
FE_DGABF_Scalar<dim, spacedim>::has_support_on_face (const unsigned int shape_index,
                                            const unsigned int face_index) const
{
  Assert (shape_index < this->dofs_per_cell,
          ExcIndexRange (shape_index, 0, this->dofs_per_cell));
  Assert (face_index < GeometryInfo<dim>::faces_per_cell,
          ExcIndexRange (face_index, 0, GeometryInfo<dim>::faces_per_cell));

  // Legendre polynomials
  return true;
}



template <int dim, int spacedim>
std::pair<Table<2,bool>, std::vector<unsigned int> >
FE_DGABF_Scalar<dim,spacedim>::get_constant_modes () const
{
  Assert (false, ExcMessage ("get_constant_modes is called by FE_DGABF_Scalar"));
  Table<2,bool> constant_modes(1, this->dofs_per_cell);
  constant_modes.fill(true);
  return std::pair<Table<2,bool>, std::vector<unsigned int> >
         (constant_modes, std::vector<unsigned int>(1, 0));
}




template <int dim, int spacedim>
std::size_t
FE_DGABF_Scalar<dim, spacedim>::memory_consumption () const
{
  Assert (false, ExcNotImplemented ());
  return 0;
}

// explicit instantiations
#include "fe_dgabf_scalar.inst"


DEAL_II_NAMESPACE_CLOSE
