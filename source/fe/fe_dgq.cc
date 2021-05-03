// ---------------------------------------------------------------------
//
// Copyright (C) 2001 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------


#include <deal.II/base/quadrature.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/fe/fe.h>
#include <deal.II/fe/fe_bernstein.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_nothing.h>
#include <deal.II/fe/fe_pyramid_p.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_q_bubbles.h>
#include <deal.II/fe/fe_q_dg0.h>
#include <deal.II/fe/fe_q_hierarchical.h>
#include <deal.II/fe/fe_q_iso_q1.h>
#include <deal.II/fe/fe_simplex_p.h>
#include <deal.II/fe/fe_simplex_p_bubbles.h>
#include <deal.II/fe/fe_tools.h>
#include <deal.II/fe/fe_wedge_p.h>

#include <deal.II/lac/vector.h>

#include <iostream>
#include <memory>
#include <sstream>


DEAL_II_NAMESPACE_OPEN


namespace internal
{
  namespace FE_DGQ
  {
    namespace
    {
      std::vector<Point<1>>
      get_QGaussLobatto_points(const unsigned int degree)
      {
        if (degree > 0)
          return QGaussLobatto<1>(degree + 1).get_points();
        else
          return std::vector<Point<1>>(1, Point<1>(0.5));
      }
    } // namespace
  }   // namespace FE_DGQ
} // namespace internal



template <int dim, int spacedim>
FE_DGQ<dim, spacedim>::FE_DGQ(const unsigned int degree)
  : FE_Poly<dim, spacedim>(
      TensorProductPolynomials<dim>(
        Polynomials::generate_complete_Lagrange_basis(
          internal::FE_DGQ::get_QGaussLobatto_points(degree))),
      FiniteElementData<dim>(get_dpo_vector(degree),
                             1,
                             degree,
                             FiniteElementData<dim>::L2),
      std::vector<bool>(
        FiniteElementData<dim>(get_dpo_vector(degree), 1, degree)
          .n_dofs_per_cell(),
        true),
      std::vector<ComponentMask>(
        FiniteElementData<dim>(get_dpo_vector(degree), 1, degree)
          .n_dofs_per_cell(),
        std::vector<bool>(1, true)))
{
  // Compute support points, which are the tensor product of the Lagrange
  // interpolation points in the constructor.
  this->unit_support_points =
    Quadrature<dim>(internal::FE_DGQ::get_QGaussLobatto_points(degree))
      .get_points();

  // do not initialize embedding and restriction here. these matrices are
  // initialized on demand in get_restriction_matrix and
  // get_prolongation_matrix

  // note: no face support points for DG elements
}



template <int dim, int spacedim>
FE_DGQ<dim, spacedim>::FE_DGQ(
  const std::vector<Polynomials::Polynomial<double>> &polynomials)
  : FE_Poly<dim, spacedim>(
      TensorProductPolynomials<dim>(polynomials),
      FiniteElementData<dim>(get_dpo_vector(polynomials.size() - 1),
                             1,
                             polynomials.size() - 1,
                             FiniteElementData<dim>::L2),
      std::vector<bool>(FiniteElementData<dim>(get_dpo_vector(
                                                 polynomials.size() - 1),
                                               1,
                                               polynomials.size() - 1)
                          .n_dofs_per_cell(),
                        true),
      std::vector<ComponentMask>(
        FiniteElementData<dim>(get_dpo_vector(polynomials.size() - 1),
                               1,
                               polynomials.size() - 1)
          .n_dofs_per_cell(),
        std::vector<bool>(1, true)))
{
  // No support points can be defined in general. Derived classes might define
  // support points like the class FE_DGQArbitraryNodes

  // do not initialize embedding and restriction here. these matrices are
  // initialized on demand in get_restriction_matrix and
  // get_prolongation_matrix
}



template <int dim, int spacedim>
std::string
FE_DGQ<dim, spacedim>::get_name() const
{
  // note that the FETools::get_fe_by_name function depends on the
  // particular format of the string this function returns, so they have to be
  // kept in sync

  std::ostringstream namebuf;
  namebuf << "FE_DGQ<" << Utilities::dim_string(dim, spacedim) << ">("
          << this->degree << ")";
  return namebuf.str();
}



template <int dim, int spacedim>
void
FE_DGQ<dim, spacedim>::convert_generalized_support_point_values_to_dof_values(
  const std::vector<Vector<double>> &support_point_values,
  std::vector<double> &              nodal_values) const
{
  AssertDimension(support_point_values.size(),
                  this->get_unit_support_points().size());
  AssertDimension(support_point_values.size(), nodal_values.size());
  AssertDimension(this->n_dofs_per_cell(), nodal_values.size());

  for (unsigned int i = 0; i < this->n_dofs_per_cell(); ++i)
    {
      AssertDimension(support_point_values[i].size(), 1);

      nodal_values[i] = support_point_values[i](0);
    }
}



template <int dim, int spacedim>
std::unique_ptr<FiniteElement<dim, spacedim>>
FE_DGQ<dim, spacedim>::clone() const
{
  return std::make_unique<FE_DGQ<dim, spacedim>>(*this);
}


//---------------------------------------------------------------------------
// Auxiliary functions
//---------------------------------------------------------------------------


template <int dim, int spacedim>
std::vector<unsigned int>
FE_DGQ<dim, spacedim>::get_dpo_vector(const unsigned int deg)
{
  std::vector<unsigned int> dpo(dim + 1, 0U);
  dpo[dim] = deg + 1;
  for (unsigned int i = 1; i < dim; ++i)
    dpo[dim] *= deg + 1;
  return dpo;
}



template <int dim, int spacedim>
void
FE_DGQ<dim, spacedim>::rotate_indices(std::vector<unsigned int> &numbers,
                                      const char direction) const
{
  const unsigned int n = this->degree + 1;
  unsigned int       s = n;
  for (unsigned int i = 1; i < dim; ++i)
    s *= n;
  numbers.resize(s);

  unsigned int l = 0;

  if (dim == 1)
    {
      // Mirror around midpoint
      for (unsigned int i = n; i > 0;)
        numbers[l++] = --i;
    }
  else
    {
      switch (direction)
        {
          // Rotate xy-plane
          // counter-clockwise
          case 'z':
            for (unsigned int iz = 0; iz < ((dim > 2) ? n : 1); ++iz)
              for (unsigned int j = 0; j < n; ++j)
                for (unsigned int i = 0; i < n; ++i)
                  {
                    unsigned int k = n * i - j + n - 1 + n * n * iz;
                    numbers[l++]   = k;
                  }
            break;
          // Rotate xy-plane
          // clockwise
          case 'Z':
            for (unsigned int iz = 0; iz < ((dim > 2) ? n : 1); ++iz)
              for (unsigned int iy = 0; iy < n; ++iy)
                for (unsigned int ix = 0; ix < n; ++ix)
                  {
                    unsigned int k = n * ix - iy + n - 1 + n * n * iz;
                    numbers[k]     = l++;
                  }
            break;
          // Rotate yz-plane
          // counter-clockwise
          case 'x':
            Assert(dim > 2, ExcDimensionMismatch(dim, 3));
            for (unsigned int iz = 0; iz < n; ++iz)
              for (unsigned int iy = 0; iy < n; ++iy)
                for (unsigned int ix = 0; ix < n; ++ix)
                  {
                    unsigned int k = n * (n * iy - iz + n - 1) + ix;
                    numbers[l++]   = k;
                  }
            break;
          // Rotate yz-plane
          // clockwise
          case 'X':
            Assert(dim > 2, ExcDimensionMismatch(dim, 3));
            for (unsigned int iz = 0; iz < n; ++iz)
              for (unsigned int iy = 0; iy < n; ++iy)
                for (unsigned int ix = 0; ix < n; ++ix)
                  {
                    unsigned int k = n * (n * iy - iz + n - 1) + ix;
                    numbers[k]     = l++;
                  }
            break;
          default:
            Assert(false, ExcNotImplemented());
        }
    }
}



template <int dim, int spacedim>
void
FE_DGQ<dim, spacedim>::get_interpolation_matrix(
  const FiniteElement<dim, spacedim> &x_source_fe,
  FullMatrix<double> &                interpolation_matrix) const
{
  // this is only implemented, if the
  // source FE is also a
  // DGQ element
  using FE = FiniteElement<dim, spacedim>;
  AssertThrow((dynamic_cast<const FE_DGQ<dim, spacedim> *>(&x_source_fe) !=
               nullptr),
              typename FE::ExcInterpolationNotImplemented());

  // ok, source is a Q element, so
  // we will be able to do the work
  const FE_DGQ<dim, spacedim> &source_fe =
    dynamic_cast<const FE_DGQ<dim, spacedim> &>(x_source_fe);

  Assert(interpolation_matrix.m() == this->n_dofs_per_cell(),
         ExcDimensionMismatch(interpolation_matrix.m(),
                              this->n_dofs_per_cell()));
  Assert(interpolation_matrix.n() == source_fe.n_dofs_per_cell(),
         ExcDimensionMismatch(interpolation_matrix.n(),
                              source_fe.n_dofs_per_cell()));


  // compute the interpolation
  // matrices in much the same way as
  // we do for the embedding matrices
  // from mother to child.
  FullMatrix<double> cell_interpolation(this->n_dofs_per_cell(),
                                        this->n_dofs_per_cell());
  FullMatrix<double> source_interpolation(this->n_dofs_per_cell(),
                                          source_fe.n_dofs_per_cell());
  FullMatrix<double> tmp(this->n_dofs_per_cell(), source_fe.n_dofs_per_cell());
  for (unsigned int j = 0; j < this->n_dofs_per_cell(); ++j)
    {
      // generate a point on this
      // cell and evaluate the
      // shape functions there
      const Point<dim> p = this->unit_support_points[j];
      for (unsigned int i = 0; i < this->n_dofs_per_cell(); ++i)
        cell_interpolation(j, i) = this->poly_space->compute_value(i, p);

      for (unsigned int i = 0; i < source_fe.n_dofs_per_cell(); ++i)
        source_interpolation(j, i) = source_fe.poly_space->compute_value(i, p);
    }

  // then compute the
  // interpolation matrix matrix
  // for this coordinate
  // direction
  cell_interpolation.gauss_jordan();
  cell_interpolation.mmult(interpolation_matrix, source_interpolation);

  // cut off very small values
  for (unsigned int i = 0; i < this->n_dofs_per_cell(); ++i)
    for (unsigned int j = 0; j < source_fe.n_dofs_per_cell(); ++j)
      if (std::fabs(interpolation_matrix(i, j)) < 1e-15)
        interpolation_matrix(i, j) = 0.;

  // make sure that the row sum of
  // each of the matrices is 1 at
  // this point. this must be so
  // since the shape functions sum up
  // to 1
  for (unsigned int i = 0; i < this->n_dofs_per_cell(); ++i)
    {
      double sum = 0.;
      for (unsigned int j = 0; j < source_fe.n_dofs_per_cell(); ++j)
        sum += interpolation_matrix(i, j);

      Assert(std::fabs(sum - 1) < 5e-14 * std::max(this->degree, 1U) * dim,
             ExcInternalError());
    }
}



template <int dim, int spacedim>
void
FE_DGQ<dim, spacedim>::get_face_interpolation_matrix(
  const FiniteElement<dim, spacedim> &x_source_fe,
  FullMatrix<double> &                interpolation_matrix,
  const unsigned int) const
{
  // this is only implemented, if the source
  // FE is also a DGQ element. in that case,
  // both elements have no dofs on their
  // faces and the face interpolation matrix
  // is necessarily empty -- i.e. there isn't
  // much we need to do here.
  (void)interpolation_matrix;
  using FE = FiniteElement<dim, spacedim>;
  AssertThrow((dynamic_cast<const FE_DGQ<dim, spacedim> *>(&x_source_fe) !=
               nullptr),
              typename FE::ExcInterpolationNotImplemented());

  Assert(interpolation_matrix.m() == 0,
         ExcDimensionMismatch(interpolation_matrix.m(), 0));
  Assert(interpolation_matrix.n() == 0,
         ExcDimensionMismatch(interpolation_matrix.m(), 0));
}



template <int dim, int spacedim>
void
FE_DGQ<dim, spacedim>::get_subface_interpolation_matrix(
  const FiniteElement<dim, spacedim> &x_source_fe,
  const unsigned int,
  FullMatrix<double> &interpolation_matrix,
  const unsigned int) const
{
  // this is only implemented, if the source
  // FE is also a DGQ element. in that case,
  // both elements have no dofs on their
  // faces and the face interpolation matrix
  // is necessarily empty -- i.e. there isn't
  // much we need to do here.
  (void)interpolation_matrix;
  using FE = FiniteElement<dim, spacedim>;
  AssertThrow((dynamic_cast<const FE_DGQ<dim, spacedim> *>(&x_source_fe) !=
               nullptr),
              typename FE::ExcInterpolationNotImplemented());

  Assert(interpolation_matrix.m() == 0,
         ExcDimensionMismatch(interpolation_matrix.m(), 0));
  Assert(interpolation_matrix.n() == 0,
         ExcDimensionMismatch(interpolation_matrix.m(), 0));
}



template <int dim, int spacedim>
const FullMatrix<double> &
FE_DGQ<dim, spacedim>::get_prolongation_matrix(
  const unsigned int         child,
  const RefinementCase<dim> &refinement_case) const
{
  AssertIndexRange(refinement_case,
                   RefinementCase<dim>::isotropic_refinement + 1);
  Assert(refinement_case != RefinementCase<dim>::no_refinement,
         ExcMessage(
           "Prolongation matrices are only available for refined cells!"));
  AssertIndexRange(child, GeometryInfo<dim>::n_children(refinement_case));

  // initialization upon first request
  if (this->prolongation[refinement_case - 1][child].n() == 0)
    {
      std::lock_guard<std::mutex> lock(this->mutex);

      // if matrix got updated while waiting for the lock
      if (this->prolongation[refinement_case - 1][child].n() ==
          this->n_dofs_per_cell())
        return this->prolongation[refinement_case - 1][child];

      // now do the work. need to get a non-const version of data in order to
      // be able to modify them inside a const function
      FE_DGQ<dim, spacedim> &this_nonconst =
        const_cast<FE_DGQ<dim, spacedim> &>(*this);
      if (refinement_case == RefinementCase<dim>::isotropic_refinement)
        {
          std::vector<std::vector<FullMatrix<double>>> isotropic_matrices(
            RefinementCase<dim>::isotropic_refinement);
          isotropic_matrices.back().resize(
            GeometryInfo<dim>::n_children(RefinementCase<dim>(refinement_case)),
            FullMatrix<double>(this->n_dofs_per_cell(),
                               this->n_dofs_per_cell()));
          if (dim == spacedim)
            FETools::compute_embedding_matrices(*this,
                                                isotropic_matrices,
                                                true);
          else
            FETools::compute_embedding_matrices(FE_DGQ<dim>(this->degree),
                                                isotropic_matrices,
                                                true);
          this_nonconst.prolongation[refinement_case - 1].swap(
            isotropic_matrices.back());
        }
      else
        {
          // must compute both restriction and prolongation matrices because
          // we only check for their size and the reinit call initializes them
          // all
          this_nonconst.reinit_restriction_and_prolongation_matrices();
          if (dim == spacedim)
            {
              FETools::compute_embedding_matrices(*this,
                                                  this_nonconst.prolongation);
              FETools::compute_projection_matrices(*this,
                                                   this_nonconst.restriction);
            }
          else
            {
              FE_DGQ<dim> tmp(this->degree);
              FETools::compute_embedding_matrices(tmp,
                                                  this_nonconst.prolongation);
              FETools::compute_projection_matrices(tmp,
                                                   this_nonconst.restriction);
            }
        }
    }

  // finally return the matrix
  return this->prolongation[refinement_case - 1][child];
}



template <int dim, int spacedim>
const FullMatrix<double> &
FE_DGQ<dim, spacedim>::get_restriction_matrix(
  const unsigned int         child,
  const RefinementCase<dim> &refinement_case) const
{
  AssertIndexRange(refinement_case,
                   RefinementCase<dim>::isotropic_refinement + 1);
  Assert(refinement_case != RefinementCase<dim>::no_refinement,
         ExcMessage(
           "Restriction matrices are only available for refined cells!"));
  AssertIndexRange(child, GeometryInfo<dim>::n_children(refinement_case));

  // initialization upon first request
  if (this->restriction[refinement_case - 1][child].n() == 0)
    {
      std::lock_guard<std::mutex> lock(this->mutex);

      // if matrix got updated while waiting for the lock...
      if (this->restriction[refinement_case - 1][child].n() ==
          this->n_dofs_per_cell())
        return this->restriction[refinement_case - 1][child];

      // now do the work. need to get a non-const version of data in order to
      // be able to modify them inside a const function
      FE_DGQ<dim, spacedim> &this_nonconst =
        const_cast<FE_DGQ<dim, spacedim> &>(*this);
      if (refinement_case == RefinementCase<dim>::isotropic_refinement)
        {
          std::vector<std::vector<FullMatrix<double>>> isotropic_matrices(
            RefinementCase<dim>::isotropic_refinement);
          isotropic_matrices.back().resize(
            GeometryInfo<dim>::n_children(RefinementCase<dim>(refinement_case)),
            FullMatrix<double>(this->n_dofs_per_cell(),
                               this->n_dofs_per_cell()));
          if (dim == spacedim)
            FETools::compute_projection_matrices(*this,
                                                 isotropic_matrices,
                                                 true);
          else
            FETools::compute_projection_matrices(FE_DGQ<dim>(this->degree),
                                                 isotropic_matrices,
                                                 true);
          this_nonconst.restriction[refinement_case - 1].swap(
            isotropic_matrices.back());
        }
      else
        {
          // must compute both restriction and prolongation matrices because
          // we only check for their size and the reinit call initializes them
          // all
          this_nonconst.reinit_restriction_and_prolongation_matrices();
          if (dim == spacedim)
            {
              FETools::compute_embedding_matrices(*this,
                                                  this_nonconst.prolongation);
              FETools::compute_projection_matrices(*this,
                                                   this_nonconst.restriction);
            }
          else
            {
              FE_DGQ<dim> tmp(this->degree);
              FETools::compute_embedding_matrices(tmp,
                                                  this_nonconst.prolongation);
              FETools::compute_projection_matrices(tmp,
                                                   this_nonconst.restriction);
            }
        }
    }

  // finally return the matrix
  return this->restriction[refinement_case - 1][child];
}



template <int dim, int spacedim>
bool
FE_DGQ<dim, spacedim>::hp_constraints_are_implemented() const
{
  return true;
}



template <int dim, int spacedim>
std::vector<std::pair<unsigned int, unsigned int>>
FE_DGQ<dim, spacedim>::hp_vertex_dof_identities(
  const FiniteElement<dim, spacedim> & /*fe_other*/) const
{
  // this element is discontinuous, so by definition there can
  // be no identities between its dofs and those of any neighbor
  // (of whichever type the neighbor may be -- after all, we have
  // no face dofs on this side to begin with)
  return std::vector<std::pair<unsigned int, unsigned int>>();
}



template <int dim, int spacedim>
std::vector<std::pair<unsigned int, unsigned int>>
FE_DGQ<dim, spacedim>::hp_line_dof_identities(
  const FiniteElement<dim, spacedim> & /*fe_other*/) const
{
  // this element is discontinuous, so by definition there can
  // be no identities between its dofs and those of any neighbor
  // (of whichever type the neighbor may be -- after all, we have
  // no face dofs on this side to begin with)
  return std::vector<std::pair<unsigned int, unsigned int>>();
}



template <int dim, int spacedim>
std::vector<std::pair<unsigned int, unsigned int>>
FE_DGQ<dim, spacedim>::hp_quad_dof_identities(
  const FiniteElement<dim, spacedim> & /*fe_other*/,
  const unsigned int) const
{
  // this element is discontinuous, so by definition there can
  // be no identities between its dofs and those of any neighbor
  // (of whichever type the neighbor may be -- after all, we have
  // no face dofs on this side to begin with)
  return std::vector<std::pair<unsigned int, unsigned int>>();
}



template <int dim, int spacedim>
FiniteElementDomination::Domination
FE_DGQ<dim, spacedim>::compare_for_domination(
  const FiniteElement<dim, spacedim> &fe_other,
  const unsigned int                  codim) const
{
  Assert(codim <= dim, ExcImpossibleInDim(dim));

  // vertex/line/face domination
  // ---------------------------
  if (codim > 0)
    // this is a discontinuous element, so by definition there will
    // be no constraints wherever this element comes together with
    // any other kind of element
    return FiniteElementDomination::no_requirements;

  // cell domination
  // ---------------
  // The following block of conditionals is rather ugly, but there is currently
  // no other way how to deal with a robust comparison of FE_DGQ elements with
  // relevant others in the current implementation.
  if (const FE_DGQ<dim, spacedim> *fe_dgq_other =
        dynamic_cast<const FE_DGQ<dim, spacedim> *>(&fe_other))
    {
      if (this->degree < fe_dgq_other->degree)
        return FiniteElementDomination::this_element_dominates;
      else if (this->degree == fe_dgq_other->degree)
        return FiniteElementDomination::either_element_can_dominate;
      else
        return FiniteElementDomination::other_element_dominates;
    }
  else if (const FE_Q<dim, spacedim> *fe_q_other =
             dynamic_cast<const FE_Q<dim, spacedim> *>(&fe_other))
    {
      if (this->degree < fe_q_other->degree)
        return FiniteElementDomination::this_element_dominates;
      else if (this->degree == fe_q_other->degree)
        return FiniteElementDomination::either_element_can_dominate;
      else
        return FiniteElementDomination::other_element_dominates;
    }
  else if (const FE_Bernstein<dim, spacedim> *fe_bernstein_other =
             dynamic_cast<const FE_Bernstein<dim, spacedim> *>(&fe_other))
    {
      if (this->degree < fe_bernstein_other->degree)
        return FiniteElementDomination::this_element_dominates;
      else if (this->degree == fe_bernstein_other->degree)
        return FiniteElementDomination::either_element_can_dominate;
      else
        return FiniteElementDomination::other_element_dominates;
    }
  else if (const FE_Q_Bubbles<dim, spacedim> *fe_bubbles_other =
             dynamic_cast<const FE_Q_Bubbles<dim, spacedim> *>(&fe_other))
    {
      if (this->degree < fe_bubbles_other->degree)
        return FiniteElementDomination::this_element_dominates;
      else if (this->degree == fe_bubbles_other->degree)
        return FiniteElementDomination::either_element_can_dominate;
      else
        return FiniteElementDomination::other_element_dominates;
    }
  else if (const FE_Q_DG0<dim, spacedim> *fe_dg0_other =
             dynamic_cast<const FE_Q_DG0<dim, spacedim> *>(&fe_other))
    {
      if (this->degree < fe_dg0_other->degree)
        return FiniteElementDomination::this_element_dominates;
      else if (this->degree == fe_dg0_other->degree)
        return FiniteElementDomination::either_element_can_dominate;
      else
        return FiniteElementDomination::other_element_dominates;
    }
  else if (const FE_Q_iso_Q1<dim, spacedim> *fe_q_iso_q1_other =
             dynamic_cast<const FE_Q_iso_Q1<dim, spacedim> *>(&fe_other))
    {
      if (this->degree < fe_q_iso_q1_other->degree)
        return FiniteElementDomination::this_element_dominates;
      else if (this->degree == fe_q_iso_q1_other->degree)
        return FiniteElementDomination::either_element_can_dominate;
      else
        return FiniteElementDomination::other_element_dominates;
    }
  else if (const FE_Q_Hierarchical<dim> *fe_hierarchical_other =
             dynamic_cast<const FE_Q_Hierarchical<dim> *>(&fe_other))
    {
      if (this->degree < fe_hierarchical_other->degree)
        return FiniteElementDomination::this_element_dominates;
      else if (this->degree == fe_hierarchical_other->degree)
        return FiniteElementDomination::either_element_can_dominate;
      else
        return FiniteElementDomination::other_element_dominates;
    }
  else if (const FE_SimplexDGP<dim, spacedim> *fe_dgp_other =
             dynamic_cast<const FE_SimplexDGP<dim, spacedim> *>(&fe_other))
    {
      if (this->degree < fe_dgp_other->degree)
        return FiniteElementDomination::this_element_dominates;
      else if (this->degree == fe_dgp_other->degree)
        return FiniteElementDomination::either_element_can_dominate;
      else
        return FiniteElementDomination::other_element_dominates;
    }
  else if (const FE_Nothing<dim> *fe_nothing =
             dynamic_cast<const FE_Nothing<dim> *>(&fe_other))
    {
      if (fe_nothing->is_dominating())
        return FiniteElementDomination::other_element_dominates;
      else
        // the FE_Nothing has no degrees of freedom and it is typically used
        // in a context where we don't require any continuity along the
        // interface
        return FiniteElementDomination::no_requirements;
    }

  Assert(false, ExcNotImplemented());
  return FiniteElementDomination::neither_element_dominates;
}



template <int dim, int spacedim>
bool
FE_DGQ<dim, spacedim>::has_support_on_face(const unsigned int shape_index,
                                           const unsigned int face_index) const
{
  AssertIndexRange(shape_index, this->n_dofs_per_cell());
  AssertIndexRange(face_index, GeometryInfo<dim>::faces_per_cell);

  unsigned int n = this->degree + 1;

  // For DGQ elements that do not define support points, we cannot define
  // whether they have support at the boundary easily, so return true in any
  // case
  if (this->unit_support_points.empty())
    return true;

  // for DGQ(0) elements or arbitrary node DGQ with support points not located
  // at the element boundary, the single shape functions is constant and
  // therefore lives on the boundary
  bool support_points_on_boundary = true;
  for (unsigned int d = 0; d < dim; ++d)
    if (std::abs(this->unit_support_points[0][d]) > 1e-13)
      support_points_on_boundary = false;
  for (unsigned int d = 0; d < dim; ++d)
    if (std::abs(this->unit_support_points.back()[d] - 1.) > 1e-13)
      support_points_on_boundary = false;
  if (support_points_on_boundary == false)
    return true;

  unsigned int n2 = n * n;

  switch (dim)
    {
      case 1:
        {
          // in 1d, things are simple. since
          // there is only one degree of
          // freedom per vertex in this
          // class, the first is on vertex 0
          // (==face 0 in some sense), the
          // second on face 1:
          return (((shape_index == 0) && (face_index == 0)) ||
                  ((shape_index == this->degree) && (face_index == 1)));
        }

      case 2:
        {
          if (face_index == 0 && (shape_index % n) == 0)
            return true;
          if (face_index == 1 && (shape_index % n) == this->degree)
            return true;
          if (face_index == 2 && shape_index < n)
            return true;
          if (face_index == 3 && shape_index >= this->n_dofs_per_cell() - n)
            return true;
          return false;
        }

      case 3:
        {
          const unsigned int in2 = shape_index % n2;

          // x=0
          if (face_index == 0 && (shape_index % n) == 0)
            return true;
          // x=1
          if (face_index == 1 && (shape_index % n) == n - 1)
            return true;
          // y=0
          if (face_index == 2 && in2 < n)
            return true;
          // y=1
          if (face_index == 3 && in2 >= n2 - n)
            return true;
          // z=0
          if (face_index == 4 && shape_index < n2)
            return true;
          // z=1
          if (face_index == 5 && shape_index >= this->n_dofs_per_cell() - n2)
            return true;
          return false;
        }

      default:
        Assert(false, ExcNotImplemented());
    }
  return true;
}



template <int dim, int spacedim>
std::pair<Table<2, bool>, std::vector<unsigned int>>
FE_DGQ<dim, spacedim>::get_constant_modes() const
{
  Table<2, bool> constant_modes(1, this->n_dofs_per_cell());
  constant_modes.fill(true);
  return std::pair<Table<2, bool>, std::vector<unsigned int>>(
    constant_modes, std::vector<unsigned int>(1, 0));
}



// ------------------------------ FE_DGQArbitraryNodes -----------------------

template <int dim, int spacedim>
FE_DGQArbitraryNodes<dim, spacedim>::FE_DGQArbitraryNodes(
  const Quadrature<1> &points)
  : FE_DGQ<dim, spacedim>(
      Polynomials::generate_complete_Lagrange_basis(points.get_points()))
{
  Assert(points.size() > 0,
         (typename FiniteElement<dim, spacedim>::ExcFEHasNoSupportPoints()));
  this->unit_support_points = Quadrature<dim>(points).get_points();
}



template <int dim, int spacedim>
std::string
FE_DGQArbitraryNodes<dim, spacedim>::get_name() const
{
  // note that the FETools::get_fe_by_name function does not work for
  // FE_DGQArbitraryNodes since there is no initialization by a degree value.
  std::ostringstream  namebuf;
  bool                equidistant = true;
  std::vector<double> points(this->degree + 1);

  auto *const polynomial_space =
    dynamic_cast<TensorProductPolynomials<dim> *>(this->poly_space.get());
  Assert(polynomial_space != nullptr, ExcInternalError());
  std::vector<unsigned int> lexicographic =
    polynomial_space->get_numbering_inverse();
  for (unsigned int j = 0; j <= this->degree; j++)
    points[j] = this->unit_support_points[lexicographic[j]][0];

  // Check whether the support points are equidistant.
  for (unsigned int j = 0; j <= this->degree; j++)
    if (std::abs(points[j] - static_cast<double>(j) / this->degree) > 1e-15)
      {
        equidistant = false;
        break;
      }
  if (this->degree == 0 && std::abs(points[0] - 0.5) < 1e-15)
    equidistant = true;

  if (equidistant == true)
    {
      if (this->degree > 2)
        namebuf << "FE_DGQArbitraryNodes<"
                << Utilities::dim_string(dim, spacedim)
                << ">(QIterated(QTrapezoid()," << this->degree << "))";
      else
        namebuf << "FE_DGQ<" << Utilities::dim_string(dim, spacedim) << ">("
                << this->degree << ")";
      return namebuf.str();
    }

  // Check whether the support points come from QGaussLobatto.
  const QGaussLobatto<1> points_gl(this->degree + 1);
  bool                   gauss_lobatto = true;
  for (unsigned int j = 0; j <= this->degree; j++)
    if (points[j] != points_gl.point(j)(0))
      {
        gauss_lobatto = false;
        break;
      }

  if (gauss_lobatto == true)
    {
      namebuf << "FE_DGQ<" << Utilities::dim_string(dim, spacedim) << ">("
              << this->degree << ")";
      return namebuf.str();
    }

  // Check whether the support points come from QGauss.
  const QGauss<1> points_g(this->degree + 1);
  bool            gauss = true;
  for (unsigned int j = 0; j <= this->degree; j++)
    if (points[j] != points_g.point(j)(0))
      {
        gauss = false;
        break;
      }

  if (gauss == true)
    {
      namebuf << "FE_DGQArbitraryNodes<" << Utilities::dim_string(dim, spacedim)
              << ">(QGauss(" << this->degree + 1 << "))";
      return namebuf.str();
    }

  // Check whether the support points come from QGauss.
  const QGaussLog<1> points_glog(this->degree + 1);
  bool               gauss_log = true;
  for (unsigned int j = 0; j <= this->degree; j++)
    if (points[j] != points_glog.point(j)(0))
      {
        gauss_log = false;
        break;
      }

  if (gauss_log == true)
    {
      namebuf << "FE_DGQArbitraryNodes<" << Utilities::dim_string(dim, spacedim)
              << ">(QGaussLog(" << this->degree + 1 << "))";
      return namebuf.str();
    }

  // All guesses exhausted
  namebuf << "FE_DGQArbitraryNodes<" << Utilities::dim_string(dim, spacedim)
          << ">(QUnknownNodes(" << this->degree + 1 << "))";
  return namebuf.str();
}



template <int dim, int spacedim>
void
FE_DGQArbitraryNodes<dim, spacedim>::
  convert_generalized_support_point_values_to_dof_values(
    const std::vector<Vector<double>> &support_point_values,
    std::vector<double> &              nodal_values) const
{
  AssertDimension(support_point_values.size(),
                  this->get_unit_support_points().size());
  AssertDimension(support_point_values.size(), nodal_values.size());
  AssertDimension(this->n_dofs_per_cell(), nodal_values.size());

  for (unsigned int i = 0; i < this->n_dofs_per_cell(); ++i)
    {
      AssertDimension(support_point_values[i].size(), 1);

      nodal_values[i] = support_point_values[i](0);
    }
}



template <int dim, int spacedim>
std::unique_ptr<FiniteElement<dim, spacedim>>
FE_DGQArbitraryNodes<dim, spacedim>::clone() const
{
  // Construct a dummy quadrature formula containing the FE's nodes:
  std::vector<Point<1>> qpoints(this->degree + 1);
  auto *const           polynomial_space =
    dynamic_cast<TensorProductPolynomials<dim> *>(this->poly_space.get());
  Assert(polynomial_space != nullptr, ExcInternalError());
  std::vector<unsigned int> lexicographic =
    polynomial_space->get_numbering_inverse();
  for (unsigned int i = 0; i <= this->degree; ++i)
    qpoints[i] = Point<1>(this->unit_support_points[lexicographic[i]][0]);
  Quadrature<1> pquadrature(qpoints);

  return std::make_unique<FE_DGQArbitraryNodes<dim, spacedim>>(pquadrature);
}



// ---------------------------------- FE_DGQLegendre -------------------------

template <int dim, int spacedim>
FE_DGQLegendre<dim, spacedim>::FE_DGQLegendre(const unsigned int degree)
  : FE_DGQ<dim, spacedim>(
      Polynomials::Legendre::generate_complete_basis(degree))
{}



template <int dim, int spacedim>
std::pair<Table<2, bool>, std::vector<unsigned int>>
FE_DGQLegendre<dim, spacedim>::get_constant_modes() const
{
  // Legendre represents a constant function by one in the first basis
  // function and zero in all others
  Table<2, bool> constant_modes(1, this->n_dofs_per_cell());
  constant_modes(0, 0) = true;
  return std::pair<Table<2, bool>, std::vector<unsigned int>>(
    constant_modes, std::vector<unsigned int>(1, 0));
}



template <int dim, int spacedim>
std::string
FE_DGQLegendre<dim, spacedim>::get_name() const
{
  return "FE_DGQLegendre<" + Utilities::dim_string(dim, spacedim) + ">(" +
         Utilities::int_to_string(this->degree) + ")";
}



template <int dim, int spacedim>
std::unique_ptr<FiniteElement<dim, spacedim>>
FE_DGQLegendre<dim, spacedim>::clone() const
{
  return std::make_unique<FE_DGQLegendre<dim, spacedim>>(this->degree);
}



// ---------------------------------- FE_DGQHermite --------------------------

template <int dim, int spacedim>
FE_DGQHermite<dim, spacedim>::FE_DGQHermite(const unsigned int degree)
  : FE_DGQ<dim, spacedim>(
      Polynomials::HermiteLikeInterpolation::generate_complete_basis(degree))
{}



template <int dim, int spacedim>
std::string
FE_DGQHermite<dim, spacedim>::get_name() const
{
  return "FE_DGQHermite<" + Utilities::dim_string(dim, spacedim) + ">(" +
         Utilities::int_to_string(this->degree) + ")";
}



template <int dim, int spacedim>
std::unique_ptr<FiniteElement<dim, spacedim>>
FE_DGQHermite<dim, spacedim>::clone() const
{
  return std::make_unique<FE_DGQHermite<dim, spacedim>>(this->degree);
}



// explicit instantiations
#include "fe_dgq.inst"


DEAL_II_NAMESPACE_CLOSE
