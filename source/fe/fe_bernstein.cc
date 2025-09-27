// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2015 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


#include <deal.II/base/polynomials_bernstein.h>
#include <deal.II/base/qprojector.h>
#include <deal.II/base/quadrature.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/fe/fe_bernstein.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_nothing.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_tools.h>

#include <memory>
#include <sstream>
#include <vector>


DEAL_II_NAMESPACE_OPEN



template <int dim, int spacedim>
FE_Bernstein<dim, spacedim>::FE_Bernstein(const unsigned int degree)
  : FE_Q_Base<dim, spacedim>(this->renumber_bases(degree),
                             FiniteElementData<dim>(this->get_dpo_vector(
                                                      degree),
                                                    1,
                                                    degree,
                                                    FiniteElementData<dim>::H1),
                             std::vector<bool>(1, false))
{}



template <int dim, int spacedim>
void
FE_Bernstein<dim, spacedim>::get_interpolation_matrix(
  const FiniteElement<dim, spacedim> &,
  FullMatrix<double> &) const
{
  // no interpolation possible. throw exception, as documentation says
  AssertThrow(
    false,
    (typename FiniteElement<dim, spacedim>::ExcInterpolationNotImplemented()));
}



template <int dim, int spacedim>
const FullMatrix<double> &
FE_Bernstein<dim, spacedim>::get_restriction_matrix(
  const unsigned int,
  const RefinementCase<dim> &) const
{
  AssertThrow(false,
              (typename FiniteElement<dim, spacedim>::ExcProjectionVoid()));
  // return dummy, nothing will happen because the base class FE_Q_Base
  // implements lazy evaluation of those matrices
  return this->restriction[0][0];
}



template <int dim, int spacedim>
const FullMatrix<double> &
FE_Bernstein<dim, spacedim>::get_prolongation_matrix(
  const unsigned int,
  const RefinementCase<dim> &) const
{
  AssertThrow(false,
              (typename FiniteElement<dim, spacedim>::ExcEmbeddingVoid()));
  // return dummy, nothing will happen because the base class FE_Q_Base
  // implements lazy evaluation of those matrices
  return this->prolongation[0][0];
}



template <int dim, int spacedim>
void
FE_Bernstein<dim, spacedim>::get_face_interpolation_matrix(
  const FiniteElement<dim, spacedim> &source_fe,
  FullMatrix<double>                 &interpolation_matrix,
  const unsigned int                  face_no) const
{
  get_subface_interpolation_matrix(source_fe,
                                   numbers::invalid_unsigned_int,
                                   interpolation_matrix,
                                   face_no);
}


template <int dim, int spacedim>
void
FE_Bernstein<dim, spacedim>::get_subface_interpolation_matrix(
  const FiniteElement<dim, spacedim> &x_source_fe,
  const unsigned int                  subface,
  FullMatrix<double>                 &interpolation_matrix,
  const unsigned int                  face_no) const
{
  Assert(interpolation_matrix.m() == x_source_fe.n_dofs_per_face(face_no),
         ExcDimensionMismatch(interpolation_matrix.m(),
                              x_source_fe.n_dofs_per_face(face_no)));

  // see if source is a Bernstein element
  if (const FE_Bernstein<dim, spacedim> *source_fe =
        dynamic_cast<const FE_Bernstein<dim, spacedim> *>(&x_source_fe))
    {
      // have this test in here since a table of size 2x0 reports its size as
      // 0x0
      Assert(interpolation_matrix.n() == this->n_dofs_per_face(face_no),
             ExcDimensionMismatch(interpolation_matrix.n(),
                                  this->n_dofs_per_face(face_no)));

      // Make sure that the element for which the DoFs should be constrained
      // is the one with the higher polynomial degree.  Actually the procedure
      // will work also if this assertion is not satisfied. But the matrices
      // produced in that case might lead to problems in the hp-procedures,
      // which use this method.
      Assert(
        this->n_dofs_per_face(face_no) <= source_fe->n_dofs_per_face(face_no),
        (typename FiniteElement<dim,
                                spacedim>::ExcInterpolationNotImplemented()));

      const Quadrature<dim - 1> quad_face_support(
        FE_Q<dim, spacedim>(QIterated<1>(QTrapezoid<1>(), source_fe->degree))
          .get_unit_face_support_points(face_no));

      // Rule of thumb for FP accuracy, that can be expected for a given
      // polynomial degree.  This value is used to cut off values close to
      // zero.
      const double eps = 2e-13 * std::max(this->degree, source_fe->degree) *
                         std::max(dim - 1, 1);

      // compute the interpolation matrix by simply taking the value at the
      // support points.
      // TODO: Verify that all faces are the same with respect to
      // these support points. Furthermore, check if something has to
      // be done for the face orientation flag in 3d.
      const Quadrature<dim> subface_quadrature =
        subface == numbers::invalid_unsigned_int ?
          QProjector<dim>::project_to_face(
            this->reference_cell(),
            quad_face_support,
            0,
            numbers::default_geometric_orientation) :
          QProjector<dim>::project_to_subface(
            this->reference_cell(),
            quad_face_support,
            0,
            subface,
            numbers::default_geometric_orientation,
            RefinementCase<dim - 1>::isotropic_refinement);

      for (unsigned int i = 0; i < source_fe->n_dofs_per_face(face_no); ++i)
        {
          const Point<dim> &p = subface_quadrature.point(i);
          for (unsigned int j = 0; j < this->n_dofs_per_face(face_no); ++j)
            {
              double matrix_entry =
                this->shape_value(this->face_to_cell_index(j, 0), p);

              // Correct the interpolated value. I.e. if it is close to 1 or
              // 0, make it exactly 1 or 0. Unfortunately, this is required to
              // avoid problems with higher order elements.
              if (std::fabs(matrix_entry - 1.0) < eps)
                matrix_entry = 1.0;
              if (std::fabs(matrix_entry) < eps)
                matrix_entry = 0.0;

              interpolation_matrix(i, j) = matrix_entry;
            }
        }

      if constexpr (running_in_debug_mode())
        {
          // make sure that the row sum of each of the matrices is 1 at this
          // point. this must be so since the shape functions sum up to 1
          for (unsigned int j = 0; j < source_fe->n_dofs_per_face(face_no); ++j)
            {
              double sum = 0.;

              for (unsigned int i = 0; i < this->n_dofs_per_face(face_no); ++i)
                sum += interpolation_matrix(j, i);

              Assert(std::fabs(sum - 1) < eps, ExcInternalError());
            }
        }
    }
  else
    {
      // When the incoming element is not FE_Bernstein we can just delegate to
      // the base class to create the interpolation matrix.
      FE_Q_Base<dim, spacedim>::get_subface_interpolation_matrix(
        x_source_fe, subface, interpolation_matrix, face_no);
    }
}



template <int dim, int spacedim>
bool
FE_Bernstein<dim, spacedim>::hp_constraints_are_implemented() const
{
  return true;
}


template <int dim, int spacedim>
std::vector<std::pair<unsigned int, unsigned int>>
FE_Bernstein<dim, spacedim>::hp_vertex_dof_identities(
  const FiniteElement<dim, spacedim> &fe_other) const
{
  // we can presently only compute these identities if both FEs are
  // FE_Bernsteins or if the other one is an FE_Nothing. in the first case,
  // there should be exactly one single DoF of each FE at a vertex, and they
  // should have identical value
  if (dynamic_cast<const FE_Bernstein<dim, spacedim> *>(&fe_other) != nullptr)
    {
      return std::vector<std::pair<unsigned int, unsigned int>>(
        1, std::make_pair(0U, 0U));
    }
  else if (dynamic_cast<const FE_Nothing<dim> *>(&fe_other) != nullptr)
    {
      // the FE_Nothing has no degrees of freedom, so there are no
      // equivalencies to be recorded
      return std::vector<std::pair<unsigned int, unsigned int>>();
    }
  else if (fe_other.n_unique_faces() == 1 && fe_other.n_dofs_per_face(0) == 0)
    {
      // if the other element has no elements on faces at all,
      // then it would be impossible to enforce any kind of
      // continuity even if we knew exactly what kind of element
      // we have -- simply because the other element declares
      // that it is discontinuous because it has no DoFs on
      // its faces. in that case, just state that we have no
      // constraints to declare
      return std::vector<std::pair<unsigned int, unsigned int>>();
    }
  else
    {
      DEAL_II_NOT_IMPLEMENTED();
      return std::vector<std::pair<unsigned int, unsigned int>>();
    }
}


template <int dim, int spacedim>
std::vector<std::pair<unsigned int, unsigned int>>
FE_Bernstein<dim, spacedim>::hp_line_dof_identities(
  const FiniteElement<dim, spacedim> &) const
{
  // Since this FE is not interpolatory but on the vertices, we can
  // not identify dofs on lines and on quads even if there are dofs
  // on lines and on quads.
  //
  // we also have nothing to say about interpolation to other finite
  // elements. consequently, we never have anything to say at all
  return std::vector<std::pair<unsigned int, unsigned int>>();
}


template <int dim, int spacedim>
std::vector<std::pair<unsigned int, unsigned int>>
FE_Bernstein<dim, spacedim>::hp_quad_dof_identities(
  const FiniteElement<dim, spacedim> &,
  const unsigned int) const
{
  // Since this FE is not interpolatory but on the vertices, we can
  // not identify dofs on lines and on quads even if there are dofs
  // on lines and on quads.
  //
  // we also have nothing to say about interpolation to other finite
  // elements. consequently, we never have anything to say at all
  return std::vector<std::pair<unsigned int, unsigned int>>();
}


template <int dim, int spacedim>
FiniteElementDomination::Domination
FE_Bernstein<dim, spacedim>::compare_for_domination(
  const FiniteElement<dim, spacedim> &fe_other,
  const unsigned int                  codim) const
{
  Assert(codim <= dim, ExcImpossibleInDim(dim));

  // vertex/line/face domination
  // (if fe_other is derived from FE_DGQ)
  // ------------------------------------
  if (codim > 0)
    if (dynamic_cast<const FE_DGQ<dim, spacedim> *>(&fe_other) != nullptr)
      // there are no requirements between continuous and discontinuous elements
      return FiniteElementDomination::no_requirements;

  // vertex/line/face domination
  // (if fe_other is not derived from FE_DGQ)
  // & cell domination
  // ----------------------------------------
  if (const FE_Bernstein<dim, spacedim> *fe_b_other =
        dynamic_cast<const FE_Bernstein<dim, spacedim> *>(&fe_other))
    {
      if (this->degree < fe_b_other->degree)
        return FiniteElementDomination::this_element_dominates;
      else if (this->degree == fe_b_other->degree)
        return FiniteElementDomination::either_element_can_dominate;
      else
        return FiniteElementDomination::other_element_dominates;
    }
  else if (const FE_Nothing<dim, spacedim> *fe_nothing =
             dynamic_cast<const FE_Nothing<dim, spacedim> *>(&fe_other))
    {
      if (fe_nothing->is_dominating())
        return FiniteElementDomination::other_element_dominates;
      else
        // the FE_Nothing has no degrees of freedom and it is typically used
        // in a context where we don't require any continuity along the
        // interface
        return FiniteElementDomination::no_requirements;
    }

  DEAL_II_NOT_IMPLEMENTED();
  return FiniteElementDomination::neither_element_dominates;
}


template <int dim, int spacedim>
std::string
FE_Bernstein<dim, spacedim>::get_name() const
{
  // note that the FETools::get_fe_by_name function depends on the
  // particular format of the string this function returns, so they have to be
  // kept in synch

  std::ostringstream namebuf;
  namebuf << "FE_Bernstein<" << Utilities::dim_string(dim, spacedim) << ">("
          << this->degree << ")";
  return namebuf.str();
}


template <int dim, int spacedim>
std::unique_ptr<FiniteElement<dim, spacedim>>
FE_Bernstein<dim, spacedim>::clone() const
{
  return std::make_unique<FE_Bernstein<dim, spacedim>>(*this);
}


/**
 * Only the assertion differs from the same function in FE_Q_Base!!
 */
template <int dim, int spacedim>
std::vector<unsigned int>
FE_Bernstein<dim, spacedim>::get_dpo_vector(const unsigned int deg)
{
  AssertThrow(deg > 0, ExcMessage("FE_Bernstein needs to be of degree > 0."));
  std::vector<unsigned int> dpo(dim + 1, 1U);
  for (unsigned int i = 1; i < dpo.size(); ++i)
    dpo[i] = dpo[i - 1] * (deg - 1);
  return dpo;
}


template <int dim, int spacedim>
TensorProductPolynomials<dim>
FE_Bernstein<dim, spacedim>::renumber_bases(const unsigned int deg)
{
  TensorProductPolynomials<dim> tpp(
    dealii::generate_complete_bernstein_basis<double>(deg));
  tpp.set_numbering(FETools::hierarchic_to_lexicographic_numbering<dim>(deg));
  return tpp;
}


// explicit instantiations
#include "fe/fe_bernstein.inst"

DEAL_II_NAMESPACE_CLOSE
