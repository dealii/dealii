// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2002 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_nothing.h>
#include <deal.II/fe/fe_q_hierarchical.h>

#include <cmath>
#include <memory>
#include <sstream>

// TODO: implement the adjust_quad_dof_index_for_face_orientation_table and
// adjust_line_dof_index_for_line_orientation_table fields, and write tests
// similar to bits/face_orientation_and_fe_q_*


DEAL_II_NAMESPACE_OPEN

namespace internal
{
  namespace FE_Q_Hierarchical
  {
    namespace
    {
      /**
       * A function which maps  in[i] to i,i.e. output[in[i]] = i;
       */
      inline std::vector<unsigned int>
      invert_numbering(const std::vector<unsigned int> &in)
      {
        std::vector<unsigned int> out(in.size());
        for (unsigned int i = 0; i < in.size(); ++i)
          {
            AssertIndexRange(in[i], out.size());
            out[in[i]] = i;
          }
        return out;
      }
    } // namespace
  }   // namespace FE_Q_Hierarchical
} // namespace internal



template <int dim>
FE_Q_Hierarchical<dim>::FE_Q_Hierarchical(const unsigned int degree)
  : FE_Poly<dim>(TensorProductPolynomials<dim>(
                   Polynomials::Hierarchical::generate_complete_basis(degree)),
                 FiniteElementData<dim>(get_dpo_vector(degree),
                                        1,
                                        degree,
                                        FiniteElementData<dim>::H1),
                 std::vector<bool>(
                   FiniteElementData<dim>(get_dpo_vector(degree), 1, degree)
                     .n_dofs_per_cell(),
                   false),
                 std::vector<ComponentMask>(
                   FiniteElementData<dim>(get_dpo_vector(degree), 1, degree)
                     .n_dofs_per_cell(),
                   ComponentMask(std::vector<bool>(1, true))))
  , face_renumber(face_fe_q_hierarchical_to_hierarchic_numbering(degree))
{
  TensorProductPolynomials<dim> *poly_space_derived_ptr =
    dynamic_cast<TensorProductPolynomials<dim> *>(this->poly_space.get());
  poly_space_derived_ptr->set_numbering(
    hierarchic_to_fe_q_hierarchical_numbering(*this));

  // The matrix @p{dofs_cell} contains the
  // values of the linear functionals of
  // the parent 1d cell applied to the
  // shape functions of the two 1d subcells.
  // The matrix @p{dofs_subcell} contains
  // the values of the linear functionals
  // on each 1d subcell applied to the
  // shape functions on the parent 1d
  // subcell.
  // We use @p{dofs_cell} and
  // @p{dofs_subcell} to compute the
  // @p{prolongation}, @p{restriction} and
  // @p{interface_constraints} matrices
  // for all dimensions.
  std::vector<FullMatrix<double>> dofs_cell(
    GeometryInfo<1>::max_children_per_cell,
    FullMatrix<double>(2 * this->n_dofs_per_vertex() + this->n_dofs_per_line(),
                       2 * this->n_dofs_per_vertex() +
                         this->n_dofs_per_line()));
  std::vector<FullMatrix<double>> dofs_subcell(
    GeometryInfo<1>::max_children_per_cell,
    FullMatrix<double>(2 * this->n_dofs_per_vertex() + this->n_dofs_per_line(),
                       2 * this->n_dofs_per_vertex() +
                         this->n_dofs_per_line()));
  // build these fields, as they are
  // needed as auxiliary fields later
  // on
  build_dofs_cell(dofs_cell, dofs_subcell);

  // then use them to initialize
  // other fields
  initialize_constraints(dofs_subcell);
  initialize_embedding_and_restriction(dofs_cell, dofs_subcell);

  // finally fill in support points
  // on cell and face
  initialize_generalized_support_points();
  initialize_generalized_face_support_points();
}



template <int dim>
std::string
FE_Q_Hierarchical<dim>::get_name() const
{
  // note that the
  // FETools::get_fe_by_name
  // function depends on the
  // particular format of the string
  // this function returns, so they
  // have to be kept in synch

  std::ostringstream namebuf;
  namebuf << "FE_Q_Hierarchical<" << dim << ">(" << this->degree << ")";

  return namebuf.str();
}



template <int dim>
std::unique_ptr<FiniteElement<dim, dim>>
FE_Q_Hierarchical<dim>::clone() const
{
  return std::make_unique<FE_Q_Hierarchical<dim>>(*this);
}



template <int dim>
void
FE_Q_Hierarchical<dim>::get_interpolation_matrix(
  const FiniteElement<dim> &source,
  FullMatrix<double>       &matrix) const
{
  // support interpolation between FE_Q_Hierarchical only.
  if (const FE_Q_Hierarchical<dim> *source_fe =
        dynamic_cast<const FE_Q_Hierarchical<dim> *>(&source))
    {
      // ok, source is a Q_Hierarchical element, so we will be able to do the
      // work
      Assert(matrix.m() == this->n_dofs_per_cell(),
             ExcDimensionMismatch(matrix.m(), this->n_dofs_per_cell()));
      Assert(matrix.n() == source.n_dofs_per_cell(),
             ExcDimensionMismatch(matrix.m(), source_fe->n_dofs_per_cell()));

      // Recall that DoFs are renumbered in the following order:
      // vertices, lines, quads, hexes.
      // As we deal with hierarchical FE, interpolation matrix is rather easy:
      // it has 1 on pairs of dofs which are the same.
      // To get those use get_embedding_dofs();

      matrix = 0.;

      // distinguish between the case when we interpolate to a richer element
      if (this->n_dofs_per_cell() >= source_fe->n_dofs_per_cell())
        {
          const std::vector<unsigned int> dof_map =
            this->get_embedding_dofs(source_fe->degree);
          for (unsigned int j = 0; j < dof_map.size(); ++j)
            matrix[dof_map[j]][j] = 1.;
        }
      // and when just truncate higher modes.
      else
        {
          const std::vector<unsigned int> dof_map =
            source_fe->get_embedding_dofs(this->degree);
          for (unsigned int j = 0; j < dof_map.size(); ++j)
            matrix[j][dof_map[j]] = 1.;
        }
    }
  else
    {
      AssertThrow(
        false,
        dealii::ExcMessage(
          "Interpolation is supported only between FE_Q_Hierarchical"));
    }
}

template <int dim>
const FullMatrix<double> &
FE_Q_Hierarchical<dim>::get_prolongation_matrix(
  const unsigned int         child,
  const RefinementCase<dim> &refinement_case) const
{
  Assert(
    refinement_case == RefinementCase<dim>::isotropic_refinement,
    ExcMessage(
      "Prolongation matrices are only available for isotropic refinement!"));

  AssertIndexRange(
    child, this->reference_cell().template n_children<dim>(refinement_case));

  return this->prolongation[refinement_case - 1][child];
}


template <int dim>
bool
FE_Q_Hierarchical<dim>::hp_constraints_are_implemented() const
{
  return true;
}


template <int dim>
std::vector<std::pair<unsigned int, unsigned int>>
FE_Q_Hierarchical<dim>::hp_vertex_dof_identities(
  const FiniteElement<dim> &fe_other) const
{
  // we can presently only compute
  // these identities if both FEs are
  // FE_Q_Hierarchicals or if the other
  // one is an FE_Nothing. in the first
  // case, there should be exactly one
  // single DoF of each FE at a
  // vertex, and they should have
  // identical value
  if (dynamic_cast<const FE_Q_Hierarchical<dim> *>(&fe_other) != nullptr)
    {
      return std::vector<std::pair<unsigned int, unsigned int>>(
        1, std::make_pair(0U, 0U));
    }
  else if (dynamic_cast<const FE_Nothing<dim> *>(&fe_other) != nullptr)
    {
      // the FE_Nothing has no
      // degrees of freedom, so there
      // are no equivalencies to be
      // recorded
      return std::vector<std::pair<unsigned int, unsigned int>>();
    }
  else
    {
      DEAL_II_NOT_IMPLEMENTED();
      return std::vector<std::pair<unsigned int, unsigned int>>();
    }
}

template <int dim>
std::vector<std::pair<unsigned int, unsigned int>>
FE_Q_Hierarchical<dim>::hp_line_dof_identities(
  const FiniteElement<dim> &fe_other) const
{
  // we can presently only compute
  // these identities if both FEs are
  // FE_Q_Hierarchicals or if the other
  // one is an FE_Nothing.
  if (dynamic_cast<const FE_Q_Hierarchical<dim> *>(&fe_other) != nullptr)
    {
      const unsigned int this_dpl  = this->n_dofs_per_line();
      const unsigned int other_dpl = fe_other.n_dofs_per_line();

      // we deal with hierarchical 1d polynomials where dofs are enumerated
      // increasingly. Thus we return a vector of pairs for the first N-1, where
      // N is minimum number of dofs_per_line for each FE_Q_Hierarchical.
      std::vector<std::pair<unsigned int, unsigned int>> res;
      res.reserve(std::min(this_dpl, other_dpl));
      for (unsigned int i = 0; i < std::min(this_dpl, other_dpl); ++i)
        res.emplace_back(i, i);

      return res;
    }
  else if (dynamic_cast<const FE_Nothing<dim> *>(&fe_other) != nullptr)
    {
      // the FE_Nothing has no
      // degrees of freedom, so there
      // are no equivalencies to be
      // recorded
      return std::vector<std::pair<unsigned int, unsigned int>>();
    }
  else
    {
      DEAL_II_NOT_IMPLEMENTED();
      return std::vector<std::pair<unsigned int, unsigned int>>();
    }
}

template <int dim>
std::vector<std::pair<unsigned int, unsigned int>>
FE_Q_Hierarchical<dim>::hp_quad_dof_identities(
  const FiniteElement<dim> &fe_other,
  const unsigned int) const
{
  // we can presently only compute
  // these identities if both FEs are
  // FE_Q_Hierarchicals or if the other
  // one is an FE_Nothing.
  if (dynamic_cast<const FE_Q_Hierarchical<dim> *>(&fe_other) != nullptr)
    {
      // TODO: the implementation makes the assumption that all faces have the
      // same number of dofs
      AssertDimension(this->n_unique_faces(), 1);
      const unsigned int face_no = 0;

      const unsigned int this_dpq  = this->n_dofs_per_quad(face_no);
      const unsigned int other_dpq = fe_other.n_dofs_per_quad(face_no);

      // we deal with hierarchical 1d polynomials where dofs are enumerated
      // increasingly. Thus we return a vector of pairs for the first N-1, where
      // N is minimum number of dofs_per_line for each FE_Q_Hierarchical.
      std::vector<std::pair<unsigned int, unsigned int>> res;
      res.reserve(std::min(this_dpq, other_dpq));
      for (unsigned int i = 0; i < std::min(this_dpq, other_dpq); ++i)
        res.emplace_back(i, i);

      return res;
    }
  else if (dynamic_cast<const FE_Nothing<dim> *>(&fe_other) != nullptr)
    {
      // the FE_Nothing has no
      // degrees of freedom, so there
      // are no equivalencies to be
      // recorded
      return std::vector<std::pair<unsigned int, unsigned int>>();
    }
  else
    {
      DEAL_II_NOT_IMPLEMENTED();
      return std::vector<std::pair<unsigned int, unsigned int>>();
    }
}

template <int dim>
FiniteElementDomination::Domination
FE_Q_Hierarchical<dim>::compare_for_domination(
  const FiniteElement<dim> &fe_other,
  const unsigned int        codim) const
{
  Assert(codim <= dim, ExcImpossibleInDim(dim));

  // vertex/line/face domination
  // (if fe_other is derived from FE_DGQ)
  // ------------------------------------
  if (codim > 0)
    if (dynamic_cast<const FE_DGQ<dim> *>(&fe_other) != nullptr)
      // there are no requirements between continuous and discontinuous elements
      return FiniteElementDomination::no_requirements;

  // vertex/line/face domination
  // (if fe_other is not derived from FE_DGQ)
  // & cell domination
  // ----------------------------------------
  if (const FE_Q_Hierarchical<dim> *fe_hierarchical_other =
        dynamic_cast<const FE_Q_Hierarchical<dim> *>(&fe_other))
    {
      if (this->degree < fe_hierarchical_other->degree)
        return FiniteElementDomination::this_element_dominates;
      else if (this->degree == fe_hierarchical_other->degree)
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

  DEAL_II_NOT_IMPLEMENTED();
  return FiniteElementDomination::neither_element_dominates;
}


//---------------------------------------------------------------------------
// Auxiliary functions
//---------------------------------------------------------------------------


template <int dim>
void
FE_Q_Hierarchical<dim>::build_dofs_cell(
  std::vector<FullMatrix<double>> &dofs_cell,
  std::vector<FullMatrix<double>> &dofs_subcell) const
{
  const unsigned int dofs_1d =
    2 * this->n_dofs_per_vertex() + this->n_dofs_per_line();

  // The dofs_subcell matrices are transposed
  // (4.19), (4.21) and (4.27),(4.28),(4.30) in
  // Demkowicz, Oden, Rachowicz, Hardy, CMAMAE 77, 79-112, 1989
  // so that
  // DoFs_c(j) = dofs_subcell[c](j,k) dofs_cell(k)

  // TODO: The dofs_subcell shall differ by a factor 2^p due to shape functions
  // defined on [0,1] instead of [-1,1]. However, that does not seem to be
  // the case. Perhaps this factor is added later on in auxiliary functions
  // which use these matrices.

  // dofs_cells[0](j,k):
  //    0  1 |  2  3  4.
  // 0  1  0 |         .
  // 1  0  0 |         .
  // -------------------
  // 2          \      .
  // 3            \  2^k * k! / (k-j)!j!
  // 4               \ .

  // dofs_subcells[0](j,k):
  //    0    1   |  2  3   4  5  6 .
  // 0  1    0   |                 .
  // 1  1/2  1/2 | -1  0  -1  0  -1.
  // -------------------------------
  // 2              \              .
  // 3                 \     (-1)^(k+j)/ 2^k * k!/(k-j)!j!
  // 4                     \       .

  // dofs_cells[1](j,k):
  //    0  1 |  2  3  4.
  // 0  0  0 |         .
  // 1  0  1 |         .
  // -------------------
  // 2          \      .
  // 3             \   (-1)^(k+j) * 2^k * k!/(k-j)!j!
  // 4                \.

  // dofs_subcells[1](j,k):
  //    0    1   |  2  3   4  5  6 .
  // 0  1/2  1/2 | -1  0  -1  0  -1.
  // 1  0    1   |                 .
  // -------------------------------
  // 2              \              .
  // 3                 \      1/ 2^k * k!/(k-j)!j!
  // 4                             .

  for (unsigned int c = 0; c < GeometryInfo<1>::max_children_per_cell; ++c)
    for (unsigned int j = 0; j < dofs_1d; ++j)
      for (unsigned int k = 0; k < dofs_1d; ++k)
        {
          // upper diagonal block
          if ((j <= 1) && (k <= 1))
            {
              if (((c == 0) && (j == 0) && (k == 0)) ||
                  ((c == 1) && (j == 1) && (k == 1)))
                dofs_cell[c](j, k) = 1.;
              else
                dofs_cell[c](j, k) = 0.;

              if (((c == 0) && (j == 1)) || ((c == 1) && (j == 0)))
                dofs_subcell[c](j, k) = .5;
              else if (((c == 0) && (k == 0)) || ((c == 1) && (k == 1)))
                dofs_subcell[c](j, k) = 1.;
              else
                dofs_subcell[c](j, k) = 0.;
            }
          // upper right block
          else if ((j <= 1) && (k >= 2))
            {
              if (((c == 0) && (j == 1) && ((k % 2) == 0)) ||
                  ((c == 1) && (j == 0) && ((k % 2) == 0)))
                dofs_subcell[c](j, k) = -1.;
            }
          // upper diagonal block
          else if ((j >= 2) && (k >= 2) && (j <= k))
            {
              double factor = 1.;
              for (unsigned int i = 1; i <= j; ++i)
                factor *=
                  (static_cast<double>(k - i + 1)) / (static_cast<double>(i));
              // factor == k * (k-1) * ... * (k-j+1) / j! = k! / (k-j)! / j!
              if (c == 0)
                {
                  dofs_subcell[c](j, k) = ((k + j) % 2 == 0) ?
                                            Utilities::pow(.5, k) * factor :
                                            -Utilities::pow(.5, k) * factor;
                  dofs_cell[c](j, k)    = Utilities::pow(2, j) * factor;
                }
              else
                {
                  dofs_subcell[c](j, k) = Utilities::pow(.5, k) * factor;
                  dofs_cell[c](j, k)    = ((k + j) % 2 == 0) ?
                                            Utilities::pow(2, j) * factor :
                                            -Utilities::pow(2, j) * factor;
                }
            }
        }
}



template <int dim>
void
FE_Q_Hierarchical<dim>::initialize_constraints(
  const std::vector<FullMatrix<double>> &dofs_subcell)
{
  const unsigned int dofs_1d =
    2 * this->n_dofs_per_vertex() + this->n_dofs_per_line();

  this->interface_constraints.TableBase<2, double>::reinit(
    this->interface_constraints_size());

  switch (dim)
    {
      case 1:
        {
          // no constraints in 1d
          break;
        }

      case 2:
        {
          // vertex node
          for (unsigned int i = 0; i < dofs_1d; ++i)
            this->interface_constraints(0, i) = dofs_subcell[0](1, i);
          // edge nodes
          for (unsigned int c = 0; c < GeometryInfo<1>::max_children_per_cell;
               ++c)
            for (unsigned int i = 0; i < dofs_1d; ++i)
              for (unsigned int j = 2; j < dofs_1d; ++j)
                this->interface_constraints(1 + c * (this->degree - 1) + j - 2,
                                            i) = dofs_subcell[c](j, i);
          break;
        }

      case 3:
        {
          for (unsigned int i = 0; i < dofs_1d * dofs_1d; ++i)
            {
              // center vertex node
              this->interface_constraints(0, face_renumber[i]) =
                dofs_subcell[0](1, i % dofs_1d) *
                dofs_subcell[0](1, (i - (i % dofs_1d)) / dofs_1d);

              // boundary vertex nodes
              this->interface_constraints(1, face_renumber[i]) =
                dofs_subcell[0](0, i % dofs_1d) *
                dofs_subcell[0](1, (i - (i % dofs_1d)) / dofs_1d);
              this->interface_constraints(2, face_renumber[i]) =
                dofs_subcell[1](1, i % dofs_1d) *
                dofs_subcell[0](1, (i - (i % dofs_1d)) / dofs_1d);
              this->interface_constraints(3, face_renumber[i]) =
                dofs_subcell[0](1, i % dofs_1d) *
                dofs_subcell[0](0, (i - (i % dofs_1d)) / dofs_1d);
              this->interface_constraints(4, face_renumber[i]) =
                dofs_subcell[1](0, i % dofs_1d) *
                dofs_subcell[1](1, (i - (i % dofs_1d)) / dofs_1d);

              // interior edges
              for (unsigned int j = 0; j < (this->degree - 1); ++j)
                {
                  this->interface_constraints(5 + j, face_renumber[i]) =
                    dofs_subcell[0](1, i % dofs_1d) *
                    dofs_subcell[0](2 + j, (i - (i % dofs_1d)) / dofs_1d);
                  this->interface_constraints(5 + (this->degree - 1) + j,
                                              face_renumber[i]) =
                    dofs_subcell[0](1, i % dofs_1d) *
                    dofs_subcell[1](2 + j, (i - (i % dofs_1d)) / dofs_1d);
                  this->interface_constraints(5 + 2 * (this->degree - 1) + j,
                                              face_renumber[i]) =
                    dofs_subcell[0](2 + j, i % dofs_1d) *
                    dofs_subcell[1](0, (i - (i % dofs_1d)) / dofs_1d);
                  this->interface_constraints(5 + 3 * (this->degree - 1) + j,
                                              face_renumber[i]) =
                    dofs_subcell[1](2 + j, i % dofs_1d) *
                    dofs_subcell[0](1, (i - (i % dofs_1d)) / dofs_1d);
                }

              // boundary edges
              for (unsigned int j = 0; j < (this->degree - 1); ++j)
                {
                  // left edge
                  this->interface_constraints(5 + 4 * (this->degree - 1) + j,
                                              face_renumber[i]) =
                    dofs_subcell[0](0, i % dofs_1d) *
                    dofs_subcell[0](2 + j, (i - (i % dofs_1d)) / dofs_1d);
                  this->interface_constraints(5 + 4 * (this->degree - 1) +
                                                (this->degree - 1) + j,
                                              face_renumber[i]) =
                    dofs_subcell[0](0, i % dofs_1d) *
                    dofs_subcell[1](2 + j, (i - (i % dofs_1d)) / dofs_1d);
                  // right edge
                  this->interface_constraints(5 + 4 * (this->degree - 1) +
                                                2 * (this->degree - 1) + j,
                                              face_renumber[i]) =
                    dofs_subcell[1](1, i % dofs_1d) *
                    dofs_subcell[0](2 + j, (i - (i % dofs_1d)) / dofs_1d);
                  this->interface_constraints(5 + 4 * (this->degree - 1) +
                                                3 * (this->degree - 1) + j,
                                              face_renumber[i]) =
                    dofs_subcell[1](1, i % dofs_1d) *
                    dofs_subcell[1](2 + j, (i - (i % dofs_1d)) / dofs_1d);
                  // bottom edge
                  this->interface_constraints(5 + 4 * (this->degree - 1) +
                                                4 * (this->degree - 1) + j,
                                              face_renumber[i]) =
                    dofs_subcell[0](2 + j, i % dofs_1d) *
                    dofs_subcell[0](0, (i - (i % dofs_1d)) / dofs_1d);
                  this->interface_constraints(5 + 4 * (this->degree - 1) +
                                                5 * (this->degree - 1) + j,
                                              face_renumber[i]) =
                    dofs_subcell[1](2 + j, i % dofs_1d) *
                    dofs_subcell[0](0, (i - (i % dofs_1d)) / dofs_1d);
                  // top edge
                  this->interface_constraints(5 + 4 * (this->degree - 1) +
                                                6 * (this->degree - 1) + j,
                                              face_renumber[i]) =
                    dofs_subcell[0](2 + j, i % dofs_1d) *
                    dofs_subcell[1](1, (i - (i % dofs_1d)) / dofs_1d);
                  this->interface_constraints(5 + 4 * (this->degree - 1) +
                                                7 * (this->degree - 1) + j,
                                              face_renumber[i]) =
                    dofs_subcell[1](2 + j, i % dofs_1d) *
                    dofs_subcell[1](1, (i - (i % dofs_1d)) / dofs_1d);
                }

              // interior faces
              for (unsigned int j = 0; j < (this->degree - 1); ++j)
                for (unsigned int k = 0; k < (this->degree - 1); ++k)
                  {
                    // subcell 0
                    this->interface_constraints(5 + 12 * (this->degree - 1) +
                                                  j + k * (this->degree - 1),
                                                face_renumber[i]) =
                      dofs_subcell[0](2 + j, i % dofs_1d) *
                      dofs_subcell[0](2 + k, (i - (i % dofs_1d)) / dofs_1d);
                    // subcell 1
                    this->interface_constraints(5 + 12 * (this->degree - 1) +
                                                  j + k * (this->degree - 1) +
                                                  (this->degree - 1) *
                                                    (this->degree - 1),
                                                face_renumber[i]) =
                      dofs_subcell[1](2 + j, i % dofs_1d) *
                      dofs_subcell[0](2 + k, (i - (i % dofs_1d)) / dofs_1d);
                    // subcell 2
                    this->interface_constraints(5 + 12 * (this->degree - 1) +
                                                  j + k * (this->degree - 1) +
                                                  2 * (this->degree - 1) *
                                                    (this->degree - 1),
                                                face_renumber[i]) =
                      dofs_subcell[0](2 + j, i % dofs_1d) *
                      dofs_subcell[1](2 + k, (i - (i % dofs_1d)) / dofs_1d);
                    // subcell 3
                    this->interface_constraints(5 + 12 * (this->degree - 1) +
                                                  j + k * (this->degree - 1) +
                                                  3 * (this->degree - 1) *
                                                    (this->degree - 1),
                                                face_renumber[i]) =
                      dofs_subcell[1](2 + j, i % dofs_1d) *
                      dofs_subcell[1](2 + k, (i - (i % dofs_1d)) / dofs_1d);
                  }
            }
          break;
        }

      default:
        DEAL_II_NOT_IMPLEMENTED();
    }
}



template <int dim>
void
FE_Q_Hierarchical<dim>::initialize_embedding_and_restriction(
  const std::vector<FullMatrix<double>> &dofs_cell,
  const std::vector<FullMatrix<double>> &dofs_subcell)
{
  unsigned int iso = RefinementCase<dim>::isotropic_refinement - 1;

  const unsigned int dofs_1d =
    2 * this->n_dofs_per_vertex() + this->n_dofs_per_line();
  TensorProductPolynomials<dim> *poly_space_derived_ptr =
    dynamic_cast<TensorProductPolynomials<dim> *>(this->poly_space.get());
  const std::vector<unsigned int> &renumber =
    poly_space_derived_ptr->get_numbering();

  for (unsigned int c = 0; c < GeometryInfo<dim>::max_children_per_cell; ++c)
    {
      this->prolongation[iso][c].reinit(this->n_dofs_per_cell(),
                                        this->n_dofs_per_cell());
      this->restriction[iso][c].reinit(this->n_dofs_per_cell(),
                                       this->n_dofs_per_cell());
    }

  // the 1d case is particularly
  // simple, so special case it:
  if (dim == 1)
    {
      for (unsigned int c = 0; c < GeometryInfo<dim>::max_children_per_cell;
           ++c)
        {
          this->prolongation[iso][c].fill(dofs_subcell[c]);
          this->restriction[iso][c].fill(dofs_cell[c]);
        }
      return;
    }

  // for higher dimensions, things
  // are a little more tricky:

  // j loops over dofs in the
  // subcell.  These are the rows in
  // the embedding matrix.
  //
  // i loops over the dofs in the
  // parent cell. These are the
  // columns in the embedding matrix.
  for (unsigned int j = 0; j < this->n_dofs_per_cell(); ++j)
    for (unsigned int i = 0; i < this->n_dofs_per_cell(); ++i)
      switch (dim)
        {
          case 2:
            {
              for (unsigned int c = 0;
                   c < GeometryInfo<2>::max_children_per_cell;
                   ++c)
                {
                  // left/right line: 0/1
                  const unsigned int c0 = c % 2;
                  // bottom/top line: 0/1
                  const unsigned int c1 = c / 2;

                  this->prolongation[iso][c](j, i) =
                    dofs_subcell[c0](renumber[j] % dofs_1d,
                                     renumber[i] % dofs_1d) *
                    dofs_subcell[c1]((renumber[j] - (renumber[j] % dofs_1d)) /
                                       dofs_1d,
                                     (renumber[i] - (renumber[i] % dofs_1d)) /
                                       dofs_1d);

                  this->restriction[iso][c](j, i) =
                    dofs_cell[c0](renumber[j] % dofs_1d,
                                  renumber[i] % dofs_1d) *
                    dofs_cell[c1]((renumber[j] - (renumber[j] % dofs_1d)) /
                                    dofs_1d,
                                  (renumber[i] - (renumber[i] % dofs_1d)) /
                                    dofs_1d);
                }
              break;
            }

          case 3:
            {
              for (unsigned int c = 0;
                   c < GeometryInfo<3>::max_children_per_cell;
                   ++c)
                {
                  // left/right face: 0/1
                  const unsigned int c0 = c % 2;
                  // front/back face: 0/1
                  const unsigned int c1 = (c % 4) / 2;
                  // bottom/top face: 0/1
                  const unsigned int c2 = c / 4;

                  this->prolongation[iso][c](j, i) =
                    dofs_subcell[c0](renumber[j] % dofs_1d,
                                     renumber[i] % dofs_1d) *
                    dofs_subcell[c1](
                      ((renumber[j] - (renumber[j] % dofs_1d)) / dofs_1d) %
                        dofs_1d,
                      ((renumber[i] - (renumber[i] % dofs_1d)) / dofs_1d) %
                        dofs_1d) *
                    dofs_subcell[c2](
                      ((renumber[j] - (renumber[j] % dofs_1d)) / dofs_1d -
                       (((renumber[j] - (renumber[j] % dofs_1d)) / dofs_1d) %
                        dofs_1d)) /
                        dofs_1d,
                      ((renumber[i] - (renumber[i] % dofs_1d)) / dofs_1d -
                       (((renumber[i] - (renumber[i] % dofs_1d)) / dofs_1d) %
                        dofs_1d)) /
                        dofs_1d);

                  this->restriction[iso][c](j, i) =
                    dofs_cell[c0](renumber[j] % dofs_1d,
                                  renumber[i] % dofs_1d) *
                    dofs_cell[c1](
                      ((renumber[j] - (renumber[j] % dofs_1d)) / dofs_1d) %
                        dofs_1d,
                      ((renumber[i] - (renumber[i] % dofs_1d)) / dofs_1d) %
                        dofs_1d) *
                    dofs_cell[c2](
                      ((renumber[j] - (renumber[j] % dofs_1d)) / dofs_1d -
                       (((renumber[j] - (renumber[j] % dofs_1d)) / dofs_1d) %
                        dofs_1d)) /
                        dofs_1d,
                      ((renumber[i] - (renumber[i] % dofs_1d)) / dofs_1d -
                       (((renumber[i] - (renumber[i] % dofs_1d)) / dofs_1d) %
                        dofs_1d)) /
                        dofs_1d);
                }
              break;
            }

          default:
            DEAL_II_NOT_IMPLEMENTED();
        }
}



template <int dim>
void
FE_Q_Hierarchical<dim>::initialize_generalized_support_points()
{
  // number of points: (degree+1)^dim
  unsigned int n = this->degree + 1;
  for (unsigned int i = 1; i < dim; ++i)
    n *= this->degree + 1;

  this->generalized_support_points.resize(n);

  TensorProductPolynomials<dim> *poly_space_derived_ptr =
    dynamic_cast<TensorProductPolynomials<dim> *>(this->poly_space.get());
  const std::vector<unsigned int> &index_map_inverse =
    poly_space_derived_ptr->get_numbering_inverse();

  Point<dim> p;
  // the method of numbering allows
  // each dof to be associated with a
  // support point. There is
  // only one support point per
  // vertex, line, quad, hex, etc.
  //
  // note, however, that the support
  // points thus associated with
  // shape functions are not unique:
  // the linear shape functions are
  // associated with the vertices,
  // but all others are associated
  // with either line, quad, or hex
  // midpoints, and there may be
  // multiple shape functions
  // associated with them. there
  // really is no other useful
  // numbering, since the
  // hierarchical shape functions do
  // not vanish at all-but-one
  // interpolation points (like the
  // Lagrange functions used in
  // FE_Q), so there's not much we
  // can do here.

  // TODO shouldn't we just at least make support points unique,
  // even though the delta property is not satisfied for this FE?
  unsigned int k = 0;
  for (unsigned int iz = 0; iz <= ((dim > 2) ? this->degree : 0); ++iz)
    for (unsigned int iy = 0; iy <= ((dim > 1) ? this->degree : 0); ++iy)
      for (unsigned int ix = 0; ix <= this->degree; ++ix)
        {
          if (ix == 0)
            p[0] = 0.;
          else if (ix == 1)
            p[0] = 1.;
          else
            p[0] = .5;
          if (dim > 1)
            {
              if (iy == 0)
                p[1] = 0.;
              else if (iy == 1)
                p[1] = 1.;
              else
                p[1] = .5;
            }
          if (dim > 2)
            {
              if (iz == 0)
                p[2] = 0.;
              else if (iz == 1)
                p[2] = 1.;
              else
                p[2] = .5;
            }
          this->generalized_support_points[index_map_inverse[k++]] = p;
        }
}



template <>
void
FE_Q_Hierarchical<1>::initialize_generalized_face_support_points()
{
  // no faces in 1d, so nothing to do
}


template <>
void
FE_Q_Hierarchical<1>::get_face_interpolation_matrix(
  const FiniteElement<1, 1> & /*x_source_fe*/,
  FullMatrix<double> & /*interpolation_matrix*/,
  const unsigned int) const
{
  Assert(false, ExcImpossibleInDim(1));
}


template <>
void
FE_Q_Hierarchical<1>::get_subface_interpolation_matrix(
  const FiniteElement<1, 1> & /*x_source_fe*/,
  const unsigned int /*subface*/,
  FullMatrix<double> & /*interpolation_matrix*/,
  const unsigned int) const
{
  Assert(false, ExcImpossibleInDim(1));
}



template <int dim>
void
FE_Q_Hierarchical<dim>::get_face_interpolation_matrix(
  const FiniteElement<dim> &x_source_fe,
  FullMatrix<double>       &interpolation_matrix,
  const unsigned int        face_no) const
{
  // this is only implemented, if the
  // source FE is also a
  // Q_Hierarchical element
  using FEQHierarchical = FE_Q_Hierarchical<dim>;
  AssertThrow((x_source_fe.get_name().find("FE_Q_Hierarchical<") == 0) ||
                (dynamic_cast<const FEQHierarchical *>(&x_source_fe) !=
                 nullptr),
              (typename FiniteElement<dim>::ExcInterpolationNotImplemented()));

  Assert(interpolation_matrix.n() == this->n_dofs_per_face(face_no),
         ExcDimensionMismatch(interpolation_matrix.n(),
                              this->n_dofs_per_face(face_no)));
  Assert(interpolation_matrix.m() == x_source_fe.n_dofs_per_face(face_no),
         ExcDimensionMismatch(interpolation_matrix.m(),
                              x_source_fe.n_dofs_per_face(face_no)));

  // ok, source is a Q_Hierarchical element, so
  // we will be able to do the work
  const FE_Q_Hierarchical<dim> &source_fe =
    dynamic_cast<const FE_Q_Hierarchical<dim> &>(x_source_fe);
  (void)source_fe;

  // Make sure, that the element,
  // for which the DoFs should be
  // constrained is the one with
  // the higher polynomial degree.
  // Actually the procedure will work
  // also if this assertion is not
  // satisfied. But the matrices
  // produced in that case might
  // lead to problems in the
  // hp-procedures, which use this
  // method.
  Assert(this->n_dofs_per_face(face_no) <= source_fe.n_dofs_per_face(face_no),
         (typename FiniteElement<dim>::ExcInterpolationNotImplemented()));
  interpolation_matrix = 0;

  switch (dim)
    {
      case 2:
        {
          // In 2-dimension the constraints are trivial.
          // First this->dofs_per_face DoFs of the constrained
          // element are made equal to the current (dominating)
          // element, which corresponds to 1 on diagonal of the matrix.
          // DoFs which correspond to higher polynomials
          // are zeroed (zero rows in the matrix).
          for (unsigned int i = 0; i < this->n_dofs_per_face(face_no); ++i)
            interpolation_matrix(i, i) = 1;

          break;
        }

      case 3:
        {
          for (unsigned int i = 0; i < GeometryInfo<3>::vertices_per_face; ++i)
            interpolation_matrix(i, i) = 1;

          for (unsigned int i = 0; i < this->degree - 1; ++i)
            {
              for (unsigned int j = 0; j < GeometryInfo<3>::lines_per_face; ++j)
                interpolation_matrix(i + j * (x_source_fe.degree - 1) +
                                       GeometryInfo<3>::vertices_per_face,
                                     i + j * (this->degree - 1) +
                                       GeometryInfo<3>::vertices_per_face) = 1;

              for (unsigned int j = 0; j < this->degree - 1; ++j)
                interpolation_matrix((i + GeometryInfo<3>::lines_per_face) *
                                         (x_source_fe.degree - 1) +
                                       j + GeometryInfo<3>::vertices_per_face,
                                     (i + GeometryInfo<3>::lines_per_face) *
                                         (this->degree - 1) +
                                       j + GeometryInfo<3>::vertices_per_face) =
                  1;
            }
        }
    }
}



template <int dim>
void
FE_Q_Hierarchical<dim>::get_subface_interpolation_matrix(
  const FiniteElement<dim> &x_source_fe,
  const unsigned int        subface,
  FullMatrix<double>       &interpolation_matrix,
  const unsigned int        face_no) const
{
  // this is only implemented, if the
  // source FE is also a
  // Q_Hierarchical element
  using FEQHierarchical = FE_Q_Hierarchical<dim>;
  AssertThrow((x_source_fe.get_name().find("FE_Q_Hierarchical<") == 0) ||
                (dynamic_cast<const FEQHierarchical *>(&x_source_fe) !=
                 nullptr),
              (typename FiniteElement<dim>::ExcInterpolationNotImplemented()));

  Assert(interpolation_matrix.n() == this->n_dofs_per_face(face_no),
         ExcDimensionMismatch(interpolation_matrix.n(),
                              this->n_dofs_per_face(face_no)));
  Assert(interpolation_matrix.m() == x_source_fe.n_dofs_per_face(face_no),
         ExcDimensionMismatch(interpolation_matrix.m(),
                              x_source_fe.n_dofs_per_face(face_no)));

  // ok, source is a Q_Hierarchical element, so
  // we will be able to do the work
  const FE_Q_Hierarchical<dim> &source_fe =
    dynamic_cast<const FE_Q_Hierarchical<dim> &>(x_source_fe);

  // Make sure, that the element,
  // for which the DoFs should be
  // constrained is the one with
  // the higher polynomial degree.
  // Actually the procedure will work
  // also if this assertion is not
  // satisfied. But the matrices
  // produced in that case might
  // lead to problems in the
  // hp-procedures, which use this
  // method.
  Assert(this->n_dofs_per_face(face_no) <= source_fe.n_dofs_per_face(face_no),
         (typename FiniteElement<dim>::ExcInterpolationNotImplemented()));

  switch (dim)
    {
      case 2:
        {
          switch (subface)
            {
              case 0:
                {
                  interpolation_matrix(0, 0) = 1.0;
                  interpolation_matrix(1, 0) = 0.5;
                  interpolation_matrix(1, 1) = 0.5;

                  for (unsigned int dof = 2;
                       dof < this->n_dofs_per_face(face_no);)
                    {
                      interpolation_matrix(1, dof) = -1.0;
                      dof                          = dof + 2;
                    }

                  int factorial_i = 1;

                  for (unsigned int i = 2; i < this->n_dofs_per_face(face_no);
                       ++i)
                    {
                      interpolation_matrix(i, i) = Utilities::pow(0.5, i);
                      factorial_i *= i;
                      int factorial_j  = factorial_i;
                      int factorial_ij = 1;

                      for (unsigned int j = i + 1;
                           j < this->n_dofs_per_face(face_no);
                           ++j)
                        {
                          factorial_ij *= j - i;
                          factorial_j *= j;

                          if (((i + j) & 1) != 0u)
                            interpolation_matrix(i, j) =
                              -1.0 * Utilities::pow(0.5, j) * factorial_j /
                              (factorial_i * factorial_ij);

                          else
                            interpolation_matrix(i, j) =
                              Utilities::pow(0.5, j) * factorial_j /
                              (factorial_i * factorial_ij);
                        }
                    }

                  break;
                }

              case 1:
                {
                  interpolation_matrix(0, 0) = 0.5;
                  interpolation_matrix(0, 1) = 0.5;

                  for (unsigned int dof = 2;
                       dof < this->n_dofs_per_face(face_no);)
                    {
                      interpolation_matrix(0, dof) = -1.0;
                      dof                          = dof + 2;
                    }

                  interpolation_matrix(1, 1) = 1.0;

                  int factorial_i = 1;

                  for (unsigned int i = 2; i < this->n_dofs_per_face(face_no);
                       ++i)
                    {
                      interpolation_matrix(i, i) = Utilities::pow(0.5, i);
                      factorial_i *= i;
                      int factorial_j  = factorial_i;
                      int factorial_ij = 1;

                      for (unsigned int j = i + 1;
                           j < this->n_dofs_per_face(face_no);
                           ++j)
                        {
                          factorial_ij *= j - i;
                          factorial_j *= j;
                          interpolation_matrix(i, j) =
                            Utilities::pow(0.5, j) * factorial_j /
                            (factorial_i * factorial_ij);
                        }
                    }
                }
            }

          break;
        }

      case 3:
        {
          switch (subface)
            {
              case 0:
                {
                  interpolation_matrix(0, 0) = 1.0;
                  interpolation_matrix(1, 0) = 0.5;
                  interpolation_matrix(1, 1) = 0.5;
                  interpolation_matrix(2, 0) = 0.5;
                  interpolation_matrix(2, 2) = 0.5;

                  for (unsigned int i = 0; i < this->degree - 1;)
                    {
                      for (unsigned int line = 0;
                           line < GeometryInfo<3>::lines_per_face;
                           ++line)
                        interpolation_matrix(3,
                                             i + line * (this->degree - 1) +
                                               4) = -0.5;

                      for (unsigned int j = 0; j < this->degree - 1;)
                        {
                          interpolation_matrix(3,
                                               i + (j + 4) * this->degree - j) =
                            1.0;
                          j = j + 2;
                        }

                      interpolation_matrix(1, i + 2 * (this->degree + 1)) =
                        -1.0;
                      interpolation_matrix(2, i + 4) = -1.0;
                      i                              = i + 2;
                    }

                  for (unsigned int vertex = 0;
                       vertex < GeometryInfo<3>::vertices_per_face;
                       ++vertex)
                    interpolation_matrix(3, vertex) = 0.25;

                  int factorial_i = 1;

                  for (unsigned int i = 2; i <= this->degree; ++i)
                    {
                      double tmp = Utilities::pow(0.5, i);
                      interpolation_matrix(i + 2, i + 2)         = tmp;
                      interpolation_matrix(i + 2 * source_fe.degree,
                                           i + 2 * this->degree) = tmp;
                      tmp *= 0.5;
                      interpolation_matrix(i + source_fe.degree + 1, i + 2) =
                        tmp;
                      interpolation_matrix(i + source_fe.degree + 1,
                                           i + this->degree + 1)     = tmp;
                      interpolation_matrix(i + 3 * source_fe.degree - 1,
                                           i + 2 * this->degree)     = tmp;
                      interpolation_matrix(i + 3 * source_fe.degree - 1,
                                           i + 3 * this->degree - 1) = tmp;
                      tmp *= -2.0;

                      for (unsigned int j = 0; j < this->degree - 1;)
                        {
                          interpolation_matrix(i + source_fe.degree + 1,
                                               (i + 2) * this->degree + j + 2 -
                                                 i) = tmp;
                          interpolation_matrix(i + 3 * source_fe.degree - 1,
                                               i + (j + 4) * this->degree - j -
                                                 2) = tmp;
                          j                         = j + 2;
                        }

                      int factorial_k = 1;

                      for (unsigned int j = 2; j <= this->degree; ++j)
                        {
                          interpolation_matrix(i + (j + 2) * source_fe.degree -
                                                 j,
                                               i + (j + 2) * this->degree - j) =
                            Utilities::pow(0.5, i + j);
                          factorial_k *= j;
                          int factorial_kl = 1;
                          int factorial_l  = factorial_k;

                          for (unsigned int k = j + 1; k < this->degree; ++k)
                            {
                              factorial_kl *= k - j;
                              factorial_l *= k;

                              if (((j + k) & 1) != 0u)
                                interpolation_matrix(
                                  i + (j + 2) * source_fe.degree - j,
                                  i + (k + 2) * this->degree - k) =
                                  -1.0 * Utilities::pow(0.5, i + k) *
                                  factorial_l / (factorial_k * factorial_kl);

                              else
                                interpolation_matrix(
                                  i + (j + 2) * source_fe.degree - j,
                                  i + (k + 2) * this->degree - k) =
                                  Utilities::pow(0.5, i + k) * factorial_l /
                                  (factorial_k * factorial_kl);
                            }
                        }

                      factorial_i *= i;
                      int factorial_j  = factorial_i;
                      int factorial_ij = 1;

                      for (unsigned int j = i + 1; j <= this->degree; ++j)
                        {
                          factorial_ij *= j - i;
                          factorial_j *= j;

                          if (((i + j) & 1) != 0u)
                            {
                              tmp = -1.0 * Utilities::pow(0.5, j) *
                                    factorial_j / (factorial_i * factorial_ij);
                              interpolation_matrix(i + 2, j + 2)         = tmp;
                              interpolation_matrix(i + 2 * source_fe.degree,
                                                   j + 2 * this->degree) = tmp;
                              factorial_k                                = 1;

                              for (unsigned int k = 2; k <= this->degree; ++k)
                                {
                                  interpolation_matrix(
                                    i + (k + 2) * source_fe.degree - k,
                                    j + (k + 2) * this->degree - k) =
                                    tmp * Utilities::pow(0.5, k);
                                  factorial_k *= k;
                                  int factorial_l  = factorial_k;
                                  int factorial_kl = 1;

                                  for (unsigned int l = k + 1;
                                       l <= this->degree;
                                       ++l)
                                    {
                                      factorial_kl *= l - k;
                                      factorial_l *= l;

                                      if (((k + l) & 1) != 0u)
                                        interpolation_matrix(
                                          i + (k + 2) * source_fe.degree - k,
                                          j + (l + 2) * this->degree - l) =
                                          -1.0 * tmp * Utilities::pow(0.5, l) *
                                          factorial_l /
                                          (factorial_k * factorial_kl);

                                      else
                                        interpolation_matrix(
                                          i + (k + 2) * source_fe.degree - k,
                                          j + (l + 2) * this->degree - l) =
                                          tmp * Utilities::pow(0.5, l) *
                                          factorial_l /
                                          (factorial_k * factorial_kl);
                                    }
                                }

                              tmp *= 0.5;
                              interpolation_matrix(i + source_fe.degree + 1,
                                                   j + 2)                = tmp;
                              interpolation_matrix(i + source_fe.degree + 1,
                                                   j + this->degree + 1) = tmp;
                              interpolation_matrix(i + 3 * source_fe.degree - 1,
                                                   j + 2 * this->degree) = tmp;
                              interpolation_matrix(i + 3 * source_fe.degree - 1,
                                                   j + 3 * this->degree - 1) =
                                tmp;
                              tmp *= -2.0;

                              for (unsigned int k = 0; k < this->degree - 1;)
                                {
                                  interpolation_matrix(i + source_fe.degree + 1,
                                                       (j + 2) * this->degree +
                                                         k + 2 - j) = tmp;
                                  interpolation_matrix(
                                    i + 3 * source_fe.degree - 1,
                                    j + (k + 4) * this->degree - k - 2) = tmp;
                                  k                                     = k + 2;
                                }
                            }
                          else
                            {
                              tmp = Utilities::pow(0.5, j) * factorial_j /
                                    (factorial_i * factorial_ij);
                              interpolation_matrix(i + 2, j + 2)         = tmp;
                              interpolation_matrix(i + 2 * source_fe.degree,
                                                   j + 2 * this->degree) = tmp;
                              factorial_k                                = 1;

                              for (unsigned int k = 2; k <= this->degree; ++k)
                                {
                                  interpolation_matrix(
                                    i + (k + 2) * source_fe.degree - k,
                                    j + (k + 2) * this->degree - k) =
                                    tmp * Utilities::pow(0.5, k);
                                  factorial_k *= k;
                                  int factorial_l  = factorial_k;
                                  int factorial_kl = 1;

                                  for (unsigned int l = k + 1;
                                       l <= this->degree;
                                       ++l)
                                    {
                                      factorial_kl *= l - k;
                                      factorial_l *= l;

                                      if (((k + l) & 1) != 0u)
                                        interpolation_matrix(
                                          i + (k + 2) * source_fe.degree - k,
                                          j + (l + 2) * this->degree - l) =
                                          -1.0 * tmp * Utilities::pow(0.5, l) *
                                          factorial_l /
                                          (factorial_k * factorial_kl);

                                      else
                                        interpolation_matrix(
                                          i + (k + 2) * source_fe.degree - k,
                                          j + (l + 2) * this->degree - l) =
                                          tmp * Utilities::pow(0.5, l) *
                                          factorial_l /
                                          (factorial_k * factorial_kl);
                                    }
                                }

                              tmp *= 0.5;
                              interpolation_matrix(i + source_fe.degree + 1,
                                                   j + 2)                = tmp;
                              interpolation_matrix(i + source_fe.degree + 1,
                                                   j + this->degree + 1) = tmp;
                              interpolation_matrix(i + 3 * source_fe.degree - 1,
                                                   j + 2 * this->degree) = tmp;
                              interpolation_matrix(i + 3 * source_fe.degree - 1,
                                                   j + 3 * this->degree - 1) =
                                tmp;
                              tmp *= -2.0;

                              for (unsigned int k = 0; k < this->degree - 1;)
                                {
                                  interpolation_matrix(i + source_fe.degree + 1,
                                                       (j + 2) * this->degree +
                                                         k + 2 - j) = tmp;
                                  interpolation_matrix(
                                    i + 3 * source_fe.degree - 1,
                                    j + (k + 4) * this->degree - k - 2) = tmp;
                                  k                                     = k + 2;
                                }
                            }
                        }
                    }

                  break;
                }

              case 1:
                {
                  interpolation_matrix(0, 0) = 0.5;
                  interpolation_matrix(0, 1) = 0.5;
                  interpolation_matrix(1, 1) = 1.0;
                  interpolation_matrix(3, 1) = 0.5;
                  interpolation_matrix(3, 3) = 0.5;

                  for (unsigned int i = 0; i < this->degree - 1;)
                    {
                      for (unsigned int line = 0;
                           line < GeometryInfo<3>::lines_per_face;
                           ++line)
                        interpolation_matrix(2,
                                             i + line * (this->degree - 1) +
                                               4) = -0.5;

                      for (unsigned int j = 0; j < this->degree - 1;)
                        {
                          interpolation_matrix(2,
                                               i + (j + 4) * this->degree - j) =
                            1.0;
                          j = j + 2;
                        }

                      interpolation_matrix(0, i + 2 * (this->degree + 1)) =
                        -1.0;
                      interpolation_matrix(3, i + this->degree + 3) = -1.0;
                      i                                             = i + 2;
                    }

                  for (unsigned int vertex = 0;
                       vertex < GeometryInfo<3>::vertices_per_face;
                       ++vertex)
                    interpolation_matrix(2, vertex) = 0.25;

                  int factorial_i = 1;

                  for (unsigned int i = 2; i <= this->degree; ++i)
                    {
                      double tmp = Utilities::pow(0.5, i + 1);
                      interpolation_matrix(i + 2, i + 2)                = tmp;
                      interpolation_matrix(i + 2, i + this->degree + 1) = tmp;
                      interpolation_matrix(i + 3 * source_fe.degree - 1,
                                           i + 2 * this->degree)        = tmp;
                      interpolation_matrix(i + 3 * source_fe.degree - 1,
                                           i + 3 * this->degree - 1)    = tmp;
                      tmp *= -2.0;

                      for (unsigned int j = 0; j < this->degree - 1;)
                        {
                          interpolation_matrix(i + 2,
                                               j + (i + 2) * this->degree + 2 -
                                                 i) = tmp;
                          interpolation_matrix(i + 3 * source_fe.degree - 1,
                                               i + (j + 4) * this->degree - j -
                                                 2) = tmp;
                          j                         = j + 2;
                        }

                      tmp *= -1.0;
                      interpolation_matrix(i + source_fe.degree + 1,
                                           i + this->degree + 1) = tmp;
                      interpolation_matrix(i + 2 * source_fe.degree,
                                           i + 2 * this->degree) = tmp;
                      factorial_i *= i;
                      int factorial_j  = factorial_i;
                      int factorial_ij = 1;

                      for (unsigned int j = i + 1; j <= this->degree; ++j)
                        {
                          factorial_ij *= j - i;
                          factorial_j *= j;
                          tmp = Utilities::pow(0.5, j) * factorial_j /
                                (factorial_i * factorial_ij);
                          interpolation_matrix(i + 2 * source_fe.degree,
                                               j + 2 * this->degree) = tmp;
                          int factorial_k                            = 1;

                          for (unsigned int k = 2; k <= this->degree; ++k)
                            {
                              interpolation_matrix(
                                i + (k + 2) * source_fe.degree - k,
                                j + (k + 2) * this->degree - k) =
                                tmp * Utilities::pow(0.5, k);
                              factorial_k *= k;
                              int factorial_l  = factorial_k;
                              int factorial_kl = 1;

                              for (unsigned int l = k + 1; l <= this->degree;
                                   ++l)
                                {
                                  factorial_kl *= l - k;
                                  factorial_l *= l;

                                  if (((k + l) & 1) != 0u)
                                    interpolation_matrix(
                                      i + (k + 2) * source_fe.degree - k,
                                      j + (l + 2) * this->degree - l) =
                                      -1.0 * tmp * Utilities::pow(0.5, l) *
                                      factorial_l /
                                      (factorial_k * factorial_kl);

                                  else
                                    interpolation_matrix(
                                      i + (k + 2) * source_fe.degree - k,
                                      j + (l + 2) * this->degree - l) =
                                      tmp * Utilities::pow(0.5, l) *
                                      factorial_l /
                                      (factorial_k * factorial_kl);
                                }
                            }

                          tmp *= -1.0;

                          for (unsigned int k = 0; k < this->degree - 1;)
                            {
                              interpolation_matrix(i + 3 * source_fe.degree - 1,
                                                   j + (k + 4) * this->degree -
                                                     k - 2) = tmp;
                              k                             = k + 2;
                            }

                          tmp *= -0.5;
                          interpolation_matrix(i + 3 * source_fe.degree - 1,
                                               j + 2 * this->degree)     = tmp;
                          interpolation_matrix(i + 3 * source_fe.degree - 1,
                                               j + 3 * this->degree - 1) = tmp;

                          if (((i + j) & 1) != 0u)
                            tmp *= -1.0;

                          interpolation_matrix(i + 2, j + 2) = tmp;
                          interpolation_matrix(i + 2, j + this->degree + 1) =
                            tmp;
                          interpolation_matrix(i + source_fe.degree + 1,
                                               j + this->degree + 1) =
                            2.0 * tmp;
                          tmp *= -2.0;

                          for (unsigned int k = 0; k < this->degree - 1;)
                            {
                              interpolation_matrix(i + 2,
                                                   k + (j + 2) * this->degree +
                                                     2 - j) = tmp;
                              k                             = k + 2;
                            }
                        }

                      int factorial_k = 1;

                      for (unsigned int j = 2; j <= this->degree; ++j)
                        {
                          interpolation_matrix(i + (j + 2) * source_fe.degree -
                                                 j,
                                               i + (j + 2) * this->degree - j) =
                            Utilities::pow(0.5, i + j);
                          factorial_k *= j;
                          int factorial_l  = factorial_k;
                          int factorial_kl = 1;

                          for (unsigned int k = j + 1; k <= this->degree; ++k)
                            {
                              factorial_kl *= k - j;
                              factorial_l *= k;

                              if (((j + k) & 1) != 0u)
                                interpolation_matrix(
                                  i + (j + 2) * source_fe.degree - j,
                                  i + (k + 2) * this->degree - k) =
                                  -1.0 * Utilities::pow(0.5, i + k) *
                                  factorial_l / (factorial_k * factorial_kl);

                              else
                                interpolation_matrix(
                                  i + (j + 2) * source_fe.degree - j,
                                  i + (k + 2) * this->degree - k) =
                                  Utilities::pow(0.5, i + k) * factorial_l /
                                  (factorial_k * factorial_kl);
                            }
                        }
                    }

                  break;
                }

              case 2:
                {
                  interpolation_matrix(0, 0) = 0.5;
                  interpolation_matrix(0, 2) = 0.5;
                  interpolation_matrix(2, 2) = 1.0;
                  interpolation_matrix(3, 2) = 0.5;
                  interpolation_matrix(3, 3) = 0.5;

                  for (unsigned int i = 0; i < this->degree - 1;)
                    {
                      for (unsigned int line = 0;
                           line < GeometryInfo<3>::lines_per_face;
                           ++line)
                        interpolation_matrix(1,
                                             i + line * (this->degree - 1) +
                                               4) = -0.5;

                      for (unsigned int j = 0; j < this->degree - 1;)
                        {
                          interpolation_matrix(1,
                                               i + (j + 4) * this->degree - j) =
                            1.0;
                          j = j + 2;
                        }

                      interpolation_matrix(0, i + 4)                    = -1.0;
                      interpolation_matrix(3, i + 3 * this->degree + 1) = -1.0;
                      i                                                 = i + 2;
                    }

                  for (unsigned int vertex = 0;
                       vertex < GeometryInfo<3>::vertices_per_face;
                       ++vertex)
                    interpolation_matrix(1, vertex) = 0.25;

                  int factorial_i = 1;

                  for (unsigned int i = 2; i <= this->degree; ++i)
                    {
                      double tmp = Utilities::pow(0.5, i);
                      interpolation_matrix(i + 2, i + 2)             = tmp;
                      interpolation_matrix(i + 3 * source_fe.degree - 1,
                                           i + 3 * this->degree - 1) = tmp;
                      tmp *= 0.5;
                      interpolation_matrix(i + source_fe.degree + 1, i + 2) =
                        tmp;
                      interpolation_matrix(i + source_fe.degree + 1,
                                           i + this->degree + 1)     = tmp;
                      interpolation_matrix(i + 2 * source_fe.degree,
                                           i + 2 * this->degree)     = tmp;
                      interpolation_matrix(i + 2 * source_fe.degree,
                                           i + 3 * this->degree - 1) = tmp;
                      tmp *= -2.0;

                      for (unsigned int j = 0; j < this->degree - 1;)
                        {
                          interpolation_matrix(i + source_fe.degree + 1,
                                               j + (i + 2) * this->degree + 2 -
                                                 i) = tmp;
                          interpolation_matrix(i + 2 * source_fe.degree,
                                               i + (j + 4) * this->degree - j -
                                                 2) = tmp;
                          j                         = j + 2;
                        }

                      int factorial_k = 1;

                      for (unsigned int j = 2; j <= this->degree; ++j)
                        {
                          interpolation_matrix(i + (j + 2) * source_fe.degree -
                                                 j,
                                               i + (j + 2) * this->degree - j) =
                            Utilities::pow(0.5, i + j);
                          factorial_k *= j;
                          int factorial_kl = 1;
                          int factorial_l  = factorial_k;

                          for (unsigned int k = j + 1; k <= this->degree; ++k)
                            {
                              factorial_kl *= k - j;
                              factorial_l *= k;
                              interpolation_matrix(
                                i + (j + 2) * source_fe.degree - j,
                                i + (k + 2) * this->degree - k) =
                                Utilities::pow(0.5, i + k) * factorial_l /
                                (factorial_k * factorial_kl);
                            }
                        }

                      factorial_i *= i;
                      int factorial_j  = factorial_i;
                      int factorial_ij = 1;

                      for (unsigned int j = i + 1; j <= this->degree; ++j)
                        {
                          factorial_ij *= j - i;
                          factorial_j *= j;
                          tmp = Utilities::pow(0.5, j) * factorial_j /
                                (factorial_i * factorial_ij);
                          interpolation_matrix(i + 2, j + 2) = tmp;
                          tmp *= -1.0;

                          for (unsigned int k = 0; k < this->degree - 1;)
                            {
                              interpolation_matrix(i + source_fe.degree + 1,
                                                   k + (j + 2) * this->degree +
                                                     2 - j) = tmp;
                              k                             = k + 2;
                            }

                          tmp *= -0.5;
                          interpolation_matrix(i + source_fe.degree + 1,
                                               j + 2)                = tmp;
                          interpolation_matrix(i + source_fe.degree + 1,
                                               j + this->degree + 1) = tmp;

                          if (((i + j) & 1) != 0u)
                            tmp *= -1.0;

                          interpolation_matrix(i + 2 * source_fe.degree,
                                               j + 2 * this->degree)     = tmp;
                          interpolation_matrix(i + 2 * source_fe.degree,
                                               j + 3 * this->degree - 1) = tmp;
                          tmp *= 2.0;
                          interpolation_matrix(i + 3 * source_fe.degree - 1,
                                               j + 3 * this->degree - 1) = tmp;
                          int factorial_k                                = 1;

                          for (unsigned int k = 2; k <= this->degree; ++k)
                            {
                              interpolation_matrix(
                                i + (k + 2) * source_fe.degree - k,
                                j + (k + 2) * this->degree - k) =
                                tmp * Utilities::pow(0.5, k);
                              factorial_k *= k;
                              int factorial_l  = factorial_k;
                              int factorial_kl = 1;

                              for (unsigned int l = k + 1; l <= this->degree;
                                   ++l)
                                {
                                  factorial_kl *= l - k;
                                  factorial_l *= l;
                                  interpolation_matrix(
                                    i + (k + 2) * source_fe.degree - k,
                                    j + (l + 2) * this->degree - l) =
                                    tmp * Utilities::pow(0.5, l) * factorial_l /
                                    (factorial_k * factorial_kl);
                                }
                            }

                          tmp *= -1.0;

                          for (unsigned int k = 0; k < this->degree - 1;)
                            {
                              interpolation_matrix(i + 2 * source_fe.degree,
                                                   j + (k + 4) * this->degree -
                                                     k - 2) = tmp;
                              k                             = k + 2;
                            }
                        }
                    }

                  break;
                }

              case 3:
                {
                  for (unsigned int vertex = 0;
                       vertex < GeometryInfo<3>::vertices_per_face;
                       ++vertex)
                    interpolation_matrix(0, vertex) = 0.25;

                  for (unsigned int i = 0; i < this->degree - 1;)
                    {
                      for (unsigned int line = 0;
                           line < GeometryInfo<3>::lines_per_face;
                           ++line)
                        interpolation_matrix(0,
                                             i + line * (this->degree - 1) +
                                               4) = -0.5;

                      for (unsigned int j = 0; j < this->degree - 1;)
                        {
                          interpolation_matrix(0,
                                               i + (j + 4) * this->degree - j) =
                            1.0;
                          j = j + 2;
                        }

                      interpolation_matrix(1, i + 4)                    = -1.0;
                      interpolation_matrix(2, i + 3 * this->degree + 1) = -1.0;
                      i                                                 = i + 2;
                    }

                  interpolation_matrix(1, 0) = 0.5;
                  interpolation_matrix(1, 1) = 0.5;
                  interpolation_matrix(2, 2) = 0.5;
                  interpolation_matrix(2, 3) = 0.5;
                  interpolation_matrix(3, 3) = 1.0;

                  int factorial_i = 1;

                  for (unsigned int i = 2; i <= this->degree; ++i)
                    {
                      double tmp = Utilities::pow(0.5, i + 1);
                      interpolation_matrix(i + 2, i + 2)                = tmp;
                      interpolation_matrix(i + 2, i + this->degree + 1) = tmp;
                      interpolation_matrix(i + 2 * source_fe.degree,
                                           i + 2 * this->degree)        = tmp;
                      interpolation_matrix(i + 2 * source_fe.degree,
                                           i + 3 * this->degree - 1)    = tmp;
                      tmp *= -2.0;

                      for (unsigned int j = 0; j < this->degree - 1;)
                        {
                          interpolation_matrix(i + 2,
                                               j + (i + 2) * this->degree + 2 -
                                                 i) = tmp;
                          interpolation_matrix(i + 2 * source_fe.degree,
                                               i + (j + 4) * this->degree - 2) =
                            tmp;
                          j = j + 2;
                        }

                      tmp *= -1.0;
                      interpolation_matrix(i + source_fe.degree + 1,
                                           i + this->degree + 1)     = tmp;
                      interpolation_matrix(i + 3 * source_fe.degree - 1,
                                           i + 3 * this->degree - 1) = tmp;
                      int factorial_k                                = 1;

                      for (unsigned int j = 2; j <= this->degree; ++j)
                        {
                          interpolation_matrix(i + (j + 2) * source_fe.degree -
                                                 j,
                                               i + (j + 2) * this->degree - j) =
                            Utilities::pow(0.5, i + j);
                          factorial_k *= j;
                          int factorial_l  = factorial_k;
                          int factorial_kl = 1;

                          for (unsigned int k = j + 1; k <= this->degree; ++k)
                            {
                              factorial_kl *= k - j;
                              factorial_l *= k;
                              interpolation_matrix(
                                i + (j + 2) * source_fe.degree - j,
                                i + (k + 2) * this->degree - k) =
                                Utilities::pow(0.5, i + k) * factorial_l /
                                (factorial_k * factorial_kl);
                            }
                        }

                      factorial_i *= i;
                      int factorial_j  = factorial_i;
                      int factorial_ij = 1;

                      for (unsigned int j = i + 1; j <= this->degree; ++j)
                        {
                          factorial_ij *= j - i;
                          factorial_j *= j;
                          tmp = Utilities::pow(0.5, j + 1) * factorial_j /
                                (factorial_i * factorial_ij);
                          interpolation_matrix(i + 2, j + 2) = tmp;
                          interpolation_matrix(i + 2, j + this->degree + 1) =
                            tmp;
                          interpolation_matrix(i + 2 * source_fe.degree,
                                               j + 2 * this->degree)     = tmp;
                          interpolation_matrix(i + 2 * source_fe.degree,
                                               j + 3 * this->degree - 1) = tmp;
                          tmp *= 2.0;
                          interpolation_matrix(i + source_fe.degree + 1,
                                               j + this->degree + 1)     = tmp;
                          interpolation_matrix(i + 3 * source_fe.degree - 1,
                                               j + 3 * this->degree - 1) = tmp;
                          int factorial_k                                = 1;

                          for (unsigned int k = 2; k <= this->degree; ++k)
                            {
                              interpolation_matrix(
                                i + (k + 2) * source_fe.degree - k,
                                j + (k + 2) * this->degree - k) =
                                tmp * Utilities::pow(0.5, k);
                              factorial_k *= k;
                              int factorial_l  = factorial_k;
                              int factorial_kl = 1;

                              for (unsigned int l = k + 1; l <= this->degree;
                                   ++l)
                                {
                                  factorial_kl *= l - k;
                                  factorial_l *= l;
                                  interpolation_matrix(
                                    i + (k + 2) * source_fe.degree - k,
                                    j + (l + 2) * this->degree - l) =
                                    tmp * Utilities::pow(0.5, l) * factorial_l /
                                    (factorial_k * factorial_kl);
                                }
                            }

                          tmp *= -1.0;

                          for (unsigned int k = 0; k < this->degree - 1;)
                            {
                              interpolation_matrix(i + 2,
                                                   k + (j + 2) * this->degree +
                                                     2 - j) = tmp;
                              interpolation_matrix(i + 2 * source_fe.degree,
                                                   j + (k + 4) * this->degree -
                                                     2)     = tmp;
                              k                             = k + 2;
                            }
                        }
                    }
                }
            }
        }
    }
}



template <int dim>
void
FE_Q_Hierarchical<dim>::initialize_generalized_face_support_points()
{
  const unsigned int codim = dim - 1;

  // TODO: the implementation makes the assumption that all faces have the
  // same number of dofs
  AssertDimension(this->n_unique_faces(), 1);
  const unsigned int face_no = 0;

  // number of points: (degree+1)^codim
  unsigned int n = this->degree + 1;
  for (unsigned int i = 1; i < codim; ++i)
    n *= this->degree + 1;

  this->generalized_face_support_points[face_no].resize(n);

  Point<codim> p;

  unsigned int k = 0;
  for (unsigned int iz = 0; iz <= ((codim > 2) ? this->degree : 0); ++iz)
    for (unsigned int iy = 0; iy <= ((codim > 1) ? this->degree : 0); ++iy)
      for (unsigned int ix = 0; ix <= this->degree; ++ix)
        {
          if (ix == 0)
            p[0] = 0.;
          else if (ix == 1)
            p[0] = 1.;
          else
            p[0] = .5;
          if (codim > 1)
            {
              if (iy == 0)
                p[1] = 0.;
              else if (iy == 1)
                p[1] = 1.;
              else
                p[1] = .5;
            }
          if (codim > 2)
            {
              if (iz == 0)
                p[2] = 0.;
              else if (iz == 1)
                p[2] = 1.;
              else
                p[2] = .5;
            }
          this->generalized_face_support_points[face_no][face_renumber[k++]] =
            p;
        }
}


// we use same dpo_vector as FE_Q
template <int dim>
std::vector<unsigned int>
FE_Q_Hierarchical<dim>::get_dpo_vector(const unsigned int deg)
{
  std::vector<unsigned int> dpo(dim + 1, 1U);
  for (unsigned int i = 1; i < dpo.size(); ++i)
    dpo[i] = dpo[i - 1] * (deg - 1);
  return dpo;
}



template <int dim>
std::vector<unsigned int>
FE_Q_Hierarchical<dim>::hierarchic_to_fe_q_hierarchical_numbering(
  const FiniteElementData<dim> &fe)
{
  Assert(fe.n_components() == 1, ExcInternalError());
  std::vector<unsigned int> h2l(fe.n_dofs_per_cell());

  // polynomial degree
  const unsigned int degree = fe.n_dofs_per_line() + 1;
  // number of grid points in each
  // direction
  const unsigned int n = degree + 1;

  // the following lines of code are
  // somewhat odd, due to the way the
  // hierarchic numbering is
  // organized. if someone would
  // really want to understand these
  // lines, you better draw some
  // pictures where you indicate the
  // indices and orders of vertices,
  // lines, etc, along with the
  // numbers of the degrees of
  // freedom in hierarchical and
  // lexicographical order
  switch (dim)
    {
      case 1:
        {
          for (unsigned int i = 0; i < fe.n_dofs_per_cell(); ++i)
            h2l[i] = i;

          break;
        }

      case 2:
        {
          // Example: degree=3
          //
          // hierarchical numbering:
          //  2 10 11  3
          //  5 14 15  7
          //  4 12 13  6
          //  0  8  9  1
          //
          // fe_q_hierarchical numbering:
          //  4  6  7  5
          // 12 14 15 13
          //  8 10 11  9
          //  0  2  3  1
          unsigned int next_index = 0;
          // first the four vertices
          h2l[next_index++] = 0;
          h2l[next_index++] = 1;
          h2l[next_index++] = n;
          h2l[next_index++] = n + 1;
          // left line
          for (unsigned int i = 0; i < fe.n_dofs_per_line(); ++i)
            h2l[next_index++] = (2 + i) * n;
          // right line
          for (unsigned int i = 0; i < fe.n_dofs_per_line(); ++i)
            h2l[next_index++] = (2 + i) * n + 1;
          // bottom line
          for (unsigned int i = 0; i < fe.n_dofs_per_line(); ++i)
            h2l[next_index++] = 2 + i;
          // top line
          for (unsigned int i = 0; i < fe.n_dofs_per_line(); ++i)
            h2l[next_index++] = n + 2 + i;
          // inside quad
          Assert(fe.n_dofs_per_quad(0 /*only one quad in 2d*/) ==
                   fe.n_dofs_per_line() * fe.n_dofs_per_line(),
                 ExcInternalError());
          for (unsigned int i = 0; i < fe.n_dofs_per_line(); ++i)
            for (unsigned int j = 0; j < fe.n_dofs_per_line(); ++j)
              h2l[next_index++] = (2 + i) * n + 2 + j;

          Assert(next_index == fe.n_dofs_per_cell(), ExcInternalError());

          break;
        }

      case 3:
        {
          unsigned int       next_index = 0;
          const unsigned int n2         = n * n;
          // first the eight vertices
          // bottom face, lexicographic
          h2l[next_index++] = 0;
          h2l[next_index++] = 1;
          h2l[next_index++] = n;
          h2l[next_index++] = n + 1;
          // top face, lexicographic
          h2l[next_index++] = n2;
          h2l[next_index++] = n2 + 1;
          h2l[next_index++] = n2 + n;
          h2l[next_index++] = n2 + n + 1;

          // now the lines
          // bottom face
          for (unsigned int i = 0; i < fe.n_dofs_per_line(); ++i)
            h2l[next_index++] = (2 + i) * n;
          for (unsigned int i = 0; i < fe.n_dofs_per_line(); ++i)
            h2l[next_index++] = (2 + i) * n + 1;
          for (unsigned int i = 0; i < fe.n_dofs_per_line(); ++i)
            h2l[next_index++] = 2 + i;
          for (unsigned int i = 0; i < fe.n_dofs_per_line(); ++i)
            h2l[next_index++] = n + 2 + i;
          // top face
          for (unsigned int i = 0; i < fe.n_dofs_per_line(); ++i)
            h2l[next_index++] = n2 + (2 + i) * n;
          for (unsigned int i = 0; i < fe.n_dofs_per_line(); ++i)
            h2l[next_index++] = n2 + (2 + i) * n + 1;
          for (unsigned int i = 0; i < fe.n_dofs_per_line(); ++i)
            h2l[next_index++] = n2 + 2 + i;
          for (unsigned int i = 0; i < fe.n_dofs_per_line(); ++i)
            h2l[next_index++] = n2 + n + 2 + i;
          // lines in z-direction
          for (unsigned int i = 0; i < fe.n_dofs_per_line(); ++i)
            h2l[next_index++] = (2 + i) * n2;
          for (unsigned int i = 0; i < fe.n_dofs_per_line(); ++i)
            h2l[next_index++] = (2 + i) * n2 + 1;
          for (unsigned int i = 0; i < fe.n_dofs_per_line(); ++i)
            h2l[next_index++] = (2 + i) * n2 + n;
          for (unsigned int i = 0; i < fe.n_dofs_per_line(); ++i)
            h2l[next_index++] = (2 + i) * n2 + n + 1;

          // TODO: the implementation makes the assumption that all faces have
          // the same number of dofs
          AssertDimension(fe.n_unique_faces(), 1);
          const unsigned int face_no = 0;
          (void)face_no;

          // inside quads
          Assert(fe.n_dofs_per_quad(face_no) ==
                   fe.n_dofs_per_line() * fe.n_dofs_per_line(),
                 ExcInternalError());
          // left face
          for (unsigned int i = 0; i < fe.n_dofs_per_line(); ++i)
            for (unsigned int j = 0; j < fe.n_dofs_per_line(); ++j)
              h2l[next_index++] = (2 + i) * n2 + (2 + j) * n;
          // right face
          for (unsigned int i = 0; i < fe.n_dofs_per_line(); ++i)
            for (unsigned int j = 0; j < fe.n_dofs_per_line(); ++j)
              h2l[next_index++] = (2 + i) * n2 + (2 + j) * n + 1;
          // front face
          for (unsigned int i = 0; i < fe.n_dofs_per_line(); ++i)
            for (unsigned int j = 0; j < fe.n_dofs_per_line(); ++j)
              h2l[next_index++] = (2 + i) * n2 + 2 + j;
          // back face
          for (unsigned int i = 0; i < fe.n_dofs_per_line(); ++i)
            for (unsigned int j = 0; j < fe.n_dofs_per_line(); ++j)
              h2l[next_index++] = (2 + i) * n2 + n + 2 + j;
          // bottom face
          for (unsigned int i = 0; i < fe.n_dofs_per_line(); ++i)
            for (unsigned int j = 0; j < fe.n_dofs_per_line(); ++j)
              h2l[next_index++] = (2 + i) * n + 2 + j;
          // top face
          for (unsigned int i = 0; i < fe.n_dofs_per_line(); ++i)
            for (unsigned int j = 0; j < fe.n_dofs_per_line(); ++j)
              h2l[next_index++] = n2 + (2 + i) * n + 2 + j;

          // inside hex
          Assert(fe.n_dofs_per_hex() ==
                   fe.n_dofs_per_quad(face_no) * fe.n_dofs_per_line(),
                 ExcInternalError());
          for (unsigned int i = 0; i < fe.n_dofs_per_line(); ++i)
            for (unsigned int j = 0; j < fe.n_dofs_per_line(); ++j)
              for (unsigned int k = 0; k < fe.n_dofs_per_line(); ++k)
                h2l[next_index++] = (2 + i) * n2 + (2 + j) * n + 2 + k;

          Assert(next_index == fe.n_dofs_per_cell(), ExcInternalError());

          break;
        }

      default:
        DEAL_II_NOT_IMPLEMENTED();
    }
  return h2l;
}


template <int dim>
std::vector<unsigned int>
FE_Q_Hierarchical<dim>::face_fe_q_hierarchical_to_hierarchic_numbering(
  const unsigned int degree)
{
  FiniteElementData<dim - 1> fe_data(
    FE_Q_Hierarchical<dim - 1>::get_dpo_vector(degree), 1, degree);
  return internal::FE_Q_Hierarchical::invert_numbering(
    FE_Q_Hierarchical<dim - 1>::hierarchic_to_fe_q_hierarchical_numbering(
      fe_data));
}



template <>
std::vector<unsigned int>
FE_Q_Hierarchical<1>::face_fe_q_hierarchical_to_hierarchic_numbering(
  const unsigned int)
{
  return std::vector<unsigned int>();
}


template <>
bool
FE_Q_Hierarchical<1>::has_support_on_face(const unsigned int shape_index,
                                          const unsigned int face_index) const
{
  AssertIndexRange(shape_index, this->n_dofs_per_cell());
  AssertIndexRange(face_index, GeometryInfo<1>::faces_per_cell);


  // in 1d, things are simple. since
  // there is only one degree of
  // freedom per vertex in this
  // class, the first is on vertex 0
  // (==face 0 in some sense), the
  // second on face 1:
  return (((shape_index == 0) && (face_index == 0)) ||
          ((shape_index == 1) && (face_index == 1)));
}



template <int dim>
bool
FE_Q_Hierarchical<dim>::has_support_on_face(const unsigned int shape_index,
                                            const unsigned int face_index) const
{
  AssertIndexRange(shape_index, this->n_dofs_per_cell());
  AssertIndexRange(face_index, GeometryInfo<dim>::faces_per_cell);

  // first, special-case interior
  // shape functions, since they
  // have no support no-where on
  // the boundary
  if (((dim == 2) && (shape_index >=
                      this->get_first_quad_index(0 /*only one quad in 2d*/))) ||
      ((dim == 3) && (shape_index >= this->get_first_hex_index())))
    return false;

  // let's see whether this is a
  // vertex
  if (shape_index < this->get_first_line_index())
    {
      // for Q elements, there is
      // one dof per vertex, so
      // shape_index==vertex_number. check
      // whether this vertex is
      // on the given face.
      const unsigned int vertex_no = shape_index;
      Assert(vertex_no < GeometryInfo<dim>::vertices_per_cell,
             ExcInternalError());
      for (unsigned int i = 0; i < GeometryInfo<dim>::vertices_per_face; ++i)
        if (GeometryInfo<dim>::face_to_cell_vertices(face_index, i) ==
            vertex_no)
          return true;
      return false;
    }
  else if (shape_index < this->get_first_quad_index(0))
    // ok, dof is on a line
    {
      const unsigned int line_index =
        (shape_index - this->get_first_line_index()) / this->n_dofs_per_line();
      Assert(line_index < GeometryInfo<dim>::lines_per_cell,
             ExcInternalError());

      for (unsigned int i = 0; i < GeometryInfo<dim>::lines_per_face; ++i)
        if (GeometryInfo<dim>::face_to_cell_lines(face_index, i) == line_index)
          return true;
      return false;
    }
  else if (shape_index < this->get_first_hex_index())
    // dof is on a quad
    {
      const unsigned int quad_index =
        (shape_index - this->get_first_quad_index(0 /*first quad*/)) /
        this->n_dofs_per_quad(face_index);
      Assert(static_cast<signed int>(quad_index) <
               static_cast<signed int>(GeometryInfo<dim>::quads_per_cell),
             ExcInternalError());

      // in 2d, cell bubble are
      // zero on all faces. but
      // we have treated this
      // case above already
      Assert(dim != 2, ExcInternalError());

      // in 3d,
      // quad_index=face_index
      if (dim == 3)
        return (quad_index == face_index);
      else
        DEAL_II_NOT_IMPLEMENTED();
    }
  else
    // dof on hex
    {
      // can only happen in 3d, but
      // this case has already been
      // covered above
      DEAL_II_NOT_IMPLEMENTED();
      return false;
    }

  // we should not have gotten here
  DEAL_II_ASSERT_UNREACHABLE();
  return false;
}



template <int dim>
std::vector<unsigned int>
FE_Q_Hierarchical<dim>::get_embedding_dofs(const unsigned int sub_degree) const
{
  Assert((sub_degree > 0) && (sub_degree <= this->degree),
         ExcIndexRange(sub_degree, 1, this->degree));

  if (dim == 1)
    {
      std::vector<unsigned int> embedding_dofs(sub_degree + 1);
      for (unsigned int i = 0; i < (sub_degree + 1); ++i)
        embedding_dofs[i] = i;

      return embedding_dofs;
    }

  if (sub_degree == 1)
    {
      std::vector<unsigned int> embedding_dofs(
        GeometryInfo<dim>::vertices_per_cell);
      for (const unsigned int i : GeometryInfo<dim>::vertex_indices())
        embedding_dofs[i] = i;

      return embedding_dofs;
    }
  else if (sub_degree == this->degree)
    {
      std::vector<unsigned int> embedding_dofs(this->n_dofs_per_cell());
      for (unsigned int i = 0; i < this->n_dofs_per_cell(); ++i)
        embedding_dofs[i] = i;

      return embedding_dofs;
    }

  if ((dim == 2) || (dim == 3))
    {
      std::vector<unsigned int> embedding_dofs(
        (dim == 2) ? (sub_degree + 1) * (sub_degree + 1) :
                     (sub_degree + 1) * (sub_degree + 1) * (sub_degree + 1));

      for (unsigned int i = 0;
           i < ((dim == 2) ?
                  (sub_degree + 1) * (sub_degree + 1) :
                  (sub_degree + 1) * (sub_degree + 1) * (sub_degree + 1));
           ++i)
        {
          // vertex
          if (i < GeometryInfo<dim>::vertices_per_cell)
            embedding_dofs[i] = i;
          // line
          else if (i < (GeometryInfo<dim>::vertices_per_cell +
                        GeometryInfo<dim>::lines_per_cell * (sub_degree - 1)))
            {
              const unsigned int j =
                (i - GeometryInfo<dim>::vertices_per_cell) % (sub_degree - 1);
              const unsigned int line =
                (i - GeometryInfo<dim>::vertices_per_cell - j) /
                (sub_degree - 1);

              embedding_dofs[i] = GeometryInfo<dim>::vertices_per_cell +
                                  line * (this->degree - 1) + j;
            }
          // quad
          else if (i < (GeometryInfo<dim>::vertices_per_cell +
                        GeometryInfo<dim>::lines_per_cell * (sub_degree - 1)) +
                         GeometryInfo<dim>::quads_per_cell * (sub_degree - 1) *
                           (sub_degree - 1))
            {
              const unsigned int j =
                (i - GeometryInfo<dim>::vertices_per_cell -
                 GeometryInfo<dim>::lines_per_cell * (sub_degree - 1)) %
                (sub_degree - 1);
              const unsigned int k =
                ((i - GeometryInfo<dim>::vertices_per_cell -
                  GeometryInfo<dim>::lines_per_cell * (sub_degree - 1) - j) /
                 (sub_degree - 1)) %
                (sub_degree - 1);
              const unsigned int face =
                (i - GeometryInfo<dim>::vertices_per_cell -
                 GeometryInfo<dim>::lines_per_cell * (sub_degree - 1) -
                 k * (sub_degree - 1) - j) /
                ((sub_degree - 1) * (sub_degree - 1));

              embedding_dofs[i] =
                GeometryInfo<dim>::vertices_per_cell +
                GeometryInfo<dim>::lines_per_cell * (this->degree - 1) +
                face * (this->degree - 1) * (this->degree - 1) +
                k * (this->degree - 1) + j;
            }
          // hex
          else if (i < (GeometryInfo<dim>::vertices_per_cell +
                        GeometryInfo<dim>::lines_per_cell * (sub_degree - 1)) +
                         GeometryInfo<dim>::quads_per_cell * (sub_degree - 1) *
                           (sub_degree - 1) +
                         GeometryInfo<dim>::hexes_per_cell * (sub_degree - 1) *
                           (sub_degree - 1) * (sub_degree - 1))
            {
              const unsigned int j =
                (i - GeometryInfo<dim>::vertices_per_cell -
                 GeometryInfo<dim>::lines_per_cell * (sub_degree - 1) -
                 GeometryInfo<dim>::quads_per_cell * (sub_degree - 1) *
                   (sub_degree - 1)) %
                (sub_degree - 1);
              const unsigned int k =
                ((i - GeometryInfo<dim>::vertices_per_cell -
                  GeometryInfo<dim>::lines_per_cell * (sub_degree - 1) -
                  GeometryInfo<dim>::quads_per_cell * (sub_degree - 1) *
                    (sub_degree - 1) -
                  j) /
                 (sub_degree - 1)) %
                (sub_degree - 1);
              const unsigned int l =
                (i - GeometryInfo<dim>::vertices_per_cell -
                 GeometryInfo<dim>::lines_per_cell * (sub_degree - 1) -
                 GeometryInfo<dim>::quads_per_cell * (sub_degree - 1) *
                   (sub_degree - 1) -
                 j - k * (sub_degree - 1)) /
                ((sub_degree - 1) * (sub_degree - 1));

              embedding_dofs[i] =
                GeometryInfo<dim>::vertices_per_cell +
                GeometryInfo<dim>::lines_per_cell * (this->degree - 1) +
                GeometryInfo<dim>::quads_per_cell * (this->degree - 1) *
                  (this->degree - 1) +
                l * (this->degree - 1) * (this->degree - 1) +
                k * (this->degree - 1) + j;
            }
        }

      return embedding_dofs;
    }
  else
    {
      DEAL_II_NOT_IMPLEMENTED();
      return std::vector<unsigned int>();
    }
}



template <int dim>
std::pair<Table<2, bool>, std::vector<unsigned int>>
FE_Q_Hierarchical<dim>::get_constant_modes() const
{
  Table<2, bool> constant_modes(1, this->n_dofs_per_cell());
  for (const unsigned int i : GeometryInfo<dim>::vertex_indices())
    constant_modes(0, i) = true;
  for (unsigned int i = GeometryInfo<dim>::vertices_per_cell;
       i < this->n_dofs_per_cell();
       ++i)
    constant_modes(0, i) = false;
  return std::pair<Table<2, bool>, std::vector<unsigned int>>(
    constant_modes, std::vector<unsigned int>(1, 0));
}



template <int dim>
std::size_t
FE_Q_Hierarchical<dim>::memory_consumption() const
{
  DEAL_II_NOT_IMPLEMENTED();
  return 0;
}



// explicit instantiations
#include "fe/fe_q_hierarchical.inst"


DEAL_II_NAMESPACE_CLOSE
