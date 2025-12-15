// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2023 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#include <deal.II/base/config.h>

#include <deal.II/base/quadrature.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/table.h>
#include <deal.II/base/template_constraints.h>

#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_face.h>
#include <deal.II/fe/fe_hermite.h>
#include <deal.II/fe/fe_nothing.h>
#include <deal.II/fe/fe_pyramid_p.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_simplex_p.h>
#include <deal.II/fe/fe_simplex_p_bubbles.h>
#include <deal.II/fe/fe_tools.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_wedge_p.h>
#include <deal.II/fe/mapping_cartesian.h>

#include <deal.II/grid/reference_cell.h>

#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_control.h>
#include <deal.II/lac/sparse_direct.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/sparsity_pattern.h>
#include <deal.II/lac/sparsity_tools.h>
#include <deal.II/lac/vector.h>

#include <deal.II/numerics/matrix_tools.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <iterator>
#include <memory>
#include <sstream>

DEAL_II_NAMESPACE_OPEN



namespace internal
{
  inline unsigned int
  get_regularity_from_degree(const unsigned int fe_degree)
  {
    Assert(fe_degree % 2 == 1,
           ExcMessage("FE_Hermite only supports odd polynomial degrees."));
    return (fe_degree == 0) ? 0 : (fe_degree - 1) / 2;
  }



  inline std::vector<unsigned int>
  get_hermite_dpo_vector(const unsigned int dim, const unsigned int regularity)
  {
    std::vector<unsigned int> result(dim + 1, 0);
    result[0] = Utilities::pow(regularity + 1, dim);

    return result;
  }



  /*
   * Renumbering function. Function needs different levels of for loop nesting
   * for different values of dim, so different definitions are used for
   * simplicity.
   */
  template <int dim>
  void
  hermite_hierarchic_to_lexicographic_numbering(const unsigned int regularity,
                                                std::vector<unsigned int> &h2l);



  template <>
  void
  hermite_hierarchic_to_lexicographic_numbering<1>(
    const unsigned int         regularity,
    std::vector<unsigned int> &h2l)
  {
    const unsigned int node_dofs_1d = regularity + 1;

    AssertDimension(h2l.size(), 2 * node_dofs_1d);

    // Assign DOFs at vertices
    for (unsigned int di = 0; di < 2; ++di)
      for (unsigned int i = 0; i < node_dofs_1d; ++i)
        h2l[i + di * node_dofs_1d] = i + di * node_dofs_1d;
  }



  template <>
  void
  hermite_hierarchic_to_lexicographic_numbering<2>(
    const unsigned int         regularity,
    std::vector<unsigned int> &h2l)
  {
    const unsigned int node_dofs_1d = regularity + 1;
    const unsigned int dim_dofs_1d  = 2 * node_dofs_1d;
    unsigned int       offset       = 0;

    AssertDimension(h2l.size(), dim_dofs_1d * dim_dofs_1d);

    // Assign DOFs at vertices
    for (unsigned int di = 0; di < 2; ++di)
      for (unsigned int dj = 0; dj < 2; ++dj)
        {
          for (unsigned int i = 0; i < node_dofs_1d; ++i)
            for (unsigned int j = 0; j < node_dofs_1d; ++j)
              h2l[j + i * node_dofs_1d + offset] =
                j + i * dim_dofs_1d + (dj + di * dim_dofs_1d) * node_dofs_1d;

          offset += node_dofs_1d * node_dofs_1d;
        }
  }



  template <>
  void
  hermite_hierarchic_to_lexicographic_numbering<3>(
    const unsigned int         regularity,
    std::vector<unsigned int> &h2l)
  {
    const unsigned int node_dofs_1d = regularity + 1;
    const unsigned int node_dofs_2d = node_dofs_1d * node_dofs_1d;

    const unsigned int dim_dofs_1d = 2 * node_dofs_1d;
    const unsigned int dim_dofs_2d = dim_dofs_1d * dim_dofs_1d;

    unsigned int offset = 0;

    AssertDimension(h2l.size(), dim_dofs_2d * dim_dofs_1d);

    // Assign DOFs at nodes
    for (unsigned int di = 0; di < 2; ++di)
      for (unsigned int dj = 0; dj < 2; ++dj)
        for (unsigned int dk = 0; dk < 2; ++dk)
          {
            for (unsigned int i = 0; i < node_dofs_1d; ++i)
              for (unsigned int j = 0; j < node_dofs_1d; ++j)
                for (unsigned int k = 0; k < node_dofs_1d; ++k)
                  h2l[k + j * node_dofs_1d + i * node_dofs_2d + offset] =
                    k + j * dim_dofs_1d + i * dim_dofs_2d +
                    node_dofs_1d * (dk + dj * dim_dofs_1d + di * dim_dofs_2d);

            offset += node_dofs_1d * node_dofs_2d;
          }
  }



  template <int dim>
  inline std::vector<unsigned int>
  hermite_hierarchic_to_lexicographic_numbering(const unsigned int regularity)
  {
    const std::vector<unsigned int> dpo =
      get_hermite_dpo_vector(dim, regularity);
    const dealii::FiniteElementData<dim> face_data(dpo, 1, 2 * regularity + 1);
    std::vector<unsigned int>            renumbering(face_data.dofs_per_cell);

    hermite_hierarchic_to_lexicographic_numbering<dim>(regularity, renumbering);

    return renumbering;
  }



  template <int dim>
  std::vector<unsigned int>
  hermite_lexicographic_to_hierarchic_numbering(const unsigned int regularity)
  {
    return Utilities::invert_permutation(
      hermite_hierarchic_to_lexicographic_numbering<dim>(regularity));
  }



  template <int dim>
  inline std::vector<unsigned int>
  hermite_face_lexicographic_to_hierarchic_numbering(
    const unsigned int regularity)
  {
    (void)regularity;
    if constexpr (dim > 1)
      return hermite_lexicographic_to_hierarchic_numbering<dim - 1>(regularity);
    else
      return std::vector<unsigned int>();
  }



  template <int dim>
  TensorProductPolynomials<dim>
  get_hermite_polynomials(const unsigned int fe_degree)
  {
    const unsigned int regularity = get_regularity_from_degree(fe_degree);

    TensorProductPolynomials<dim> polynomial_basis(
      Polynomials::PolynomialsHermite::generate_complete_basis(regularity));

    std::vector<unsigned int> renumber =
      internal::hermite_hierarchic_to_lexicographic_numbering<dim>(regularity);
    polynomial_basis.set_numbering(renumber);

    return polynomial_basis;
  }



  /**
   * The @p Rescaler class implements the re-scaling of individual shape
   * functions required by Hermite bases on non-uniform meshes. The three
   * cases for different element dimensions are all defined separately
   * due to the requirement for different levels of nesting of for loops.
   */
  class Rescaler
  {
  public:
    template <int spacedim, typename Number>
    void
    rescale_fe_hermite_values(
      const FE_Hermite<1, spacedim>                         &fe_herm,
      const typename Mapping<1, spacedim>::InternalDataBase &mapping_data,
      Table<2, Number>                                      &value_list)
    {
      double cell_extent = 1.0;

      // Check mapping_data is associated with a compatible mapping class
      if (dynamic_cast<const typename MappingCartesian<1>::InternalData *>(
            &mapping_data) != nullptr)
        {
          const typename MappingCartesian<1>::InternalData *data =
            dynamic_cast<const typename MappingCartesian<1>::InternalData *>(
              &mapping_data);
          cell_extent = data->cell_extents[0];
        }
      else
        DEAL_II_ASSERT_UNREACHABLE();

      const unsigned int regularity      = fe_herm.get_regularity();
      const unsigned int n_dofs_per_cell = fe_herm.n_dofs_per_cell();
      const unsigned int n_q_points_out  = value_list.size(1);
      (void)n_dofs_per_cell;

      AssertDimension(value_list.size(0), n_dofs_per_cell);
      AssertDimension(n_dofs_per_cell, 2 * regularity + 2);

      std::vector<unsigned int> l2h =
        hermite_lexicographic_to_hierarchic_numbering<1>(regularity);

      for (unsigned int q = 0; q < n_q_points_out; ++q)
        {
          double factor_1 = 1.0;

          for (unsigned int d1 = 0, d2 = regularity + 1; d2 < n_dofs_per_cell;
               ++d1, ++d2)
            {
              /*
               * d1 is used to count over indices on the left and d2 counts
               * over indices on the right. These variables are used
               * to avoid the need to loop over vertices.
               */
              value_list(l2h[d1], q) *= factor_1;
              value_list(l2h[d2], q) *= factor_1;

              factor_1 *= cell_extent;
            }
        }
    }



    template <int spacedim, typename Number>
    void
    rescale_fe_hermite_values(
      const FE_Hermite<2, spacedim>                         &fe_herm,
      const typename Mapping<2, spacedim>::InternalDataBase &mapping_data,
      Table<2, Number>                                      &value_list)
    {
      Tensor<1, 2> cell_extents;

      // Check mapping_data is associated with a compatible mapping class
      if (dynamic_cast<const typename MappingCartesian<2>::InternalData *>(
            &mapping_data) != nullptr)
        {
          const typename MappingCartesian<2>::InternalData *data =
            dynamic_cast<const typename MappingCartesian<2>::InternalData *>(
              &mapping_data);
          cell_extents = data->cell_extents;
        }
      else
        DEAL_II_ASSERT_UNREACHABLE();

      const unsigned int regularity      = fe_herm.get_regularity();
      const unsigned int n_dofs_per_cell = fe_herm.n_dofs_per_cell();
      const unsigned int n_dofs_per_dim  = 2 * regularity + 2;
      const unsigned int n_q_points_out  = value_list.size(1);
      (void)n_dofs_per_cell;

      AssertDimension(value_list.size(0), n_dofs_per_cell);
      AssertDimension(n_dofs_per_dim * n_dofs_per_dim, n_dofs_per_cell);

      std::vector<unsigned int> l2h =
        hermite_lexicographic_to_hierarchic_numbering<2>(regularity);

      AssertDimension(l2h.size(), n_dofs_per_cell);

      for (unsigned int q = 0; q < n_q_points_out; ++q)
        {
          double factor_2 = 1.0;

          for (unsigned int d3 = 0, d4 = regularity + 1; d4 < n_dofs_per_dim;
               ++d3, ++d4)
            {
              double factor_1 = factor_2;

              for (unsigned int d1 = 0, d2 = regularity + 1;
                   d2 < n_dofs_per_dim;
                   ++d1, ++d2)
                {
                  /*
                   * d1 and d2 represent "left" and "right" in the
                   * x-direction, d3 and d4 represent "bottom" and "top"
                   * in the y-direction. As before, this is to avoid looping
                   * over vertices.
                   */
                  value_list(l2h[d1 + d3 * n_dofs_per_dim], q) *= factor_1;
                  value_list(l2h[d2 + d3 * n_dofs_per_dim], q) *= factor_1;
                  value_list(l2h[d1 + d4 * n_dofs_per_dim], q) *= factor_1;
                  value_list(l2h[d2 + d4 * n_dofs_per_dim], q) *= factor_1;

                  factor_1 *= cell_extents[0];
                }

              factor_2 *= cell_extents[1];
            }
        }
    }



    template <int spacedim, typename Number>
    void
    rescale_fe_hermite_values(
      const FE_Hermite<3, spacedim>                         &fe_herm,
      const typename Mapping<3, spacedim>::InternalDataBase &mapping_data,
      Table<2, Number>                                      &value_list)
    {
      Tensor<1, 3> cell_extents;

      // Check mapping_data is associated with a compatible mapping class
      if (dynamic_cast<const typename MappingCartesian<3>::InternalData *>(
            &mapping_data) != nullptr)
        {
          const typename MappingCartesian<3>::InternalData *data =
            dynamic_cast<const typename MappingCartesian<3>::InternalData *>(
              &mapping_data);
          cell_extents = data->cell_extents;
        }
      else
        DEAL_II_ASSERT_UNREACHABLE();

      const unsigned int regularity      = fe_herm.get_regularity();
      const unsigned int n_dofs_per_cell = fe_herm.n_dofs_per_cell();
      const unsigned int n_dofs_per_dim  = 2 * regularity + 2;
      const unsigned int n_dofs_per_quad = n_dofs_per_dim * n_dofs_per_dim;
      const unsigned int n_q_points_out  = value_list.size(1);
      (void)n_dofs_per_cell;

      AssertDimension(value_list.size(0), n_dofs_per_cell);
      AssertDimension(Utilities::pow(n_dofs_per_dim, 3), n_dofs_per_cell);

      std::vector<unsigned int> l2h =
        hermite_lexicographic_to_hierarchic_numbering<3>(regularity);

      for (unsigned int q = 0; q < n_q_points_out; ++q)
        {
          double factor_3 = 1.0;

          for (unsigned int d5 = 0, d6 = regularity + 1; d6 < n_dofs_per_dim;
               ++d5, ++d6)
            {
              double factor_2 = factor_3;

              for (unsigned int d3 = 0, d4 = regularity + 1;
                   d4 < n_dofs_per_dim;
                   ++d3, ++d4)
                {
                  double factor_1 = factor_2;

                  for (unsigned int d1 = 0, d2 = regularity + 1;
                       d2 < n_dofs_per_dim;
                       ++d1, ++d2)
                    {
                      /*
                       * d1, d2: "left" and "right" (x-direction)
                       * d3, d4: "bottom" and "top" (y-direction)
                       * d5, d6: "down" and "up"    (z-direction)
                       * This avoids looping over vertices
                       */
                      value_list(
                        l2h[d1 + d3 * n_dofs_per_dim + d5 * n_dofs_per_quad],
                        q) *= factor_1;
                      value_list(
                        l2h[d2 + d3 * n_dofs_per_dim + d5 * n_dofs_per_quad],
                        q) *= factor_1;
                      value_list(
                        l2h[d1 + d4 * n_dofs_per_dim + d5 * n_dofs_per_quad],
                        q) *= factor_1;
                      value_list(
                        l2h[d2 + d4 * n_dofs_per_dim + d5 * n_dofs_per_quad],
                        q) *= factor_1;
                      value_list(
                        l2h[d1 + d3 * n_dofs_per_dim + d6 * n_dofs_per_quad],
                        q) *= factor_1;
                      value_list(
                        l2h[d2 + d3 * n_dofs_per_dim + d6 * n_dofs_per_quad],
                        q) *= factor_1;
                      value_list(
                        l2h[d1 + d4 * n_dofs_per_dim + d6 * n_dofs_per_quad],
                        q) *= factor_1;
                      value_list(
                        l2h[d2 + d4 * n_dofs_per_dim + d6 * n_dofs_per_quad],
                        q) *= factor_1;

                      factor_1 *= cell_extents[0];
                    }

                  factor_2 *= cell_extents[1];
                }

              factor_3 *= cell_extents[2];
            }
        }
    }
  }; // class Rescaler
} // namespace internal



// Constructors
template <int dim, int spacedim>
FE_Hermite<dim, spacedim>::FE_Hermite(const unsigned int fe_degree)
  : FE_Poly<dim, spacedim>(
      internal::get_hermite_polynomials<dim>(fe_degree),
      FiniteElementData<dim>(internal::get_hermite_dpo_vector(
                               dim,
                               internal::get_regularity_from_degree(fe_degree)),
                             1,
                             std::max(1U, fe_degree),
                             ((fe_degree > 2) ? FiniteElementData<dim>::H2 :
                                                FiniteElementData<dim>::H1)),
      std::vector<bool>(Utilities::pow(std::max(2U, fe_degree + 1), dim),
                        false),
      std::vector<ComponentMask>(Utilities::pow(std::max(2U, fe_degree + 1),
                                                dim),
                                 ComponentMask(1, true)))
  , regularity(internal::get_regularity_from_degree(fe_degree))
{
  Assert((fe_degree % 2 == 1),
         ExcMessage(
           "ERROR: The current implementation of Hermite interpolation "
           "polynomials is only defined for odd polynomial degrees. Running "
           "in release mode will use a polynomial degree of max(1,fe_degree-1) "
           "to protect against unexpected internal bugs."));
}



template <int dim, int spacedim>
std::string
FE_Hermite<dim, spacedim>::get_name() const
{
  std::ostringstream name_buffer;
  name_buffer << "FE_Hermite<" << Utilities::dim_string(dim, spacedim) << ">("
              << this->degree << ")";
  return name_buffer.str();
}



template <int dim, int spacedim>
std::unique_ptr<FiniteElement<dim, spacedim>>
FE_Hermite<dim, spacedim>::clone() const
{
  return std::make_unique<FE_Hermite<dim, spacedim>>(*this);
}



template <int dim, int spacedim>
UpdateFlags
FE_Hermite<dim, spacedim>::requires_update_flags(const UpdateFlags flags) const
{
  UpdateFlags out = FE_Poly<dim, spacedim>::requires_update_flags(flags);
  if (flags & (update_values | update_gradients | update_hessians |
               update_3rd_derivatives))
    out |= update_rescale; // since we need to rescale values, gradients, ...
  return out;
}



/**
 * A large part of the following function is copied from FE_Q_Base. At present
 * the case of two Hermite bases meeting is not implemented as it is unlikely
 * that two different Hermite bases would be used in an hp method. This may be
 * implemented in a later update.
 */
template <int dim, int spacedim>
std::vector<std::pair<unsigned int, unsigned int>>
FE_Hermite<dim, spacedim>::hp_vertex_dof_identities(
  const FiniteElement<dim, spacedim> &fe_other) const
{
  if (dynamic_cast<const FE_Q_Base<dim, spacedim> *>(&fe_other) != nullptr)
    {
      // there should be exactly one single DoF of FE_Q_Base at a vertex, and it
      // should have an identical value to the first Hermite DoF
      return {{0U, 0U}};
    }
  else if (dynamic_cast<const FE_SimplexP<dim, spacedim> *>(&fe_other) !=
           nullptr)
    {
      // there should be exactly one single DoF of FE_Q_Base at a vertex, and it
      // should have an identical value to the first Hermite DoF
      return {{0U, 0U}};
    }
  else if (dynamic_cast<const FE_Nothing<dim> *>(&fe_other) != nullptr)
    {
      // the FE_Nothing has no degrees of freedom, so there are no
      // equivalencies to be recorded
      return {};
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
      return {};
    }
  else if (const FE_Hermite<dim, spacedim> *fe_herm_other =
             dynamic_cast<const FE_Hermite<dim, spacedim> *>(&fe_other))
    {
      DEAL_II_NOT_IMPLEMENTED();
      return {};
    }
  else
    {
      DEAL_II_NOT_IMPLEMENTED();
      return {};
    }
}



/**
 * This function only supplies empty lists of pairs, since Hermite
 * stores all DoFs on vertices meaning there is no continuity that
 * could be enforced.
 */
template <int dim, int spacedim>
std::vector<std::pair<unsigned int, unsigned int>>
FE_Hermite<dim, spacedim>::hp_line_dof_identities(
  const FiniteElement<dim, spacedim> &fe_other) const
{
  (void)fe_other;
  return {};
}



/**
 * Similar to above, no continuity can be enforced on quads
 * for Hermite.
 */
template <int dim, int spacedim>
std::vector<std::pair<unsigned int, unsigned int>>
FE_Hermite<dim, spacedim>::hp_quad_dof_identities(
  const FiniteElement<dim, spacedim> &fe_other,
  const unsigned int                  face_no) const
{
  (void)fe_other;
  (void)face_no;
  return {};
}



/*
 * The layout of this function is largely copied directly from FE_Q,
 * however FE_Hermite can behave significantly differently in terms
 * of domination due to how the function space is defined */
template <int dim, int spacedim>
FiniteElementDomination::Domination
FE_Hermite<dim, spacedim>::compare_for_domination(
  const FiniteElement<dim, spacedim> &fe_other,
  const unsigned int                  codim) const
{
  Assert(codim <= dim, ExcImpossibleInDim(dim));

  if (codim > 0)
    if (dynamic_cast<const FE_DGQ<dim, spacedim> *>(&fe_other) != nullptr)
      // there are no requirements between continuous and discontinuous elements
      return FiniteElementDomination::no_requirements;


  // vertex/line/face domination
  // (if fe_other is not derived from FE_DGQ)
  // & cell domination
  // ----------------------------------------
  if (const FE_Hermite<dim, spacedim> *fe_hermite_other =
        dynamic_cast<const FE_Hermite<dim, spacedim> *>(&fe_other))
    {
      if (this->degree < fe_hermite_other->degree)
        return FiniteElementDomination::this_element_dominates;
      else if (this->degree == fe_hermite_other->degree)
        return FiniteElementDomination::either_element_can_dominate;
      else
        return FiniteElementDomination::other_element_dominates;
    }
  if (const FE_Q<dim, spacedim> *fe_q_other =
        dynamic_cast<const FE_Q<dim, spacedim> *>(&fe_other))
    {
      if (fe_q_other->degree == 1)
        {
          if (this->degree == 1)
            return FiniteElementDomination::either_element_can_dominate;
          else
            return FiniteElementDomination::other_element_dominates;
        }
      else if (this->degree <= fe_q_other->degree)
        return FiniteElementDomination::this_element_dominates;
      else
        return FiniteElementDomination::neither_element_dominates;
    }
  else if (const FE_SimplexP<dim, spacedim> *fe_p_other =
             dynamic_cast<const FE_SimplexP<dim, spacedim> *>(&fe_other))
    {
      if (fe_p_other->degree == 1)
        {
          if (this->degree == 1)
            return FiniteElementDomination::either_element_can_dominate;
          else
            return FiniteElementDomination::other_element_dominates;
        }
      else if (this->degree <= fe_p_other->degree)
        return FiniteElementDomination::this_element_dominates;
      else
        return FiniteElementDomination::neither_element_dominates;
    }
  else if (const FE_WedgeP<dim, spacedim> *fe_wp_other =
             dynamic_cast<const FE_WedgeP<dim, spacedim> *>(&fe_other))
    {
      if (fe_wp_other->degree == 1)
        {
          if (this->degree == 1)
            return FiniteElementDomination::either_element_can_dominate;
          else
            return FiniteElementDomination::other_element_dominates;
        }
      else if (this->degree <= fe_wp_other->degree)
        return FiniteElementDomination::this_element_dominates;
      else
        return FiniteElementDomination::neither_element_dominates;
    }
  else if (const FE_PyramidP<dim, spacedim> *fe_pp_other =
             dynamic_cast<const FE_PyramidP<dim, spacedim> *>(&fe_other))
    {
      if (fe_pp_other->degree == 1)
        {
          if (this->degree == 1)
            return FiniteElementDomination::either_element_can_dominate;
          else
            return FiniteElementDomination::other_element_dominates;
        }
      else if (this->degree <= fe_pp_other->degree)
        return FiniteElementDomination::this_element_dominates;
      else
        return FiniteElementDomination::neither_element_dominates;
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
std::vector<unsigned int>
FE_Hermite<dim, spacedim>::get_lexicographic_to_hierarchic_numbering() const
{
  return internal::hermite_lexicographic_to_hierarchic_numbering<dim>(
    this->regularity);
}



template <int dim, int spacedim>
Table<2, unsigned int>
FE_Hermite<dim, spacedim>::get_dofs_corresponding_to_outward_normal_derivatives(
  const unsigned int derivative_order) const
{
  /*
   * Create a look-up table for finding relevant dofs on all
   * 2*dim faces of reference cell
   */
  const unsigned int degree        = this->degree;
  const unsigned int regularity    = this->get_regularity();
  const unsigned int dofs_per_face = this->n_dofs_per_face();
  AssertIndexRange(derivative_order, regularity + 1);
  AssertDimension(dofs_per_face,
                  (regularity + 1) * Utilities::pow(degree + 1, dim - 1));

  const unsigned int relevant_dofs_per_face = dofs_per_face / (regularity + 1);
  Table<2, unsigned int> dofs_on_each_face(2 * dim, relevant_dofs_per_face);

  /*
   * Use knowledge of the local degree numbering for this version,
   * saving expensive calls to reinit().
   */
  const std::vector<unsigned int> l2h =
    get_lexicographic_to_hierarchic_numbering();
  const unsigned int dofs_per_cell = Utilities::pow(degree + 1, dim);
  AssertDimension(dofs_per_cell, l2h.size());
  (void)dofs_per_cell;

  /*
   * The following loop uses the variables batch_size, batch_index
   * and local_index to simplify calculations. The idea is to find
   * relevant DoFs in batches, with each batch representing a
   * set of DoFs of interest on a given face that occur consecutively
   * in the ordering of all DoFs on the reference cell.
   * To quickly summarise what the variables mean:
   * sublist_index: index of a DoF in the list of relevant DoF indices
   * index: index of a DoF in the list of all DoFs on the cell
   * batch_size: Number of consecutive DoFs in the ordering that are
   *             all of interest,
   * batch index: Index of the current batch in the list of batches
   * local_index: Index of the current DoF within a batch
   *
   * The variable correction is used because the pattern of relevant
   * DoFs on opposite face pairs is always the same, just separated by
   * a constant offset value in the indices, so it's easier to calculate
   * the pattern once and find this offset value.
   */
  for (unsigned int d = 0, batch_size = 1; d < dim;
       ++d, batch_size *= degree + 1)
    for (unsigned int sublist_index = 0; sublist_index < relevant_dofs_per_face;
         ++sublist_index)
      {
        const unsigned int local_index = sublist_index % batch_size;
        const unsigned int batch_index = sublist_index / batch_size;

        unsigned int index =
          local_index +
          (batch_index * (degree + 1) + derivative_order) * batch_size;
        unsigned int correction = batch_size * (regularity + 1);
        Assert(index + correction < dofs_per_cell,
               ExcDimensionMismatch(index + correction, dofs_per_cell));

        dofs_on_each_face(2 * d, sublist_index)     = l2h[index];
        dofs_on_each_face(2 * d + 1, sublist_index) = l2h[index + correction];
      }

  return dofs_on_each_face;
}



template <int dim, int spacedim>
void
FE_Hermite<dim, spacedim>::fill_fe_values(
  const typename Triangulation<dim, spacedim>::cell_iterator &,
  const CellSimilarity::Similarity cell_similarity,
  const Quadrature<dim> & /*quadrature*/,
  const Mapping<dim, spacedim>                            &mapping,
  const typename Mapping<dim, spacedim>::InternalDataBase &mapping_internal,
  const dealii::internal::FEValuesImplementation::MappingRelatedData<dim,
                                                                     spacedim>
    & /*mapping_data*/,
  const typename FiniteElement<dim, spacedim>::InternalDataBase &fe_internal,
  dealii::internal::FEValuesImplementation::FiniteElementRelatedData<dim,
                                                                     spacedim>
    &output_data) const
{
  // Convert data object to internal data for this class.
  // Fails with an exception if that is not possible.
  Assert(
    (dynamic_cast<const typename FE_Hermite<dim, spacedim>::InternalData *>(
       &fe_internal) != nullptr),
    ExcInternalError());
  const typename FE_Hermite<dim, spacedim>::InternalData &fe_data =
    static_cast<const typename FE_Hermite<dim, spacedim>::InternalData &>(
      fe_internal);

  const UpdateFlags flags(fe_data.update_each);

  // Transform values, gradients and higher derivatives. Values also need to
  // be rescaled according the the nodal derivative they correspond to.
  if ((flags & update_values) &&
      (cell_similarity != CellSimilarity::translation))
    {
      internal::Rescaler shape_fix;
      for (unsigned int i = 0; i < output_data.shape_values.size(0); ++i)
        for (unsigned int q = 0; q < output_data.shape_values.size(1); ++q)
          output_data.shape_values(i, q) = fe_data.shape_values(i, q);
      shape_fix.rescale_fe_hermite_values(*this,
                                          mapping_internal,
                                          output_data.shape_values);
    }

  if ((flags & update_gradients) &&
      (cell_similarity != CellSimilarity::translation))
    {
      for (unsigned int k = 0; k < this->n_dofs_per_cell(); ++k)
        mapping.transform(make_array_view(fe_data.shape_gradients, k),
                          mapping_covariant,
                          mapping_internal,
                          make_array_view(output_data.shape_gradients, k));

      internal::Rescaler grad_fix;
      grad_fix.rescale_fe_hermite_values(*this,
                                         mapping_internal,
                                         output_data.shape_gradients);
    }

  if ((flags & update_hessians) &&
      (cell_similarity != CellSimilarity::translation))
    {
      for (unsigned int k = 0; k < this->n_dofs_per_cell(); ++k)
        mapping.transform(make_array_view(fe_data.shape_hessians, k),
                          mapping_covariant_gradient,
                          mapping_internal,
                          make_array_view(output_data.shape_hessians, k));

      internal::Rescaler hessian_fix;
      hessian_fix.rescale_fe_hermite_values(*this,
                                            mapping_internal,
                                            output_data.shape_hessians);
    }

  if ((flags & update_3rd_derivatives) &&
      (cell_similarity != CellSimilarity::translation))
    {
      for (unsigned int k = 0; k < this->n_dofs_per_cell(); ++k)
        mapping.transform(make_array_view(fe_data.shape_3rd_derivatives, k),
                          mapping_covariant_hessian,
                          mapping_internal,
                          make_array_view(output_data.shape_3rd_derivatives,
                                          k));

      internal::Rescaler third_dev_fix;
      third_dev_fix.rescale_fe_hermite_values(
        *this, mapping_internal, output_data.shape_3rd_derivatives);
    }
}



template <int dim, int spacedim>
void
FE_Hermite<dim, spacedim>::fill_fe_face_values(
  const typename Triangulation<dim, spacedim>::cell_iterator &cell,
  const unsigned int                                          face_no,
  const hp::QCollection<dim - 1>                             &quadrature,
  const Mapping<dim, spacedim>                               &mapping,
  const typename Mapping<dim, spacedim>::InternalDataBase    &mapping_internal,
  const dealii::internal::FEValuesImplementation::MappingRelatedData<dim,
                                                                     spacedim>
    &,
  const typename FiniteElement<dim, spacedim>::InternalDataBase &fe_internal,
  dealii::internal::FEValuesImplementation::FiniteElementRelatedData<dim,
                                                                     spacedim>
    &output_data) const
{
  /*
   * Convert data object to internal data for this class. Fails with
   * an exception if that is not possible.
   */
  Assert(
    (dynamic_cast<const typename FE_Hermite<dim, spacedim>::InternalData *>(
       &fe_internal) != nullptr),
    ExcInternalError());
  const typename FE_Hermite<dim, spacedim>::InternalData &fe_data =
    static_cast<const typename FE_Hermite<dim, spacedim>::InternalData &>(
      fe_internal);

  Assert((dynamic_cast<
            const typename MappingCartesian<dim, spacedim>::InternalData *>(
            &mapping_internal) != nullptr),
         ExcInternalError());

  AssertDimension(quadrature.size(), 1U);

  /*
   * offset determines which data set to take (all data sets for all
   * faces are stored contiguously)
   */
  const typename QProjector<dim>::DataSetDescriptor offset =
    QProjector<dim>::DataSetDescriptor::face(
      ReferenceCells::get_hypercube<dim>(),
      face_no,
      cell->combined_face_orientation(face_no),
      quadrature[0].size());

  const UpdateFlags flags(fe_data.update_each);

  // Transform values, gradients and higher derivatives.
  if (flags & update_values)
    {
      for (unsigned int k = 0; k < this->dofs_per_cell; ++k)
        for (unsigned int i = 0; i < quadrature[0].size(); ++i)
          output_data.shape_values(k, i) = fe_data.shape_values[k][i + offset];

      internal::Rescaler shape_face_fix;
      shape_face_fix.rescale_fe_hermite_values(*this,
                                               mapping_internal,
                                               output_data.shape_values);
    }

  if (flags & update_gradients)
    {
      for (unsigned int k = 0; k < this->dofs_per_cell; ++k)
        mapping.transform(make_array_view(fe_data.shape_gradients,
                                          k,
                                          offset,
                                          quadrature[0].size()),
                          mapping_covariant,
                          mapping_internal,
                          make_array_view(output_data.shape_gradients, k));

      internal::Rescaler grad_face_fix;
      grad_face_fix.rescale_fe_hermite_values(*this,
                                              mapping_internal,
                                              output_data.shape_gradients);
    }

  if (flags & update_hessians)
    {
      for (unsigned int k = 0; k < this->dofs_per_cell; ++k)
        mapping.transform(make_array_view(fe_data.shape_hessians,
                                          k,
                                          offset,
                                          quadrature[0].size()),
                          mapping_covariant_gradient,
                          mapping_internal,
                          make_array_view(output_data.shape_hessians, k));

      internal::Rescaler hessian_face_fix;
      hessian_face_fix.rescale_fe_hermite_values(*this,
                                                 mapping_internal,
                                                 output_data.shape_hessians);
    }

  if (flags & update_3rd_derivatives)
    {
      for (unsigned int k = 0; k < this->dofs_per_cell; ++k)
        mapping.transform(make_array_view(fe_data.shape_3rd_derivatives,
                                          k,
                                          offset,
                                          quadrature[0].size()),
                          mapping_covariant_hessian,
                          mapping_internal,
                          make_array_view(output_data.shape_3rd_derivatives,
                                          k));

      internal::Rescaler shape_3rd_face_fix;
      shape_3rd_face_fix.rescale_fe_hermite_values(
        *this, mapping_internal, output_data.shape_3rd_derivatives);
    }
}



// Explicit instantiations
#include "fe/fe_hermite.inst"

DEAL_II_NAMESPACE_CLOSE
