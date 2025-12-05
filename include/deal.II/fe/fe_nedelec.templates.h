// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_fe_nedelec_templates_h
#define dealii_fe_nedelec_templates_h

#include <deal.II/fe/fe_nedelec.h>

DEAL_II_NAMESPACE_OPEN

template <int dim>
template <typename Number, int dim_q>
void
FE_NedelecNodal<dim>::fill_shape_info(
  internal::MatrixFreeFunctions::ShapeInfo<Number> *shape_info,
  const Quadrature<dim_q>                          &quad_in,
  const unsigned int                                base_element_number) const
{
  const auto quad = quad_in.get_tensor_basis()[0];
  shape_info->element_type =
    internal::MatrixFreeFunctions::ElementType::tensor_nedelec;

  const FiniteElement<dim, dim> &fe = this->base_element(base_element_number);
  shape_info->n_dimensions          = dim;
  shape_info->n_components          = this->n_components();

  shape_info->data.resize(2);
  const unsigned int n_q_points_1d = quad.size();

  shape_info->n_q_points      = Utilities::fixed_power<dim>(n_q_points_1d);
  shape_info->n_q_points_face = Utilities::fixed_power<dim - 1>(n_q_points_1d);

  shape_info->dofs_per_component_on_cell =
    this->n_dofs_per_cell() / this->n_components();

  // NOTE dofs_per_component_on_face is in tangential direction!
  shape_info->dofs_per_component_on_face =
    this->n_dofs_per_face() + Utilities::pow(this->degree, dim - 2);
  const unsigned int dofs_per_face = this->n_dofs_per_face();

  const unsigned int dofs_per_line = this->n_dofs_per_line();

  // degrees of freedom on (dim-2)-dimensional elements
  const unsigned int dofs_dim_2 =
    dim == 3 ? this->n_dofs_per_line() : this->n_dofs_per_vertex();


  shape_info->lexicographic_numbering =
    get_lexicographic_numbering(this->degree - 1);


  // To get the right shape_values of the Nedelec element
  std::vector<unsigned int> lex_normal, lex_tangent;
  for (unsigned int i = dofs_per_line * 2; i < dofs_per_line * 2 + fe.degree;
       ++i)
    lex_normal.push_back(i);
  lex_tangent.push_back(lex_normal[0]);
  for (unsigned int i = 4 * dofs_per_face - 4 * dofs_dim_2;
       i < 4 * dofs_per_face - 4 * dofs_dim_2 + (fe.degree - 1) * fe.degree;
       i += fe.degree)
    lex_tangent.push_back(i);
  lex_tangent.push_back(dofs_per_line * 3);


  // 'direction' distinguishes between normal=0 and tangential=1 direction
  for (unsigned int direction = 0; direction < 2; ++direction)
    {
      shape_info->data[direction].element_type =
        internal::MatrixFreeFunctions::ElementType::tensor_nedelec;
      shape_info->data[direction].quadrature    = quad;
      shape_info->data[direction].n_q_points_1d = n_q_points_1d;
      shape_info->data[direction].fe_degree     = fe.degree - (1 - direction);
      const std::vector<unsigned int> &lexicographic =
        direction == 0 ? lex_normal : lex_tangent;


      shape_info->data[direction].evaluate_shape_functions(fe,
                                                           quad,
                                                           lexicographic,
                                                           direction);
      shape_info->data[direction].evaluate_collocation_space(fe,
                                                             quad,
                                                             lexicographic,
                                                             direction);

      shape_info->data[direction].check_and_set_shapes_symmetric();
    }
}

DEAL_II_NAMESPACE_CLOSE

#endif
