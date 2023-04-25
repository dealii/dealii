// ---------------------------------------------------------------------
//
// Copyright (C) 2015 - 2022 by the deal.II authors
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


#include <deal.II/fe/fe_nedelec_sz.h>

#include <memory>

DEAL_II_NAMESPACE_OPEN

// Constructor:
template <int dim, int spacedim>
FE_NedelecSZ<dim, spacedim>::FE_NedelecSZ(const unsigned int order)
  : FiniteElement<dim, dim>(
      FiniteElementData<dim>(get_dpo_vector(order),
                             dim,
                             order + 1,
                             FiniteElementData<dim>::Hcurl),
      std::vector<bool>(compute_num_dofs(order), true),
      std::vector<ComponentMask>(compute_num_dofs(order),
                                 std::vector<bool>(dim, true)))
{
  Assert(dim >= 2, ExcImpossibleInDim(dim));

  this->mapping_kind = mapping_nedelec;
  // Set up the table converting components to base components. Since we have
  // only one base element, everything remains zero except the component in the
  // base, which is the component itself.
  for (unsigned int comp = 0; comp < this->n_components(); ++comp)
    {
      this->component_to_base_table[comp].first.second = comp;
    }

  // Generate the 1-D polynomial basis.
  create_polynomials(order);
}



// Shape functions:
template <int dim, int spacedim>
double
FE_NedelecSZ<dim, spacedim>::shape_value(const unsigned int /*i*/,
                                         const Point<dim> & /*p*/) const
{
  Assert(false, (typename FiniteElement<dim, spacedim>::ExcFENotPrimitive()));
  return 0.;
}



template <int dim, int spacedim>
double
FE_NedelecSZ<dim, spacedim>::shape_value_component(
  const unsigned int /*i*/,
  const Point<dim> & /*p*/,
  const unsigned int /*component*/) const
{
  // Not implemented yet:
  Assert(false, ExcNotImplemented());
  return 0.;
}



template <int dim, int spacedim>
Tensor<1, dim>
FE_NedelecSZ<dim, spacedim>::shape_grad(const unsigned int /*i*/,
                                        const Point<dim> & /*p*/) const
{
  Assert(false, (typename FiniteElement<dim, spacedim>::ExcFENotPrimitive()));
  return Tensor<1, dim>();
}



template <int dim, int spacedim>
Tensor<1, dim>
FE_NedelecSZ<dim, spacedim>::shape_grad_component(
  const unsigned int /*i*/,
  const Point<dim> & /*p*/,
  const unsigned int /*component*/) const
{
  Assert(false, ExcNotImplemented());
  return Tensor<1, dim>();
}



template <int dim, int spacedim>
Tensor<2, dim>
FE_NedelecSZ<dim, spacedim>::shape_grad_grad(const unsigned int /*i*/,
                                             const Point<dim> & /*p*/) const
{
  Assert(false, (typename FiniteElement<dim, spacedim>::ExcFENotPrimitive()));
  return Tensor<2, dim>();
}



template <int dim, int spacedim>
Tensor<2, dim>
FE_NedelecSZ<dim, spacedim>::shape_grad_grad_component(
  const unsigned int /*i*/,
  const Point<dim> & /*p*/,
  const unsigned int /*component*/) const
{
  Assert(false, ExcNotImplemented());
  return Tensor<2, dim>();
}



template <int dim, int spacedim>
std::unique_ptr<typename dealii::FiniteElement<dim, spacedim>::InternalDataBase>
FE_NedelecSZ<dim, spacedim>::get_data(
  const UpdateFlags update_flags,
  const Mapping<dim, spacedim> & /*mapping*/,
  const Quadrature<dim> &quadrature,
  dealii::internal::FEValuesImplementation::FiniteElementRelatedData<dim,
                                                                     spacedim>
    & /*output_data*/) const
{
  std::unique_ptr<
    typename dealii::FiniteElement<dim, spacedim>::InternalDataBase>
        data_ptr   = std::make_unique<InternalData>();
  auto &data       = dynamic_cast<InternalData &>(*data_ptr);
  data.update_each = requires_update_flags(update_flags);

  // Useful quantities:
  const unsigned int degree(this->degree - 1); // Note: FE holds input degree+1

  const unsigned int vertices_per_cell = GeometryInfo<dim>::vertices_per_cell;
  const unsigned int lines_per_cell    = GeometryInfo<dim>::lines_per_cell;
  const unsigned int faces_per_cell    = GeometryInfo<dim>::faces_per_cell;

  const unsigned int n_line_dofs = this->n_dofs_per_line() * lines_per_cell;

  // we assume that all quads have the same number of dofs
  const unsigned int n_face_dofs = this->n_dofs_per_quad(0) * faces_per_cell;

  const UpdateFlags  flags(data.update_each);
  const unsigned int n_q_points = quadrature.size();

  // Resize the internal data storage:
  data.sigma_imj_values.resize(
    n_q_points,
    std::vector<std::vector<double>>(vertices_per_cell,
                                     std::vector<double>(vertices_per_cell)));

  data.sigma_imj_grads.resize(vertices_per_cell,
                              std::vector<std::vector<double>>(
                                vertices_per_cell, std::vector<double>(dim)));

  // Resize shape function arrays according to update flags:
  if (flags & update_values)
    {
      data.shape_values.resize(this->n_dofs_per_cell(),
                               std::vector<Tensor<1, dim>>(n_q_points));
    }

  if (flags & update_gradients)
    {
      data.shape_grads.resize(this->n_dofs_per_cell(),
                              std::vector<DerivativeForm<1, dim, dim>>(
                                n_q_points));
    }

  if (flags & update_hessians)
    {
      data.shape_hessians.resize(this->n_dofs_per_cell(),
                                 std::vector<DerivativeForm<2, dim, dim>>(
                                   n_q_points));
    }

  std::vector<Point<dim>> p_list(n_q_points);
  p_list = quadrature.get_points();


  switch (dim)
    {
      case 2:
        {
          // Compute values of sigma & lambda and the sigma differences and
          // lambda additions.
          std::vector<std::vector<double>> sigma(
            n_q_points, std::vector<double>(lines_per_cell));
          std::vector<std::vector<double>> lambda(
            n_q_points, std::vector<double>(lines_per_cell));

          for (unsigned int q = 0; q < n_q_points; ++q)
            {
              sigma[q][0] = (1.0 - p_list[q][0]) + (1.0 - p_list[q][1]);
              sigma[q][1] = p_list[q][0] + (1.0 - p_list[q][1]);
              sigma[q][2] = (1.0 - p_list[q][0]) + p_list[q][1];
              sigma[q][3] = p_list[q][0] + p_list[q][1];

              lambda[q][0] = (1.0 - p_list[q][0]) * (1.0 - p_list[q][1]);
              lambda[q][1] = p_list[q][0] * (1.0 - p_list[q][1]);
              lambda[q][2] = (1.0 - p_list[q][0]) * p_list[q][1];
              lambda[q][3] = p_list[q][0] * p_list[q][1];
              for (unsigned int i = 0; i < vertices_per_cell; ++i)
                {
                  for (unsigned int j = 0; j < vertices_per_cell; ++j)
                    {
                      data.sigma_imj_values[q][i][j] =
                        sigma[q][i] - sigma[q][j];
                    }
                }
            }

          // Calculate the gradient of sigma_imj_values[q][i][j] =
          // sigma[q][i]-sigma[q][j]
          //   - this depends on the component and the direction of the
          //   corresponding edge.
          //   - the direction of the edge is determined by
          //   sigma_imj_sign[i][j].
          // Helper arrays:
          const int sigma_comp_signs[GeometryInfo<2>::vertices_per_cell][2] = {
            {-1, -1}, {1, -1}, {-1, 1}, {1, 1}};
          int          sigma_imj_sign[vertices_per_cell][vertices_per_cell];
          unsigned int sigma_imj_component[vertices_per_cell]
                                          [vertices_per_cell];

          for (unsigned int i = 0; i < vertices_per_cell; ++i)
            {
              for (unsigned int j = 0; j < vertices_per_cell; ++j)
                {
                  // sigma_imj_sign is the sign (+/-) of the coefficient of
                  // x/y/z in sigma_imj_values Due to the numbering of vertices
                  // on the reference element it is easy to find edges in the
                  // positive direction are from smaller to higher local vertex
                  // numbering.
                  sigma_imj_sign[i][j] = (i < j) ? -1 : 1;
                  sigma_imj_sign[i][j] = (i == j) ? 0 : sigma_imj_sign[i][j];

                  // Now store the component which the sigma_i - sigma_j
                  // corresponds to:
                  sigma_imj_component[i][j] = 0;
                  for (unsigned int d = 0; d < dim; ++d)
                    {
                      int temp_imj =
                        sigma_comp_signs[i][d] - sigma_comp_signs[j][d];
                      // Only interested in the first non-zero
                      // as if there is a second, it can not be a valid edge.
                      if (temp_imj != 0)
                        {
                          sigma_imj_component[i][j] = d;
                          break;
                        }
                    }
                  // Can now calculate the gradient, only non-zero in the
                  // component given: Note some i,j combinations will be
                  // incorrect, but only on invalid edges.
                  data.sigma_imj_grads[i][j][sigma_imj_component[i][j]] =
                    2.0 * sigma_imj_sign[i][j];
                }
            }

          // Now compute the edge parameterisations for a single element
          // with global numbering matching that of the reference element:

          // Resize the edge parameterisations
          data.edge_sigma_values.resize(lines_per_cell);
          data.edge_sigma_grads.resize(lines_per_cell);
          for (unsigned int m = 0; m < lines_per_cell; ++m)
            {
              data.edge_sigma_values[m].resize(n_q_points);

              // sigma grads are constant in a cell (no need for quad points)
              data.edge_sigma_grads[m].resize(dim);
            }

          // Fill the values for edge lambda and edge sigma:
          const unsigned int
            edge_sigma_direction[GeometryInfo<2>::lines_per_cell] = {1,
                                                                     1,
                                                                     0,
                                                                     0};

          data.edge_lambda_values.resize(lines_per_cell,
                                         std::vector<double>(n_q_points));
          data.edge_lambda_grads_2d.resize(lines_per_cell,
                                           std::vector<double>(dim));
          for (unsigned int m = 0; m < lines_per_cell; ++m)
            {
              // e1=max(reference vertex numbering on this edge)
              // e2=min(reference vertex numbering on this edge)
              // Which is guaranteed to be:
              const unsigned int e1(
                GeometryInfo<dim>::line_to_cell_vertices(m, 1));
              const unsigned int e2(
                GeometryInfo<dim>::line_to_cell_vertices(m, 0));
              for (unsigned int q = 0; q < n_q_points; ++q)
                {
                  data.edge_sigma_values[m][q] =
                    data.sigma_imj_values[q][e2][e1];
                  data.edge_lambda_values[m][q] = lambda[q][e1] + lambda[q][e2];
                }

              data.edge_sigma_grads[m][edge_sigma_direction[m]] = -2.0;
            }
          data.edge_lambda_grads_2d[0] = {-1.0, 0.0};
          data.edge_lambda_grads_2d[1] = {1.0, 0.0};
          data.edge_lambda_grads_2d[2] = {0.0, -1.0};
          data.edge_lambda_grads_2d[3] = {0.0, 1.0};

          // If the polynomial order is 0, then no more work to do:
          if (degree < 1)
            {
              break;
            }

          // Otherwise, we can compute the non-cell dependent shape functions.
          //
          // Note: the local dof numberings follow the usual order of lines ->
          // faces -> cells
          //       (we have no vertex-based DoFs in this element).
          // For a given cell we have:
          //      n_line_dofs = dofs_per_line*lines_per_cell.
          //      n_face_dofs = dofs_per_face*faces_per_cell.
          //      n_cell_dofs = dofs_per_quad (2d)
          //                  = dofs_per_hex (3d)
          //
          // i.e. For the local dof numbering:
          //      the first line dof is 0,
          //      the first face dof is n_line_dofs,
          //      the first cell dof is n_line_dofs + n_face_dofs.
          //
          // On a line, DoFs are ordered first by line_dof and then line_index:
          // i.e. line_dof_index = line_dof + line_index*(dofs_per_line)
          //
          // and similarly for faces:
          // i.e. face_dof_index = face_dof + face_index*(dofs_per_face).
          //
          // HOWEVER, we have different types of DoFs on a line/face/cell.
          // On a line we have two types, lowest order and higher order
          // gradients.
          //    - The numbering is such the lowest order is first, then higher
          //    order.
          //      This is simple enough as there is only 1 lowest order and
          //      degree higher orders DoFs per line.
          //
          // On a 2d cell, we have 3 types: Type 1/2/3:
          //    - The ordering done by type:
          //      - Type 1: 0 <= i1,j1 < degree. degree^2 in total.
          //        Numbered: ij1 = i1 + j1*(degree).        i.e. cell_dof_index
          //        = ij1.
          //      - Type 2: 0 <= i2,j2 < degree. degree^2 in total.
          //        Numbered: ij2 = i2 + j2*(degree).        i.e. cell_dof_index
          //        = degree^2 + ij2
          //      - Type 3: 0 <= i3 < 2*degree. 2*degree in total.
          //        Numbered: ij3 = i3.                      i.e. cell_dof_index
          //        =  2*(degree^2) + ij3.
          //
          // These then fit into the local dof numbering described above:
          // - local dof numberings are:
          //   line_dofs: local_dof = line_dof_index.    0 <= local_dof <
          //   dofs_per_line*lines_per_cell face_dofs: local_dof =
          //   n_line_dofs*lines_per_cell + face_dof_index. cell dofs: local_dof
          //   = n_lines_dof + n_face_dofs + cell_dof_index.
          //
          // The cell-based shape functions are:
          //
          // Type 1 (gradients):
          // \phi^{C_{1}}_{ij) = grad( L_{i+2}(2x-1)L_{j+2}(2y-1) ),
          //
          // 0 <= i,j < degree.
          //
          // NOTE: The derivative produced by IntegratedLegendrePolynomials does
          // not account for the
          //       (2*x-1) or (2*y-1) so we must take this into account when
          //       taking derivatives.
          const unsigned int cell_type1_offset = n_line_dofs;

          // Type 2:
          // \phi^{C_{2}}_{ij) = L'_{i+2}(2x-1) L_{j+2}(2y-1) \mathbf{e}_{x}
          //                     - L_{i+2}(2x-1) L'_{j+2}(2y-1) \mathbf{e}_{y},
          //
          // 0 <= i,j < degree.
          const unsigned int cell_type2_offset =
            cell_type1_offset + degree * degree;

          // Type 3 (two subtypes):
          // \phi^{C_{3}}_{j)        = L_{j+2}(2y-1) \mathbf{e}_{x}
          //
          // \phi^{C_{3}}_{j+degree) = L_{j+2}(2x-1) \mathbf{e}_{y}
          //
          // 0 <= j < degree
          const unsigned int cell_type3_offset1 =
            cell_type2_offset + degree * degree;
          const unsigned int cell_type3_offset2 = cell_type3_offset1 + degree;

          if (flags & (update_values | update_gradients | update_hessians))
            {
              // compute all points we must evaluate the 1d polynomials at:
              std::vector<Point<dim>> cell_points(n_q_points);
              for (unsigned int q = 0; q < n_q_points; ++q)
                {
                  for (unsigned int d = 0; d < dim; ++d)
                    {
                      cell_points[q][d] = 2.0 * p_list[q][d] - 1.0;
                    }
                }

              // Loop through quad points:
              for (unsigned int q = 0; q < n_q_points; ++q)
                {
                  // pre-compute values & required derivatives at this quad
                  // point (x,y): polyx = L_{i+2}(2x-1), polyy = L_{j+2}(2y-1),
                  //
                  // for each polyc[d], c=x,y, contains the d-th derivative with
                  // respect to the coordinate c.

                  // We only need poly values and 1st derivative for
                  // update_values, but need the 2nd derivative too for
                  // update_gradients. For update_hessians we also need the 3rd
                  // derivatives.
                  const unsigned int poly_length =
                    (flags & update_hessians) ?
                      4 :
                      ((flags & update_gradients) ? 3 : 2);

                  std::vector<std::vector<double>> polyx(
                    degree, std::vector<double>(poly_length));
                  std::vector<std::vector<double>> polyy(
                    degree, std::vector<double>(poly_length));
                  for (unsigned int i = 0; i < degree; ++i)
                    {
                      // Compute all required 1d polynomials and their
                      // derivatives, starting at degree 2. e.g. to access
                      // L'_{3}(2x-1) use polyx[1][1].
                      IntegratedLegendrePolynomials[i + 2].value(
                        cell_points[q][0], polyx[i]);
                      IntegratedLegendrePolynomials[i + 2].value(
                        cell_points[q][1], polyy[i]);
                    }
                  // Now use these to compute the shape functions:
                  if (flags & update_values)
                    {
                      for (unsigned int j = 0; j < degree; ++j)
                        {
                          const unsigned int shift_j(j * degree);
                          for (unsigned int i = 0; i < degree; ++i)
                            {
                              const unsigned int shift_ij(i + shift_j);

                              // Type 1:
                              const unsigned int dof_index1(cell_type1_offset +
                                                            shift_ij);
                              data.shape_values[dof_index1][q][0] =
                                2.0 * polyx[i][1] * polyy[j][0];
                              data.shape_values[dof_index1][q][1] =
                                2.0 * polyx[i][0] * polyy[j][1];

                              // Type 2:
                              const unsigned int dof_index2(cell_type2_offset +
                                                            shift_ij);
                              data.shape_values[dof_index2][q][0] =
                                data.shape_values[dof_index1][q][0];
                              data.shape_values[dof_index2][q][1] =
                                -1.0 * data.shape_values[dof_index1][q][1];
                            }
                          // Type 3:
                          const unsigned int dof_index3_1(cell_type3_offset1 +
                                                          j);
                          data.shape_values[dof_index3_1][q][0] = polyy[j][0];
                          data.shape_values[dof_index3_1][q][1] = 0.0;

                          const unsigned int dof_index3_2(cell_type3_offset2 +
                                                          j);
                          data.shape_values[dof_index3_2][q][0] = 0.0;
                          data.shape_values[dof_index3_2][q][1] = polyx[j][0];
                        }
                    }
                  if (flags & update_gradients)
                    {
                      for (unsigned int j = 0; j < degree; ++j)
                        {
                          const unsigned int shift_j(j * degree);
                          for (unsigned int i = 0; i < degree; ++i)
                            {
                              const unsigned int shift_ij(i + shift_j);

                              // Type 1:
                              const unsigned int dof_index1(cell_type1_offset +
                                                            shift_ij);
                              data.shape_grads[dof_index1][q][0][0] =
                                4.0 * polyx[i][2] * polyy[j][0];
                              data.shape_grads[dof_index1][q][0][1] =
                                4.0 * polyx[i][1] * polyy[j][1];
                              data.shape_grads[dof_index1][q][1][0] =
                                data.shape_grads[dof_index1][q][0][1];
                              data.shape_grads[dof_index1][q][1][1] =
                                4.0 * polyx[i][0] * polyy[j][2];

                              // Type 2:
                              const unsigned int dof_index2(cell_type2_offset +
                                                            shift_ij);
                              data.shape_grads[dof_index2][q][0][0] =
                                data.shape_grads[dof_index1][q][0][0];
                              data.shape_grads[dof_index2][q][0][1] =
                                data.shape_grads[dof_index1][q][0][1];
                              data.shape_grads[dof_index2][q][1][0] =
                                -1.0 * data.shape_grads[dof_index1][q][1][0];
                              data.shape_grads[dof_index2][q][1][1] =
                                -1.0 * data.shape_grads[dof_index1][q][1][1];
                            }
                          // Type 3:
                          const unsigned int dof_index3_1(cell_type3_offset1 +
                                                          j);
                          data.shape_grads[dof_index3_1][q][0][0] = 0.0;
                          data.shape_grads[dof_index3_1][q][0][1] =
                            2.0 * polyy[j][1];
                          data.shape_grads[dof_index3_1][q][1][0] = 0.0;
                          data.shape_grads[dof_index3_1][q][1][1] = 0.0;

                          const unsigned int dof_index3_2(cell_type3_offset2 +
                                                          j);
                          data.shape_grads[dof_index3_2][q][0][0] = 0.0;
                          data.shape_grads[dof_index3_2][q][0][1] = 0.0;
                          data.shape_grads[dof_index3_2][q][1][0] =
                            2.0 * polyx[j][1];
                          data.shape_grads[dof_index3_2][q][1][1] = 0.0;
                        }
                    }
                  if (flags & update_hessians)
                    {
                      for (unsigned int j = 0; j < degree; ++j)
                        {
                          const unsigned int shift_j(j * degree);
                          for (unsigned int i = 0; i < degree; ++i)
                            {
                              const unsigned int shift_ij(i + shift_j);

                              // Type 1:
                              const unsigned int dof_index1(cell_type1_offset +
                                                            shift_ij);
                              data.shape_hessians[dof_index1][q][0][0][0] =
                                8.0 * polyx[i][3] * polyy[j][0];
                              data.shape_hessians[dof_index1][q][1][0][0] =
                                8.0 * polyx[i][2] * polyy[j][1];

                              data.shape_hessians[dof_index1][q][0][1][0] =
                                data.shape_hessians[dof_index1][q][1][0][0];
                              data.shape_hessians[dof_index1][q][1][1][0] =
                                8.0 * polyx[i][1] * polyy[j][2];

                              data.shape_hessians[dof_index1][q][0][0][1] =
                                data.shape_hessians[dof_index1][q][1][0][0];
                              data.shape_hessians[dof_index1][q][1][0][1] =
                                data.shape_hessians[dof_index1][q][1][1][0];

                              data.shape_hessians[dof_index1][q][0][1][1] =
                                data.shape_hessians[dof_index1][q][1][1][0];
                              data.shape_hessians[dof_index1][q][1][1][1] =
                                8.0 * polyx[i][0] * polyy[j][3];



                              // Type 2:
                              const unsigned int dof_index2(cell_type2_offset +
                                                            shift_ij);
                              for (unsigned int d = 0; d < dim; ++d)
                                {
                                  data.shape_hessians[dof_index2][q][0][0][d] =
                                    data.shape_hessians[dof_index1][q][0][0][d];
                                  data.shape_hessians[dof_index2][q][0][1][d] =
                                    data.shape_hessians[dof_index1][q][0][1][d];
                                  data.shape_hessians[dof_index2][q][1][0][d] =
                                    -1.0 *
                                    data.shape_hessians[dof_index1][q][1][0][d];
                                  data.shape_hessians[dof_index2][q][1][1][d] =
                                    -1.0 *
                                    data.shape_hessians[dof_index1][q][1][1][d];
                                }
                            }
                          // Type 3:
                          const unsigned int dof_index3_1(cell_type3_offset1 +
                                                          j);
                          data.shape_hessians[dof_index3_1][q][0][0][0] = 0.0;
                          data.shape_hessians[dof_index3_1][q][0][0][1] = 0.0;
                          data.shape_hessians[dof_index3_1][q][0][1][0] = 0.0;
                          data.shape_hessians[dof_index3_1][q][0][1][1] =
                            4.0 * polyy[j][2];
                          data.shape_hessians[dof_index3_1][q][1][0][0] = 0.0;
                          data.shape_hessians[dof_index3_1][q][1][0][1] = 0.0;
                          data.shape_hessians[dof_index3_1][q][1][1][0] = 0.0;
                          data.shape_hessians[dof_index3_1][q][1][1][1] = 0.0;

                          const unsigned int dof_index3_2(cell_type3_offset2 +
                                                          j);
                          data.shape_hessians[dof_index3_2][q][0][0][0] = 0.0;
                          data.shape_hessians[dof_index3_2][q][0][0][1] = 0.0;
                          data.shape_hessians[dof_index3_2][q][0][1][0] = 0.0;
                          data.shape_hessians[dof_index3_2][q][0][1][1] = 0.0;
                          data.shape_hessians[dof_index3_2][q][1][0][0] =
                            4.0 * polyx[j][2];
                          data.shape_hessians[dof_index3_2][q][1][0][1] = 0.0;
                          data.shape_hessians[dof_index3_2][q][1][1][0] = 0.0;
                          data.shape_hessians[dof_index3_2][q][1][1][1] = 0.0;
                        }
                    }
                }
            }
          break;
        }
      case 3:
        {
          // Compute values of sigma & lambda and the sigma differences and
          // lambda additions.
          std::vector<std::vector<double>> sigma(
            n_q_points, std::vector<double>(lines_per_cell));
          std::vector<std::vector<double>> lambda(
            n_q_points, std::vector<double>(lines_per_cell));
          for (unsigned int q = 0; q < n_q_points; ++q)
            {
              sigma[q][0] = (1.0 - p_list[q][0]) + (1.0 - p_list[q][1]) +
                            (1 - p_list[q][2]);
              sigma[q][1] =
                p_list[q][0] + (1.0 - p_list[q][1]) + (1 - p_list[q][2]);
              sigma[q][2] =
                (1.0 - p_list[q][0]) + p_list[q][1] + (1 - p_list[q][2]);
              sigma[q][3] = p_list[q][0] + p_list[q][1] + (1 - p_list[q][2]);
              sigma[q][4] =
                (1.0 - p_list[q][0]) + (1.0 - p_list[q][1]) + p_list[q][2];
              sigma[q][5] = p_list[q][0] + (1.0 - p_list[q][1]) + p_list[q][2];
              sigma[q][6] = (1.0 - p_list[q][0]) + p_list[q][1] + p_list[q][2];
              sigma[q][7] = p_list[q][0] + p_list[q][1] + p_list[q][2];

              lambda[q][0] = (1.0 - p_list[q][0]) * (1.0 - p_list[q][1]) *
                             (1.0 - p_list[q][2]);
              lambda[q][1] =
                p_list[q][0] * (1.0 - p_list[q][1]) * (1.0 - p_list[q][2]);
              lambda[q][2] =
                (1.0 - p_list[q][0]) * p_list[q][1] * (1.0 - p_list[q][2]);
              lambda[q][3] = p_list[q][0] * p_list[q][1] * (1.0 - p_list[q][2]);
              lambda[q][4] =
                (1.0 - p_list[q][0]) * (1.0 - p_list[q][1]) * p_list[q][2];
              lambda[q][5] = p_list[q][0] * (1.0 - p_list[q][1]) * p_list[q][2];
              lambda[q][6] = (1.0 - p_list[q][0]) * p_list[q][1] * p_list[q][2];
              lambda[q][7] = p_list[q][0] * p_list[q][1] * p_list[q][2];

              // Compute values of sigma_imj = \sigma_{i} - \sigma_{j}
              // and lambda_ipj = \lambda_{i} + \lambda_{j}.
              for (unsigned int i = 0; i < vertices_per_cell; ++i)
                {
                  for (unsigned int j = 0; j < vertices_per_cell; ++j)
                    {
                      data.sigma_imj_values[q][i][j] =
                        sigma[q][i] - sigma[q][j];
                    }
                }
            }

          // We now want some additional information about
          // sigma_imj_values[q][i][j] = sigma[q][i]-sigma[q][j] In order to
          // calculate values & derivatives of the shape functions we need to
          // know:
          // - The component the sigma_imj value corresponds to - this varies
          // with i & j.
          // - The gradient of the sigma_imj value
          //   - this depends on the component and the direction of the
          //   corresponding edge.
          //   - the direction of the edge is determined by
          //   sigma_imj_sign[i][j].
          //
          // Note that not every i,j combination is a valid edge (there are only
          // 12 valid edges in 3d), but we compute them all as it simplifies
          // things.

          // store the sign of each component x, y, z in the sigma list.
          // can use this to fill in the sigma_imj_component data.
          const int sigma_comp_signs[GeometryInfo<3>::vertices_per_cell][3] = {
            {-1, -1, -1},
            {1, -1, -1},
            {-1, 1, -1},
            {1, 1, -1},
            {-1, -1, 1},
            {1, -1, 1},
            {-1, 1, 1},
            {1, 1, 1}};

          int          sigma_imj_sign[vertices_per_cell][vertices_per_cell];
          unsigned int sigma_imj_component[vertices_per_cell]
                                          [vertices_per_cell];

          for (unsigned int i = 0; i < vertices_per_cell; ++i)
            {
              for (unsigned int j = 0; j < vertices_per_cell; ++j)
                {
                  // sigma_imj_sign is the sign (+/-) of the coefficient of
                  // x/y/z in sigma_imj. Due to the numbering of vertices on the
                  // reference element this is easy to work out because edges in
                  // the positive direction go from smaller to higher local
                  // vertex numbering.
                  sigma_imj_sign[i][j] = (i < j) ? -1 : 1;
                  sigma_imj_sign[i][j] = (i == j) ? 0 : sigma_imj_sign[i][j];

                  // Now store the component which the sigma_i - sigma_j
                  // corresponds to:
                  sigma_imj_component[i][j] = 0;
                  for (unsigned int d = 0; d < dim; ++d)
                    {
                      int temp_imj =
                        sigma_comp_signs[i][d] - sigma_comp_signs[j][d];
                      // Only interested in the first non-zero
                      // as if there is a second, it will not be a valid edge.
                      if (temp_imj != 0)
                        {
                          sigma_imj_component[i][j] = d;
                          break;
                        }
                    }
                  // Can now calculate the gradient, only non-zero in the
                  // component given: Note some i,j combinations will be
                  // incorrect, but only on invalid edges.
                  data.sigma_imj_grads[i][j][sigma_imj_component[i][j]] =
                    2.0 * sigma_imj_sign[i][j];
                }
            }
          // Now compute the edge parameterisations for a single element
          // with global numbering matching that of the reference element:

          // resize the edge parameterisations
          data.edge_sigma_values.resize(lines_per_cell);
          data.edge_lambda_values.resize(lines_per_cell);
          data.edge_sigma_grads.resize(lines_per_cell);
          data.edge_lambda_grads_3d.resize(lines_per_cell);
          data.edge_lambda_gradgrads_3d.resize(lines_per_cell);
          for (unsigned int m = 0; m < lines_per_cell; ++m)
            {
              data.edge_sigma_values[m].resize(n_q_points);
              data.edge_lambda_values[m].resize(n_q_points);

              // sigma grads are constant in a cell (no need for quad points)
              data.edge_sigma_grads[m].resize(dim);

              data.edge_lambda_grads_3d[m].resize(n_q_points);
              for (unsigned int q = 0; q < n_q_points; ++q)
                {
                  data.edge_lambda_grads_3d[m][q].resize(dim);
                }
              // lambda_gradgrads are constant in a cell (no need for quad
              // points)
              data.edge_lambda_gradgrads_3d[m].resize(dim);
              for (unsigned int d = 0; d < dim; ++d)
                {
                  data.edge_lambda_gradgrads_3d[m][d].resize(dim);
                }
            }

          // Fill the values:
          const unsigned int
            edge_sigma_direction[GeometryInfo<3>::lines_per_cell] = {
              1, 1, 0, 0, 1, 1, 0, 0, 2, 2, 2, 2};

          for (unsigned int m = 0; m < lines_per_cell; ++m)
            {
              // e1=max(reference vertex numbering on this edge)
              // e2=min(reference vertex numbering on this edge)
              // Which is guaranteed to be:
              const unsigned int e1(
                GeometryInfo<dim>::line_to_cell_vertices(m, 1));
              const unsigned int e2(
                GeometryInfo<dim>::line_to_cell_vertices(m, 0));

              for (unsigned int q = 0; q < n_q_points; ++q)
                {
                  data.edge_sigma_values[m][q] =
                    data.sigma_imj_values[q][e2][e1];
                  data.edge_lambda_values[m][q] = lambda[q][e1] + lambda[q][e2];
                }

              data.edge_sigma_grads[m][edge_sigma_direction[m]] = -2.0;
            }
          // edge_lambda_grads
          for (unsigned int q = 0; q < n_q_points; ++q)
            {
              double x(p_list[q][0]);
              double y(p_list[q][1]);
              double z(p_list[q][2]);
              data.edge_lambda_grads_3d[0][q]  = {z - 1.0, 0.0, x - 1.0};
              data.edge_lambda_grads_3d[1][q]  = {1.0 - z, 0.0, -x};
              data.edge_lambda_grads_3d[2][q]  = {0.0, z - 1.0, y - 1.0};
              data.edge_lambda_grads_3d[3][q]  = {0.0, 1.0 - z, -y};
              data.edge_lambda_grads_3d[4][q]  = {-z, 0.0, 1.0 - x};
              data.edge_lambda_grads_3d[5][q]  = {z, 0.0, x};
              data.edge_lambda_grads_3d[6][q]  = {0.0, -z, 1.0 - y};
              data.edge_lambda_grads_3d[7][q]  = {0.0, z, y};
              data.edge_lambda_grads_3d[8][q]  = {y - 1.0, x - 1.0, 0.0};
              data.edge_lambda_grads_3d[9][q]  = {1.0 - y, -x, 0.0};
              data.edge_lambda_grads_3d[10][q] = {-y, 1.0 - x, 0.0};
              data.edge_lambda_grads_3d[11][q] = {y, x, 0.0};
            }
          // edge_lambda gradgrads:
          const int edge_lambda_sign[GeometryInfo<3>::lines_per_cell] = {
            1,
            -1,
            1,
            -1,
            -1,
            1,
            -1,
            1,
            1,
            -1,
            -1,
            1}; // sign of the 2nd derivative for each edge.
          const unsigned int
            edge_lambda_directions[GeometryInfo<3>::lines_per_cell][2] = {
              {0, 2},
              {0, 2},
              {1, 2},
              {1, 2},
              {0, 2},
              {0, 2},
              {1, 2},
              {1, 2},
              {0, 1},
              {0, 1},
              {0, 1},
              {0, 1}}; // component which edge_lambda[m] depends on.
          for (unsigned int m = 0; m < lines_per_cell; ++m)
            {
              data.edge_lambda_gradgrads_3d[m][edge_lambda_directions[m][0]]
                                           [edge_lambda_directions[m][1]] =
                edge_lambda_sign[m];
              data.edge_lambda_gradgrads_3d[m][edge_lambda_directions[m][1]]
                                           [edge_lambda_directions[m][0]] =
                edge_lambda_sign[m];
            }
          // Precomputation for higher order shape functions,
          // and the face parameterisation.
          if (degree > 0)
            {
              // resize required data:
              data.face_lambda_values.resize(faces_per_cell);
              data.face_lambda_grads.resize(faces_per_cell);
              // for face-based shape functions:
              for (unsigned int m = 0; m < faces_per_cell; ++m)
                {
                  data.face_lambda_values[m].resize(n_q_points);
                  data.face_lambda_grads[m].resize(3);
                }
              // Fill in the values (these don't change between cells).
              for (unsigned int q = 0; q < n_q_points; ++q)
                {
                  double x(p_list[q][0]);
                  double y(p_list[q][1]);
                  double z(p_list[q][2]);
                  data.face_lambda_values[0][q] = 1.0 - x;
                  data.face_lambda_values[1][q] = x;
                  data.face_lambda_values[2][q] = 1.0 - y;
                  data.face_lambda_values[3][q] = y;
                  data.face_lambda_values[4][q] = 1.0 - z;
                  data.face_lambda_values[5][q] = z;
                }
              // gradients are constant:
              data.face_lambda_grads[0] = {-1.0, 0.0, 0.0};
              data.face_lambda_grads[1] = {1.0, 0.0, 0.0};
              data.face_lambda_grads[2] = {0.0, -1.0, 0.0};
              data.face_lambda_grads[3] = {0.0, 1.0, 0.0};
              data.face_lambda_grads[4] = {0.0, 0.0, -1.0};
              data.face_lambda_grads[5] = {0.0, 0.0, 1.0};

              // for cell-based shape functions:
              // these don't depend on the cell, so can precompute all here:
              if (flags & (update_values | update_gradients | update_hessians))
                {
                  // Cell-based shape functions:
                  //
                  // Type-1 (gradients):
                  // \phi^{C_{1}}_{ijk} = grad(
                  // L_{i+2}(2x-1)L_{j+2}(2y-1)L_{k+2}(2z-1) ),
                  //
                  // 0 <= i,j,k < degree. (in a group of degree*degree*degree)
                  const unsigned int cell_type1_offset(n_line_dofs +
                                                       n_face_dofs);
                  // Type-2:
                  //
                  // \phi^{C_{2}}_{ijk} = diag(1, -1, 1)\phi^{C_{1}}_{ijk}
                  // \phi^{C_{2}}_{ijk + p^3} = diag(1, -1,
                  // -1)\phi^{C_{1}}_{ijk}
                  //
                  // 0 <= i,j,k < degree. (subtypes in groups of
                  // degree*degree*degree)
                  //
                  // here we order so that all of subtype 1 comes first, then
                  // subtype 2.
                  const unsigned int cell_type2_offset1(
                    cell_type1_offset + degree * degree * degree);
                  const unsigned int cell_type2_offset2(
                    cell_type2_offset1 + degree * degree * degree);
                  // Type-3
                  // \phi^{C_{3}}_{jk} = L_{j+2}(2y-1)L_{k+2}(2z-1)e_{x}
                  // \phi^{C_{3}}_{ik} = L_{i+2}(2x-1)L_{k+2}(2z-1)e_{y}
                  // \phi^{C_{3}}_{ij} = L_{i+2}(2x-1)L_{j+2}(2y-1)e_{z}
                  //
                  // 0 <= i,j,k < degree. (subtypes in groups of degree*degree)
                  //
                  // again we order so we compute all of subtype 1 first, then
                  // subtype 2, etc.
                  const unsigned int cell_type3_offset1(
                    cell_type2_offset2 + degree * degree * degree);
                  const unsigned int cell_type3_offset2(cell_type3_offset1 +
                                                        degree * degree);
                  const unsigned int cell_type3_offset3(cell_type3_offset2 +
                                                        degree * degree);

                  // compute all points we must evaluate the 1d polynomials at:
                  std::vector<Point<dim>> cell_points(n_q_points);
                  for (unsigned int q = 0; q < n_q_points; ++q)
                    {
                      for (unsigned int d = 0; d < dim; ++d)
                        {
                          cell_points[q][d] = 2.0 * p_list[q][d] - 1.0;
                        }
                    }

                  // We only need poly values and 1st derivative for
                  // update_values, but need the 2nd derivative too for
                  // update_gradients. For update_hessians we also need 3rd
                  // derivative.
                  const unsigned int poly_length =
                    (flags & update_hessians) ?
                      4 :
                      ((flags & update_gradients) ? 3 : 2);

                  // Loop through quad points:
                  for (unsigned int q = 0; q < n_q_points; ++q)
                    {
                      // pre-compute values & required derivatives at this quad
                      // point, (x,y,z): polyx = L_{i+2}(2x-1), polyy =
                      // L_{j+2}(2y-1), polyz = L_{k+2}(2z-1).
                      //
                      // for each polyc[d], c=x,y,z, contains the d-th
                      // derivative with respect to the coordinate c.
                      std::vector<std::vector<double>> polyx(
                        degree, std::vector<double>(poly_length));
                      std::vector<std::vector<double>> polyy(
                        degree, std::vector<double>(poly_length));
                      std::vector<std::vector<double>> polyz(
                        degree, std::vector<double>(poly_length));
                      for (unsigned int i = 0; i < degree; ++i)
                        {
                          // compute all required 1d polynomials for i
                          IntegratedLegendrePolynomials[i + 2].value(
                            cell_points[q][0], polyx[i]);
                          IntegratedLegendrePolynomials[i + 2].value(
                            cell_points[q][1], polyy[i]);
                          IntegratedLegendrePolynomials[i + 2].value(
                            cell_points[q][2], polyz[i]);
                        }
                      // Now use these to compute the shape functions:
                      if (flags & update_values)
                        {
                          for (unsigned int k = 0; k < degree; ++k)
                            {
                              const unsigned int shift_k(k * degree * degree);
                              const unsigned int shift_j(
                                k * degree); // Used below when subbing k for j
                                             // (type 3)
                              for (unsigned int j = 0; j < degree; ++j)
                                {
                                  const unsigned int shift_jk(j * degree +
                                                              shift_k);
                                  for (unsigned int i = 0; i < degree; ++i)
                                    {
                                      const unsigned int shift_ijk(shift_jk +
                                                                   i);

                                      // Type 1:
                                      const unsigned int dof_index1(
                                        cell_type1_offset + shift_ijk);

                                      data.shape_values[dof_index1][q][0] =
                                        2.0 * polyx[i][1] * polyy[j][0] *
                                        polyz[k][0];
                                      data.shape_values[dof_index1][q][1] =
                                        2.0 * polyx[i][0] * polyy[j][1] *
                                        polyz[k][0];
                                      data.shape_values[dof_index1][q][2] =
                                        2.0 * polyx[i][0] * polyy[j][0] *
                                        polyz[k][1];

                                      // Type 2:
                                      const unsigned int dof_index2_1(
                                        cell_type2_offset1 + shift_ijk);
                                      const unsigned int dof_index2_2(
                                        cell_type2_offset2 + shift_ijk);

                                      data.shape_values[dof_index2_1][q][0] =
                                        data.shape_values[dof_index1][q][0];
                                      data.shape_values[dof_index2_1][q][1] =
                                        -1.0 *
                                        data.shape_values[dof_index1][q][1];
                                      data.shape_values[dof_index2_1][q][2] =
                                        data.shape_values[dof_index1][q][2];

                                      data.shape_values[dof_index2_2][q][0] =
                                        data.shape_values[dof_index1][q][0];
                                      data.shape_values[dof_index2_2][q][1] =
                                        -1.0 *
                                        data.shape_values[dof_index1][q][1];
                                      data.shape_values[dof_index2_2][q][2] =
                                        -1.0 *
                                        data.shape_values[dof_index1][q][2];
                                    }
                                  // Type 3: (note we re-use k and j for
                                  // convenience):
                                  const unsigned int shift_ij(
                                    j + shift_j); // here we've subbed j for i,
                                                  // k for j.
                                  const unsigned int dof_index3_1(
                                    cell_type3_offset1 + shift_ij);
                                  const unsigned int dof_index3_2(
                                    cell_type3_offset2 + shift_ij);
                                  const unsigned int dof_index3_3(
                                    cell_type3_offset3 + shift_ij);

                                  data.shape_values[dof_index3_1][q][0] =
                                    polyy[j][0] * polyz[k][0];
                                  data.shape_values[dof_index3_1][q][1] = 0.0;
                                  data.shape_values[dof_index3_1][q][2] = 0.0;

                                  data.shape_values[dof_index3_2][q][0] = 0.0;
                                  data.shape_values[dof_index3_2][q][1] =
                                    polyx[j][0] * polyz[k][0];
                                  data.shape_values[dof_index3_2][q][2] = 0.0;

                                  data.shape_values[dof_index3_3][q][0] = 0.0;
                                  data.shape_values[dof_index3_3][q][1] = 0.0;
                                  data.shape_values[dof_index3_3][q][2] =
                                    polyx[j][0] * polyy[k][0];
                                }
                            }
                        }
                      if (flags & update_gradients)
                        {
                          for (unsigned int k = 0; k < degree; ++k)
                            {
                              const unsigned int shift_k(k * degree * degree);
                              const unsigned int shift_j(
                                k * degree); // Used below when subbing k for j
                                             // (type 3)
                              for (unsigned int j = 0; j < degree; ++j)
                                {
                                  const unsigned int shift_jk(j * degree +
                                                              shift_k);
                                  for (unsigned int i = 0; i < degree; ++i)
                                    {
                                      const unsigned int shift_ijk(shift_jk +
                                                                   i);

                                      // Type 1:
                                      const unsigned int dof_index1(
                                        cell_type1_offset + shift_ijk);

                                      data.shape_grads[dof_index1][q][0][0] =
                                        4.0 * polyx[i][2] * polyy[j][0] *
                                        polyz[k][0];
                                      data.shape_grads[dof_index1][q][0][1] =
                                        4.0 * polyx[i][1] * polyy[j][1] *
                                        polyz[k][0];
                                      data.shape_grads[dof_index1][q][0][2] =
                                        4.0 * polyx[i][1] * polyy[j][0] *
                                        polyz[k][1];

                                      data.shape_grads[dof_index1][q][1][0] =
                                        data.shape_grads[dof_index1][q][0][1];
                                      data.shape_grads[dof_index1][q][1][1] =
                                        4.0 * polyx[i][0] * polyy[j][2] *
                                        polyz[k][0];
                                      data.shape_grads[dof_index1][q][1][2] =
                                        4.0 * polyx[i][0] * polyy[j][1] *
                                        polyz[k][1];

                                      data.shape_grads[dof_index1][q][2][0] =
                                        data.shape_grads[dof_index1][q][0][2];
                                      data.shape_grads[dof_index1][q][2][1] =
                                        data.shape_grads[dof_index1][q][1][2];
                                      data.shape_grads[dof_index1][q][2][2] =
                                        4.0 * polyx[i][0] * polyy[j][0] *
                                        polyz[k][2];

                                      // Type 2:
                                      const unsigned int dof_index2_1(
                                        cell_type2_offset1 + shift_ijk);
                                      const unsigned int dof_index2_2(
                                        cell_type2_offset2 + shift_ijk);

                                      for (unsigned int d = 0; d < dim; ++d)
                                        {
                                          data.shape_grads[dof_index2_1][q][0]
                                                          [d] =
                                            data
                                              .shape_grads[dof_index1][q][0][d];
                                          data.shape_grads[dof_index2_1][q][1]
                                                          [d] =
                                            -1.0 *
                                            data
                                              .shape_grads[dof_index1][q][1][d];
                                          data.shape_grads[dof_index2_1][q][2]
                                                          [d] =
                                            data
                                              .shape_grads[dof_index1][q][2][d];

                                          data.shape_grads[dof_index2_2][q][0]
                                                          [d] =
                                            data
                                              .shape_grads[dof_index1][q][0][d];
                                          data.shape_grads[dof_index2_2][q][1]
                                                          [d] =
                                            -1.0 *
                                            data
                                              .shape_grads[dof_index1][q][1][d];
                                          data.shape_grads[dof_index2_2][q][2]
                                                          [d] =
                                            -1.0 *
                                            data
                                              .shape_grads[dof_index1][q][2][d];
                                        }
                                    }
                                  // Type 3: (note we re-use k and j for
                                  // convenience):
                                  const unsigned int shift_ij(
                                    j + shift_j); // here we've subbed j for i,
                                                  // k for j.
                                  const unsigned int dof_index3_1(
                                    cell_type3_offset1 + shift_ij);
                                  const unsigned int dof_index3_2(
                                    cell_type3_offset2 + shift_ij);
                                  const unsigned int dof_index3_3(
                                    cell_type3_offset3 + shift_ij);
                                  for (unsigned int d1 = 0; d1 < dim; ++d1)
                                    {
                                      for (unsigned int d2 = 0; d2 < dim; ++d2)
                                        {
                                          data.shape_grads[dof_index3_1][q][d1]
                                                          [d2] = 0.0;
                                          data.shape_grads[dof_index3_2][q][d1]
                                                          [d2] = 0.0;
                                          data.shape_grads[dof_index3_3][q][d1]
                                                          [d2] = 0.0;
                                        }
                                    }
                                  data.shape_grads[dof_index3_1][q][0][1] =
                                    2.0 * polyy[j][1] * polyz[k][0];
                                  data.shape_grads[dof_index3_1][q][0][2] =
                                    2.0 * polyy[j][0] * polyz[k][1];

                                  data.shape_grads[dof_index3_2][q][1][0] =
                                    2.0 * polyx[j][1] * polyz[k][0];
                                  data.shape_grads[dof_index3_2][q][1][2] =
                                    2.0 * polyx[j][0] * polyz[k][1];

                                  data.shape_grads[dof_index3_3][q][2][0] =
                                    2.0 * polyx[j][1] * polyy[k][0];
                                  data.shape_grads[dof_index3_3][q][2][1] =
                                    2.0 * polyx[j][0] * polyy[k][1];
                                }
                            }
                        }
                      if (flags & update_hessians)
                        {
                          for (unsigned int k = 0; k < degree; ++k)
                            {
                              const unsigned int shift_k(k * degree * degree);
                              const unsigned int shift_j(
                                k * degree); // Used below when subbing k for j
                                             // type 3

                              for (unsigned int j = 0; j < degree; ++j)
                                {
                                  const unsigned int shift_jk(j * degree +
                                                              shift_k);
                                  for (unsigned int i = 0; i < degree; ++i)
                                    {
                                      const unsigned int shift_ijk(shift_jk +
                                                                   i);

                                      // Type 1:
                                      const unsigned int dof_index1(
                                        cell_type1_offset + shift_ijk);

                                      data.shape_hessians[dof_index1][q][0][0]
                                                         [0] =
                                        8.0 * polyx[i][3] * polyy[j][0] *
                                        polyz[k][0];
                                      data.shape_hessians[dof_index1][q][1][0]
                                                         [0] =
                                        8.0 * polyx[i][2] * polyy[j][1] *
                                        polyz[k][0];
                                      data.shape_hessians[dof_index1][q][2][0]
                                                         [0] =
                                        8.0 * polyx[i][2] * polyy[j][0] *
                                        polyz[k][1];

                                      data.shape_hessians[dof_index1][q][0][1]
                                                         [0] =
                                        data.shape_hessians[dof_index1][q][1][0]
                                                           [0];
                                      data.shape_hessians[dof_index1][q][1][1]
                                                         [0] =
                                        8.0 * polyx[i][1] * polyy[j][2] *
                                        polyz[k][0];
                                      data.shape_hessians[dof_index1][q][2][1]
                                                         [0] =
                                        8.0 * polyx[i][1] * polyy[j][1] *
                                        polyz[k][1];

                                      data.shape_hessians[dof_index1][q][0][2]
                                                         [0] =
                                        data.shape_hessians[dof_index1][q][2][0]
                                                           [0];
                                      data.shape_hessians[dof_index1][q][1][2]
                                                         [0] =
                                        data.shape_hessians[dof_index1][q][2][1]
                                                           [0];
                                      data.shape_hessians[dof_index1][q][2][2]
                                                         [0] =
                                        8.0 * polyx[i][1] * polyy[j][0] *
                                        polyz[k][2];


                                      data.shape_hessians[dof_index1][q][0][0]
                                                         [1] =
                                        data.shape_hessians[dof_index1][q][1][0]
                                                           [0];
                                      data.shape_hessians[dof_index1][q][1][0]
                                                         [1] =
                                        data.shape_hessians[dof_index1][q][1][1]
                                                           [0];
                                      data.shape_hessians[dof_index1][q][2][0]
                                                         [1] =
                                        data.shape_hessians[dof_index1][q][2][1]
                                                           [0];

                                      data.shape_hessians[dof_index1][q][0][1]
                                                         [1] =
                                        data.shape_hessians[dof_index1][q][1][1]
                                                           [0];
                                      data.shape_hessians[dof_index1][q][1][1]
                                                         [1] =
                                        8.0 * polyx[i][0] * polyy[j][3] *
                                        polyz[k][0];
                                      data.shape_hessians[dof_index1][q][2][1]
                                                         [1] =
                                        8.0 * polyx[i][0] * polyy[j][2] *
                                        polyz[k][1];

                                      data.shape_hessians[dof_index1][q][0][2]
                                                         [1] =
                                        data.shape_hessians[dof_index1][q][2][1]
                                                           [0];
                                      data.shape_hessians[dof_index1][q][1][2]
                                                         [1] =
                                        data.shape_hessians[dof_index1][q][2][1]
                                                           [1];
                                      data.shape_hessians[dof_index1][q][2][2]
                                                         [1] =
                                        8.0 * polyx[i][0] * polyy[j][1] *
                                        polyz[k][2];


                                      data.shape_hessians[dof_index1][q][0][0]
                                                         [2] =
                                        data.shape_hessians[dof_index1][q][2][0]
                                                           [0];
                                      data.shape_hessians[dof_index1][q][1][0]
                                                         [2] =
                                        data.shape_hessians[dof_index1][q][2][1]
                                                           [0];
                                      data.shape_hessians[dof_index1][q][2][0]
                                                         [2] =
                                        data.shape_hessians[dof_index1][q][2][2]
                                                           [0];

                                      data.shape_hessians[dof_index1][q][0][1]
                                                         [2] =
                                        data.shape_hessians[dof_index1][q][2][1]
                                                           [0];
                                      data.shape_hessians[dof_index1][q][1][1]
                                                         [2] =
                                        data.shape_hessians[dof_index1][q][2][1]
                                                           [1];
                                      data.shape_hessians[dof_index1][q][2][1]
                                                         [2] =
                                        data.shape_hessians[dof_index1][q][2][2]
                                                           [1];

                                      data.shape_hessians[dof_index1][q][0][2]
                                                         [2] =
                                        data.shape_hessians[dof_index1][q][2][2]
                                                           [0];
                                      data.shape_hessians[dof_index1][q][1][2]
                                                         [2] =
                                        data.shape_hessians[dof_index1][q][2][2]
                                                           [1];
                                      data.shape_hessians[dof_index1][q][2][2]
                                                         [2] =
                                        8.0 * polyx[i][0] * polyy[j][0] *
                                        polyz[k][3];


                                      // Type 2:
                                      const unsigned int dof_index2_1(
                                        cell_type2_offset1 + shift_ijk);
                                      const unsigned int dof_index2_2(
                                        cell_type2_offset2 + shift_ijk);

                                      for (unsigned int d1 = 0; d1 < dim; ++d1)
                                        {
                                          for (unsigned int d2 = 0; d2 < dim;
                                               ++d2)
                                            {
                                              data
                                                .shape_hessians[dof_index2_1][q]
                                                               [0][d1][d2] =
                                                data
                                                  .shape_hessians[dof_index1][q]
                                                                 [0][d1][d2];
                                              data
                                                .shape_hessians[dof_index2_1][q]
                                                               [1][d1][d2] =
                                                -1.0 *
                                                data
                                                  .shape_hessians[dof_index1][q]
                                                                 [1][d1][d2];
                                              data
                                                .shape_hessians[dof_index2_1][q]
                                                               [2][d1][d2] =
                                                data
                                                  .shape_hessians[dof_index1][q]
                                                                 [2][d1][d2];

                                              data
                                                .shape_hessians[dof_index2_2][q]
                                                               [0][d1][d2] =
                                                data
                                                  .shape_hessians[dof_index1][q]
                                                                 [0][d1][d2];
                                              data
                                                .shape_hessians[dof_index2_2][q]
                                                               [1][d1][d2] =
                                                -1.0 *
                                                data
                                                  .shape_hessians[dof_index1][q]
                                                                 [1][d1][d2];
                                              data
                                                .shape_hessians[dof_index2_2][q]
                                                               [2][d1][d2] =
                                                -1.0 *
                                                data
                                                  .shape_hessians[dof_index1][q]
                                                                 [2][d1][d2];
                                            }
                                        }
                                    }
                                  // Type 3: (note we re-use k and j for
                                  // convenience):
                                  const unsigned int shift_ij(
                                    j + shift_j); // here we've subbed j for i,
                                                  // k for j.
                                  const unsigned int dof_index3_1(
                                    cell_type3_offset1 + shift_ij);
                                  const unsigned int dof_index3_2(
                                    cell_type3_offset2 + shift_ij);
                                  const unsigned int dof_index3_3(
                                    cell_type3_offset3 + shift_ij);
                                  for (unsigned int d1 = 0; d1 < dim; ++d1)
                                    {
                                      for (unsigned int d2 = 0; d2 < dim; ++d2)
                                        {
                                          for (unsigned int d3 = 0; d3 < dim;
                                               ++d3)
                                            {
                                              data
                                                .shape_hessians[dof_index3_1][q]
                                                               [d1][d2][d3] =
                                                0.0;
                                              data
                                                .shape_hessians[dof_index3_2][q]
                                                               [d1][d2][d3] =
                                                0.0;
                                              data
                                                .shape_hessians[dof_index3_3][q]
                                                               [d1][d2][d3] =
                                                0.0;
                                            }
                                        }
                                    }
                                  data
                                    .shape_hessians[dof_index3_1][q][0][1][1] =
                                    4.0 * polyy[j][2] * polyz[k][0];
                                  data
                                    .shape_hessians[dof_index3_1][q][0][1][2] =
                                    4.0 * polyy[j][1] * polyz[k][1];

                                  data
                                    .shape_hessians[dof_index3_1][q][0][2][1] =
                                    data
                                      .shape_hessians[dof_index3_1][q][0][1][2];
                                  data
                                    .shape_hessians[dof_index3_1][q][0][2][2] =
                                    4.0 * polyy[j][0] * polyz[k][2];


                                  data
                                    .shape_hessians[dof_index3_2][q][1][0][0] =
                                    4.0 * polyx[j][2] * polyz[k][0];
                                  data
                                    .shape_hessians[dof_index3_2][q][1][0][2] =
                                    4.0 * polyx[j][1] * polyz[k][1];

                                  data
                                    .shape_hessians[dof_index3_2][q][1][2][0] =
                                    data
                                      .shape_hessians[dof_index3_2][q][1][0][2];
                                  data
                                    .shape_hessians[dof_index3_2][q][1][2][2] =
                                    4.0 * polyx[j][0] * polyz[k][2];


                                  data
                                    .shape_hessians[dof_index3_3][q][2][0][0] =
                                    4.0 * polyx[j][2] * polyy[k][0];
                                  data
                                    .shape_hessians[dof_index3_3][q][2][0][1] =
                                    4.0 * polyx[j][1] * polyy[k][1];

                                  data
                                    .shape_hessians[dof_index3_3][q][2][1][0] =
                                    data
                                      .shape_hessians[dof_index3_3][q][2][0][1];
                                  data
                                    .shape_hessians[dof_index3_3][q][2][1][1] =
                                    4.0 * polyx[j][0] * polyy[k][2];
                                }
                            }
                        }
                    }
                }
            }
          break;
        }
      default:
        {
          Assert(false, ExcNotImplemented());
        }
    }
  return data_ptr;
}



template <int dim, int spacedim>
void
FE_NedelecSZ<dim, spacedim>::fill_edge_values(
  const typename Triangulation<dim, dim>::cell_iterator &cell,
  const Quadrature<dim> &                                quadrature,
  const InternalData &                                   fe_data) const
{
  // This function handles the cell-dependent construction of the EDGE-based
  // shape functions.
  //
  // Note it will handle both 2d and 3d, in 2d, the edges are faces, but we
  // handle them here.
  //
  // It will fill in the missing parts of fe_data which were not possible to
  // fill in the get_data routine, with respect to the edge-based shape
  // functions.
  //
  // It should be called by the fill_fe_*_values routines in order to complete
  // the basis set at quadrature points on the current cell for each edge.

  const UpdateFlags  flags(fe_data.update_each);
  const unsigned int n_q_points = quadrature.size();

  Assert(!(flags & update_values) ||
           fe_data.shape_values.size() == this->n_dofs_per_cell(),
         ExcDimensionMismatch(fe_data.shape_values.size(),
                              this->n_dofs_per_cell()));
  Assert(!(flags & update_values) ||
           fe_data.shape_values[0].size() == n_q_points,
         ExcDimensionMismatch(fe_data.shape_values[0].size(), n_q_points));

  // Useful constants:
  const unsigned int degree(
    this->degree -
    1); // Note: constructor takes input degree + 1, so need to knock 1 off.

  // Useful geometry info:
  const unsigned int vertices_per_line(2);
  const unsigned int lines_per_cell(GeometryInfo<dim>::lines_per_cell);

  // Calculate edge orderings:
  std::vector<std::vector<unsigned int>> edge_order(
    lines_per_cell, std::vector<unsigned int>(vertices_per_line));


  switch (dim)
    {
      case 2:
        {
          if (flags & (update_values | update_gradients | update_hessians))
            {
              // Define an edge numbering so that each edge, E_{m} = [e^{m}_{1},
              // e^{m}_{2}] e1 = higher global numbering of the two local
              // vertices e2 = lower global numbering of the two local vertices
              std::vector<int> edge_sign(lines_per_cell);
              for (unsigned int m = 0; m < lines_per_cell; ++m)
                {
                  unsigned int v0_loc =
                    GeometryInfo<dim>::line_to_cell_vertices(m, 0);
                  unsigned int v1_loc =
                    GeometryInfo<dim>::line_to_cell_vertices(m, 1);
                  unsigned int v0_glob = cell->vertex_index(v0_loc);
                  unsigned int v1_glob = cell->vertex_index(v1_loc);

                  if (v0_glob > v1_glob)
                    {
                      // Opposite to global numbering on our reference element
                      edge_sign[m] = -1.0;
                    }
                  else
                    {
                      // Aligns with global numbering on our reference element.
                      edge_sign[m] = 1.0;
                    }
                }

              // Define \sigma_{m} = sigma_{e^{m}_{2}} - sigma_{e^{m}_{1}}
              //        \lambda_{m} = \lambda_{e^{m}_{1}} + \lambda_{e^{m}_{2}}
              //
              // To help things, in fe_data, we have precomputed (sigma_{i} -
              // sigma_{j}) and (lambda_{i} + lambda_{j}) for 0<= i,j <
              // lines_per_cell.
              //
              // There are two types:
              // - lower order (1 per edge, m):
              //   \phi_{m}^{\mathcal{N}_{0}} = 1/2 grad(\sigma_{m})\lambda_{m}
              //
              // - higher order (degree per edge, m):
              //   \phi_{i}^{E_{m}} = grad( L_{i+2}(\sigma_{m}) (\lambda_{m}) ).
              //
              //   NOTE: sigma_{m} and lambda_{m} are either a function of x OR
              //   y
              //         and if sigma is of x, then lambda is of y, and vice
              //         versa. This means that grad(\sigma) requires
              //         multiplication by d(sigma)/dx_{i} for the i^th comp of
              //         grad(sigma) and similarly when taking derivatives of
              //         lambda.
              //
              // First handle the lowest order edges (dofs 0 to 3)
              // 0 and 1 are the edges in the y dir. (sigma is function of y,
              // lambda is function of x). 2 and 3 are the edges in the x dir.
              // (sigma is function of x, lambda is function of y).
              //
              // More more info: see GeometryInfo for picture of the standard
              // element.
              //
              // Fill edge-based points:
              //      std::vector<std::vector< Point<dim> > >
              //      edge_points(lines_per_cell, std::vector<Point<dim>>
              //      (n_q_points));

              std::vector<std::vector<double>> edge_sigma_values(
                fe_data.edge_sigma_values);
              std::vector<std::vector<double>> edge_sigma_grads(
                fe_data.edge_sigma_grads);

              std::vector<std::vector<double>> edge_lambda_values(
                fe_data.edge_lambda_values);
              std::vector<std::vector<double>> edge_lambda_grads(
                fe_data.edge_lambda_grads_2d);

              // Adjust the edge_sigma_* for the current cell:
              for (unsigned int m = 0; m < lines_per_cell; ++m)
                {
                  std::transform(edge_sigma_values[m].begin(),
                                 edge_sigma_values[m].end(),
                                 edge_sigma_values[m].begin(),
                                 [&](const double edge_sigma_value) {
                                   return edge_sign[m] * edge_sigma_value;
                                 });

                  std::transform(edge_sigma_grads[m].begin(),
                                 edge_sigma_grads[m].end(),
                                 edge_sigma_grads[m].begin(),
                                 [&](const double edge_sigma_grad) {
                                   return edge_sign[m] * edge_sigma_grad;
                                 });
                }

              // If we want to generate shape gradients then we need second
              // derivatives of the 1d polynomials, but only first derivatives
              // for the shape values.
              const unsigned int poly_length =
                (flags & update_hessians) ?
                  4 :
                  ((flags & update_gradients) ? 3 : 2);


              for (unsigned int m = 0; m < lines_per_cell; ++m)
                {
                  const unsigned int shift_m(m * this->n_dofs_per_line());
                  for (unsigned int q = 0; q < n_q_points; ++q)
                    {
                      // Only compute 1d polynomials if degree>0.
                      std::vector<std::vector<double>> poly(
                        degree, std::vector<double>(poly_length));
                      for (unsigned int i = 1; i < degree + 1; ++i)
                        {
                          // Compute all required 1d polynomials and their
                          // derivatives, starting at degree 2. e.g. to access
                          // L'_{i+2}(edge_sigma) use polyx[i][1].
                          IntegratedLegendrePolynomials[i + 1].value(
                            edge_sigma_values[m][q], poly[i - 1]);
                        }
                      if (flags & update_values)
                        {
                          // Lowest order edge shape functions:
                          for (unsigned int d = 0; d < dim; ++d)
                            {
                              fe_data.shape_values[shift_m][q][d] =
                                0.5 * edge_sigma_grads[m][d] *
                                edge_lambda_values[m][q];
                            }
                          // Higher order edge shape functions:
                          for (unsigned int i = 1; i < degree + 1; ++i)
                            {
                              const unsigned int poly_index = i - 1;
                              const unsigned int dof_index(i + shift_m);
                              for (unsigned int d = 0; d < dim; ++d)
                                {
                                  fe_data.shape_values[dof_index][q][d] =
                                    edge_sigma_grads[m][d] *
                                      poly[poly_index][1] *
                                      edge_lambda_values[m][q] +
                                    poly[poly_index][0] *
                                      edge_lambda_grads[m][d];
                                }
                            }
                        }
                      if (flags & update_gradients)
                        {
                          // Lowest order edge shape functions:
                          for (unsigned int d1 = 0; d1 < dim; ++d1)
                            {
                              for (unsigned int d2 = 0; d2 < dim; ++d2)
                                {
                                  // Note: gradient is constant for a given
                                  // edge.
                                  fe_data.shape_grads[shift_m][q][d1][d2] =
                                    0.5 * edge_sigma_grads[m][d1] *
                                    edge_lambda_grads[m][d2];
                                }
                            }
                          // Higher order edge shape functions:
                          for (unsigned int i = 1; i < degree + 1; ++i)
                            {
                              const unsigned int poly_index = i - 1;
                              const unsigned int dof_index(i + shift_m);

                              fe_data.shape_grads[dof_index][q][0][0] =
                                edge_sigma_grads[m][0] *
                                edge_sigma_grads[m][0] *
                                edge_lambda_values[m][q] * poly[poly_index][2];

                              fe_data.shape_grads[dof_index][q][0][1] =
                                (edge_sigma_grads[m][0] *
                                   edge_lambda_grads[m][1] +
                                 edge_sigma_grads[m][1] *
                                   edge_lambda_grads[m][0]) *
                                poly[poly_index][1];

                              fe_data.shape_grads[dof_index][q][1][0] =
                                fe_data.shape_grads[dof_index][q][0][1];

                              fe_data.shape_grads[dof_index][q][1][1] =
                                edge_sigma_grads[m][1] *
                                edge_sigma_grads[m][1] *
                                edge_lambda_values[m][q] * poly[poly_index][2];
                            }
                        }
                      if (flags & update_hessians)
                        {
                          // Lowest order edge shape function
                          for (unsigned int d1 = 0; d1 < dim; ++d1)
                            {
                              for (unsigned int d2 = 0; d2 < dim; ++d2)
                                {
                                  for (unsigned int d3 = 0; d3 < dim; ++d3)
                                    {
                                      fe_data.shape_hessians[shift_m][q][d1][d2]
                                                            [d3] = 0;
                                    }
                                }
                            }

                          // Higher order edge shape function
                          for (unsigned int i = 0; i < degree; ++i)
                            {
                              const unsigned int dof_index(i + 1 + shift_m);

                              for (unsigned int d1 = 0; d1 < dim; ++d1)
                                {
                                  for (unsigned int d2 = 0; d2 < dim; ++d2)
                                    {
                                      for (unsigned int d3 = 0; d3 < dim; ++d3)
                                        {
                                          fe_data.shape_hessians[dof_index][q]
                                                                [d1][d2][d3] =
                                            edge_sigma_grads[m][d1] *
                                              edge_sigma_grads[m][d2] *
                                              edge_sigma_grads[m][d3] *
                                              poly[i][3] *
                                              edge_lambda_values[m][q] +
                                            poly[i][2] *
                                              (edge_sigma_grads[m][d1] *
                                                 edge_sigma_grads[m][d2] *
                                                 edge_lambda_grads[m][d3] +
                                               edge_sigma_grads[m][d3] *
                                                 edge_sigma_grads[m][d1] *
                                                 edge_lambda_grads[m][d2] +
                                               edge_sigma_grads[m][d3] *
                                                 edge_sigma_grads[m][d2] *
                                                 edge_lambda_grads[m][d1]);
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
          break;
        }
      case 3:
        {
          if (flags & (update_values | update_gradients | update_hessians))
            {
              // Define an edge numbering so that each edge, E_{m} = [e^{m}_{1},
              // e^{m}_{2}] e1 = higher global numbering of the two local
              // vertices e2 = lower global numbering of the two local vertices
              std::vector<int> edge_sign(lines_per_cell);
              for (unsigned int m = 0; m < lines_per_cell; ++m)
                {
                  unsigned int v0_loc =
                    GeometryInfo<dim>::line_to_cell_vertices(m, 0);
                  unsigned int v1_loc =
                    GeometryInfo<dim>::line_to_cell_vertices(m, 1);
                  unsigned int v0_glob = cell->vertex_index(v0_loc);
                  unsigned int v1_glob = cell->vertex_index(v1_loc);

                  if (v0_glob > v1_glob)
                    {
                      // Opposite to global numbering on our reference element
                      edge_sign[m] = -1.0;
                    }
                  else
                    {
                      // Aligns with global numbering on our reference element.
                      edge_sign[m] = 1.0;
                    }
                }

              // Define \sigma_{m} = sigma_{e^{m}_{1}} - sigma_{e^{m}_{2}}
              //        \lambda_{m} = \lambda_{e^{m}_{1}} + \lambda_{e^{m}_{2}}
              //
              // To help things, in fe_data, we have precomputed (sigma_{i} -
              // sigma_{j}) and (lambda_{i} + lambda_{j}) for 0<= i,j <
              // lines_per_cell.
              //
              // There are two types:
              // - lower order (1 per edge, m):
              //   \phi_{m}^{\mathcal{N}_{0}} = 1/2 grad(\sigma_{m})\lambda_{m}
              //
              // - higher order (degree per edge, m):
              //   \phi_{i}^{E_{m}} = grad( L_{i+2}(\sigma_{m}) (\lambda_{m}) ).
              //
              //   NOTE: In the ref cell, sigma_{m} is a function of x OR y OR Z
              //   and lambda_{m} a function of the remaining co-ords.
              //         for example, if sigma is of x, then lambda is of y AND
              //         z, and so on. This means that grad(\sigma) requires
              //         multiplication by d(sigma)/dx_{i} for the i^th comp of
              //         grad(sigma) and similarly when taking derivatives of
              //         lambda.
              //
              // First handle the lowest order edges (dofs 0 to 11)
              // 0 and 1 are the edges in the y dir at z=0. (sigma is a fn of y,
              // lambda is a fn of x & z). 2 and 3 are the edges in the x dir at
              // z=0. (sigma is a fn of x, lambda is a fn of y & z). 4 and 5 are
              // the edges in the y dir at z=1. (sigma is a fn of y, lambda is a
              // fn of x & z). 6 and 7 are the edges in the x dir at z=1. (sigma
              // is a fn of x, lambda is a fn of y & z). 8 and 9 are the edges
              // in the z dir at y=0. (sigma is a fn of z, lambda is a fn of x &
              // y). 10 and 11 are the edges in the z dir at y=1. (sigma is a fn
              // of z, lambda is a fn of x & y).
              //
              // For more info: see GeometryInfo for picture of the standard
              // element.

              // Copy over required edge-based data:
              std::vector<std::vector<double>> edge_sigma_values(
                fe_data.edge_sigma_values);
              std::vector<std::vector<double>> edge_lambda_values(
                fe_data.edge_lambda_values);
              std::vector<std::vector<double>> edge_sigma_grads(
                fe_data.edge_sigma_grads);
              std::vector<std::vector<std::vector<double>>> edge_lambda_grads(
                fe_data.edge_lambda_grads_3d);
              std::vector<std::vector<std::vector<double>>>
                edge_lambda_gradgrads_3d(fe_data.edge_lambda_gradgrads_3d);

              // Adjust the edge_sigma_* for the current cell:
              for (unsigned int m = 0; m < lines_per_cell; ++m)
                {
                  std::transform(edge_sigma_values[m].begin(),
                                 edge_sigma_values[m].end(),
                                 edge_sigma_values[m].begin(),
                                 [&](const double edge_sigma_value) {
                                   return edge_sign[m] * edge_sigma_value;
                                 });
                  std::transform(edge_sigma_grads[m].begin(),
                                 edge_sigma_grads[m].end(),
                                 edge_sigma_grads[m].begin(),
                                 [&](const double edge_sigma_grad) {
                                   return edge_sign[m] * edge_sigma_grad;
                                 });
                }

              // Now calculate the edge-based shape functions:
              // If we want to generate shape gradients then we need second
              // derivatives of the 1d polynomials, but only first derivatives
              // for the shape values.
              const unsigned int poly_length =
                (flags & update_hessians) ?
                  4 :
                  ((flags & update_gradients) ? 3 : 2);

              std::vector<std::vector<double>> poly(
                degree, std::vector<double>(poly_length));
              for (unsigned int m = 0; m < lines_per_cell; ++m)
                {
                  const unsigned int shift_m(m * this->n_dofs_per_line());
                  for (unsigned int q = 0; q < n_q_points; ++q)
                    {
                      // precompute values of all 1d polynomials required:
                      // for example poly[i][1] = L'_{i+2}(edge_sigma_values)
                      if (degree > 0)
                        {
                          for (unsigned int i = 0; i < degree; ++i)
                            {
                              IntegratedLegendrePolynomials[i + 2].value(
                                edge_sigma_values[m][q], poly[i]);
                            }
                        }
                      if (flags & update_values)
                        {
                          // Lowest order edge shape functions:
                          for (unsigned int d = 0; d < dim; ++d)
                            {
                              fe_data.shape_values[shift_m][q][d] =
                                0.5 * edge_sigma_grads[m][d] *
                                edge_lambda_values[m][q];
                            }
                          // Higher order edge shape functions
                          if (degree > 0)
                            {
                              for (unsigned int i = 0; i < degree; ++i)
                                {
                                  const unsigned int dof_index(i + 1 + shift_m);
                                  for (unsigned int d = 0; d < dim; ++d)
                                    {
                                      fe_data.shape_values[dof_index][q][d] =
                                        edge_sigma_grads[m][d] * poly[i][1] *
                                          edge_lambda_values[m][q] +
                                        poly[i][0] * edge_lambda_grads[m][q][d];
                                    }
                                }
                            }
                        }
                      if (flags & update_gradients)
                        {
                          // Lowest order edge shape functions:
                          for (unsigned int d1 = 0; d1 < dim; ++d1)
                            {
                              for (unsigned int d2 = 0; d2 < dim; ++d2)
                                {
                                  fe_data.shape_grads[shift_m][q][d1][d2] =
                                    0.5 * edge_sigma_grads[m][d1] *
                                    edge_lambda_grads[m][q][d2];
                                }
                            }
                          // Higher order edge shape functions
                          if (degree > 0)
                            {
                              for (unsigned int i = 0; i < degree; ++i)
                                {
                                  const unsigned int dof_index(i + 1 + shift_m);

                                  for (unsigned int d1 = 0; d1 < dim; ++d1)
                                    {
                                      for (unsigned int d2 = 0; d2 < dim; ++d2)
                                        {
                                          fe_data
                                            .shape_grads[dof_index][q][d1][d2] =
                                            edge_sigma_grads[m][d1] *
                                              edge_sigma_grads[m][d2] *
                                              poly[i][2] *
                                              edge_lambda_values[m][q] +
                                            (edge_sigma_grads[m][d1] *
                                               edge_lambda_grads[m][q][d2] +
                                             edge_sigma_grads[m][d2] *
                                               edge_lambda_grads[m][q][d1]) *
                                              poly[i][1] +
                                            edge_lambda_gradgrads_3d[m][d1]
                                                                    [d2] *
                                              poly[i][0];
                                        }
                                    }
                                }
                            }
                        }
                      if (flags & update_hessians)
                        {
                          // Lowest order edge shape functions:
                          for (unsigned int d1 = 0; d1 < dim; ++d1)
                            {
                              for (unsigned int d2 = 0; d2 < dim; ++d2)
                                {
                                  for (unsigned int d3 = 0; d3 < dim; ++d3)
                                    {
                                      fe_data.shape_hessians[shift_m][q][d1][d2]
                                                            [d3] =
                                        0.5 * edge_sigma_grads[m][d1] *
                                        edge_lambda_gradgrads_3d[m][d3][d2];
                                    }
                                }
                            }

                          // Higher order edge shape functions
                          if (degree > 0)
                            {
                              for (unsigned int i = 0; i < degree; ++i)
                                {
                                  const unsigned int dof_index(i + 1 + shift_m);

                                  for (unsigned int d1 = 0; d1 < dim; ++d1)
                                    {
                                      for (unsigned int d2 = 0; d2 < dim; ++d2)
                                        {
                                          for (unsigned int d3 = 0; d3 < dim;
                                               ++d3)
                                            {
                                              fe_data
                                                .shape_hessians[dof_index][q]
                                                               [d1][d2][d3] =
                                                edge_sigma_grads[m][d1] *
                                                  edge_sigma_grads[m][d2] *
                                                  edge_sigma_grads[m][d3] *
                                                  poly[i][3] *
                                                  edge_lambda_values[m][q] +
                                                poly[i][2] *
                                                  (edge_sigma_grads[m][d1] *
                                                     edge_sigma_grads[m][d2] *
                                                     edge_lambda_grads[m][q]
                                                                      [d3] +
                                                   edge_sigma_grads[m][d3] *
                                                     edge_sigma_grads[m][d1] *
                                                     edge_lambda_grads[m][q]
                                                                      [d2] +
                                                   edge_sigma_grads[m][d3] *
                                                     edge_sigma_grads[m][d2] *
                                                     edge_lambda_grads[m][q]
                                                                      [d1]) +
                                                poly[i][1] *
                                                  (edge_sigma_grads[m][d1] *
                                                     edge_lambda_gradgrads_3d
                                                       [m][d3][d2] +
                                                   edge_sigma_grads[m][d2] *
                                                     edge_lambda_gradgrads_3d
                                                       [m][d3][d1] +
                                                   edge_sigma_grads[m][d3] *
                                                     edge_lambda_gradgrads_3d
                                                       [m][d2][d1]);
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
          break;
        }
      default:
        {
          Assert(false, ExcNotImplemented());
        }
    }
}



template <int dim, int spacedim>
void
FE_NedelecSZ<dim, spacedim>::fill_face_values(
  const typename Triangulation<dim, dim>::cell_iterator &cell,
  const Quadrature<dim> &                                quadrature,
  const InternalData &                                   fe_data) const
{
  // This function handles the cell-dependent construction of the FACE-based
  // shape functions.
  //
  // Note that it should only be called in 3d.
  Assert(dim == 3, ExcDimensionMismatch(dim, 3));
  //
  // It will fill in the missing parts of fe_data which were not possible to
  // fill in the get_data routine, with respect to face-based shape functions.
  //
  // It should be called by the fill_fe_*_values routines in order to complete
  // the basis set at quadrature points on the current cell for each face.

  // Useful constants:
  const unsigned int degree(
    this->degree -
    1); // Note: constructor takes input degree + 1, so need to knock 1 off.

  // Do nothing if FE degree is 0.
  if (degree > 0)
    {
      const UpdateFlags flags(fe_data.update_each);

      if (flags & (update_values | update_gradients | update_hessians))
        {
          const unsigned int n_q_points = quadrature.size();

          Assert(!(flags & update_values) ||
                   fe_data.shape_values.size() == this->n_dofs_per_cell(),
                 ExcDimensionMismatch(fe_data.shape_values.size(),
                                      this->n_dofs_per_cell()));
          Assert(!(flags & update_values) ||
                   fe_data.shape_values[0].size() == n_q_points,
                 ExcDimensionMismatch(fe_data.shape_values[0].size(),
                                      n_q_points));

          // Useful geometry info:
          const unsigned int vertices_per_face(
            GeometryInfo<dim>::vertices_per_face);
          const unsigned int faces_per_cell(GeometryInfo<dim>::faces_per_cell);

          // DoF info:
          const unsigned int n_line_dofs =
            this->n_dofs_per_line() * GeometryInfo<dim>::lines_per_cell;

          // First we find the global face orientations on the current cell.
          std::vector<std::vector<unsigned int>> face_orientation(
            faces_per_cell, std::vector<unsigned int>(vertices_per_face));

          const unsigned int
            vertex_opposite_on_face[GeometryInfo<3>::vertices_per_face] = {3,
                                                                           2,
                                                                           1,
                                                                           0};

          const unsigned int
            vertices_adjacent_on_face[GeometryInfo<3>::vertices_per_face][2] = {
              {1, 2}, {0, 3}, {0, 3}, {1, 2}};

          for (unsigned int m = 0; m < faces_per_cell; ++m)
            {
              // Find the local vertex on this face with the highest global
              // numbering. This is f^m_0.
              unsigned int current_max  = 0;
              unsigned int current_glob = cell->vertex_index(
                GeometryInfo<dim>::face_to_cell_vertices(m, 0));
              for (unsigned int v = 1; v < vertices_per_face; ++v)
                {
                  if (current_glob <
                      cell->vertex_index(
                        GeometryInfo<dim>::face_to_cell_vertices(m, v)))
                    {
                      current_max  = v;
                      current_glob = cell->vertex_index(
                        GeometryInfo<dim>::face_to_cell_vertices(m, v));
                    }
                }
              face_orientation[m][0] =
                GeometryInfo<dim>::face_to_cell_vertices(m, current_max);

              // f^m_2 is the vertex opposite f^m_0.
              face_orientation[m][2] = GeometryInfo<dim>::face_to_cell_vertices(
                m, vertex_opposite_on_face[current_max]);

              // Finally, f^m_1 is the vertex with the greater global numbering
              // of the remaining two local vertices. Then, f^m_3 is the other.
              if (cell->vertex_index(GeometryInfo<dim>::face_to_cell_vertices(
                    m, vertices_adjacent_on_face[current_max][0])) >
                  cell->vertex_index(GeometryInfo<dim>::face_to_cell_vertices(
                    m, vertices_adjacent_on_face[current_max][1])))
                {
                  face_orientation[m][1] =
                    GeometryInfo<dim>::face_to_cell_vertices(
                      m, vertices_adjacent_on_face[current_max][0]);
                  face_orientation[m][3] =
                    GeometryInfo<dim>::face_to_cell_vertices(
                      m, vertices_adjacent_on_face[current_max][1]);
                }
              else
                {
                  face_orientation[m][1] =
                    GeometryInfo<dim>::face_to_cell_vertices(
                      m, vertices_adjacent_on_face[current_max][1]);
                  face_orientation[m][3] =
                    GeometryInfo<dim>::face_to_cell_vertices(
                      m, vertices_adjacent_on_face[current_max][0]);
                }
            }

          // Now we know the face orientation on the current cell, we can
          // generate the parameterisation:
          std::vector<std::vector<double>> face_xi_values(
            faces_per_cell, std::vector<double>(n_q_points));
          std::vector<std::vector<double>> face_xi_grads(
            faces_per_cell, std::vector<double>(dim));
          std::vector<std::vector<double>> face_eta_values(
            faces_per_cell, std::vector<double>(n_q_points));
          std::vector<std::vector<double>> face_eta_grads(
            faces_per_cell, std::vector<double>(dim));

          std::vector<std::vector<double>> face_lambda_values(
            faces_per_cell, std::vector<double>(n_q_points));
          std::vector<std::vector<double>> face_lambda_grads(
            faces_per_cell, std::vector<double>(dim));
          for (unsigned int m = 0; m < faces_per_cell; ++m)
            {
              for (unsigned int q = 0; q < n_q_points; ++q)
                {
                  face_xi_values[m][q] =
                    fe_data.sigma_imj_values[q][face_orientation[m][0]]
                                            [face_orientation[m][1]];
                  face_eta_values[m][q] =
                    fe_data.sigma_imj_values[q][face_orientation[m][0]]
                                            [face_orientation[m][3]];
                  face_lambda_values[m][q] = fe_data.face_lambda_values[m][q];
                }
              for (unsigned int d = 0; d < dim; ++d)
                {
                  face_xi_grads[m][d] =
                    fe_data.sigma_imj_grads[face_orientation[m][0]]
                                           [face_orientation[m][1]][d];
                  face_eta_grads[m][d] =
                    fe_data.sigma_imj_grads[face_orientation[m][0]]
                                           [face_orientation[m][3]][d];

                  face_lambda_grads[m][d] = fe_data.face_lambda_grads[m][d];
                }
            }
          // Now can generate the basis
          const unsigned int poly_length =
            (flags & update_hessians) ? 4 :
                                        ((flags & update_gradients) ? 3 : 2);


          std::vector<std::vector<double>> polyxi(
            degree, std::vector<double>(poly_length));
          std::vector<std::vector<double>> polyeta(
            degree, std::vector<double>(poly_length));

          // Loop through quad points:
          for (unsigned int m = 0; m < faces_per_cell; ++m)
            {
              // we assume that all quads have the same number of dofs
              const unsigned int shift_m(m * this->n_dofs_per_quad(0));
              // Calculate the offsets for each face-based shape function:
              //
              // Type-1 (gradients)
              // \phi^{F_m,1}_{ij} = \nabla( L_{i+2}(\xi_{F_{m}})
              // L_{j+2}(\eta_{F_{m}}) \lambda_{F_{m}} )
              //
              // 0 <= i,j < degree (in a group of degree*degree)
              const unsigned int face_type1_offset(n_line_dofs + shift_m);
              // Type-2:
              //
              // \phi^{F_m,2}_{ij} = ( L'_{i+2}(\xi_{F_{m}})
              // L_{j+2}(\eta_{F_{m}}) \nabla\xi_{F_{m}}
              //                       - L_{i+2}(\xi_{F_{m}})
              //                       L'_{j+2}(\eta_{F_{m}}) \nabla\eta_{F_{m}}
              //                       ) \lambda_{F_{m}}
              //
              // 0 <= i,j < degree (in a group of degree*degree)
              const unsigned int face_type2_offset(face_type1_offset +
                                                   degree * degree);
              // Type-3:
              //
              // \phi^{F_m,3}_{i} = L_{i+2}(\eta_{F_{m}}) \lambda_{F_{m}}
              // \nabla\xi_{F_{m}} \phi^{F_m,3}_{i+p} = L_{i+2}(\xi_{F_{m}})
              // \lambda_{F_{m}} \nabla\eta_{F_{m}}
              //
              // 0 <= i < degree.
              //
              // here we order so that all of subtype 1 comes first, then
              // subtype 2.
              const unsigned int face_type3_offset1(face_type2_offset +
                                                    degree * degree);
              const unsigned int face_type3_offset2(face_type3_offset1 +
                                                    degree);

              // Loop over all faces:
              for (unsigned int q = 0; q < n_q_points; ++q)
                {
                  // pre-compute values & required derivatives at this quad
                  // point: polyxi = L_{i+2}(\xi_{F_{m}}), polyeta =
                  // L_{j+2}(\eta_{F_{m}}),
                  //
                  // each polypoint[k][d], contains the dth derivative of
                  // L_{k+2} at the point \xi or \eta. Note that this doesn't
                  // include the derivative of xi/eta via the chain rule.
                  for (unsigned int i = 0; i < degree; ++i)
                    {
                      // compute all required 1d polynomials:
                      IntegratedLegendrePolynomials[i + 2].value(
                        face_xi_values[m][q], polyxi[i]);
                      IntegratedLegendrePolynomials[i + 2].value(
                        face_eta_values[m][q], polyeta[i]);
                    }
                  // Now use these to compute the shape functions:
                  if (flags & update_values)
                    {
                      for (unsigned int j = 0; j < degree; ++j)
                        {
                          const unsigned int shift_j(j * degree);
                          for (unsigned int i = 0; i < degree; ++i)
                            {
                              const unsigned int shift_ij(shift_j + i);
                              // Type 1:
                              const unsigned int dof_index1(face_type1_offset +
                                                            shift_ij);
                              for (unsigned int d = 0; d < dim; ++d)
                                {
                                  fe_data.shape_values[dof_index1][q][d] =
                                    (face_xi_grads[m][d] * polyxi[i][1] *
                                       polyeta[j][0] +
                                     face_eta_grads[m][d] * polyxi[i][0] *
                                       polyeta[j][1]) *
                                      face_lambda_values[m][q] +
                                    face_lambda_grads[m][d] * polyxi[i][0] *
                                      polyeta[j][0];
                                }
                              // Type 2:
                              const unsigned int dof_index2(face_type2_offset +
                                                            shift_ij);
                              for (unsigned int d = 0; d < dim; ++d)
                                {
                                  fe_data.shape_values[dof_index2][q][d] =
                                    (face_xi_grads[m][d] * polyxi[i][1] *
                                       polyeta[j][0] -
                                     face_eta_grads[m][d] * polyxi[i][0] *
                                       polyeta[j][1]) *
                                    face_lambda_values[m][q];
                                }
                            }
                          // Type 3:
                          const unsigned int dof_index3_1(face_type3_offset1 +
                                                          j);
                          const unsigned int dof_index3_2(face_type3_offset2 +
                                                          j);
                          for (unsigned int d = 0; d < dim; ++d)
                            {
                              fe_data.shape_values[dof_index3_1][q][d] =
                                face_xi_grads[m][d] * polyeta[j][0] *
                                face_lambda_values[m][q];

                              fe_data.shape_values[dof_index3_2][q][d] =
                                face_eta_grads[m][d] * polyxi[j][0] *
                                face_lambda_values[m][q];
                            }
                        }
                    }
                  if (flags & update_gradients)
                    {
                      for (unsigned int j = 0; j < degree; ++j)
                        {
                          const unsigned int shift_j(j * degree);
                          for (unsigned int i = 0; i < degree; ++i)
                            {
                              const unsigned int shift_ij(shift_j + i);
                              // Type 1:
                              const unsigned int dof_index1(face_type1_offset +
                                                            shift_ij);
                              for (unsigned int d1 = 0; d1 < dim; ++d1)
                                {
                                  for (unsigned int d2 = 0; d2 < dim; ++d2)
                                    {
                                      fe_data
                                        .shape_grads[dof_index1][q][d1][d2] =
                                        (face_xi_grads[m][d1] *
                                           face_xi_grads[m][d2] * polyxi[i][2] *
                                           polyeta[j][0] +
                                         (face_xi_grads[m][d1] *
                                            face_eta_grads[m][d2] +
                                          face_xi_grads[m][d2] *
                                            face_eta_grads[m][d1]) *
                                           polyxi[i][1] * polyeta[j][1] +
                                         face_eta_grads[m][d1] *
                                           face_eta_grads[m][d2] *
                                           polyxi[i][0] * polyeta[j][2]) *
                                          face_lambda_values[m][q] +
                                        (face_xi_grads[m][d2] * polyxi[i][1] *
                                           polyeta[j][0] +
                                         face_eta_grads[m][d2] * polyxi[i][0] *
                                           polyeta[j][1]) *
                                          face_lambda_grads[m][d1] +
                                        (face_xi_grads[m][d1] * polyxi[i][1] *
                                           polyeta[j][0] +
                                         face_eta_grads[m][d1] * polyxi[i][0] *
                                           polyeta[j][1]) *
                                          face_lambda_grads[m][d2];
                                    }
                                }
                              // Type 2:
                              const unsigned int dof_index2(face_type2_offset +
                                                            shift_ij);
                              for (unsigned int d1 = 0; d1 < dim; ++d1)
                                {
                                  for (unsigned int d2 = 0; d2 < dim; ++d2)
                                    {
                                      fe_data
                                        .shape_grads[dof_index2][q][d1][d2] =
                                        (face_xi_grads[m][d1] *
                                           face_xi_grads[m][d2] * polyxi[i][2] *
                                           polyeta[j][0] +
                                         (face_xi_grads[m][d1] *
                                            face_eta_grads[m][d2] -
                                          face_xi_grads[m][d2] *
                                            face_eta_grads[m][d1]) *
                                           polyxi[i][1] * polyeta[j][1] -
                                         face_eta_grads[m][d1] *
                                           face_eta_grads[m][d2] *
                                           polyxi[i][0] * polyeta[j][2]) *
                                          face_lambda_values[m][q] +
                                        (face_xi_grads[m][d1] * polyxi[i][1] *
                                           polyeta[j][0] -
                                         face_eta_grads[m][d1] * polyxi[i][0] *
                                           polyeta[j][1]) *
                                          face_lambda_grads[m][d2];
                                    }
                                }
                            }
                          // Type 3:
                          const unsigned int dof_index3_1(face_type3_offset1 +
                                                          j);
                          const unsigned int dof_index3_2(face_type3_offset2 +
                                                          j);
                          for (unsigned int d1 = 0; d1 < dim; ++d1)
                            {
                              for (unsigned int d2 = 0; d2 < dim; ++d2)
                                {
                                  fe_data.shape_grads[dof_index3_1][q][d1][d2] =
                                    face_xi_grads[m][d1] *
                                    (face_eta_grads[m][d2] * polyeta[j][1] *
                                       face_lambda_values[m][q] +
                                     face_lambda_grads[m][d2] * polyeta[j][0]);

                                  fe_data.shape_grads[dof_index3_2][q][d1][d2] =
                                    face_eta_grads[m][d1] *
                                    (face_xi_grads[m][d2] * polyxi[j][1] *
                                       face_lambda_values[m][q] +
                                     face_lambda_grads[m][d2] * polyxi[j][0]);
                                }
                            }
                        }
                    }
                  if (flags & update_hessians)
                    {
                      for (unsigned int j = 0; j < degree; ++j)
                        {
                          const unsigned int shift_j(j * degree);
                          for (unsigned int i = 0; i < degree; ++i)
                            {
                              const unsigned int shift_ij(shift_j + i);

                              // Type 1:
                              const unsigned int dof_index1(face_type1_offset +
                                                            shift_ij);
                              for (unsigned int d1 = 0; d1 < dim; ++d1)
                                {
                                  for (unsigned int d2 = 0; d2 < dim; ++d2)
                                    {
                                      for (unsigned int d3 = 0; d3 < dim; ++d3)
                                        {
                                          fe_data.shape_hessians[dof_index1][q]
                                                                [d1][d2][d3] =
                                            polyxi[i][1] *
                                              face_xi_grads[m][d3] *
                                              (face_eta_grads[m][d1] *
                                                 (polyeta[j][2] *
                                                    face_eta_grads[m][d2] *
                                                    face_lambda_values[m][q] +
                                                  polyeta[j][1] *
                                                    face_lambda_grads[m][d2]) +
                                               polyeta[j][1] *
                                                 face_eta_grads[m][d2] *
                                                 face_lambda_grads[m][d1]) +
                                            polyxi[i][0] *
                                              (polyeta[j][3] *
                                                 face_eta_grads[m][d1] *
                                                 face_eta_grads[m][d2] *
                                                 face_eta_grads[m][d3] *
                                                 face_lambda_values[m][q] +
                                               polyeta[j][2] *
                                                 (face_eta_grads[m][d1] *
                                                    face_eta_grads[m][d2] *
                                                    face_lambda_grads[m][d3] +
                                                  face_eta_grads[m][d3] *
                                                    (face_eta_grads[m][d1] *
                                                       face_lambda_grads[m]
                                                                        [d2] +
                                                     face_eta_grads[m][d2] *
                                                       face_lambda_grads
                                                         [m][d1]))) +
                                            (polyxi[i][1] * polyeta[j][1] *
                                               face_eta_grads[m][d3] +
                                             polyxi[i][2] * polyeta[j][0] *
                                               face_xi_grads[m][d3]) *
                                              (face_xi_grads[m][d1] *
                                                 face_lambda_grads[m][d2] +
                                               face_xi_grads[m][d2] *
                                                 face_lambda_grads[m][d1]) +
                                            face_lambda_grads[m][d3] *
                                              (polyxi[i][2] * polyeta[j][0] *
                                                 face_xi_grads[m][d1] *
                                                 face_xi_grads[m][d2] +
                                               polyxi[i][1] * polyeta[j][1] *
                                                 (face_xi_grads[m][d1] *
                                                    face_eta_grads[m][d2] +
                                                  face_xi_grads[m][d2] *
                                                    face_eta_grads[m][d1])) +
                                            face_lambda_values[m][q] *
                                              (polyxi[i][3] * polyeta[j][0] *
                                                 face_xi_grads[m][d1] *
                                                 face_xi_grads[m][d2] *
                                                 face_xi_grads[m][d3] +
                                               polyxi[i][1] * polyeta[j][2] *
                                                 face_eta_grads[m][d3] *
                                                 (face_xi_grads[m][d1] *
                                                    face_eta_grads[m][d2] +
                                                  face_xi_grads[m][d2] *
                                                    face_eta_grads[m][d1]) +
                                               polyxi[i][2] * polyeta[j][1] *
                                                 (face_xi_grads[m][d3] *
                                                    face_xi_grads[m][d2] *
                                                    face_eta_grads[m][d1] +
                                                  face_xi_grads[m][d1] *
                                                    (face_xi_grads[m][d2] *
                                                       face_eta_grads[m][d3] +
                                                     face_xi_grads[m][d3] *
                                                       face_eta_grads[m][d2])));
                                        }
                                    }
                                }

                              // Type 2:
                              const unsigned int dof_index2(face_type2_offset +
                                                            shift_ij);
                              for (unsigned int d1 = 0; d1 < dim; ++d1)
                                {
                                  for (unsigned int d2 = 0; d2 < dim; ++d2)
                                    {
                                      for (unsigned int d3 = 0; d3 < dim; ++d3)
                                        {
                                          fe_data.shape_hessians[dof_index2][q]
                                                                [d1][d2][d3] =
                                            face_xi_grads[m][d1] *
                                              (polyxi[i][1] * polyeta[j][1] *
                                                 (face_eta_grads[m][d2] *
                                                    face_lambda_grads[m][d3] +
                                                  face_eta_grads[m][d3] *
                                                    face_lambda_grads[m][d2]) +
                                               polyxi[i][2] * polyeta[j][0] *
                                                 (face_xi_grads[m][d2] *
                                                    face_lambda_grads[m][d3] +
                                                  face_xi_grads[m][d3] *
                                                    face_lambda_grads[m][d2]) +
                                               face_lambda_values[m][q] *
                                                 (face_eta_grads[m][d2] *
                                                    (polyxi[i][1] *
                                                       polyeta[j][2] *
                                                       face_eta_grads[m][d3] +
                                                     polyxi[i][2] *
                                                       polyeta[j][1] *
                                                       face_xi_grads[m][d3]) +
                                                  face_xi_grads[m][d2] *
                                                    (polyxi[i][2] *
                                                       polyeta[j][1] *
                                                       face_eta_grads[m][d3] +
                                                     polyxi[i][3] *
                                                       polyeta[j][0] *
                                                       face_xi_grads[m][d3]))) -
                                            polyxi[i][0] *
                                              face_eta_grads[m][d1] *
                                              (face_eta_grads[m][d2] *
                                                 (polyeta[j][3] *
                                                    face_eta_grads[m][d3] *
                                                    face_lambda_values[m][q] +
                                                  polyeta[j][2] *
                                                    face_lambda_grads[m][d3]) +
                                               polyeta[j][2] *
                                                 face_eta_grads[m][d3] *
                                                 face_lambda_grads[m][d2]) -
                                            face_eta_grads[m][d1] *
                                              (polyxi[i][1] *
                                                 face_xi_grads[m][d3] *
                                                 (polyeta[j][2] *
                                                    face_eta_grads[m][d2] *
                                                    face_lambda_values[m][q] +
                                                  polyeta[j][1] *
                                                    face_lambda_grads[m][d2]) +
                                               face_xi_grads[m][d2] *
                                                 (polyxi[i][1] *
                                                    (polyeta[j][2] *
                                                       face_eta_grads[m][d3] *
                                                       face_lambda_values[m]
                                                                         [q] +
                                                     polyeta[j][1] *
                                                       face_lambda_grads[m]
                                                                        [d3]) +
                                                  polyxi[i][2] * polyeta[j][1] *
                                                    face_xi_grads[m][d3] *
                                                    face_lambda_values[m][q]));
                                        }
                                    }
                                }
                            }
                          // Type 3:
                          const unsigned int dof_index3_1(face_type3_offset1 +
                                                          j);
                          const unsigned int dof_index3_2(face_type3_offset2 +
                                                          j);
                          for (unsigned int d1 = 0; d1 < dim; ++d1)
                            {
                              for (unsigned int d2 = 0; d2 < dim; ++d2)
                                {
                                  for (unsigned int d3 = 0; d3 < dim; ++d3)
                                    {
                                      fe_data.shape_hessians[dof_index3_1][q]
                                                            [d1][d2][d3] =
                                        face_xi_grads[m][d1] *
                                        (face_eta_grads[m][d2] *
                                           (polyeta[j][2] *
                                              face_eta_grads[m][d3] *
                                              face_lambda_values[m][q] +
                                            polyeta[j][1] *
                                              face_lambda_grads[m][d3]) +
                                         face_lambda_grads[m][d2] *
                                           polyeta[j][1] *
                                           face_eta_grads[m][d3]);

                                      fe_data.shape_hessians[dof_index3_2][q]
                                                            [d1][d2][d3] =
                                        face_eta_grads[m][d1] *
                                        (face_xi_grads[m][d2] *
                                           (polyxi[j][2] *
                                              face_xi_grads[m][d3] *
                                              face_lambda_values[m][q] +
                                            polyxi[j][1] *
                                              face_lambda_grads[m][d3]) +
                                         face_lambda_grads[m][d2] *
                                           polyxi[j][1] * face_xi_grads[m][d3]);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}



template <int dim, int spacedim>
void
FE_NedelecSZ<dim, spacedim>::fill_fe_values(
  const typename Triangulation<dim, dim>::cell_iterator &cell,
  const CellSimilarity::Similarity /*cell_similarity*/,
  const Quadrature<dim> &                             quadrature,
  const Mapping<dim, dim> &                           mapping,
  const typename Mapping<dim, dim>::InternalDataBase &mapping_internal,
  const internal::FEValuesImplementation::MappingRelatedData<dim, dim>
    &                                                       mapping_data,
  const typename FiniteElement<dim, dim>::InternalDataBase &fe_internal,
  dealii::internal::FEValuesImplementation::FiniteElementRelatedData<dim, dim>
    &data) const
{
  // Convert to the correct internal data class for this FE class.
  Assert(dynamic_cast<const InternalData *>(&fe_internal) != nullptr,
         ExcInternalError());
  const InternalData &fe_data = static_cast<const InternalData &>(fe_internal);

  // Now update the edge-based DoFs, which depend on the cell.
  // This will fill in the missing items in the InternalData
  // (fe_internal/fe_data) which was not filled in by get_data.
  fill_edge_values(cell, quadrature, fe_data);
  if (dim == 3 && this->degree > 1)
    {
      fill_face_values(cell, quadrature, fe_data);
    }

  const UpdateFlags  flags(fe_data.update_each);
  const unsigned int n_q_points = quadrature.size();

  Assert(!(flags & update_values) ||
           fe_data.shape_values.size() == this->n_dofs_per_cell(),
         ExcDimensionMismatch(fe_data.shape_values.size(),
                              this->n_dofs_per_cell()));
  Assert(!(flags & update_values) ||
           fe_data.shape_values[0].size() == n_q_points,
         ExcDimensionMismatch(fe_data.shape_values[0].size(), n_q_points));

  if (flags & update_values)
    {
      // Now have all shape_values stored on the reference cell.
      // Must now transform to the physical cell.
      std::vector<Tensor<1, dim>> transformed_shape_values(n_q_points);
      for (unsigned int dof = 0; dof < this->n_dofs_per_cell(); ++dof)
        {
          const unsigned int first =
            data.shape_function_to_row_table[dof * this->n_components() +
                                             this->get_nonzero_components(dof)
                                               .first_selected_component()];

          mapping.transform(make_array_view(fe_data.shape_values[dof]),
                            mapping_covariant,
                            mapping_internal,
                            make_array_view(transformed_shape_values));
          for (unsigned int q = 0; q < n_q_points; ++q)
            {
              for (unsigned int d = 0; d < dim; ++d)
                {
                  data.shape_values(first + d, q) =
                    transformed_shape_values[q][d];
                }
            }
        }
    }

  if (flags & update_gradients)
    {
      // Now have all shape_grads stored on the reference cell.
      // Must now transform to the physical cell.
      std::vector<Tensor<2, dim>> input(n_q_points);
      std::vector<Tensor<2, dim>> transformed_shape_grads(n_q_points);
      for (unsigned int dof = 0; dof < this->n_dofs_per_cell(); ++dof)
        {
          for (unsigned int q = 0; q < n_q_points; ++q)
            {
              input[q] = fe_data.shape_grads[dof][q];
            }
          mapping.transform(make_array_view(input),
                            mapping_covariant_gradient,
                            mapping_internal,
                            make_array_view(transformed_shape_grads));

          const unsigned int first =
            data.shape_function_to_row_table[dof * this->n_components() +
                                             this->get_nonzero_components(dof)
                                               .first_selected_component()];

          for (unsigned int q = 0; q < n_q_points; ++q)
            {
              for (unsigned int d1 = 0; d1 < dim; ++d1)
                {
                  for (unsigned int d2 = 0; d2 < dim; ++d2)
                    {
                      transformed_shape_grads[q][d1] -=
                        data.shape_values(first + d2, q) *
                        mapping_data.jacobian_pushed_forward_grads[q][d2][d1];
                    }
                }
            }

          for (unsigned int q = 0; q < n_q_points; ++q)
            {
              for (unsigned int d = 0; d < dim; ++d)
                {
                  data.shape_gradients[first + d][q] =
                    transformed_shape_grads[q][d];
                }
            }
        }
    }

  if (flags & update_hessians)
    {
      // Now have all shape_grads stored on the reference cell.
      // Must now transform to the physical cell.
      std::vector<Tensor<3, dim>> input(n_q_points);
      std::vector<Tensor<3, dim>> transformed_shape_hessians(n_q_points);
      for (unsigned int dof = 0; dof < this->n_dofs_per_cell(); ++dof)
        {
          for (unsigned int q = 0; q < n_q_points; ++q)
            {
              input[q] = fe_data.shape_hessians[dof][q];
            }
          mapping.transform(make_array_view(input),
                            mapping_covariant_hessian,
                            mapping_internal,
                            make_array_view(transformed_shape_hessians));

          const unsigned int first =
            data.shape_function_to_row_table[dof * this->n_components() +
                                             this->get_nonzero_components(dof)
                                               .first_selected_component()];

          for (unsigned int q = 0; q < n_q_points; ++q)
            {
              for (unsigned int d1 = 0; d1 < dim; ++d1)
                {
                  for (unsigned int d2 = 0; d2 < dim; ++d2)
                    {
                      for (unsigned int d3 = 0; d3 < dim; ++d3)
                        {
                          for (unsigned int d4 = 0; d4 < dim; ++d4)
                            {
                              transformed_shape_hessians[q][d1][d3][d4] -=
                                (data.shape_values(first + d2, q) *
                                 mapping_data
                                   .jacobian_pushed_forward_2nd_derivatives
                                     [q][d2][d1][d3][d4]) +
                                (data.shape_gradients[first + d1][q][d2] *
                                 mapping_data
                                   .jacobian_pushed_forward_grads[q][d2][d3]
                                                                 [d4]) +
                                (data.shape_gradients[first + d2][q][d3] *
                                 mapping_data
                                   .jacobian_pushed_forward_grads[q][d2][d1]
                                                                 [d4]) +
                                (data.shape_gradients[first + d2][q][d4] *
                                 mapping_data
                                   .jacobian_pushed_forward_grads[q][d2][d3]
                                                                 [d1]);
                            }
                        }
                    }
                }
            }

          for (unsigned int q = 0; q < n_q_points; ++q)
            {
              for (unsigned int d = 0; d < dim; ++d)
                {
                  data.shape_hessians[first + d][q] =
                    transformed_shape_hessians[q][d];
                }
            }
        }
    }
}



template <int dim, int spacedim>
void
FE_NedelecSZ<dim, spacedim>::fill_fe_face_values(
  const typename Triangulation<dim, dim>::cell_iterator &cell,
  const unsigned int                                     face_no,
  const hp::QCollection<dim - 1> &                       quadrature,
  const Mapping<dim, dim> &                              mapping,
  const typename Mapping<dim, dim>::InternalDataBase &   mapping_internal,
  const internal::FEValuesImplementation::MappingRelatedData<dim, dim>
    &                                                       mapping_data,
  const typename FiniteElement<dim, dim>::InternalDataBase &fe_internal,
  dealii::internal::FEValuesImplementation::FiniteElementRelatedData<dim, dim>
    &data) const
{
  AssertDimension(quadrature.size(), 1);

  // Note for future improvement:
  // We don't have the full quadrature - should use QProjector to create the 2d
  // quadrature.
  //
  // For now I am effectively generating all of the shape function vals/grads,
  // etc. On all quad points on all faces and then only using them for one face.
  // This is obviously inefficient. I should cache the cell number and cache
  // all of the shape_values/gradients etc and then reuse them for each face.

  // convert data object to internal
  // data for this class. fails with
  // an exception if that is not
  // possible
  Assert(dynamic_cast<const InternalData *>(&fe_internal) != nullptr,
         ExcInternalError());
  const InternalData &fe_data = static_cast<const InternalData &>(fe_internal);

  // Now update the edge-based DoFs, which depend on the cell.
  // This will fill in the missing items in the InternalData
  // (fe_internal/fe_data) which was not filled in by get_data.
  fill_edge_values(cell,
                   QProjector<dim>::project_to_all_faces(this->reference_cell(),
                                                         quadrature[0]),
                   fe_data);
  if (dim == 3 && this->degree > 1)
    {
      fill_face_values(cell,
                       QProjector<dim>::project_to_all_faces(
                         this->reference_cell(), quadrature[0]),
                       fe_data);
    }

  const UpdateFlags  flags(fe_data.update_each);
  const unsigned int n_q_points = quadrature[0].size();
  const auto         offset =
    QProjector<dim>::DataSetDescriptor::face(this->reference_cell(),
                                             face_no,
                                             cell->face_orientation(face_no),
                                             cell->face_flip(face_no),
                                             cell->face_rotation(face_no),
                                             n_q_points);

  if (flags & update_values)
    {
      // Now have all shape_values stored on the reference cell.
      // Must now transform to the physical cell.
      std::vector<Tensor<1, dim>> transformed_shape_values(n_q_points);
      for (unsigned int dof = 0; dof < this->n_dofs_per_cell(); ++dof)
        {
          mapping.transform(make_array_view(fe_data.shape_values[dof],
                                            offset,
                                            n_q_points),
                            mapping_covariant,
                            mapping_internal,
                            make_array_view(transformed_shape_values));

          const unsigned int first =
            data.shape_function_to_row_table[dof * this->n_components() +
                                             this->get_nonzero_components(dof)
                                               .first_selected_component()];

          for (unsigned int q = 0; q < n_q_points; ++q)
            {
              for (unsigned int d = 0; d < dim; ++d)
                {
                  data.shape_values(first + d, q) =
                    transformed_shape_values[q][d];
                }
            }
        }
    }
  if (flags & update_gradients)
    {
      // Now have all shape_grads stored on the reference cell.
      // Must now transform to the physical cell.
      std::vector<Tensor<2, dim>> input(n_q_points);
      std::vector<Tensor<2, dim>> transformed_shape_grads(n_q_points);
      for (unsigned int dof = 0; dof < this->n_dofs_per_cell(); ++dof)
        {
          for (unsigned int q = 0; q < n_q_points; ++q)
            {
              input[q] = fe_data.shape_grads[dof][offset + q];
            }
          mapping.transform(input,
                            mapping_covariant_gradient,
                            mapping_internal,
                            make_array_view(transformed_shape_grads));

          const unsigned int first =
            data.shape_function_to_row_table[dof * this->n_components() +
                                             this->get_nonzero_components(dof)
                                               .first_selected_component()];

          for (unsigned int q = 0; q < n_q_points; ++q)
            {
              for (unsigned int d1 = 0; d1 < dim; ++d1)
                {
                  for (unsigned int d2 = 0; d2 < dim; ++d2)
                    {
                      transformed_shape_grads[q][d1] -=
                        data.shape_values(first + d2, q) *
                        mapping_data.jacobian_pushed_forward_grads[q][d2][d1];
                    }
                }
            }

          for (unsigned int q = 0; q < n_q_points; ++q)
            {
              for (unsigned int d = 0; d < dim; ++d)
                {
                  data.shape_gradients[first + d][q] =
                    transformed_shape_grads[q][d];
                }
            }
        }
    }
  if (flags & update_hessians)
    {
      // Now have all shape_grads stored on the reference cell.
      // Must now transform to the physical cell.
      std::vector<Tensor<3, dim>> input(n_q_points);
      std::vector<Tensor<3, dim>> transformed_shape_hessians(n_q_points);
      for (unsigned int dof = 0; dof < this->n_dofs_per_cell(); ++dof)
        {
          for (unsigned int q = 0; q < n_q_points; ++q)
            input[q] = fe_data.shape_hessians[dof][offset + q];

          mapping.transform(input,
                            mapping_covariant_hessian,
                            mapping_internal,
                            make_array_view(transformed_shape_hessians));

          const unsigned int first =
            data.shape_function_to_row_table[dof * this->n_components() +
                                             this->get_nonzero_components(dof)
                                               .first_selected_component()];

          for (unsigned int q = 0; q < n_q_points; ++q)
            {
              for (unsigned int d1 = 0; d1 < dim; ++d1)
                {
                  for (unsigned int d2 = 0; d2 < dim; ++d2)
                    {
                      for (unsigned int d3 = 0; d3 < dim; ++d3)
                        {
                          for (unsigned int d4 = 0; d4 < dim; ++d4)
                            {
                              transformed_shape_hessians[q][d1][d3][d4] -=
                                (data.shape_values(first + d2, q) *
                                 mapping_data
                                   .jacobian_pushed_forward_2nd_derivatives
                                     [q][d2][d1][d3][d4]) +
                                (data.shape_gradients[first + d1][q][d2] *
                                 mapping_data
                                   .jacobian_pushed_forward_grads[q][d2][d3]
                                                                 [d4]) +
                                (data.shape_gradients[first + d2][q][d3] *
                                 mapping_data
                                   .jacobian_pushed_forward_grads[q][d2][d1]
                                                                 [d4]) +
                                (data.shape_gradients[first + d2][q][d4] *
                                 mapping_data
                                   .jacobian_pushed_forward_grads[q][d2][d3]
                                                                 [d1]);
                            }
                        }
                    }
                }
            }

          for (unsigned int q = 0; q < n_q_points; ++q)
            {
              for (unsigned int d = 0; d < dim; ++d)
                {
                  data.shape_hessians[first + d][q] =
                    transformed_shape_hessians[q][d];
                }
            }
        }
    }
}



template <int dim, int spacedim>
void
FE_NedelecSZ<dim, spacedim>::fill_fe_subface_values(
  const typename Triangulation<dim, dim>::cell_iterator & /*cell*/,
  const unsigned int /*face_no*/,
  const unsigned int /*sub_no*/,
  const Quadrature<dim - 1> & /*quadrature*/,
  const Mapping<dim, dim> & /*mapping*/,
  const typename Mapping<dim, dim>::InternalDataBase & /*mapping_internal*/,
  const internal::FEValuesImplementation::MappingRelatedData<dim, dim>
    & /*mapping_data*/,
  const typename FiniteElement<dim, dim>::InternalDataBase & /*fe_internal*/,
  dealii::internal::FEValuesImplementation::FiniteElementRelatedData<dim, dim>
    & /*data*/) const
{
  Assert(false, ExcNotImplemented());
}



template <int dim, int spacedim>
UpdateFlags
FE_NedelecSZ<dim, spacedim>::requires_update_flags(
  const UpdateFlags flags) const
{
  UpdateFlags out = update_default;

  if (flags & update_values)
    out |= update_values | update_covariant_transformation;

  if (flags & update_gradients)
    out |= update_gradients | update_values |
           update_jacobian_pushed_forward_grads |
           update_covariant_transformation;

  if (flags & update_hessians)
    out |= update_hessians | update_values | update_gradients |
           update_jacobian_pushed_forward_grads |
           update_jacobian_pushed_forward_2nd_derivatives |
           update_covariant_transformation;

  return out;
}



template <int dim, int spacedim>
std::string
FE_NedelecSZ<dim, spacedim>::get_name() const
{
  // note that the FETools::get_fe_by_name function depends on the particular
  // format of the string this function returns, so they have to be kept in sync
  std::ostringstream namebuf;
  namebuf << "FE_NedelecSZ<" << dim << ">(" << this->degree - 1 << ")";

  return namebuf.str();
}



template <int dim, int spacedim>
std::unique_ptr<FiniteElement<dim, dim>>
FE_NedelecSZ<dim, spacedim>::clone() const
{
  return std::make_unique<FE_NedelecSZ<dim, spacedim>>(*this);
}



template <int dim, int spacedim>
std::vector<unsigned int>
FE_NedelecSZ<dim, spacedim>::get_dpo_vector(const unsigned int degree)
{
  // internal function to return a vector of "dofs per object"
  // where the objects inside the vector refer to:
  // 0 = vertex
  // 1 = edge
  // 2 = face (which is a cell in 2d)
  // 3 = cell
  std::vector<unsigned int> dpo;

  dpo.push_back(0);
  dpo.push_back(degree + 1);
  if (dim > 1)
    dpo.push_back(2 * degree * (degree + 1));
  if (dim > 2)
    dpo.push_back(3 * degree * degree * (degree + 1));

  return dpo;
}



template <int dim, int spacedim>
unsigned int
FE_NedelecSZ<dim, spacedim>::compute_num_dofs(const unsigned int degree) const
{
  // Internal function to compute the number of DoFs
  // for a given dimension & polynomial order.
  switch (dim)
    {
      case 2:
        return 2 * (degree + 1) * (degree + 2);

      case 3:
        return 3 * (degree + 1) * (degree + 2) * (degree + 2);

      default:
        {
          Assert(false, ExcNotImplemented());
          return 0;
        }
    }
}



template <int dim, int spacedim>
void
FE_NedelecSZ<dim, spacedim>::create_polynomials(const unsigned int degree)
{
  // fill the 1d polynomials vector:
  IntegratedLegendrePolynomials =
    IntegratedLegendreSZ::generate_complete_basis(degree + 1);
}



// explicit instantiations
#include "fe_nedelec_sz.inst"

DEAL_II_NAMESPACE_CLOSE
