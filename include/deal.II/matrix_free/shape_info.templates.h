// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2012 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_matrix_free_shape_info_templates_h
#define dealii_matrix_free_shape_info_templates_h


#include <deal.II/base/config.h>

#include <deal.II/base/memory_consumption.h>
#include <deal.II/base/polynomial.h>
#include <deal.II/base/polynomials_piecewise.h>
#include <deal.II/base/qprojector.h>
#include <deal.II/base/tensor_product_polynomials.h>
#include <deal.II/base/utilities.h>

#include <deal.II/fe/fe.h>
#include <deal.II/fe/fe_dgp.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_nothing.h>
#include <deal.II/fe/fe_poly.h>
#include <deal.II/fe/fe_pyramid_p.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_q_dg0.h>
#include <deal.II/fe/fe_q_iso_q1.h>
#include <deal.II/fe/fe_raviart_thomas.h>
#include <deal.II/fe/fe_simplex_p.h>
#include <deal.II/fe/fe_simplex_p_bubbles.h>
#include <deal.II/fe/fe_tools.h>
#include <deal.II/fe/fe_wedge_p.h>

#include <deal.II/grid/reference_cell.h>

#include <deal.II/lac/householder.h>

#include <deal.II/matrix_free/shape_info.h>
#include <deal.II/matrix_free/util.h>


DEAL_II_NAMESPACE_OPEN


namespace internal
{
  namespace MatrixFreeFunctions
  {
    template <typename Number>
    UnivariateShapeData<Number>::UnivariateShapeData()
      : element_type(tensor_general)
      , fe_degree(0)
      , n_q_points_1d(0)
      , nodal_at_cell_boundaries(false)
    {}



    template <typename Number>
    Number
    get_first_array_element(const Number a)
    {
      return a;
    }



    template <typename Number, std::size_t width>
    Number
    get_first_array_element(const VectorizedArray<Number, width> a)
    {
      return a[0];
    }



    template <int dim, int spacedim>
    void
    get_element_type_specific_information(
      const FiniteElement<dim, spacedim> &fe_in,
      const FiniteElement<dim, spacedim> &fe,
      const unsigned int                  base_element_number,
      ElementType                        &element_type,
      std::vector<unsigned int>          &scalar_lexicographic,
      std::vector<unsigned int>          &lexicographic_numbering)
    {
      element_type = tensor_general;

      const auto fe_poly = dynamic_cast<const FE_Poly<dim, spacedim> *>(&fe);

      if (dynamic_cast<const FE_SimplexPoly<dim, spacedim> *>(&fe) != nullptr ||
          dynamic_cast<const FE_WedgePoly<dim, spacedim> *>(&fe) != nullptr ||
          dynamic_cast<const FE_PyramidPoly<dim, spacedim> *>(&fe) != nullptr)
        {
          scalar_lexicographic.resize(fe.n_dofs_per_cell());
          for (unsigned int i = 0; i < scalar_lexicographic.size(); ++i)
            scalar_lexicographic[i] = i;
          element_type = tensor_none;
        }
      else if (fe_poly != nullptr &&
               (dynamic_cast<const TensorProductPolynomials<dim> *>(
                  &fe_poly->get_poly_space()) != nullptr ||
                dynamic_cast<const TensorProductPolynomials<
                    dim,
                    Polynomials::PiecewisePolynomial<double>> *>(
                  &fe_poly->get_poly_space()) != nullptr))
        scalar_lexicographic = fe_poly->get_poly_space_numbering_inverse();
      else if (const auto fe_dgp =
                 dynamic_cast<const FE_DGP<dim, spacedim> *>(&fe))
        {
          scalar_lexicographic.resize(fe_dgp->n_dofs_per_cell());
          for (unsigned int i = 0; i < fe_dgp->n_dofs_per_cell(); ++i)
            scalar_lexicographic[i] = i;
          element_type = truncated_tensor;
        }
      else if (const auto fe_q_dg0 =
                 dynamic_cast<const FE_Q_DG0<dim, spacedim> *>(&fe))
        {
          scalar_lexicographic = fe_q_dg0->get_poly_space_numbering_inverse();
          element_type         = tensor_symmetric_plus_dg0;
        }
      else if (fe.n_dofs_per_cell() == 0)
        {
          // FE_Nothing case -> nothing to do here
        }
      else
        DEAL_II_NOT_IMPLEMENTED();

      // Finally store the renumbering into the respective field
      if (fe_in.n_components() == 1)
        lexicographic_numbering = scalar_lexicographic;
      else
        {
          // have more than one component, get the inverse permutation, invert
          // it, sort the components one by one, and invert back
          std::vector<unsigned int> scalar_inv =
            Utilities::invert_permutation(scalar_lexicographic);
          std::vector<unsigned int> lexicographic(
            fe_in.n_dofs_per_cell(), numbers::invalid_unsigned_int);
          unsigned int components_before = 0;
          for (unsigned int e = 0; e < base_element_number; ++e)
            components_before += fe_in.element_multiplicity(e);
          for (unsigned int comp = 0;
               comp < fe_in.element_multiplicity(base_element_number);
               ++comp)
            for (unsigned int i = 0; i < scalar_inv.size(); ++i)
              lexicographic[fe_in.component_to_system_index(
                comp + components_before, i)] =
                scalar_inv.size() * comp + scalar_inv[i];

          // invert numbering again. Need to do it manually because we might
          // have undefined blocks
          lexicographic_numbering.resize(fe_in.element_multiplicity(
                                           base_element_number) *
                                           fe.n_dofs_per_cell(),
                                         numbers::invalid_unsigned_int);
          for (unsigned int i = 0; i < lexicographic.size(); ++i)
            if (lexicographic[i] != numbers::invalid_unsigned_int)
              {
                AssertIndexRange(lexicographic[i],
                                 lexicographic_numbering.size());
                lexicographic_numbering[lexicographic[i]] = i;
              }
        }
    }



    template <int dim_to, int dim, int spacedim>
    std::unique_ptr<FiniteElement<dim_to, dim_to>>
    create_fe(const FiniteElement<dim, spacedim> &fe)
    {
      std::string fe_name = fe.get_name();

      Assert(
        fe_name.find("FESystem") == std::string::npos,
        ExcMessage(
          "This function can not accept FESystem but only base elements."));

      {
        const std::size_t template_starts = fe_name.find_first_of('<');
        Assert(fe_name[template_starts + 1] ==
                 (dim == 1 ? '1' : (dim == 2 ? '2' : '3')),
               ExcInternalError());
        fe_name[template_starts + 1] = std::to_string(dim_to)[0];
      }
      return FETools::get_fe_by_name<dim_to, dim_to>(fe_name);
    }


    // ----------------- actual ShapeInfo implementation --------------------

    template <typename Number>
    ShapeInfo<Number>::ShapeInfo()
      : element_type(tensor_general)
      , n_dimensions(0)
      , n_components(0)
      , n_q_points(0)
      , dofs_per_component_on_cell(0)
      , n_q_points_face(0)
      , dofs_per_component_on_face(0)
    {}



    template <typename Number>
    template <int dim, int spacedim, int dim_q>
    inline ShapeInfo<Number>::ShapeInfo(
      const Quadrature<dim_q>            &quad,
      const FiniteElement<dim, spacedim> &fe_in,
      const unsigned int                  base_element_number)
      : element_type(tensor_general)
      , n_dimensions(0)
      , n_components(0)
      , n_q_points(0)
      , dofs_per_component_on_cell(0)
      , n_q_points_face(0)
      , dofs_per_component_on_face(0)
    {
      reinit(quad, fe_in, base_element_number);
    }



    template <typename Number>
    template <int dim, int spacedim, int dim_q>
    void
    ShapeInfo<Number>::reinit(const Quadrature<dim_q>            &quad_in,
                              const FiniteElement<dim, spacedim> &fe_in,
                              const unsigned int base_element_number)
    {
      // ShapeInfo for RT elements. Here, data is of size 2 instead of 1.
      // data[0] is univariate_shape_data in normal direction and
      // data[1] is univariate_shape_data in tangential direction
      //
      if (dynamic_cast<const FE_RaviartThomasNodal<dim> *>(
            &fe_in.base_element(base_element_number)))
        {
          element_type = tensor_raviart_thomas;

          const auto quad = quad_in.get_tensor_basis()[0];

          const FiniteElement<dim, spacedim> &fe =
            fe_in.base_element(base_element_number);
          n_dimensions = dim;
          n_components = fe_in.n_components();

          data.resize(2);
          const unsigned int n_q_points_1d = quad.size();

          n_q_points      = Utilities::fixed_power<dim>(n_q_points_1d);
          n_q_points_face = Utilities::fixed_power<dim - 1>(n_q_points_1d);

          dofs_per_component_on_cell = fe_in.n_dofs_per_cell() / n_components;

          // NOTE dofs_per_component_on_face is in tangential direction!
          dofs_per_component_on_face =
            fe_in.n_dofs_per_face() + Utilities::pow(fe_in.degree, dim - 2);
          const unsigned int dofs_per_face_normal = fe_in.n_dofs_per_face();

          lexicographic_numbering =
            FE_RaviartThomas<dim>::get_lexicographic_numbering(fe_in.degree -
                                                               1);

          // To get the right shape_values of the RT element
          std::vector<unsigned int> lex_normal, lex_tangent;
          for (unsigned int i = 0; i < fe.degree; ++i)
            lex_tangent.push_back(i);

          lex_normal.push_back(0);
          for (unsigned int i = dofs_per_face_normal * 2 * dim;
               i < dofs_per_face_normal * 2 * dim + fe.degree - 1;
               ++i)
            lex_normal.push_back(i);
          lex_normal.push_back(dofs_per_face_normal);

          // 'direction' distinguishes between normal and tangential direction
          for (unsigned int direction = 0; direction < 2; ++direction)
            {
              data[direction].element_type  = tensor_raviart_thomas;
              data[direction].quadrature    = quad;
              data[direction].n_q_points_1d = n_q_points_1d;
              data[direction].fe_degree     = fe.degree - direction;
              const std::vector<unsigned int> &lexicographic =
                direction == 0 ? lex_normal : lex_tangent;

              data[direction].evaluate_shape_functions(fe,
                                                       quad,
                                                       lexicographic,
                                                       direction);
              data[direction].evaluate_collocation_space(fe,
                                                         quad,
                                                         lexicographic,
                                                         direction);
              data[direction].check_and_set_shapes_symmetric();
            }

          if (dim == 3)
            face_orientations_quad = compute_orientation_table(n_q_points_1d);

          return;
        }
      else if (quad_in.is_tensor_product() == false ||
               dynamic_cast<const FE_SimplexPoly<dim, spacedim> *>(
                 &fe_in.base_element(base_element_number)) != nullptr ||
               dynamic_cast<const FE_WedgePoly<dim, spacedim> *>(
                 &fe_in.base_element(base_element_number)) != nullptr ||
               dynamic_cast<const FE_PyramidPoly<dim, spacedim> *>(
                 &fe_in.base_element(base_element_number)) != nullptr)
        {
          // specialization for arbitrary finite elements and quadrature rules
          // as needed in the context, e.g., of simplices

          AssertDimension(dim, dim_q);

          const auto  quad           = Quadrature<dim>(quad_in);
          const auto &fe             = fe_in.base_element(base_element_number);
          n_dimensions               = dim;
          n_components               = fe_in.n_components();
          n_q_points                 = quad.size();
          dofs_per_component_on_cell = fe.n_dofs_per_cell();
          n_q_points_face            = 0; // not implemented yet
          dofs_per_component_on_face = 0; //

          Assert(fe.n_components() == 1,
                 ExcMessage(
                   "FEEvaluation only works for scalar finite elements."));

          data.resize(1);
          UnivariateShapeData<Number> &univariate_shape_data = data.front();
          data_access.reinit(n_dimensions, n_components);
          data_access.fill(&univariate_shape_data);

          // note: we cannot write `univariate_shape_data.quadrature = quad`,
          // since the quadrature rule within UnivariateShapeData expects
          // a 1d quadrature rule. However, in this case we are not able to
          // define that rule anyway so other code cannot use this information.

          univariate_shape_data.fe_degree     = fe.degree;
          univariate_shape_data.n_q_points_1d = quad.size();

          if ((fe.n_dofs_per_cell() == 0) || (quad.empty()))
            return;

          // grant write access to common univariate shape data
          auto &shape_values      = univariate_shape_data.shape_values;
          auto &shape_values_face = univariate_shape_data.shape_values_face;
          auto &shape_gradients   = univariate_shape_data.shape_gradients;
          auto &shape_gradients_face =
            univariate_shape_data.shape_gradients_face;

          const unsigned int n_dofs = fe.n_dofs_per_cell();

          const unsigned int array_size = n_dofs * n_q_points;

          shape_values.resize_fast(array_size);
          shape_gradients.resize_fast(array_size * dim);

          for (unsigned int i = 0; i < n_dofs; ++i)
            for (unsigned int q = 0; q < n_q_points; ++q)
              {
                shape_values[i * n_q_points + q] =
                  fe.shape_value(i, quad.point(q));

                const auto grad = fe.shape_grad(i, quad.point(q));

                for (unsigned int d = 0; d < dim; ++d)
                  shape_gradients[i * dim * n_q_points + q * dim + d] = grad[d];
              }

          {
            const auto reference_cell = fe.reference_cell();

            const auto  temp      = get_face_quadrature_collection(quad, false);
            const auto &quad_face = temp.second;

            if (reference_cell != temp.first)
              {
                // TODO: this might happen if the quadrature rule and the
                // FE do not match
                this->n_q_points_face = 0;
              }
            else
              {
                this->n_q_points_face = quad_face[0].size();

                const unsigned int n_faces = temp.first.n_faces();

                n_q_points_faces.resize(n_faces);
                for (unsigned int i = 0; i < n_faces; ++i)
                  n_q_points_faces[i] =
                    quad_face[quad_face.size() == 1 ? 0 : i].size();

                unsigned int n_q_points_face_max = 0;

                for (unsigned int i = 0; i < quad_face.size(); ++i)
                  n_q_points_face_max =
                    std::max(n_q_points_face_max, quad_face[i].size());

                unsigned int n_max_vertices = 0;

                for (unsigned int face_no = 0; face_no < n_faces; ++face_no)
                  n_max_vertices = std::max(
                    n_max_vertices,
                    reference_cell.face_reference_cell(face_no).n_vertices());

                const auto projected_quad_face =
                  QProjector<dim>::project_to_all_faces(reference_cell,
                                                        quad_face);

                const unsigned int n_max_face_orientations =
                  dim == 2 ? 2 : (2 * n_max_vertices);

                shape_values_face.reinit({n_faces,
                                          n_max_face_orientations,
                                          n_dofs * n_q_points_face_max});

                shape_gradients_face.reinit(
                  {n_faces,
                   n_max_face_orientations,
                   dim * n_dofs * n_q_points_face_max});

                for (unsigned int f = 0; f < n_faces; ++f)
                  {
                    const unsigned int n_face_orientations =
                      reference_cell.n_face_orientations(f);

                    const unsigned int n_q_points_face =
                      quad_face[quad_face.size() == 1 ? 0 : f].size();

                    for (types::geometric_orientation orientation = 0;
                         orientation < n_face_orientations;
                         ++orientation)
                      {
                        const auto offset =
                          QProjector<dim>::DataSetDescriptor::face(
                            reference_cell, f, orientation, quad_face);

                        for (unsigned int i = 0; i < n_dofs; ++i)
                          for (unsigned int q = 0; q < n_q_points_face; ++q)
                            {
                              const auto &point =
                                projected_quad_face.point(q + offset);

                              shape_values_face(f,
                                                orientation,
                                                i * n_q_points_face + q) =
                                fe.shape_value(i, point);

                              const auto grad = fe.shape_grad(i, point);

                              for (unsigned int d = 0; d < dim; ++d)
                                shape_gradients_face(f,
                                                     orientation,
                                                     i * dim * n_q_points_face +
                                                       q * dim + d) = grad[d];
                            }
                      }
                  }
              }
          }

          // TODO: also fill shape_hessians, inverse_shape_values,
          //   shape_data_on_face, quadrature_data_on_face,
          //   values_within_subface, gradients_within_subface,
          //   hessians_within_subface

          // note: shape_gradients_collocation, shape_hessians_collocation,
          //  shape_values_eo, shape_gradients_eo, shape_hessians_eo,
          //  shape_gradients_collocation_eo, shape_hessians_collocation_eo,
          //  inverse_shape_values_eo cannot be filled

          std::vector<unsigned int> scalar_lexicographic;
          get_element_type_specific_information(fe_in,
                                                fe,
                                                base_element_number,
                                                element_type,
                                                scalar_lexicographic,
                                                lexicographic_numbering);

          univariate_shape_data.element_type = this->element_type;

          univariate_shape_data.nodal_at_cell_boundaries = true;

          const ReferenceCell reference_cell = fe.reference_cell();
          if (reference_cell.is_simplex())
            {
              if (dim == 2)
                dofs_per_component_on_face = fe.degree + 1;
              else
                dofs_per_component_on_face =
                  (fe.degree + 1) * (fe.degree + 2) / 2;

              face_to_cell_index_nodal.reinit(reference_cell.n_faces(),
                                              dofs_per_component_on_face);


              for (unsigned int face = 0; face < reference_cell.n_faces();
                   ++face)
                {
                  // first get info from reference cell, i.e. the linear case
                  unsigned int d = 0;
                  for (; d < dim; ++d)
                    face_to_cell_index_nodal[face][d] =
                      reference_cell.face_to_cell_vertices(
                        face, d, numbers::default_geometric_orientation);

                  // now fill the rest of the indices, start with the lines
                  if (fe.degree == 2)
                    for (; d < dofs_per_component_on_face; ++d)
                      face_to_cell_index_nodal[face][d] =
                        reference_cell.n_vertices() +
                        reference_cell.face_to_cell_lines(
                          face,
                          d - dim,
                          numbers::default_geometric_orientation);

                  // in the cubic case it is more complicated as more DoFs are
                  // on the lines
                  else if (fe.degree == 3)
                    {
                      for (unsigned int line = 0;
                           d < dofs_per_component_on_face - 1;
                           ++line, d += 2)
                        {
                          const unsigned int face_to_cell_lines =
                            reference_cell.face_to_cell_lines(
                              face,
                              line,
                              numbers::default_geometric_orientation);
                          // check the direction of the line
                          // is it 0 -> 1 or 1 -> 0
                          // as DoFs on the line are ordered differently
                          if (reference_cell.line_to_cell_vertices(
                                face_to_cell_lines, 0) ==
                              reference_cell.face_to_cell_vertices(
                                face,
                                line,
                                numbers::default_geometric_orientation))
                            {
                              face_to_cell_index_nodal[face][d] =
                                reference_cell.n_vertices() +
                                2 * face_to_cell_lines;
                              face_to_cell_index_nodal[face][d + 1] =
                                reference_cell.n_vertices() +
                                2 * face_to_cell_lines + 1;
                            }
                          else
                            {
                              face_to_cell_index_nodal[face][d + 1] =
                                reference_cell.n_vertices() +
                                2 * face_to_cell_lines;
                              face_to_cell_index_nodal[face][d] =
                                reference_cell.n_vertices() +
                                2 * face_to_cell_lines + 1;
                            }
                        }
                      //  in 3D we also need the DoFs on the quads
                      if (dim == 3)
                        {
                          face_to_cell_index_nodal
                            [face][dofs_per_component_on_face - 1] =
                              reference_cell.n_vertices() +
                              2 * reference_cell.n_lines() + face;
                        }
                    }
                  else if (fe.degree > 3)
                    DEAL_II_NOT_IMPLEMENTED();
                }
            }
          // TODO: set up face_to_cell_index_nodal, face_to_cell_index_hermite,
          //  face_orientations

          return;
        }

      const auto quad = quad_in.get_tensor_basis()[0];

      const FiniteElement<dim, spacedim> &fe =
        fe_in.base_element(base_element_number);
      n_dimensions = dim;
      n_components = fe_in.n_components();

      Assert(fe.n_components() == 1,
             ExcMessage("FEEvaluation only works for scalar finite elements."));

      // assuming isotropy of dimensions and components
      data.resize(1);
      UnivariateShapeData<Number> &univariate_shape_data = data.front();
      data_access.reinit(n_dimensions, n_components);
      data_access.fill(&univariate_shape_data);
      univariate_shape_data.quadrature    = quad;
      univariate_shape_data.fe_degree     = fe.degree;
      univariate_shape_data.n_q_points_1d = quad.size();

      if ((fe.n_dofs_per_cell() == 0) || (quad.empty()))
        return;

      const unsigned int fe_degree     = fe.degree;
      const unsigned int n_q_points_1d = quad.size();
      const unsigned int n_dofs_1d =
        std::min(fe.n_dofs_per_cell(), fe_degree + 1);

      // renumber (this is necessary for FE_Q, for example, since there the
      // vertex DoFs come first, which is incompatible with the lexicographic
      // ordering necessary to apply tensor products efficiently)
      std::vector<unsigned int> scalar_lexicographic;
      Assert(fe.n_components() == 1, ExcMessage("Expected a scalar element"));

      get_element_type_specific_information(fe_in,
                                            fe,
                                            base_element_number,
                                            element_type,
                                            scalar_lexicographic,
                                            lexicographic_numbering);

      n_q_points = Utilities::fixed_power<dim>(n_q_points_1d);
      n_q_points_face =
        (dim > 1 ? Utilities::fixed_power<dim - 1>(n_q_points_1d) : 1);
      dofs_per_component_on_cell = fe.n_dofs_per_cell();
      dofs_per_component_on_face =
        (dim > 1 ? Utilities::fixed_power<dim - 1>(fe_degree + 1) : 1);

      univariate_shape_data.evaluate_shape_functions(fe,
                                                     quad,
                                                     scalar_lexicographic,
                                                     0);

      if (dim > 1 && (dynamic_cast<const FE_Q<dim, spacedim> *>(&fe) ||
                      dynamic_cast<const FE_Q_iso_Q1<dim, spacedim> *>(&fe)))
        {
          auto &subface_interpolation_matrix_0 =
            univariate_shape_data.subface_interpolation_matrices[0];
          auto &subface_interpolation_matrix_1 =
            univariate_shape_data.subface_interpolation_matrices[1];
          auto &subface_interpolation_matrix_scalar_0 =
            univariate_shape_data.subface_interpolation_matrices_scalar[0];
          auto &subface_interpolation_matrix_scalar_1 =
            univariate_shape_data.subface_interpolation_matrices_scalar[1];

          const unsigned int nn = fe_degree + 1;
          subface_interpolation_matrix_0.resize(nn * nn);
          subface_interpolation_matrix_1.resize(nn * nn);
          subface_interpolation_matrix_scalar_0.resize(nn * nn);
          subface_interpolation_matrix_scalar_1.resize(nn * nn);

          const bool is_feq =
            dynamic_cast<const FE_Q<dim, spacedim> *>(&fe) != nullptr;

          std::vector<Point<1>> fe_q_points =
            is_feq ? QGaussLobatto<1>(nn).get_points() :
                     QIterated<1>(QTrapezoid<1>(), nn - 1).get_points();

          const std::vector<Polynomials::Polynomial<double>> poly_feq =
            Polynomials::generate_complete_Lagrange_basis(fe_q_points);

          const std::vector<Polynomials::PiecewisePolynomial<double>>
            poly_feq_iso_q1 =
              Polynomials::generate_complete_Lagrange_basis_on_subdivisions(nn -
                                                                              1,
                                                                            1);

          for (unsigned int i = 0, c = 0; i < nn; ++i)
            for (unsigned int j = 0; j < nn; ++j, ++c)
              {
                subface_interpolation_matrix_scalar_0[c] =
                  is_feq ? poly_feq[j].value(0.5 * fe_q_points[i][0]) :
                           poly_feq_iso_q1[j].value(0.5 * fe_q_points[i][0]);
                subface_interpolation_matrix_0[c] =
                  subface_interpolation_matrix_scalar_0[c];
                subface_interpolation_matrix_scalar_1[c] =
                  is_feq ?
                    poly_feq[j].value(0.5 + 0.5 * fe_q_points[i][0]) :
                    poly_feq_iso_q1[j].value(0.5 + 0.5 * fe_q_points[i][0]);
                subface_interpolation_matrix_1[c] =
                  subface_interpolation_matrix_scalar_1[c];
              }
        }

      univariate_shape_data.evaluate_collocation_space(fe,
                                                       quad,
                                                       scalar_lexicographic,
                                                       0);

      const auto &shape_data_on_face = univariate_shape_data.shape_data_on_face;

      if (element_type == tensor_general &&
          univariate_shape_data.check_and_set_shapes_symmetric())
        {
          if (dynamic_cast<const FE_Q_iso_Q1<dim, spacedim> *>(&fe) &&
              fe.tensor_degree() > 1)
            element_type = tensor_symmetric_no_collocation;
          else if (univariate_shape_data.check_shapes_collocation())
            element_type = tensor_symmetric_collocation;
          else
            element_type = tensor_symmetric;

          if (n_dofs_1d > 2 && element_type == tensor_symmetric)
            {
              // check if we are a Hermite type
              element_type = tensor_symmetric_hermite;
              for (unsigned int i = 1; i < n_dofs_1d; ++i)
                if (std::abs(get_first_array_element(
                      shape_data_on_face[0][i])) > 1e-12)
                  element_type = tensor_symmetric;
              for (unsigned int i = 2; i < n_dofs_1d; ++i)
                if (std::abs(get_first_array_element(
                      shape_data_on_face[0][n_dofs_1d + i])) > 1e-12)
                  element_type = tensor_symmetric;
            }
        }
      else if (element_type == tensor_symmetric_plus_dg0)
        univariate_shape_data.check_and_set_shapes_symmetric();

      univariate_shape_data.nodal_at_cell_boundaries = true;
      for (unsigned int i = 1; i < n_dofs_1d; ++i)
        if (std::abs(get_first_array_element(shape_data_on_face[0][i])) >
              1e-13 ||
            std::abs(get_first_array_element(shape_data_on_face[1][i - 1])) >
              1e-13)
          univariate_shape_data.nodal_at_cell_boundaries = false;

      if (univariate_shape_data.nodal_at_cell_boundaries == true)
        {
          face_to_cell_index_nodal.reinit(GeometryInfo<dim>::faces_per_cell,
                                          dofs_per_component_on_face);
          for (const auto f : GeometryInfo<dim>::face_indices())
            {
              const unsigned int direction = f / 2;
              const unsigned int stride =
                direction < dim - 1 ? (fe_degree + 1) : 1;
              int shift = 1;
              for (unsigned int d = 0; d < direction; ++d)
                shift *= fe_degree + 1;
              const unsigned int offset = (f % 2) * fe_degree * shift;

              if (direction == 0 || direction == dim - 1)
                for (unsigned int i = 0; i < dofs_per_component_on_face; ++i)
                  face_to_cell_index_nodal(f, i) = offset + i * stride;
              else
                // local coordinate system on faces 2 and 3 is zx in
                // deal.II, not xz as expected for tensor products -> adjust
                // that here
                for (unsigned int j = 0; j <= fe_degree; ++j)
                  for (unsigned int i = 0; i <= fe_degree; ++i)
                    {
                      const unsigned int ind =
                        offset + j * dofs_per_component_on_face + i;
                      AssertIndexRange(ind, dofs_per_component_on_cell);
                      const unsigned int l           = i * (fe_degree + 1) + j;
                      face_to_cell_index_nodal(f, l) = ind;
                    }
            }

          // face orientation for faces in 3d
          // (similar to MappingInfoStorage::QuadratureDescriptor::initialize)
          if (dim == 3)
            {
              face_orientations_dofs = compute_orientation_table(fe_degree + 1);
              face_orientations_quad = compute_orientation_table(n_q_points_1d);
            }
          else
            {
              face_orientations_dofs.reinit(1, 1);
              face_orientations_quad.reinit(1, 1);
            }
        }

      if (element_type == tensor_symmetric_hermite)
        {
          face_to_cell_index_hermite.reinit(GeometryInfo<dim>::faces_per_cell,
                                            2 * dofs_per_component_on_face);
          for (const auto f : GeometryInfo<dim>::face_indices())
            {
              const unsigned int direction = f / 2;
              const unsigned int stride =
                direction < dim - 1 ? (fe_degree + 1) : 1;
              int shift = 1;
              for (unsigned int d = 0; d < direction; ++d)
                shift *= fe_degree + 1;
              const unsigned int offset = (f % 2) * fe_degree * shift;
              if (f % 2 == 1)
                shift = -shift;

              if (direction == 0 || direction == dim - 1)
                for (unsigned int i = 0; i < dofs_per_component_on_face; ++i)
                  {
                    face_to_cell_index_hermite(f, 2 * i) = offset + i * stride;
                    face_to_cell_index_hermite(f, 2 * i + 1) =
                      offset + i * stride + shift;
                  }
              else
                // local coordinate system on faces 2 and 3 is zx in
                // deal.II, not xz as expected for tensor products -> adjust
                // that here
                for (unsigned int j = 0; j <= fe_degree; ++j)
                  for (unsigned int i = 0; i <= fe_degree; ++i)
                    {
                      const unsigned int ind =
                        offset + j * dofs_per_component_on_face + i;
                      AssertIndexRange(ind, dofs_per_component_on_cell);
                      const unsigned int l = i * (fe_degree + 1) + j;
                      face_to_cell_index_hermite(f, 2 * l)     = ind;
                      face_to_cell_index_hermite(f, 2 * l + 1) = ind + shift;
                    }
            }
        }

      univariate_shape_data.element_type = this->element_type;
    }



    template <int dim, int spacedim>
    Point<dim>
    get_unit_point(const FiniteElement<dim, spacedim> &fe,
                   const std::vector<unsigned int>    &lexicographic)
    {
      Point<dim> unit_point;
      // to evaluate 1d polynomials, evaluate along the line with the first
      // unit support point, assuming that fe.shape_value(0,unit_point) ==
      // 1. otherwise, need other entry point (e.g. generating a 1d element
      // by reading the name, as done before r29356)
      if (fe.has_support_points())
        unit_point = fe.get_unit_support_points()[lexicographic[0]];
      Assert(fe.n_dofs_per_cell() == 0 ||
               std::abs(
                 fe.shape_value_component(lexicographic[0], unit_point, 0) -
                 1) < 1e-13,
             ExcInternalError("Could not decode 1d shape functions for the "
                              "element " +
                              fe.get_name()));
      return unit_point;
    }



    template <typename Number>
    template <int dim, int spacedim>
    void
    UnivariateShapeData<Number>::evaluate_shape_functions(
      const FiniteElement<dim, spacedim> &fe,
      const Quadrature<1>                &quad,
      const std::vector<unsigned int>    &lexicographic,
      const unsigned int                  direction)
    {
      const unsigned int n_dofs_1d =
        std::min(fe.n_dofs_per_cell(), fe_degree + 1);

      const unsigned int array_size = n_dofs_1d * n_q_points_1d;
      shape_gradients.resize_fast(array_size);
      shape_values.resize_fast(array_size);
      shape_hessians.resize_fast(array_size);

      shape_data_on_face[0].resize(3 * n_dofs_1d);
      shape_data_on_face[1].resize(3 * n_dofs_1d);
      values_within_subface[0].resize(array_size);
      values_within_subface[1].resize(array_size);
      gradients_within_subface[0].resize(array_size);
      gradients_within_subface[1].resize(array_size);
      hessians_within_subface[0].resize(array_size);
      hessians_within_subface[1].resize(array_size);

      for (unsigned int i = 0; i < n_dofs_1d; ++i)
        {
          // need to reorder from hierarchical to lexicographic to get the
          // DoFs correct
          const unsigned int my_i = lexicographic[i];
          for (unsigned int q = 0; q < n_q_points_1d; ++q)
            {
              Point<dim> q_point = get_unit_point(fe, lexicographic);
              q_point[direction] = quad.get_points()[q][0];

              shape_values[i * n_q_points_1d + q] =
                fe.shape_value_component(my_i, q_point, 0);
              shape_gradients[i * n_q_points_1d + q] =
                fe.shape_grad_component(my_i, q_point, 0)[direction];
              shape_hessians[i * n_q_points_1d + q] =
                fe.shape_grad_grad_component(my_i,
                                             q_point,
                                             0)[direction][direction];

              // evaluate basis functions on the two 1d subfaces (i.e., at the
              // positions divided by one half and shifted by one half,
              // respectively)
              q_point[direction] *= 0.5;
              values_within_subface[0][i * n_q_points_1d + q] =
                fe.shape_value_component(my_i, q_point, 0);
              gradients_within_subface[0][i * n_q_points_1d + q] =
                fe.shape_grad_component(my_i, q_point, 0)[direction];
              hessians_within_subface[0][i * n_q_points_1d + q] =
                fe.shape_grad_grad_component(my_i,
                                             q_point,
                                             0)[direction][direction];
              q_point[direction] += 0.5;
              values_within_subface[1][i * n_q_points_1d + q] =
                fe.shape_value_component(my_i, q_point, 0);
              gradients_within_subface[1][i * n_q_points_1d + q] =
                fe.shape_grad_component(my_i, q_point, 0)[direction];
              hessians_within_subface[1][i * n_q_points_1d + q] =
                fe.shape_grad_grad_component(my_i,
                                             q_point,
                                             0)[direction][direction];
            }

          // evaluate basis functions on the 1d faces, i.e., in zero and one
          Point<dim> q_point       = get_unit_point(fe, lexicographic);
          q_point[direction]       = 0;
          shape_data_on_face[0][i] = fe.shape_value_component(my_i, q_point, 0);
          shape_data_on_face[0][i + n_dofs_1d] =
            fe.shape_grad_component(my_i, q_point, 0)[direction];
          shape_data_on_face[0][i + 2 * n_dofs_1d] =
            fe.shape_grad_grad_component(my_i,
                                         q_point,
                                         0)[direction][direction];
          q_point[direction]       = 1;
          shape_data_on_face[1][i] = fe.shape_value_component(my_i, q_point, 0);
          shape_data_on_face[1][i + n_dofs_1d] =
            fe.shape_grad_component(my_i, q_point, 0)[direction];
          shape_data_on_face[1][i + 2 * n_dofs_1d] =
            fe.shape_grad_grad_component(my_i,
                                         q_point,
                                         0)[direction][direction];
        }
    }



    template <typename Number>
    template <int dim, int spacedim>
    void
    UnivariateShapeData<Number>::evaluate_collocation_space(
      const FiniteElement<dim, spacedim> &fe,
      const Quadrature<1>                &quad,
      const std::vector<unsigned int>    &lexicographic,
      const unsigned int                  direction)
    {
      const unsigned int n_dofs_1d =
        std::min(fe.n_dofs_per_cell(), fe_degree + 1);

      // get gradient and Hessian transformation matrix for the polynomial
      // space associated with the quadrature rule (collocation space). We
      // need to avoid the case with more than a few hundreds of quadrature
      // points when the Lagrange polynomials might underflow. Note that 200
      // is not an exact value, as different quadrature formulas behave
      // slightly differently, but 200 has been observed to be low enough for
      // all common quadrature formula types. For QGauss, the actual limit is
      // 517 points, for example.
      if (n_q_points_1d >= 200)
        return;

      shape_gradients_collocation.resize(n_q_points_1d * n_q_points_1d);
      shape_hessians_collocation.resize(n_q_points_1d * n_q_points_1d);
      const std::vector<Polynomials::Polynomial<double>> poly_coll =
        Polynomials::generate_complete_Lagrange_basis(quad.get_points());
      std::array<double, 3> values;
      for (unsigned int i = 0; i < n_q_points_1d; ++i)
        for (unsigned int q = 0; q < n_q_points_1d; ++q)
          {
            poly_coll[i].value(quad.get_points()[q][0], 2, values.data());
            shape_gradients_collocation[i * n_q_points_1d + q] = values[1];
            shape_hessians_collocation[i * n_q_points_1d + q]  = values[2];
          }

      // compute the inverse shape functions in three steps: we first
      // change from the given quadrature formula and the associated
      // Lagrange polynomials to the Lagrange polynomials at quadrature
      // points. in this basis, we can then perform the second step, which
      // is the computation of a projection matrix from the potentially
      // higher polynomial degree associated to the quadrature points to a
      // polynomial space of degree equal to the degree of the given
      // elements. in the third step, we change from the Lagrange
      // polynomials in the Gauss quadrature points to the polynomial
      // space of the given element

      // step 1: change basis from the Lagrange polynomials at the given
      // quadrature points to the Lagrange basis at Gauss quadrature
      // points. this is often the identity operation as we often compute
      // with Gaussian quadrature, but not necessarily so
      QGauss<1>          quad_gauss(n_q_points_1d);
      FullMatrix<double> transform_to_gauss(n_q_points_1d, n_q_points_1d);
      for (unsigned int i = 0; i < n_q_points_1d; ++i)
        for (unsigned int j = 0; j < n_q_points_1d; ++j)
          transform_to_gauss(i, j) = poly_coll[j].value(quad_gauss.point(i)[0]);

      // step 2: computation for the projection (in reference coordinates)
      // from higher to lower polynomial degree
      //
      // loop over quadrature points, multiply by q-weight on high degree
      // integrate loop going from high degree to low degree loop over new
      // points, multiply by inverse q-weight on low degree
      //
      // This projection step is for the special case of Lagrange
      // polynomials where most of the interpolation matrices are unit
      // matrices when applying the inverse mass matrix, so we do not need
      // to compute much.
      QGauss<1> quad_project(n_dofs_1d);
      const std::vector<Polynomials::Polynomial<double>> poly_project =
        Polynomials::generate_complete_Lagrange_basis(
          quad_project.get_points());

      FullMatrix<double> project_gauss(n_dofs_1d, n_q_points_1d);

      for (unsigned int i = 0; i < n_dofs_1d; ++i)
        for (unsigned int q = 0; q < n_q_points_1d; ++q)
          project_gauss(i, q) =
            poly_project[i].value(quad_gauss.get_points()[q][0]) *
            (quad_gauss.weight(q) / quad_project.weight(i));
      FullMatrix<double> project_to_dof_space(n_dofs_1d, n_q_points_1d);
      project_gauss.mmult(project_to_dof_space, transform_to_gauss);

      // step 3: change the basis back to the given finite element
      // space. we can use a shortcut for elements that define support
      // points, in which case we can evaluate the Lagrange polynomials of
      // the Gauss quadrature in those points. this will give more
      // accurate results than the inversion of a matrix. for more general
      // polynomial spaces, we must invert a matrix of a Vandermonde type,
      // which we do by a Householder transformation to keep roundoff
      // errors low.
      inverse_shape_values.resize_fast(shape_values.size());
      FullMatrix<double> transform_from_gauss(n_dofs_1d, n_dofs_1d);
      if (fe.has_support_points())
        {
          for (unsigned int i = 0; i < n_dofs_1d; ++i)
            for (unsigned int j = 0; j < n_dofs_1d; ++j)
              transform_from_gauss(i, j) = poly_project[j].value(
                fe.get_unit_support_points()[lexicographic[i]][0]);
          FullMatrix<double> result(n_dofs_1d, n_q_points_1d);
          transform_from_gauss.mmult(result, project_to_dof_space);

          // set very small entries to zero - we are in reference space
          // with normalized numbers, so this is straight-forward to check
          // here
          for (unsigned int i = 0; i < n_dofs_1d; ++i)
            for (unsigned int q = 0; q < n_q_points_1d; ++q)
              inverse_shape_values[i * n_q_points_1d + q] =
                std::abs(result(i, q)) < 1e-15 ? 0 : result(i, q);
        }
      else
        {
          for (unsigned int i = 0; i < n_dofs_1d; ++i)
            for (unsigned int j = 0; j < n_dofs_1d; ++j)
              {
                Point<dim> q_point = get_unit_point(fe, lexicographic);
                q_point[direction] = quad_project.point(i)[0];

                transform_from_gauss(i, j) =
                  fe.shape_value_component(lexicographic[j], q_point, 0);
              }
          Householder<double> H(transform_from_gauss);
          Vector<double>      in(n_dofs_1d), out(n_dofs_1d);
          for (unsigned int q = 0; q < n_q_points_1d; ++q)
            {
              for (unsigned int i = 0; i < n_dofs_1d; ++i)
                in(i) = project_to_dof_space(i, q);
              H.least_squares(out, in);
              for (unsigned int i = 0; i < n_dofs_1d; ++i)
                inverse_shape_values[i * n_q_points_1d + q] =
                  std::abs(out(i)) < 1e-15 ? 0. : out(i);
            }
        }
      quadrature_data_on_face[0].resize(quad.size() * 3);
      quadrature_data_on_face[1].resize(quad.size() * 3);

      for (unsigned int i = 0; i < quad.size(); ++i)
        {
          std::array<double, 3> values;
          poly_coll[i].value(0.0, 2, values.data());
          for (unsigned int d = 0; d < 3; ++d)
            quadrature_data_on_face[0][i + d * quad.size()] = values[d];
          poly_coll[i].value(1.0, 2, values.data());
          for (unsigned int d = 0; d < 3; ++d)
            quadrature_data_on_face[1][i + d * quad.size()] = values[d];
        }
    }



    template <typename Number>
    bool
    UnivariateShapeData<Number>::check_and_set_shapes_symmetric()
    {
      const double zero_tol =
        std::is_same_v<Number, double> == true ? 1e-12 : 1e-7;
      // symmetry for values
      const unsigned int n_dofs_1d = fe_degree + 1;
      for (unsigned int i = 0; i < (n_dofs_1d + 1) / 2; ++i)
        for (unsigned int j = 0; j < n_q_points_1d; ++j)
          if (std::abs(get_first_array_element(
                shape_values[i * n_q_points_1d + j] -
                shape_values[(n_dofs_1d - i) * n_q_points_1d - j - 1])) >
              std::max(zero_tol,
                       zero_tol * std::abs(get_first_array_element(
                                    shape_values[i * n_q_points_1d + j]))))
            return false;

      // shape values should be zero at x=0.5 for all basis functions except
      // for the middle one for degrees of 4 and higher
      if (n_dofs_1d > 3 && n_q_points_1d % 2 == 1 && n_dofs_1d % 2 == 1)
        {
          for (unsigned int i = 0; i < n_dofs_1d / 2; ++i)
            if (std::abs(get_first_array_element(
                  shape_values[i * n_q_points_1d + n_q_points_1d / 2])) >
                zero_tol)
              return false;
        }

      // skew-symmetry for gradient, zero of middle basis function in middle
      // quadrature point. Multiply tolerance by degree of the element to
      // the power of 1.5 to get a suitable gradient scaling
      const double zero_tol_gradient =
        zero_tol * std::sqrt(fe_degree + 1.) * (fe_degree + 1);
      for (unsigned int i = 0; i < (n_dofs_1d + 1) / 2; ++i)
        for (unsigned int j = 0; j < n_q_points_1d; ++j)
          if (std::abs(get_first_array_element(
                shape_gradients[i * n_q_points_1d + j] +
                shape_gradients[(n_dofs_1d - i) * n_q_points_1d - j - 1])) >
              zero_tol_gradient)
            return false;
      if (n_dofs_1d % 2 == 1 && n_q_points_1d % 2 == 1)
        if (std::abs(get_first_array_element(
              shape_gradients[(n_dofs_1d / 2) * n_q_points_1d +
                              (n_q_points_1d / 2)])) > zero_tol_gradient)
          return false;

      // symmetry for Hessian. Multiply tolerance by degree^3 of the element
      // to get a suitable Hessian scaling
      const double zero_tol_hessian =
        zero_tol * (fe_degree + 1) * (fe_degree + 1) * (fe_degree + 1);
      for (unsigned int i = 0; i < (n_dofs_1d + 1) / 2; ++i)
        for (unsigned int j = 0; j < n_q_points_1d; ++j)
          if (std::abs(get_first_array_element(
                shape_hessians[i * n_q_points_1d + j] -
                shape_hessians[(n_dofs_1d - i) * n_q_points_1d - j - 1])) >
              zero_tol_hessian)
            return false;

      auto convert_to_eo = [](const AlignedVector<Number> &array,
                              const unsigned               n_rows,
                              const unsigned               n_cols) {
        const unsigned int    stride = (n_cols + 1) / 2;
        AlignedVector<Number> array_eo(n_rows * stride);
        for (unsigned int i = 0; i < n_rows / 2; ++i)
          for (unsigned int q = 0; q < stride; ++q)
            {
              array_eo[i * stride + q] =
                0.5 *
                (array[i * n_cols + q] + array[i * n_cols + n_cols - 1 - q]);
              array_eo[(n_rows - 1 - i) * stride + q] =
                0.5 *
                (array[i * n_cols + q] - array[i * n_cols + n_cols - 1 - q]);
            }
        if ((n_rows - 1) % 2 == 0)
          for (unsigned int q = 0; q < stride; ++q)
            {
              array_eo[(n_rows - 1) / 2 * stride + q] =
                array[((n_rows - 1) / 2) * n_cols + q];
            }

        return array_eo;
      };

      shape_values_eo =
        convert_to_eo(shape_values, fe_degree + 1, n_q_points_1d);
      shape_gradients_eo =
        convert_to_eo(shape_gradients, fe_degree + 1, n_q_points_1d);
      shape_hessians_eo =
        convert_to_eo(shape_hessians, fe_degree + 1, n_q_points_1d);

      // Avoid underflow of Lagrange polynomials on typical quadrature
      // formulas (see also above where shape_gradients_collocation and
      // shape_hessians_collocation is set up).
      if (n_q_points_1d < 200)
        {
          shape_gradients_collocation_eo =
            convert_to_eo(shape_gradients_collocation,
                          n_q_points_1d,
                          n_q_points_1d);
          shape_hessians_collocation_eo =
            convert_to_eo(shape_hessians_collocation,
                          n_q_points_1d,
                          n_q_points_1d);
          inverse_shape_values_eo =
            convert_to_eo(inverse_shape_values, fe_degree + 1, n_q_points_1d);
        }

      return true;
    }



    template <typename Number>
    bool
    UnivariateShapeData<Number>::check_shapes_collocation() const
    {
      if (fe_degree + 1 != n_q_points_1d)
        return false;

      const double zero_tol =
        std::is_same_v<Number, double> == true ? 1e-12 : 1e-7;
      // check: identity operation for shape values
      const unsigned int n_points_1d = fe_degree + 1;
      for (unsigned int i = 0; i < n_points_1d; ++i)
        for (unsigned int j = 0; j < n_points_1d; ++j)
          if (i != j)
            {
              if (std::abs(get_first_array_element(
                    shape_values[i * n_points_1d + j])) > zero_tol)
                return false;
            }
          else
            {
              if (std::abs(
                    get_first_array_element(shape_values[i * n_points_1d + j]) -
                    1.) > zero_tol)
                return false;
            }
      return true;
    }



    template <typename Number>
    template <int dim, int spacedim>
    bool
    ShapeInfo<Number>::is_supported(const FiniteElement<dim, spacedim> &fe)
    {
      if (dynamic_cast<const FE_RaviartThomasNodal<dim> *>(&fe))
        return true;

      for (unsigned int base = 0; base < fe.n_base_elements(); ++base)
        {
          const FiniteElement<dim, spacedim> *fe_ptr = &(fe.base_element(base));
          if (fe_ptr->n_components() != 1)
            return false;

          // then check if the base element is supported or not
          if (dynamic_cast<const FE_Poly<dim, spacedim> *>(fe_ptr) != nullptr)
            {
              const FE_Poly<dim, spacedim> *fe_poly_ptr =
                dynamic_cast<const FE_Poly<dim, spacedim> *>(fe_ptr);
              // Simplices are a special case since the polynomial family is not
              // indicative of their support
              if (dynamic_cast<const FE_SimplexPoly<dim, spacedim> *>(
                    fe_poly_ptr) != nullptr ||
                  dynamic_cast<const FE_WedgePoly<dim, spacedim> *>(
                    fe_poly_ptr) != nullptr ||
                  dynamic_cast<const FE_PyramidPoly<dim, spacedim> *>(
                    fe_poly_ptr) != nullptr)
                return true;

              if (dynamic_cast<const TensorProductPolynomials<
                      dim,
                      Polynomials::Polynomial<double>> *>(
                    &fe_poly_ptr->get_poly_space()) == nullptr &&
                  dynamic_cast<const TensorProductPolynomials<
                      dim,
                      Polynomials::PiecewisePolynomial<double>> *>(
                    &fe_poly_ptr->get_poly_space()) == nullptr &&
                  dynamic_cast<const FE_DGP<dim, spacedim> *>(fe_ptr) ==
                    nullptr &&
                  dynamic_cast<const FE_Q_DG0<dim, spacedim> *>(fe_ptr) ==
                    nullptr)
                return false;
            }
          else if (dynamic_cast<const FE_Nothing<dim, spacedim> *>(fe_ptr) !=
                   nullptr)
            return true;
          else
            return false;
        }

      // if we arrived here, all base elements were supported so we can
      // support the present element
      return true;
    }



    template <typename Number>
    Table<2, unsigned int>
    ShapeInfo<Number>::compute_orientation_table(const unsigned int n)
    {
      Table<2, unsigned int> face_orientations(8, n * n);
      for (unsigned int j = 0, i = 0; j < n; ++j)
        for (unsigned int k = 0; k < n; ++k, ++i)
          {
            // face_orientation=true,  face_rotation=false, face_flip=false
            face_orientations[0][i] = i;
            // face_orientation=false, face_rotation=false, face_flip=false
            face_orientations[1][i] = j + k * n;
            // face_orientation=true,  face_rotation=true, face_flip=false
            face_orientations[2][i] = j + (n - 1 - k) * n;
            // face_orientation=false, face_rotation=true, face_flip=false
            face_orientations[3][i] = k + (n - 1 - j) * n;
            // face_orientation=true,  face_rotation=false, face_flip=true
            face_orientations[4][i] = (n - 1 - k) + (n - 1 - j) * n;
            // face_orientation=false, face_rotation=false, face_flip=true
            face_orientations[5][i] = (n - 1 - j) + (n - 1 - k) * n;
            // face_orientation=true,  face_rotation=true, face_flip=true
            face_orientations[6][i] = (n - 1 - j) + k * n;
            // face_orientation=false, face_rotation=true, face_flip=true
            face_orientations[7][i] = (n - 1 - k) + j * n;
          }
      return face_orientations;
    }



    template <typename Number>
    std::size_t
    ShapeInfo<Number>::memory_consumption() const
    {
      std::size_t memory = sizeof(*this);
      for (const auto &univariate_shape_data : data)
        memory += univariate_shape_data.memory_consumption();
      return memory;
    }

    template <typename Number>
    std::size_t
    UnivariateShapeData<Number>::memory_consumption() const
    {
      std::size_t memory = sizeof(*this);
      memory += MemoryConsumption::memory_consumption(shape_values);
      memory += MemoryConsumption::memory_consumption(shape_gradients);
      memory += MemoryConsumption::memory_consumption(shape_hessians);
      memory +=
        MemoryConsumption::memory_consumption(shape_gradients_collocation);
      memory +=
        MemoryConsumption::memory_consumption(shape_hessians_collocation);
      memory += MemoryConsumption::memory_consumption(shape_values_eo);
      memory += MemoryConsumption::memory_consumption(shape_gradients_eo);
      memory += MemoryConsumption::memory_consumption(shape_hessians_eo);
      memory +=
        MemoryConsumption::memory_consumption(shape_gradients_collocation_eo);
      memory +=
        MemoryConsumption::memory_consumption(shape_hessians_collocation_eo);
      for (unsigned int i = 0; i < 2; ++i)
        {
          memory +=
            MemoryConsumption::memory_consumption(shape_data_on_face[i]);
          memory +=
            MemoryConsumption::memory_consumption(quadrature_data_on_face[i]);
          memory +=
            MemoryConsumption::memory_consumption(values_within_subface[i]);
          memory +=
            MemoryConsumption::memory_consumption(gradients_within_subface[i]);
        }
      return memory;
    }

  } // namespace MatrixFreeFunctions

} // namespace internal


DEAL_II_NAMESPACE_CLOSE

#endif
