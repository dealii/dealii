// ---------------------------------------------------------------------
//
// Copyright (C) 2013 - 2020 by the deal.II authors
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

#include <deal.II/base/logstream.h>
#include <deal.II/base/qprojector.h>
#include <deal.II/base/quadrature.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/utilities.h>

#include <deal.II/dofs/dof_accessor.h>

#include <deal.II/fe/fe_nedelec.h>
#include <deal.II/fe/fe_nothing.h>
#include <deal.II/fe/fe_tools.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/vector.h>

#include <iostream>
#include <memory>
#include <sstream>

DEAL_II_NAMESPACE_OPEN

//#define DEBUG_NEDELEC

namespace internal
{
  namespace FE_Nedelec
  {
    namespace
    {
      double
      get_embedding_computation_tolerance(const unsigned int p)
      {
        // This heuristic was computed by monitoring the worst residual
        // resulting from the least squares computation when computing
        // the face embedding matrices in the FE_Nedelec constructor.
        // The residual growth is exponential, but is bounded by this
        // function up to degree 12.
        return 1.e-15 * std::exp(std::pow(p, 1.075));
      }
    } // namespace
  }   // namespace FE_Nedelec
} // namespace internal


// TODO: implement the adjust_quad_dof_index_for_face_orientation_table and
// adjust_line_dof_index_for_line_orientation_table fields, and write tests
// similar to bits/face_orientation_and_fe_q_*

template <int dim>
FE_Nedelec<dim>::FE_Nedelec(const unsigned int order)
  : FE_PolyTensor<dim>(
      PolynomialsNedelec<dim>(order),
      FiniteElementData<dim>(get_dpo_vector(order),
                             dim,
                             order + 1,
                             FiniteElementData<dim>::Hcurl),
      std::vector<bool>(PolynomialsNedelec<dim>::n_polynomials(order), true),
      std::vector<ComponentMask>(PolynomialsNedelec<dim>::n_polynomials(order),
                                 std::vector<bool>(dim, true)))
{
#ifdef DEBUG_NEDELEC
  deallog << get_name() << std::endl;
#endif

  Assert(dim >= 2, ExcImpossibleInDim(dim));

  const unsigned int n_dofs = this->n_dofs_per_cell();

  this->mapping_kind = {mapping_nedelec};
  // First, initialize the
  // generalized support points and
  // quadrature weights, since they
  // are required for interpolation.
  initialize_support_points(order);
  this->inverse_node_matrix.reinit(n_dofs, n_dofs);
  this->inverse_node_matrix.fill(FullMatrix<double>(IdentityMatrix(n_dofs)));
  // From now on, the shape functions
  // will be the correct ones, not
  // the raw shape functions anymore.

  // do not initialize embedding and restriction here. these matrices are
  // initialized on demand in get_restriction_matrix and
  // get_prolongation_matrix

#ifdef DEBUG_NEDELEC
  deallog << "Face Embedding" << std::endl;
#endif
  FullMatrix<double> face_embeddings[GeometryInfo<dim>::max_children_per_face];

  // TODO: the implementation makes the assumption that all faces have the
  // same number of dofs
  AssertDimension(this->n_unique_faces(), 1);
  const unsigned int face_no = 0;

  for (unsigned int i = 0; i < GeometryInfo<dim>::max_children_per_face; ++i)
    face_embeddings[i].reinit(this->n_dofs_per_face(face_no),
                              this->n_dofs_per_face(face_no));

  FETools::compute_face_embedding_matrices<dim, double>(
    *this,
    face_embeddings,
    0,
    0,
    internal::FE_Nedelec::get_embedding_computation_tolerance(order));

  switch (dim)
    {
      case 1:
        {
          this->interface_constraints.reinit(0, 0);
          break;
        }

      case 2:
        {
          this->interface_constraints.reinit(2 * this->n_dofs_per_face(face_no),
                                             this->n_dofs_per_face(face_no));

          for (unsigned int i = 0; i < GeometryInfo<2>::max_children_per_face;
               ++i)
            for (unsigned int j = 0; j < this->n_dofs_per_face(face_no); ++j)
              for (unsigned int k = 0; k < this->n_dofs_per_face(face_no); ++k)
                this->interface_constraints(i * this->n_dofs_per_face(face_no) +
                                              j,
                                            k) = face_embeddings[i](j, k);

          break;
        }

      case 3:
        {
          this->interface_constraints.reinit(
            4 * (this->n_dofs_per_face(face_no) - this->degree),
            this->n_dofs_per_face(face_no));

          unsigned int target_row = 0;

          for (unsigned int i = 0; i < 2; ++i)
            for (unsigned int j = this->degree; j < 2 * this->degree;
                 ++j, ++target_row)
              for (unsigned int k = 0; k < this->n_dofs_per_face(face_no); ++k)
                this->interface_constraints(target_row, k) =
                  face_embeddings[2 * i](j, k);

          for (unsigned int i = 0; i < 2; ++i)
            for (unsigned int j = 3 * this->degree;
                 j < GeometryInfo<3>::lines_per_face * this->degree;
                 ++j, ++target_row)
              for (unsigned int k = 0; k < this->n_dofs_per_face(face_no); ++k)
                this->interface_constraints(target_row, k) =
                  face_embeddings[i](j, k);

          for (unsigned int i = 0; i < 2; ++i)
            for (unsigned int j = 0; j < 2; ++j)
              for (unsigned int k = i * this->degree;
                   k < (i + 1) * this->degree;
                   ++k, ++target_row)
                for (unsigned int l = 0; l < this->n_dofs_per_face(face_no);
                     ++l)
                  this->interface_constraints(target_row, l) =
                    face_embeddings[i + 2 * j](k, l);

          for (unsigned int i = 0; i < 2; ++i)
            for (unsigned int j = 0; j < 2; ++j)
              for (unsigned int k = (i + 2) * this->degree;
                   k < (i + 3) * this->degree;
                   ++k, ++target_row)
                for (unsigned int l = 0; l < this->n_dofs_per_face(face_no);
                     ++l)
                  this->interface_constraints(target_row, l) =
                    face_embeddings[2 * i + j](k, l);

          for (unsigned int i = 0; i < GeometryInfo<3>::max_children_per_face;
               ++i)
            for (unsigned int j =
                   GeometryInfo<3>::lines_per_face * this->degree;
                 j < this->n_dofs_per_face(face_no);
                 ++j, ++target_row)
              for (unsigned int k = 0; k < this->n_dofs_per_face(face_no); ++k)
                this->interface_constraints(target_row, k) =
                  face_embeddings[i](j, k);

          break;
        }

      default:
        Assert(false, ExcNotImplemented());
    }

  // We need to initialize the dof permuation table and the one for the sign
  // change.
  initialize_quad_dof_index_permutation_and_sign_change();
}


template <int dim>
void
FE_Nedelec<dim>::initialize_quad_dof_index_permutation_and_sign_change()
{
  // for 1D and 2D, do nothing
  if (dim < 3)
    return;

  // TODO: Implement this for this class
  return;
}


template <int dim>
std::string
FE_Nedelec<dim>::get_name() const
{
  // note that the
  // FETools::get_fe_by_name
  // function depends on the
  // particular format of the string
  // this function returns, so they
  // have to be kept in synch

  std::ostringstream namebuf;
  namebuf << "FE_Nedelec<" << dim << ">(" << this->degree - 1 << ")";

  return namebuf.str();
}


template <int dim>
std::unique_ptr<FiniteElement<dim, dim>>
FE_Nedelec<dim>::clone() const
{
  return std::make_unique<FE_Nedelec<dim>>(*this);
}

//---------------------------------------------------------------------------
// Auxiliary and internal functions
//---------------------------------------------------------------------------



// Set the generalized support
// points and precompute the
// parts of the projection-based
// interpolation, which does
// not depend on the interpolated
// function.
template <>
void
FE_Nedelec<1>::initialize_support_points(const unsigned int)
{
  Assert(false, ExcNotImplemented());
}



template <>
void
FE_Nedelec<2>::initialize_support_points(const unsigned int order)
{
  const int dim = 2;

  // TODO: the implementation makes the assumption that all faces have the
  // same number of dofs
  AssertDimension(this->n_unique_faces(), 1);
  const unsigned int face_no = 0;

  // Create polynomial basis.
  const std::vector<Polynomials::Polynomial<double>> &lobatto_polynomials =
    Polynomials::Lobatto::generate_complete_basis(order + 1);
  std::vector<Polynomials::Polynomial<double>> lobatto_polynomials_grad(order +
                                                                        1);

  for (unsigned int i = 0; i < lobatto_polynomials_grad.size(); ++i)
    lobatto_polynomials_grad[i] = lobatto_polynomials[i + 1].derivative();

  // Initialize quadratures to obtain
  // quadrature points later on.
  const QGauss<dim - 1> reference_edge_quadrature(order + 1);
  const unsigned int    n_edge_points = reference_edge_quadrature.size();
  const unsigned int    n_boundary_points =
    GeometryInfo<dim>::lines_per_cell * n_edge_points;
  const Quadrature<dim> edge_quadrature =
    QProjector<dim>::project_to_all_faces(this->reference_cell(),
                                          reference_edge_quadrature);

  this->generalized_face_support_points[face_no].resize(n_edge_points);

  // Create face support points.
  for (unsigned int q_point = 0; q_point < n_edge_points; ++q_point)
    this->generalized_face_support_points[face_no][q_point] =
      reference_edge_quadrature.point(q_point);

  if (order > 0)
    {
      // If the polynomial degree is positive
      // we have support points on the faces
      // and in the interior of a cell.
      const QGauss<dim>  quadrature(order + 1);
      const unsigned int n_interior_points = quadrature.size();

      this->generalized_support_points.resize(n_boundary_points +
                                              n_interior_points);
      boundary_weights.reinit(n_edge_points, order);

      for (unsigned int q_point = 0; q_point < n_edge_points; ++q_point)
        {
          for (unsigned int line = 0; line < GeometryInfo<dim>::lines_per_cell;
               ++line)
            this->generalized_support_points[line * n_edge_points + q_point] =
              edge_quadrature.point(
                QProjector<dim>::DataSetDescriptor::face(this->reference_cell(),
                                                         line,
                                                         true,
                                                         false,
                                                         false,
                                                         n_edge_points) +
                q_point);

          for (unsigned int i = 0; i < order; ++i)
            boundary_weights(q_point, i) =
              reference_edge_quadrature.weight(q_point) *
              lobatto_polynomials_grad[i + 1].value(
                this->generalized_face_support_points[face_no][q_point](0));
        }

      for (unsigned int q_point = 0; q_point < n_interior_points; ++q_point)
        this->generalized_support_points[q_point + n_boundary_points] =
          quadrature.point(q_point);
    }

  else
    {
      // In this case we only need support points
      // on the faces of a cell.
      this->generalized_support_points.resize(n_boundary_points);

      for (unsigned int line = 0; line < GeometryInfo<dim>::lines_per_cell;
           ++line)
        for (unsigned int q_point = 0; q_point < n_edge_points; ++q_point)
          this->generalized_support_points[line * n_edge_points + q_point] =
            edge_quadrature.point(
              QProjector<dim>::DataSetDescriptor::face(this->reference_cell(),
                                                       line,
                                                       true,
                                                       false,
                                                       false,
                                                       n_edge_points) +
              q_point);
    }
}



template <>
void
FE_Nedelec<3>::initialize_support_points(const unsigned int order)
{
  const int dim = 3;

  // TODO: the implementation makes the assumption that all faces have the
  // same number of dofs
  AssertDimension(this->n_unique_faces(), 1);
  const unsigned int face_no = 0;

  // Create polynomial basis.
  const std::vector<Polynomials::Polynomial<double>> &lobatto_polynomials =
    Polynomials::Lobatto::generate_complete_basis(order + 1);
  std::vector<Polynomials::Polynomial<double>> lobatto_polynomials_grad(order +
                                                                        1);

  for (unsigned int i = 0; i < lobatto_polynomials_grad.size(); ++i)
    lobatto_polynomials_grad[i] = lobatto_polynomials[i + 1].derivative();

  // Initialize quadratures to obtain
  // quadrature points later on.
  const QGauss<1>            reference_edge_quadrature(order + 1);
  const unsigned int         n_edge_points = reference_edge_quadrature.size();
  const Quadrature<dim - 1> &edge_quadrature =
    QProjector<dim - 1>::project_to_all_faces(
      ReferenceCells::get_hypercube<dim - 1>(), reference_edge_quadrature);

  if (order > 0)
    {
      // If the polynomial order is positive
      // we have support points on the edges,
      // faces and in the interior of a cell.
      const QGauss<dim - 1> reference_face_quadrature(order + 1);
      const unsigned int    n_face_points = reference_face_quadrature.size();
      const unsigned int    n_boundary_points =
        GeometryInfo<dim>::lines_per_cell * n_edge_points +
        GeometryInfo<dim>::faces_per_cell * n_face_points;
      const QGauss<dim>  quadrature(order + 1);
      const unsigned int n_interior_points = quadrature.size();

      boundary_weights.reinit(n_edge_points + n_face_points,
                              2 * (order + 1) * order);
      this->generalized_face_support_points[face_no].resize(4 * n_edge_points +
                                                            n_face_points);
      this->generalized_support_points.resize(n_boundary_points +
                                              n_interior_points);

      // Create support points on edges.
      for (unsigned int q_point = 0; q_point < n_edge_points; ++q_point)
        {
          for (unsigned int line = 0;
               line < GeometryInfo<dim - 1>::lines_per_cell;
               ++line)
            this
              ->generalized_face_support_points[face_no][line * n_edge_points +
                                                         q_point] =
              edge_quadrature.point(
                QProjector<dim - 1>::DataSetDescriptor::face(
                  ReferenceCells::get_hypercube<dim - 1>(),
                  line,
                  true,
                  false,
                  false,
                  n_edge_points) +
                q_point);

          for (unsigned int i = 0; i < 2; ++i)
            for (unsigned int j = 0; j < 2; ++j)
              {
                this->generalized_support_points[q_point +
                                                 (i + 4 * j) * n_edge_points] =
                  Point<dim>(i, reference_edge_quadrature.point(q_point)(0), j);
                this->generalized_support_points[q_point + (i + 4 * j + 2) *
                                                             n_edge_points] =
                  Point<dim>(reference_edge_quadrature.point(q_point)(0), i, j);
                this->generalized_support_points[q_point + (i + 2 * (j + 4)) *
                                                             n_edge_points] =
                  Point<dim>(i, j, reference_edge_quadrature.point(q_point)(0));
              }

          for (unsigned int i = 0; i < order; ++i)
            boundary_weights(q_point, i) =
              reference_edge_quadrature.weight(q_point) *
              lobatto_polynomials_grad[i + 1].value(
                this->generalized_face_support_points[face_no][q_point](1));
        }

      // Create support points on faces.
      for (unsigned int q_point = 0; q_point < n_face_points; ++q_point)
        {
          this->generalized_face_support_points[face_no]
                                               [q_point + 4 * n_edge_points] =
            reference_face_quadrature.point(q_point);

          for (unsigned int i = 0; i <= order; ++i)
            for (unsigned int j = 0; j < order; ++j)
              {
                boundary_weights(q_point + n_edge_points, 2 * (i * order + j)) =
                  reference_face_quadrature.weight(q_point) *
                  lobatto_polynomials_grad[i].value(
                    this->generalized_face_support_points
                      [face_no][q_point + 4 * n_edge_points](0)) *
                  lobatto_polynomials[j + 2].value(
                    this->generalized_face_support_points
                      [face_no][q_point + 4 * n_edge_points](1));
                boundary_weights(q_point + n_edge_points,
                                 2 * (i * order + j) + 1) =
                  reference_face_quadrature.weight(q_point) *
                  lobatto_polynomials_grad[i].value(
                    this->generalized_face_support_points
                      [face_no][q_point + 4 * n_edge_points](1)) *
                  lobatto_polynomials[j + 2].value(
                    this->generalized_face_support_points
                      [face_no][q_point + 4 * n_edge_points](0));
              }
        }

      const Quadrature<dim> &face_quadrature =
        QProjector<dim>::project_to_all_faces(this->reference_cell(),
                                              reference_face_quadrature);

      for (const unsigned int face : GeometryInfo<dim>::face_indices())
        for (unsigned int q_point = 0; q_point < n_face_points; ++q_point)
          {
            this->generalized_support_points[face * n_face_points + q_point +
                                             GeometryInfo<dim>::lines_per_cell *
                                               n_edge_points] =
              face_quadrature.point(
                QProjector<dim>::DataSetDescriptor::face(this->reference_cell(),
                                                         face,
                                                         true,
                                                         false,
                                                         false,
                                                         n_face_points) +
                q_point);
          }

      // Create support points in the interior.
      for (unsigned int q_point = 0; q_point < n_interior_points; ++q_point)
        this->generalized_support_points[q_point + n_boundary_points] =
          quadrature.point(q_point);
    }

  else
    {
      this->generalized_face_support_points[face_no].resize(4 * n_edge_points);
      this->generalized_support_points.resize(
        GeometryInfo<dim>::lines_per_cell * n_edge_points);

      for (unsigned int q_point = 0; q_point < n_edge_points; ++q_point)
        {
          for (unsigned int line = 0;
               line < GeometryInfo<dim - 1>::lines_per_cell;
               ++line)
            this
              ->generalized_face_support_points[face_no][line * n_edge_points +
                                                         q_point] =
              edge_quadrature.point(
                QProjector<dim - 1>::DataSetDescriptor::face(
                  ReferenceCells::get_hypercube<dim - 1>(),
                  line,
                  true,
                  false,
                  false,
                  n_edge_points) +
                q_point);

          for (unsigned int i = 0; i < 2; ++i)
            for (unsigned int j = 0; j < 2; ++j)
              {
                this->generalized_support_points[q_point +
                                                 (i + 4 * j) * n_edge_points] =
                  Point<dim>(i, reference_edge_quadrature.point(q_point)(0), j);
                this->generalized_support_points[q_point + (i + 4 * j + 2) *
                                                             n_edge_points] =
                  Point<dim>(reference_edge_quadrature.point(q_point)(0), i, j);
                this->generalized_support_points[q_point + (i + 2 * (j + 4)) *
                                                             n_edge_points] =
                  Point<dim>(i, j, reference_edge_quadrature.point(q_point)(0));
              }
        }
    }
}



// Set the restriction matrices.
template <>
void
FE_Nedelec<1>::initialize_restriction()
{
  // there is only one refinement case in 1d,
  // which is the isotropic one
  for (unsigned int i = 0; i < GeometryInfo<1>::max_children_per_cell; ++i)
    this->restriction[0][i].reinit(0, 0);
}



// Restriction operator
template <int dim>
void
FE_Nedelec<dim>::initialize_restriction()
{
  // This function does the same as the
  // function interpolate further below.
  // But since the functions, which we
  // interpolate here, are discontinuous
  // we have to use more quadrature
  // points as in interpolate.
  const QGauss<1>              edge_quadrature(2 * this->degree);
  const std::vector<Point<1>> &edge_quadrature_points =
    edge_quadrature.get_points();
  const unsigned int n_edge_quadrature_points = edge_quadrature.size();
  const unsigned int index = RefinementCase<dim>::isotropic_refinement - 1;

  switch (dim)
    {
      case 2:
        {
          // First interpolate the shape
          // functions of the child cells
          // to the lowest order shape
          // functions of the parent cell.
          for (unsigned int dof = 0; dof < this->n_dofs_per_cell(); ++dof)
            for (unsigned int q_point = 0; q_point < n_edge_quadrature_points;
                 ++q_point)
              {
                const double weight = 2.0 * edge_quadrature.weight(q_point);

                if (edge_quadrature_points[q_point](0) < 0.5)
                  {
                    Point<dim> quadrature_point(
                      0.0, 2.0 * edge_quadrature_points[q_point](0));

                    this->restriction[index][0](0, dof) +=
                      weight *
                      this->shape_value_component(dof, quadrature_point, 1);
                    quadrature_point(0) = 1.0;
                    this->restriction[index][1](this->degree, dof) +=
                      weight *
                      this->shape_value_component(dof, quadrature_point, 1);
                    quadrature_point(0) = quadrature_point(1);
                    quadrature_point(1) = 0.0;
                    this->restriction[index][0](2 * this->degree, dof) +=
                      weight *
                      this->shape_value_component(dof, quadrature_point, 0);
                    quadrature_point(1) = 1.0;
                    this->restriction[index][2](3 * this->degree, dof) +=
                      weight *
                      this->shape_value_component(dof, quadrature_point, 0);
                  }

                else
                  {
                    Point<dim> quadrature_point(
                      0.0, 2.0 * edge_quadrature_points[q_point](0) - 1.0);

                    this->restriction[index][2](0, dof) +=
                      weight *
                      this->shape_value_component(dof, quadrature_point, 1);
                    quadrature_point(0) = 1.0;
                    this->restriction[index][3](this->degree, dof) +=
                      weight *
                      this->shape_value_component(dof, quadrature_point, 1);
                    quadrature_point(0) = quadrature_point(1);
                    quadrature_point(1) = 0.0;
                    this->restriction[index][1](2 * this->degree, dof) +=
                      weight *
                      this->shape_value_component(dof, quadrature_point, 0);
                    quadrature_point(1) = 1.0;
                    this->restriction[index][3](3 * this->degree, dof) +=
                      weight *
                      this->shape_value_component(dof, quadrature_point, 0);
                  }
              }

          // Then project the shape functions
          // of the child cells to the higher
          // order shape functions of the
          // parent cell.
          if (this->degree > 1)
            {
              const unsigned int deg = this->degree - 1;
              const std::vector<Polynomials::Polynomial<double>>
                &legendre_polynomials =
                  Polynomials::Legendre::generate_complete_basis(deg);
              FullMatrix<double> system_matrix_inv(deg, deg);

              {
                FullMatrix<double> assembling_matrix(deg,
                                                     n_edge_quadrature_points);

                for (unsigned int q_point = 0;
                     q_point < n_edge_quadrature_points;
                     ++q_point)
                  {
                    const double weight =
                      std::sqrt(edge_quadrature.weight(q_point));

                    for (unsigned int i = 0; i < deg; ++i)
                      assembling_matrix(i, q_point) =
                        weight * legendre_polynomials[i + 1].value(
                                   edge_quadrature_points[q_point](0));
                  }

                FullMatrix<double> system_matrix(deg, deg);

                assembling_matrix.mTmult(system_matrix, assembling_matrix);
                system_matrix_inv.invert(system_matrix);
              }

              FullMatrix<double> solution(this->degree - 1, 4);
              FullMatrix<double> system_rhs(this->degree - 1, 4);
              Vector<double>     tmp(4);

              for (unsigned int dof = 0; dof < this->n_dofs_per_cell(); ++dof)
                for (unsigned int i = 0; i < 2; ++i)
                  {
                    system_rhs = 0.0;

                    for (unsigned int q_point = 0;
                         q_point < n_edge_quadrature_points;
                         ++q_point)
                      {
                        const double weight = edge_quadrature.weight(q_point);
                        const Point<dim> quadrature_point_0(
                          i, edge_quadrature_points[q_point](0));
                        const Point<dim> quadrature_point_1(
                          edge_quadrature_points[q_point](0), i);

                        if (edge_quadrature_points[q_point](0) < 0.5)
                          {
                            Point<dim> quadrature_point_2(
                              i, 2.0 * edge_quadrature_points[q_point](0));

                            tmp(0) =
                              weight *
                              (2.0 * this->shape_value_component(
                                       dof, quadrature_point_2, 1) -
                               this->restriction[index][i](i * this->degree,
                                                           dof) *
                                 this->shape_value_component(i * this->degree,
                                                             quadrature_point_0,
                                                             1));
                            tmp(1) =
                              -1.0 * weight *
                              this->restriction[index][i + 2](i * this->degree,
                                                              dof) *
                              this->shape_value_component(i * this->degree,
                                                          quadrature_point_0,
                                                          1);
                            quadrature_point_2 = Point<dim>(
                              2.0 * edge_quadrature_points[q_point](0), i);
                            tmp(2) =
                              weight *
                              (2.0 * this->shape_value_component(
                                       dof, quadrature_point_2, 0) -
                               this->restriction[index][2 * i]((i + 2) *
                                                                 this->degree,
                                                               dof) *
                                 this->shape_value_component((i + 2) *
                                                               this->degree,
                                                             quadrature_point_1,
                                                             0));
                            tmp(3) =
                              -1.0 * weight *
                              this->restriction[index][2 * i + 1](
                                (i + 2) * this->degree, dof) *
                              this->shape_value_component(
                                (i + 2) * this->degree, quadrature_point_1, 0);
                          }

                        else
                          {
                            tmp(0) =
                              -1.0 * weight *
                              this->restriction[index][i](i * this->degree,
                                                          dof) *
                              this->shape_value_component(i * this->degree,
                                                          quadrature_point_0,
                                                          1);

                            Point<dim> quadrature_point_2(
                              i,
                              2.0 * edge_quadrature_points[q_point](0) - 1.0);

                            tmp(1) =
                              weight *
                              (2.0 * this->shape_value_component(
                                       dof, quadrature_point_2, 1) -
                               this->restriction[index][i + 2](i * this->degree,
                                                               dof) *
                                 this->shape_value_component(i * this->degree,
                                                             quadrature_point_0,
                                                             1));
                            tmp(2) =
                              -1.0 * weight *
                              this->restriction[index][2 * i]((i + 2) *
                                                                this->degree,
                                                              dof) *
                              this->shape_value_component(
                                (i + 2) * this->degree, quadrature_point_1, 0);
                            quadrature_point_2 = Point<dim>(
                              2.0 * edge_quadrature_points[q_point](0) - 1.0,
                              i);
                            tmp(3) =
                              weight *
                              (2.0 * this->shape_value_component(
                                       dof, quadrature_point_2, 0) -
                               this->restriction[index][2 * i + 1](
                                 (i + 2) * this->degree, dof) *
                                 this->shape_value_component((i + 2) *
                                                               this->degree,
                                                             quadrature_point_1,
                                                             0));
                          }

                        for (unsigned int j = 0; j < this->degree - 1; ++j)
                          {
                            const double L_j =
                              legendre_polynomials[j + 1].value(
                                edge_quadrature_points[q_point](0));

                            for (unsigned int k = 0; k < tmp.size(); ++k)
                              system_rhs(j, k) += tmp(k) * L_j;
                          }
                      }

                    system_matrix_inv.mmult(solution, system_rhs);

                    for (unsigned int j = 0; j < this->degree - 1; ++j)
                      for (unsigned int k = 0; k < 2; ++k)
                        {
                          if (std::abs(solution(j, k)) > 1e-14)
                            this->restriction[index][i + 2 * k](
                              i * this->degree + j + 1, dof) = solution(j, k);

                          if (std::abs(solution(j, k + 2)) > 1e-14)
                            this->restriction[index][2 * i + k](
                              (i + 2) * this->degree + j + 1, dof) =
                              solution(j, k + 2);
                        }
                  }

              const QGauss<dim>              quadrature(2 * this->degree);
              const std::vector<Point<dim>> &quadrature_points =
                quadrature.get_points();
              const std::vector<Polynomials::Polynomial<double>>
                &lobatto_polynomials =
                  Polynomials::Lobatto::generate_complete_basis(this->degree);
              const unsigned int n_boundary_dofs =
                GeometryInfo<dim>::faces_per_cell * this->degree;
              const unsigned int n_quadrature_points = quadrature.size();

              {
                FullMatrix<double> assembling_matrix((this->degree - 1) *
                                                       this->degree,
                                                     n_quadrature_points);

                for (unsigned int q_point = 0; q_point < n_quadrature_points;
                     ++q_point)
                  {
                    const double weight = std::sqrt(quadrature.weight(q_point));

                    for (unsigned int i = 0; i < this->degree; ++i)
                      {
                        const double L_i =
                          weight * legendre_polynomials[i].value(
                                     quadrature_points[q_point](0));

                        for (unsigned int j = 0; j < this->degree - 1; ++j)
                          assembling_matrix(i * (this->degree - 1) + j,
                                            q_point) =
                            L_i * lobatto_polynomials[j + 2].value(
                                    quadrature_points[q_point](1));
                      }
                  }

                FullMatrix<double> system_matrix(assembling_matrix.m(),
                                                 assembling_matrix.m());

                assembling_matrix.mTmult(system_matrix, assembling_matrix);
                system_matrix_inv.reinit(system_matrix.m(), system_matrix.m());
                system_matrix_inv.invert(system_matrix);
              }

              solution.reinit(system_matrix_inv.m(), 8);
              system_rhs.reinit(system_matrix_inv.m(), 8);
              tmp.reinit(8);

              for (unsigned int dof = 0; dof < this->n_dofs_per_cell(); ++dof)
                {
                  system_rhs = 0.0;

                  for (unsigned int q_point = 0; q_point < n_quadrature_points;
                       ++q_point)
                    {
                      tmp = 0.0;

                      if (quadrature_points[q_point](0) < 0.5)
                        {
                          if (quadrature_points[q_point](1) < 0.5)
                            {
                              const Point<dim> quadrature_point(
                                2.0 * quadrature_points[q_point](0),
                                2.0 * quadrature_points[q_point](1));

                              tmp(0) += 2.0 * this->shape_value_component(
                                                dof, quadrature_point, 0);
                              tmp(1) += 2.0 * this->shape_value_component(
                                                dof, quadrature_point, 1);
                            }

                          else
                            {
                              const Point<dim> quadrature_point(
                                2.0 * quadrature_points[q_point](0),
                                2.0 * quadrature_points[q_point](1) - 1.0);

                              tmp(4) += 2.0 * this->shape_value_component(
                                                dof, quadrature_point, 0);
                              tmp(5) += 2.0 * this->shape_value_component(
                                                dof, quadrature_point, 1);
                            }
                        }

                      else if (quadrature_points[q_point](1) < 0.5)
                        {
                          const Point<dim> quadrature_point(
                            2.0 * quadrature_points[q_point](0) - 1.0,
                            2.0 * quadrature_points[q_point](1));

                          tmp(2) +=
                            2.0 * this->shape_value_component(dof,
                                                              quadrature_point,
                                                              0);
                          tmp(3) +=
                            2.0 * this->shape_value_component(dof,
                                                              quadrature_point,
                                                              1);
                        }

                      else
                        {
                          const Point<dim> quadrature_point(
                            2.0 * quadrature_points[q_point](0) - 1.0,
                            2.0 * quadrature_points[q_point](1) - 1.0);

                          tmp(6) +=
                            2.0 * this->shape_value_component(dof,
                                                              quadrature_point,
                                                              0);
                          tmp(7) +=
                            2.0 * this->shape_value_component(dof,
                                                              quadrature_point,
                                                              1);
                        }

                      for (unsigned int i = 0; i < 2; ++i)
                        for (unsigned int j = 0; j < this->degree; ++j)
                          {
                            tmp(2 * i) -=
                              this->restriction[index][i](j + 2 * this->degree,
                                                          dof) *
                              this->shape_value_component(
                                j + 2 * this->degree,
                                quadrature_points[q_point],
                                0);
                            tmp(2 * i + 1) -=
                              this->restriction[index][i](i * this->degree + j,
                                                          dof) *
                              this->shape_value_component(
                                i * this->degree + j,
                                quadrature_points[q_point],
                                1);
                            tmp(2 * (i + 2)) -= this->restriction[index][i + 2](
                                                  j + 3 * this->degree, dof) *
                                                this->shape_value_component(
                                                  j + 3 * this->degree,
                                                  quadrature_points[q_point],
                                                  0);
                            tmp(2 * i + 5) -= this->restriction[index][i + 2](
                                                i * this->degree + j, dof) *
                                              this->shape_value_component(
                                                i * this->degree + j,
                                                quadrature_points[q_point],
                                                1);
                          }

                      tmp *= quadrature.weight(q_point);

                      for (unsigned int i = 0; i < this->degree; ++i)
                        {
                          const double L_i_0 = legendre_polynomials[i].value(
                            quadrature_points[q_point](0));
                          const double L_i_1 = legendre_polynomials[i].value(
                            quadrature_points[q_point](1));

                          for (unsigned int j = 0; j < this->degree - 1; ++j)
                            {
                              const double l_j_0 =
                                L_i_0 * lobatto_polynomials[j + 2].value(
                                          quadrature_points[q_point](1));
                              const double l_j_1 =
                                L_i_1 * lobatto_polynomials[j + 2].value(
                                          quadrature_points[q_point](0));

                              for (unsigned int k = 0; k < 4; ++k)
                                {
                                  system_rhs(i * (this->degree - 1) + j,
                                             2 * k) += tmp(2 * k) * l_j_0;
                                  system_rhs(i * (this->degree - 1) + j,
                                             2 * k + 1) +=
                                    tmp(2 * k + 1) * l_j_1;
                                }
                            }
                        }
                    }

                  system_matrix_inv.mmult(solution, system_rhs);

                  for (unsigned int i = 0; i < this->degree; ++i)
                    for (unsigned int j = 0; j < this->degree - 1; ++j)
                      for (unsigned int k = 0; k < 4; ++k)
                        {
                          if (std::abs(solution(i * (this->degree - 1) + j,
                                                2 * k)) > 1e-14)
                            this->restriction[index][k](i * (this->degree - 1) +
                                                          j + n_boundary_dofs,
                                                        dof) =
                              solution(i * (this->degree - 1) + j, 2 * k);

                          if (std::abs(solution(i * (this->degree - 1) + j,
                                                2 * k + 1)) > 1e-14)
                            this->restriction[index][k](
                              i + (this->degree - 1 + j) * this->degree +
                                n_boundary_dofs,
                              dof) =
                              solution(i * (this->degree - 1) + j, 2 * k + 1);
                        }
                }
            }

          break;
        }

      case 3:
        {
          // First interpolate the shape
          // functions of the child cells
          // to the lowest order shape
          // functions of the parent cell.
          for (unsigned int dof = 0; dof < this->n_dofs_per_cell(); ++dof)
            for (unsigned int q_point = 0; q_point < n_edge_quadrature_points;
                 ++q_point)
              {
                const double weight = 2.0 * edge_quadrature.weight(q_point);

                if (edge_quadrature_points[q_point](0) < 0.5)
                  for (unsigned int i = 0; i < 2; ++i)
                    for (unsigned int j = 0; j < 2; ++j)
                      {
                        Point<dim> quadrature_point(
                          i, 2.0 * edge_quadrature_points[q_point](0), j);

                        this->restriction[index][i + 4 * j]((i + 4 * j) *
                                                              this->degree,
                                                            dof) +=
                          weight *
                          this->shape_value_component(dof, quadrature_point, 1);
                        quadrature_point =
                          Point<dim>(2.0 * edge_quadrature_points[q_point](0),
                                     i,
                                     j);
                        this->restriction[index][2 * (i + 2 * j)](
                          (i + 4 * j + 2) * this->degree, dof) +=
                          weight *
                          this->shape_value_component(dof, quadrature_point, 0);
                        quadrature_point =
                          Point<dim>(i,
                                     j,
                                     2.0 * edge_quadrature_points[q_point](0));
                        this->restriction[index][i + 2 * j]((i + 2 * (j + 4)) *
                                                              this->degree,
                                                            dof) +=
                          weight *
                          this->shape_value_component(dof, quadrature_point, 2);
                      }

                else
                  for (unsigned int i = 0; i < 2; ++i)
                    for (unsigned int j = 0; j < 2; ++j)
                      {
                        Point<dim> quadrature_point(
                          i, 2.0 * edge_quadrature_points[q_point](0) - 1.0, j);

                        this->restriction[index][i + 4 * j + 2]((i + 4 * j) *
                                                                  this->degree,
                                                                dof) +=
                          weight *
                          this->shape_value_component(dof, quadrature_point, 1);
                        quadrature_point = Point<dim>(
                          2.0 * edge_quadrature_points[q_point](0) - 1.0, i, j);
                        this->restriction[index][2 * (i + 2 * j) + 1](
                          (i + 4 * j + 2) * this->degree, dof) +=
                          weight *
                          this->shape_value_component(dof, quadrature_point, 0);
                        quadrature_point = Point<dim>(
                          i, j, 2.0 * edge_quadrature_points[q_point](0) - 1.0);
                        this->restriction[index][i + 2 * (j + 2)](
                          (i + 2 * (j + 4)) * this->degree, dof) +=
                          weight *
                          this->shape_value_component(dof, quadrature_point, 2);
                      }
              }

          // Then project the shape functions
          // of the child cells to the higher
          // order shape functions of the
          // parent cell.
          if (this->degree > 1)
            {
              const unsigned int deg = this->degree - 1;
              const std::vector<Polynomials::Polynomial<double>>
                &legendre_polynomials =
                  Polynomials::Legendre::generate_complete_basis(deg);
              FullMatrix<double> system_matrix_inv(deg, deg);

              {
                FullMatrix<double> assembling_matrix(deg,
                                                     n_edge_quadrature_points);

                for (unsigned int q_point = 0;
                     q_point < n_edge_quadrature_points;
                     ++q_point)
                  {
                    const double weight =
                      std::sqrt(edge_quadrature.weight(q_point));

                    for (unsigned int i = 0; i < deg; ++i)
                      assembling_matrix(i, q_point) =
                        weight * legendre_polynomials[i + 1].value(
                                   edge_quadrature_points[q_point](0));
                  }

                FullMatrix<double> system_matrix(deg, deg);

                assembling_matrix.mTmult(system_matrix, assembling_matrix);
                system_matrix_inv.invert(system_matrix);
              }

              FullMatrix<double> solution(deg, 6);
              FullMatrix<double> system_rhs(deg, 6);
              Vector<double>     tmp(6);

              for (unsigned int i = 0; i < 2; ++i)
                for (unsigned int j = 0; j < 2; ++j)
                  for (unsigned int dof = 0; dof < this->n_dofs_per_cell();
                       ++dof)
                    {
                      system_rhs = 0.0;

                      for (unsigned int q_point = 0;
                           q_point < n_edge_quadrature_points;
                           ++q_point)
                        {
                          const double weight = edge_quadrature.weight(q_point);
                          const Point<dim> quadrature_point_0(
                            i, edge_quadrature_points[q_point](0), j);
                          const Point<dim> quadrature_point_1(
                            edge_quadrature_points[q_point](0), i, j);
                          const Point<dim> quadrature_point_2(
                            i, j, edge_quadrature_points[q_point](0));

                          if (edge_quadrature_points[q_point](0) < 0.5)
                            {
                              Point<dim> quadrature_point_3(
                                i, 2.0 * edge_quadrature_points[q_point](0), j);

                              tmp(0) =
                                weight * (2.0 * this->shape_value_component(
                                                  dof, quadrature_point_3, 1) -
                                          this->restriction[index][i + 4 * j](
                                            (i + 4 * j) * this->degree, dof) *
                                            this->shape_value_component(
                                              (i + 4 * j) * this->degree,
                                              quadrature_point_0,
                                              1));
                              tmp(1) =
                                -1.0 * weight *
                                this->restriction[index][i + 4 * j + 2](
                                  (i + 4 * j) * this->degree, dof) *
                                this->shape_value_component((i + 4 * j) *
                                                              this->degree,
                                                            quadrature_point_0,
                                                            1);
                              quadrature_point_3 = Point<dim>(
                                2.0 * edge_quadrature_points[q_point](0), i, j);
                              tmp(2) =
                                weight *
                                (2.0 * this->shape_value_component(
                                         dof, quadrature_point_3, 0) -
                                 this->restriction[index][2 * (i + 2 * j)](
                                   (i + 4 * j + 2) * this->degree, dof) *
                                   this->shape_value_component(
                                     (i + 4 * j + 2) * this->degree,
                                     quadrature_point_1,
                                     0));
                              tmp(3) =
                                -1.0 * weight *
                                this->restriction[index][2 * (i + 2 * j) + 1](
                                  (i + 4 * j + 2) * this->degree, dof) *
                                this->shape_value_component((i + 4 * j + 2) *
                                                              this->degree,
                                                            quadrature_point_1,
                                                            0);
                              quadrature_point_3 = Point<dim>(
                                i, j, 2.0 * edge_quadrature_points[q_point](0));
                              tmp(4) =
                                weight *
                                (2.0 * this->shape_value_component(
                                         dof, quadrature_point_3, 2) -
                                 this->restriction[index][i + 2 * j](
                                   (i + 2 * (j + 4)) * this->degree, dof) *
                                   this->shape_value_component(
                                     (i + 2 * (j + 4)) * this->degree,
                                     quadrature_point_2,
                                     2));
                              tmp(5) =
                                -1.0 * weight *
                                this->restriction[index][i + 2 * (j + 2)](
                                  (i + 2 * (j + 4)) * this->degree, dof) *
                                this->shape_value_component((i + 2 * (j + 4)) *
                                                              this->degree,
                                                            quadrature_point_2,
                                                            2);
                            }

                          else
                            {
                              tmp(0) =
                                -1.0 * weight *
                                this->restriction[index][i + 4 * j](
                                  (i + 4 * j) * this->degree, dof) *
                                this->shape_value_component((i + 4 * j) *
                                                              this->degree,
                                                            quadrature_point_0,
                                                            1);

                              Point<dim> quadrature_point_3(
                                i,
                                2.0 * edge_quadrature_points[q_point](0) - 1.0,
                                j);

                              tmp(1) = weight *
                                       (2.0 * this->shape_value_component(
                                                dof, quadrature_point_3, 1) -
                                        this->restriction[index][i + 4 * j + 2](
                                          (i + 4 * j) * this->degree, dof) *
                                          this->shape_value_component(
                                            (i + 4 * j) * this->degree,
                                            quadrature_point_0,
                                            1));
                              tmp(2) =
                                -1.0 * weight *
                                this->restriction[index][2 * (i + 2 * j)](
                                  (i + 4 * j + 2) * this->degree, dof) *
                                this->shape_value_component((i + 4 * j + 2) *
                                                              this->degree,
                                                            quadrature_point_1,
                                                            0);
                              quadrature_point_3 = Point<dim>(
                                2.0 * edge_quadrature_points[q_point](0) - 1.0,
                                i,
                                j);
                              tmp(3) =
                                weight *
                                (2.0 * this->shape_value_component(
                                         dof, quadrature_point_3, 0) -
                                 this->restriction[index][2 * (i + 2 * j) + 1](
                                   (i + 4 * j + 2) * this->degree, dof) *
                                   this->shape_value_component(
                                     (i + 4 * j + 2) * this->degree,
                                     quadrature_point_1,
                                     0));
                              tmp(4) =
                                -1.0 * weight *
                                this->restriction[index][i + 2 * j](
                                  (i + 2 * (j + 4)) * this->degree, dof) *
                                this->shape_value_component((i + 2 * (j + 4)) *
                                                              this->degree,
                                                            quadrature_point_2,
                                                            2);
                              quadrature_point_3 = Point<dim>(
                                i,
                                j,
                                2.0 * edge_quadrature_points[q_point](0) - 1.0);
                              tmp(5) =
                                weight *
                                (2.0 * this->shape_value_component(
                                         dof, quadrature_point_3, 2) -
                                 this->restriction[index][i + 2 * (j + 2)](
                                   (i + 2 * (j + 4)) * this->degree, dof) *
                                   this->shape_value_component(
                                     (i + 2 * (j + 4)) * this->degree,
                                     quadrature_point_2,
                                     2));
                            }

                          for (unsigned int k = 0; k < deg; ++k)
                            {
                              const double L_k =
                                legendre_polynomials[k + 1].value(
                                  edge_quadrature_points[q_point](0));

                              for (unsigned int l = 0; l < tmp.size(); ++l)
                                system_rhs(k, l) += tmp(l) * L_k;
                            }
                        }

                      system_matrix_inv.mmult(solution, system_rhs);

                      for (unsigned int k = 0; k < 2; ++k)
                        for (unsigned int l = 0; l < deg; ++l)
                          {
                            if (std::abs(solution(l, k)) > 1e-14)
                              this->restriction[index][i + 2 * (2 * j + k)](
                                (i + 4 * j) * this->degree + l + 1, dof) =
                                solution(l, k);

                            if (std::abs(solution(l, k + 2)) > 1e-14)
                              this->restriction[index][2 * (i + 2 * j) + k](
                                (i + 4 * j + 2) * this->degree + l + 1, dof) =
                                solution(l, k + 2);

                            if (std::abs(solution(l, k + 4)) > 1e-14)
                              this->restriction[index][i + 2 * (j + 2 * k)](
                                (i + 2 * (j + 4)) * this->degree + l + 1, dof) =
                                solution(l, k + 4);
                          }
                    }

              const QGauss<2>              face_quadrature(2 * this->degree);
              const std::vector<Point<2>> &face_quadrature_points =
                face_quadrature.get_points();
              const std::vector<Polynomials::Polynomial<double>>
                &lobatto_polynomials =
                  Polynomials::Lobatto::generate_complete_basis(this->degree);
              const unsigned int n_edge_dofs =
                GeometryInfo<dim>::lines_per_cell * this->degree;
              const unsigned int n_face_quadrature_points =
                face_quadrature.size();

              {
                FullMatrix<double> assembling_matrix(deg * this->degree,
                                                     n_face_quadrature_points);

                for (unsigned int q_point = 0;
                     q_point < n_face_quadrature_points;
                     ++q_point)
                  {
                    const double weight =
                      std::sqrt(face_quadrature.weight(q_point));

                    for (unsigned int i = 0; i <= deg; ++i)
                      {
                        const double L_i =
                          weight * legendre_polynomials[i].value(
                                     face_quadrature_points[q_point](0));

                        for (unsigned int j = 0; j < deg; ++j)
                          assembling_matrix(i * deg + j, q_point) =
                            L_i * lobatto_polynomials[j + 2].value(
                                    face_quadrature_points[q_point](1));
                      }
                  }

                FullMatrix<double> system_matrix(assembling_matrix.m(),
                                                 assembling_matrix.m());

                assembling_matrix.mTmult(system_matrix, assembling_matrix);
                system_matrix_inv.reinit(system_matrix.m(), system_matrix.m());
                system_matrix_inv.invert(system_matrix);
              }

              solution.reinit(system_matrix_inv.m(), 24);
              system_rhs.reinit(system_matrix_inv.m(), 24);
              tmp.reinit(24);

              for (unsigned int i = 0; i < 2; ++i)
                for (unsigned int dof = 0; dof < this->n_dofs_per_cell(); ++dof)
                  {
                    system_rhs = 0.0;

                    for (unsigned int q_point = 0;
                         q_point < n_face_quadrature_points;
                         ++q_point)
                      {
                        tmp = 0.0;

                        if (face_quadrature_points[q_point](0) < 0.5)
                          {
                            if (face_quadrature_points[q_point](1) < 0.5)
                              {
                                Point<dim> quadrature_point_0(
                                  i,
                                  2.0 * face_quadrature_points[q_point](0),
                                  2.0 * face_quadrature_points[q_point](1));

                                tmp(0) += 2.0 * this->shape_value_component(
                                                  dof, quadrature_point_0, 1);
                                tmp(1) += 2.0 * this->shape_value_component(
                                                  dof, quadrature_point_0, 2);
                                quadrature_point_0 = Point<dim>(
                                  2.0 * face_quadrature_points[q_point](0),
                                  i,
                                  2.0 * face_quadrature_points[q_point](1));
                                tmp(8) += 2.0 * this->shape_value_component(
                                                  dof, quadrature_point_0, 2);
                                tmp(9) += 2.0 * this->shape_value_component(
                                                  dof, quadrature_point_0, 0);
                                quadrature_point_0 = Point<dim>(
                                  2.0 * face_quadrature_points[q_point](0),
                                  2.0 * face_quadrature_points[q_point](1),
                                  i);
                                tmp(16) += 2.0 * this->shape_value_component(
                                                   dof, quadrature_point_0, 0);
                                tmp(17) += 2.0 * this->shape_value_component(
                                                   dof, quadrature_point_0, 1);
                              }

                            else
                              {
                                Point<dim> quadrature_point_0(
                                  i,
                                  2.0 * face_quadrature_points[q_point](0),
                                  2.0 * face_quadrature_points[q_point](1) -
                                    1.0);

                                tmp(2) += 2.0 * this->shape_value_component(
                                                  dof, quadrature_point_0, 1);
                                tmp(3) += 2.0 * this->shape_value_component(
                                                  dof, quadrature_point_0, 2);
                                quadrature_point_0 = Point<dim>(
                                  2.0 * face_quadrature_points[q_point](0),
                                  i,
                                  2.0 * face_quadrature_points[q_point](1) -
                                    1.0);
                                tmp(10) += 2.0 * this->shape_value_component(
                                                   dof, quadrature_point_0, 2);
                                tmp(11) += 2.0 * this->shape_value_component(
                                                   dof, quadrature_point_0, 0);
                                quadrature_point_0 = Point<dim>(
                                  2.0 * face_quadrature_points[q_point](0),
                                  2.0 * face_quadrature_points[q_point](1) -
                                    1.0,
                                  i);
                                tmp(18) += 2.0 * this->shape_value_component(
                                                   dof, quadrature_point_0, 0);
                                tmp(19) += 2.0 * this->shape_value_component(
                                                   dof, quadrature_point_0, 1);
                              }
                          }

                        else if (face_quadrature_points[q_point](1) < 0.5)
                          {
                            Point<dim> quadrature_point_0(
                              i,
                              2.0 * face_quadrature_points[q_point](0) - 1.0,
                              2.0 * face_quadrature_points[q_point](1));

                            tmp(4) += 2.0 * this->shape_value_component(
                                              dof, quadrature_point_0, 1);
                            tmp(5) += 2.0 * this->shape_value_component(
                                              dof, quadrature_point_0, 2);
                            quadrature_point_0 = Point<dim>(
                              2.0 * face_quadrature_points[q_point](0) - 1.0,
                              i,
                              2.0 * face_quadrature_points[q_point](1));
                            tmp(12) += 2.0 * this->shape_value_component(
                                               dof, quadrature_point_0, 2);
                            tmp(13) += 2.0 * this->shape_value_component(
                                               dof, quadrature_point_0, 0);
                            quadrature_point_0 = Point<dim>(
                              2.0 * face_quadrature_points[q_point](0) - 1.0,
                              2.0 * face_quadrature_points[q_point](1),
                              i);
                            tmp(20) += 2.0 * this->shape_value_component(
                                               dof, quadrature_point_0, 0);
                            tmp(21) += 2.0 * this->shape_value_component(
                                               dof, quadrature_point_0, 1);
                          }

                        else
                          {
                            Point<dim> quadrature_point_0(
                              i,
                              2.0 * face_quadrature_points[q_point](0) - 1.0,
                              2.0 * face_quadrature_points[q_point](1) - 1.0);

                            tmp(6) += 2.0 * this->shape_value_component(
                                              dof, quadrature_point_0, 1);
                            tmp(7) += 2.0 * this->shape_value_component(
                                              dof, quadrature_point_0, 2);
                            quadrature_point_0 = Point<dim>(
                              2.0 * face_quadrature_points[q_point](0) - 1.0,
                              i,
                              2.0 * face_quadrature_points[q_point](1) - 1.0);
                            tmp(14) += 2.0 * this->shape_value_component(
                                               dof, quadrature_point_0, 2);
                            tmp(15) += 2.0 * this->shape_value_component(
                                               dof, quadrature_point_0, 0);
                            quadrature_point_0 = Point<dim>(
                              2.0 * face_quadrature_points[q_point](0) - 1.0,
                              2.0 * face_quadrature_points[q_point](1) - 1.0,
                              i);
                            tmp(22) += 2.0 * this->shape_value_component(
                                               dof, quadrature_point_0, 0);
                            tmp(23) += 2.0 * this->shape_value_component(
                                               dof, quadrature_point_0, 1);
                          }

                        const Point<dim> quadrature_point_0(
                          i,
                          face_quadrature_points[q_point](0),
                          face_quadrature_points[q_point](1));
                        const Point<dim> quadrature_point_1(
                          face_quadrature_points[q_point](0),
                          i,
                          face_quadrature_points[q_point](1));
                        const Point<dim> quadrature_point_2(
                          face_quadrature_points[q_point](0),
                          face_quadrature_points[q_point](1),
                          i);

                        for (unsigned int j = 0; j < 2; ++j)
                          for (unsigned int k = 0; k < 2; ++k)
                            for (unsigned int l = 0; l <= deg; ++l)
                              {
                                tmp(2 * (j + 2 * k)) -=
                                  this->restriction[index][i + 2 * (2 * j + k)](
                                    (i + 4 * j) * this->degree + l, dof) *
                                  this->shape_value_component(
                                    (i + 4 * j) * this->degree + l,
                                    quadrature_point_0,
                                    1);
                                tmp(2 * (j + 2 * k) + 1) -=
                                  this->restriction[index][i + 2 * (2 * j + k)](
                                    (i + 2 * (k + 4)) * this->degree + l, dof) *
                                  this->shape_value_component(
                                    (i + 2 * (k + 4)) * this->degree + l,
                                    quadrature_point_0,
                                    2);
                                tmp(2 * (j + 2 * (k + 2))) -=
                                  this->restriction[index][2 * (i + 2 * j) + k](
                                    (2 * (i + 4) + k) * this->degree + l, dof) *
                                  this->shape_value_component(
                                    (2 * (i + 4) + k) * this->degree + l,
                                    quadrature_point_1,
                                    2);
                                tmp(2 * (j + 2 * k) + 9) -=
                                  this->restriction[index][2 * (i + 2 * j) + k](
                                    (i + 4 * j + 2) * this->degree + l, dof) *
                                  this->shape_value_component(
                                    (i + 4 * j + 2) * this->degree + l,
                                    quadrature_point_1,
                                    0);
                                tmp(2 * (j + 2 * (k + 4))) -=
                                  this->restriction[index][2 * (2 * i + j) + k](
                                    (4 * i + j + 2) * this->degree + l, dof) *
                                  this->shape_value_component(
                                    (4 * i + j + 2) * this->degree + l,
                                    quadrature_point_2,
                                    0);
                                tmp(2 * (j + 2 * k) + 17) -=
                                  this->restriction[index][2 * (2 * i + j) + k](
                                    (4 * i + k) * this->degree + l, dof) *
                                  this->shape_value_component(
                                    (4 * i + k) * this->degree + l,
                                    quadrature_point_2,
                                    1);
                              }

                        tmp *= face_quadrature.weight(q_point);

                        for (unsigned int j = 0; j <= deg; ++j)
                          {
                            const double L_j_0 = legendre_polynomials[j].value(
                              face_quadrature_points[q_point](0));
                            const double L_j_1 = legendre_polynomials[j].value(
                              face_quadrature_points[q_point](1));

                            for (unsigned int k = 0; k < deg; ++k)
                              {
                                const double l_k_0 =
                                  L_j_0 * lobatto_polynomials[k + 2].value(
                                            face_quadrature_points[q_point](1));
                                const double l_k_1 =
                                  L_j_1 * lobatto_polynomials[k + 2].value(
                                            face_quadrature_points[q_point](0));

                                for (unsigned int l = 0; l < 4; ++l)
                                  {
                                    system_rhs(j * deg + k, 2 * l) +=
                                      tmp(2 * l) * l_k_0;
                                    system_rhs(j * deg + k, 2 * l + 1) +=
                                      tmp(2 * l + 1) * l_k_1;
                                    system_rhs(j * deg + k, 2 * (l + 4)) +=
                                      tmp(2 * (l + 4)) * l_k_1;
                                    system_rhs(j * deg + k, 2 * l + 9) +=
                                      tmp(2 * l + 9) * l_k_0;
                                    system_rhs(j * deg + k, 2 * (l + 8)) +=
                                      tmp(2 * (l + 8)) * l_k_0;
                                    system_rhs(j * deg + k, 2 * l + 17) +=
                                      tmp(2 * l + 17) * l_k_1;
                                  }
                              }
                          }
                      }

                    system_matrix_inv.mmult(solution, system_rhs);

                    for (unsigned int j = 0; j < 2; ++j)
                      for (unsigned int k = 0; k < 2; ++k)
                        for (unsigned int l = 0; l <= deg; ++l)
                          for (unsigned int m = 0; m < deg; ++m)
                            {
                              if (std::abs(solution(l * deg + m,
                                                    2 * (j + 2 * k))) > 1e-14)
                                this->restriction[index][i + 2 * (2 * j + k)](
                                  (2 * i * this->degree + l) * deg + m +
                                    n_edge_dofs,
                                  dof) = solution(l * deg + m, 2 * (j + 2 * k));

                              if (std::abs(solution(l * deg + m,
                                                    2 * (j + 2 * k) + 1)) >
                                  1e-14)
                                this->restriction[index][i + 2 * (2 * j + k)](
                                  ((2 * i + 1) * deg + m) * this->degree + l +
                                    n_edge_dofs,
                                  dof) =
                                  solution(l * deg + m, 2 * (j + 2 * k) + 1);

                              if (std::abs(solution(l * deg + m,
                                                    2 * (j + 2 * (k + 2)))) >
                                  1e-14)
                                this->restriction[index][2 * (i + 2 * j) + k](
                                  (2 * (i + 2) * this->degree + l) * deg + m +
                                    n_edge_dofs,
                                  dof) =
                                  solution(l * deg + m, 2 * (j + 2 * (k + 2)));

                              if (std::abs(solution(l * deg + m,
                                                    2 * (j + 2 * k) + 9)) >
                                  1e-14)
                                this->restriction[index][2 * (i + 2 * j) + k](
                                  ((2 * i + 5) * deg + m) * this->degree + l +
                                    n_edge_dofs,
                                  dof) =
                                  solution(l * deg + m, 2 * (j + 2 * k) + 9);

                              if (std::abs(solution(l * deg + m,
                                                    2 * (j + 2 * (k + 4)))) >
                                  1e-14)
                                this->restriction[index][2 * (2 * i + j) + k](
                                  (2 * (i + 4) * this->degree + l) * deg + m +
                                    n_edge_dofs,
                                  dof) =
                                  solution(l * deg + m, 2 * (j + 2 * (k + 4)));

                              if (std::abs(solution(l * deg + m,
                                                    2 * (j + 2 * k) + 17)) >
                                  1e-14)
                                this->restriction[index][2 * (2 * i + j) + k](
                                  ((2 * i + 9) * deg + m) * this->degree + l +
                                    n_edge_dofs,
                                  dof) =
                                  solution(l * deg + m, 2 * (j + 2 * k) + 17);
                            }
                  }

              const QGauss<dim>              quadrature(2 * this->degree);
              const std::vector<Point<dim>> &quadrature_points =
                quadrature.get_points();
              const unsigned int n_boundary_dofs =
                2 * GeometryInfo<dim>::faces_per_cell * deg * this->degree +
                n_edge_dofs;
              const unsigned int n_quadrature_points = quadrature.size();

              {
                FullMatrix<double> assembling_matrix(deg * deg * this->degree,
                                                     n_quadrature_points);

                for (unsigned int q_point = 0; q_point < n_quadrature_points;
                     ++q_point)
                  {
                    const double weight = std::sqrt(quadrature.weight(q_point));

                    for (unsigned int i = 0; i <= deg; ++i)
                      {
                        const double L_i =
                          weight * legendre_polynomials[i].value(
                                     quadrature_points[q_point](0));

                        for (unsigned int j = 0; j < deg; ++j)
                          {
                            const double l_j =
                              L_i * lobatto_polynomials[j + 2].value(
                                      quadrature_points[q_point](1));

                            for (unsigned int k = 0; k < deg; ++k)
                              assembling_matrix((i * deg + j) * deg + k,
                                                q_point) =
                                l_j * lobatto_polynomials[k + 2].value(
                                        quadrature_points[q_point](2));
                          }
                      }
                  }

                FullMatrix<double> system_matrix(assembling_matrix.m(),
                                                 assembling_matrix.m());

                assembling_matrix.mTmult(system_matrix, assembling_matrix);
                system_matrix_inv.reinit(system_matrix.m(), system_matrix.m());
                system_matrix_inv.invert(system_matrix);
              }

              solution.reinit(system_matrix_inv.m(), 24);
              system_rhs.reinit(system_matrix_inv.m(), 24);
              tmp.reinit(24);

              for (unsigned int dof = 0; dof < this->n_dofs_per_cell(); ++dof)
                {
                  system_rhs = 0.0;

                  for (unsigned int q_point = 0; q_point < n_quadrature_points;
                       ++q_point)
                    {
                      tmp = 0.0;

                      if (quadrature_points[q_point](0) < 0.5)
                        {
                          if (quadrature_points[q_point](1) < 0.5)
                            {
                              if (quadrature_points[q_point](2) < 0.5)
                                {
                                  const Point<dim> quadrature_point(
                                    2.0 * quadrature_points[q_point](0),
                                    2.0 * quadrature_points[q_point](1),
                                    2.0 * quadrature_points[q_point](2));

                                  tmp(0) += 2.0 * this->shape_value_component(
                                                    dof, quadrature_point, 0);
                                  tmp(1) += 2.0 * this->shape_value_component(
                                                    dof, quadrature_point, 1);
                                  tmp(2) += 2.0 * this->shape_value_component(
                                                    dof, quadrature_point, 2);
                                }

                              else
                                {
                                  const Point<dim> quadrature_point(
                                    2.0 * quadrature_points[q_point](0),
                                    2.0 * quadrature_points[q_point](1),
                                    2.0 * quadrature_points[q_point](2) - 1.0);

                                  tmp(3) += 2.0 * this->shape_value_component(
                                                    dof, quadrature_point, 0);
                                  tmp(4) += 2.0 * this->shape_value_component(
                                                    dof, quadrature_point, 1);
                                  tmp(5) += 2.0 * this->shape_value_component(
                                                    dof, quadrature_point, 2);
                                }
                            }

                          else if (quadrature_points[q_point](2) < 0.5)
                            {
                              const Point<dim> quadrature_point(
                                2.0 * quadrature_points[q_point](0),
                                2.0 * quadrature_points[q_point](1) - 1.0,
                                2.0 * quadrature_points[q_point](2));

                              tmp(6) += 2.0 * this->shape_value_component(
                                                dof, quadrature_point, 0);
                              tmp(7) += 2.0 * this->shape_value_component(
                                                dof, quadrature_point, 1);
                              tmp(8) += 2.0 * this->shape_value_component(
                                                dof, quadrature_point, 2);
                            }

                          else
                            {
                              const Point<dim> quadrature_point(
                                2.0 * quadrature_points[q_point](0),
                                2.0 * quadrature_points[q_point](1) - 1.0,
                                2.0 * quadrature_points[q_point](2) - 1.0);

                              tmp(9) += 2.0 * this->shape_value_component(
                                                dof, quadrature_point, 0);
                              tmp(10) += 2.0 * this->shape_value_component(
                                                 dof, quadrature_point, 1);
                              tmp(11) += 2.0 * this->shape_value_component(
                                                 dof, quadrature_point, 2);
                            }
                        }

                      else if (quadrature_points[q_point](1) < 0.5)
                        {
                          if (quadrature_points[q_point](2) < 0.5)
                            {
                              const Point<dim> quadrature_point(
                                2.0 * quadrature_points[q_point](0) - 1.0,
                                2.0 * quadrature_points[q_point](1),
                                2.0 * quadrature_points[q_point](2));

                              tmp(12) += 2.0 * this->shape_value_component(
                                                 dof, quadrature_point, 0);
                              tmp(13) += 2.0 * this->shape_value_component(
                                                 dof, quadrature_point, 1);
                              tmp(14) += 2.0 * this->shape_value_component(
                                                 dof, quadrature_point, 2);
                            }

                          else
                            {
                              const Point<dim> quadrature_point(
                                2.0 * quadrature_points[q_point](0) - 1.0,
                                2.0 * quadrature_points[q_point](1),
                                2.0 * quadrature_points[q_point](2) - 1.0);

                              tmp(15) += 2.0 * this->shape_value_component(
                                                 dof, quadrature_point, 0);
                              tmp(16) += 2.0 * this->shape_value_component(
                                                 dof, quadrature_point, 1);
                              tmp(17) += 2.0 * this->shape_value_component(
                                                 dof, quadrature_point, 2);
                            }
                        }

                      else if (quadrature_points[q_point](2) < 0.5)
                        {
                          const Point<dim> quadrature_point(
                            2.0 * quadrature_points[q_point](0) - 1.0,
                            2.0 * quadrature_points[q_point](1) - 1.0,
                            2.0 * quadrature_points[q_point](2));

                          tmp(18) +=
                            2.0 * this->shape_value_component(dof,
                                                              quadrature_point,
                                                              0);
                          tmp(19) +=
                            2.0 * this->shape_value_component(dof,
                                                              quadrature_point,
                                                              1);
                          tmp(20) +=
                            2.0 * this->shape_value_component(dof,
                                                              quadrature_point,
                                                              2);
                        }

                      else
                        {
                          const Point<dim> quadrature_point(
                            2.0 * quadrature_points[q_point](0) - 1.0,
                            2.0 * quadrature_points[q_point](1) - 1.0,
                            2.0 * quadrature_points[q_point](2) - 1.0);

                          tmp(21) +=
                            2.0 * this->shape_value_component(dof,
                                                              quadrature_point,
                                                              0);
                          tmp(22) +=
                            2.0 * this->shape_value_component(dof,
                                                              quadrature_point,
                                                              1);
                          tmp(23) +=
                            2.0 * this->shape_value_component(dof,
                                                              quadrature_point,
                                                              2);
                        }

                      for (unsigned int i = 0; i < 2; ++i)
                        for (unsigned int j = 0; j < 2; ++j)
                          for (unsigned int k = 0; k < 2; ++k)
                            for (unsigned int l = 0; l <= deg; ++l)
                              {
                                tmp(3 * (i + 2 * (j + 2 * k))) -=
                                  this->restriction[index][2 * (2 * i + j) + k](
                                    (4 * i + j + 2) * this->degree + l, dof) *
                                  this->shape_value_component(
                                    (4 * i + j + 2) * this->degree + l,
                                    quadrature_points[q_point],
                                    0);
                                tmp(3 * (i + 2 * (j + 2 * k)) + 1) -=
                                  this->restriction[index][2 * (2 * i + j) + k](
                                    (4 * i + k) * this->degree + l, dof) *
                                  this->shape_value_component(
                                    (4 * i + k) * this->degree + l,
                                    quadrature_points[q_point],
                                    1);
                                tmp(3 * (i + 2 * (j + 2 * k)) + 2) -=
                                  this->restriction[index][2 * (2 * i + j) + k](
                                    (2 * (j + 4) + k) * this->degree + l, dof) *
                                  this->shape_value_component(
                                    (2 * (j + 4) + k) * this->degree + l,
                                    quadrature_points[q_point],
                                    2);

                                for (unsigned int m = 0; m < deg; ++m)
                                  {
                                    tmp(3 * (i + 2 * (j + 2 * k))) -=
                                      this->restriction[index][2 * (2 * i + j) +
                                                               k](
                                        ((2 * j + 5) * deg + m) * this->degree +
                                          l + n_edge_dofs,
                                        dof) *
                                      this->shape_value_component(
                                        ((2 * j + 5) * deg + m) * this->degree +
                                          l + n_edge_dofs,
                                        quadrature_points[q_point],
                                        0);
                                    tmp(3 * (i + 2 * (j + 2 * k))) -=
                                      this->restriction[index][2 * (2 * i + j) +
                                                               k](
                                        (2 * (i + 4) * this->degree + l) * deg +
                                          m + n_edge_dofs,
                                        dof) *
                                      this->shape_value_component(
                                        (2 * (i + 4) * this->degree + l) * deg +
                                          m + n_edge_dofs,
                                        quadrature_points[q_point],
                                        0);
                                    tmp(3 * (i + 2 * (j + 2 * k)) + 1) -=
                                      this->restriction[index][2 * (2 * i + j) +
                                                               k](
                                        (2 * k * this->degree + l) * deg + m +
                                          n_edge_dofs,
                                        dof) *
                                      this->shape_value_component(
                                        (2 * k * this->degree + l) * deg + m +
                                          n_edge_dofs,
                                        quadrature_points[q_point],
                                        1);
                                    tmp(3 * (i + 2 * (j + 2 * k)) + 1) -=
                                      this->restriction[index][2 * (2 * i + j) +
                                                               k](
                                        ((2 * i + 9) * deg + m) * this->degree +
                                          l + n_edge_dofs,
                                        dof) *
                                      this->shape_value_component(
                                        ((2 * i + 9) * deg + m) * this->degree +
                                          l + n_edge_dofs,
                                        quadrature_points[q_point],
                                        1);
                                    tmp(3 * (i + 2 * (j + 2 * k)) + 2) -=
                                      this->restriction[index][2 * (2 * i + j) +
                                                               k](
                                        ((2 * k + 1) * deg + m) * this->degree +
                                          l + n_edge_dofs,
                                        dof) *
                                      this->shape_value_component(
                                        ((2 * k + 1) * deg + m) * this->degree +
                                          l + n_edge_dofs,
                                        quadrature_points[q_point],
                                        2);
                                    tmp(3 * (i + 2 * (j + 2 * k)) + 2) -=
                                      this->restriction[index][2 * (2 * i + j) +
                                                               k](
                                        (2 * (j + 2) * this->degree + l) * deg +
                                          m + n_edge_dofs,
                                        dof) *
                                      this->shape_value_component(
                                        (2 * (j + 2) * this->degree + l) * deg +
                                          m + n_edge_dofs,
                                        quadrature_points[q_point],
                                        2);
                                  }
                              }

                      tmp *= quadrature.weight(q_point);

                      for (unsigned int i = 0; i <= deg; ++i)
                        {
                          const double L_i_0 = legendre_polynomials[i].value(
                            quadrature_points[q_point](0));
                          const double L_i_1 = legendre_polynomials[i].value(
                            quadrature_points[q_point](1));
                          const double L_i_2 = legendre_polynomials[i].value(
                            quadrature_points[q_point](2));

                          for (unsigned int j = 0; j < deg; ++j)
                            {
                              const double l_j_0 =
                                L_i_0 * lobatto_polynomials[j + 2].value(
                                          quadrature_points[q_point](1));
                              const double l_j_1 =
                                L_i_1 * lobatto_polynomials[j + 2].value(
                                          quadrature_points[q_point](0));
                              const double l_j_2 =
                                L_i_2 * lobatto_polynomials[j + 2].value(
                                          quadrature_points[q_point](0));

                              for (unsigned int k = 0; k < deg; ++k)
                                {
                                  const double l_k_0 =
                                    l_j_0 * lobatto_polynomials[k + 2].value(
                                              quadrature_points[q_point](2));
                                  const double l_k_1 =
                                    l_j_1 * lobatto_polynomials[k + 2].value(
                                              quadrature_points[q_point](2));
                                  const double l_k_2 =
                                    l_j_2 * lobatto_polynomials[k + 2].value(
                                              quadrature_points[q_point](1));

                                  for (unsigned int l = 0; l < 8; ++l)
                                    {
                                      system_rhs((i * deg + j) * deg + k,
                                                 3 * l) += tmp(3 * l) * l_k_0;
                                      system_rhs((i * deg + j) * deg + k,
                                                 3 * l + 1) +=
                                        tmp(3 * l + 1) * l_k_1;
                                      system_rhs((i * deg + j) * deg + k,
                                                 3 * l + 2) +=
                                        tmp(3 * l + 2) * l_k_2;
                                    }
                                }
                            }
                        }
                    }

                  system_matrix_inv.mmult(solution, system_rhs);

                  for (unsigned int i = 0; i < 2; ++i)
                    for (unsigned int j = 0; j < 2; ++j)
                      for (unsigned int k = 0; k < 2; ++k)
                        for (unsigned int l = 0; l <= deg; ++l)
                          for (unsigned int m = 0; m < deg; ++m)
                            for (unsigned int n = 0; n < deg; ++n)
                              {
                                if (std::abs(
                                      solution((l * deg + m) * deg + n,
                                               3 * (i + 2 * (j + 2 * k)))) >
                                    1e-14)
                                  this->restriction[index][2 * (2 * i + j) + k](
                                    (l * deg + m) * deg + n + n_boundary_dofs,
                                    dof) = solution((l * deg + m) * deg + n,
                                                    3 * (i + 2 * (j + 2 * k)));

                                if (std::abs(
                                      solution((l * deg + m) * deg + n,
                                               3 * (i + 2 * (j + 2 * k)) + 1)) >
                                    1e-14)
                                  this->restriction[index][2 * (2 * i + j) + k](
                                    (l + (m + deg) * this->degree) * deg + n +
                                      n_boundary_dofs,
                                    dof) =
                                    solution((l * deg + m) * deg + n,
                                             3 * (i + 2 * (j + 2 * k)) + 1);

                                if (std::abs(
                                      solution((l * deg + m) * deg + n,
                                               3 * (i + 2 * (j + 2 * k)) + 2)) >
                                    1e-14)
                                  this->restriction[index][2 * (2 * i + j) + k](
                                    l +
                                      ((m + 2 * deg) * deg + n) * this->degree +
                                      n_boundary_dofs,
                                    dof) =
                                    solution((l * deg + m) * deg + n,
                                             3 * (i + 2 * (j + 2 * k)) + 2);
                              }
                }
            }

          break;
        }

      default:
        Assert(false, ExcNotImplemented());
    }
}



template <int dim>
std::vector<unsigned int>
FE_Nedelec<dim>::get_dpo_vector(const unsigned int degree, bool dg)
{
  std::vector<unsigned int> dpo(dim + 1);

  if (dg)
    {
      dpo[dim] = PolynomialsNedelec<dim>::n_polynomials(degree);
    }
  else
    {
      dpo[0] = 0;
      dpo[1] = degree + 1;
      dpo[2] = 2 * degree * (degree + 1);

      if (dim == 3)
        dpo[3] = 3 * degree * degree * (degree + 1);
    }

  return dpo;
}

//---------------------------------------------------------------------------
// Data field initialization
//---------------------------------------------------------------------------

// Check whether a given shape
// function has support on a
// given face.

// We just switch through the
// faces of the cell and return
// true, if the shape function
// has support on the face
// and false otherwise.
template <int dim>
bool
FE_Nedelec<dim>::has_support_on_face(const unsigned int shape_index,
                                     const unsigned int face_index) const
{
  AssertIndexRange(shape_index, this->n_dofs_per_cell());
  AssertIndexRange(face_index, GeometryInfo<dim>::faces_per_cell);

  const unsigned int deg = this->degree - 1;
  switch (dim)
    {
      case 2:
        switch (face_index)
          {
            case 0:
              if (!((shape_index > deg) && (shape_index < 2 * this->degree)))
                return true;

              else
                return false;

            case 1:
              if ((shape_index > deg) &&
                  (shape_index <
                   GeometryInfo<2>::lines_per_cell * this->degree))
                return true;

              else
                return false;

            case 2:
              if (shape_index < 3 * this->degree)
                return true;

              else
                return false;

            case 3:
              if (!((shape_index >= 2 * this->degree) &&
                    (shape_index < 3 * this->degree)))
                return true;

              else
                return false;

            default:
              {
                Assert(false, ExcNotImplemented());
                return false;
              }
          }

      case 3:
        switch (face_index)
          {
            case 0:
              if (((shape_index > deg) && (shape_index < 2 * this->degree)) ||
                  ((shape_index >= 5 * this->degree) &&
                   (shape_index < 6 * this->degree)) ||
                  ((shape_index >= 9 * this->degree) &&
                   (shape_index < 10 * this->degree)) ||
                  ((shape_index >= 11 * this->degree) &&
                   (shape_index <
                    GeometryInfo<3>::lines_per_cell * this->degree)) ||
                  ((shape_index >= (GeometryInfo<3>::lines_per_cell + 2 * deg) *
                                     this->degree) &&
                   (shape_index < (GeometryInfo<3>::lines_per_cell + 4 * deg) *
                                    this->degree)) ||
                  ((shape_index >= (GeometryInfo<3>::lines_per_cell + 5 * deg) *
                                     this->degree) &&
                   (shape_index < (GeometryInfo<3>::lines_per_cell + 6 * deg) *
                                    this->degree)) ||
                  ((shape_index >= (GeometryInfo<3>::lines_per_cell + 7 * deg) *
                                     this->degree) &&
                   (shape_index < (GeometryInfo<3>::lines_per_cell + 9 * deg) *
                                    this->degree)) ||
                  ((shape_index >=
                    (GeometryInfo<3>::lines_per_cell + 10 * deg) *
                      this->degree) &&
                   (shape_index < (GeometryInfo<3>::lines_per_cell + 11 * deg) *
                                    this->degree)))
                return false;

              else
                return true;

            case 1:
              if (((shape_index > deg) && (shape_index < 4 * this->degree)) ||
                  ((shape_index >= 5 * this->degree) &&
                   (shape_index < 8 * this->degree)) ||
                  ((shape_index >= 9 * this->degree) &&
                   (shape_index < 10 * this->degree)) ||
                  ((shape_index >= 11 * this->degree) &&
                   (shape_index <
                    GeometryInfo<3>::lines_per_cell * this->degree)) ||
                  ((shape_index >= (GeometryInfo<3>::lines_per_cell + 2 * deg) *
                                     this->degree) &&
                   (shape_index < (GeometryInfo<3>::lines_per_cell + 5 * deg) *
                                    this->degree)) ||
                  ((shape_index >= (GeometryInfo<3>::lines_per_cell + 6 * deg) *
                                     this->degree) &&
                   (shape_index < (GeometryInfo<3>::lines_per_cell + 7 * deg) *
                                    this->degree)) ||
                  ((shape_index >= (GeometryInfo<3>::lines_per_cell + 9 * deg) *
                                     this->degree) &&
                   (shape_index < (GeometryInfo<3>::lines_per_cell + 10 * deg) *
                                    this->degree)) ||
                  ((shape_index >=
                    (GeometryInfo<3>::lines_per_cell + 11 * deg) *
                      this->degree) &&
                   (shape_index < (GeometryInfo<3>::lines_per_cell + 12 * deg) *
                                    this->degree)))
                return true;

              else
                return false;

            case 2:
              if ((shape_index < 3 * this->degree) ||
                  ((shape_index >= 4 * this->degree) &&
                   (shape_index < 7 * this->degree)) ||
                  ((shape_index >= 8 * this->degree) &&
                   (shape_index < 10 * this->degree)) ||
                  ((shape_index >=
                    (GeometryInfo<3>::lines_per_cell + deg) * this->degree) &&
                   (shape_index < (GeometryInfo<3>::lines_per_cell + 2 * deg) *
                                    this->degree)) ||
                  ((shape_index >= (GeometryInfo<3>::lines_per_cell + 3 * deg) *
                                     this->degree) &&
                   (shape_index < (GeometryInfo<3>::lines_per_cell + 6 * deg) *
                                    this->degree)) ||
                  ((shape_index >= (GeometryInfo<3>::lines_per_cell + 8 * deg) *
                                     this->degree) &&
                   (shape_index < (GeometryInfo<3>::lines_per_cell + 9 * deg) *
                                    this->degree)) ||
                  ((shape_index >=
                    (GeometryInfo<3>::lines_per_cell + 10 * deg) *
                      this->degree) &&
                   (shape_index < (GeometryInfo<3>::lines_per_cell + 11 * deg) *
                                    this->degree)))
                return true;

              else
                return false;

            case 3:
              if ((shape_index < 2 * this->degree) ||
                  ((shape_index >= 3 * this->degree) &&
                   (shape_index < 6 * this->degree)) ||
                  ((shape_index >= 7 * this->degree) &&
                   (shape_index < 8 * this->degree)) ||
                  ((shape_index >= 10 * this->degree) &&
                   (shape_index <
                    GeometryInfo<3>::lines_per_cell * this->degree)) ||
                  ((shape_index >=
                    (GeometryInfo<3>::lines_per_cell + deg) * this->degree) &&
                   (shape_index < (GeometryInfo<3>::lines_per_cell + 2 * deg) *
                                    this->degree)) ||
                  ((shape_index >= (GeometryInfo<3>::lines_per_cell + 3 * deg) *
                                     this->degree) &&
                   (shape_index < (GeometryInfo<3>::lines_per_cell + 4 * deg) *
                                    this->degree)) ||
                  ((shape_index >= (GeometryInfo<3>::lines_per_cell + 6 * deg) *
                                     this->degree) &&
                   (shape_index < (GeometryInfo<3>::lines_per_cell + 9 * deg) *
                                    this->degree)) ||
                  ((shape_index >=
                    (GeometryInfo<3>::lines_per_cell + 10 * deg) *
                      this->degree) &&
                   (shape_index < (GeometryInfo<3>::lines_per_cell + 11 * deg) *
                                    this->degree)))
                return true;

              else
                return false;

            case 4:
              if ((shape_index < 4 * this->degree) ||
                  ((shape_index >= 8 * this->degree) &&
                   (shape_index <
                    (GeometryInfo<3>::lines_per_cell + deg) * this->degree)) ||
                  ((shape_index >= (GeometryInfo<3>::lines_per_cell + 2 * deg) *
                                     this->degree) &&
                   (shape_index < (GeometryInfo<3>::lines_per_cell + 3 * deg) *
                                    this->degree)) ||
                  ((shape_index >= (GeometryInfo<3>::lines_per_cell + 5 * deg) *
                                     this->degree) &&
                   (shape_index < (GeometryInfo<3>::lines_per_cell + 6 * deg) *
                                    this->degree)) ||
                  ((shape_index >= (GeometryInfo<3>::lines_per_cell + 7 * deg) *
                                     this->degree) &&
                   (shape_index < (GeometryInfo<3>::lines_per_cell + 10 * deg) *
                                    this->degree)))
                return true;

              else
                return false;

            case 5:
              if (((shape_index >= 4 * this->degree) &&
                   (shape_index <
                    (GeometryInfo<3>::lines_per_cell + deg) * this->degree)) ||
                  ((shape_index >= (GeometryInfo<3>::lines_per_cell + 2 * deg) *
                                     this->degree) &&
                   (shape_index < (GeometryInfo<3>::lines_per_cell + 3 * deg) *
                                    this->degree)) ||
                  ((shape_index >= (GeometryInfo<3>::lines_per_cell + 5 * deg) *
                                     this->degree) &&
                   (shape_index < (GeometryInfo<3>::lines_per_cell + 6 * deg) *
                                    this->degree)) ||
                  ((shape_index >= (GeometryInfo<3>::lines_per_cell + 7 * deg) *
                                     this->degree) &&
                   (shape_index < (GeometryInfo<3>::lines_per_cell + 8 * deg) *
                                    this->degree)) ||
                  ((shape_index >=
                    (GeometryInfo<3>::lines_per_cell + 10 * deg) *
                      this->degree) &&
                   (shape_index < (GeometryInfo<3>::lines_per_cell + 12 * deg) *
                                    this->degree)))
                return true;

              else
                return false;

            default:
              {
                Assert(false, ExcNotImplemented());
                return false;
              }
          }

      default:
        {
          Assert(false, ExcNotImplemented());
          return false;
        }
    }
}

template <int dim>
FiniteElementDomination::Domination
FE_Nedelec<dim>::compare_for_domination(const FiniteElement<dim> &fe_other,
                                        const unsigned int        codim) const
{
  Assert(codim <= dim, ExcImpossibleInDim(dim));
  (void)codim;

  // vertex/line/face/cell domination
  // --------------------------------
  if (const FE_Nedelec<dim> *fe_nedelec_other =
        dynamic_cast<const FE_Nedelec<dim> *>(&fe_other))
    {
      if (this->degree < fe_nedelec_other->degree)
        return FiniteElementDomination::this_element_dominates;
      else if (this->degree == fe_nedelec_other->degree)
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

template <int dim>
bool
FE_Nedelec<dim>::hp_constraints_are_implemented() const
{
  return true;
}

template <int dim>
std::vector<std::pair<unsigned int, unsigned int>>
FE_Nedelec<dim>::hp_vertex_dof_identities(const FiniteElement<dim> &) const
{
  // Nedelec elements do not have any dofs
  // on vertices, hence return an empty vector.
  return std::vector<std::pair<unsigned int, unsigned int>>();
}

template <int dim>
std::vector<std::pair<unsigned int, unsigned int>>
FE_Nedelec<dim>::hp_line_dof_identities(
  const FiniteElement<dim> &fe_other) const
{
  // we can presently only compute these
  // identities if both FEs are
  // FE_Nedelec or if the other one is an
  // FE_Nothing
  if (const FE_Nedelec<dim> *fe_nedelec_other =
        dynamic_cast<const FE_Nedelec<dim> *>(&fe_other))
    {
      // dofs are located on lines, so
      // two dofs are identical, if their
      // edge shape functions have the
      // same polynomial degree.
      std::vector<std::pair<unsigned int, unsigned int>> identities;

      for (unsigned int i = 0;
           i < std::min(fe_nedelec_other->degree, this->degree);
           ++i)
        identities.emplace_back(i, i);

      return identities;
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
      Assert(false, ExcNotImplemented());
      return std::vector<std::pair<unsigned int, unsigned int>>();
    }
}

template <int dim>
std::vector<std::pair<unsigned int, unsigned int>>
FE_Nedelec<dim>::hp_quad_dof_identities(const FiniteElement<dim> &fe_other,
                                        const unsigned int) const
{
  // we can presently only compute
  // these identities if both FEs are
  // FE_Nedelec or if the other one is an
  // FE_Nothing
  if (const FE_Nedelec<dim> *fe_nedelec_other =
        dynamic_cast<const FE_Nedelec<dim> *>(&fe_other))
    {
      // dofs are located on the interior
      // of faces, so two dofs are identical,
      // if their face shape functions have
      // the same polynomial degree.
      const unsigned int p     = fe_nedelec_other->degree;
      const unsigned int q     = this->degree;
      const unsigned int p_min = std::min(p, q);
      std::vector<std::pair<unsigned int, unsigned int>> identities;

      for (unsigned int i = 0; i < p_min; ++i)
        for (unsigned int j = 0; j < p_min - 1; ++j)
          {
            identities.emplace_back(i * (q - 1) + j, i * (p - 1) + j);
            identities.emplace_back(i + (j + q - 1) * q, i + (j + p - 1) * p);
          }

      return identities;
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
      Assert(false, ExcNotImplemented());
      return std::vector<std::pair<unsigned int, unsigned int>>();
    }
}

// In this function we compute the face
// interpolation matrix. This is usually
// done by projection-based interpolation,
// but, since one can compute the entries
// easy per hand, we save some computation
// time at this point and just fill in the
// correct values.
template <int dim>
void
FE_Nedelec<dim>::get_face_interpolation_matrix(
  const FiniteElement<dim> &source,
  FullMatrix<double> &      interpolation_matrix,
  const unsigned int        face_no) const
{
  (void)face_no;
  // this is only implemented, if the
  // source FE is also a
  // Nedelec element
  AssertThrow((source.get_name().find("FE_Nedelec<") == 0) ||
                (dynamic_cast<const FE_Nedelec<dim> *>(&source) != nullptr),
              (typename FiniteElement<dim>::ExcInterpolationNotImplemented()));
  Assert(interpolation_matrix.m() == source.n_dofs_per_face(face_no),
         ExcDimensionMismatch(interpolation_matrix.m(),
                              source.n_dofs_per_face(face_no)));
  Assert(interpolation_matrix.n() == this->n_dofs_per_face(face_no),
         ExcDimensionMismatch(interpolation_matrix.n(),
                              this->n_dofs_per_face(face_no)));

  // ok, source is a Nedelec element, so
  // we will be able to do the work
  const FE_Nedelec<dim> &source_fe =
    dynamic_cast<const FE_Nedelec<dim> &>(source);

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

  // On lines we can just identify
  // all degrees of freedom.
  for (unsigned int i = 0; i < this->degree; ++i)
    interpolation_matrix(i, i) = 1.0;

  // In 3d we have some lines more
  // and a face. The procedure stays
  // the same as above, but we have
  // to take a bit more care of the
  // indices of the degrees of
  // freedom.
  if (dim == 3)
    {
      const unsigned int p = source_fe.degree;
      const unsigned int q = this->degree;

      for (unsigned int i = 0; i < q; ++i)
        {
          for (unsigned int j = 1; j < GeometryInfo<dim>::lines_per_face; ++j)
            interpolation_matrix(j * p + i, j * q + i) = 1.0;

          for (unsigned int j = 0; j < q - 1; ++j)
            {
              interpolation_matrix(GeometryInfo<dim>::lines_per_face * p +
                                     i * (p - 1) + j,
                                   GeometryInfo<dim>::lines_per_face * q +
                                     i * (q - 1) + j) = 1.0;
              interpolation_matrix(GeometryInfo<dim>::lines_per_face * p + i +
                                     (j + p - 1) * p,
                                   GeometryInfo<dim>::lines_per_face * q + i +
                                     (j + q - 1) * q) = 1.0;
            }
        }
    }
}



template <>
void
FE_Nedelec<1>::get_subface_interpolation_matrix(const FiniteElement<1, 1> &,
                                                const unsigned int,
                                                FullMatrix<double> &,
                                                const unsigned int) const
{
  Assert(false, ExcNotImplemented());
}



// In this function we compute the
// subface interpolation matrix.
// This is done by a projection-
// based interpolation. Therefore
// we first interpolate the
// shape functions of the higher
// order element on the lowest
// order edge shape functions.
// Then the remaining part of
// the interpolated shape
// functions is projected on the
// higher order edge shape
// functions, the face shape
// functions and the interior
// shape functions (if they all
// exist).
template <int dim>
void
FE_Nedelec<dim>::get_subface_interpolation_matrix(
  const FiniteElement<dim> &source,
  const unsigned int        subface,
  FullMatrix<double> &      interpolation_matrix,
  const unsigned int        face_no) const
{
  // this is only implemented, if the
  // source FE is also a
  // Nedelec element
  AssertThrow((source.get_name().find("FE_Nedelec<") == 0) ||
                (dynamic_cast<const FE_Nedelec<dim> *>(&source) != nullptr),
              typename FiniteElement<dim>::ExcInterpolationNotImplemented());
  Assert(interpolation_matrix.m() == source.n_dofs_per_face(face_no),
         ExcDimensionMismatch(interpolation_matrix.m(),
                              source.n_dofs_per_face(face_no)));
  Assert(interpolation_matrix.n() == this->n_dofs_per_face(face_no),
         ExcDimensionMismatch(interpolation_matrix.n(),
                              this->n_dofs_per_face(face_no)));

  // ok, source is a Nedelec element, so
  // we will be able to do the work
  const FE_Nedelec<dim> &source_fe =
    dynamic_cast<const FE_Nedelec<dim> &>(source);

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
  interpolation_matrix = 0.0;
  // Perform projection-based interpolation
  // as usual.
  const QGauss<1>              edge_quadrature(source_fe.degree);
  const std::vector<Point<1>> &edge_quadrature_points =
    edge_quadrature.get_points();
  const unsigned int n_edge_quadrature_points = edge_quadrature.size();

  switch (dim)
    {
      case 2:
        {
          for (unsigned int dof = 0; dof < this->n_dofs_per_face(face_no);
               ++dof)
            for (unsigned int q_point = 0; q_point < n_edge_quadrature_points;
                 ++q_point)
              {
                const Point<dim> quadrature_point(
                  0.0, 0.5 * (edge_quadrature_points[q_point](0) + subface));

                interpolation_matrix(0, dof) +=
                  0.5 * edge_quadrature.weight(q_point) *
                  this->shape_value_component(dof, quadrature_point, 1);
              }

          if (source_fe.degree > 1)
            {
              const std::vector<Polynomials::Polynomial<double>>
                &legendre_polynomials =
                  Polynomials::Legendre::generate_complete_basis(
                    source_fe.degree - 1);
              FullMatrix<double> system_matrix_inv(source_fe.degree - 1,
                                                   source_fe.degree - 1);

              {
                FullMatrix<double> assembling_matrix(source_fe.degree - 1,
                                                     n_edge_quadrature_points);

                for (unsigned int q_point = 0;
                     q_point < n_edge_quadrature_points;
                     ++q_point)
                  {
                    const double weight =
                      std::sqrt(edge_quadrature.weight(q_point));

                    for (unsigned int i = 0; i < source_fe.degree - 1; ++i)
                      assembling_matrix(i, q_point) =
                        weight * legendre_polynomials[i + 1].value(
                                   edge_quadrature_points[q_point](0));
                  }

                FullMatrix<double> system_matrix(source_fe.degree - 1,
                                                 source_fe.degree - 1);

                assembling_matrix.mTmult(system_matrix, assembling_matrix);
                system_matrix_inv.invert(system_matrix);
              }

              Vector<double> solution(source_fe.degree - 1);
              Vector<double> system_rhs(source_fe.degree - 1);

              for (unsigned int dof = 0; dof < this->n_dofs_per_face(face_no);
                   ++dof)
                {
                  system_rhs = 0.0;

                  for (unsigned int q_point = 0;
                       q_point < n_edge_quadrature_points;
                       ++q_point)
                    {
                      const Point<dim> quadrature_point_0(
                        0.0,
                        0.5 * (edge_quadrature_points[q_point](0) + subface));
                      const Point<dim> quadrature_point_1(
                        0.0, edge_quadrature_points[q_point](0));
                      const double tmp =
                        edge_quadrature.weight(q_point) *
                        (0.5 * this->shape_value_component(dof,
                                                           quadrature_point_0,
                                                           1) -
                         interpolation_matrix(0, dof) *
                           source_fe.shape_value_component(0,
                                                           quadrature_point_1,
                                                           1));

                      for (unsigned int i = 0; i < source_fe.degree - 1; ++i)
                        system_rhs(i) +=
                          tmp * legendre_polynomials[i + 1].value(
                                  edge_quadrature_points[q_point](0));
                    }

                  system_matrix_inv.vmult(solution, system_rhs);

                  for (unsigned int i = 0; i < source_fe.degree - 1; ++i)
                    if (std::abs(solution(i)) > 1e-14)
                      interpolation_matrix(i + 1, dof) = solution(i);
                }
            }

          break;
        }

      case 3:
        {
          const double shifts[4][2] = {{0.0, 0.0},
                                       {1.0, 0.0},
                                       {0.0, 1.0},
                                       {1.0, 1.0}};

          for (unsigned int dof = 0; dof < this->n_dofs_per_face(face_no);
               ++dof)
            for (unsigned int q_point = 0; q_point < n_edge_quadrature_points;
                 ++q_point)
              {
                const double weight = 0.5 * edge_quadrature.weight(q_point);

                for (unsigned int i = 0; i < 2; ++i)
                  {
                    Point<dim> quadrature_point(
                      0.5 * (i + shifts[subface][0]),
                      0.5 * (edge_quadrature_points[q_point](0) +
                             shifts[subface][1]),
                      0.0);

                    interpolation_matrix(i * source_fe.degree, dof) +=
                      weight *
                      this->shape_value_component(
                        this->face_to_cell_index(dof, 4), quadrature_point, 1);
                    quadrature_point =
                      Point<dim>(0.5 * (edge_quadrature_points[q_point](0) +
                                        shifts[subface][0]),
                                 0.5 * (i + shifts[subface][1]),
                                 0.0);
                    interpolation_matrix((i + 2) * source_fe.degree, dof) +=
                      weight *
                      this->shape_value_component(
                        this->face_to_cell_index(dof, 4), quadrature_point, 0);
                  }
              }

          if (source_fe.degree > 1)
            {
              const std::vector<Polynomials::Polynomial<double>>
                &legendre_polynomials =
                  Polynomials::Legendre::generate_complete_basis(
                    source_fe.degree - 1);
              FullMatrix<double> system_matrix_inv(source_fe.degree - 1,
                                                   source_fe.degree - 1);

              {
                FullMatrix<double> assembling_matrix(source_fe.degree - 1,
                                                     n_edge_quadrature_points);

                for (unsigned int q_point = 0;
                     q_point < n_edge_quadrature_points;
                     ++q_point)
                  {
                    const double weight =
                      std::sqrt(edge_quadrature.weight(q_point));

                    for (unsigned int i = 0; i < source_fe.degree - 1; ++i)
                      assembling_matrix(i, q_point) =
                        weight * legendre_polynomials[i + 1].value(
                                   edge_quadrature_points[q_point](0));
                  }

                FullMatrix<double> system_matrix(source_fe.degree - 1,
                                                 source_fe.degree - 1);

                assembling_matrix.mTmult(system_matrix, assembling_matrix);
                system_matrix_inv.invert(system_matrix);
              }

              FullMatrix<double> solution(source_fe.degree - 1,
                                          GeometryInfo<dim>::lines_per_face);
              FullMatrix<double> system_rhs(source_fe.degree - 1,
                                            GeometryInfo<dim>::lines_per_face);
              Vector<double>     tmp(GeometryInfo<dim>::lines_per_face);

              for (unsigned int dof = 0; dof < this->n_dofs_per_face(face_no);
                   ++dof)
                {
                  system_rhs = 0.0;

                  for (unsigned int q_point = 0;
                       q_point < n_edge_quadrature_points;
                       ++q_point)
                    {
                      const double weight = edge_quadrature.weight(q_point);

                      for (unsigned int i = 0; i < 2; ++i)
                        {
                          Point<dim> quadrature_point_0(
                            0.5 * (i + shifts[subface][0]),
                            0.5 * (edge_quadrature_points[q_point](0) +
                                   shifts[subface][1]),
                            0.0);
                          Point<dim> quadrature_point_1(
                            i, edge_quadrature_points[q_point](0), 0.0);

                          tmp(i) =
                            weight *
                            (0.5 * this->shape_value_component(
                                     this->face_to_cell_index(dof, 4),
                                     quadrature_point_0,
                                     1) -
                             interpolation_matrix(i * source_fe.degree, dof) *
                               source_fe.shape_value_component(
                                 i * source_fe.degree, quadrature_point_1, 1));
                          quadrature_point_0 =
                            Point<dim>(0.5 *
                                         (edge_quadrature_points[q_point](0) +
                                          shifts[subface][0]),
                                       0.5 * (i + shifts[subface][1]),
                                       0.0);
                          quadrature_point_1 =
                            Point<dim>(edge_quadrature_points[q_point](0),
                                       i,
                                       0.0);
                          tmp(i + 2) =
                            weight *
                            (0.5 * this->shape_value_component(
                                     this->face_to_cell_index(dof, 4),
                                     quadrature_point_0,
                                     0) -
                             interpolation_matrix((i + 2) * source_fe.degree,
                                                  dof) *
                               source_fe.shape_value_component(
                                 (i + 2) * source_fe.degree,
                                 quadrature_point_1,
                                 0));
                        }

                      for (unsigned int i = 0; i < source_fe.degree - 1; ++i)
                        {
                          const double L_i = legendre_polynomials[i + 1].value(
                            edge_quadrature_points[q_point](0));

                          for (unsigned int j = 0;
                               j < GeometryInfo<dim>::lines_per_face;
                               ++j)
                            system_rhs(i, j) += tmp(j) * L_i;
                        }
                    }

                  system_matrix_inv.mmult(solution, system_rhs);

                  for (unsigned int i = 0;
                       i < GeometryInfo<dim>::lines_per_face;
                       ++i)
                    for (unsigned int j = 0; j < source_fe.degree - 1; ++j)
                      if (std::abs(solution(j, i)) > 1e-14)
                        interpolation_matrix(i * source_fe.degree + j + 1,
                                             dof) = solution(j, i);
                }

              const QGauss<2>              quadrature(source_fe.degree);
              const std::vector<Point<2>> &quadrature_points =
                quadrature.get_points();
              const std::vector<Polynomials::Polynomial<double>>
                &lobatto_polynomials =
                  Polynomials::Lobatto::generate_complete_basis(
                    source_fe.degree);
              const unsigned int n_boundary_dofs =
                GeometryInfo<dim>::lines_per_face * source_fe.degree;
              const unsigned int n_quadrature_points = quadrature.size();

              {
                FullMatrix<double> assembling_matrix(source_fe.degree *
                                                       (source_fe.degree - 1),
                                                     n_quadrature_points);

                for (unsigned int q_point = 0; q_point < n_quadrature_points;
                     ++q_point)
                  {
                    const double weight = std::sqrt(quadrature.weight(q_point));

                    for (unsigned int i = 0; i < source_fe.degree; ++i)
                      {
                        const double L_i =
                          weight * legendre_polynomials[i].value(
                                     quadrature_points[q_point](0));

                        for (unsigned int j = 0; j < source_fe.degree - 1; ++j)
                          assembling_matrix(i * (source_fe.degree - 1) + j,
                                            q_point) =
                            L_i * lobatto_polynomials[j + 2].value(
                                    quadrature_points[q_point](1));
                      }
                  }

                FullMatrix<double> system_matrix(assembling_matrix.m(),
                                                 assembling_matrix.m());

                assembling_matrix.mTmult(system_matrix, assembling_matrix);
                system_matrix_inv.reinit(system_matrix.m(), system_matrix.m());
                system_matrix_inv.invert(system_matrix);
              }

              solution.reinit(system_matrix_inv.m(), 2);
              system_rhs.reinit(system_matrix_inv.m(), 2);
              tmp.reinit(2);

              for (unsigned int dof = 0; dof < this->n_dofs_per_face(face_no);
                   ++dof)
                {
                  system_rhs = 0.0;

                  for (unsigned int q_point = 0; q_point < n_quadrature_points;
                       ++q_point)
                    {
                      Point<dim> quadrature_point(
                        0.5 *
                          (quadrature_points[q_point](0) + shifts[subface][0]),
                        0.5 *
                          (quadrature_points[q_point](1) + shifts[subface][1]),
                        0.0);
                      tmp(0) = 0.5 * this->shape_value_component(
                                       this->face_to_cell_index(dof, 4),
                                       quadrature_point,
                                       0);
                      tmp(1) = 0.5 * this->shape_value_component(
                                       this->face_to_cell_index(dof, 4),
                                       quadrature_point,
                                       1);
                      quadrature_point =
                        Point<dim>(quadrature_points[q_point](0),
                                   quadrature_points[q_point](1),
                                   0.0);

                      for (unsigned int i = 0; i < 2; ++i)
                        for (unsigned int j = 0; j < source_fe.degree; ++j)
                          {
                            tmp(0) -= interpolation_matrix(
                                        (i + 2) * source_fe.degree + j, dof) *
                                      source_fe.shape_value_component(
                                        (i + 2) * source_fe.degree + j,
                                        quadrature_point,
                                        0);
                            tmp(1) -=
                              interpolation_matrix(i * source_fe.degree + j,
                                                   dof) *
                              source_fe.shape_value_component(
                                i * source_fe.degree + j, quadrature_point, 1);
                          }

                      tmp *= quadrature.weight(q_point);

                      for (unsigned int i = 0; i < source_fe.degree; ++i)
                        {
                          const double L_i_0 = legendre_polynomials[i].value(
                            quadrature_points[q_point](0));
                          const double L_i_1 = legendre_polynomials[i].value(
                            quadrature_points[q_point](1));

                          for (unsigned int j = 0; j < source_fe.degree - 1;
                               ++j)
                            {
                              system_rhs(i * (source_fe.degree - 1) + j, 0) +=
                                tmp(0) * L_i_0 *
                                lobatto_polynomials[j + 2].value(
                                  quadrature_points[q_point](1));
                              system_rhs(i * (source_fe.degree - 1) + j, 1) +=
                                tmp(1) * L_i_1 *
                                lobatto_polynomials[j + 2].value(
                                  quadrature_points[q_point](0));
                            }
                        }
                    }

                  system_matrix_inv.mmult(solution, system_rhs);

                  for (unsigned int i = 0; i < source_fe.degree; ++i)
                    for (unsigned int j = 0; j < source_fe.degree - 1; ++j)
                      {
                        if (std::abs(solution(i * (source_fe.degree - 1) + j,
                                              0)) > 1e-14)
                          interpolation_matrix(i * (source_fe.degree - 1) + j +
                                                 n_boundary_dofs,
                                               dof) =
                            solution(i * (source_fe.degree - 1) + j, 0);

                        if (std::abs(solution(i * (source_fe.degree - 1) + j,
                                              1)) > 1e-14)
                          interpolation_matrix(
                            i + (j + source_fe.degree - 1) * source_fe.degree +
                              n_boundary_dofs,
                            dof) = solution(i * (source_fe.degree - 1) + j, 1);
                      }
                }
            }

          break;
        }

      default:
        Assert(false, ExcNotImplemented());
    }
}

template <int dim>
const FullMatrix<double> &
FE_Nedelec<dim>::get_prolongation_matrix(
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
      FE_Nedelec<dim> &this_nonconst = const_cast<FE_Nedelec<dim> &>(*this);

      // Reinit the vectors of
      // restriction and prolongation
      // matrices to the right sizes.
      // Restriction only for isotropic
      // refinement
#ifdef DEBUG_NEDELEC
      deallog << "Embedding" << std::endl;
#endif
      this_nonconst.reinit_restriction_and_prolongation_matrices();
      // Fill prolongation matrices with embedding operators
      FETools::compute_embedding_matrices(
        this_nonconst,
        this_nonconst.prolongation,
        true,
        internal::FE_Nedelec::get_embedding_computation_tolerance(
          this->degree));
#ifdef DEBUG_NEDELEC
      deallog << "Restriction" << std::endl;
#endif
      this_nonconst.initialize_restriction();
    }

  // we use refinement_case-1 here. the -1 takes care of the origin of the
  // vector, as for RefinementCase<dim>::no_refinement (=0) there is no data
  // available and so the vector indices are shifted
  return this->prolongation[refinement_case - 1][child];
}

template <int dim>
const FullMatrix<double> &
FE_Nedelec<dim>::get_restriction_matrix(
  const unsigned int         child,
  const RefinementCase<dim> &refinement_case) const
{
  AssertIndexRange(refinement_case,
                   RefinementCase<dim>::isotropic_refinement + 1);
  Assert(refinement_case != RefinementCase<dim>::no_refinement,
         ExcMessage(
           "Restriction matrices are only available for refined cells!"));
  AssertIndexRange(
    child, GeometryInfo<dim>::n_children(RefinementCase<dim>(refinement_case)));

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
      FE_Nedelec<dim> &this_nonconst = const_cast<FE_Nedelec<dim> &>(*this);

      // Reinit the vectors of
      // restriction and prolongation
      // matrices to the right sizes.
      // Restriction only for isotropic
      // refinement
#ifdef DEBUG_NEDELEC
      deallog << "Embedding" << std::endl;
#endif
      this_nonconst.reinit_restriction_and_prolongation_matrices();
      // Fill prolongation matrices with embedding operators
      FETools::compute_embedding_matrices(
        this_nonconst,
        this_nonconst.prolongation,
        true,
        internal::FE_Nedelec::get_embedding_computation_tolerance(
          this->degree));
#ifdef DEBUG_NEDELEC
      deallog << "Restriction" << std::endl;
#endif
      this_nonconst.initialize_restriction();
    }

  // we use refinement_case-1 here. the -1 takes care of the origin of the
  // vector, as for RefinementCase<dim>::no_refinement (=0) there is no data
  // available and so the vector indices are shifted
  return this->restriction[refinement_case - 1][child];
}


// Interpolate a function, which is given by
// its values at the generalized support
// points in the finite element space on the
// reference cell.
// This is done as usual by projection-based
// interpolation.
template <int dim>
void
FE_Nedelec<dim>::convert_generalized_support_point_values_to_dof_values(
  const std::vector<Vector<double>> &support_point_values,
  std::vector<double> &              nodal_values) const
{
  // TODO: the implementation makes the assumption that all faces have the
  // same number of dofs
  AssertDimension(this->n_unique_faces(), 1);
  const unsigned int face_no = 0;

  const unsigned int deg = this->degree - 1;
  Assert(support_point_values.size() == this->generalized_support_points.size(),
         ExcDimensionMismatch(support_point_values.size(),
                              this->generalized_support_points.size()));
  Assert(support_point_values[0].size() == this->n_components(),
         ExcDimensionMismatch(support_point_values[0].size(),
                              this->n_components()));
  Assert(nodal_values.size() == this->n_dofs_per_cell(),
         ExcDimensionMismatch(nodal_values.size(), this->n_dofs_per_cell()));
  std::fill(nodal_values.begin(), nodal_values.end(), 0.0);

  switch (dim)
    {
      case 2:
        {
          // Let us begin with the
          // interpolation part.
          const QGauss<1>    reference_edge_quadrature(this->degree);
          const unsigned int n_edge_points = reference_edge_quadrature.size();

          for (unsigned int i = 0; i < 2; ++i)
            for (unsigned int j = 0; j < 2; ++j)
              {
                for (unsigned int q_point = 0; q_point < n_edge_points;
                     ++q_point)
                  nodal_values[(i + 2 * j) * this->degree] +=
                    reference_edge_quadrature.weight(q_point) *
                    support_point_values[q_point + (i + 2 * j) * n_edge_points]
                                        [1 - j];

                // Add the computed support_point_values to the resulting vector
                // only, if they are not
                // too small.
                if (std::abs(nodal_values[(i + 2 * j) * this->degree]) < 1e-14)
                  nodal_values[(i + 2 * j) * this->degree] = 0.0;
              }

          // If the Nedelec element degree is greater
          // than 0 (i.e., the polynomial degree is greater than 1),
          // then we have still some higher order edge
          // shape functions to consider.
          // Note that this->degree returns the polynomial
          // degree.
          // Here the projection part starts.
          // The dof support_point_values are obtained by solving
          // a linear system of
          // equations.
          if (this->degree > 1)
            {
              // We start with projection
              // on the higher order edge
              // shape function.
              const std::vector<Polynomials::Polynomial<double>>
                &lobatto_polynomials =
                  Polynomials::Lobatto::generate_complete_basis(this->degree);
              FullMatrix<double> system_matrix(this->degree - 1,
                                               this->degree - 1);
              std::vector<Polynomials::Polynomial<double>>
                lobatto_polynomials_grad(this->degree);

              for (unsigned int i = 0; i < lobatto_polynomials_grad.size(); ++i)
                lobatto_polynomials_grad[i] =
                  lobatto_polynomials[i + 1].derivative();

              // Set up the system matrix.
              // This can be used for all
              // edges.
              for (unsigned int i = 0; i < system_matrix.m(); ++i)
                for (unsigned int j = 0; j < system_matrix.n(); ++j)
                  for (unsigned int q_point = 0; q_point < n_edge_points;
                       ++q_point)
                    system_matrix(i, j) +=
                      boundary_weights(q_point, j) *
                      lobatto_polynomials_grad[i + 1].value(
                        this->generalized_face_support_points[face_no][q_point](
                          0));

              FullMatrix<double> system_matrix_inv(this->degree - 1,
                                                   this->degree - 1);

              system_matrix_inv.invert(system_matrix);

              const unsigned int
                             line_coordinate[GeometryInfo<2>::lines_per_cell] = {1, 1, 0, 0};
              Vector<double> system_rhs(system_matrix.m());
              Vector<double> solution(system_rhs.size());

              for (unsigned int line = 0;
                   line < GeometryInfo<dim>::lines_per_cell;
                   ++line)
                {
                  // Set up the right hand side.
                  system_rhs = 0;

                  for (unsigned int q_point = 0; q_point < n_edge_points;
                       ++q_point)
                    {
                      const double tmp =
                        support_point_values[line * n_edge_points + q_point]
                                            [line_coordinate[line]] -
                        nodal_values[line * this->degree] *
                          this->shape_value_component(
                            line * this->degree,
                            this->generalized_support_points[line *
                                                               n_edge_points +
                                                             q_point],
                            line_coordinate[line]);

                      for (unsigned int i = 0; i < system_rhs.size(); ++i)
                        system_rhs(i) += boundary_weights(q_point, i) * tmp;
                    }

                  system_matrix_inv.vmult(solution, system_rhs);

                  // Add the computed support_point_values
                  // to the resulting vector
                  // only, if they are not
                  // too small.
                  for (unsigned int i = 0; i < solution.size(); ++i)
                    if (std::abs(solution(i)) > 1e-14)
                      nodal_values[line * this->degree + i + 1] = solution(i);
                }

              // Then we go on to the
              // interior shape
              // functions. Again we
              // set up the system
              // matrix and use it
              // for both, the
              // horizontal and the
              // vertical, interior
              // shape functions.
              const QGauss<dim>  reference_quadrature(this->degree);
              const unsigned int n_interior_points =
                reference_quadrature.size();
              const std::vector<Polynomials::Polynomial<double>>
                &legendre_polynomials =
                  Polynomials::Legendre::generate_complete_basis(this->degree -
                                                                 1);

              system_matrix.reinit((this->degree - 1) * this->degree,
                                   (this->degree - 1) * this->degree);
              system_matrix = 0;

              for (unsigned int i = 0; i < this->degree; ++i)
                for (unsigned int j = 0; j < this->degree - 1; ++j)
                  for (unsigned int k = 0; k < this->degree; ++k)
                    for (unsigned int l = 0; l < this->degree - 1; ++l)
                      for (unsigned int q_point = 0;
                           q_point < n_interior_points;
                           ++q_point)
                        system_matrix(i * (this->degree - 1) + j,
                                      k * (this->degree - 1) + l) +=
                          reference_quadrature.weight(q_point) *
                          legendre_polynomials[i].value(
                            this->generalized_support_points
                              [q_point + GeometryInfo<dim>::lines_per_cell *
                                           n_edge_points](0)) *
                          lobatto_polynomials[j + 2].value(
                            this->generalized_support_points
                              [q_point + GeometryInfo<dim>::lines_per_cell *
                                           n_edge_points](1)) *
                          lobatto_polynomials_grad[k].value(
                            this->generalized_support_points
                              [q_point + GeometryInfo<dim>::lines_per_cell *
                                           n_edge_points](0)) *
                          lobatto_polynomials[l + 2].value(
                            this->generalized_support_points
                              [q_point + GeometryInfo<dim>::lines_per_cell *
                                           n_edge_points](1));

              system_matrix_inv.reinit(system_matrix.m(), system_matrix.m());
              system_matrix_inv.invert(system_matrix);
              // Set up the right hand side
              // for the horizontal shape
              // functions.
              system_rhs.reinit(system_matrix_inv.m());
              system_rhs = 0;

              for (unsigned int q_point = 0; q_point < n_interior_points;
                   ++q_point)
                {
                  double tmp =
                    support_point_values[q_point +
                                         GeometryInfo<dim>::lines_per_cell *
                                           n_edge_points][0];

                  for (unsigned int i = 0; i < 2; ++i)
                    for (unsigned int j = 0; j <= deg; ++j)
                      tmp -= nodal_values[(i + 2) * this->degree + j] *
                             this->shape_value_component(
                               (i + 2) * this->degree + j,
                               this->generalized_support_points
                                 [q_point + GeometryInfo<dim>::lines_per_cell *
                                              n_edge_points],
                               0);

                  for (unsigned int i = 0; i <= deg; ++i)
                    for (unsigned int j = 0; j < deg; ++j)
                      system_rhs(i * deg + j) +=
                        reference_quadrature.weight(q_point) * tmp *
                        lobatto_polynomials_grad[i].value(
                          this->generalized_support_points
                            [q_point + GeometryInfo<dim>::lines_per_cell *
                                         n_edge_points](0)) *
                        lobatto_polynomials[j + 2].value(
                          this->generalized_support_points
                            [q_point + GeometryInfo<dim>::lines_per_cell *
                                         n_edge_points](1));
                }

              solution.reinit(system_matrix.m());
              system_matrix_inv.vmult(solution, system_rhs);

              // Add the computed support_point_values
              // to the resulting vector
              // only, if they are not
              // too small.
              for (unsigned int i = 0; i <= deg; ++i)
                for (unsigned int j = 0; j < deg; ++j)
                  if (std::abs(solution(i * deg + j)) > 1e-14)
                    nodal_values[(i + GeometryInfo<dim>::lines_per_cell) * deg +
                                 j + GeometryInfo<dim>::lines_per_cell] =
                      solution(i * deg + j);

              system_rhs = 0;
              // Set up the right hand side
              // for the vertical shape
              // functions.

              for (unsigned int q_point = 0; q_point < n_interior_points;
                   ++q_point)
                {
                  double tmp =
                    support_point_values[q_point +
                                         GeometryInfo<dim>::lines_per_cell *
                                           n_edge_points][1];

                  for (unsigned int i = 0; i < 2; ++i)
                    for (unsigned int j = 0; j <= deg; ++j)
                      tmp -= nodal_values[i * this->degree + j] *
                             this->shape_value_component(
                               i * this->degree + j,
                               this->generalized_support_points
                                 [q_point + GeometryInfo<dim>::lines_per_cell *
                                              n_edge_points],
                               1);

                  for (unsigned int i = 0; i <= deg; ++i)
                    for (unsigned int j = 0; j < deg; ++j)
                      system_rhs(i * deg + j) +=
                        reference_quadrature.weight(q_point) * tmp *
                        lobatto_polynomials_grad[i].value(
                          this->generalized_support_points
                            [q_point + GeometryInfo<dim>::lines_per_cell *
                                         n_edge_points](1)) *
                        lobatto_polynomials[j + 2].value(
                          this->generalized_support_points
                            [q_point + GeometryInfo<dim>::lines_per_cell *
                                         n_edge_points](0));
                }

              system_matrix_inv.vmult(solution, system_rhs);

              // Add the computed support_point_values
              // to the resulting vector
              // only, if they are not
              // too small.
              for (unsigned int i = 0; i <= deg; ++i)
                for (unsigned int j = 0; j < deg; ++j)
                  if (std::abs(solution(i * deg + j)) > 1e-14)
                    nodal_values[i +
                                 (j + GeometryInfo<dim>::lines_per_cell + deg) *
                                   this->degree] = solution(i * deg + j);
            }

          break;
        }

      case 3:
        {
          // Let us begin with the
          // interpolation part.
          const QGauss<1>    reference_edge_quadrature(this->degree);
          const unsigned int n_edge_points = reference_edge_quadrature.size();

          for (unsigned int q_point = 0; q_point < n_edge_points; ++q_point)
            {
              for (unsigned int i = 0; i < 4; ++i)
                nodal_values[(i + 8) * this->degree] +=
                  reference_edge_quadrature.weight(q_point) *
                  support_point_values[q_point + (i + 8) * n_edge_points][2];

              for (unsigned int i = 0; i < 2; ++i)
                for (unsigned int j = 0; j < 2; ++j)
                  for (unsigned int k = 0; k < 2; ++k)
                    nodal_values[(i + 2 * (2 * j + k)) * this->degree] +=
                      reference_edge_quadrature.weight(q_point) *
                      support_point_values[q_point + (i + 2 * (2 * j + k)) *
                                                       n_edge_points][1 - k];
            }

          // Add the computed support_point_values
          // to the resulting vector
          // only, if they are not
          // too small.
          for (unsigned int i = 0; i < 4; ++i)
            if (std::abs(nodal_values[(i + 8) * this->degree]) < 1e-14)
              nodal_values[(i + 8) * this->degree] = 0.0;

          for (unsigned int i = 0; i < 2; ++i)
            for (unsigned int j = 0; j < 2; ++j)
              for (unsigned int k = 0; k < 2; ++k)
                if (std::abs(
                      nodal_values[(i + 2 * (2 * j + k)) * this->degree]) <
                    1e-14)
                  nodal_values[(i + 2 * (2 * j + k)) * this->degree] = 0.0;

          // If the degree is greater
          // than 0, then we have still
          // some higher order shape
          // functions to consider.
          // Here the projection part
          // starts. The dof support_point_values
          // are obtained by solving
          // a linear system of
          // equations.
          if (this->degree > 1)
            {
              // We start with projection
              // on the higher order edge
              // shape function.
              const std::vector<Polynomials::Polynomial<double>>
                &lobatto_polynomials =
                  Polynomials::Lobatto::generate_complete_basis(this->degree);
              FullMatrix<double> system_matrix(this->degree - 1,
                                               this->degree - 1);
              std::vector<Polynomials::Polynomial<double>>
                lobatto_polynomials_grad(this->degree);

              for (unsigned int i = 0; i < lobatto_polynomials_grad.size(); ++i)
                lobatto_polynomials_grad[i] =
                  lobatto_polynomials[i + 1].derivative();

              // Set up the system matrix.
              // This can be used for all
              // edges.
              for (unsigned int i = 0; i < system_matrix.m(); ++i)
                for (unsigned int j = 0; j < system_matrix.n(); ++j)
                  for (unsigned int q_point = 0; q_point < n_edge_points;
                       ++q_point)
                    system_matrix(i, j) +=
                      boundary_weights(q_point, j) *
                      lobatto_polynomials_grad[i + 1].value(
                        this->generalized_face_support_points[face_no][q_point](
                          1));

              FullMatrix<double> system_matrix_inv(this->degree - 1,
                                                   this->degree - 1);

              system_matrix_inv.invert(system_matrix);

              const unsigned int
                line_coordinate[GeometryInfo<3>::lines_per_cell] = {
                  1, 1, 0, 0, 1, 1, 0, 0, 2, 2, 2, 2};
              Vector<double> system_rhs(system_matrix.m());
              Vector<double> solution(system_rhs.size());

              for (unsigned int line = 0;
                   line < GeometryInfo<dim>::lines_per_cell;
                   ++line)
                {
                  // Set up the right hand side.
                  system_rhs = 0;

                  for (unsigned int q_point = 0; q_point < this->degree;
                       ++q_point)
                    {
                      const double tmp =
                        support_point_values[line * this->degree + q_point]
                                            [line_coordinate[line]] -
                        nodal_values[line * this->degree] *
                          this->shape_value_component(
                            line * this->degree,
                            this
                              ->generalized_support_points[line * this->degree +
                                                           q_point],
                            line_coordinate[line]);

                      for (unsigned int i = 0; i < system_rhs.size(); ++i)
                        system_rhs(i) += boundary_weights(q_point, i) * tmp;
                    }

                  system_matrix_inv.vmult(solution, system_rhs);

                  // Add the computed values
                  // to the resulting vector
                  // only, if they are not
                  // too small.
                  for (unsigned int i = 0; i < solution.size(); ++i)
                    if (std::abs(solution(i)) > 1e-14)
                      nodal_values[line * this->degree + i + 1] = solution(i);
                }

              // Then we go on to the
              // face shape functions.
              // Again we set up the
              // system matrix and
              // use it for both, the
              // horizontal and the
              // vertical, shape
              // functions.
              const std::vector<Polynomials::Polynomial<double>>
                &legendre_polynomials =
                  Polynomials::Legendre::generate_complete_basis(this->degree -
                                                                 1);
              const unsigned int n_face_points = n_edge_points * n_edge_points;

              system_matrix.reinit((this->degree - 1) * this->degree,
                                   (this->degree - 1) * this->degree);
              system_matrix = 0;

              for (unsigned int i = 0; i < this->degree; ++i)
                for (unsigned int j = 0; j < this->degree - 1; ++j)
                  for (unsigned int k = 0; k < this->degree; ++k)
                    for (unsigned int l = 0; l < this->degree - 1; ++l)
                      for (unsigned int q_point = 0; q_point < n_face_points;
                           ++q_point)
                        system_matrix(i * (this->degree - 1) + j,
                                      k * (this->degree - 1) + l) +=
                          boundary_weights(q_point + n_edge_points,
                                           2 * (k * (this->degree - 1) + l)) *
                          legendre_polynomials[i].value(
                            this->generalized_face_support_points
                              [face_no][q_point + 4 * n_edge_points](0)) *
                          lobatto_polynomials[j + 2].value(
                            this->generalized_face_support_points
                              [face_no][q_point + 4 * n_edge_points](1));

              system_matrix_inv.reinit(system_matrix.m(), system_matrix.m());
              system_matrix_inv.invert(system_matrix);
              solution.reinit(system_matrix.m());
              system_rhs.reinit(system_matrix.m());

              const unsigned int
                face_coordinates[GeometryInfo<3>::faces_per_cell][2] = {
                  {1, 2}, {1, 2}, {2, 0}, {2, 0}, {0, 1}, {0, 1}};
              const unsigned int edge_indices[GeometryInfo<3>::faces_per_cell]
                                             [GeometryInfo<3>::lines_per_face] =
                                               {{0, 4, 8, 10},
                                                {1, 5, 9, 11},
                                                {8, 9, 2, 6},
                                                {10, 11, 3, 7},
                                                {2, 3, 0, 1},
                                                {6, 7, 4, 5}};

              for (const unsigned int face : GeometryInfo<dim>::face_indices())
                {
                  // Set up the right hand side
                  // for the horizontal shape
                  // functions.
                  system_rhs = 0;

                  for (unsigned int q_point = 0; q_point < n_face_points;
                       ++q_point)
                    {
                      double tmp =
                        support_point_values[q_point +
                                             GeometryInfo<dim>::lines_per_cell *
                                               n_edge_points]
                                            [face_coordinates[face][0]];

                      for (unsigned int i = 0; i < 2; ++i)
                        for (unsigned int j = 0; j <= deg; ++j)
                          tmp -=
                            nodal_values[edge_indices[face][i] * this->degree +
                                         j] *
                            this->shape_value_component(
                              edge_indices[face][i] * this->degree + j,
                              this->generalized_support_points
                                [q_point + GeometryInfo<dim>::lines_per_cell *
                                             n_edge_points],
                              face_coordinates[face][0]);

                      for (unsigned int i = 0; i <= deg; ++i)
                        for (unsigned int j = 0; j < deg; ++j)
                          system_rhs(i * deg + j) +=
                            boundary_weights(q_point + n_edge_points,
                                             2 * (i * deg + j)) *
                            tmp;
                    }

                  system_matrix_inv.vmult(solution, system_rhs);

                  // Add the computed support_point_values
                  // to the resulting vector
                  // only, if they are not
                  // too small.
                  for (unsigned int i = 0; i <= deg; ++i)
                    for (unsigned int j = 0; j < deg; ++j)
                      if (std::abs(solution(i * deg + j)) > 1e-14)
                        nodal_values[(2 * face * this->degree + i +
                                      GeometryInfo<dim>::lines_per_cell) *
                                       deg +
                                     j + GeometryInfo<dim>::lines_per_cell] =
                          solution(i * deg + j);

                  // Set up the right hand side
                  // for the vertical shape
                  // functions.
                  system_rhs = 0;

                  for (unsigned int q_point = 0; q_point < n_face_points;
                       ++q_point)
                    {
                      double tmp =
                        support_point_values[q_point +
                                             GeometryInfo<dim>::lines_per_cell *
                                               n_edge_points]
                                            [face_coordinates[face][1]];

                      for (unsigned int i = 2;
                           i < GeometryInfo<dim>::lines_per_face;
                           ++i)
                        for (unsigned int j = 0; j <= deg; ++j)
                          tmp -=
                            nodal_values[edge_indices[face][i] * this->degree +
                                         j] *
                            this->shape_value_component(
                              edge_indices[face][i] * this->degree + j,
                              this->generalized_support_points
                                [q_point + GeometryInfo<dim>::lines_per_cell *
                                             n_edge_points],
                              face_coordinates[face][1]);

                      for (unsigned int i = 0; i <= deg; ++i)
                        for (unsigned int j = 0; j < deg; ++j)
                          system_rhs(i * deg + j) +=
                            boundary_weights(q_point + n_edge_points,
                                             2 * (i * deg + j) + 1) *
                            tmp;
                    }

                  system_matrix_inv.vmult(solution, system_rhs);

                  // Add the computed support_point_values
                  // to the resulting vector
                  // only, if they are not
                  // too small.
                  for (unsigned int i = 0; i <= deg; ++i)
                    for (unsigned int j = 0; j < deg; ++j)
                      if (std::abs(solution(i * deg + j)) > 1e-14)
                        nodal_values[((2 * face + 1) * deg + j +
                                      GeometryInfo<dim>::lines_per_cell) *
                                       this->degree +
                                     i] = solution(i * deg + j);
                }

              // Finally we project
              // the remaining parts
              // of the function on
              // the interior shape
              // functions.
              const QGauss<dim>  reference_quadrature(this->degree);
              const unsigned int n_interior_points =
                reference_quadrature.size();

              // We create the
              // system matrix.
              system_matrix.reinit(this->degree * deg * deg,
                                   this->degree * deg * deg);
              system_matrix = 0;

              for (unsigned int i = 0; i <= deg; ++i)
                for (unsigned int j = 0; j < deg; ++j)
                  for (unsigned int k = 0; k < deg; ++k)
                    for (unsigned int l = 0; l <= deg; ++l)
                      for (unsigned int m = 0; m < deg; ++m)
                        for (unsigned int n = 0; n < deg; ++n)
                          for (unsigned int q_point = 0;
                               q_point < n_interior_points;
                               ++q_point)
                            system_matrix((i * deg + j) * deg + k,
                                          (l * deg + m) * deg + n) +=
                              reference_quadrature.weight(q_point) *
                              legendre_polynomials[i].value(
                                this->generalized_support_points
                                  [q_point +
                                   GeometryInfo<dim>::lines_per_cell *
                                     n_edge_points +
                                   GeometryInfo<dim>::faces_per_cell *
                                     n_face_points](0)) *
                              lobatto_polynomials[j + 2].value(
                                this->generalized_support_points
                                  [q_point +
                                   GeometryInfo<dim>::lines_per_cell *
                                     n_edge_points +
                                   GeometryInfo<dim>::faces_per_cell *
                                     n_face_points](1)) *
                              lobatto_polynomials[k + 2].value(
                                this->generalized_support_points
                                  [q_point +
                                   GeometryInfo<dim>::lines_per_cell *
                                     n_edge_points +
                                   GeometryInfo<dim>::faces_per_cell *
                                     n_face_points](2)) *
                              lobatto_polynomials_grad[l].value(
                                this->generalized_support_points
                                  [q_point +
                                   GeometryInfo<dim>::lines_per_cell *
                                     n_edge_points +
                                   GeometryInfo<dim>::faces_per_cell *
                                     n_face_points](0)) *
                              lobatto_polynomials[m + 2].value(
                                this->generalized_support_points
                                  [q_point +
                                   GeometryInfo<dim>::lines_per_cell *
                                     n_edge_points +
                                   GeometryInfo<dim>::faces_per_cell *
                                     n_face_points](1)) *
                              lobatto_polynomials[n + 2].value(
                                this->generalized_support_points
                                  [q_point +
                                   GeometryInfo<dim>::lines_per_cell *
                                     n_edge_points +
                                   GeometryInfo<dim>::faces_per_cell *
                                     n_face_points](2));

              system_matrix_inv.reinit(system_matrix.m(), system_matrix.m());
              system_matrix_inv.invert(system_matrix);
              // Set up the right hand side.
              system_rhs.reinit(system_matrix.m());
              system_rhs = 0;

              for (unsigned int q_point = 0; q_point < n_interior_points;
                   ++q_point)
                {
                  double tmp =
                    support_point_values[q_point +
                                         GeometryInfo<dim>::lines_per_cell *
                                           n_edge_points +
                                         GeometryInfo<dim>::faces_per_cell *
                                           n_face_points][0];

                  for (unsigned int i = 0; i <= deg; ++i)
                    {
                      for (unsigned int j = 0; j < 2; ++j)
                        for (unsigned int k = 0; k < 2; ++k)
                          tmp -=
                            nodal_values[i + (j + 4 * k + 2) * this->degree] *
                            this->shape_value_component(
                              i + (j + 4 * k + 2) * this->degree,
                              this->generalized_support_points
                                [q_point +
                                 GeometryInfo<dim>::lines_per_cell *
                                   n_edge_points +
                                 GeometryInfo<dim>::faces_per_cell *
                                   n_face_points],
                              0);

                      for (unsigned int j = 0; j < deg; ++j)
                        for (unsigned int k = 0; k < 4; ++k)
                          tmp -=
                            nodal_values[(i + 2 * (k + 2) * this->degree +
                                          GeometryInfo<dim>::lines_per_cell) *
                                           deg +
                                         j +
                                         GeometryInfo<dim>::lines_per_cell] *
                            this->shape_value_component(
                              (i + 2 * (k + 2) * this->degree +
                               GeometryInfo<dim>::lines_per_cell) *
                                  deg +
                                j + GeometryInfo<dim>::lines_per_cell,
                              this->generalized_support_points
                                [q_point +
                                 GeometryInfo<dim>::lines_per_cell *
                                   n_edge_points +
                                 GeometryInfo<dim>::faces_per_cell *
                                   n_face_points],
                              0);
                    }

                  for (unsigned int i = 0; i <= deg; ++i)
                    for (unsigned int j = 0; j < deg; ++j)
                      for (unsigned int k = 0; k < deg; ++k)
                        system_rhs((i * deg + j) * deg + k) +=
                          reference_quadrature.weight(q_point) * tmp *
                          lobatto_polynomials_grad[i].value(
                            this->generalized_support_points
                              [q_point +
                               GeometryInfo<dim>::lines_per_cell *
                                 n_edge_points +
                               GeometryInfo<dim>::faces_per_cell *
                                 n_face_points](0)) *
                          lobatto_polynomials[j + 2].value(
                            this->generalized_support_points
                              [q_point +
                               GeometryInfo<dim>::lines_per_cell *
                                 n_edge_points +
                               GeometryInfo<dim>::faces_per_cell *
                                 n_face_points](1)) *
                          lobatto_polynomials[k + 2].value(
                            this->generalized_support_points
                              [q_point +
                               GeometryInfo<dim>::lines_per_cell *
                                 n_edge_points +
                               GeometryInfo<dim>::faces_per_cell *
                                 n_face_points](2));
                }

              solution.reinit(system_rhs.size());
              system_matrix_inv.vmult(solution, system_rhs);

              // Add the computed values
              // to the resulting vector
              // only, if they are not
              // too small.
              for (unsigned int i = 0; i <= deg; ++i)
                for (unsigned int j = 0; j < deg; ++j)
                  for (unsigned int k = 0; k < deg; ++k)
                    if (std::abs(solution((i * deg + j) * deg + k)) > 1e-14)
                      nodal_values
                        [((i + 2 * GeometryInfo<dim>::faces_per_cell) * deg +
                          j + GeometryInfo<dim>::lines_per_cell +
                          2 * GeometryInfo<dim>::faces_per_cell) *
                           deg +
                         k + GeometryInfo<dim>::lines_per_cell] =
                          solution((i * deg + j) * deg + k);

              // Set up the right hand side.
              system_rhs = 0;

              for (unsigned int q_point = 0; q_point < n_interior_points;
                   ++q_point)
                {
                  double tmp =
                    support_point_values[q_point +
                                         GeometryInfo<dim>::lines_per_cell *
                                           n_edge_points +
                                         GeometryInfo<dim>::faces_per_cell *
                                           n_face_points][1];

                  for (unsigned int i = 0; i <= deg; ++i)
                    for (unsigned int j = 0; j < 2; ++j)
                      {
                        for (unsigned int k = 0; k < 2; ++k)
                          tmp -= nodal_values[i + (4 * j + k) * this->degree] *
                                 this->shape_value_component(
                                   i + (4 * j + k) * this->degree,
                                   this->generalized_support_points
                                     [q_point +
                                      GeometryInfo<dim>::lines_per_cell *
                                        n_edge_points +
                                      GeometryInfo<dim>::faces_per_cell *
                                        n_face_points],
                                   1);

                        for (unsigned int k = 0; k < deg; ++k)
                          tmp -=
                            nodal_values[(i + 2 * j * this->degree +
                                          GeometryInfo<dim>::lines_per_cell) *
                                           deg +
                                         k +
                                         GeometryInfo<dim>::lines_per_cell] *
                              this->shape_value_component(
                                (i + 2 * j * this->degree +
                                 GeometryInfo<dim>::lines_per_cell) *
                                    deg +
                                  k + GeometryInfo<dim>::lines_per_cell,
                                this->generalized_support_points
                                  [q_point +
                                   GeometryInfo<dim>::lines_per_cell *
                                     n_edge_points +
                                   GeometryInfo<dim>::faces_per_cell *
                                     n_face_points],
                                1) +
                            nodal_values[i +
                                         ((2 * j + 9) * deg + k +
                                          GeometryInfo<dim>::lines_per_cell) *
                                           this->degree] *
                              this->shape_value_component(
                                i + ((2 * j + 9) * deg + k +
                                     GeometryInfo<dim>::lines_per_cell) *
                                      this->degree,
                                this->generalized_support_points
                                  [q_point +
                                   GeometryInfo<dim>::lines_per_cell *
                                     n_edge_points +
                                   GeometryInfo<dim>::faces_per_cell *
                                     n_face_points],
                                1);
                      }

                  for (unsigned int i = 0; i <= deg; ++i)
                    for (unsigned int j = 0; j < deg; ++j)
                      for (unsigned int k = 0; k < deg; ++k)
                        system_rhs((i * deg + j) * deg + k) +=
                          reference_quadrature.weight(q_point) * tmp *
                          lobatto_polynomials_grad[i].value(
                            this->generalized_support_points
                              [q_point +
                               GeometryInfo<dim>::lines_per_cell *
                                 n_edge_points +
                               GeometryInfo<dim>::faces_per_cell *
                                 n_face_points](1)) *
                          lobatto_polynomials[j + 2].value(
                            this->generalized_support_points
                              [q_point +
                               GeometryInfo<dim>::lines_per_cell *
                                 n_edge_points +
                               GeometryInfo<dim>::faces_per_cell *
                                 n_face_points](0)) *
                          lobatto_polynomials[k + 2].value(
                            this->generalized_support_points
                              [q_point +
                               GeometryInfo<dim>::lines_per_cell *
                                 n_edge_points +
                               GeometryInfo<dim>::faces_per_cell *
                                 n_face_points](2));
                }

              system_matrix_inv.vmult(solution, system_rhs);

              // Add the computed support_point_values
              // to the resulting vector
              // only, if they are not
              // too small.
              for (unsigned int i = 0; i <= deg; ++i)
                for (unsigned int j = 0; j < deg; ++j)
                  for (unsigned int k = 0; k < deg; ++k)
                    if (std::abs(solution((i * deg + j) * deg + k)) > 1e-14)
                      nodal_values[((i + this->degree +
                                     2 * GeometryInfo<dim>::faces_per_cell) *
                                      deg +
                                    j + GeometryInfo<dim>::lines_per_cell +
                                    2 * GeometryInfo<dim>::faces_per_cell) *
                                     deg +
                                   k + GeometryInfo<dim>::lines_per_cell] =
                        solution((i * deg + j) * deg + k);

              // Set up the right hand side.
              system_rhs = 0;

              for (unsigned int q_point = 0; q_point < n_interior_points;
                   ++q_point)
                {
                  double tmp =
                    support_point_values[q_point +
                                         GeometryInfo<dim>::lines_per_cell *
                                           n_edge_points +
                                         GeometryInfo<dim>::faces_per_cell *
                                           n_face_points][2];

                  for (unsigned int i = 0; i <= deg; ++i)
                    for (unsigned int j = 0; j < 4; ++j)
                      {
                        tmp -= nodal_values[i + (j + 8) * this->degree] *
                               this->shape_value_component(
                                 i + (j + 8) * this->degree,
                                 this->generalized_support_points
                                   [q_point +
                                    GeometryInfo<dim>::lines_per_cell *
                                      n_edge_points +
                                    GeometryInfo<dim>::faces_per_cell *
                                      n_face_points],
                                 2);

                        for (unsigned int k = 0; k < deg; ++k)
                          tmp -=
                            nodal_values[i +
                                         ((2 * j + 1) * deg + k +
                                          GeometryInfo<dim>::lines_per_cell) *
                                           this->degree] *
                            this->shape_value_component(
                              i + ((2 * j + 1) * deg + k +
                                   GeometryInfo<dim>::lines_per_cell) *
                                    this->degree,
                              this->generalized_support_points
                                [q_point +
                                 GeometryInfo<dim>::lines_per_cell *
                                   n_edge_points +
                                 GeometryInfo<dim>::faces_per_cell *
                                   n_face_points],
                              2);
                      }

                  for (unsigned int i = 0; i <= deg; ++i)
                    for (unsigned int j = 0; j < deg; ++j)
                      for (unsigned int k = 0; k < deg; ++k)
                        system_rhs((i * deg + j) * deg + k) +=
                          reference_quadrature.weight(q_point) * tmp *
                          lobatto_polynomials_grad[i].value(
                            this->generalized_support_points
                              [q_point +
                               GeometryInfo<dim>::lines_per_cell *
                                 n_edge_points +
                               GeometryInfo<dim>::faces_per_cell *
                                 n_face_points](2)) *
                          lobatto_polynomials[j + 2].value(
                            this->generalized_support_points
                              [q_point +
                               GeometryInfo<dim>::lines_per_cell *
                                 n_edge_points +
                               GeometryInfo<dim>::faces_per_cell *
                                 n_face_points](0)) *
                          lobatto_polynomials[k + 2].value(
                            this->generalized_support_points
                              [q_point +
                               GeometryInfo<dim>::lines_per_cell *
                                 n_edge_points +
                               GeometryInfo<dim>::faces_per_cell *
                                 n_face_points](1));
                }

              system_matrix_inv.vmult(solution, system_rhs);

              // Add the computed support_point_values
              // to the resulting vector
              // only, if they are not
              // too small.
              for (unsigned int i = 0; i <= deg; ++i)
                for (unsigned int j = 0; j < deg; ++j)
                  for (unsigned int k = 0; k < deg; ++k)
                    if (std::abs(solution((i * deg + j) * deg + k)) > 1e-14)
                      nodal_values
                        [i +
                         ((j + 2 * (deg + GeometryInfo<dim>::faces_per_cell)) *
                            deg +
                          k + GeometryInfo<dim>::lines_per_cell) *
                           this->degree] = solution((i * deg + j) * deg + k);
            }

          break;
        }

      default:
        Assert(false, ExcNotImplemented());
    }
}



template <int dim>
std::pair<Table<2, bool>, std::vector<unsigned int>>
FE_Nedelec<dim>::get_constant_modes() const
{
  Table<2, bool> constant_modes(dim, this->n_dofs_per_cell());
  for (unsigned int d = 0; d < dim; ++d)
    for (unsigned int i = 0; i < this->n_dofs_per_cell(); ++i)
      constant_modes(d, i) = true;
  std::vector<unsigned int> components;
  for (unsigned int d = 0; d < dim; ++d)
    components.push_back(d);
  return std::pair<Table<2, bool>, std::vector<unsigned int>>(constant_modes,
                                                              components);
}


template <int dim>
std::size_t
FE_Nedelec<dim>::memory_consumption() const
{
  Assert(false, ExcNotImplemented());
  return 0;
}


//----------------------------------------------------------------------//


// explicit instantiations
#include "fe_nedelec.inst"


DEAL_II_NAMESPACE_CLOSE
