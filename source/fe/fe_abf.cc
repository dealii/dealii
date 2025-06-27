// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2003 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


#include <deal.II/base/qprojector.h>
#include <deal.II/base/quadrature.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/table.h>
#include <deal.II/base/utilities.h>

#include <deal.II/dofs/dof_accessor.h>

#include <deal.II/fe/fe.h>
#include <deal.II/fe/fe_abf.h>
#include <deal.II/fe/fe_tools.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>

#include <iostream>
#include <memory>
#include <sstream>

DEAL_II_NAMESPACE_OPEN

// TODO: implement the adjust_quad_dof_index_for_face_orientation_table and
// adjust_line_dof_index_for_line_orientation_table fields, and write tests
// similar to bits/face_orientation_and_fe_q_*

template <int dim>
FE_ABF<dim>::FE_ABF(const unsigned int deg)
  : FE_PolyTensor<dim>(
      PolynomialsABF<dim>(deg),
      FiniteElementData<dim>(get_dpo_vector(deg),
                             dim,
                             deg + 2,
                             FiniteElementData<dim>::Hdiv),
      std::vector<bool>(PolynomialsABF<dim>::n_polynomials(deg), true),
      std::vector<ComponentMask>(PolynomialsABF<dim>::n_polynomials(deg),
                                 ComponentMask(std::vector<bool>(dim, true))))
  , rt_order(deg)
{
  Assert(dim >= 2, ExcImpossibleInDim(dim));
  const unsigned int n_dofs = this->n_dofs_per_cell();

  this->mapping_kind = {mapping_raviart_thomas};
  // First, initialize the
  // generalized support points and
  // quadrature weights, since they
  // are required for interpolation.
  initialize_support_points(deg);

  // Now compute the inverse node matrix, generating the correct
  // basis functions from the raw ones. For a discussion of what
  // exactly happens here, see FETools::compute_node_matrix.
  const FullMatrix<double> M = FETools::compute_node_matrix(*this);
  this->inverse_node_matrix.reinit(n_dofs, n_dofs);
  this->inverse_node_matrix.invert(M);
  // From now on, the shape functions provided by FiniteElement::shape_value
  // and similar functions will be the correct ones, not
  // the raw shape functions from the polynomial space anymore.

  // Reinit the vectors of
  // restriction and prolongation
  // matrices to the right sizes.
  // Restriction only for isotropic
  // refinement
  this->reinit_restriction_and_prolongation_matrices(true);
  // Fill prolongation matrices with embedding operators
  FETools::compute_embedding_matrices(*this, this->prolongation, false, 1.e-10);

  initialize_restriction();

  // TODO: the implementation makes the assumption that all faces have the
  // same number of dofs
  AssertDimension(this->n_unique_faces(), 1);
  const unsigned int face_no = 0;

  // TODO[TL]: for anisotropic refinement we will probably need a table of
  // submatrices with an array for each refine case
  std::vector<FullMatrix<double>> face_embeddings(
    1 << (dim - 1),
    FullMatrix<double>(this->n_dofs_per_face(face_no),
                       this->n_dofs_per_face(face_no)));
  // TODO: Something goes wrong there. The error of the least squares fit
  // is to large ...
  // FETools::compute_face_embedding_matrices(*this, face_embeddings.data(), 0,
  // 0);
  this->interface_constraints.reinit((1 << (dim - 1)) *
                                       this->n_dofs_per_face(face_no),
                                     this->n_dofs_per_face(face_no));
  unsigned int target_row = 0;
  for (const auto &face_embedding : face_embeddings)
    for (unsigned int i = 0; i < face_embedding.m(); ++i)
      {
        for (unsigned int j = 0; j < face_embedding.n(); ++j)
          this->interface_constraints(target_row, j) = face_embedding(i, j);
        ++target_row;
      }

  // We need to initialize the dof permutation table and the one for the sign
  // change.
  initialize_quad_dof_index_permutation_and_sign_change();
}


template <int dim>
void
FE_ABF<dim>::initialize_quad_dof_index_permutation_and_sign_change()
{
  // for 1d and 2d, do nothing
  if (dim < 3)
    return;

  // TODO: Implement this for this class
  return;
}


template <int dim>
std::string
FE_ABF<dim>::get_name() const
{
  // note that the
  // FETools::get_fe_by_name
  // function depends on the
  // particular format of the string
  // this function returns, so they
  // have to be kept in synch

  std::ostringstream namebuf;

  namebuf << "FE_ABF<" << dim << ">(" << rt_order << ")";

  return namebuf.str();
}



template <int dim>
std::unique_ptr<FiniteElement<dim, dim>>
FE_ABF<dim>::clone() const
{
  return std::make_unique<FE_ABF<dim>>(rt_order);
}


//---------------------------------------------------------------------------
// Auxiliary and internal functions
//---------------------------------------------------------------------------



// Version for 2d and higher. See above for 1d version
template <int dim>
void
FE_ABF<dim>::initialize_support_points(const unsigned int deg)
{
  const QGauss<dim>  cell_quadrature(deg + 2);
  const unsigned int n_interior_points = cell_quadrature.size();

  // TODO: the implementation makes the assumption that all faces have the
  // same number of dofs
  AssertDimension(this->n_unique_faces(), 1);
  const unsigned int face_no = 0;

  unsigned int n_face_points = (dim > 1) ? 1 : 0;
  // compute (deg+1)^(dim-1)
  for (unsigned int d = 1; d < dim; ++d)
    n_face_points *= deg + 1;

  this->generalized_support_points.resize(
    GeometryInfo<dim>::faces_per_cell * n_face_points + n_interior_points);
  this->generalized_face_support_points[face_no].resize(n_face_points);


  // These might be required when the faces contribution is computed
  // Therefore they will be initialized at this point.
  std::array<std::unique_ptr<AnisotropicPolynomials<dim>>, dim> polynomials_abf;

  // Generate x_1^{i} x_2^{r+1} ...
  for (unsigned int dd = 0; dd < dim; ++dd)
    {
      std::vector<std::vector<Polynomials::Polynomial<double>>> poly(dim);
      for (unsigned int d = 0; d < dim; ++d)
        poly[d].push_back(Polynomials::Monomial<double>(deg + 1));
      poly[dd] = Polynomials::Monomial<double>::generate_complete_basis(deg);

      polynomials_abf[dd] = std::make_unique<AnisotropicPolynomials<dim>>(poly);
    }

  // Number of the point being entered
  unsigned int current = 0;

  if (dim > 1)
    {
      const QGauss<dim - 1>             face_points(deg + 1);
      TensorProductPolynomials<dim - 1> legendre =
        Polynomials::Legendre::generate_complete_basis(deg);

      boundary_weights.reinit(n_face_points, legendre.n());

      //       Assert (face_points.size() == this->n_dofs_per_face(),
      //            ExcInternalError());

      for (unsigned int k = 0; k < n_face_points; ++k)
        {
          this->generalized_face_support_points[face_no][k] =
            face_points.point(k);
          // Compute its quadrature
          // contribution for each
          // moment.
          for (unsigned int i = 0; i < legendre.n(); ++i)
            {
              boundary_weights(k, i) =
                face_points.weight(k) *
                legendre.compute_value(i, face_points.point(k));
            }
        }

      Quadrature<dim> faces =
        QProjector<dim>::project_to_all_faces(this->reference_cell(),
                                              face_points);
      for (unsigned int face_no = 0;
           face_no < GeometryInfo<dim>::faces_per_cell;
           ++face_no)
        {
          const auto offset = QProjector<dim>::DataSetDescriptor::face(
            this->reference_cell(),
            face_no,
            numbers::default_geometric_orientation,
            n_face_points);
          for (unsigned int face_point = 0; face_point < n_face_points;
               ++face_point)
            {
              // Enter the support point into the vector
              this->generalized_support_points[current] =
                faces.point(offset + face_point);
              ++current;
            }
        }


      // Now initialize edge interior weights for the ABF elements.
      // These are completely independent from the usual edge moments. They
      // stem from applying the Gauss theorem to the nodal values, which
      // was necessary to cast the ABF elements into the deal.II framework
      // for vector valued elements.
      boundary_weights_abf.reinit(faces.size(), polynomials_abf[0]->n() * dim);
      for (unsigned int k = 0; k < faces.size(); ++k)
        {
          for (unsigned int i = 0; i < polynomials_abf[0]->n() * dim; ++i)
            {
              boundary_weights_abf(k, i) =
                polynomials_abf[i % dim]->compute_value(i / dim,
                                                        faces.point(k)) *
                faces.weight(k);
            }
        }
    }

  // Create Legendre basis for the
  // space D_xi Q_k
  if (deg > 0)
    {
      std::array<std::unique_ptr<AnisotropicPolynomials<dim>>, dim> polynomials;

      for (unsigned int dd = 0; dd < dim; ++dd)
        {
          std::vector<std::vector<Polynomials::Polynomial<double>>> poly(dim);
          for (unsigned int d = 0; d < dim; ++d)
            poly[d] = Polynomials::Legendre::generate_complete_basis(deg);
          poly[dd] = Polynomials::Legendre::generate_complete_basis(deg - 1);

          polynomials[dd] = std::make_unique<AnisotropicPolynomials<dim>>(poly);
        }

      interior_weights.reinit(
        TableIndices<3>(n_interior_points, polynomials[0]->n(), dim));

      for (unsigned int k = 0; k < cell_quadrature.size(); ++k)
        {
          for (unsigned int i = 0; i < polynomials[0]->n(); ++i)
            for (unsigned int d = 0; d < dim; ++d)
              interior_weights(k, i, d) =
                cell_quadrature.weight(k) *
                polynomials[d]->compute_value(i, cell_quadrature.point(k));
        }
    }


  // Decouple the creation of the generalized support points
  // from computation of interior weights.
  for (unsigned int k = 0; k < cell_quadrature.size(); ++k)
    this->generalized_support_points[current++] = cell_quadrature.point(k);

  // Additional functionality for the ABF elements
  // TODO: Here the canonical extension of the principle
  // behind the ABF elements is implemented. It is unclear,
  // if this really leads to the ABF spaces in 3d!
  interior_weights_abf.reinit(TableIndices<3>(cell_quadrature.size(),
                                              polynomials_abf[0]->n() * dim,
                                              dim));
  Tensor<1, dim> poly_grad;

  for (unsigned int k = 0; k < cell_quadrature.size(); ++k)
    {
      for (unsigned int i = 0; i < polynomials_abf[0]->n() * dim; ++i)
        {
          poly_grad =
            polynomials_abf[i % dim]->compute_grad(i / dim,
                                                   cell_quadrature.point(k)) *
            cell_quadrature.weight(k);
          // The minus sign comes from the use of the Gauss theorem to replace
          // the divergence.
          for (unsigned int d = 0; d < dim; ++d)
            interior_weights_abf(k, i, d) = -poly_grad[d];
        }
    }

  Assert(current == this->generalized_support_points.size(),
         ExcInternalError());
}



// This function is the same Raviart-Thomas interpolation performed by
// interpolate. Still, we cannot use interpolate, since it was written
// for smooth functions. The functions interpolated here are not
// smooth, maybe even not continuous. Therefore, we must double the
// number of quadrature points in each direction in order to integrate
// only smooth functions.

// Then again, the interpolated function is chosen such that the
// moments coincide with the function to be interpolated.

template <int dim>
void
FE_ABF<dim>::initialize_restriction()
{
  if (dim == 1)
    {
      unsigned int iso = RefinementCase<dim>::isotropic_refinement - 1;
      for (unsigned int i = 0; i < GeometryInfo<dim>::max_children_per_cell;
           ++i)
        this->restriction[iso][i].reinit(0, 0);
      return;
    }
  unsigned int          iso = RefinementCase<dim>::isotropic_refinement - 1;
  const QGauss<dim - 1> q_base(rt_order + 1);
  const unsigned int    n_face_points = q_base.size();
  // First, compute interpolation on
  // subfaces
  for (const unsigned int face : GeometryInfo<dim>::face_indices())
    {
      // The shape functions of the
      // child cell are evaluated
      // in the quadrature points
      // of a full face.
      Quadrature<dim> q_face = QProjector<dim>::project_to_face(
        this->reference_cell(),
        q_base,
        face,
        numbers::default_geometric_orientation);
      // Store shape values, since the
      // evaluation suffers if not
      // ordered by point
      Table<2, double> cached_values_face(this->n_dofs_per_cell(),
                                          q_face.size());
      for (unsigned int k = 0; k < q_face.size(); ++k)
        for (unsigned int i = 0; i < this->n_dofs_per_cell(); ++i)
          cached_values_face(i, k) = this->shape_value_component(
            i, q_face.point(k), GeometryInfo<dim>::unit_normal_direction[face]);

      for (unsigned int sub = 0; sub < GeometryInfo<dim>::max_children_per_face;
           ++sub)
        {
          // The weight functions for
          // the coarse face are
          // evaluated on the subface
          // only.
          Quadrature<dim> q_sub = QProjector<dim>::project_to_subface(
            this->reference_cell(),
            q_base,
            face,
            sub,
            numbers::default_geometric_orientation,
            RefinementCase<dim - 1>::isotropic_refinement);
          const unsigned int child = GeometryInfo<dim>::child_cell_on_face(
            RefinementCase<dim>::isotropic_refinement, face, sub);

          // On a certain face, we must
          // compute the moments of ALL
          // fine level functions with
          // the coarse level weight
          // functions belonging to
          // that face. Due to the
          // orthogonalization process
          // when building the shape
          // functions, these weights
          // are equal to the
          // corresponding shape
          // functions.
          for (unsigned int k = 0; k < n_face_points; ++k)
            for (unsigned int i_child = 0; i_child < this->n_dofs_per_cell();
                 ++i_child)
              for (unsigned int i_face = 0;
                   i_face < this->n_dofs_per_face(face);
                   ++i_face)
                {
                  // The quadrature
                  // weights on the
                  // subcell are NOT
                  // transformed, so we
                  // have to do it here.
                  this->restriction[iso][child](
                    face * this->n_dofs_per_face(face) + i_face, i_child) +=
                    Utilities::fixed_power<dim - 1>(.5) * q_sub.weight(k) *
                    cached_values_face(i_child, k) *
                    this->shape_value_component(
                      face * this->n_dofs_per_face(face) + i_face,
                      q_sub.point(k),
                      GeometryInfo<dim>::unit_normal_direction[face]);
                }
        }
    }

  if (rt_order == 0)
    return;

  // Create Legendre basis for the
  // space D_xi Q_k. Here, we cannot
  // use the shape functions
  std::array<std::unique_ptr<AnisotropicPolynomials<dim>>, dim> polynomials;
  for (unsigned int dd = 0; dd < dim; ++dd)
    {
      std::vector<std::vector<Polynomials::Polynomial<double>>> poly(dim);
      for (unsigned int d = 0; d < dim; ++d)
        poly[d] = Polynomials::Legendre::generate_complete_basis(rt_order);
      poly[dd] = Polynomials::Legendre::generate_complete_basis(rt_order - 1);

      polynomials[dd] = std::make_unique<AnisotropicPolynomials<dim>>(poly);
    }

  // TODO: the implementation makes the assumption that all faces have the
  // same number of dofs
  AssertDimension(this->n_unique_faces(), 1);
  const unsigned int face_no = 0;

  const QGauss<dim>  q_cell(rt_order + 1);
  const unsigned int start_cell_dofs =
    GeometryInfo<dim>::faces_per_cell * this->n_dofs_per_face(face_no);

  // Store shape values, since the
  // evaluation suffers if not
  // ordered by point
  Table<3, double> cached_values_cell(this->n_dofs_per_cell(),
                                      q_cell.size(),
                                      dim);
  for (unsigned int k = 0; k < q_cell.size(); ++k)
    for (unsigned int i = 0; i < this->n_dofs_per_cell(); ++i)
      for (unsigned int d = 0; d < dim; ++d)
        cached_values_cell(i, k, d) =
          this->shape_value_component(i, q_cell.point(k), d);

  for (unsigned int child = 0; child < GeometryInfo<dim>::max_children_per_cell;
       ++child)
    {
      Quadrature<dim> q_sub =
        QProjector<dim>::project_to_child(this->reference_cell(),
                                          q_cell,
                                          child);

      for (unsigned int k = 0; k < q_sub.size(); ++k)
        for (unsigned int i_child = 0; i_child < this->n_dofs_per_cell();
             ++i_child)
          for (unsigned int d = 0; d < dim; ++d)
            for (unsigned int i_weight = 0; i_weight < polynomials[d]->n();
                 ++i_weight)
              {
                this->restriction[iso][child](start_cell_dofs + i_weight * dim +
                                                d,
                                              i_child) +=
                  q_sub.weight(k) * cached_values_cell(i_child, k, d) *
                  polynomials[d]->compute_value(i_weight, q_sub.point(k));
              }
    }
}



template <int dim>
std::vector<unsigned int>
FE_ABF<dim>::get_dpo_vector(const unsigned int rt_order)
{
  if (dim == 1)
    {
      Assert(false, ExcImpossibleInDim(1));
      return std::vector<unsigned int>();
    }

  // the element is face-based (not
  // to be confused with George
  // W. Bush's Faith Based
  // Initiative...), and we have
  // (rt_order+1)^(dim-1) DoFs per face
  unsigned int dofs_per_face = 1;
  for (unsigned int d = 0; d < dim - 1; ++d)
    dofs_per_face *= rt_order + 1;

  // and then there are interior dofs
  const unsigned int interior_dofs = dim * (rt_order + 1) * dofs_per_face;

  std::vector<unsigned int> dpo(dim + 1);
  dpo[dim - 1] = dofs_per_face;
  dpo[dim]     = interior_dofs;

  return dpo;
}

//---------------------------------------------------------------------------
// Data field initialization
//---------------------------------------------------------------------------

template <int dim>
bool
FE_ABF<dim>::has_support_on_face(const unsigned int shape_index,
                                 const unsigned int face_index) const
{
  AssertIndexRange(shape_index, this->n_dofs_per_cell());
  AssertIndexRange(face_index, GeometryInfo<dim>::faces_per_cell);

  // Return computed values if we
  // know them easily. Otherwise, it
  // is always safe to return true.
  switch (rt_order)
    {
      case 0:
        {
          switch (dim)
            {
              case 2:
                {
                  // only on the one
                  // non-adjacent face
                  // are the values
                  // actually zero. list
                  // these in a table
                  return (face_index !=
                          GeometryInfo<dim>::opposite_face[shape_index]);
                }

              default:
                return true;
            }
        }

      default: // other rt_order
        return true;
    }

  return true;
}



template <int dim>
void
FE_ABF<dim>::convert_generalized_support_point_values_to_dof_values(
  const std::vector<Vector<double>> &support_point_values,
  std::vector<double>               &nodal_values) const
{
  Assert(support_point_values.size() == this->generalized_support_points.size(),
         ExcDimensionMismatch(support_point_values.size(),
                              this->generalized_support_points.size()));
  Assert(support_point_values[0].size() == this->n_components(),
         ExcDimensionMismatch(support_point_values[0].size(),
                              this->n_components()));
  Assert(nodal_values.size() == this->n_dofs_per_cell(),
         ExcDimensionMismatch(nodal_values.size(), this->n_dofs_per_cell()));

  std::fill(nodal_values.begin(), nodal_values.end(), 0.);

  const unsigned int n_face_points = boundary_weights.size(0);
  for (const unsigned int face : GeometryInfo<dim>::face_indices())
    for (unsigned int k = 0; k < n_face_points; ++k)
      for (unsigned int i = 0; i < boundary_weights.size(1); ++i)
        {
          nodal_values[i + face * this->n_dofs_per_face(face)] +=
            boundary_weights(k, i) *
            support_point_values[face * n_face_points + k][GeometryInfo<
              dim>::unit_normal_direction[face]];
        }

  // TODO: the implementation makes the assumption that all faces have the
  // same number of dofs
  AssertDimension(this->n_unique_faces(), 1);
  const unsigned int face_no = 0;

  const unsigned int start_cell_dofs =
    GeometryInfo<dim>::faces_per_cell * this->n_dofs_per_face(face_no);
  const unsigned int start_cell_points =
    GeometryInfo<dim>::faces_per_cell * n_face_points;

  for (unsigned int k = 0; k < interior_weights.size(0); ++k)
    for (unsigned int i = 0; i < interior_weights.size(1); ++i)
      for (unsigned int d = 0; d < dim; ++d)
        nodal_values[start_cell_dofs + i * dim + d] +=
          interior_weights(k, i, d) *
          support_point_values[k + start_cell_points][d];

  const unsigned int start_abf_dofs =
    start_cell_dofs + interior_weights.size(1) * dim;

  // Cell integral of ABF terms
  for (unsigned int k = 0; k < interior_weights_abf.size(0); ++k)
    for (unsigned int i = 0; i < interior_weights_abf.size(1); ++i)
      for (unsigned int d = 0; d < dim; ++d)
        nodal_values[start_abf_dofs + i] +=
          interior_weights_abf(k, i, d) *
          support_point_values[k + start_cell_points][d];

  // Face integral of ABF terms
  for (const unsigned int face : GeometryInfo<dim>::face_indices())
    {
      const double n_orient = GeometryInfo<dim>::unit_normal_orientation[face];
      for (unsigned int fp = 0; fp < n_face_points; ++fp)
        {
          // TODO: Check what the face_orientation, face_flip and face_rotation
          // have to be in 3d
          unsigned int k = QProjector<dim>::DataSetDescriptor::face(
            this->reference_cell(),
            face,
            numbers::default_geometric_orientation,
            n_face_points);
          for (unsigned int i = 0; i < boundary_weights_abf.size(1); ++i)
            nodal_values[start_abf_dofs + i] +=
              n_orient * boundary_weights_abf(k + fp, i) *
              support_point_values[face * n_face_points + fp][GeometryInfo<
                dim>::unit_normal_direction[face]];
        }
    }

  // TODO: Check if this "correction" can be removed.
  for (unsigned int i = 0; i < boundary_weights_abf.size(1); ++i)
    if (std::fabs(nodal_values[start_abf_dofs + i]) < 1.0e-16)
      nodal_values[start_abf_dofs + i] = 0.0;
}



template <int dim>
std::size_t
FE_ABF<dim>::memory_consumption() const
{
  DEAL_II_NOT_IMPLEMENTED();
  return 0;
}



/*-------------- Explicit Instantiations -------------------------------*/
#include "fe/fe_abf.inst"

DEAL_II_NAMESPACE_CLOSE
