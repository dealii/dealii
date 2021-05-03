// ---------------------------------------------------------------------
//
// Copyright (C) 2003 - 2020 by the deal.II authors
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


#include <deal.II/base/qprojector.h>
#include <deal.II/base/quadrature.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/table.h>
#include <deal.II/base/utilities.h>

#include <deal.II/dofs/dof_accessor.h>

#include <deal.II/fe/fe.h>
#include <deal.II/fe/fe_raviart_thomas.h>
#include <deal.II/fe/fe_tools.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>

#include <iostream>
#include <memory>
#include <sstream>


DEAL_II_NAMESPACE_OPEN


template <int dim>
FE_RaviartThomas<dim>::FE_RaviartThomas(const unsigned int deg)
  : FE_PolyTensor<dim>(
      PolynomialsRaviartThomas<dim>(deg),
      FiniteElementData<dim>(get_dpo_vector(deg),
                             dim,
                             deg + 1,
                             FiniteElementData<dim>::Hdiv),
      std::vector<bool>(PolynomialsRaviartThomas<dim>::n_polynomials(deg),
                        true),
      std::vector<ComponentMask>(PolynomialsRaviartThomas<dim>::n_polynomials(
                                   deg),
                                 std::vector<bool>(dim, true)))
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
  FETools::compute_embedding_matrices(*this, this->prolongation);
  initialize_restriction();

  // TODO: the implementation makes the assumption that all faces have the
  // same number of dofs
  AssertDimension(this->n_unique_faces(), 1);
  const unsigned int face_no = 0;

  // TODO[TL]: for anisotropic refinement we will probably need a table of
  // submatrices with an array for each refine case
  FullMatrix<double> face_embeddings[GeometryInfo<dim>::max_children_per_face];
  for (unsigned int i = 0; i < GeometryInfo<dim>::max_children_per_face; ++i)
    face_embeddings[i].reinit(this->n_dofs_per_face(face_no),
                              this->n_dofs_per_face(face_no));
  FETools::compute_face_embedding_matrices<dim, double>(*this,
                                                        face_embeddings,
                                                        0,
                                                        0);
  this->interface_constraints.reinit((1 << (dim - 1)) *
                                       this->n_dofs_per_face(face_no),
                                     this->n_dofs_per_face(face_no));
  unsigned int target_row = 0;
  for (unsigned int d = 0; d < GeometryInfo<dim>::max_children_per_face; ++d)
    for (unsigned int i = 0; i < face_embeddings[d].m(); ++i)
      {
        for (unsigned int j = 0; j < face_embeddings[d].n(); ++j)
          this->interface_constraints(target_row, j) = face_embeddings[d](i, j);
        ++target_row;
      }

  // We need to initialize the dof permuation table and the one for the sign
  // change.
  initialize_quad_dof_index_permutation_and_sign_change();
}



template <int dim>
std::string
FE_RaviartThomas<dim>::get_name() const
{
  // note that the
  // FETools::get_fe_by_name
  // function depends on the
  // particular format of the string
  // this function returns, so they
  // have to be kept in synch

  // note that this->degree is the maximal
  // polynomial degree and is thus one higher
  // than the argument given to the
  // constructor
  std::ostringstream namebuf;
  namebuf << "FE_RaviartThomas<" << dim << ">(" << this->degree - 1 << ")";

  return namebuf.str();
}


template <int dim>
std::unique_ptr<FiniteElement<dim, dim>>
FE_RaviartThomas<dim>::clone() const
{
  return std::make_unique<FE_RaviartThomas<dim>>(*this);
}


//---------------------------------------------------------------------------
// Auxiliary and internal functions
//---------------------------------------------------------------------------


template <int dim>
void
FE_RaviartThomas<dim>::initialize_support_points(const unsigned int deg)
{
  QGauss<dim>        cell_quadrature(deg + 1);
  const unsigned int n_interior_points = (deg > 0) ? cell_quadrature.size() : 0;

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

  // Number of the point being entered
  unsigned int current = 0;

  if (dim > 1)
    {
      QGauss<dim - 1>                   face_points(deg + 1);
      TensorProductPolynomials<dim - 1> legendre =
        Polynomials::Legendre::generate_complete_basis(deg);

      boundary_weights.reinit(n_face_points, legendre.n());

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
      for (; current < GeometryInfo<dim>::faces_per_cell * n_face_points;
           ++current)
        {
          // Enter the support point
          // into the vector
          this->generalized_support_points[current] = faces.point(
            current +
            QProjector<dim>::DataSetDescriptor::face(
              this->reference_cell(), 0, true, false, false, n_face_points));
        }
    }

  if (deg == 0)
    return;

  // Create Legendre basis for the space D_xi Q_k
  std::unique_ptr<AnisotropicPolynomials<dim>> polynomials[dim];
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
      this->generalized_support_points[current++] = cell_quadrature.point(k);
      for (unsigned int i = 0; i < polynomials[0]->n(); ++i)
        for (unsigned int d = 0; d < dim; ++d)
          interior_weights(k, i, d) =
            cell_quadrature.weight(k) *
            polynomials[d]->compute_value(i, cell_quadrature.point(k));
    }

  Assert(current == this->generalized_support_points.size(),
         ExcInternalError());
}


template <int dim>
void
FE_RaviartThomas<dim>::initialize_quad_dof_index_permutation_and_sign_change()
{
  // For 1D do nothing.
  //
  // TODO: For 2D we simply keep the legacy behavior for now. This should be
  // changed in the future and can be taken care of by similar means as the 3D
  // case below. The legacy behavior can be found in fe_poly_tensor.cc in the
  // function internal::FE_PolyTensor::get_dof_sign_change_h_div(...)
  if (dim < 3)
    return;

  // TODO: the implementation makes the assumption that all faces have the
  // same number of dofs
  AssertDimension(this->n_unique_faces(), 1);
  const unsigned int face_no = 0;

  Assert(this->adjust_quad_dof_index_for_face_orientation_table[0]
             .n_elements() == 8 * this->n_dofs_per_quad(face_no),
         ExcInternalError());

  // The 3D RaviartThomas space has tensor_degree*tensor_degree face dofs
  const unsigned int n = this->tensor_degree();
  Assert(n * n == this->n_dofs_per_quad(face_no), ExcInternalError());

  // The vector of tables adjust_quad_dof_index_for_face_orientation_table
  // contains offsets for local face_dofs and is being filled with zeros in
  // fe.cc. We need to fill it with the correct values in case of non-standard,
  // flipped (rotated by +180 degrees) or rotated (rotated by +90 degrees) faces

  // The dofs on a face are connected to a n x n matrix. for example, for
  // tensor_degree==3 we have the following dofs on a quad:

  //  ___________
  // |           |
  // |  6  7  8  |
  // |           |
  // |  3  4  5  |
  // |           |
  // |  0  1  2  |
  // |___________|
  //
  // We have dof_index=i+n*j with index i in x-direction and index j in
  // y-direction running from 0 to n-1.  to extract i and j we can use
  // i=dof_index%n and j=dof_index/n. The indices i and j can then be used to
  // compute the offset.

  // Example: if the switches are (true | true | true) that means we rotate the
  // face first by + 90 degree(counterclockwise) then by another +180
  // degrees but we do not flip it since the face has standard
  // orientation. The flip axis is the diagonal from the lower left to the upper
  // right corner of the face. With these flags the configuration above becomes:
  //  ___________
  // |           |
  // |  2  5  8  |
  // |           |
  // |  1  4  7  |
  // |           |
  // |  0  3  6  |
  // |___________|
  //
  // This is exactly what is realized by the formulas implemented below. Note
  // that the necessity of a permuattion depends on the three flags.

  // There is also a pattern for the sign change of the mapped face_dofs
  // depending on the switches. In the above example it would be
  //  ___________
  // |           |
  // |  +  -  +  |
  // |           |
  // |  +  -  +  |
  // |           |
  // |  +  -  +  |
  // |___________|
  //

  for (unsigned int local_face_dof = 0;
       local_face_dof < this->n_dofs_per_quad(face_no);
       ++local_face_dof)
    {
      // Row and column
      unsigned int i = local_face_dof % n;
      unsigned int j = local_face_dof / n;

      // We have 8 cases that are all treated the same way. Note that the
      // corresponding case to case_no is just its binary representation.
      // The above example of (false | true | true) would be case_no=3
      for (unsigned int case_no = 0; case_no < 8; ++case_no)
        {
          // Get the binary representation of the case
          const bool face_orientation = Utilities::get_bit(case_no, 2);
          const bool face_flip        = Utilities::get_bit(case_no, 1);
          const bool face_rotation    = Utilities::get_bit(case_no, 0);

          if (((!face_orientation) && (!face_rotation)) ||
              ((face_orientation) && (face_rotation)))
            {
              // We flip across the diagonal
              // This is the local face dof offset
              this->adjust_quad_dof_index_for_face_orientation_table[face_no](
                local_face_dof, case_no) = j + i * n - local_face_dof;
            }
          else
            {
              // Offset is zero
              this->adjust_quad_dof_index_for_face_orientation_table[face_no](
                local_face_dof, case_no) = 0;
            } // if face needs dof permutation

          // Get new local face_dof by adding offset
          const unsigned int new_local_face_dof =
            local_face_dof +
            this->adjust_quad_dof_index_for_face_orientation_table[face_no](
              local_face_dof, case_no);
          // compute new row and column index
          i = new_local_face_dof % n;
          j = new_local_face_dof / n;

          /*
           * Now compute if a sign change is necessary. This is done for the
           * case of face_orientation==true
           */
          // flip = false, rotation=true
          if (!face_flip && face_rotation)
            {
              this->adjust_quad_dof_sign_for_face_orientation_table[face_no](
                local_face_dof, case_no) = ((j % 2) == 1);
            }
          // flip = true, rotation=false
          else if (face_flip && !face_rotation)
            {
              // This case is symmetric (although row and column may be
              // switched)
              this->adjust_quad_dof_sign_for_face_orientation_table[face_no](
                local_face_dof, case_no) = ((j % 2) == 1) != ((i % 2) == 1);
            }
          // flip = true, rotation=true
          else if (face_flip && face_rotation)
            {
              this->adjust_quad_dof_sign_for_face_orientation_table[face_no](
                local_face_dof, case_no) = ((i % 2) == 1);
            }
          /*
           * flip = false, rotation=false => nothing to do
           */

          /*
           * If face_orientation==false the sign flip is exactly the opposite.
           */
          if (!face_orientation)
            this->adjust_quad_dof_sign_for_face_orientation_table[face_no](
              local_face_dof, case_no) =
              !this->adjust_quad_dof_sign_for_face_orientation_table[face_no](
                local_face_dof, case_no);
        } // case_no
    }     // local_face_dof
}


template <>
void
FE_RaviartThomas<1>::initialize_restriction()
{
  // there is only one refinement case in 1d,
  // which is the isotropic one (first index of
  // the matrix array has to be 0)
  for (unsigned int i = 0; i < GeometryInfo<1>::max_children_per_cell; ++i)
    this->restriction[0][i].reinit(0, 0);
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
FE_RaviartThomas<dim>::initialize_restriction()
{
  const unsigned int iso = RefinementCase<dim>::isotropic_refinement - 1;

  QGauss<dim - 1>    q_base(this->degree);
  const unsigned int n_face_points = q_base.size();
  // First, compute interpolation on
  // subfaces
  for (unsigned int face : GeometryInfo<dim>::face_indices())
    {
      // The shape functions of the
      // child cell are evaluated
      // in the quadrature points
      // of a full face.
      Quadrature<dim> q_face =
        QProjector<dim>::project_to_face(this->reference_cell(), q_base, face);
      // Store shape values, since the
      // evaluation suffers if not
      // ordered by point
      Table<2, double> cached_values_on_face(this->n_dofs_per_cell(),
                                             q_face.size());
      for (unsigned int k = 0; k < q_face.size(); ++k)
        for (unsigned int i = 0; i < this->n_dofs_per_cell(); ++i)
          cached_values_on_face(i, k) = this->shape_value_component(
            i, q_face.point(k), GeometryInfo<dim>::unit_normal_direction[face]);

      for (unsigned int sub = 0; sub < GeometryInfo<dim>::max_children_per_face;
           ++sub)
        {
          // The weight functions for
          // the coarse face are
          // evaluated on the subface
          // only.
          Quadrature<dim> q_sub = QProjector<dim>::project_to_subface(
            this->reference_cell(), q_base, face, sub);
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
                    cached_values_on_face(i_child, k) *
                    this->shape_value_component(
                      face * this->n_dofs_per_face(face) + i_face,
                      q_sub.point(k),
                      GeometryInfo<dim>::unit_normal_direction[face]);
                }
        }
    }

  if (this->degree == 1)
    return;

  // Create Legendre basis for the space D_xi Q_k. Here, we cannot
  // use the shape functions
  std::unique_ptr<AnisotropicPolynomials<dim>> polynomials[dim];
  for (unsigned int dd = 0; dd < dim; ++dd)
    {
      std::vector<std::vector<Polynomials::Polynomial<double>>> poly(dim);
      for (unsigned int d = 0; d < dim; ++d)
        poly[d] =
          Polynomials::Legendre::generate_complete_basis(this->degree - 1);
      poly[dd] =
        Polynomials::Legendre::generate_complete_basis(this->degree - 2);

      polynomials[dd] = std::make_unique<AnisotropicPolynomials<dim>>(poly);
    }

  // TODO: the implementation makes the assumption that all faces have the
  // same number of dofs
  AssertDimension(this->n_unique_faces(), 1);
  const unsigned int face_no = 0;

  QGauss<dim>        q_cell(this->degree);
  const unsigned int start_cell_dofs =
    GeometryInfo<dim>::faces_per_cell * this->n_dofs_per_face(face_no);

  // Store shape values, since the
  // evaluation suffers if not
  // ordered by point
  Table<3, double> cached_values_on_cell(this->n_dofs_per_cell(),
                                         q_cell.size(),
                                         dim);
  for (unsigned int k = 0; k < q_cell.size(); ++k)
    for (unsigned int i = 0; i < this->n_dofs_per_cell(); ++i)
      for (unsigned int d = 0; d < dim; ++d)
        cached_values_on_cell(i, k, d) =
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
                  q_sub.weight(k) * cached_values_on_cell(i_child, k, d) *
                  polynomials[d]->compute_value(i_weight, q_sub.point(k));
              }
    }
}



template <int dim>
std::vector<unsigned int>
FE_RaviartThomas<dim>::get_dpo_vector(const unsigned int deg)
{
  // the element is face-based and we have
  // (deg+1)^(dim-1) DoFs per face
  unsigned int dofs_per_face = 1;
  for (unsigned int d = 1; d < dim; ++d)
    dofs_per_face *= deg + 1;

  // and then there are interior dofs
  const unsigned int interior_dofs = dim * deg * dofs_per_face;

  std::vector<unsigned int> dpo(dim + 1);
  dpo[dim - 1] = dofs_per_face;
  dpo[dim]     = interior_dofs;

  return dpo;
}



template <int dim>
std::pair<Table<2, bool>, std::vector<unsigned int>>
FE_RaviartThomas<dim>::get_constant_modes() const
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



//---------------------------------------------------------------------------
// Data field initialization
//---------------------------------------------------------------------------


template <int dim>
bool
FE_RaviartThomas<dim>::has_support_on_face(const unsigned int shape_index,
                                           const unsigned int face_index) const
{
  AssertIndexRange(shape_index, this->n_dofs_per_cell());
  AssertIndexRange(face_index, GeometryInfo<dim>::faces_per_cell);

  // Return computed values if we
  // know them easily. Otherwise, it
  // is always safe to return true.
  switch (this->degree)
    {
      case 1:
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
FE_RaviartThomas<dim>::convert_generalized_support_point_values_to_dof_values(
  const std::vector<Vector<double>> &support_point_values,
  std::vector<double> &              nodal_values) const
{
  Assert(support_point_values.size() == this->generalized_support_points.size(),
         ExcDimensionMismatch(support_point_values.size(),
                              this->generalized_support_points.size()));
  Assert(nodal_values.size() == this->n_dofs_per_cell(),
         ExcDimensionMismatch(nodal_values.size(), this->n_dofs_per_cell()));
  Assert(support_point_values[0].size() == this->n_components(),
         ExcDimensionMismatch(support_point_values[0].size(),
                              this->n_components()));

  std::fill(nodal_values.begin(), nodal_values.end(), 0.);

  const unsigned int n_face_points = boundary_weights.size(0);
  for (unsigned int face : GeometryInfo<dim>::face_indices())
    for (unsigned int k = 0; k < n_face_points; ++k)
      for (unsigned int i = 0; i < boundary_weights.size(1); ++i)
        {
          nodal_values[i + face * this->n_dofs_per_face(face)] +=
            boundary_weights(k, i) *
            support_point_values[face * n_face_points + k](
              GeometryInfo<dim>::unit_normal_direction[face]);
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
          support_point_values[k + start_cell_points](d);
}



template <int dim>
std::size_t
FE_RaviartThomas<dim>::memory_consumption() const
{
  Assert(false, ExcNotImplemented());
  return 0;
}



// explicit instantiations
#include "fe_raviart_thomas.inst"


DEAL_II_NAMESPACE_CLOSE
