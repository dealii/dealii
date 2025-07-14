// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2018 - 2025 by the deal.II authors
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
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_accessor.h>

#include <deal.II/fe/fe.h>
#include <deal.II/fe/fe_rt_bubbles.h>
#include <deal.II/fe/fe_tools.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>

#include <memory>
#include <sstream>


DEAL_II_NAMESPACE_OPEN

// TODO: implement the adjust_quad_dof_index_for_face_orientation_table and
// adjust_line_dof_index_for_line_orientation_table fields, and write tests
// similar to bits/face_orientation_and_fe_q_*

template <int dim>
FE_RT_Bubbles<dim>::FE_RT_Bubbles(const unsigned int deg)
  : FE_PolyTensor<dim>(
      PolynomialsRT_Bubbles<dim>(deg),
      FiniteElementData<dim>(get_dpo_vector(deg),
                             dim,
                             deg + 1,
                             FiniteElementData<dim>::Hdiv),
      get_ria_vector(deg),
      std::vector<ComponentMask>(PolynomialsRT_Bubbles<dim>::n_polynomials(deg),
                                 ComponentMask(std::vector<bool>(dim, true))))
{
  Assert(dim >= 2, ExcImpossibleInDim(dim));
  Assert(
    deg >= 1,
    ExcMessage(
      "Lowest order RT_Bubbles element is degree 1, but you requested for degree 0"));
  const unsigned int n_dofs = this->n_dofs_per_cell();

  this->mapping_kind = {mapping_raviart_thomas};
  // Initialize support points and quadrature weights
  initialize_support_points(deg);
  // Compute the inverse node matrix to get
  // the correct basis functions
  FullMatrix<double> M = FETools::compute_node_matrix(*this);
  this->inverse_node_matrix.reinit(n_dofs, n_dofs);
  this->inverse_node_matrix.invert(M);

  // Reinit the vectors of prolongation matrices to the
  // right sizes. There are no restriction matrices implemented
  for (const unsigned int ref_case :
       RefinementCase<dim>::all_refinement_cases())
    if (ref_case != RefinementCase<dim>::no_refinement)
      {
        const unsigned int nc = this->reference_cell().template n_children<dim>(
          RefinementCase<dim>(ref_case));

        for (unsigned int i = 0; i < nc; ++i)
          this->prolongation[ref_case - 1][i].reinit(n_dofs, n_dofs);
      }

  // TODO: the implementation makes the assumption that all faces have the
  // same number of dofs
  AssertDimension(this->n_unique_faces(), 1);
  const unsigned int face_no = 0;

  // Fill prolongation matrices with embedding operators
  // set tolerance to 1, as embedding error accumulate quickly
  FETools::compute_embedding_matrices(*this, this->prolongation, true, 1.0);
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

  // We need to initialize the dof permutation table and the one for the sign
  // change.
  initialize_quad_dof_index_permutation_and_sign_change();
}


template <int dim>
void
FE_RT_Bubbles<dim>::initialize_quad_dof_index_permutation_and_sign_change()
{
  // for 1d and 2d, do nothing
  if (dim < 3)
    return;

  // TODO: Implement this for this class
  return;
}



template <int dim>
std::string
FE_RT_Bubbles<dim>::get_name() const
{
  // Note: this->degree is the maximal polynomial degree and is thus one higher
  // than the argument given to the constructor
  std::ostringstream namebuf;
  namebuf << "FE_RT_Bubbles<" << dim << ">(" << this->degree << ")";

  return namebuf.str();
}



template <int dim>
std::unique_ptr<FiniteElement<dim, dim>>
FE_RT_Bubbles<dim>::clone() const
{
  return std::make_unique<FE_RT_Bubbles<dim>>(*this);
}


//---------------------------------------------------------------------------
// Auxiliary and internal functions
//---------------------------------------------------------------------------



template <int dim>
void
FE_RT_Bubbles<dim>::initialize_support_points(const unsigned int deg)
{
  // TODO: the implementation makes the assumption that all faces have the
  // same number of dofs
  AssertDimension(this->n_unique_faces(), 1);
  const unsigned int face_no = 0;

  this->generalized_support_points.resize(this->n_dofs_per_cell());
  this->generalized_face_support_points[face_no].resize(
    this->n_dofs_per_face(face_no));

  // Index of the point being entered
  unsigned int current = 0;

  // On the faces, we choose as many Gauss-Lobatto points
  // as required to determine the normal component uniquely.
  // This is the deg of the RT_Bubble element plus one.
  if (dim > 1)
    {
      const QGaussLobatto<dim - 1> face_points(deg + 1);
      Assert(face_points.size() == this->n_dofs_per_face(face_no),
             ExcInternalError());
      for (unsigned int k = 0; k < this->n_dofs_per_face(face_no); ++k)
        this->generalized_face_support_points[face_no][k] =
          face_points.point(k);
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
            face_points.size());
          for (unsigned int face_point = 0; face_point < face_points.size();
               ++face_point)
            {
              // Enter the support point into the vector
              this->generalized_support_points[current] =
                faces.point(offset + face_point);
              ++current;
            }
        }
    }

  if (deg == 1)
    return;

  // In the interior, we need anisotropic Gauss-Lobatto quadratures,
  // one for each direction
  const QGaussLobatto<1> high(deg + 1);
  std::vector<Point<1>>  pts = high.get_points();
  if (pts.size() > 2)
    {
      pts.erase(pts.begin());
      pts.erase(pts.end() - 1);
    }

  std::vector<double> wts(pts.size(), 1);
  const Quadrature<1> low(pts, wts);

  for (unsigned int d = 0; d < dim; ++d)
    {
      std::unique_ptr<QAnisotropic<dim>> quadrature;
      switch (dim)
        {
          case 1:
            quadrature = std::make_unique<QAnisotropic<dim>>(high);
            break;
          case 2:
            quadrature =
              std::make_unique<QAnisotropic<dim>>(((d == 0) ? low : high),
                                                  ((d == 1) ? low : high));
            break;
          case 3:
            quadrature =
              std::make_unique<QAnisotropic<dim>>(((d == 0) ? low : high),
                                                  ((d == 1) ? low : high),
                                                  ((d == 2) ? low : high));
            break;
          default:
            DEAL_II_NOT_IMPLEMENTED();
        }

      for (unsigned int k = 0; k < quadrature->size(); ++k)
        this->generalized_support_points[current++] = quadrature->point(k);
    }
  Assert(current == this->n_dofs_per_cell(), ExcInternalError());
}



template <int dim>
std::vector<unsigned int>
FE_RT_Bubbles<dim>::get_dpo_vector(const unsigned int deg)
{
  // We have (deg+1)^(dim-1) DoFs per face...
  unsigned int dofs_per_face = 1;
  for (unsigned int d = 1; d < dim; ++d)
    dofs_per_face *= deg + 1;

  // ...plus the interior DoFs for the total of dim*(deg+1)^dim
  const unsigned int interior_dofs =
    dim * (deg - 1) * Utilities::pow(deg + 1, dim - 1);

  std::vector<unsigned int> dpo(dim + 1);
  dpo[dim - 1] = dofs_per_face;
  dpo[dim]     = interior_dofs;

  return dpo;
}



template <>
std::vector<bool>
FE_RT_Bubbles<1>::get_ria_vector(const unsigned int)
{
  Assert(false, ExcImpossibleInDim(1));
  return std::vector<bool>();
}



template <int dim>
std::vector<bool>
FE_RT_Bubbles<dim>::get_ria_vector(const unsigned int deg)
{
  const unsigned int dofs_per_cell =
    PolynomialsRT_Bubbles<dim>::n_polynomials(deg);
  unsigned int dofs_per_face = deg + 1;
  for (unsigned int d = 2; d < dim; ++d)
    dofs_per_face *= deg + 1;
  // All face dofs need to be non-additive, since they have
  // continuity requirements. The interior dofs are
  // made additive.
  std::vector<bool> ret_val(dofs_per_cell, false);
  for (unsigned int i = GeometryInfo<dim>::faces_per_cell * dofs_per_face;
       i < dofs_per_cell;
       ++i)
    ret_val[i] = true;

  return ret_val;
}



template <int dim>
void
FE_RT_Bubbles<dim>::convert_generalized_support_point_values_to_dof_values(
  const std::vector<Vector<double>> &support_point_values,
  std::vector<double>               &nodal_values) const
{
  Assert(support_point_values.size() == this->generalized_support_points.size(),
         ExcDimensionMismatch(support_point_values.size(),
                              this->generalized_support_points.size()));
  Assert(nodal_values.size() == this->n_dofs_per_cell(),
         ExcDimensionMismatch(nodal_values.size(), this->n_dofs_per_cell()));
  Assert(support_point_values[0].size() == this->n_components(),
         ExcDimensionMismatch(support_point_values[0].size(),
                              this->n_components()));

  // First do interpolation on faces. There, the component
  // evaluated depends on the face direction and orientation.
  unsigned int fbase = 0;
  unsigned int f     = 0;
  for (; f < GeometryInfo<dim>::faces_per_cell;
       ++f, fbase += this->n_dofs_per_face(f))
    {
      for (unsigned int i = 0; i < this->n_dofs_per_face(f); ++i)
        {
          nodal_values[fbase + i] = support_point_values[fbase + i](
            GeometryInfo<dim>::unit_normal_direction[f]);
        }
    }

  // The remaining points form dim chunks, one for each component.
  const unsigned int istep = (this->n_dofs_per_cell() - fbase) / dim;
  Assert((this->n_dofs_per_cell() - fbase) % dim == 0, ExcInternalError());

  f = 0;
  while (fbase < this->n_dofs_per_cell())
    {
      for (unsigned int i = 0; i < istep; ++i)
        {
          nodal_values[fbase + i] = support_point_values[fbase + i](f);
        }
      fbase += istep;
      ++f;
    }
  Assert(fbase == this->n_dofs_per_cell(), ExcInternalError());
}



// explicit instantiations
#include "fe/fe_rt_bubbles.inst"


DEAL_II_NAMESPACE_CLOSE
