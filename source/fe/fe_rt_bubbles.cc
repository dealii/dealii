// ---------------------------------------------------------------------
//
// Copyright (C) 2018 - 2019 by the deal.II authors
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
                                 std::vector<bool>(dim, true)))
{
  Assert(dim >= 2, ExcImpossibleInDim(dim));
  Assert(
    deg >= 1,
    ExcMessage(
      "Lowest order RT_Bubbles element is degree 1, but you requested for degree 0"));
  const unsigned int n_dofs = this->dofs_per_cell;

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
  for (unsigned int ref_case = RefinementCase<dim>::cut_x;
       ref_case < RefinementCase<dim>::isotropic_refinement + 1;
       ++ref_case)
    {
      const unsigned int nc =
        GeometryInfo<dim>::n_children(RefinementCase<dim>(ref_case));

      for (unsigned int i = 0; i < nc; ++i)
        this->prolongation[ref_case - 1][i].reinit(n_dofs, n_dofs);
    }
  // Fill prolongation matrices with embedding operators
  // set tolerance to 1, as embedding error accumulate quickly
  FETools::compute_embedding_matrices(*this, this->prolongation, true, 1.0);
  FullMatrix<double> face_embeddings[GeometryInfo<dim>::max_children_per_face];
  for (unsigned int i = 0; i < GeometryInfo<dim>::max_children_per_face; ++i)
    face_embeddings[i].reinit(this->dofs_per_face, this->dofs_per_face);
  FETools::compute_face_embedding_matrices<dim, double>(*this,
                                                        face_embeddings,
                                                        0,
                                                        0);
  this->interface_constraints.reinit((1 << (dim - 1)) * this->dofs_per_face,
                                     this->dofs_per_face);
  unsigned int target_row = 0;
  for (unsigned int d = 0; d < GeometryInfo<dim>::max_children_per_face; ++d)
    for (unsigned int i = 0; i < face_embeddings[d].m(); ++i)
      {
        for (unsigned int j = 0; j < face_embeddings[d].n(); ++j)
          this->interface_constraints(target_row, j) = face_embeddings[d](i, j);
        ++target_row;
      }
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
  this->generalized_support_points.resize(this->dofs_per_cell);
  this->generalized_face_support_points.resize(this->dofs_per_face);

  // Index of the point being entered
  unsigned int current = 0;

  // On the faces, we choose as many Gauss-Lobatto points
  // as required to determine the normal component uniquely.
  // This is the deg of the RT_Bubble element plus one.
  if (dim > 1)
    {
      QGaussLobatto<dim - 1> face_points(deg + 1);
      Assert(face_points.size() == this->dofs_per_face, ExcInternalError());
      for (unsigned int k = 0; k < this->dofs_per_face; ++k)
        this->generalized_face_support_points[k] = face_points.point(k);
      Quadrature<dim> faces =
        QProjector<dim>::project_to_all_faces(face_points);
      for (unsigned int k = 0;
           k < this->dofs_per_face * GeometryInfo<dim>::faces_per_cell;
           ++k)
        this->generalized_support_points[k] =
          faces.point(k + QProjector<dim>::DataSetDescriptor::face(
                            0, true, false, false, this->dofs_per_face));

      current = this->dofs_per_face * GeometryInfo<dim>::faces_per_cell;
    }

  if (deg == 1)
    return;

  // In the interior, we need anisotropic Gauss-Lobatto quadratures,
  // one for each direction
  QGaussLobatto<1>      high(deg + 1);
  std::vector<Point<1>> pts = high.get_points();
  pts.erase(pts.begin());
  pts.erase(pts.end() - 1);

  std::vector<double> wts(pts.size(), 1);
  Quadrature<1>       low(pts, wts);

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
            Assert(false, ExcNotImplemented());
        }

      for (unsigned int k = 0; k < quadrature->size(); ++k)
        this->generalized_support_points[current++] = quadrature->point(k);
    }
  Assert(current == this->dofs_per_cell, ExcInternalError());
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
  std::vector<double> &              nodal_values) const
{
  Assert(support_point_values.size() == this->generalized_support_points.size(),
         ExcDimensionMismatch(support_point_values.size(),
                              this->generalized_support_points.size()));
  Assert(nodal_values.size() == this->dofs_per_cell,
         ExcDimensionMismatch(nodal_values.size(), this->dofs_per_cell));
  Assert(support_point_values[0].size() == this->n_components(),
         ExcDimensionMismatch(support_point_values[0].size(),
                              this->n_components()));

  // First do interpolation on faces. There, the component
  // evaluated depends on the face direction and orientation.
  unsigned int fbase = 0;
  unsigned int f     = 0;
  for (; f < GeometryInfo<dim>::faces_per_cell;
       ++f, fbase += this->dofs_per_face)
    {
      for (unsigned int i = 0; i < this->dofs_per_face; ++i)
        {
          nodal_values[fbase + i] = support_point_values[fbase + i](
            GeometryInfo<dim>::unit_normal_direction[f]);
        }
    }

  // The remaining points form dim chunks, one for each component.
  const unsigned int istep = (this->dofs_per_cell - fbase) / dim;
  Assert((this->dofs_per_cell - fbase) % dim == 0, ExcInternalError());

  f = 0;
  while (fbase < this->dofs_per_cell)
    {
      for (unsigned int i = 0; i < istep; ++i)
        {
          nodal_values[fbase + i] = support_point_values[fbase + i](f);
        }
      fbase += istep;
      ++f;
    }
  Assert(fbase == this->dofs_per_cell, ExcInternalError());
}



// explicit instantiations
#include "fe_rt_bubbles.inst"


DEAL_II_NAMESPACE_CLOSE
