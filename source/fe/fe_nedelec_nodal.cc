// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2005 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


#include <deal.II/base/polynomial.h>
#include <deal.II/base/qprojector.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_accessor.h>

#include <deal.II/fe/fe.h>
#include <deal.II/fe/fe_nedelec.h>
#include <deal.II/fe/fe_nedelec.templates.h>
#include <deal.II/fe/fe_nothing.h>
#include <deal.II/fe/fe_tools.h>
#include <deal.II/fe/mapping.h>

#include <deal.II/matrix_free/shape_info.templates.h>

#include <memory>
#include <sstream>


DEAL_II_NAMESPACE_OPEN


// ---------------- polynomial class for FE_NedelecNodal ---------------

namespace
{
  // Return a vector of "dofs per object" where the components of the returned
  // vector refer to:
  // 0 = vertex
  // 1 = edge
  // 2 = face (which is a cell in 2d)
  // 3 = cell.
  std::vector<unsigned int>
  get_nedelec_dpo_vector(const unsigned int dim, const unsigned int degree)
  {
    std::vector<unsigned int> dpo(dim + 1);
    dpo[0] = 0;
    dpo[1] = degree + 1;
    dpo[2] = 2 * degree * (degree + 1);

    if (dim == 3)
      dpo[3] = 3 * degree * degree * (degree + 1);

    return dpo;
  }

} // namespace



// --------------------- actual implementation of element --------------------

template <int dim>
FE_NedelecNodal<dim>::FE_NedelecNodal(const unsigned int degree)
  : FE_PolyTensor<dim>(
      PolynomialsVectorAnisotropic<dim>(degree,
                                        degree + 1,
                                        get_lexicographic_numbering(degree)),
      FiniteElementData<dim>(get_nedelec_dpo_vector(dim, degree),
                             dim,
                             degree + 1,
                             FiniteElementData<dim>::Hcurl),
      std::vector<bool>(
        PolynomialsVectorAnisotropic<dim>::n_polynomials(degree, degree + 1),
        true),
      std::vector<ComponentMask>(
        PolynomialsVectorAnisotropic<dim>::n_polynomials(degree, degree + 1),
        std::vector<bool>(dim, true)))
{
  Assert(dim >= 2, ExcImpossibleInDim(dim));

  const std::vector<unsigned int> &renumber =
    get_lexicographic_numbering(degree);
  this->mapping_kind = {mapping_nedelec};
  this->generalized_support_points =
    PolynomialsVectorAnisotropic<dim>(degree, degree + 1, renumber)
      .get_polynomial_support_points();

  AssertDimension(this->generalized_support_points.size(),
                  this->n_dofs_per_cell());

  PolynomialsVectorAnisotropic<dim> polynomials(degree, degree + 1, renumber);
  QMidpoint<1>                      gauss;
  QGaussLobatto<1>                  gl(degree + 2);
  this->unit_support_points.resize(this->dofs_per_cell);

  const unsigned int n_pols = polynomials.n() / dim;
  for (unsigned int k = 0; k < (dim > 2 ? degree + 2 : 1); ++k)
    for (unsigned int j = 0; j < (dim > 1 ? degree + 2 : 1); ++j)
      for (unsigned int i = 0; i < degree + 1; ++i)
        {
          this->unit_support_points
            [renumber[(k * (degree + 2) + j) * (degree + 1) + i]][0] =
            (degree == 0 ? gauss.point(i)[0] : gl.point(i)[0]);
          if (dim > 1)
            this->unit_support_points
              [renumber[(k * (degree + 2) + j) * (degree + 1) + i]][1] =
              gl.point(j)[0];
          if (dim > 2)
            this->unit_support_points
              [renumber[(k * (degree + 2) + j) * (degree + 1) + i]][2] =
              gl.point(k)[0];
        }
  if (dim > 1)
    for (unsigned int k = 0; k < (dim > 2 ? degree + 2 : 1); ++k)
      for (unsigned int j = 0; j < degree + 1; ++j)
        for (unsigned int i = 0; i < degree + 2; ++i)
          {
            this->unit_support_points
              [renumber[n_pols + (k * (degree + 1) + j) * (degree + 2) + i]]
              [0] = gl.point(i)[0];
            this->unit_support_points
              [renumber[n_pols + (k * (degree + 1) + j) * (degree + 2) + i]]
              [1] = (degree == 0 ? gauss.point(j)[0] : gl.point(j)[0]);
            if (dim > 2)
              this->unit_support_points
                [renumber[n_pols + (k * (degree + 1) + j) * (degree + 2) + i]]
                [2] = gl.point(k)[0];
          }
  if (dim > 2)
    for (unsigned int k = 0; k < degree + 1; ++k)
      for (unsigned int j = 0; j < degree + 2; ++j)
        for (unsigned int i = 0; i < degree + 2; ++i)
          {
            this->unit_support_points
              [renumber[2 * n_pols + (k * (degree + 2) + j) * (degree + 2) + i]]
              [0] = gl.point(i)[0];
            this->unit_support_points
              [renumber[2 * n_pols + (k * (degree + 2) + j) * (degree + 2) + i]]
              [1] = gl.point(j)[0];
            this->unit_support_points
              [renumber[2 * n_pols + (k * (degree + 2) + j) * (degree + 2) + i]]
              [2] = (degree == 0 ? gauss.point(k)[0] : gl.point(k)[0]);
          }
}



template <int dim>
std::string
FE_NedelecNodal<dim>::get_name() const
{
  /*Note, that the FETools::get_fe_by_name function depends on the particular
   format of the string this function returns. So, they have to be kept in
   sync.*/

  // Note, that this->degree is the maximal polynomial degree and is thus one
  // higher than the argument given to the constructor.
  return "FE_NedelecNodal<" + std::to_string(dim) + ">(" +
         std::to_string(this->degree - 1) + ")";
}


template <int dim>
std::unique_ptr<FiniteElement<dim, dim>>
FE_NedelecNodal<dim>::clone() const
{
  return std::make_unique<FE_NedelecNodal<dim>>(*this);
}

template <int dim>
void
FE_NedelecNodal<dim>::convert_generalized_support_point_values_to_dof_values(
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

  // First do interpolation on lines. The component evaluated there depends
  // on the face direction and orientation.
  unsigned int base = 0;
  for (unsigned int l = 0; l < GeometryInfo<dim>::lines_per_cell;
       ++l, base += this->n_dofs_per_line())
    {
      for (unsigned int i = 0; i < this->n_dofs_per_line(); ++i)
        {
          nodal_values[base + i] =
            support_point_values[base + i](l < 8 ? ((l % 4) / 2) ^ 1 : 2);
        }
    }

  if (dim == 3)
    {
      std::vector<unsigned int> componentFreq(dim,
                                              this->n_dofs_per_quad(0) / 2);
      unsigned int              index = -1;
      for (unsigned int f = 0; f < GeometryInfo<dim>::faces_per_cell; ++f)
        {
          index                = GeometryInfo<dim>::unit_normal_direction[f];
          componentFreq[index] = 0;
          for (unsigned int compInd = 0; compInd < dim; ++compInd)
            {
              for (unsigned int compQuant = 0;
                   compQuant < componentFreq[compInd];
                   ++compQuant)
                {
                  nodal_values[base + compQuant] =
                    support_point_values[base + compQuant](compInd);
                }
              base += componentFreq[compInd];
            }
          componentFreq[index] = this->n_dofs_per_quad(f) / 2;
        }
    }
  // The remaining points form dim chunks, one for each component.
  const unsigned int istep = (this->n_dofs_per_cell() - base) / dim;
  Assert((this->n_dofs_per_cell() - base) % dim == 0, ExcInternalError());

  int f = 0;
  while (base < this->n_dofs_per_cell())
    {
      for (unsigned int i = 0; i < istep; ++i)
        {
          nodal_values[base + i] = support_point_values[base + i](f);
        }
      base += istep;
      ++f;
    }
  Assert(base == this->n_dofs_per_cell(), ExcInternalError());
}



template <int dim>
std::vector<unsigned int>
FE_NedelecNodal<dim>::get_lexicographic_numbering(
  const unsigned int degree) const
{
  const unsigned int n_pols = (dim == 2) ?
                                (degree + 1) * (degree + 2) :
                                (degree + 1) * (degree + 2) * (degree + 2);

  std::vector<unsigned int> renumber_hierarchic_to_lexicographic;

  // Edges 0-3
  for (unsigned int i = 0; i < degree + 1; ++i)
    renumber_hierarchic_to_lexicographic.push_back(n_pols + i * (degree + 2));
  for (unsigned int i = 0; i < degree + 1; ++i)
    renumber_hierarchic_to_lexicographic.push_back(n_pols + degree + 1 +
                                                   i * (degree + 2));
  for (unsigned int i = 0; i < degree + 1; ++i)
    renumber_hierarchic_to_lexicographic.push_back(i);
  for (unsigned int i = 0; i < degree + 1; ++i)
    renumber_hierarchic_to_lexicographic.push_back((degree + 1) * (degree + 1) +
                                                   i);
  if (dim == 2)
    {
      // Quads for 2D
      for (unsigned int j = 1; j < degree + 1; ++j)
        for (unsigned int i = 0; i < degree + 1; ++i)
          renumber_hierarchic_to_lexicographic.push_back(j * (degree + 1) + i);
      for (unsigned int j = 0; j < degree + 1; ++j)
        for (unsigned int i = 1; i < degree + 1; ++i)
          renumber_hierarchic_to_lexicographic.push_back(n_pols +
                                                         j * (degree + 2) + i);
    }
  else if (dim == 3)
    {
      // Edges 4-7
      const unsigned int offset_z = (degree + 1) * (degree + 1) * (degree + 2);
      for (unsigned int i = 0; i < degree + 1; ++i)
        renumber_hierarchic_to_lexicographic.push_back(
          n_pols + i * (degree + 2) + offset_z);
      for (unsigned int i = 0; i < degree + 1; ++i)
        renumber_hierarchic_to_lexicographic.push_back(
          n_pols + degree + 1 + i * (degree + 2) + offset_z);
      for (unsigned int i = 0; i < degree + 1; ++i)
        renumber_hierarchic_to_lexicographic.push_back(i + offset_z);
      for (unsigned int i = 0; i < degree + 1; ++i)
        renumber_hierarchic_to_lexicographic.push_back(
          (degree + 1) * (degree + 1) + i + offset_z);

      // Edges 8-11
      for (unsigned int i = 0; i < degree + 1; ++i)
        renumber_hierarchic_to_lexicographic.push_back(
          2 * n_pols + i * (degree + 2) * (degree + 2));
      for (unsigned int i = 0; i < degree + 1; ++i)
        renumber_hierarchic_to_lexicographic.push_back(
          2 * n_pols + degree + 1 + i * (degree + 2) * (degree + 2));
      for (unsigned int i = 0; i < degree + 1; ++i)
        renumber_hierarchic_to_lexicographic.push_back(
          2 * n_pols + (degree + 1) * (degree + 2) +
          i * (degree + 2) * (degree + 2));
      for (unsigned int i = 0; i < degree + 1; ++i)
        renumber_hierarchic_to_lexicographic.push_back(
          2 * n_pols + (degree + 2) * (degree + 2) - 1 +
          i * (degree + 2) * (degree + 2));

      // Quad 0
      for (unsigned int j = 1; j < degree + 1; ++j)
        for (unsigned int i = 0; i < degree + 1; ++i)
          renumber_hierarchic_to_lexicographic.push_back(
            n_pols + i * (degree + 2) + j * (degree + 2) * (degree + 1));
      for (unsigned int j = 0; j < degree + 1; ++j)
        for (unsigned int i = 1; i < degree + 1; ++i)
          renumber_hierarchic_to_lexicographic.push_back(
            2 * n_pols + i * (degree + 2) + j * (degree + 2) * (degree + 2));
      // Quad 1
      for (unsigned int j = 1; j < degree + 1; ++j)
        for (unsigned int i = 0; i < degree + 1; ++i)
          renumber_hierarchic_to_lexicographic.push_back(
            n_pols + degree + 1 + i * (degree + 2) +
            j * (degree + 2) * (degree + 1));
      for (unsigned int j = 0; j < degree + 1; ++j)
        for (unsigned int i = 1; i < degree + 1; ++i)
          renumber_hierarchic_to_lexicographic.push_back(
            2 * n_pols + degree + 1 + i * (degree + 2) +
            j * (degree + 2) * (degree + 2));
      // Quad 2
      for (unsigned int j = 1; j < degree + 1; ++j)
        for (unsigned int i = 0; i < degree + 1; ++i)
          renumber_hierarchic_to_lexicographic.push_back(i + j * (degree + 2) *
                                                               (degree + 1));
      for (unsigned int j = 0; j < degree + 1; ++j)
        for (unsigned int i = 1; i < degree + 1; ++i)
          renumber_hierarchic_to_lexicographic.push_back(
            2 * n_pols + i + j * (degree + 2) * (degree + 2));
      // Quad 3
      for (unsigned int j = 1; j < degree + 1; ++j)
        for (unsigned int i = 0; i < degree + 1; ++i)
          renumber_hierarchic_to_lexicographic.push_back(
            (degree + 1) * (degree + 1) + i + j * (degree + 2) * (degree + 1));
      for (unsigned int j = 0; j < degree + 1; ++j)
        for (unsigned int i = 1; i < degree + 1; ++i)
          renumber_hierarchic_to_lexicographic.push_back(
            2 * n_pols + (degree + 1) * (degree + 2) + i +
            j * (degree + 2) * (degree + 2));
      // Quad 4
      for (unsigned int j = 1; j < degree + 1; ++j)
        for (unsigned int i = 0; i < degree + 1; ++i)
          renumber_hierarchic_to_lexicographic.push_back(i + j * (degree + 1));
      for (unsigned int j = 0; j < degree + 1; ++j)
        for (unsigned int i = 1; i < degree + 1; ++i)
          renumber_hierarchic_to_lexicographic.push_back(n_pols + i +
                                                         j * (degree + 2));
      // Quad 5
      for (unsigned int j = 1; j < degree + 1; ++j)
        for (unsigned int i = 0; i < degree + 1; ++i)
          renumber_hierarchic_to_lexicographic.push_back(
            (degree + 1) * (degree + 1) * (degree + 2) + i + j * (degree + 1));
      for (unsigned int j = 0; j < degree + 1; ++j)
        for (unsigned int i = 1; i < degree + 1; ++i)
          renumber_hierarchic_to_lexicographic.push_back(
            n_pols + (degree + 1) * (degree + 1) * (degree + 2) + i +
            j * (degree + 2));

      // Hexes
      for (unsigned int k = 1; k < degree + 1; ++k)
        for (unsigned int j = 1; j < degree + 1; ++j)
          for (unsigned int i = 0; i < degree + 1; ++i)
            renumber_hierarchic_to_lexicographic.push_back(
              k * (degree + 1) * (degree + 2) + j * (degree + 1) + i);
      for (unsigned int k = 1; k < degree + 1; ++k)
        for (unsigned int j = 0; j < degree + 1; ++j)
          for (unsigned int i = 1; i < degree + 1; ++i)
            renumber_hierarchic_to_lexicographic.push_back(
              n_pols + k * (degree + 2) * (degree + 1) + j * (degree + 2) + i);
      for (unsigned int k = 0; k < degree + 1; ++k)
        for (unsigned int j = 1; j < degree + 1; ++j)
          for (unsigned int i = 1; i < degree + 1; ++i)
            renumber_hierarchic_to_lexicographic.push_back(
              2 * n_pols + k * (degree + 2) * (degree + 2) + j * (degree + 2) +
              i);
    }

  return Utilities::invert_permutation(renumber_hierarchic_to_lexicographic);
}



template <int dim>
FiniteElementDomination::Domination
FE_NedelecNodal<dim>::compare_for_domination(const FiniteElement<dim> &fe_other,
                                             const unsigned int codim) const
{
  Assert(codim <= dim, ExcImpossibleInDim(dim));
  (void)codim;

  // vertex/line/face/cell domination
  // --------------------------------
  if (const FE_NedelecNodal<dim> *fe_nedelec_nodal_other =
        dynamic_cast<const FE_NedelecNodal<dim> *>(&fe_other))
    {
      if (this->degree < fe_nedelec_nodal_other->degree)
        return FiniteElementDomination::this_element_dominates;
      else if (this->degree == fe_nedelec_nodal_other->degree)
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


// Required for UnivariateShapeData from shape_info
#include "matrix_free/shape_info.inst"
// explicit instantiations
#include "fe/fe_nedelec_nodal.inst"



DEAL_II_NAMESPACE_CLOSE
