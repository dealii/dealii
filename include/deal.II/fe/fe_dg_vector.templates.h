// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2010 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_fe_dg_vector_templates_h
#define dealii_fe_dg_vector_templates_h


#include <deal.II/base/config.h>

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/fe/fe_dg_vector.h>
#include <deal.II/fe/fe_tools.h>

#include <memory>

DEAL_II_NAMESPACE_OPEN


// TODO:[GK] deg+1 is wrong here and should be fixed after FiniteElementData was
// cleaned up

template <typename PolynomialType, int dim, int spacedim>
FE_DGVector<PolynomialType, dim, spacedim>::FE_DGVector(const unsigned int deg,
                                                        MappingKind        map)
  : FE_PolyTensor<dim, spacedim>(
      PolynomialType(deg),
      FiniteElementData<dim>(get_dpo_vector(deg),
                             dim,
                             deg + 1,
                             FiniteElementData<dim>::L2),
      std::vector<bool>(PolynomialType::n_polynomials(deg), true),
      std::vector<ComponentMask>(PolynomialType::n_polynomials(deg),
                                 ComponentMask(dim, true)))
{
  this->mapping_kind                   = {map};
  const unsigned int polynomial_degree = this->tensor_degree();

  QGauss<dim> quadrature(polynomial_degree + 1);
  this->generalized_support_points = quadrature.get_points();

  this->reinit_restriction_and_prolongation_matrices(true, true);
  FETools::compute_projection_matrices(*this, this->restriction, true);
  FETools::compute_embedding_matrices(*this, this->prolongation, true);
}


template <typename PolynomialType, int dim, int spacedim>
std::unique_ptr<FiniteElement<dim, spacedim>>
FE_DGVector<PolynomialType, dim, spacedim>::clone() const
{
  return std::make_unique<FE_DGVector<PolynomialType, dim, spacedim>>(*this);
}


template <typename PolynomialType, int dim, int spacedim>
std::string
FE_DGVector<PolynomialType, dim, spacedim>::get_name() const
{
  std::ostringstream namebuf;
  namebuf << "FE_DGVector_" << this->poly_space->name() << '<' << dim << ">("
          << this->degree - 1 << ')';
  return namebuf.str();
}


template <typename PolynomialType, int dim, int spacedim>
std::vector<unsigned int>
FE_DGVector<PolynomialType, dim, spacedim>::get_dpo_vector(
  const unsigned int deg)
{
  std::vector<unsigned int> dpo(dim + 1);
  dpo[dim] = PolynomialType::n_polynomials(deg);

  return dpo;
}


template <typename PolynomialType, int dim, int spacedim>
bool
FE_DGVector<PolynomialType, dim, spacedim>::has_support_on_face(
  const unsigned int,
  const unsigned int) const
{
  return true;
}


template <typename PolynomialType, int dim, int spacedim>
std::size_t
FE_DGVector<PolynomialType, dim, spacedim>::memory_consumption() const
{
  DEAL_II_NOT_IMPLEMENTED();
  return 0;
}

DEAL_II_NAMESPACE_CLOSE

#endif
