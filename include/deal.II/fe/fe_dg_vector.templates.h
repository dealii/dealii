// ---------------------------------------------------------------------
//
// Copyright (C) 2006 - 2014 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------


#include <deal.II/fe/fe_dg_vector.h>
#include <deal.II/fe/fe_tools.h>
#include <deal.II/base/quadrature_lib.h>

DEAL_II_NAMESPACE_OPEN


//TODO:[GK] deg+1 is wrong here and should bew fixed after FiniteElementData was cleaned up

template <class POLY, int dim, int spacedim>
FE_DGVector<POLY,dim,spacedim>::FE_DGVector (
  const unsigned int deg, MappingType map)
  :
  FE_PolyTensor<POLY, dim, spacedim>(
    deg,
    FiniteElementData<dim>(
      get_dpo_vector(deg), dim, deg+1, FiniteElementData<dim>::L2, 1),
    std::vector<bool>(POLY::compute_n_pols(deg), true),
    std::vector<ComponentMask>(POLY::compute_n_pols(deg),
                               ComponentMask(dim,true)))
{
  this->mapping_type = map;
  const unsigned int polynomial_degree = this->tensor_degree();

  QGauss<dim> quadrature(polynomial_degree+1);
  this->generalized_support_points = quadrature.get_points();

  this->reinit_restriction_and_prolongation_matrices(true, true);
  FETools::compute_projection_matrices (*this, this->restriction, true);
  FETools::compute_embedding_matrices (*this, this->prolongation, true);
}


template <class POLY, int dim, int spacedim>
FiniteElement<dim, spacedim> *
FE_DGVector<POLY,dim,spacedim>::clone() const
{
  return new FE_DGVector<POLY, dim, spacedim>(*this);
}


template <class POLY, int dim, int spacedim>
std::string
FE_DGVector<POLY,dim,spacedim>::get_name() const
{
  std::ostringstream namebuf;
  namebuf << "FE_DGVector_" << this->poly_space.name()
          << "<" << dim << ">(" << this->degree-1 << ")";
  return namebuf.str();
}


template <class POLY, int dim, int spacedim>
std::vector<unsigned int>
FE_DGVector<POLY,dim,spacedim>::get_dpo_vector (const unsigned int deg)
{
  std::vector<unsigned int> dpo(dim+1);
  dpo[dim] = POLY::compute_n_pols(deg);

  return dpo;
}


template <class POLY, int dim, int spacedim>
bool
FE_DGVector<POLY,dim,spacedim>::has_support_on_face (const unsigned int,
                                                     const unsigned int) const
{
  return true;
}


template <class POLY, int dim, int spacedim>
void
FE_DGVector<POLY,dim,spacedim>::interpolate(
  std::vector<double> &,
  const std::vector<double> &) const
{
  Assert(false, ExcNotImplemented());
}


template <class POLY, int dim, int spacedim>
void
FE_DGVector<POLY,dim,spacedim>::interpolate(
  std::vector<double> & /*local_dofs*/,
  const std::vector<Vector<double> > & /*values*/,
  unsigned int /*offset*/) const
{
  Assert(false, ExcNotImplemented());
}

template <class POLY, int dim, int spacedim>
void
FE_DGVector<POLY,dim,spacedim>::interpolate(
  std::vector<double> & /*local_dofs*/,
  const VectorSlice<const std::vector<std::vector<double> > > & /*values*/) const
{
  Assert(false, ExcNotImplemented());
}


template <class POLY, int dim, int spacedim>
std::size_t
FE_DGVector<POLY,dim,spacedim>::memory_consumption() const
{
  Assert(false, ExcNotImplemented());
  return 0;
}

DEAL_II_NAMESPACE_CLOSE
