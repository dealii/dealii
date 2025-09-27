// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2009 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


#include <deal.II/lac/block_indices.h>

#include <deal.II/meshworker/local_integrator.h>
#include <deal.II/meshworker/local_results.h>

DEAL_II_NAMESPACE_OPEN

namespace MeshWorker
{
  template <typename number>
  void
  LocalResults<number>::reinit(const BlockIndices &bi)
  {
    for (unsigned int i = 0; i < J.size(); ++i)
      J[i] = 0.;
    for (unsigned int i = 0; i < R.size(); ++i)
      R[i].reinit(bi);
    for (unsigned int i = 0; i < M1.size(); ++i)
      M1[i].matrix.reinit(bi.block_size(M1[i].row),
                          bi.block_size(M1[i].column));
    for (unsigned int i = 0; i < M2.size(); ++i)
      M2[i].matrix.reinit(bi.block_size(M2[i].row),
                          bi.block_size(M2[i].column));
    quadrature_data.reset_values();
  }


  template <typename number>
  std::size_t
  LocalResults<number>::memory_consumption() const
  {
    std::size_t mem = sizeof(*this) + MemoryConsumption::memory_consumption(J) +
                      MemoryConsumption::memory_consumption(R) +
                      MemoryConsumption::memory_consumption(M1) +
                      MemoryConsumption::memory_consumption(M2) +
                      MemoryConsumption::memory_consumption(quadrature_data);
    return mem;
  }


  template class LocalResults<float>;
  template class LocalResults<double>;

  template <int dim, int spacedim, typename number>
  LocalIntegrator<dim, spacedim, number>::LocalIntegrator()
    : use_cell(true)
    , use_boundary(true)
    , use_face(true)
  {}


  template <int dim, int spacedim, typename number>
  LocalIntegrator<dim, spacedim, number>::LocalIntegrator(bool c,
                                                          bool b,
                                                          bool f)
    : use_cell(c)
    , use_boundary(b)
    , use_face(f)
  {}



  template <int dim, int spacedim, typename number>
  void
  LocalIntegrator<dim, spacedim, number>::cell(
    DoFInfo<dim, spacedim, number> &,
    IntegrationInfo<dim, spacedim> &) const
  {
    Assert(false, ExcPureFunction());
  }


  template <int dim, int spacedim, typename number>
  void
  LocalIntegrator<dim, spacedim, number>::boundary(
    DoFInfo<dim, spacedim, number> &,
    IntegrationInfo<dim, spacedim> &) const
  {
    Assert(false, ExcPureFunction());
  }


  template <int dim, int spacedim, typename number>
  void
  LocalIntegrator<dim, spacedim, number>::face(
    DoFInfo<dim, spacedim, number> &,
    DoFInfo<dim, spacedim, number> &,
    IntegrationInfo<dim, spacedim> &,
    IntegrationInfo<dim, spacedim> &) const
  {
    Assert(false, ExcPureFunction());
  }


  template class LocalIntegrator<1, 1, float>;
  template class LocalIntegrator<1, 1, double>;
  template class LocalIntegrator<1, 2, float>;
  template class LocalIntegrator<1, 2, double>;
  template class LocalIntegrator<1, 3, float>;
  template class LocalIntegrator<1, 3, double>;
  template class LocalIntegrator<2, 2, float>;
  template class LocalIntegrator<2, 2, double>;
  template class LocalIntegrator<2, 3, float>;
  template class LocalIntegrator<2, 3, double>;
  template class LocalIntegrator<3, 3, float>;
  template class LocalIntegrator<3, 3, double>;
} // namespace MeshWorker


DEAL_II_NAMESPACE_CLOSE
