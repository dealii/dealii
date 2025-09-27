// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 1999 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


#ifndef dealii_mg_base_templates_h
#define dealii_mg_base_templates_h

#include <deal.II/base/config.h>

#include <deal.II/multigrid/mg_base.h>

DEAL_II_NAMESPACE_OPEN

template <typename VectorType>
void
MGSmootherBase<VectorType>::apply(const unsigned int level,
                                  VectorType        &u,
                                  const VectorType  &rhs) const
{
  u = typename VectorType::value_type(0.);
  smooth(level, u, rhs);
}



template <typename VectorType>
void
MGTransferBase<VectorType>::prolongate_and_add(const unsigned int to_level,
                                               VectorType        &dst,
                                               const VectorType  &src) const
{
  VectorType temp;
  temp.reinit(dst, true);

  this->prolongate(to_level, temp, src);

  dst += temp;
}

DEAL_II_NAMESPACE_CLOSE

#endif
