// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2011 - 2021 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Provide scratch data objects with no contents for the MeshWorker::loop()

#include <deal.II/meshworker/dof_info.h>

class EmptyInfo
{
public:
  template <class DOFINFO>
  void
  reinit(const DOFINFO &)
  {}
};


class EmptyInfoBox
{
public:
  using CellInfo = EmptyInfo;
  template <int dim, class DOFINFO>
  void
  post_cell(const MeshWorker::DoFInfoBox<dim, DOFINFO> &)
  {}

  template <int dim, class DOFINFO>
  void
  post_faces(const MeshWorker::DoFInfoBox<dim, DOFINFO> &)
  {}

  EmptyInfo cell;
  EmptyInfo boundary;
  EmptyInfo face;
  EmptyInfo subface;
  EmptyInfo neighbor;
};
