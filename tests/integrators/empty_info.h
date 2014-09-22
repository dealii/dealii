// ---------------------------------------------------------------------
//
// Copyright (C) 2000 - 2013 by the deal.II authors
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

// Provide scratch data objects with no contents for the MeshWorker::loop()

#include <deal.II/meshworker/dof_info.h>

class EmptyInfo
{
public:
  template <class DOFINFO>
  void reinit(const DOFINFO &)
  {}
};


class EmptyInfoBox
{
public:
  typedef EmptyInfo CellInfo;
  template <int dim, class DOFINFO>
  void post_cell(const MeshWorker::DoFInfoBox<dim, DOFINFO> &)
  {}

  template <int dim, class DOFINFO>
  void post_faces(const MeshWorker::DoFInfoBox<dim, DOFINFO> &)
  {}

  EmptyInfo cell;
  EmptyInfo boundary;
  EmptyInfo face;
  EmptyInfo subface;
  EmptyInfo neighbor;
};


