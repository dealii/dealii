// ---------------------------------------------------------------------
// $Id$
//
// Copyright (C) 2010 - 2013 by the deal.II authors
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


#ifndef __deal2__dof_output_operator_h
#define __deal2__dof_output_operator_h

#include <deal.II/base/config.h>
#include <deal.II/base/named_data.h>
#include <deal.II/base/event.h>
#include <deal.II/algorithms/operator.h>
#include <deal.II/dofs/dof_handler.h>

#include <fstream>

DEAL_II_NAMESPACE_OPEN

namespace Algorithms
{
  template <class VECTOR, int dim, int spacedim=dim>
  class DoFOutputOperator : public OutputOperator<VECTOR>
  {
  public:
    void initialize (DoFHandler<dim, spacedim> &dof_handler);

    virtual OutputOperator<VECTOR> &
    operator<<(const NamedData<VECTOR *> &vectors);

  private:
    SmartPointer<DoFHandler<dim, spacedim>,
                 DoFOutputOperator<VECTOR, dim, spacedim> > dof;
  };

  template <class VECTOR, int dim, int spacedim>
  void DoFOutputOperator<VECTOR, dim, spacedim>::initialize(
    DoFHandler<dim, spacedim> &dof_handler)
  {
    dof = &dof_handler;
  }
}


DEAL_II_NAMESPACE_CLOSE

#endif
