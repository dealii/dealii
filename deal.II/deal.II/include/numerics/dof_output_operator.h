//---------------------------------------------------------------------------
//    $Id$
//
//    Copyright (C) 2010 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------

#ifndef __deal2__dof_output_operator_h
#define __deal2__dof_output_operator_h

#include <base/config.h>
#include <base/named_data.h>
#include <base/event.h>
#include <algorithms/operator.h>
#include <dofs/dof_handler.h>

#include <fstream>

DEAL_II_NAMESPACE_OPEN

namespace Algorithms
{
  template <class VECTOR, int dim, int spacedim=dim>
  class DoFOutputOperator : public OutputOperator<VECTOR>
  {
    public:
      void initialize (DoFHandler<dim, spacedim>& dof_handler);

      virtual OutputOperator<VECTOR>& 
        operator<<(const NamedData<VECTOR*>& vectors);

    private:
      SmartPointer<DoFHandler<dim, spacedim>,
        DoFOutputOperator<VECTOR, dim, spacedim> > dof;
  };

  template <class VECTOR, int dim, int spacedim>
  void DoFOutputOperator<VECTOR, dim, spacedim>::initialize(
      DoFHandler<dim, spacedim>& dof_handler)
  {
    dof = &dof_handler;  
  }
}


DEAL_II_NAMESPACE_CLOSE

#endif
