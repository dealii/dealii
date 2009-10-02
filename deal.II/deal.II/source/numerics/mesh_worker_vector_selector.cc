//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2009 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------

#include <lac/vector.h>
#include <lac/block_vector.h>

#include <numerics/mesh_worker_vector_selector.templates.h>

DEAL_II_NAMESPACE_OPEN

namespace MeshWorker
{
  template class VectorDataBase<deal_II_dimension>;

#include "mesh_worker_vector_selector.inst"
}

DEAL_II_NAMESPACE_CLOSE
