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
#include <lac/petsc_vector.h>
#include <lac/petsc_block_vector.h>
#include <lac/trilinos_vector.h>
#include <lac/trilinos_block_vector.h>

#include <numerics/mesh_worker_info.templates.h>

DEAL_II_NAMESPACE_OPEN

#if deal_II_dimension > 1

namespace MeshWorker
{
  template class LocalResults<double>;
  template class DoFInfo<deal_II_dimension,deal_II_dimension>;
  
#include "mesh_worker_info.inst"
  
  template void IntegrationInfo<deal_II_dimension, FEValuesBase<deal_II_dimension> >
  ::initialize<FEValues<deal_II_dimension> >(
    const FiniteElement<deal_II_dimension>&, const Mapping<deal_II_dimension>&,
    const Quadrature<FEValues<deal_II_dimension>::integral_dimension>&, const UpdateFlags);
  template void IntegrationInfo<deal_II_dimension, FEFaceValuesBase<deal_II_dimension> >
  ::initialize<FEFaceValues<deal_II_dimension> >(
    const FiniteElement<deal_II_dimension>&, const Mapping<deal_II_dimension>&,
    const Quadrature<FEFaceValues<deal_II_dimension>::integral_dimension>&, const UpdateFlags);
  template void IntegrationInfo<deal_II_dimension, FEFaceValuesBase<deal_II_dimension> >
  ::initialize<FESubfaceValues<deal_II_dimension> >(
    const FiniteElement<deal_II_dimension>&, const Mapping<deal_II_dimension>&,
    const Quadrature<FESubfaceValues<deal_II_dimension>::integral_dimension>&, const UpdateFlags);
}

#endif

DEAL_II_NAMESPACE_CLOSE

