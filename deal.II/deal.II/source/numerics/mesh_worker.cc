//---------------------------------------------------------------------------
//    $Id$
//
//    Copyright (C) 2006, 2007, 2008, 2009, 2010 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------

#include <numerics/mesh_worker.h>
#include <lac/block_indices.h>

DEAL_II_NAMESPACE_OPEN

using namespace MeshWorker;

template <typename number>
void
LocalResults<number>::reinit(const BlockIndices& bi)
{
  for (unsigned int i=0;i<J.size();++i)
      J[i] = 0.;
  for (unsigned int i=0;i<R.size();++i)
    R[i].reinit(bi);
  for (unsigned int i=0;i<M1.size();++i)
    M1[i].matrix.reinit(bi.block_size(M1[i].row),
			bi.block_size(M1[i].column));
  for (unsigned int i=0;i<M2.size();++i)
    M2[i].matrix.reinit(bi.block_size(M2[i].row),
			bi.block_size(M2[i].column));
  for (unsigned int k=0;k<quadrature_values.size();++k)
    for (unsigned int i=0;i<quadrature_values[k].size();++i)
      quadrature_values[k][i] = 0.;
}


template <typename number>
unsigned int
LocalResults<number>::memory_consumption () const
{
  unsigned int mem = sizeof(*this)
		     + MemoryConsumption::memory_consumption(J)
		     + MemoryConsumption::memory_consumption(R)
		     + MemoryConsumption::memory_consumption(M1)
		     + MemoryConsumption::memory_consumption(M2)
		     + MemoryConsumption::memory_consumption(quadrature_values);
  return mem;
}


template class LocalResults<float>;
template class LocalResults<double>;

DEAL_II_NAMESPACE_CLOSE
