//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2007 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------


#include <lac/vector_memory.h>
#include <lac/vector.h>
#include <lac/block_vector.h>

DEAL_II_NAMESPACE_OPEN

  
namespace
{
  GrowingVectorMemory<Vector<double> > default_pool_Vector_double;
  GrowingVectorMemory<Vector<float> > default_pool_Vector_float;
  GrowingVectorMemory<BlockVector<double> > default_pool_BlockVector_double;
  GrowingVectorMemory<BlockVector<float> > default_pool_BlockVector_float;
  
  template<class VECTOR>
  inline
  VectorMemory<VECTOR>*
  default_pool_select()
  {
    Assert(false, typename VectorMemory<VECTOR>::ExcNoDefaultMemoryForThisVectorClass());
    return 0;
  }


  template <>
  inline
  VectorMemory<Vector<double> >*
  default_pool_select<Vector<double> >()
  {
    return &default_pool_Vector_double;
  }
  
  
  template <>
  inline
  VectorMemory<Vector<float> >*
  default_pool_select<Vector<float> >()
  {
    return &default_pool_Vector_float;
  }
  
  
  
  template <>
  inline
  VectorMemory<BlockVector<double> >*
  default_pool_select<BlockVector<double> >()
  {
    return &default_pool_BlockVector_double;
  }
  
  
  template <>
  inline
  VectorMemory<BlockVector<float> >*
  default_pool_select<BlockVector<float> >()
  {
    return &default_pool_BlockVector_float;
  }
  
}


template <class VECTOR>
VectorMemory<VECTOR>&
VectorMemory<VECTOR>::default_pool()
{
  return *default_pool_select<VECTOR>();
}


template class VectorMemory<Vector<double> >;
template class VectorMemory<Vector<float> >;
template class VectorMemory<BlockVector<double> >;
template class VectorMemory<BlockVector<float> >;


DEAL_II_NAMESPACE_CLOSE
