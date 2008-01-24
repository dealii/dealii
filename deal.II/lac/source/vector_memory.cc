//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2007, 2008 by the deal.II authors
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


template <typename VECTOR>
GrowingVectorMemory<VECTOR>::Pool::Pool()
		:
		data(0)
{}



template <typename VECTOR>
GrowingVectorMemory<VECTOR>::Pool::~Pool()
{
				   // Nothing to do if memory was unused.
  if (data == 0) return;

				   // First, delete all remaining
				   // vectors. Actually, there should
				   // be none, if there is no memory
				   // leak
  unsigned int n=0;
  for (typename std::vector<entry_type>::iterator i=data->begin();
       i != data->end();
       ++i)
    {
      if (i->first == true)
	++n;
      delete i->second;
    }
  delete data;
}


template <typename VECTOR>
void
GrowingVectorMemory<VECTOR>::Pool::initialize(const unsigned int size)
{
  if (data == 0)
    {
      data = new std::vector<entry_type>(size);
      
      for (typename std::vector<entry_type>::iterator i= data->begin();
	   i != data->end();
	   ++i)
	{
	  i->first = false;
	  i->second = new VECTOR;
	}
    }
}


template <typename VECTOR>
typename GrowingVectorMemory<VECTOR>::Pool GrowingVectorMemory<VECTOR>::pool;

template <typename VECTOR>
Threads::ThreadMutex GrowingVectorMemory<VECTOR>::mutex;





#include "vector_memory.inst"

DEAL_II_NAMESPACE_CLOSE
