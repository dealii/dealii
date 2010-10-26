//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2007, 2008, 2010 by the deal.II authors
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
#include <lac/petsc_vector.h>
#include <lac/petsc_block_vector.h>
#include <lac/petsc_parallel_vector.h>
#include <lac/petsc_parallel_block_vector.h>
#include <lac/trilinos_vector.h>
#include <lac/trilinos_block_vector.h>

DEAL_II_NAMESPACE_OPEN


template <typename VECTOR>
typename GrowingVectorMemory<VECTOR>::Pool GrowingVectorMemory<VECTOR>::pool;

template <typename VECTOR>
Threads::ThreadMutex GrowingVectorMemory<VECTOR>::mutex;

template <typename VECTOR>
inline
GrowingVectorMemory<VECTOR>::Pool::Pool()
		:
		data(0)
{}



template <typename VECTOR>
inline
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
inline
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
inline
GrowingVectorMemory<VECTOR>::GrowingVectorMemory (const unsigned int initial_size,
						  const bool log_statistics)

		:
		total_alloc(0),
		current_alloc(0),
		log_statistics(log_statistics)
{
  Threads::ThreadMutex::ScopedLock lock(mutex);
  pool.initialize(initial_size);
}


template<typename VECTOR>
inline
GrowingVectorMemory<VECTOR>::~GrowingVectorMemory()
{
  AssertThrow(current_alloc == 0,
	      StandardExceptions::ExcMemoryLeak(current_alloc));
  if (log_statistics)
    {
      deallog << "GrowingVectorMemory:Overall allocated vectors: "
	      << total_alloc << std::endl;
      deallog << "GrowingVectorMemory:Maximum allocated vectors: "
	      << pool.data->size() << std::endl;
    }
}



template<typename VECTOR>
inline
VECTOR *
GrowingVectorMemory<VECTOR>::alloc ()
{
  Threads::ThreadMutex::ScopedLock lock(mutex);
  ++total_alloc;
  ++current_alloc;
				   // see if there is a free vector
				   // available in our list
  for (typename std::vector<entry_type>::iterator i=pool.data->begin();
       i != pool.data->end(); ++i)
    {
      if (i->first == false)
	{
	  i->first = true;
	  return (i->second);
	}
    }

				   // no free vector found, so let's
				   // just allocate a new one
  const entry_type t (true, new VECTOR);
  pool.data->push_back(t);

  return t.second;
}



template<typename VECTOR>
inline
void
GrowingVectorMemory<VECTOR>::free(const VECTOR* const v)
{
  Threads::ThreadMutex::ScopedLock lock(mutex);
  for (typename std::vector<entry_type>::iterator i=pool.data->begin();
       i != pool.data->end(); ++i)
    {
      if (v == (i->second))
	{
	  i->first = false;
	  --current_alloc;
	  return;
	}
    }
  Assert(false, typename VectorMemory<VECTOR>::ExcNotAllocatedHere());
}



template<typename VECTOR>
inline
void
GrowingVectorMemory<VECTOR>::release_unused_memory ()
{
  Threads::ThreadMutex::ScopedLock lock(mutex);

  std::vector<entry_type> new_data;

  if (pool.data != 0)
    {
      const typename std::vector<entry_type>::const_iterator
	end = pool.data->end();
      for (typename std::vector<entry_type>::const_iterator
	     i = pool.data->begin(); i != end ; ++i)
	if (i->first == false)
	  delete i->second;
	else
	  new_data.push_back (*i);

      *pool.data = new_data;
    }
}



template<typename VECTOR>
inline
unsigned int
GrowingVectorMemory<VECTOR>::memory_consumption () const
{
  Threads::ThreadMutex::ScopedLock lock(mutex);

  unsigned int result = sizeof (*this);
  const typename std::vector<entry_type>::const_iterator
    end = pool.data->end();
  for (typename std::vector<entry_type>::const_iterator
         i = pool.data->begin(); i != end ; ++i)
    result += sizeof (*i) + i->second->memory_consumption();

  return result;
}


// -------------------------------------------------------------
// explicit instantiations


#include "vector_memory.inst"

//TODO: Fold this into the list of vectors to be instantiated
#ifdef DEAL_II_USE_PETSC
    template class VectorMemory<PETScWrappers::Vector>;
    template class GrowingVectorMemory<PETScWrappers::Vector>;

    template class VectorMemory<PETScWrappers::BlockVector>;
    template class GrowingVectorMemory<PETScWrappers::BlockVector>;

    template class VectorMemory<PETScWrappers::MPI::Vector>;
    template class GrowingVectorMemory<PETScWrappers::MPI::Vector>;

    template class VectorMemory<PETScWrappers::MPI::BlockVector>;
    template class GrowingVectorMemory<PETScWrappers::MPI::BlockVector>;
#endif

#ifdef DEAL_II_USE_TRILINOS
    template class VectorMemory<TrilinosWrappers::Vector>;
    template class GrowingVectorMemory<TrilinosWrappers::Vector>;

    template class VectorMemory<TrilinosWrappers::BlockVector>;
    template class GrowingVectorMemory<TrilinosWrappers::BlockVector>;

    template class VectorMemory<TrilinosWrappers::MPI::Vector>;
    template class GrowingVectorMemory<TrilinosWrappers::MPI::Vector>;

    template class VectorMemory<TrilinosWrappers::MPI::BlockVector>;
    template class GrowingVectorMemory<TrilinosWrappers::MPI::BlockVector>;
#endif

DEAL_II_NAMESPACE_CLOSE
