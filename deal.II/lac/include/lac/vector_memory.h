//----------------------------  vector_memory.h  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001, 2002, 2003 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  vector_memory.h  ---------------------------
#ifndef __deal2__vector_memory_h
#define __deal2__vector_memory_h


#include <base/config.h>
#include <base/subscriptor.h>
#include <base/logstream.h>
#include <base/thread_management.h>
#include <lac/vector.h>

#include <vector>

/*!@addtogroup VMemory */
/*@{*/

//! Base class for vector memory management
/**
 * Memory management base class for vectors. This is an abstract base
 * class used by all iterative methods to allocate space for
 * auxilliary vectors. An object of a derived class with an actual
 * implementation of the memory management must be created by the
 * application program to be used in iterative methods.
 *
 * The purpose of this class is avoiding interaction with the
 * operating system whenever a vector is needed. Especially, when an
 * iterative method is invoked as part of an outer iteration, this
 * would lead to runtime overhead and memory fragmentation.
 *
 * @author Guido Kanschat, 1998-2003
*/
template<class VECTOR = ::Vector<double> >
class VectorMemory : public Subscriptor
{
  public:

				     /**
				      * Virtual destructor is needed
				      * as there are virtual functions
				      * in this class.
				      */
    virtual ~VectorMemory() {};

				     /**
				      * Return new vector from the pool.
				      */
    virtual VECTOR* alloc() = 0;
    
				     /**
				      * Return a vector into the pool
				      * for later use.
				      */
    virtual void free (const VECTOR * const) = 0;
				     /**
				      * No more available vectors.
				      */
    DeclException0(ExcNoMoreVectors);
				     /**
				      * Vector was not allocated from
				      * this memory pool.
				      */
    DeclException0(ExcNotAllocatedHere);
};

//! Sample implementation using system memory management.
/**
 * Simple memory management.  This memory class is just made for
 * tests. It just allocates and deletes
 * vectors as needed from the global heap, i.e. performs no
 * specially adapted actions to the purpose of this class.
 */
template<class VECTOR = ::Vector<double> >
class PrimitiveVectorMemory : public VectorMemory<VECTOR>
{
  public:
				     /**
				      * Constructor.
				      */
    PrimitiveVectorMemory () {}

				     /**
				      * Allocate a vector
				      * from the global heap.
				      */
    virtual VECTOR* alloc()
      {
	return new VECTOR();
      }
    
				     /**
				      * Return a vector to the global heap.
				      */
    virtual void free (const VECTOR * const v)
      {
	delete v;
      }
};

//! Keeps all vectors and avoids reallocation.
/**
 * Memory keeping allocated vectors.  This @p{VectorMemory} is able to
 * grow infinitely (according to computer memory).  A vector once
 * allocated will stay in the memory pool, though, and it will be
 * reused in later allocation.
 * 
 * @author Guido Kanschat, 1999
 */
template<class VECTOR = ::Vector<double> >
class GrowingVectorMemory : public VectorMemory<VECTOR>
{
  public:
				     /**
				      * Constructor.
				      * The argument allows to
				      * preallocate a certain number
				      * of vectors.
				      */
    GrowingVectorMemory (const unsigned int initial_size = 0);

				     /**
				      * Destructor.
				      * Release all vectors.
				      * This destructor also offers
				      * some statistic on the number
				      * of allocated vectors.
				      *
				      * The log file will also contain
				      * a warning message, if there
				      * are allocated vectors left.
				      */
    ~GrowingVectorMemory();
    
				     /**
				      * Return new vector from the pool.
				      */
    virtual VECTOR* alloc();
    
				     /**
				      * Return a vector into the pool
				      * for later use.
				      */
    virtual void free (const VECTOR * const);

  private:
				     /**
				      * Type to enter into the
				      * array. First component will be
				      * a flag telling whether the
				      * vector is used, second the
				      * vector itself.
				      */
    typedef std::pair<bool, VECTOR*> entry_type;

				     /**
				      * Array of allocated vectors.
				      */
    std::vector<entry_type> pool;
    
				     /**
				      * Overall number of allocations.
				      */
    unsigned int n_alloc;

				     /**
				      * Memory consumed by this class
				      * and all currently allocated
				      * vectors.
				      */
    unsigned int memory_consumption() const;
    
				     /**
				      * Mutex to synchronise access to
				      * internal data of this object
				      * from multiple threads.
				      */
    Threads::ThreadMutex mutex;
};

/*@}*/

/* --------------------- inline functions ---------------------- */


template <typename VECTOR>
GrowingVectorMemory<VECTOR>::GrowingVectorMemory(const unsigned int initial_size)
		: pool(initial_size)
{
  Threads::ThreadMutex::ScopedLock lock(mutex);
  for (typename std::vector<entry_type>::iterator i=pool.begin();
       i != pool.end();
       ++i)
    {
      i->first = false;
      i->second = new VECTOR;
    };

				   // no vectors yet claimed
  n_alloc = 0;
}



template<typename VECTOR>
GrowingVectorMemory<VECTOR>::~GrowingVectorMemory()
{
  Threads::ThreadMutex::ScopedLock lock(mutex);

				   // deallocate all vectors and count
				   // number of vectors that is still
				   // claimed by other users
  unsigned int n = 0;
  for (typename std::vector<entry_type>::iterator i=pool.begin();
       i != pool.end();
       ++i)
    {
      if (i->first == true)
	++n;
      delete i->second;
    }
  deallog << "GrowingVectorMemory:Overall allocated vectors: "
	  << n_alloc << std::endl;
  deallog << "GrowingVectorMemory:Maximum allocated vectors: "
	  << pool.size() << std::endl;
  pool.clear ();

				   // write out warning if memory leak
  if (n!=0)
    deallog << "GrowingVectorMemory:Vectors not deallocated: "
	    << n << ". Memory leak!" << std::endl;
}



template<typename VECTOR>
VECTOR *
GrowingVectorMemory<VECTOR>::alloc()
{
  Threads::ThreadMutex::ScopedLock lock(mutex);
  ++n_alloc;
  for (typename std::vector<entry_type>::iterator i=pool.begin();
       i != pool.end();
       ++i)
    {
      if (i->first == false)
	{
	  i->first = true;
	  return (i->second);
	}
    }
  
  const entry_type t (true, new VECTOR);
  pool.push_back(t);
  
  return t.second;
}



template<typename VECTOR>
void
GrowingVectorMemory<VECTOR>::free(const VECTOR* const v)
{
  Threads::ThreadMutex::ScopedLock lock(mutex);
  for (typename std::vector<entry_type>::iterator i=pool.begin();i != pool.end() ;++i)
    {
      if (v == (i->second))
	{
	  i->first = false;
	  return;
	}
    }
  Assert(false, typename VectorMemory<VECTOR>::ExcNotAllocatedHere());
}


template<typename VECTOR>
unsigned int
GrowingVectorMemory<VECTOR>::memory_consumption () const
{
  unsigned int result = sizeof (*this);
  const typename std::vector<entry_type>::const_iterator
    end = pool.end();
  for (typename std::vector<entry_type>::const_iterator i = pool.begin()
						     ; i != end ; ++i)
    result += i->second->memory_consumption();
}

#endif
