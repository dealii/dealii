//----------------------------  vector_memory.h  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001, 2002 by the deal.II authors
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


/**
 * Memory management for vectors. This class is used by all
 * iterative methods to allocate space for auxilliary
 * vectors. This class is used to avoid interaction with the
 * operating system whenever a vector is needed. Especially, when
 * an iterative method is invoked as part of an outer iteration,
 * this would lead to runtime overhead and memory fragmentation.
 *
 * Classes derived from this class implement a more or less
 * sophisticated management of vectors. One of these has to be
 * applied by the user according to his needs.
 */
template<class Vector = ::Vector<double> >
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
    virtual Vector* alloc() = 0;
    
				     /**
				      * Return a vector into the pool
				      * for later use.
				      */
    virtual void free (const Vector * const) = 0;
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


/**
 * Simple memory management.  This memory class is just made for
 * tests. It just allocates and deletes
 * vectors as needed from the global heap, i.e. performs no
 * specially adapted actions to the purpose of this class.
 */
template<class Vector = ::Vector<double> >
class PrimitiveVectorMemory : public VectorMemory<Vector>
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
    virtual Vector* alloc()
      {
	return new Vector();
      }
    
				     /**
				      * Return a vector to the global heap.
				      */
    virtual void free (const Vector * const v)
      {
	delete v;
      }
};


/**
 * Memory keeping allocated vectors.  This @p{VectorMemory} is able to
 * grow infinitely (according to computer memory).  A vector once
 * allocated will stay in the memory pool, though, and it will be
 * reused in later allocation.
 * 
 * @author Guido Kanschat, 1999
 */
template<class Vector = ::Vector<double> >
class GrowingVectorMemory : public VectorMemory<Vector>
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
    virtual Vector* alloc();
    
				     /**
				      * Return a vector into the pool
				      * for later use.
				      */
    virtual void free (const Vector * const);

  private:
				     /**
				      * Type to enter into the
				      * array. First component will be
				      * a flag telling whether the
				      * vector is used, second the
				      * vector itself.
				      */
    typedef typename std::pair<bool, Vector*> entry_type;

				     /**
				      * Array of allocated vectors.
				      */
    typename std::vector<entry_type> pool;
    
				     /**
				      * Overall number of allocations.
				      */
    unsigned int n_alloc;

				     /**
				      * Mutex to synchronise access to
				      * internal data of this object
				      * from multiple threads.
				      */
    Threads::ThreadMutex mutex;
};



/* --------------------- inline functions ---------------------- */


template <typename Vector>
GrowingVectorMemory<Vector>::GrowingVectorMemory(const unsigned int initial_size)
		: pool(initial_size)
{
  mutex.acquire ();
  for (typename std::vector<entry_type>::iterator i=pool.begin();
       i != pool.end();
       ++i)
    {
      i->first = false;
      i->second = new Vector;
    };

				   // no vectors yet claimed
  n_alloc = 0;
  mutex.release ();
}



template<typename Vector>
GrowingVectorMemory<Vector>::~GrowingVectorMemory()
{
  mutex.acquire ();

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
  mutex.release ();

				   // write out warning if memory leak
  if (n!=0)
    deallog << "GrowingVectorMemory:Vectors not deallocated: "
	    << n << ". Memory leak!" << std::endl;
}



template<typename Vector>
Vector*
GrowingVectorMemory<Vector>::alloc()
{
  mutex.acquire ();
  ++n_alloc;
  for (typename std::vector<entry_type>::iterator i=pool.begin();
       i != pool.end();
       ++i)
    {
      if (i->first == false)
	{
	  i->first = true;
	  mutex.release ();
	  return (i->second);
	}
    }
  entry_type t;
  t.first = true;
  t.second = new Vector;
  pool.push_back(t);
  mutex.release ();
  
  return t.second;
}



template<typename Vector>
void
GrowingVectorMemory<Vector>::free(const Vector* const v)
{
  mutex.acquire ();
  for (typename std::vector<entry_type>::iterator i=pool.begin();i != pool.end() ;++i)
    {
      if (v == (i->second))
	{
	  i->first = false;
	  mutex.release ();
	  return;
	}
    }
  mutex.release ();
  
  Assert(false, typename VectorMemory<Vector>::ExcNotAllocatedHere());
}


#endif
