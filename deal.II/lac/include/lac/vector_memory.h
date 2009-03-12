//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008, 2009 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------
#ifndef __deal2__vector_memory_h
#define __deal2__vector_memory_h


#include <base/config.h>
#include <base/smartpointer.h>
#include <base/logstream.h>
#include <base/thread_management.h>
#include <lac/vector.h>

#include <vector>
#include <iostream>

DEAL_II_NAMESPACE_OPEN


/*!@addtogroup VMemory */
/*@{*/

/**
 * Memory management base class for vectors. This is an abstract base
 * class used, among other places, by all iterative methods to
 * allocate space for auxiliary vectors.
 *
 * The purpose of this class is as follows: in iterative solvers and
 * other places, one needs to allocate temporary storage for vectors,
 * for example for auxiliary vectors. One could allocate and release
 * them anew every time, but this may be expensive in some situations
 * if it has to happen very frequently. A common case for this is when
 * an iterative method is used to invert a matrix in each iteration an
 * outer solver, such as when inverting a matrix block for a Schur
 * complement solver.
 *
 * In such situations, allocating and deallocating vectors anew in
 * each call to the inner solver is expensive and leads to memory
 * fragmentation. The present class allows to avoid this by offering
 * an interface that other classes can use to allocate and deallocate
 * vectors. Different derived classes then implement different
 * strategies to provide temporary storage vectors to using
 * classes.
 *
 * For example, the PrimitiveVectorMemory class simply allocates and
 * deallocated vectors each time it is asked for a vector. It is an
 * appropriate implementation to use for iterative solvers that are
 * called only once, or very infrequently.
 *
 * On the other hand, the GrowingVectorMemory class never returns
 * memory space to the operating system memory management subsystem
 * during its lifetime; it only marks them as unused and allows them
 * to be reused next time a vector is requested.
 *
 * Yet other classes, when implemented, could provide even other
 * strategies for memory management.
 *
 * @author Guido Kanschat, 1998-2003
*/
template<class VECTOR = dealii::Vector<double> >
class VectorMemory : public Subscriptor
{
  public:

				     /**
				      * Virtual destructor is needed
				      * as there are virtual functions
				      * in this class.
				      */
    virtual ~VectorMemory() {}

				     /**
				      * Return a pointer to a new
				      * vector. The number of elements
				      * or their subdivision into
				      * blocks (if applicable) is
				      * unspecified and users of this
				      * function should reset vectors
				      * to their proper size. The same
				      * holds for the contents of
				      * vectors: they are unspecified.
				      */
    virtual VECTOR* alloc () = 0;
    
				     /**
				      * Return a vector and indicate
				      * that it is not going to be
				      * used any further by the
				      * instance that called alloc()
				      * to get a pointer to it.
				      */
    virtual void free (const VECTOR * const) = 0;

				     /** @addtogroup Exceptions
				      * @{ */
    
				     /**
				      * No more available vectors.
				      */
    DeclException0(ExcNoMoreVectors);
				     /**
				      * Vector was not allocated from
				      * this memory pool.
				      */
    DeclException0(ExcNotAllocatedHere);

				     //@}
/**
 * Pointer to vectors allocated from VectorMemory objects. This
 * pointer is safe in the sense that it automatically calls free()
 * when it is destroyed, thus relieving the user from using vector
 * management functions at all.
 *
 * @author Guido Kanschat, 2009
 */
    class Pointer
    {
      public:
					 /**
					  * Constructor, automatically
					  * allocating a vector from
					  * @p mem.
					  */
	Pointer(VectorMemory<VECTOR>& mem);
					 /**
					  * Destructor, automatically
					  * releasing the vector from
					  * the memory #pool.
					  */
	~Pointer();

					 /**
					  * Conversion to regular pointer.
					  */
	operator VECTOR* () const;
	
					 /**
					  * Dereferencing operator.
					  */
	VECTOR& operator * () const;
	
					 /**
					  * Dereferencing operator.
					  */
	VECTOR * operator -> () const;
      private:
					 /**
					  * The memory pool used.
					  */
	SmartPointer<VectorMemory<VECTOR> > pool;
					 /**
					  * The pointer to the vector.
					  */
	VECTOR* v;
    };    
};



/**
 * Simple memory management. See the documentation of the base class
 * for a description of its purpose.
 *
 * This class allocates and deletes
 * vectors as needed from the global heap, i.e. performs no
 * specially adapted actions for memory management.
 */
template<class VECTOR = dealii::Vector<double> >
class PrimitiveVectorMemory : public VectorMemory<VECTOR>
{
  public:
				     /**
				      * Constructor.
				      */
    PrimitiveVectorMemory () {}

				     /**
				      * Return a pointer to a new
				      * vector. The number of elements
				      * or their subdivision into
				      * blocks (if applicable) is
				      * unspecified and users of this
				      * function should reset vectors
				      * to their proper size. The same
				      * holds for the contents of
				      * vectors: they are unspecified.
				      *
				      * For the present class, calling
				      * this function will allocate a
				      * new vector on the heap and
				      * returning a pointer to it.
				      */
    virtual VECTOR* alloc ()
      {
	return new VECTOR();
      }
    
				     /**
				      * Return a vector and indicate
				      * that it is not going to be
				      * used any further by the
				      * instance that called alloc()
				      * to get a pointer to it.
				      *
				      *
				      * For the present class, this
				      * means that the vector is
				      * returned to the global heap.
				      */
    virtual void free (const VECTOR * const v)
      {
	delete v;
      }
};



/**
 * A pool based memory management class. See the documentation of the
 * base class for a description of its purpose.
 *
 * Each time a vector is requested from this class, it checks if it
 * has one available and returns its address, or allocates a new one
 * on the heap. If a vector is returned, through the free() member
 * function, it doesn't return it to the operating system memory
 * subsystem, but keeps it around unused for later use if alloc() is
 * called again, or until the object is destroyed. The class therefore
 * avoid the overhead of repeatedly allocating memory on the heap if
 * temporary vectors are required and released frequently; on the
 * other hand, it doesn't release once-allocated memory at the
 * earliest possible time and may therefore lead to an increased
 * overall memory consumption.
 *
 * All GrowingVectorMemory objects of the same vector type use the
 * same memory Pool. Therefore, functions can create such a
 * VectorMemory object whenever needed without preformance penalty. A
 * drawback of this policy might be that vectors once allocated are
 * only released at the end of the program run. Nebertheless, the
 * since they are reused, this should be of no concern. Additionally,
 * the destructor of the Pool warns about memory leaks.
 *
 * @note Instantiations for this class are provided for the types
 * Vector and BlockVector with number types <tt>float</tt>,
 * <tt>double</tt>, <tt>long double</tt>, <tt>@<std::complex@<float@>@></tt>,
 * <tt>@<std::complex@<double@>@></tt> and <tt>@<std::complex@<long
 * double@>@></tt>. In order to generate additional instances, it is
 * sufficient to define the Pool variable somewhere by
 * @code
 * #include <lac/vector_memory.h>
 *
 * template GrowingVectorMemory<MyVector>::Pool
 *   GrowingVectorMemory<MyVector>::pool;
 * @endcode
 * 
 * @author Guido Kanschat, 1999, 2007
 */
template<class VECTOR = dealii::Vector<double> >
class GrowingVectorMemory : public VectorMemory<VECTOR>
{
  public:
				     /**
				      * Constructor.  The argument
				      * allows to preallocate a
				      * certain number of vectors. The
				      * default is not to do this.
				      */
    GrowingVectorMemory (const unsigned int initial_size = 0,
			 const bool log_statistics = false);

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
				      * Return a pointer to a new
				      * vector. The number of elements
				      * or their subdivision into
				      * blocks (if applicable) is
				      * unspecified and users of this
				      * function should reset vectors
				      * to their proper size. The same
				      * holds for the contents of
				      * vectors: they are unspecified.
				      */
    virtual VECTOR* alloc ();
    
				     /**
				      * Return a vector and indicate
				      * that it is not going to be
				      * used any further by the
				      * instance that called alloc()
				      * to get a pointer to it.
				      *
				      * For the present class, this
				      * means retaining the vector for
				      * later reuse by the alloc()
				      * method.
				      */
    virtual void free (const VECTOR * const);

				     /**
				      * Memory consumed by this class
				      * and all currently allocated
				      * vectors.
				      */
    unsigned int memory_consumption() const;

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
				      * The class providing the actual
				      * storage for the memory pool.
				      *
				      * This is where the actual
				      * storage for GrowingVectorMemory
				      * is provided. Only one of these
				      * pools is used for each vector
				      * type, thus allocating all
				      * vectors from the same storage.
				      *
				      * @author Guido Kanschat, 2007
				      */
    struct Pool
    {
/// Standard constructor creating an empty pool
	Pool();
/// Destructor. Frees memory and warns about memory leaks
	~Pool();
/// Create data vector; does nothing after first initialization
	void initialize(const unsigned int size);
/// Pointer to the storage object
	std::vector<entry_type>* data;
    };
    
				     /**
				      * Array of allocated vectors.
				      */
    static Pool pool;
    
				     /**
				      * Overall number of
				      * allocations. Only used for
				      * bookkeeping and to generate
				      * output at the end of an
				      * object's lifetime.
				      */
    unsigned int total_alloc;
				     /**
				      * Number of vectors currently
				      * allocated in this object; used
				      * for detecting memory leaks.
				      */
    unsigned int current_alloc;
    
				     /**
				      * A flag controlling the logging
				      * of statistics by the
				      * destructor.
				      */
    bool log_statistics;
    
				     /**
				      * Mutex to synchronise access to
				      * internal data of this object
				      * from multiple threads.
				      */
    static Threads::ThreadMutex mutex;
};

/*@}*/

#ifndef DOXYGEN
/* --------------------- inline functions ---------------------- */


template <typename VECTOR>
inline
VectorMemory<VECTOR>::Pointer::Pointer(VectorMemory<VECTOR>& mem)
		:
		pool(&mem), v(0)
{
  v = pool->alloc();
}


template <typename VECTOR>
inline
VectorMemory<VECTOR>::Pointer::~Pointer()
{
  pool->free(v);
}


template <typename VECTOR>
inline
VectorMemory<VECTOR>::Pointer::operator VECTOR* () const
{
  return v;
}


template <typename VECTOR>
inline
VECTOR & VectorMemory<VECTOR>::Pointer::operator * () const
{
  return *v;
}


template <typename VECTOR>
inline
VECTOR * VectorMemory<VECTOR>::Pointer::operator -> () const
{
  return v;
}


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

#endif // DOXYGEN

DEAL_II_NAMESPACE_CLOSE

#endif
