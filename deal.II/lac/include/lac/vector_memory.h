/*----------------------------   vector_memory.h     ---------------------------*/
/*      $Id$                 */
#ifndef __vector_memory_H
#define __vector_memory_H
/*----------------------------   vector_memory.h     ---------------------------*/


#include <base/subscriptor.h>

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
template<class Vector>
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
    virtual void free(const Vector*) = 0;
};





/**
 * Simple memory management.  This memory class is just made for
 * tests. It just allocates and deletes
 * vectors as needed from the global heap, i.e. performs no
 * specially adapted actions to the purpose of this class.
 */
template<class Vector>
class PrimitiveVectorMemory : public VectorMemory<Vector> {
  public:
				     /**
				      * Constructor.
				      */
    PrimitiveVectorMemory () {};

				     /**
				      * Allocate a vector
				      * from the global heap.
				      */
    virtual Vector* alloc() {
      return new Vector();
    };
    
				     /**
				      * Return a vector to the global heap.
				      */
    virtual void free (const Vector* v) {
      delete v;
    };
};



/*----------------------------   vector_memory.h     ---------------------------*/
/* end of #ifndef __vector_memory_H */
#endif
/*----------------------------   vector_memory.h     ---------------------------*/

