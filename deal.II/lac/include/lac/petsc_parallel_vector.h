//----------------------------  petsc_parallel_vector.h  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2004 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  petsc_parallel_vector.h  ---------------------------
#ifndef __deal2__petsc_parallel_vector_h
#define __deal2__petsc_parallel_vector_h

#include <base/config.h>
#include <base/exceptions.h>
#include <base/subscriptor.h>

#include <lac/vector.h>

#ifdef DEAL_II_USE_PETSC

#include <lac/petsc_vector_base.h>


namespace PETScWrappers
{
/**
 * Namespace for PETSc classes that work in parallel over MPI, such as
 * distributed vectors and matrices.
 *
 * @author Wolfgang Bangerth, 2004
 */
  namespace MPI
  {
    
/**
 * Implementation of a parallel vector class based on PETSC and using MPI
 * communication to synchronise distributed operations. All the functionality
 * is actually in the base class, except for the calls to generate a
 * parallel vector. This is possible since PETSc only works on an abstract
 * vector type and internally distributes to functions that do the actual work
 * depending on the actual vector type (much like using virtual
 * functions). Only the functions creating a vector of specific type differ,
 * and are implemented in this particular class.
 *
 * @author Wolfgang Bangerth, 2004
 */
    class Vector : public VectorBase
    {
      public:
                                         /**
                                          * Default constructor. Initialize the
                                          * vector as empty.
                                          */
        Vector ();
      
                                         /**
                                          * Constructor. Set dimension to
                                          * @p{n} and initialize all
                                          * elements with zero.
                                          *
                                          * The constructor is made explicit to
                                          * avoid accidents like this:
                                          * @p{v=0;}. Presumably, the user wants
                                          * to set every element of the vector to
                                          * zero, but instead, what happens is
                                          * this call: @p{v=Vector<number>(0);},
                                          * i.e. the vector is replaced by one of
                                          * length zero.
                                          */
        explicit Vector (const unsigned int n,
                         const unsigned int local_size,
                         const MPI_Comm    &communicator);
    
                                         /**
                                          * Copy-constructor from deal.II
                                          * vectors. Sets the dimension to that
                                          * of the given vector, and copies all
                                          * elements.
                                          */
        template <typename Number>
        explicit Vector (const ::Vector<Number> &v,
                         const unsigned int      local_size,
                         const MPI_Comm         &communicator);

                                         /**
                                          * Copy-constructor the
                                          * values from a PETSc wrapper vector
                                          * class.
                                          */
        explicit Vector (const VectorBase &v,
                         const unsigned int local_size,
                         const MPI_Comm    &communicator);

                                         /**
                                          * Set all components of the vector to
                                          * the given number @p{s}. Simply pass
                                          * this down to the base class, but we
                                          * still need to declare this function
                                          * to make the example given in the
                                          * discussion about making the
                                          * constructor explicit work.
                                          */
        Vector & operator = (const PetscScalar s);

                                         /**
                                          * Copy the values of a deal.II vector
                                          * (as opposed to those of the PETSc
                                          * vector wrapper class) into this
                                          * object.
                                          */
        template <typename number>
        Vector & operator = (const ::Vector<number> &v);
      
      protected:
                                         /**
                                          * Create a vector of length @p{n}. For
                                          * this class, we create a parallel
                                          * vector. @arg n denotes the total
                                          * size of the vector to be
                                          * created. @arg local_size denotes how
                                          * many of these elements shall be
                                          * stored locally. The last argument is
                                          * ignored for sequential vectors.
                                          */
        virtual void create_vector (const unsigned int n,
                                    const unsigned int local_size);

      private:
                                         /**
                                          * Copy of the communicator object to
                                          * be used for this parallel vector.
                                          */
        MPI_Comm communicator;
    };



// ------------------ template and inline functions -------------


    template <typename number>
    Vector::Vector (const ::Vector<number> &v,
                    const unsigned int      local_size,
                    const MPI_Comm         &communicator)
                    :
                    communicator (communicator)
    {
      Vector::create_vector (v.size(), local_size);

      VectorBase::operator = (v);
    }

  
  
    inline
    Vector &
    Vector::operator = (const PetscScalar s)
    {
      VectorBase::operator = (s);

      return *this;
    }
  

    template <typename number>
    inline
    Vector &
    Vector::operator = (const ::Vector<number> &v)
    {
      VectorBase::operator = (v);

      return *this;
    }
    
  }
  
}


#endif // DEAL_II_USE_PETSC

/*----------------------------   petsc_parallel_vector.h     ---------------------------*/

#endif
/*----------------------------   petsc_parallel_vector.h     ---------------------------*/
