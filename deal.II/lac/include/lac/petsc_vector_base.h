//----------------------------  petsc_vector_base.h  ---------------------------
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
//----------------------------  petsc_vector_base.h  ---------------------------
#ifndef __deal2__petsc_vector_base_h
#define __deal2__petsc_vector_base_h


#include <base/config.h>
#include <base/exceptions.h>
#include <base/subscriptor.h>

#ifdef DEAL_II_USE_PETSC

#include <petscvec.h>


                                 // forward declaration
template <typename number> class Vector;


/**
 * A namespace in which wrapper classes for PETSc objects reside.
 *
 * @author Wolfgang Bangerth, 2004
 */
namespace PETScWrappers
{
                                   // forward declaration
  class VectorBase;

  
/**
 * A namespace for internal implementation details of the PETScWrapper
 * members.
 */
  namespace internal
  {
                                     /**
                                      * Since access to PETSc vectors only
                                      * goes through functions, rather than by
                                      * obtaining a reference to a vector
                                      * element, we need a wrapper class that
                                      * acts as if it was a reference, and
                                      * basically redirects all accesses (read
                                      * and write) to member functions of this
                                      * class.
                                      *
                                      * This class implements such a wrapper:
                                      * it is initialized with a vector and an
                                      * element within it, and has a
                                      * conversion operator to extract the
                                      * scalar value of this element. It also
                                      * has a variety of assignment operator
                                      * for writing to this one element.
                                      */
    class VectorReference
    {
      private:
                                         /**
                                          * Constructor. It is made private so
                                          * as to only allow the actual vector
                                          * class to create it.
                                          */
        VectorReference (const VectorBase  &vector,
                         const unsigned int index);
        
      public:
                                         /**
                                          * Set the referenced element of the
                                          * vector to <tt>s</tt>.
                                          */
        const VectorReference & operator = (const PetscScalar &s) const;

                                         /**
                                          * Add <tt>s</tt> to the referenced
                                          * element of the vector.
                                          */
        const VectorReference & operator += (const PetscScalar &s) const;

                                         /**
                                          * Subtract <tt>s</tt> from the referenced
                                          * element of the vector.
                                          */
        const VectorReference & operator -= (const PetscScalar &s) const;

                                         /**
                                          * Multiply the referenced element of
                                          * the vector by <tt>s</tt>.
                                          */
        const VectorReference & operator *= (const PetscScalar &s) const;

                                         /**
                                          * Divide the referenced element of
                                          * the vector by <tt>s</tt>.
                                          */
        const VectorReference & operator /= (const PetscScalar &s) const;

                                         /**
                                          * Convert the reference to an actual
                                          * value, i.e. return the value of
                                          * the referenced element of the
                                          * vector.
                                          */
        operator PetscScalar () const;

                                         /**
                                          * Exception
                                          */
        DeclException1 (ExcPETScError,
                        int,
                        << "An error with error number " << arg1
                        << " occured while calling a PETSc function");
      private:
                                         /**
                                          * Point to the vector we are
                                          * referencing.
                                          */
        const VectorBase   &vector;

                                         /**
                                          * Index of the referenced element of
                                          * the vector.
                                          */
        const unsigned int  index;

                                         /**
                                          * Make the vector class a friend, so
                                          * that it can create objects of the
                                          * present type.
                                          */
        friend class ::PETScWrappers::VectorBase;
    };
  }
  

/**
 * Base class for all vector classes that are implemented on top of the PETSc
 * vector types. Since in PETSc all vector types (i.e. sequential and parallel
 * ones) are built by filling the contents of an abstract object that is only
 * referenced through a pointer of a type that is independent of the actual
 * vector type, we can implement almost all functionality of vectors in this
 * base class. Derived classes will then only have to provide the
 * functionality to create one or the other kind of vector.
 *
 * The interface of this class is modeled after the existing @ref{Vector}
 * class in deal.II. It has almost the same member functions, and is often
 * exchangable. However, since PETSc only supports a single scalar type
 * (either double, float, or a complex data type), it is not templated, and
 * only works with whatever your PETSc installation has defined the data type
 * @p{PetscScalar} to.
 *
 * Note that PETSc only guarantees that operations do what you expect if the
 * functions @p{VecAssemblyBegin} and @p{VecAssemblyEnd} have been called
 * after vector assembly. Therefore, you need to call @ref{Vector::compress}
 * before you actually use the vector.
 *
 * @author Wolfgang Bangerth, 2004
 */
  class VectorBase
  {
    public:
                                       /**
                                        * Declare some of the standard types
                                        * used in all containers. These types
                                        * parallel those in the @p{C++}
                                        * standard libraries @p{vector<...>}
                                        * class.
                                        */
      typedef PetscScalar       value_type;
      typedef size_t            size_type;

                                       /**
                                        * Default constructor. It doesn't do
                                        * anything, derived classes will have
                                        * to initialize the data.
                                        */
      VectorBase ();
      
                                       /**
                                        * Copy constructor. Sets the dimension
                                        * to that of the given vector, and
                                        * copies all elements.
                                        */
      VectorBase (const VectorBase &v);

                                       /**
                                        * Destructor
                                        */
      virtual ~VectorBase ();

                                       /**
                                        * Compress the underlying
                                        * representation of the PETSc object,
                                        * i.e. flush the buffers of the vector
                                        * object if it has any. This function
                                        * is necessary after writing into a
                                        * vector element-by-element and before
                                        * anything else can be done on it.
                                        */
      void compress ();

                                       /**
                                        * Set all entries to zero. Equivalent
                                        * to @p{v = 0}, but more obvious and
                                        * faster.  Note that this function
                                        * does not change the size of the
                                        * vector, unlike the STL's
                                        * @p{vector<>::clear} function.
                                        */
      void clear ();    
      
                                       /**
                                        * Set all components of the vector to
                                        * the given number @p{s}.
                                        */
      VectorBase & operator = (const PetscScalar s);    
      
                                       /**
                                        * Test for equality. This function
                                        * assumes that the present vector and
                                        * the one to compare with have the same
                                        * size already, since comparing vectors
                                        * of different sizes makes not much
                                        * sense anyway.
                                        */
      bool operator == (const VectorBase &v) const;
    
                                       /**
                                        * Test for inequality. This function
                                        * assumes that the present vector and
                                        * the one to compare with have the same
                                        * size already, since comparing vectors
                                        * of different sizes makes not much
                                        * sense anyway.
                                        */
      bool operator != (const VectorBase &v) const;

                                       /**
                                        * Return the global dimension of the vector.
                                        */
      unsigned int size () const;

                                       /**
                                        * Return the local dimension of the
                                        * vector, i.e. the number of elements
                                        * stored on the present MPI
                                        * process. For sequential vectors,
                                        * this number is the same as size(),
                                        * but for parallel vectors it may be
                                        * smaller.
                                        */
      unsigned int local_size () const;
      
                                       /**
                                        * Provide access to a given element,
                                        * both read and write.
                                        */
      internal::VectorReference
      operator () (const unsigned int index);

                                       /**
                                        * Provide read-only access to an
                                        * element.
                                        */
      PetscScalar
      operator () (const unsigned int index) const;
      
                                       /**
                                        * Return the scalar product of two
                                        * vectors. The vectors must have the
                                        * same size.
                                        */
      PetscScalar operator * (const VectorBase &vec) const;

                                       /**
                                        * Return square of the $l_2$-norm.
                                        */
      PetscScalar norm_sqr () const;

                                       /**
                                        * Mean value of the elements of
                                        * this vector.
                                        */
      PetscScalar mean_value () const;

                                       /**
                                        * $l_1$-norm of the vector.
                                        * The sum of the absolute values.
                                        */
      PetscScalar l1_norm () const;

                                       /**
                                        * $l_2$-norm of the vector.  The
                                        * square root of the sum of the
                                        * squares of the elements.
                                        */
      PetscScalar l2_norm () const;

                                       /**
                                        * $l_p$-norm of the vector. The
                                        * pth root of the sum of the pth
                                        * powers of the absolute values
                                        * of the elements.
                                        */
      PetscScalar lp_norm (const PetscScalar p) const;

                                       /**
                                        * Maximum absolute value of the
                                        * elements.
                                        */
      PetscScalar linfty_norm () const;

                                       /**
                                        * Return whether the vector contains
                                        * only elements with value zero. This
                                        * function is mainly for internal
                                        * consistency checks and should
                                        * seldomly be used when not in debug
                                        * mode since it uses quite some time.
                                        */
      bool all_zero () const;

                                       /**
                                        * Return @p{true} if the vector has no
                                        * negative entries, i.e. all entries
                                        * are zero or positive. This function
                                        * is used, for example, to check
                                        * whether refinement indicators are
                                        * really all positive (or zero).
                                        */
      bool is_non_negative () const;
      
                                       /**
                                        * Multiply the entire vector by a
                                        * fixed factor.
                                        */
      VectorBase & operator *= (const PetscScalar factor);
    
                                       /**
                                        * Divide the entire vector by a
                                        * fixed factor.
                                        */
      VectorBase & operator /= (const PetscScalar factor);

                                       /**
                                        * Add the given vector to the present
                                        * one.
                                        */
      VectorBase & operator += (const VectorBase &V);

                                       /**
                                        * Subtract the given vector from the
                                        * present one.
                                        */
      VectorBase & operator -= (const VectorBase &V);

                                       /**
                                        * Addition of @p{s} to all
                                        * components. Note that @p{s} is a
                                        * scalar and not a vector.
                                        */
      void add (const PetscScalar s);
    
                                       /**
                                        * Simple vector addition, equal to the
                                        * @p{operator +=}.
                                        */
      void add (const VectorBase &V);
    
                                       /**
                                        * Simple addition of a multiple of a
                                        * vector, i.e. @p{*this += a*V}.
                                        */
      void add (const PetscScalar a, const VectorBase &V);
    
                                       /**
                                        * Multiple addition of scaled vectors,
                                        * i.e. @p{*this += a*V+b*W}.
                                        */
      void add (const PetscScalar a, const VectorBase &V,
                const PetscScalar b, const VectorBase &W);
    
                                       /**
                                        * Scaling and simple vector addition,
                                        * i.e.
                                        * @p{*this = s*(*this)+V}.
                                        */
      void sadd (const PetscScalar s,
                 const VectorBase     &V);
    
                                       /**
                                        * Scaling and simple addition, i.e.
                                        * @p{*this = s*(*this)+a*V}.
                                        */
      void sadd (const PetscScalar s,
                 const PetscScalar a,
                 const VectorBase     &V);
    
                                       /**
                                        * Scaling and multiple addition.
                                        */
      void sadd (const PetscScalar s,
                 const PetscScalar a,
                 const VectorBase     &V,
                 const PetscScalar b,
                 const VectorBase     &W);
    
                                       /**
                                        * Scaling and multiple addition.
                                        * @p{*this = s*(*this)+a*V + b*W + c*X}.
                                        */
      void sadd (const PetscScalar s,
                 const PetscScalar a,
                 const VectorBase     &V,
                 const PetscScalar b,
                 const VectorBase     &W, 
                 const PetscScalar c,
                 const VectorBase     &X);
    
                                       /**
                                        * Scale each element of this
                                        * vector by the corresponding
                                        * element in the argument. This
                                        * function is mostly meant to
                                        * simulate multiplication (and
                                        * immediate re-assignment) by a
                                        * diagonal scaling matrix.
                                        */
      void scale (const VectorBase &scaling_factors);
    
                                       /**
                                        * Assignment @p{*this = a*V}.
                                        */
      void equ (const PetscScalar a, const VectorBase &V);
    
                                       /**
                                        * Assignment @p{*this = a*V + b*W}.
                                        */
      void equ (const PetscScalar a, const VectorBase &V,
                const PetscScalar b, const VectorBase &W);

                                       /**
                                        * Compute the elementwise ratio of the
                                        * two given vectors, that is let
                                        * @p{this[i] = a[i]/b[i]}. This is
                                        * useful for example if you want to
                                        * compute the cellwise ratio of true to
                                        * estimated error.
                                        *
                                        * This vector is appropriately
                                        * scaled to hold the result.
                                        *
                                        * If any of the @p{b[i]} is
                                        * zero, the result is
                                        * undefined. No attempt is made
                                        * to catch such situations.
                                        */
      void ratio (const VectorBase &a,
                  const VectorBase &b);

                                       /**
                                        * Print to a
                                        * stream. @p{precision} denotes
                                        * the desired precision with
                                        * which values shall be printed,
                                        * @p{scientific} whether
                                        * scientific notation shall be
                                        * used. If @p{across} is
                                        * @p{true} then the vector is
                                        * printed in a line, while if
                                        * @p{false} then the elements
                                        * are printed on a separate line
                                        * each.
                                        */
      void print (std::ostream       &out,
                  const unsigned int  precision  = 3,
                  const bool          scientific = true,
                  const bool          across     = true) const;

                                       /**
                                        * Swap the contents of this
                                        * vector and the other vector
                                        * @p{v}. One could do this
                                        * operation with a temporary
                                        * variable and copying over the
                                        * data elements, but this
                                        * function is significantly more
                                        * efficient since it only swaps
                                        * the pointers to the data of
                                        * the two vectors and therefore
                                        * does not need to allocate
                                        * temporary storage and move
                                        * data around.
                                        *
                                        * This function is analog to the
                                        * the @p{swap} function of all C++
                                        * standard containers. Also,
                                        * there is a global function
                                        * @p{swap(u,v)} that simply calls
                                        * @p{u.swap(v)}, again in analogy
                                        * to standard functions.
                                        */
      void swap (VectorBase &v);
      
      
                                       /**
                                        * Conversion operator to gain access
                                        * to the underlying PETSc type. If you
                                        * do this, you cut this class off some
                                        * information it may need, so this
                                        * conversion operator should only be
                                        * used if you know what you do. In
                                        * particular, it should only be used
                                        * for read-only operations into the
                                        * vector.
                                        */
      operator const Vec & () const;
      
                                       /**
                                        * Exception
                                        */
      DeclException1 (ExcPETScError,
                      int,
                      << "An error with error number " << arg1
                      << " occured while calling a PETSc function");
                                       /**
                                        * Exception
                                        */
      DeclException2 (ExcNonMatchingSizes,
                      int, int,
                      << "The sizes " << arg1 << " and " << arg2
                      << " are supposed to be equal, but are not.");
      
    protected:
                                       /**
                                        * A generic vector object in
                                        * PETSc. The actual type, a sequential
                                        * vector, is set in the constructor.
                                        */
      Vec vector;


                                       /**
                                        * PETSc doesn't allow to mix additions
                                        * to matrix entries and overwriting
                                        * them (to make synchronisation of
                                        * parallel computations
                                        * simpler). Since the interface of the
                                        * existing classes don't support the
                                        * notion of not interleaving things,
                                        * we have to emulate this
                                        * ourselves. The way we do it is to,
                                        * for each access operation, store
                                        * whether it is an insertion or an
                                        * addition. If the previous one was of
                                        * different type, then we first have
                                        * to flush the PETSc buffers;
                                        * otherwise, we can simply go on.
                                        *
                                        * The following structure and variable
                                        * declare and store the previous
                                        * state.
                                        */
      struct LastAction
      {
          enum Values { none, insert, add };
      };


                                       /**
                                        * Store whether the last action was a
                                        * write or add operation. This
                                        * variable is @p{mutable} so that the
                                        * accessor classes can write to it,
                                        * even though the vector object they
                                        * refer to is constant.
                                        */
      mutable LastAction::Values last_action;

                                       /**
                                        * Make the reference class a friend.
                                        */
      friend class internal::VectorReference;
  };



// ------------------- inline and template functions --------------  

/**
 * Global function @p{swap} which overloads the default implementation
 * of the C++ standard library which uses a temporary object. The
 * function simply exchanges the data of the two vectors.
 *
 * @author Wolfgang Bangerth, 2004
 */
  inline
  void swap (VectorBase &u, VectorBase &v)
  {
    u.swap (v);
  }


  namespace internal
  {
    inline
    VectorReference::VectorReference (const VectorBase  &vector,
                                      const unsigned int index)
                    :
                    vector (vector),
                    index (index)
    {}


    inline
    const VectorReference &
    VectorReference::operator = (const PetscScalar &value) const
    {
      if (vector.last_action != VectorBase::LastAction::insert)
        {
          int ierr;
          ierr = VecAssemblyBegin (vector);
          AssertThrow (ierr == 0, ExcPETScError(ierr));

          ierr = VecAssemblyEnd (vector);
          AssertThrow (ierr == 0, ExcPETScError(ierr));
        }
      
      const signed int petsc_i = index;
      
      const int ierr
        = VecSetValues (vector, 1, &petsc_i, &value, INSERT_VALUES);
      AssertThrow (ierr == 0, ExcPETScError(ierr));

      vector.last_action = VectorBase::LastAction::insert;
      
      return *this;
    }



    inline
    const VectorReference &
    VectorReference::operator += (const PetscScalar &value) const
    {
      if (vector.last_action != VectorBase::LastAction::add)
        {
          int ierr;
          ierr = VecAssemblyBegin (vector);
          AssertThrow (ierr == 0, ExcPETScError(ierr));

          ierr = VecAssemblyEnd (vector);
          AssertThrow (ierr == 0, ExcPETScError(ierr));
        }
      
      const signed int petsc_i = index;

                                       // use the PETSc function to add something
      const int ierr
        = VecSetValues (vector, 1, &petsc_i, &value, ADD_VALUES);
      AssertThrow (ierr == 0, ExcPETScError(ierr));

      vector.last_action = VectorBase::LastAction::add;
      
      return *this;
    }



    inline
    const VectorReference &
    VectorReference::operator -= (const PetscScalar &value) const
    {
      if (vector.last_action != VectorBase::LastAction::add)
        {
          int ierr;
          ierr = VecAssemblyBegin (vector);
          AssertThrow (ierr == 0, ExcPETScError(ierr));

          ierr = VecAssemblyEnd (vector);
          AssertThrow (ierr == 0, ExcPETScError(ierr));
        }
      
      const signed int petsc_i = index;

                                       // use the PETSc function to add something
      const PetscScalar subtractand = -value;
      const int ierr
        = VecSetValues (vector, 1, &petsc_i, &subtractand, ADD_VALUES);
      AssertThrow (ierr == 0, ExcPETScError(ierr));

      vector.last_action = VectorBase::LastAction::add;
      
      return *this;
    }



    inline
    const VectorReference &
    VectorReference::operator *= (const PetscScalar &value) const
    {
      if (vector.last_action != VectorBase::LastAction::insert)
        {
          int ierr;
          ierr = VecAssemblyBegin (vector);
          AssertThrow (ierr == 0, ExcPETScError(ierr));

          ierr = VecAssemblyEnd (vector);
          AssertThrow (ierr == 0, ExcPETScError(ierr));
        }
      
      const signed int petsc_i = index;

      const PetscScalar new_value
        = static_cast<PetscScalar>(*this) * value;
      
      const int ierr
        = VecSetValues (vector, 1, &petsc_i, &new_value, INSERT_VALUES);
      AssertThrow (ierr == 0, ExcPETScError(ierr));

      vector.last_action = VectorBase::LastAction::insert;
      
      return *this;
    }



    inline
    const VectorReference &
    VectorReference::operator /= (const PetscScalar &value) const
    {
      if (vector.last_action != VectorBase::LastAction::insert)
        {
          int ierr;
          ierr = VecAssemblyBegin (vector);
          AssertThrow (ierr == 0, ExcPETScError(ierr));

          ierr = VecAssemblyEnd (vector);
          AssertThrow (ierr == 0, ExcPETScError(ierr));
        }
      
      const signed int petsc_i = index;

      const PetscScalar new_value
        = static_cast<PetscScalar>(*this) / value;
      
      const int ierr
        = VecSetValues (vector, 1, &petsc_i, &new_value, INSERT_VALUES);
      AssertThrow (ierr == 0, ExcPETScError(ierr));

      vector.last_action = VectorBase::LastAction::insert;
      
      return *this;
    }



    inline
    VectorReference::operator PetscScalar () const
    {
                                       // this is clumsy: there is no simple
                                       // way in PETSc read an element from a
                                       // vector, i.e. there is no function
                                       // VecGetValue or so. The only way is
                                       // to obtain a pointer to a contiguous
                                       // representation of the vector and
                                       // read from it. Subsequently, the
                                       // vector representation has to be
                                       // restored. If the vector has some
                                       // kind of non-standard format, such as
                                       // for parallel vectors, then this is a
                                       // costly operation, just for a single
                                       // read access..
      PetscScalar *ptr;
      int ierr
        = VecGetArray (static_cast<const Vec &>(vector), &ptr);
      AssertThrow (ierr == 0, ExcPETScError(ierr));

      const PetscScalar value = *(ptr+index);

      ierr = VecRestoreArray (static_cast<const Vec &>(vector), &ptr);
      AssertThrow (ierr == 0, ExcPETScError(ierr));
      
      return value;
    }
  }


  inline
  internal::VectorReference
  VectorBase::operator () (const unsigned int index)
  {
    return internal::VectorReference (*this, index);
  }



  inline
  PetscScalar
  VectorBase::operator () (const unsigned int index) const
  {
    return static_cast<PetscScalar>(internal::VectorReference (*this, index));
  }
  


}

#endif // DEAL_II_USE_PETSC

/*----------------------------   petsc_vector_base.h     ---------------------------*/

#endif
/*----------------------------   petsc_vector_base.h     ---------------------------*/
