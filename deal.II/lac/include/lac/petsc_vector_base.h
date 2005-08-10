//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2004, 2005 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------
#ifndef __deal2__petsc_vector_base_h
#define __deal2__petsc_vector_base_h


#include <base/config.h>
#include <base/subscriptor.h>
#include <lac/exceptions.h>

#ifdef DEAL_II_USE_PETSC

#include <petscvec.h>

#include <vector>
#include <utility>


                                 // forward declaration
template <typename number> class Vector;


/**
 * A namespace in which wrapper classes for PETSc objects reside.
 *
 * @ingroup PETScWrappers
 * @author Wolfgang Bangerth, 2004
 */
namespace PETScWrappers
{
                                   // forward declaration
  class VectorBase;

  
/**
 * A namespace for internal implementation details of the PETScWrapper
 * members.
 * @ingroup PETScWrappers
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
				      * @ingroup PETScWrappers
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
                                          * This looks like a copy operator,
                                          * but does something different than
                                          * usual. In particular, it does not
                                          * copy the member variables of this
                                          * reference. Rather, it handles the
                                          * situation where we have two
                                          * vectors @p v and @p w, and assign
                                          * elements like in
                                          * <tt>v(i)=w(i)</tt>. Here, both
                                          * left and right hand side of the
                                          * assignment have data type
                                          * VectorReference, but what we
                                          * really mean is to assign the
                                          * vector elements represented by the
                                          * two references. This operator
                                          * implements this operation. Note
                                          * also that this allows us to make
                                          * the assignment operator const.
                                          */
        const VectorReference & operator = (const VectorReference &r) const;
        
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
                                          * Subtract <tt>s</tt> from the
                                          * referenced element of the vector.
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
                                         /**
                                          * Exception
                                          */
        DeclException3 (ExcAccessToNonlocalElement,
                        int, int, int,
                        << "You tried to access element " << arg1
                        << " of a distributed vector, but only elements "
                        << arg2 << " through " << arg3
                        << " are stored locally and can be accessed.");
        
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
 * The interface of this class is modeled after the existing Vector
 * class in deal.II. It has almost the same member functions, and is often
 * exchangable. However, since PETSc only supports a single scalar type
 * (either double, float, or a complex data type), it is not templated, and
 * only works with whatever your PETSc installation has defined the data type
 * @p PetscScalar to.
 *
 * Note that PETSc only guarantees that operations do what you expect if the
 * functions @p VecAssemblyBegin and @p VecAssemblyEnd have been called
 * after vector assembly. Therefore, you need to call Vector::compress()
 * before you actually use the vector.
 *
 * @ingroup PETScWrappers
 * @author Wolfgang Bangerth, 2004
 */
  class VectorBase
  {
    public:
                                       /**
                                        * Declare some of the standard types
                                        * used in all containers. These types
                                        * parallel those in the <tt>C++</tt>
                                        * standard libraries <tt>vector<...></tt>
                                        * class.
                                        */
      typedef PetscScalar               value_type;
      typedef size_t                    size_type;
      typedef internal::VectorReference reference;
      typedef const internal::VectorReference const_reference;

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
                                        * Set all components of the vector to
                                        * the given number @p s. Simply pass
                                        * this down to the individual block
                                        * objects, but we still need to declare
                                        * this function to make the example
                                        * given in the discussion about making
                                        * the constructor explicit work.
                                        *
                                        *
                                        * Since the semantics of assigning a
                                        * scalar to a vector are not
                                        * immediately clear, this operator
                                        * should really only be used if you
                                        * want to set the entire vector to
                                        * zero. This allows the intuitive
                                        * notation <tt>v=0</tt>. Assigning
                                        * other values is deprecated and may
                                        * be disallowed in the future.
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
                                        * Return the global dimension of the
                                        * vector.
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
					*
					* To figure out which elements
					* exactly are stored locally,
					* use local_range().
                                        */
      unsigned int local_size () const;

                                       /**
					* Return a pair of indices
					* indicating which elements of
					* this vector are stored
					* locally. The first number is
					* the index of the first
					* element stored, the second
					* the index of the one past
					* the last one that is stored
					* locally. If this is a
					* sequential vector, then the
					* result will be the pair
					* (0,N), otherwise it will be
					* a pair (i,i+n), where
					* <tt>n=local_size()</tt>.
					*/
      std::pair<unsigned int, unsigned int>
      local_range () const;

				       /**
					* Return whether @p index is
					* in the local range or not,
					* see also local_range().
					*/
      bool in_local_range (const unsigned int index) const;

                                       /**
                                        * Provide access to a given element,
                                        * both read and write.
                                        */
      reference
      operator () (const unsigned int index);

                                       /**
                                        * Provide read-only access to an
                                        * element.
                                        */
      PetscScalar
      operator () (const unsigned int index) const;

                                       /**
                                        * A collective set operation: instead
                                        * of setting individual elements of a
                                        * vector, this function allows to set
                                        * a whole set of elements at once. The
                                        * indices of the elements to be set
                                        * are stated in the first argument,
                                        * the corresponding values in the
                                        * second.
                                        */
      void set (const std::vector<unsigned int> &indices,
                const std::vector<PetscScalar>  &values);
      
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
                                        * Return @p true if the vector has no
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
                                        * Addition of @p s to all
                                        * components. Note that @p s is a
                                        * scalar and not a vector.
                                        */
      void add (const PetscScalar s);
    
                                       /**
                                        * Simple vector addition, equal to the
                                        * <tt>operator +=</tt>.
                                        */
      void add (const VectorBase &V);
    
                                       /**
                                        * Simple addition of a multiple of a
                                        * vector, i.e. <tt>*this += a*V</tt>.
                                        */
      void add (const PetscScalar a, const VectorBase &V);
    
                                       /**
                                        * Multiple addition of scaled vectors,
                                        * i.e. <tt>*this += a*V+b*W</tt>.
                                        */
      void add (const PetscScalar a, const VectorBase &V,
                const PetscScalar b, const VectorBase &W);
    
                                       /**
                                        * Scaling and simple vector addition,
                                        * i.e.
                                        * <tt>*this = s*(*this)+V</tt>.
                                        */
      void sadd (const PetscScalar s,
                 const VectorBase     &V);
    
                                       /**
                                        * Scaling and simple addition, i.e.
                                        * <tt>*this = s*(*this)+a*V</tt>.
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
                                        * <tt>*this = s*(*this)+a*V + b*W + c*X</tt>.
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
                                        * Assignment <tt>*this = a*V</tt>.
                                        */
      void equ (const PetscScalar a, const VectorBase &V);
    
                                       /**
                                        * Assignment <tt>*this = a*V + b*W</tt>.
                                        */
      void equ (const PetscScalar a, const VectorBase &V,
                const PetscScalar b, const VectorBase &W);

                                       /**
                                        * Compute the elementwise ratio of the
                                        * two given vectors, that is let
                                        * <tt>this[i] = a[i]/b[i]</tt>. This is
                                        * useful for example if you want to
                                        * compute the cellwise ratio of true to
                                        * estimated error.
                                        *
                                        * This vector is appropriately
                                        * scaled to hold the result.
                                        *
                                        * If any of the <tt>b[i]</tt> is
                                        * zero, the result is
                                        * undefined. No attempt is made
                                        * to catch such situations.
                                        */
      void ratio (const VectorBase &a,
                  const VectorBase &b);

                                       /**
                                        * Print to a
                                        * stream. @p precision denotes
                                        * the desired precision with
                                        * which values shall be printed,
                                        * @p scientific whether
                                        * scientific notation shall be
                                        * used. If @p across is
                                        * @p true then the vector is
                                        * printed in a line, while if
                                        * @p false then the elements
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
                                        * @p v. One could do this
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
                                        * the @p swap function of all C++
                                        * standard containers. Also,
                                        * there is a global function
                                        * <tt>swap(u,v)</tt> that simply calls
                                        * <tt>u.swap(v)</tt>, again in analogy
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
                                        * variable is @p mutable so that the
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
 * Global function @p swap which overloads the default implementation
 * of the C++ standard library which uses a temporary object. The
 * function simply exchanges the data of the two vectors.
 *
 * @ingroup PETScWrappers
 * @author Wolfgang Bangerth, 2004
 */
  inline
  void swap (VectorBase &u, VectorBase &v)
  {
    u.swap (v);
  }

///@if NoDoc
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
    VectorReference::operator = (const VectorReference &r) const
    {
                                       // as explained in the class
                                       // documentation, this is not the copy
                                       // operator. so simply pass on to the
                                       // "correct" assignment operator
      *this = static_cast<PetscScalar> (r);

      return *this;
    }


    
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

          vector.last_action = VectorBase::LastAction::add;
        }

                                       // we have to do above actions in any
                                       // case to be consistent with the MPI
                                       // communication model (see the
                                       // comments in the documentation of
                                       // PETScWrappers::MPI::Vector), but we
                                       // can save some work if the addend is
                                       // zero
      if (value == 0)
        return *this;
      
                                       // use the PETSc function to add something
      const signed int petsc_i = index;
      const int ierr
        = VecSetValues (vector, 1, &petsc_i, &value, ADD_VALUES);
      AssertThrow (ierr == 0, ExcPETScError(ierr));

      
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

          vector.last_action = VectorBase::LastAction::add;
        }
      
                                       // we have to do above actions in any
                                       // case to be consistent with the MPI
                                       // communication model (see the
                                       // comments in the documentation of
                                       // PETScWrappers::MPI::Vector), but we
                                       // can save some work if the addend is
                                       // zero
      if (value == 0)
        return *this;

                                       // use the PETSc function to add something
      const signed int petsc_i = index;
      const PetscScalar subtractand = -value;
      const int ierr
        = VecSetValues (vector, 1, &petsc_i, &subtractand, ADD_VALUES);
      AssertThrow (ierr == 0, ExcPETScError(ierr));
      
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

          vector.last_action = VectorBase::LastAction::insert;
        }
      
                                       // we have to do above actions in any
                                       // case to be consistent with the MPI
                                       // communication model (see the
                                       // comments in the documentation of
                                       // PETScWrappers::MPI::Vector), but we
                                       // can save some work if the factor is
                                       // one
      if (value == 1.)
        return *this;

      const signed int petsc_i = index;

      const PetscScalar new_value
        = static_cast<PetscScalar>(*this) * value;
      
      const int ierr
        = VecSetValues (vector, 1, &petsc_i, &new_value, INSERT_VALUES);
      AssertThrow (ierr == 0, ExcPETScError(ierr));
      
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

          vector.last_action = VectorBase::LastAction::insert;
        }
      
                                       // we have to do above actions in any
                                       // case to be consistent with the MPI
                                       // communication model (see the
                                       // comments in the documentation of
                                       // PETScWrappers::MPI::Vector), but we
                                       // can save some work if the factor is
                                       // one
      if (value == 1.)
        return *this;

      const signed int petsc_i = index;

      const PetscScalar new_value
        = static_cast<PetscScalar>(*this) / value;
      
      const int ierr
        = VecSetValues (vector, 1, &petsc_i, &new_value, INSERT_VALUES);
      AssertThrow (ierr == 0, ExcPETScError(ierr));
      
      return *this;
    }
  }
  


  inline
  bool
  VectorBase::in_local_range (const unsigned int index) const
  {
    int begin, end;
    const int ierr = VecGetOwnershipRange (static_cast<const Vec &>(vector),
					   &begin, &end);
    AssertThrow (ierr == 0, ExcPETScError(ierr));
    
    return ((index >= static_cast<unsigned int>(begin)) &&
            (index < static_cast<unsigned int>(end)));
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
  

///@endif
}

#endif // DEAL_II_USE_PETSC

/*----------------------------   petsc_vector_base.h     ---------------------------*/

#endif
/*----------------------------   petsc_vector_base.h     ---------------------------*/
