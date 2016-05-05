// ---------------------------------------------------------------------
//
// Copyright (C) 2015 - 2016 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------

#ifndef dealii__read_write_vector_h
#define dealii__read_write_vector_h

#include <deal.II/base/config.h>
#include <deal.II/base/index_set.h>
#include <deal.II/base/mpi.h>
#include <deal.II/base/template_constraints.h>
#include <deal.II/base/types.h>
#include <deal.II/base/utilities.h>
#include <deal.II/base/memory_consumption.h>
#include <deal.II/base/thread_management.h>
#include <deal.II/lac/vector_view.h>

#include <cstring>
#include <iomanip>

#ifdef DEAL_II_WITH_TRILINOS
#include <deal.II/lac/trilinos_epetra_communication_pattern.h>
#include "Epetra_MultiVector.h"
#endif

DEAL_II_NAMESPACE_OPEN

namespace LinearAlgebra
{
  class CommunicationPatternBase;
}

#ifdef DEAL_II_WITH_PETSC
namespace PETScWrappers
{
  namespace MPI
  {
    class Vector;
  }
}
#endif

#ifdef DEAL_II_WITH_TRILINOS
namespace TrilinosWrappers
{
  namespace MPI
  {
    class Vector;
  }
}

namespace LinearAlgebra
{
  namespace EpetraWrappers
  {
    class Vector;
  }
}
#endif

namespace LinearAlgebra
{
  /*! @addtogroup Vectors
   *@{
   */

  /**
   * ReadWriteVector is intended to represent vectors in ${\mathbb R}^N$ for
   * which it stores all or a subset of elements. The latter case in important
   * in parallel computations, where $N$ may be so large that no processor can
   * actually all elements of a solution vector, but where this is also not
   * necessary: one typically only has to store the values of degrees of
   * freedom that live on cells that are locally owned plus potentially those
   * degrees of freedom that live on ghost cells.
   *
   * This class allows to access individual elements to be read or written.
   * However, it does not allow global operations such as taking the norm.
   * ReadWriteVector can be used to read and write elements in vectors derived
   * from VectorSpaceVector such as TrilinosWrappers::MPI::Vector and
   * PETScWrappers::MPI::Vector.
   *
   * <h3>Storing elements</h3> Most of the time, one will simply read from or
   * write into a vector of the current class using the global numbers of
   * these degrees of freedom. This is done using operator() or operator[]
   * which call global_to_local() to transform the <i>global</i> index into a
   * <i>local</i> one. In such cases, it is clear that one can only access
   * elements of the vector that the current object indeed stores.
   *
   * However, it is also possible to access elements in the order in which
   * they are stored by the current object. In other words, one is not
   * interested in accessing elements with their <i>global</i> indices, but
   * instead using an enumeration that only takes into account the elements
   * that are actually stored. This is facilitated by the local_element()
   * function. To this end, it is necessary to know <i>in which order</i> the
   * current class stores its element. The elements of all the consecutive
   * ranges are stored in ascending order of the first index of each range.
   * The function largest_range_starting_index() of IndexSet can be used to
   * get the first index of the largest range.
   *
   * @author Bruno Turcksin, 2015.
   */
  template <typename Number>
  class ReadWriteVector : public Subscriptor
  {
  public:
    /**
     * Declare standard types used in all containers. These types parallel
     * those in the <tt>C++</tt> standard libraries <tt>vector<...></tt>
     * class.
     */
    typedef Number                                            value_type;
    typedef value_type                                       *pointer;
    typedef const value_type                                 *const_pointer;
    typedef value_type                                       *iterator;
    typedef const value_type                                 *const_iterator;
    typedef value_type                                       &reference;
    typedef const value_type                                 &const_reference;
    typedef types::global_dof_index                           size_type;
    typedef typename numbers::NumberTraits<Number>::real_type real_type;

    /**
     * @name 1: Basic Object-handling
     */
    //@{
    /**
     * Empty constructor.
     */
    ReadWriteVector ();

    /**
     * Copy constructor.
     */
    ReadWriteVector (const ReadWriteVector<Number> &in_vector);

    /**
     * Constructs a vector given the size, the stored elements have their
     * index in [0,size).
     */
    explicit ReadWriteVector (const size_type size);

    /**
     * Constructs a vector whose stored elements indices are given by the
     * IndexSet @p locally_stored_indices.
     */
    explicit ReadWriteVector (const IndexSet &locally_stored_indices);

    /**
     * Destructor.
     */
    ~ReadWriteVector ();

    /**
     * Sets the global size of the vector to @p size. The stored elements have
     * their index in [0,size).
     *
     * If the flag @p omit_zeroing_entries is set to false, the memory will be
     * initialized with zero, otherwise the memory will be untouched (and the
     * user must make sure to fill it with reasonable data before using it).
     */
    void reinit (const size_type size,
                 const bool      omit_zeroing_entries = false);

    /**
     * Uses the same IndexSet as the one of the input vector @p in_vector and
     * allocates memory for this vector.
     *
     * If the flag @p omit_zeroing_entries is set to false, the memory will be
     * initialized with zero, otherwise the memory will be untouched (and the
     * user must make sure to fill it with reasonable data before using it).
     */
    template <typename Number2>
    void reinit(const ReadWriteVector<Number2> &in_vector,
                const bool                      omit_zeroing_entries = false);

    /**
     * Initializes the vector. The indices are specified by @p
     * locally_stored_indices.
     *
     * If the flag @p omit_zeroing_entries is set to false, the memory will be
     * initialized with zero, otherwise the memory will be untouched (and the
     * user must make sure to fill it with reasonable data before using it).
     * locally_stored_indices.
     */
    void reinit (const IndexSet &locally_stored_indices,
                 const bool      omit_zeroing_entries = false);

    /**
     * Swap the contents of this vector and the other vector @p v. One could
     * do this operation with a temporary variable and copying over the data
     * elements, but this function is significantly more efficient since it
     * only swaps the pointers to the data of the two vectors and therefore
     * does not need to allocate temporary storage and move data around.
     *
     * This function is analog to the the @p swap function of all C++ standard
     * containers. Also, there is a global function <tt>swap(u,v)</tt> that
     * simply calls <tt>u.swap(v)</tt>, again in analogy to standard
     * functions.
     */
    void swap (ReadWriteVector<Number> &v);

    /**
     * Copies the data and the IndexSet of the input vector @p in_vector.
     */
    ReadWriteVector<Number> &
    operator= (const ReadWriteVector<Number> &in_vector);

    /**
     * Copies the data and the IndexSet of the input vector @p in_vector.
     */
    template <typename Number2>
    ReadWriteVector<Number> &
    operator= (const ReadWriteVector<Number2> &in_vector);

    /**
     * Sets all elements of the vector to the scalar @p s. This operation is
     * only allowed if @p s is equal to zero.
     */
    ReadWriteVector<Number> &operator = (const Number s);

#ifdef DEAL_II_WITH_PETSC
    /**
     * Imports all the elements present in the vector's IndexSet from the input
     * vector @p petsc_vec. VectorOperation::values @p operation is used to decide
     * if the elements in @p V should be added to the current vector or replace
     * the current elements. The last parameter can be used if the same
     * communication pattern is used multiple times. This can be used to improve
     * performance.
     */
    void import(const PETScWrappers::MPI::Vector &petsc_vec,
                VectorOperation::values operation,
                std_cxx11::shared_ptr<const CommunicationPatternBase> communication_pattern =
                  std_cxx11::shared_ptr<const CommunicationPatternBase> ());
#endif

#ifdef DEAL_II_WITH_TRILINOS
    /**
     * Imports all the elements present in the vector's IndexSet from the input
     * vector @p trilinos_vec. VectorOperation::values @p operation is used to
     * decide if the elements in @p V should be added to the current vector or
     * replace the current elements. The last parameter can be used if the same
     * communication pattern is used multiple times. This can be used to improve
     * performance.
     */
    void import(const TrilinosWrappers::MPI::Vector &trilinos_vec,
                VectorOperation::values operation,
                std_cxx11::shared_ptr<const CommunicationPatternBase> communication_pattern =
                  std_cxx11::shared_ptr<const CommunicationPatternBase> ());

    /**
     * Imports all the elements present in the vector's IndexSet from the input
     * vector @p epetra_vec. VectorOperation::values @p operation is used to
     * decide if the elements in @p V should be added to the current vector or
     * replace the current elements. The last parameter can be used if the same
     * communication pattern is used multiple times. This can be used to improve
     * performance.
     */
    void import(const EpetraWrappers::Vector &epetra_vec,
                VectorOperation::values operation,
                std_cxx11::shared_ptr<const CommunicationPatternBase> communication_pattern =
                  std_cxx11::shared_ptr<const CommunicationPatternBase> ());
#endif

    /**
     * The value returned by this function denotes the dimension of the vector
     * spaces that are modeled by objects of this kind. However, objects of
     * the current class do not actually stores all elements of vectors of
     * this space but may, in fact store only a subset. The number of elements
     * stored is returned by n_elements() and is smaller or equal to the
     * number returned by the current function.
     */
    size_type size() const;

    /**
     * This function returns the number of elements stored. It is smaller or
     * equal to the dimension of the vector space that is modeled by an object
     * of this kind. This dimension is return by size().
     */
    size_type n_elements() const;

    /**
     * Return the IndexSet that represents the indices of the elements stored.
     */
    const IndexSet &get_stored_elements () const;

    /**
     * Make the @p ReadWriteVector class a bit like the <tt>vector<></tt>
     * class of the C++ standard library by returning iterators to the start
     * and end of the <i>locally stored</i> elements of this vector.
     */
    iterator begin ();

    /**
     * Returns constant iterator to the start of the locally stored elements
     * of the vector.
     */
    const_iterator begin () const;

    /**
     * Returns an iterator pointing to the element past the end of the array
     * of locally stored entries.
     */
    iterator end ();

    /**
     * Returns a constant iterator pointing to the element past the end of the
     * array of the locally stored entries.
     */
    const_iterator end () const;
    //@}


    /**
     * @name 2: Data-Access
     */
    //@{

    /**
     * Read access to the data in the position corresponding to @p
     * global_index. An exception is thrown if @p global_index is not stored
     * by the current object.
     */
    Number operator () (const size_type global_index) const;

    /**
     * Read and write access to the data in the position corresponding to @p
     * global_index. An exception is thrown if @p global_index is not stored
     * by the current object.
     */
    Number &operator () (const size_type global_index);

    /**
     * Read access to the data in the position corresponding to @p
     * global_index. An exception is thrown if @p global_index is not stored
     * by the current object.
     *
     * This function does the same thing as operator().
     */
    Number operator [] (const size_type global_index) const;

    /**
     * Read and write access to the data in the position corresponding to @p
     * global_index. An exception is thrown if @p global_index is not stored
     * by the current object.
     *
     * This function does the same thing as operator().
     */
    Number &operator [] (const size_type global_index);

    /**
     * Instead of getting individual elements of a vector, this function
     * allows to get a whole set of elements at once. The indices of the
     * elements to be read are stated in the first argument, the corresponding
     * values are returned in the second.
     */
    template <typename Number2>
    void extract_subvector_to (const std::vector<size_type> &indices,
                               std::vector<Number2> &values) const;

    /**
     * Just as the above, but with pointers. Useful in minimizing copying of
     * data around.
     */
    template <typename ForwardIterator, typename OutputIterator>
    void extract_subvector_to (ForwardIterator          indices_begin,
                               const ForwardIterator    indices_end,
                               OutputIterator           values_begin) const;

    /**
     * Read access to the data field specified by @p local_index. When you
     * access elements in the order in which they are stored, it is necessary
     * that you know in which they are stored. In other words, you need to
     * know the map between the global indices of the elements this class
     * stores, and the local indices into the contiguous array of these global
     * elements. For this, see the general documentation of this class.
     *
     * Performance: Direct array access (fast).
     */
    Number local_element (const size_type local_index) const;

    /**
     * Read and write access to the data field specified by @p local_index.
     * When you access elements in the order in which they are stored, it is
     * necessary that you know in which they are stored. In other words, you
     * need to know the map between the global indices of the elements this
     * class stores, and the local indices into the contiguous array of these
     * global elements. For this, see the general documentation of this class.
     *
     * Performance: Direct array access (fast).
     */
    Number &local_element (const size_type local_index);
    //@}


    /**
     * @name 3: Modification of vectors
     */
    //@{

    /**
     * This function adds a whole set of values stored in @p values to the
     * vector components specified by @p indices.
     */
    template <typename Number2>
    void add (const std::vector<size_type>  &indices,
              const std::vector<Number2>    &values);

    /**
     * This function is similar to the previous one but takes a
     * ReadWriteVector of values.
     */
    template <typename Number2>
    void add (const std::vector<size_type>   &indices,
              const ReadWriteVector<Number2> &values);

    /**
     * Take an address where <tt>n_elements</tt> are stored contiguously and
     * add them into the vector. Handles all cases which are not covered by
     * the other two <tt>add()</tt> functions above.
     */
    template <typename Number2>
    void add (const size_type  n_elements,
              const size_type *indices,
              const Number2   *values);

    /**
     * Prints the vector to the output stream @p out.
     */
    void print (std::ostream       &out,
                const unsigned int  precision  = 3,
                const bool          scientific = true) const;

    /**
     * Returns the memory consumption of this class in bytes.
     */
    std::size_t memory_consumption () const;
    //@}

  protected:
#ifdef DEAL_II_WITH_TRILINOS
    /**
     * Import all the elements present in the vector's IndexSet from the input
     * vector @p multivector. This is an helper function and it should not be
     * used directly.
     */
    void import(const Epetra_MultiVector                       &multivector,
                const IndexSet                                 &locally_owned_elements,
                VectorOperation::values                         operation,
                const MPI_Comm                                 &mpi_comm,
                std_cxx11::shared_ptr<const CommunicationPatternBase> communication_pattern);
#endif

    /**
     * Return the local position of @p global_index.
     */
    unsigned int
    global_to_local (const types::global_dof_index global_index) const
    {
      // the following will throw an exception if the global_index is not
      // in the remaining_elements
      return static_cast<unsigned int>(stored_elements.index_within_set(global_index));
    }

    /**
     * A helper function that is used to resize the val array.
     */
    void resize_val (const size_type new_allocated_size);

#if defined(DEAL_II_WITH_TRILINOS) && defined(DEAL_II_WITH_MPI)
    /**
     * Return a EpetraWrappers::Communication pattern and store it for future
     * use.
     */
    EpetraWrappers::CommunicationPattern
    create_epetra_comm_pattern(const IndexSet &source_index_set,
                               const MPI_Comm &mpi_comm);
#endif

    /**
     * Indices of the elements stored.
     */
    IndexSet stored_elements;

    /**
     * IndexSet of the elements of the last imported vector;
     */
    IndexSet source_stored_elements;

    /**
     * CommunicationPattern for the communication between the
     * source_stored_elements IndexSet and the current vector.
     */
    std_cxx11::shared_ptr<CommunicationPatternBase> comm_pattern;

    /**
     * Pointer to the array of local elements of this vector.
     */
    Number *val;

    /**
     * For parallel loops with TBB, this member variable stores the affinity
     * information of loops.
     */
    mutable std_cxx11::shared_ptr<parallel::internal::TBBPartitioner> thread_loop_partitioner;

    /**
     * Make all other ReadWriteVector types friends.
     */
    template <typename Number2> friend class ReadWriteVector;
  };

  /*@}*/


  /*----------------------- Inline functions ----------------------------------*/

#ifndef DOXYGEN

  template <typename Number>
  inline
  ReadWriteVector<Number>::ReadWriteVector ()
    :
    val(NULL)
  {
    reinit(0, true);
  }



  template <typename Number>
  inline
  ReadWriteVector<Number>::ReadWriteVector (const ReadWriteVector<Number> &v)
    :
    Subscriptor(),
    val(NULL)
  {
    this->operator=(v);
  }



  template <typename Number>
  inline
  ReadWriteVector<Number>::ReadWriteVector (const size_type size)
    :
    val(NULL)
  {
    reinit (size, false);
  }



  template <typename Number>
  inline
  ReadWriteVector<Number>::ReadWriteVector (const IndexSet &locally_stored_indices)
    :
    val(NULL)
  {
    reinit (locally_stored_indices);
  }



  template <typename Number>
  inline
  ReadWriteVector<Number>::~ReadWriteVector ()
  {
    resize_val(0);
  }



  template <typename Number>
  inline
  typename ReadWriteVector<Number>::size_type
  ReadWriteVector<Number>::size() const
  {
    return stored_elements.size();
  }



  template <typename Number>
  inline
  typename ReadWriteVector<Number>::size_type
  ReadWriteVector<Number>::n_elements() const
  {
    return stored_elements.n_elements();
  }



  template <typename Number>
  inline
  const IndexSet &
  ReadWriteVector<Number>::get_stored_elements () const
  {
    return stored_elements;
  }



  template <typename Number>
  inline
  typename ReadWriteVector<Number>::iterator
  ReadWriteVector<Number>::begin ()
  {
    return &val[0];
  }



  template <typename Number>
  inline
  typename ReadWriteVector<Number>::const_iterator
  ReadWriteVector<Number>::begin () const
  {
    return &val[0];
  }



  template <typename Number>
  inline
  typename ReadWriteVector<Number>::iterator
  ReadWriteVector<Number>::end ()
  {
    return &val[this->n_elements()];
  }



  template <typename Number>
  inline
  typename ReadWriteVector<Number>::const_iterator
  ReadWriteVector<Number>::end () const
  {
    return &val[this->n_elements()];
  }



  template <typename Number>
  inline
  Number
  ReadWriteVector<Number>::operator() (const size_type global_index) const
  {
    return val[global_to_local(global_index)];
  }



  template <typename Number>
  inline
  Number &
  ReadWriteVector<Number>::operator() (const size_type global_index)
  {
    return val[global_to_local (global_index)];
  }



  template <typename Number>
  inline
  Number
  ReadWriteVector<Number>::operator[] (const size_type global_index) const
  {
    return operator()(global_index);
  }



  template <typename Number>
  inline
  Number &
  ReadWriteVector<Number>::operator[] (const size_type global_index)
  {
    return operator()(global_index);
  }



  template <typename Number>
  template <typename Number2>
  inline
  void ReadWriteVector<Number>::extract_subvector_to (const std::vector<size_type> &indices,
                                                      std::vector<Number2> &values) const
  {
    for (size_type i = 0; i < indices.size(); ++i)
      values[i] = operator()(indices[i]);
  }



  template <typename Number>
  template <typename ForwardIterator, typename OutputIterator>
  inline
  void ReadWriteVector<Number>::extract_subvector_to (ForwardIterator          indices_begin,
                                                      const ForwardIterator    indices_end,
                                                      OutputIterator           values_begin) const
  {
    while (indices_begin != indices_end)
      {
        *values_begin = operator()(*indices_begin);
        indices_begin++;
        values_begin++;
      }
  }



  template <typename Number>
  inline
  Number
  ReadWriteVector<Number>::local_element (const size_type local_index) const
  {
    AssertIndexRange (local_index, this->n_elements());

    return val[local_index];
  }



  template <typename Number>
  inline
  Number &
  ReadWriteVector<Number>::local_element (const size_type local_index)
  {
    AssertIndexRange (local_index, this->n_elements());

    return val[local_index];
  }



  template <typename Number>
  template <typename Number2>
  inline
  void
  ReadWriteVector<Number>::add (const std::vector<size_type> &indices,
                                const std::vector<Number2>   &values)
  {
    AssertDimension (indices.size(), values.size());
    add (indices.size(), &indices[0], &values[0]);
  }



  template <typename Number>
  template <typename Number2>
  inline
  void
  ReadWriteVector<Number>::add (const std::vector<size_type>   &indices,
                                const ReadWriteVector<Number2> &values)
  {
    const size_type size = indices.size();
    for (size_type i=0; i<size; ++i)
      {
        Assert (numbers::is_finite(values[i]),
                ExcMessage("The given value is not finite but either infinite or Not A Number (NaN)"));
        this->operator()(indices[i]) += values[indices[i]];
      }
  }



  template <typename Number>
  template <typename Number2>
  inline
  void
  ReadWriteVector<Number>::add (const size_type    n_indices,
                                const size_type   *indices,
                                const Number2     *values)
  {
    for (size_type i=0; i<n_indices; ++i)
      {
        Assert (numbers::is_finite(values[i]),
                ExcMessage("The given value is not finite but either infinite or Not A Number (NaN)"));
        this->operator()(indices[i]) += values[i];
      }
  }

#endif  // ifndef DOXYGEN

} // end of namespace LinearAlgebra




/**
 * Global function @p swap which overloads the default implementation of the
 * C++ standard library which uses a temporary object. The function simply
 * exchanges the data of the two vectors.
 *
 * @relates Vector
 */
template <typename Number>
inline
void swap (LinearAlgebra::ReadWriteVector<Number> &u,
           LinearAlgebra::ReadWriteVector<Number> &v)
{
  u.swap (v);
}


DEAL_II_NAMESPACE_CLOSE

#endif
