// ---------------------------------------------------------------------
//
// Copyright (C) 1999 - 2013 by the deal.II authors
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

#ifndef __deal2__parallel_block_vector_h
#define __deal2__parallel_block_vector_h


#include <deal.II/base/config.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/lac/block_indices.h>
#include <deal.II/lac/block_vector_base.h>
#include <deal.II/lac/parallel_vector.h>

#include <deal.II/lac/petsc_parallel_block_vector.h>
#include <deal.II/lac/trilinos_parallel_block_vector.h>


#include <cstdio>
#include <vector>

DEAL_II_NAMESPACE_OPEN


namespace parallel
{
  namespace distributed
  {

    /*! @addtogroup Vectors
     *@{
     */


    /**
     * An implementation of block vectors based on distribued deal.II
     * vectors. While the base class provides for most of the interface, this
     * class handles the actual allocation of vectors and provides functions that
     * are specific to the underlying vector type.
     *
     * @note Instantiations for this template are provided for <tt>@<float@> and
     * @<double@></tt>; others can be generated in application programs (see the
     * section on @ref Instantiations in the manual).
     *
     * @see @ref GlossBlockLA "Block (linear algebra)"
     * @author Katharina Kormann, Martin Kronbichler, 2011
     */
    template <typename Number>
    class BlockVector : public BlockVectorBase<Vector<Number> >
    {
    public:
      /**
       * Typedef the base class for simpler access to its own typedefs.
       */
      typedef BlockVectorBase<Vector<Number> > BaseClass;

      /**
       * Typedef the type of the underlying vector.
       */
      typedef typename BaseClass::BlockType  BlockType;

      /**
       * Import the typedefs from the base class.
       */
      typedef typename BaseClass::value_type      value_type;
      typedef typename BaseClass::real_type       real_type;
      typedef typename BaseClass::pointer         pointer;
      typedef typename BaseClass::const_pointer   const_pointer;
      typedef typename BaseClass::reference       reference;
      typedef typename BaseClass::const_reference const_reference;
      typedef typename BaseClass::size_type       size_type;
      typedef typename BaseClass::iterator        iterator;
      typedef typename BaseClass::const_iterator  const_iterator;

      /**
       *  Constructor. There are three ways to use this constructor. First,
       *  without any arguments, it generates an object with no blocks. Given
       *  one argument, it initializes <tt>num_blocks</tt> blocks, but these
       *  blocks have size zero. The third variant finally initializes all
       *  blocks to the same size <tt>block_size</tt>.
       *
       *  Confer the other constructor further down if you intend to use
       *  blocks of different sizes.
       */
      explicit BlockVector (const size_type num_blocks = 0,
                            const size_type block_size = 0);

      /**
       * Copy-Constructor. Dimension set to that of V, all components are
       * copied from V
       */
      BlockVector (const BlockVector<Number> &V);


#ifndef DEAL_II_EXPLICIT_CONSTRUCTOR_BUG
      /**
       * Copy constructor taking a BlockVector of another data type. This will
       * fail if there is no conversion path from <tt>OtherNumber</tt> to
       * <tt>Number</tt>. Note that you may lose accuracy when copying to a
       * BlockVector with data elements with less accuracy.
       *
       * Older versions of gcc did not honor the @p explicit keyword on
       * template constructors. In such cases, it is easy to accidentally
       * write code that can be very inefficient, since the compiler starts
       * performing hidden conversions. To avoid this, this function is
       * disabled if we have detected a broken compiler during configuration.
       */
      template <typename OtherNumber>
      explicit
      BlockVector (const BlockVector<OtherNumber> &v);
#endif

      /**
       * Constructor. Set the number of blocks to <tt>block_sizes.size()</tt>
       * and initialize each block with <tt>block_sizes[i]</tt> zero elements.
       */
      BlockVector (const std::vector<size_type> &block_sizes);

      /**
       * Construct a block vector with an IndexSet for the local range
       * and ghost entries for each block.
       */
      BlockVector (const std::vector<IndexSet> &local_ranges,
                   const std::vector<IndexSet> &ghost_indices,
                   const MPI_Comm  communicator);

      /**
       * Same as above but the ghost indicies are assumed to be empty.
       */
      BlockVector (const std::vector<IndexSet> &local_ranges,
                   const MPI_Comm  communicator);

      /**
       * Destructor. Clears memory.
       */
      ~BlockVector ();

      /**
       * Copy operator: fill all components of the vector with the given
       * scalar value.
       */
      BlockVector &operator = (const value_type s);

      /**
       * Copy operator for arguments of the same type. Resize the present
       * vector if necessary.
       */
      BlockVector &
      operator= (const BlockVector &V);

      /**
       * Copy operator for template arguments of different types. Resize the
       * present vector if necessary.
       */
      template <class Number2>
      BlockVector &
      operator= (const BlockVector<Number2> &V);

      /**
       * Copy a regular vector into a block vector.
       */
      BlockVector &
      operator= (const Vector<Number> &V);

#ifdef DEAL_II_WITH_PETSC
      /**
       * Copy the content of a PETSc vector into the calling vector. This
       * function assumes that the vectors layouts have already been
       * initialized to match.
       *
       * This operator is only available if deal.II was configured with PETSc.
       */
      BlockVector<Number> &
      operator = (const PETScWrappers::MPI::BlockVector &petsc_vec);
#endif

#ifdef DEAL_II_WITH_TRILINOS
      /**
       * Copy the content of a Trilinos vector into the calling vector. This
       * function assumes that the vectors layouts have already been
       * initialized to match.
       *
       * This operator is only available if deal.II was configured with
       * Trilinos.
       */
      BlockVector<Number> &
      operator = (const TrilinosWrappers::MPI::BlockVector &trilinos_vec);
#endif

      /**
       * Reinitialize the BlockVector to contain <tt>num_blocks</tt> blocks of
       * size <tt>block_size</tt> each.
       *
       * If the second argument is left at its default value, then the block
       * vector allocates the specified number of blocks but leaves them at
       * zero size. You then need to later reinitialize the individual blocks,
       * and call collect_sizes() to update the block system's knowledge of
       * its individual block's sizes.
       *
       * If <tt>fast==false</tt>, the vector is filled with zeros.
       */
      void reinit (const size_type num_blocks,
                   const size_type block_size = 0,
                   const bool fast = false);

      /**
       * Reinitialize the BlockVector such that it contains
       * <tt>block_sizes.size()</tt> blocks. Each block is reinitialized to
       * dimension <tt>block_sizes[i]</tt>.
       *
       * If the number of blocks is the same as before this function was
       * called, all vectors remain the same and reinit() is called for each
       * vector.
       *
       * If <tt>fast==false</tt>, the vector is filled with zeros.
       *
       * Note that you must call this (or the other reinit() functions)
       * function, rather than calling the reinit() functions of an individual
       * block, to allow the block vector to update its caches of vector
       * sizes. If you call reinit() on one of the blocks, then subsequent
       * actions on this object may yield unpredictable results since they may
       * be routed to the wrong block.
       */
      void reinit (const std::vector<size_type> &N,
                   const bool                    fast=false);

      /**
       * Change the dimension to that of the vector <tt>V</tt>. The same
       * applies as for the other reinit() function.
       *
       * The elements of <tt>V</tt> are not copied, i.e.  this function is the
       * same as calling <tt>reinit (V.size(), fast)</tt>.
       *
       * Note that you must call this (or the other reinit() functions)
       * function, rather than calling the reinit() functions of an individual
       * block, to allow the block vector to update its caches of vector
       * sizes. If you call reinit() of one of the blocks, then subsequent
       * actions of this object may yield unpredictable results since they may
       * be routed to the wrong block.
       */
      template <typename Number2>
      void reinit (const BlockVector<Number2> &V,
                   const bool                 fast=false);

      /**
       * This function copies the data that has accumulated in the data buffer
       * for ghost indices to the owning processor. For the meaning of the
       * argument @p operation, see the entry on @ref GlossCompress
       * "Compressing distributed vectors and matrices" in the glossary.
       *
       * There are two variants for this function. If called with argument @p
       * VectorOperation::add adds all the data accumulated in ghost elements
       * to the respective elements on the owning processor and clears the
       * ghost array afterwards. If called with argument @p
       * VectorOperation::insert, a set operation is performed. Since setting
       * elements in a vector with ghost elements is ambiguous (as one can set
       * both the element on the ghost site as well as the owning site), this
       * operation makes the assumption that all data is set correctly on the
       * owning processor. Upon call of compress(VectorOperation::insert), all
       * ghost entries are therefore simply zeroed out (using
       * zero_ghost_values()). In debug mode, a check is performed that makes
       * sure that the data set is actually consistent between processors,
       * i.e., whenever a non-zero ghost element is found, it is compared to
       * the value on the owning processor and an exception is thrown if these
       * elements do not agree.
       *
       */
      void compress (::dealii::VectorOperation::values operation);

      /**
       * Fills the data field for ghost indices with the values stored in the
       * respective positions of the owning processor. This function is needed
       * before reading from ghosts. The function is @p const even though
       * ghost data is changed. This is needed to allow functions with a @p
       * const vector to perform the data exchange without creating
       * temporaries.
       */
      void update_ghost_values () const;

      /**
       * This method zeros the entries on ghost dofs, but does not touch
       * locally owned DoFs.
       *
       * After calling this method, read access to ghost elements of the
       * vector is forbidden and an exception is thrown. Only write access to
       * ghost elements is allowed in this state.
       */
      void zero_out_ghosts ();

      /**
       * Returns if this Vector contains ghost elements.
       */
      bool has_ghost_elements() const;

      /**
       * Return whether the vector contains only elements with value
       * zero. This function is mainly for internal consistency checks and
       * should seldom be used when not in debug mode since it uses quite some
       * time.
       */
      bool all_zero () const;

      /**
       * Return @p true if the vector has no negative entries, i.e. all
       * entries are zero or positive. This function is used, for example, to
       * check whether refinement indicators are really all positive (or
       * zero).
       *
       * The function obviously only makes sense if the template argument of
       * this class is a real type. If it is a complex type, then an exception
       * is thrown.
       */
      bool is_non_negative () const;

      /**
       * Checks for equality of the two vectors.
       */
      template <typename Number2>
      bool operator == (const BlockVector<Number2> &v) const;

      /**
       * Checks for inequality of the two vectors.
       */
      template <typename Number2>
      bool operator != (const BlockVector<Number2> &v) const;

      /**
       * Perform the inner product of two vectors.
       */
      template <typename Number2>
      Number operator * (const BlockVector<Number2> &V) const;

      /**
       * Computes the square of the l<sub>2</sub> norm of the vector (i.e.,
       * the sum of the squares of all entries among all processors).
       */
      real_type norm_sqr () const;

      /**
       * Computes the mean value of all the entries in the vector.
       */
      Number mean_value () const;

      /**
       * Returns the l<sub>1</sub> norm of the vector (i.e., the sum of the
       * absolute values of all entries among all processors).
       */
      real_type l1_norm () const;

      /**
       * Returns the l<sub>2</sub> norm of the vector (i.e., square root of
       * the sum of the square of all entries among all processors).
       */
      real_type l2_norm () const;

      /**
       * Returns the l<sub>p</sub> norm with real @p p of the vector (i.e.,
       * the pth root of sum of the pth power of all entries among all
       * processors).
       */
      real_type lp_norm (const real_type p) const;

      /**
       * Returns the maximum norm of the vector (i.e., maximum absolute value
       * among all entries among all processors).
       */
      real_type linfty_norm () const;

      /**
       * Scale each element of the vector by the given factor.
       *
       * This function is deprecated and will be removed in a future
       * version. Use <tt>operator *=</tt> and <tt>operator /=</tt> instead.
       *
       * @deprecated Use <tt>operator*=</tt> instead.
       */
      void scale (const value_type factor) DEAL_II_DEPRECATED;

      /**
       * Multiply each element of this vector by the corresponding element of
       * <tt>v</tt>.
       */
      template <class BlockVector2>
      void scale (const BlockVector2 &v);

      /**
       * Swap the contents of this vector and the other vector <tt>v</tt>. One
       * could do this operation with a temporary variable and copying over
       * the data elements, but this function is significantly more efficient
       * since it only swaps the pointers to the data of the two vectors and
       * therefore does not need to allocate temporary storage and move data
       * around.
       *
       * Limitation: right now this function only works if both vectors have
       * the same number of blocks. If needed, the numbers of blocks should be
       * exchanged, too.
       *
       * This function is analog to the the swap() function of all C++
       * standard containers. Also, there is a global function swap(u,v) that
       * simply calls <tt>u.swap(v)</tt>, again in analogy to standard
       * functions.
       */
      void swap (BlockVector<Number> &v);

      /** @addtogroup Exceptions
       * @{ */

      /**
       * Exception
       */
      DeclException0 (ExcIteratorRangeDoesNotMatchVectorSize);
      //@}
    };

    /*@}*/

#ifndef DOXYGEN
    /*----------------------- Inline functions ----------------------------------*/


    template <typename Number>
    inline
    BlockVector<Number>::BlockVector (const size_type n_blocks,
                                      const size_type block_size)
    {
      reinit (n_blocks, block_size);
    }



    template <typename Number>
    inline
    BlockVector<Number>::BlockVector (const std::vector<size_type> &n)
    {
      reinit (n, false);
    }


    template <typename Number>
    inline
    BlockVector<Number>::BlockVector (const std::vector<IndexSet> &local_ranges,
                                      const std::vector<IndexSet> &ghost_indices,
                                      const MPI_Comm  communicator)
    {
      std::vector<size_type> sizes(local_ranges.size());
      for (unsigned int i=0; i<local_ranges.size(); ++i)
        sizes[i] = local_ranges[i].size();

      this->block_indices.reinit(sizes);
      this->components.resize(this->n_blocks());

      for (unsigned int i=0; i<this->n_blocks(); ++i)
        this->block(i).reinit(local_ranges[i], ghost_indices[i], communicator);
    }


    template <typename Number>
    inline
    BlockVector<Number>::BlockVector (const std::vector<IndexSet> &local_ranges,
                                      const MPI_Comm  communicator)
    {
      std::vector<size_type> sizes(local_ranges.size());
      for (unsigned int i=0; i<local_ranges.size(); ++i)
        sizes[i] = local_ranges[i].size();

      this->block_indices.reinit(sizes);
      this->components.resize(this->n_blocks());

      for (unsigned int i=0; i<this->n_blocks(); ++i)
        this->block(i).reinit(local_ranges[i], communicator);
    }



    template <typename Number>
    inline
    BlockVector<Number>::BlockVector (const BlockVector<Number> &v)
      :
      BlockVectorBase<Vector<Number> > ()
    {
      this->components.resize (v.n_blocks());
      this->block_indices = v.block_indices;

      for (size_type i=0; i<this->n_blocks(); ++i)
        this->components[i] = v.components[i];
    }


#ifndef DEAL_II_EXPLICIT_CONSTRUCTOR_BUG

    template <typename Number>
    template <typename OtherNumber>
    inline
    BlockVector<Number>::BlockVector (const BlockVector<OtherNumber> &v)
    {
      reinit (v, true);
      *this = v;
    }

#endif



    template <typename Number>
    inline
    void BlockVector<Number>::reinit (const size_type n_bl,
                                      const size_type bl_sz,
                                      const bool         fast)
    {
      std::vector<size_type> n(n_bl, bl_sz);
      reinit(n, fast);
    }


    template <typename Number>
    inline
    void BlockVector<Number>::reinit (const std::vector<size_type> &n,
                                      const bool                    fast)
    {
      this->block_indices.reinit (n);
      if (this->components.size() != this->n_blocks())
        this->components.resize(this->n_blocks());

      for (size_type i=0; i<this->n_blocks(); ++i)
        this->components[i].reinit(n[i], fast);
    }



    template <typename Number>
    template <typename Number2>
    inline
    void BlockVector<Number>::reinit (const BlockVector<Number2> &v,
                                      const bool fast)
    {
      this->block_indices = v.get_block_indices();
      if (this->components.size() != this->n_blocks())
        this->components.resize(this->n_blocks());

      for (unsigned int i=0; i<this->n_blocks(); ++i)
        this->block(i).reinit(v.block(i), fast);
    }



    template <typename Number>
    inline
    BlockVector<Number>::~BlockVector ()
    {}



    template <typename Number>
    inline
    BlockVector<Number> &
    BlockVector<Number>::operator = (const value_type s)
    {

      Assert (numbers::is_finite(s), ExcNumberNotFinite());

      BaseClass::operator = (s);
      return *this;
    }



    template <typename Number>
    inline
    BlockVector<Number> &
    BlockVector<Number>::operator = (const BlockVector &v)
    {
      // we only allow assignment to vectors with the same number of blocks
      // or to an empty BlockVector
      Assert (this->n_blocks() == 0 || this->n_blocks() == v.n_blocks(),
              ExcDimensionMismatch(this->n_blocks(), v.n_blocks()));

      if (this->n_blocks() != v.n_blocks())
        reinit(v.n_blocks(), true);

      for (size_type i=0; i<this->n_blocks(); ++i)
        this->components[i] = v.block(i);

      this->collect_sizes();
      return *this;
    }



    template <typename Number>
    inline
    BlockVector<Number> &
    BlockVector<Number>::operator = (const Vector<Number> &v)
    {
      BaseClass::operator = (v);
      return *this;
    }



    template <typename Number>
    template <typename Number2>
    inline
    BlockVector<Number> &
    BlockVector<Number>::operator = (const BlockVector<Number2> &v)
    {
      reinit (v, true);
      BaseClass::operator = (v);
      return *this;
    }



#ifdef DEAL_II_WITH_PETSC

    template <typename Number>
    inline
    BlockVector<Number> &
    BlockVector<Number>::operator = (const PETScWrappers::MPI::BlockVector &petsc_vec)
    {
      AssertDimension(this->n_blocks(), petsc_vec.n_blocks());
      for (unsigned int i=0; i<this->n_blocks(); ++i)
        this->block(i) = petsc_vec.block(i);

      return *this;
    }

#endif



#ifdef DEAL_II_WITH_TRILINOS

    template <typename Number>
    inline
    BlockVector<Number> &
    BlockVector<Number>::operator = (const TrilinosWrappers::MPI::BlockVector &trilinos_vec)
    {
      AssertDimension(this->n_blocks(), trilinos_vec.n_blocks());
      for (unsigned int i=0; i<this->n_blocks(); ++i)
        this->block(i) = trilinos_vec.block(i);

      return *this;
    }

#endif



    template <typename Number>
    inline
    void
    BlockVector<Number>::compress (::dealii::VectorOperation::values operation)
    {
      // start all requests for all blocks before finishing the transfers as
      // this saves repeated synchronizations
      for (unsigned int block=0; block<this->n_blocks(); ++block)
        this->block(block).compress_start(block*10 + 8273, operation);
      for (unsigned int block=0; block<this->n_blocks(); ++block)
        this->block(block).compress_finish(operation);
    }



    template <typename Number>
    inline
    void
    BlockVector<Number>::update_ghost_values () const
    {
      for (unsigned int block=0; block<this->n_blocks(); ++block)
        this->block(block).update_ghost_values_start(block*10 + 9923);
      for (unsigned int block=0; block<this->n_blocks(); ++block)
        this->block(block).update_ghost_values_finish();
    }



    template <typename Number>
    inline
    void
    BlockVector<Number>::zero_out_ghosts ()
    {
      for (unsigned int block=0; block<this->n_blocks(); ++block)
        this->block(block).zero_out_ghosts();
    }



    template <typename Number>
    inline
    bool
    BlockVector<Number>::has_ghost_elements () const
    {
      bool has_ghost_elements = false;
      for (unsigned int block=0; block<this->n_blocks(); ++block)
        if (this->block(block).has_ghost_elements() == true)
          has_ghost_elements = true;
      return has_ghost_elements;
    }



    template <typename Number>
    inline
    bool
    BlockVector<Number>::all_zero () const
    {
      Assert (this->n_blocks() > 0, ExcEmptyObject());

      // use int instead of bool. in order to make global reduction operations
      // work also when MPI_Init was not called, only call MPI_Allreduce
      // commands when there is more than one processor (note that reinit()
      // functions handle this case correctly through the job_supports_mpi()
      // query). this is the same in all the functions below
      int local_result = -1;
      for (unsigned int i=0; i<this->n_blocks(); ++i)
        local_result = std::max(local_result,
                                -static_cast<int>(this->block(i).all_zero_local()));

      if (this->block(0).partitioner->n_mpi_processes() > 1)
        return -Utilities::MPI::max(local_result,
                                    this->block(0).partitioner->get_communicator());
      else
        return local_result;
    }



    template <typename Number>
    inline
    bool
    BlockVector<Number>::is_non_negative () const
    {
      Assert (this->n_blocks() > 0, ExcEmptyObject());
      int local_result = -1;
      for (unsigned int i=0; i<this->n_blocks(); ++i)
        local_result = std::max(local_result,
                                -static_cast<int>(this->block(i).is_non_negative_local()));
      if (this->block(0).partitioner->n_mpi_processes() > 1)
        return Utilities::MPI::max(local_result,
                                   this->block(0).partitioner->get_communicator());
      else
        return local_result;
    }



    template <typename Number>
    template <typename Number2>
    inline
    bool
    BlockVector<Number>::operator == (const BlockVector<Number2> &v) const
    {
      Assert (this->n_blocks() > 0, ExcEmptyObject());
      AssertDimension (this->n_blocks(), v.n_blocks());

      // MPI does not support bools, so use unsigned int instead. Two vectors
      // are equal if the check for non-equal fails on all processors
      unsigned int local_result = 0;
      for (unsigned int i=0; i<this->n_blocks(); ++i)
        local_result = std::max(local_result,
                                static_cast<unsigned int>(!this->block(i).vectors_equal_local(v.block(i))));
      unsigned int result =
        this->block(0).partitioner->n_mpi_processes() > 1
        ?
        Utilities::MPI::max(local_result, this->block(0).partitioner->get_communicator())
        :
        local_result;
      return result==0;
    }



    template <typename Number>
    template <typename Number2>
    inline
    bool
    BlockVector<Number>::operator != (const BlockVector<Number2> &v) const
    {
      return !(operator == (v));
    }



    template <typename Number>
    template <typename Number2>
    inline
    Number
    BlockVector<Number>::operator * (const BlockVector<Number2> &v) const
    {
      Assert (this->n_blocks() > 0, ExcEmptyObject());
      AssertDimension (this->n_blocks(), v.n_blocks());

      Number local_result = Number();
      for (unsigned int i=0; i<this->n_blocks(); ++i)
        local_result += this->block(i).inner_product_local(v.block(i));

      if (this->block(0).partitioner->n_mpi_processes() > 1)
        return Utilities::MPI::sum (local_result,
                                    this->block(0).partitioner->get_communicator());
      else
        return local_result;
    }



    template <typename Number>
    inline
    typename BlockVector<Number>::real_type
    BlockVector<Number>::norm_sqr () const
    {
      Assert (this->n_blocks() > 0, ExcEmptyObject());

      real_type local_result = real_type();
      for (unsigned int i=0; i<this->n_blocks(); ++i)
        local_result += this->block(i).norm_sqr_local();

      if (this->block(0).partitioner->n_mpi_processes() > 1)
        return Utilities::MPI::sum (local_result,
                                    this->block(0).partitioner->get_communicator());
      else
        return local_result;
    }



    template <typename Number>
    inline
    Number
    BlockVector<Number>::mean_value () const
    {
      Assert (this->n_blocks() > 0, ExcEmptyObject());

      Number local_result = Number();
      for (unsigned int i=0; i<this->n_blocks(); ++i)
        local_result += this->block(i).mean_value_local()*(real_type)this->block(i).partitioner->local_size();

      if (this->block(0).partitioner->n_mpi_processes() > 1)
        return Utilities::MPI::sum (local_result,
                                    this->block(0).partitioner->get_communicator())/
               (real_type)this->size();
      else
        return local_result/(real_type)this->size();
    }



    template <typename Number>
    inline
    typename BlockVector<Number>::real_type
    BlockVector<Number>::l1_norm () const
    {
      Assert (this->n_blocks() > 0, ExcEmptyObject());

      real_type local_result = real_type();
      for (unsigned int i=0; i<this->n_blocks(); ++i)
        local_result += this->block(i).l1_norm_local();

      if (this->block(0).partitioner->n_mpi_processes() > 1)
        return Utilities::MPI::sum (local_result,
                                    this->block(0).partitioner->get_communicator());
      else
        return local_result;
    }



    template <typename Number>
    inline
    typename BlockVector<Number>::real_type
    BlockVector<Number>::l2_norm () const
    {
      return std::sqrt(norm_sqr());
    }



    template <typename Number>
    inline
    typename BlockVector<Number>::real_type
    BlockVector<Number>::lp_norm (const real_type p) const
    {
      Assert (this->n_blocks() > 0, ExcEmptyObject());

      real_type local_result = real_type();
      for (unsigned int i=0; i<this->n_blocks(); ++i)
        local_result += std::pow(this->block(i).lp_norm_local(p), p);

      if (this->block(0).partitioner->n_mpi_processes() > 1)
        return std::pow (Utilities::MPI::sum(local_result,
                                             this->block(0).partitioner->get_communicator()),
                         static_cast<real_type>(1.0/p));
      else
        return std::pow (local_result, static_cast<real_type>(1.0/p));
    }



    template <typename Number>
    inline
    typename BlockVector<Number>::real_type
    BlockVector<Number>::linfty_norm () const
    {
      Assert (this->n_blocks() > 0, ExcEmptyObject());

      real_type local_result = real_type();
      for (unsigned int i=0; i<this->n_blocks(); ++i)
        local_result = std::max(local_result, this->block(i).linfty_norm_local());

      if (this->block(0).partitioner->n_mpi_processes() > 1)
        return Utilities::MPI::max (local_result,
                                    this->block(0).partitioner->get_communicator());
      else
        return local_result;
    }



    template <typename Number>
    inline
    void BlockVector<Number>::swap (BlockVector<Number> &v)
    {
      Assert (this->n_blocks() == v.n_blocks(),
              ExcDimensionMismatch(this->n_blocks(), v.n_blocks()));

      for (size_type i=0; i<this->n_blocks(); ++i)
        dealii::swap (this->components[i], v.components[i]);
      dealii::swap (this->block_indices, v.block_indices);
    }



    template <typename Number>
    void BlockVector<Number>::scale (const value_type factor)
    {

      Assert (numbers::is_finite(factor), ExcNumberNotFinite());

      for (size_type i=0; i<this->n_blocks(); ++i)
        this->components[i].scale(factor);
    }



    template <typename Number>
    template <class BlockVector2>
    void BlockVector<Number>::scale (const BlockVector2 &v)
    {
      BaseClass::scale (v);
    }

#endif // DOXYGEN

  } // end of namespace distributed

} // end of namespace parallel

/**
 * Global function which overloads the default implementation
 * of the C++ standard library which uses a temporary object. The
 * function simply exchanges the data of the two vectors.
 *
 * @relates BlockVector
 * @author Katharina Kormann, Martin Kronbichler, 2011
 */
template <typename Number>
inline
void swap (parallel::distributed::BlockVector<Number> &u,
           parallel::distributed::BlockVector<Number> &v)
{
  u.swap (v);
}

DEAL_II_NAMESPACE_CLOSE

#endif
