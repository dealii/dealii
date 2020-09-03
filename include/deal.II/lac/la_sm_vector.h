// ---------------------------------------------------------------------
//
// Copyright (C) 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------

#ifndef dealii_la_sm_vector_h
#define dealii_la_sm_vector_h

#include <deal.II/base/config.h>

#include <deal.II/base/memory_space.h>
#include <deal.II/base/mpi.h>
#include <deal.II/base/numbers.h>
#include <deal.II/base/partitioner.h>
#include <deal.II/base/thread_management.h>

#include <deal.II/lac/la_sm_partitioner.h>
#include <deal.II/lac/vector_operation.h>
#include <deal.II/lac/vector_space_vector.h>
#include <deal.II/lac/vector_type_traits.h>

#include <iomanip>
#include <memory>

DEAL_II_NAMESPACE_OPEN

namespace LinearAlgebra
{
  namespace SharedMPI
  {
    /**
     * Data for MPI shared memory.
     */
    template <typename Number>
    struct MemorySpaceData
      : public ::dealii::MemorySpace::MemorySpaceData<Number, MemorySpace::Host>
    {
      /**
       * Copy the active data (values) to @p begin.
       *
       * @note Not implemented yet.
       */
      void
      copy_to(Number *begin, std::size_t n_elements) override
      {
        Assert(false, ExcNotImplemented());
        (void)begin;
        (void)n_elements;
      }

      /**
       * Copy the data in @p begin to the active data of the structure (values).
       *
       * @note Not implemented yet.
       */
      void
      copy_from(Number *begin, std::size_t n_elements) override
      {
        Assert(false, ExcNotImplemented());
        (void)begin;
        (void)n_elements;
      }

      /**
       * Return memory consumption.
       *
       * @note It also includes added memory for alignment.
       */
      std::size_t
      memory_consumption() const
      {
        // note: values_win is not accounted for
        return MemoryConsumption::memory_consumption(others) +
               MemoryConsumption::memory_consumption(
                 memory_constumption_values) +
               memory_constumption_values;
      }

      /**
       * Memory consumption of the allocated values (reserved during
       * MPI_Win_allocate_shared).
       *
       * @note It also includes added memory for alignment.
       */
      std::size_t memory_constumption_values;

#ifdef DEAL_II_WITH_MPI
      /**
       * MPI window. It is connected to the destructor of the values; its
       * destruction (MPI_Win_free) leads actually to the deallocation of
       * values.
       */
      MPI_Win *values_win = nullptr;
#endif

      /**
       * Pointer to all the
       */
      std::vector<Number *> others;
    };


    template <typename Number, typename MemorySpace = MemorySpace::Host>
    class Vector : public ::dealii::LinearAlgebra::VectorSpaceVector<Number>,
                   public Subscriptor
    {
    public:
      using memory_space    = MemorySpace;
      using value_type      = Number;
      using pointer         = value_type *;
      using const_pointer   = const value_type *;
      using iterator        = value_type *;
      using const_iterator  = const value_type *;
      using reference       = value_type &;
      using const_reference = const value_type &;
      using size_type       = types::global_dof_index;
      using real_type       = typename numbers::NumberTraits<Number>::real_type;

      static_assert(
        std::is_same<MemorySpace, ::dealii::MemorySpace::Host>::value ||
          std::is_same<MemorySpace, ::dealii::MemorySpace::CUDA>::value,
        "MemorySpace should be Host or CUDA");

      /**
       * Empty constructor. To be used together with a reinit() function call.
       */
      Vector();

      /**
       * Destructor. Clear all MPI_Requests.
       */
      virtual ~Vector() override;

      /**
       * Change the dimension to that of the vector @p V. The elements of V are not
       * copied.
       *
       * @note Not implemented yet.
       */
      virtual void
      reinit(const VectorSpaceVector<Number> &V,
             const bool omit_zeroing_entries = false) override;

      /**
       * Change the dimension to that of the vector @p V. The elements of V are not
       * copied.
       */
      template <typename Number2>
      void
      reinit(const Vector<Number2, MemorySpace> &v,
             const bool                          omit_zeroing_entries = false);

      /**
       * Set up internal data structures based on the given partitioner.
       */
      void
      reinit(const std::shared_ptr<const Utilities::MPI::Partitioner> &,
             const std::shared_ptr<const PartitionerBase> &partitioner,
             const bool                                    setup_ghosts = true);

      /**
       * Get pointers to the beginning of the values of the other
       * processes of the same shared-memory domain.
       *
       * TODO: name of the function?
       */
      const std::vector<Number *> &
      other_values() const;

      /**
       * Swap the contents of this vector and the other vector @p v.
       *
       * @note Not implemented yet.
       */
      void
      swap(Vector<Number, MemorySpace> &v);

      /**
       * Copy assignment.
       */
      Vector<Number, MemorySpace> &
      operator=(const Vector<Number, MemorySpace> &in_vector);

      /**
       * Copy assignment for different underlying value types.
       */
      template <typename Number2>
      Vector<Number, MemorySpace> &
      operator=(const Vector<Number2, MemorySpace> &in_vector);

      /**
       * Update ghost values.
       */
      void
      update_ghost_values() const;

      /**
       * Start updating ghost values.
       */
      void
      update_ghost_values_start(
        const unsigned int communication_channel = 0) const;

      /**
       * Finish updating ghost values.
       */
      void
      update_ghost_values_finish() const;

      /**
       * Perform compression.
       */
      virtual void
      compress(::dealii::VectorOperation::values operation) override;

      /**
       * Start compression.
       */
      void
      compress_start(
        const unsigned int                communication_channel = 0,
        ::dealii::VectorOperation::values operation = VectorOperation::add);

      /**
       * Finish compression.
       */
      void
      compress_finish(::dealii::VectorOperation::values operation);

      /**
       * This method zeros the entries on ghost dofs, but does not touch
       * locally owned DoFs.
       */
      void
      zero_out_ghosts() const;

      /**
       * Return whether the vector currently is in a state where ghost values
       * can be read or not.
       */
      bool
      has_ghost_elements() const;

      /**
       * This method copies the data in the locally owned range from another
       * distributed vector @p src into the calling vector.
       *
       * @note Not implemented yet.
       */
      template <typename Number2>
      void
      copy_locally_owned_data_from(const Vector<Number2, MemorySpace> &src);

      /**
       * Import all the elements present in the distributed vector @p src.
       *
       * @note Not implemented yet.
       */
      template <typename MemorySpace2>
      void
      import(const Vector<Number, MemorySpace2> &src,
             VectorOperation::values             operation);

      /**
       * Multiply the entire vector by a fixed factor.
       *
       * @note Not implemented yet.
       */
      virtual Vector<Number, MemorySpace> &
      operator*=(const Number factor) override;

      /**
       * Divide the entire vector by a fixed factor.
       *
       * @note Not implemented yet.
       */
      virtual Vector<Number, MemorySpace> &
      operator/=(const Number factor) override;

      /**
       * Add the vector @p V to the present one.
       *
       * @note Not implemented yet.
       */
      virtual Vector<Number, MemorySpace> &
      operator+=(const VectorSpaceVector<Number> &V) override;

      /**
       * Subtract the vector @p V from the present one.
       *
       * @note Not implemented yet.
       */
      virtual Vector<Number, MemorySpace> &
      operator-=(const VectorSpaceVector<Number> &V) override;

      /**
       * mport all the elements present in the vector's IndexSet from the input
       * vector @p V.
       *
       * @note Not implemented yet.
       */
      virtual void
      import(
        const LinearAlgebra::ReadWriteVector<Number> &  V,
        VectorOperation::values                         operation,
        std::shared_ptr<const CommunicationPatternBase> communication_pattern =
          std::shared_ptr<const CommunicationPatternBase>()) override;

      /**
       * Return the scalar product of two vectors.
       */
      virtual Number
      operator*(const VectorSpaceVector<Number> &V) const override;

      /**
       * Add @p a to all components. Note that @p a is a scalar not a vector.
       *
       * @note Not implemented yet.
       */
      virtual void
      add(const Number a) override;

      /**
       * Simple addition of a multiple of a vector, i.e. <tt>*this += a*V</tt>.
       *
       * @note Not implemented yet.
       */
      virtual void
      add(const Number a, const VectorSpaceVector<Number> &V) override;

      /**
       * Multiple addition of scaled vectors, i.e. <tt>*this += a*V+b*W</tt>.
       *
       * @note Not implemented yet.
       */
      virtual void
      add(const Number                     a,
          const VectorSpaceVector<Number> &V,
          const Number                     b,
          const VectorSpaceVector<Number> &W) override;

      /**
       * A collective add operation: This function adds a whole set of values
       * stored in @p values to the vector components specified by @p indices.
       *
       * @note Not implemented yet.
       */
      virtual void
      add(const std::vector<size_type> &indices,
          const std::vector<Number> &   values);

      /**
       * Scaling and simple addition of a multiple of a vector, i.e. <tt>*this =
       * s*(*this)+a*V</tt>.
       */
      virtual void
      sadd(const Number                     s,
           const Number                     a,
           const VectorSpaceVector<Number> &V) override;

      /**
       * Scale each element of this vector by the corresponding element in the
       * argument.
       */
      virtual void
      scale(const VectorSpaceVector<Number> &scaling_factors) override;

      /**
       * Assignment <tt>*this = a*V</tt>.
       */
      virtual void
      equ(const Number a, const VectorSpaceVector<Number> &V) override;

      /**
       * Return the l<sub>1</sub> norm of the vector (i.e., the sum of the
       * absolute values of all entries among all processors).
       *
       * @note Not implemented yet.
       */
      virtual real_type
      l1_norm() const override;

      /**
       * Return the $l_2$ norm of the vector (i.e., the square root of
       * the sum of the square of all entries among all processors).
       */
      virtual real_type
      l2_norm() const override;

      /**
       * Return the square of the $l_2$ norm of the vector.
       */
      real_type
      norm_sqr() const;

      /**
       * Return the maximum norm of the vector (i.e., the maximum absolute value
       * among all entries and among all processors).
       */
      virtual real_type
      linfty_norm() const override;

      /**
       * Perform a combined operation of a vector addition and a subsequent
       * inner product, returning the value of the inner product. In other
       * words, the result of this function is the same as if the user called
       * @code
       * this->add(a, V);
       * return_value = *this * W;
       * @endcode
       */
      virtual Number
      add_and_dot(const Number                     a,
                  const VectorSpaceVector<Number> &V,
                  const VectorSpaceVector<Number> &W) override;

      /**
       * Return the global size of the vector.
       */
      virtual size_type
      size() const override;

      /**
       * Return index set that describes which elements of this vector are
       * owned by the current processor.
       *
       * @note Not implemented since the underlying partitioner might not be
       *   built around indices.
       */
      virtual dealii::IndexSet
      locally_owned_elements() const override;

      /**
       * Print the vector to the output stream @p out.
       *
       * @note Not implemented yet.
       */
      virtual void
      print(std::ostream &     out,
            const unsigned int precision  = 3,
            const bool         scientific = true,
            const bool         across     = true) const override;

      /**
       * Return the memory consumption of this class in bytes.
       */
      virtual std::size_t
      memory_consumption() const override;

      /**
       * Sets all elements of the vector to the scalar @p s.
       */
      virtual Vector<Number, MemorySpace> &
      operator=(const Number s) override;

      /**
       * This is a collective add operation that adds a whole set of values
       * stored in @p values to the vector components specified by @p indices.
       */
      template <typename OtherNumber>
      void
      add(const std::vector<size_type> &       indices,
          const ::dealii::Vector<OtherNumber> &values);

      /**
       * Take an address where n_elements are stored contiguously and add them
       * into the vector.
       */
      template <typename OtherNumber>
      void
      add(const size_type    n_elements,
          const size_type *  indices,
          const OtherNumber *values);

      /**
       * Scaling and simple vector addition, i.e.  <tt>*this =
       * s*(*this)+V</tt>.
       */
      void
      sadd(const Number s, const Vector<Number, MemorySpace> &V);

      /**
       * Return the local size of the vector
       */
      size_type
      local_size() const;

      /**
       * Return iterator to the start of the locally owned elements
       * of the vector.
       */
      iterator
      begin();

      /**
       * Return constant iterator to the start of the locally owned elements
       * of the vector.
       */
      const_iterator
      begin() const;

      /**
       * Return an iterator pointing to the element past the end of the array
       * of locally owned entries.
       */
      iterator
      end();

      /**
       * Return a constant iterator pointing to the element past the end of
       * the array of the locally owned entries.
       */
      const_iterator
      end() const;

      /**
       * Indirect array read access via global indices.
       *
       * @note Not implemented yet.
       */
      Number
      operator()(const size_type global_index) const;

      /**
       * Indirect array write access via global indices.
       *
       * @note Not implemented yet.
       */
      Number &
      operator()(const size_type global_index);

      /**
       * Indirect array read access via global indices.
       *
       * @note Not implemented yet.
       */
      Number operator[](const size_type global_index) const;

      /**
       * Indirect array write access via global indices.
       *
       * @note Not implemented yet.
       */
      Number &operator[](const size_type global_index);

      /**
       * Direct array access.
       */
      Number
      local_element(const size_type local_index) const;

      /**
       * Direct array access.
       */
      Number &
      local_element(const size_type local_index);

      /**
       * Return the pointer to the underlying raw array.
       */
      Number *
      get_values() const;

      /**
       * Instead of getting individual elements of a vector via operator(),
       * this function allows getting a whole set of elements at once.
       *
       * @note Not implemented yet.
       */
      template <typename OtherNumber>
      void
      extract_subvector_to(const std::vector<size_type> &indices,
                           std::vector<OtherNumber> &    values) const;

      /**
       * Instead of getting individual elements of a vector via operator(),
       * this function allows getting a whole set of elements at once.
       *
       * @note Not implemented yet.
       */
      template <typename ForwardIterator, typename OutputIterator>
      void
      extract_subvector_to(ForwardIterator       indices_begin,
                           const ForwardIterator indices_end,
                           OutputIterator        values_begin) const;

      /**
       * Return whether the vector contains only elements with value zero.
       * This is a collective operation. This function is expensive, because
       * potentially all elements have to be checked.
       */
      virtual bool
      all_zero() const override;

      /**
       * Compute the mean value of all the entries in the vector.
       *
       * @note Not implemented yet.
       */
      virtual Number
      mean_value() const override;

      /**
       * $l_p$-norm of the vector. The pth root of the sum of the pth powers
       * of the absolute values of the elements.
       *
       * @note Not implemented yet.
       */
      real_type
      lp_norm(const real_type p) const;

      /**
       * Return a reference to the MPI communicator object in use with this
       * vector.
       */
      const MPI_Comm &
      get_mpi_communicator() const;

      /**
       * Return underlying partitioner.
       *
       * @note Not implemented yet.
       */
      const std::shared_ptr<const Utilities::MPI::Partitioner> &
      get_partitioner() const;

      /**
       * Check whether the given partitioner is compatible with the
       * partitioner used for this vector.
       *
       * @note Not implemented yet.
       */
      bool
      partitioners_are_compatible(
        const Utilities::MPI::Partitioner &part) const;

      /**
       * Check whether the given partitioner is compatible with the
       * partitioner used for this vector.
       *
       * @note Not implemented yet.
       */
      bool
      partitioners_are_globally_compatible(
        const Utilities::MPI::Partitioner &part) const;

      /**
       * Change the ghost state of this vector to @p ghosted.
       */
      void
      set_ghost_state(const bool ghosted) const;

      DeclException0(ExcVectorTypeNotCompatible);

    private:
      /**
       * Simple addition of a multiple of a vector, i.e. <tt>*this += a*V</tt>
       * without MPI communication.
       */
      void
      add_local(const Number a, const VectorSpaceVector<Number> &V);

      /**
       * Scaling and simple addition of a multiple of a vector, i.e. <tt>*this =
       * s*(*this)+a*V</tt> without MPI communication.
       */
      void
      sadd_local(const Number                     s,
                 const Number                     a,
                 const VectorSpaceVector<Number> &V);

      /**
       * Local part of the inner product of two vectors.
       */
      template <typename Number2>
      Number
      inner_product_local(const Vector<Number2, MemorySpace> &V) const;

      /**
       * Local part of norm_sqr().
       */
      real_type
      norm_sqr_local() const;

      /**
       * Local part of mean_value().
       */
      Number
      mean_value_local() const;

      /**
       * Local part of l1_norm().
       */
      real_type
      l1_norm_local() const;

      /**
       * Local part of lp_norm().
       */
      real_type
      lp_norm_local(const real_type p) const;


      /**
       * Local part of linfty_norm().
       */
      real_type
      linfty_norm_local() const;

      /**
       * Local part of the addition followed by an inner product of two
       * vectors. The same applies for complex-valued vectors as for
       * the add_and_dot() function.
       */
      Number
      add_and_dot_local(const Number                       a,
                        const Vector<Number, MemorySpace> &V,
                        const Vector<Number, MemorySpace> &W);

      /**
       * Dummy Utilities::MPI::Partitioner for compatibility purposes.
       */
      std::shared_ptr<const Utilities::MPI::Partitioner> partitioner_old;

      /**
       * Partitioner.
       */
      std::shared_ptr<const PartitionerBase> partitioner;

      /**
       * The size that is currently allocated in the val array.
       */
      bool setup_ghosts = true;

      /**
       * The size that is currently allocated in the val array.
       */
      size_type allocated_size;

      /**
       * Underlying data structure storing the local elements of this vector.
       */
      mutable MemorySpaceData<Number> data;

      /**
       * For parallel loops with TBB, this member variable stores the affinity
       * information of loops.
       */
      mutable std::shared_ptr<::dealii::parallel::internal::TBBPartitioner>
        thread_loop_partitioner;

      /**
       * Temporary storage that holds the data that is sent to this processor
       * in compress() or sent from this processor in update_ghost_values().
       */
      mutable dealii::AlignedVector<Number> import_data;

      /**
       * Stores whether the vector currently allows for reading ghost elements
       * or not.
       */
      mutable bool vector_is_ghosted;

      /**
       * A vector that collects all requests from compress().
       */
      mutable std::vector<MPI_Request> compress_requests;

      /**
       * A vector that collects all requests from update_ghost_values().
       */
      mutable std::vector<MPI_Request> update_ghost_values_requests;

      /**
       * A lock that makes sure that the compress() and update_ghost_values()
       * functions give reasonable results also when used with several threads.
       */
      mutable std::mutex mutex;

      /**
       * A helper function that clears the compress_requests and
       * update_ghost_values_requests field. Used in reinit() functions.
       */
      void
      clear_mpi_requests();

      /**
       * A helper function that is used to resize the val array.
       */
      void
      resize_val(const size_type new_allocated_size, const MPI_Comm &comm_sm);

      template <typename Number2, typename MemorySpace2>
      friend class Vector;

      template <typename Number2>
      friend class BlockVector;
    };


    /*-------------------- Inline functions ---------------------------------*/

#ifndef DOXYGEN

    namespace internal
    {
      template <typename Number, typename MemorySpace>
      struct Policy
      {
        static inline typename Vector<Number, MemorySpace>::iterator
        begin(MemorySpaceData<Number> &data)
        {
          return data.values.get();
        }

        static inline typename Vector<Number, MemorySpace>::const_iterator
        begin(const MemorySpaceData<Number> &data)
        {
          return data.values.get();
        }

        static inline Number *
        get_values(MemorySpaceData<Number> &)
        {
          Assert(false, ExcNotImplemented());
          return nullptr;
        }
      };
    } // namespace internal


    template <typename Number, typename MemorySpace>
    inline bool
    Vector<Number, MemorySpace>::has_ghost_elements() const
    {
      return vector_is_ghosted && setup_ghosts;
    }



    template <typename Number, typename MemorySpace>
    inline typename Vector<Number, MemorySpace>::size_type
    Vector<Number, MemorySpace>::size() const
    {
      Assert(false, ExcNotImplemented());
      return 0;
    }



    template <typename Number, typename MemorySpace>
    inline typename Vector<Number, MemorySpace>::size_type
    Vector<Number, MemorySpace>::local_size() const
    {
      return partitioner->local_size();
    }



    template <typename Number, typename MemorySpace>
    inline IndexSet
    Vector<Number, MemorySpace>::locally_owned_elements() const
    {
      Assert(false, ExcNotImplemented());
      IndexSet is;

      return is;
    }



    template <typename Number, typename MemorySpace>
    inline typename Vector<Number, MemorySpace>::iterator
    Vector<Number, MemorySpace>::begin()
    {
      return internal::Policy<Number, MemorySpace>::begin(data);
    }



    template <typename Number, typename MemorySpace>
    inline typename Vector<Number, MemorySpace>::const_iterator
    Vector<Number, MemorySpace>::begin() const
    {
      return internal::Policy<Number, MemorySpace>::begin(data);
    }



    template <typename Number, typename MemorySpace>
    inline typename Vector<Number, MemorySpace>::iterator
    Vector<Number, MemorySpace>::end()
    {
      return internal::Policy<Number, MemorySpace>::begin(data) +
             partitioner->local_size();
    }



    template <typename Number, typename MemorySpace>
    inline typename Vector<Number, MemorySpace>::const_iterator
    Vector<Number, MemorySpace>::end() const
    {
      return internal::Policy<Number, MemorySpace>::begin(data) +
             partitioner->local_size();
    }



    template <typename Number, typename MemorySpace>
    inline Number
    Vector<Number, MemorySpace>::operator()(const size_type global_index) const
    {
      Assert(false, ExcNotImplemented());
      return data.values[partitioner_old->global_to_local(global_index)];
    }



    template <typename Number, typename MemorySpace>
    inline Number &
    Vector<Number, MemorySpace>::operator()(const size_type global_index)
    {
      Assert(false, ExcNotImplemented());
      return data.values[partitioner_old->global_to_local(global_index)];
    }



    template <typename Number, typename MemorySpace>
    inline Number Vector<Number, MemorySpace>::
                  operator[](const size_type global_index) const
    {
      Assert(false, ExcNotImplemented());
      return operator()(global_index);
    }



    template <typename Number, typename MemorySpace>
    inline Number &Vector<Number, MemorySpace>::
                   operator[](const size_type global_index)
    {
      Assert(false, ExcNotImplemented());
      return operator()(global_index);
    }



    template <typename Number, typename MemorySpace>
    inline Number
    Vector<Number, MemorySpace>::local_element(
      const size_type local_index) const
    {
      Assert((std::is_same<MemorySpace, ::dealii::MemorySpace::Host>::value),
             ExcMessage(
               "This function is only implemented for the Host memory space"));
      AssertIndexRange(local_index,
                       partitioner->local_size() +
                         partitioner->n_ghost_indices());
      // do not allow reading a vector which is not in ghost mode
      Assert(local_index < local_size() || vector_is_ghosted == true,
             ExcMessage("You tried to read a ghost element of this vector, "
                        "but it has not imported its ghost values."));

      return data.values[local_index];
    }



    template <typename Number, typename MemorySpace>
    inline Number &
    Vector<Number, MemorySpace>::local_element(const size_type local_index)
    {
      Assert((std::is_same<MemorySpace, ::dealii::MemorySpace::Host>::value),
             ExcMessage(
               "This function is only implemented for the Host memory space"));

      AssertIndexRange(local_index,
                       partitioner->local_size() +
                         partitioner->n_ghost_indices());

      return data.values[local_index];
    }



    template <typename Number, typename MemorySpace>
    inline Number *
    Vector<Number, MemorySpace>::get_values() const
    {
      Assert(false, ExcNotImplemented());

      return internal::Policy<Number, MemorySpace>::get_values(data);
    }



    template <typename Number, typename MemorySpace>
    template <typename OtherNumber>
    inline void
    Vector<Number, MemorySpace>::extract_subvector_to(
      const std::vector<size_type> &indices,
      std::vector<OtherNumber> &    values) const
    {
      Assert(false, ExcNotImplemented());
      (void)indices;
      (void)values;
    }



    template <typename Number, typename MemorySpace>
    template <typename ForwardIterator, typename OutputIterator>
    inline void
    Vector<Number, MemorySpace>::extract_subvector_to(
      ForwardIterator       indices_begin,
      const ForwardIterator indices_end,
      OutputIterator        values_begin) const
    {
      Assert(false, ExcNotImplemented());
      (void)indices_begin;
      (void)indices_end;
      (void)values_begin;
    }



    template <typename Number, typename MemorySpace>
    template <typename OtherNumber>
    inline void
    Vector<Number, MemorySpace>::add(
      const std::vector<size_type> &       indices,
      const ::dealii::Vector<OtherNumber> &values)
    {
      Assert(false, ExcNotImplemented());
      (void)indices;
      (void)values;
    }



    template <typename Number, typename MemorySpace>
    template <typename OtherNumber>
    inline void
    Vector<Number, MemorySpace>::add(const size_type    n_elements,
                                     const size_type *  indices,
                                     const OtherNumber *values)
    {
      Assert(false, ExcNotImplemented());
      (void)n_elements;
      (void)indices;
      (void)values;
    }



    template <typename Number, typename MemorySpace>
    inline const MPI_Comm &
    Vector<Number, MemorySpace>::get_mpi_communicator() const
    {
      Assert(false, ExcNotImplemented());
      return partitioner->get_mpi_communicator();
    }



    template <typename Number, typename MemorySpace>
    inline const std::shared_ptr<const Utilities::MPI::Partitioner> &
    Vector<Number, MemorySpace>::get_partitioner() const
    {
      return partitioner_old;
    }



    template <typename Number, typename MemorySpace>
    inline void
    Vector<Number, MemorySpace>::set_ghost_state(const bool ghosted) const
    {
      vector_is_ghosted = ghosted;
    }



    template <typename Number, typename MemorySpace>
    const std::vector<Number *> &
    Vector<Number, MemorySpace>::other_values() const
    {
      return data.others;
    }

#endif

  } // namespace SharedMPI
} // namespace LinearAlgebra


template <typename Number, typename MemorySpace>
inline void
swap(LinearAlgebra::SharedMPI::Vector<Number, MemorySpace> &u,
     LinearAlgebra::SharedMPI::Vector<Number, MemorySpace> &v)
{
  Assert(false, ExcNotImplemented());
  (void)u;
  (void)v;
}


template <typename Number, typename MemorySpace>
struct is_serial_vector<LinearAlgebra::SharedMPI::Vector<Number, MemorySpace>>
  : std::false_type
{};



namespace internal
{
  namespace LinearOperatorImplementation
  {
    template <typename>
    class ReinitHelper;

    template <typename Number>
    class ReinitHelper<LinearAlgebra::SharedMPI::Vector<Number>>
    {
    public:
      template <typename Matrix>
      static void
      reinit_range_vector(const Matrix &                            matrix,
                          LinearAlgebra::SharedMPI::Vector<Number> &v,
                          bool omit_zeroing_entries)
      {
        Assert(false, ExcNotImplemented());
        (void)matrix;
        (void)v;
        (void)omit_zeroing_entries;
      }

      template <typename Matrix>
      static void
      reinit_domain_vector(const Matrix &                            matrix,
                           LinearAlgebra::SharedMPI::Vector<Number> &v,
                           bool omit_zeroing_entries)
      {
        Assert(false, ExcNotImplemented());
        (void)matrix;
        (void)v;
        (void)omit_zeroing_entries;
      }
    };

  } // namespace LinearOperatorImplementation
} /* namespace internal */


DEAL_II_NAMESPACE_CLOSE

#endif
