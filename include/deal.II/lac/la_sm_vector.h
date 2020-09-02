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

      Vector(const size_type size);

      Vector(const IndexSet &local_range,
             const IndexSet &ghost_indices,
             const MPI_Comm  communicator);

      Vector(const IndexSet &local_range, const MPI_Comm communicator);

      Vector(const std::shared_ptr<const Utilities::MPI::Partitioner>
               &partitioner_old);

      /**
       * Destructor. Clear all MPI_Requests.
       */
      virtual ~Vector() override;

      void
      reinit(const size_type size, const bool omit_zeroing_entries = false);

      template <typename Number2>
      void
      reinit(const Vector<Number2, MemorySpace> &in_vector,
             const bool                          omit_zeroing_entries = false);

      void
      reinit(const IndexSet &local_range,
             const IndexSet &ghost_indices,
             const MPI_Comm  communicator);

      void
      reinit(const IndexSet &local_range, const MPI_Comm communicator);

      void
      reinit(const std::shared_ptr<const Utilities::MPI::Partitioner>
               &                                           partitioner_old,
             const std::shared_ptr<const PartitionerBase> &partitioner,
             const bool                                    setup_ghosts = true);

      void
      swap(Vector<Number, MemorySpace> &v);

      Vector<Number, MemorySpace> &
      operator=(const Vector<Number, MemorySpace> &in_vector);

      template <typename Number2>
      Vector<Number, MemorySpace> &
      operator=(const Vector<Number2, MemorySpace> &in_vector);

      virtual void
      compress(::dealii::VectorOperation::values operation) override;

      void
      update_ghost_values() const;

      void
      compress_start(
        const unsigned int                communication_channel = 0,
        ::dealii::VectorOperation::values operation = VectorOperation::add);

      void
      compress_finish(::dealii::VectorOperation::values operation);

      void
      update_ghost_values_start(
        const unsigned int communication_channel = 0) const;

      void
      update_ghost_values_finish() const;

      void
      zero_out_ghosts() const;

      bool
      has_ghost_elements() const;

      template <typename Number2>
      void
      copy_locally_owned_data_from(const Vector<Number2, MemorySpace> &src);

      template <typename MemorySpace2>
      void
      import(const Vector<Number, MemorySpace2> &src,
             VectorOperation::values             operation);

      virtual void
      reinit(const VectorSpaceVector<Number> &V,
             const bool omit_zeroing_entries = false) override;

      virtual Vector<Number, MemorySpace> &
      operator*=(const Number factor) override;

      virtual Vector<Number, MemorySpace> &
      operator/=(const Number factor) override;

      virtual Vector<Number, MemorySpace> &
      operator+=(const VectorSpaceVector<Number> &V) override;

      virtual Vector<Number, MemorySpace> &
      operator-=(const VectorSpaceVector<Number> &V) override;

      virtual void
      import(
        const LinearAlgebra::ReadWriteVector<Number> &  V,
        VectorOperation::values                         operation,
        std::shared_ptr<const CommunicationPatternBase> communication_pattern =
          std::shared_ptr<const CommunicationPatternBase>()) override;

      virtual Number
      operator*(const VectorSpaceVector<Number> &V) const override;

      virtual void
      add(const Number a) override;

      virtual void
      add(const Number a, const VectorSpaceVector<Number> &V) override;

      virtual void
      add(const Number                     a,
          const VectorSpaceVector<Number> &V,
          const Number                     b,
          const VectorSpaceVector<Number> &W) override;

      virtual void
      add(const std::vector<size_type> &indices,
          const std::vector<Number> &   values);

      virtual void
      sadd(const Number                     s,
           const Number                     a,
           const VectorSpaceVector<Number> &V) override;

      virtual void
      scale(const VectorSpaceVector<Number> &scaling_factors) override;

      virtual void
      equ(const Number a, const VectorSpaceVector<Number> &V) override;

      virtual real_type
      l1_norm() const override;

      virtual real_type
      l2_norm() const override;

      real_type
      norm_sqr() const;

      virtual real_type
      linfty_norm() const override;

      virtual Number
      add_and_dot(const Number                     a,
                  const VectorSpaceVector<Number> &V,
                  const VectorSpaceVector<Number> &W) override;

      virtual size_type
      size() const override;

      virtual dealii::IndexSet
      locally_owned_elements() const override;

      virtual void
      print(std::ostream &     out,
            const unsigned int precision  = 3,
            const bool         scientific = true,
            const bool         across     = true) const override;

      virtual std::size_t
      memory_consumption() const override;

      virtual Vector<Number, MemorySpace> &
      operator=(const Number s) override;

      template <typename OtherNumber>
      void
      add(const std::vector<size_type> &       indices,
          const ::dealii::Vector<OtherNumber> &values);

      template <typename OtherNumber>
      void
      add(const size_type    n_elements,
          const size_type *  indices,
          const OtherNumber *values);

      void
      sadd(const Number s, const Vector<Number, MemorySpace> &V);

      size_type
      local_size() const;

      Number *
      begin_sm();

      Number *
      begin_sm() const;

      /**
       * Get const pointers to the beginning of the values of the other
       * processes of the same shared-memory domain.
       *
       * TODO: name of the function?
       */
      std::vector<Number *> &
      other_values();

      /**
       * Get pointers to the beginning of the values of the other
       * processes of the same shared-memory domain.
       *
       * TODO: name of the function?
       */
      const std::vector<Number *> &
      other_values() const;

      iterator
      begin();

      const_iterator
      begin() const;

      iterator
      end();

      const_iterator
      end() const;

      Number
      operator()(const size_type global_index) const;

      Number &
      operator()(const size_type global_index);

      Number operator[](const size_type global_index) const;

      Number &operator[](const size_type global_index);

      Number
      local_element(const size_type local_index) const;

      Number &
      local_element(const size_type local_index);

      Number *
      get_values() const;

      template <typename OtherNumber>
      void
      extract_subvector_to(const std::vector<size_type> &indices,
                           std::vector<OtherNumber> &    values) const;

      template <typename ForwardIterator, typename OutputIterator>
      void
      extract_subvector_to(ForwardIterator       indices_begin,
                           const ForwardIterator indices_end,
                           OutputIterator        values_begin) const;

      virtual bool
      all_zero() const override;

      virtual Number
      mean_value() const override;

      real_type
      lp_norm(const real_type p) const;

      const MPI_Comm &
      get_mpi_communicator() const;

      const std::shared_ptr<const Utilities::MPI::Partitioner> &
      get_partitioner() const;

      bool
      partitioners_are_compatible(
        const Utilities::MPI::Partitioner &part) const;

      bool
      partitioners_are_globally_compatible(
        const Utilities::MPI::Partitioner &part) const;

      void
      set_ghost_state(const bool ghosted) const;

      DeclException0(ExcVectorTypeNotCompatible);

      DeclException0(ExcNotAllowedForCuda);

      DeclException3(ExcNonMatchingElements,
                     Number,
                     Number,
                     unsigned int,
                     << "Called compress(VectorOperation::insert), but"
                     << " the element received from a remote processor, value "
                     << std::setprecision(16) << arg1
                     << ", does not match with the value "
                     << std::setprecision(16) << arg2
                     << " on the owner processor " << arg3);

      DeclException4(ExcAccessToNonLocalElement,
                     size_type,
                     size_type,
                     size_type,
                     size_type,
                     << "You tried to access element " << arg1
                     << " of a SharedMPI vector, but this element is not "
                     << "stored on the current processor. Note: The range of "
                     << "locally owned elements is " << arg2 << " to " << arg3
                     << ", and there are " << arg4 << " ghost elements "
                     << "that this vector can access.");

    private:
      void
      add_local(const Number a, const VectorSpaceVector<Number> &V);

      void
      sadd_local(const Number                     s,
                 const Number                     a,
                 const VectorSpaceVector<Number> &V);

      template <typename Number2>
      Number
      inner_product_local(const Vector<Number2, MemorySpace> &V) const;

      real_type
      norm_sqr_local() const;

      Number
      mean_value_local() const;

      real_type
      l1_norm_local() const;

      real_type
      lp_norm_local(const real_type p) const;

      real_type
      linfty_norm_local() const;

      Number
      add_and_dot_local(const Number                       a,
                        const Vector<Number, MemorySpace> &V,
                        const Vector<Number, MemorySpace> &W);

      std::shared_ptr<const Utilities::MPI::Partitioner> partitioner_old;

      std::shared_ptr<const PartitionerBase> partitioner;

      bool setup_ghosts = true;

      size_type allocated_size;

      mutable MemorySpaceData<Number> data;

      /**
       * For parallel loops with TBB, this member variable stores the affinity
       * information of loops.
       */
      mutable std::shared_ptr<::dealii::parallel::internal::TBBPartitioner>
        thread_loop_partitioner;

      // needed?
      mutable dealii::AlignedVector<Number> import_data;

      mutable bool vector_is_ghosted;

      mutable std::vector<MPI_Request> requests;

      mutable std::mutex mutex;

      void
      clear_mpi_requests();

      void
      resize_val(const size_type new_allocated_size, const MPI_Comm &comm_sm);

      template <typename Number2, typename MemorySpace2>
      friend class Vector;

      template <typename Number2>
      friend class BlockVector;

      // std::vector<Number *> data_others;
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
        begin_sm(MemorySpaceData<Number> &data)
        {
          return data.others[0];
        }

        static inline Number *
        begin_sm(const MemorySpaceData<Number> &data)
        {
          return data.others[0];
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
    Number *
    Vector<Number, MemorySpace>::begin_sm()
    {
      return internal::Policy<Number, MemorySpace>::begin_sm(data);
    }

    template <typename Number, typename MemorySpace>
    Number *
    Vector<Number, MemorySpace>::begin_sm() const
    {
      return internal::Policy<Number, MemorySpace>::begin_sm(data);
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
    std::vector<Number *> &
    Vector<Number, MemorySpace>::other_values()
    {
      return data.others;
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
