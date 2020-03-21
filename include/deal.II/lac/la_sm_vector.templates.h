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

#ifndef dealii_la_sm_vector_templates_h
#define dealii_la_sm_vector_templates_h


#include <deal.II/base/config.h>

#include <deal.II/base/cuda.h>
#include <deal.II/base/cuda_size.h>
#include <deal.II/base/std_cxx14/memory.h>

#include <deal.II/lac/exceptions.h>
#include <deal.II/lac/la_sm_vector.h>
#include <deal.II/lac/petsc_vector.h>
#include <deal.II/lac/read_write_vector.h>
#include <deal.II/lac/trilinos_vector.h>
#include <deal.II/lac/vector_operations_internal.h>


DEAL_II_NAMESPACE_OPEN


namespace LinearAlgebra
{
  namespace SharedMPI
  {
    namespace internal
    {
      template <typename Number, typename MemorySpaceType>
      struct la_parallel_vector_templates_functions
      {
        using size_type = types::global_dof_index;

        static void
        resize_val(const types::global_dof_index new_alloc_size,
                   types::global_dof_index &     allocated_size,
                   MemorySpaceData<Number> &     data,
                   const MPI_Comm &              comm_shared)
        {
          Assert(allocated_size == 0, ExcNotImplemented());
          Assert(data.values == nullptr, ExcNotImplemented());

          allocated_size = new_alloc_size;

          // TODO
          std::vector<Number *> data_others;

          data.values_win   = new MPI_Win;
          Number *data_this = (Number *)malloc(0);
          data_others.resize(Utilities::MPI::n_mpi_processes(comm_shared));


          MPI_Info info;
          MPI_Info_create(&info);
          MPI_Info_set(info, "alloc_shared_noncontig", "true");

          MPI_Win_allocate_shared(new_alloc_size * sizeof(Number),
                                  sizeof(Number),
                                  info,
                                  comm_shared,
                                  data_this,
                                  data.values_win);

          for (unsigned int i = 0;
               i < Utilities::MPI::n_mpi_processes(comm_shared);
               i++)
            {
              int      disp_unit;
              MPI_Aint ssize;
              MPI_Win_shared_query(
                *data.values_win, i, &ssize, &disp_unit, &data_others[i]);
            }

          data.values = {
            data_others[Utilities::MPI::this_mpi_process(comm_shared)],
            [&data](Number *&) { MPI_Win_free(data.values_win); }};
        }
      };
    } // namespace internal


    template <typename Number, typename MemorySpaceType>
    void
    Vector<Number, MemorySpaceType>::clear_mpi_requests()
    {
      // TODO: what is the use of this?
      // it should probably delegated to partitioner
    }



    template <typename Number, typename MemorySpaceType>
    void
    Vector<Number, MemorySpaceType>::resize_val(const size_type new_alloc_size,
                                                const MPI_Comm &comm_sm)
    {
      internal::la_parallel_vector_templates_functions<
        Number,
        MemorySpaceType>::resize_val(new_alloc_size,
                                     allocated_size,
                                     data,
                                     comm_sm);
    }



    template <typename Number, typename MemorySpaceType>
    void
    Vector<Number, MemorySpaceType>::reinit(const size_type size,
                                            const bool omit_zeroing_entries)
    {
      Assert(false, ExcNotImplemented());
      (void)size;
      (void)omit_zeroing_entries;
    }



    template <typename Number, typename MemorySpaceType>
    template <typename Number2>
    void
    Vector<Number, MemorySpaceType>::reinit(
      const Vector<Number2, MemorySpaceType> &v,
      const bool                              omit_zeroing_entries)
    {
      Assert(false, ExcNotImplemented());
      (void)v;
      (void)omit_zeroing_entries;
    }



    template <typename Number, typename MemorySpaceType>
    void
    Vector<Number, MemorySpaceType>::reinit(
      const IndexSet &locally_owned_indices,
      const IndexSet &ghost_indices,
      const MPI_Comm  communicator)
    {
      Assert(false, ExcNotImplemented());
      (void)locally_owned_indices;
      (void)ghost_indices;
      (void)communicator;
    }



    template <typename Number, typename MemorySpaceType>
    void
    Vector<Number, MemorySpaceType>::reinit(
      const IndexSet &locally_owned_indices,
      const MPI_Comm  communicator)
    {
      Assert(false, ExcNotImplemented());
      (void)locally_owned_indices;
      (void)communicator;
    }



    template <typename Number, typename MemorySpaceType>
    void
    Vector<Number, MemorySpaceType>::reinit(
      const std::shared_ptr<const Utilities::MPI::Partitioner> &partitioner,
      const std::shared_ptr<const Partitioner<Number>> &        partitioner_sm)
    {
      clear_mpi_requests();
      this->partitioner    = partitioner;
      this->partitioner_sm = partitioner_sm;

      // set vector size and allocate memory
      const size_type new_allocated_size =
        this->partitioner->local_size() + this->partitioner->n_ghost_indices();
      resize_val(new_allocated_size, partitioner_sm->get_sm_mpi_communicator());

      // initialize to zero
      this->operator=(Number());


      // do not reallocate import_data directly, but only upon request. It
      // is only used as temporary storage for compress() and
      // update_ghost_values, and we might have vectors where we never
      // call these methods and hence do not need to have the storage.
      import_data.values.reset();
      import_data.values_dev.reset();

      vector_is_ghosted = false;
    }



    template <typename Number, typename MemorySpaceType>
    Vector<Number, MemorySpaceType>::Vector()
      : partitioner(new Utilities::MPI::Partitioner())
      , allocated_size(0)
    {}



    template <typename Number, typename MemorySpaceType>
    Vector<Number, MemorySpaceType>::Vector(
      const Vector<Number, MemorySpaceType> &v)
      : Subscriptor()
      , allocated_size(0)
      , vector_is_ghosted(false)
    {
      Assert(false, ExcNotImplemented());
      (void)v;
    }



    template <typename Number, typename MemorySpaceType>
    Vector<Number, MemorySpaceType>::Vector(const IndexSet &local_range,
                                            const IndexSet &ghost_indices,
                                            const MPI_Comm  communicator)
      : allocated_size(0)
      , vector_is_ghosted(false)
    {
      Assert(false, ExcNotImplemented());
      (void)local_range;
      (void)ghost_indices;
      (void)communicator;
    }



    template <typename Number, typename MemorySpaceType>
    Vector<Number, MemorySpaceType>::Vector(const IndexSet &local_range,
                                            const MPI_Comm  communicator)
      : allocated_size(0)
      , vector_is_ghosted(false)
    {
      Assert(false, ExcNotImplemented());
      (void)local_range;
      (void)communicator;
    }



    template <typename Number, typename MemorySpaceType>
    Vector<Number, MemorySpaceType>::Vector(const size_type size)
      : allocated_size(0)
      , vector_is_ghosted(false)
    {
      Assert(false, ExcNotImplemented());
      (void)size;
    }



    template <typename Number, typename MemorySpaceType>
    Vector<Number, MemorySpaceType>::Vector(
      const std::shared_ptr<const Utilities::MPI::Partitioner> &partitioner)
      : allocated_size(0)
      , vector_is_ghosted(false)
    {
      Assert(false, ExcNotImplemented());
      (void)partitioner;
    }



    template <typename Number, typename MemorySpaceType>
    inline Vector<Number, MemorySpaceType>::~Vector()
    {
      try
        {
          clear_mpi_requests();
        }
      catch (...)
        {}
    }



    template <typename Number, typename MemorySpaceType>
    inline Vector<Number, MemorySpaceType> &
    Vector<Number, MemorySpaceType>::
    operator=(const Vector<Number, MemorySpaceType> &c)
    {
      Assert(false, ExcNotImplemented());
      (void)c;
      return *this;
    }



    template <typename Number, typename MemorySpaceType>
    template <typename Number2>
    inline Vector<Number, MemorySpaceType> &
    Vector<Number, MemorySpaceType>::
    operator=(const Vector<Number2, MemorySpaceType> &c)
    {
      Assert(false, ExcNotImplemented());
      (void)c;
      return *this;
    }



    template <typename Number, typename MemorySpaceType>
    template <typename Number2>
    void
    Vector<Number, MemorySpaceType>::copy_locally_owned_data_from(
      const Vector<Number2, MemorySpaceType> &src)
    {
      Assert(false, ExcNotImplemented());
      (void)src;
    }



    template <typename Number, typename MemorySpaceType>
    template <typename MemorySpaceType2>
    void
    Vector<Number, MemorySpaceType>::import(
      const Vector<Number, MemorySpaceType2> &src,
      VectorOperation::values                 operation)
    {
      Assert(false, ExcNotImplemented());
      (void)src;
      (void)operation;
    }



    template <typename Number, typename MemorySpaceType>
    void
    Vector<Number, MemorySpaceType>::compress(
      ::dealii::VectorOperation::values operation)
    {
      Assert(false, ExcNotImplemented());
      compress_start(0, operation);
      compress_finish(operation);
    }



    template <typename Number, typename MemorySpaceType>
    void
    Vector<Number, MemorySpaceType>::update_ghost_values() const
    {
      Assert(false, ExcNotImplemented());
      update_ghost_values_start();
      update_ghost_values_finish();
    }



    template <typename Number, typename MemorySpaceType>
    void
    Vector<Number, MemorySpaceType>::zero_out_ghosts() const
    {
      if (data.values != nullptr)
        std::fill_n(data.values.get() + partitioner->local_size(),
                    partitioner->n_ghost_indices(),
                    Number());
    }



    template <typename Number, typename MemorySpaceType>
    void
    Vector<Number, MemorySpaceType>::compress_start(
      const unsigned int                communication_channel,
      ::dealii::VectorOperation::values operation)
    {
      Assert(false, ExcNotImplemented());
      (void)communication_channel;
      (void)operation;
    }



    template <typename Number, typename MemorySpaceType>
    void
    Vector<Number, MemorySpaceType>::compress_finish(
      ::dealii::VectorOperation::values operation)
    {
      Assert(false, ExcNotImplemented());
      (void)operation;
    }



    template <typename Number, typename MemorySpaceType>
    void
    Vector<Number, MemorySpaceType>::update_ghost_values_start(
      const unsigned int communication_channel) const
    {
      Assert(false, ExcNotImplemented());
      (void)communication_channel;
    }



    template <typename Number, typename MemorySpaceType>
    void
    Vector<Number, MemorySpaceType>::update_ghost_values_finish() const
    {
      Assert(false, ExcNotImplemented());
    }



    template <typename Number, typename MemorySpaceType>
    void
    Vector<Number, MemorySpaceType>::import(
      const ReadWriteVector<Number> &                 V,
      VectorOperation::values                         operation,
      std::shared_ptr<const CommunicationPatternBase> communication_pattern)
    {
      Assert(false, ExcNotImplemented());
      (void)V;
      (void)operation;
      (void)communication_pattern;
    }

    template <typename Number, typename MemorySpaceType>
    void
    Vector<Number, MemorySpaceType>::swap(Vector<Number, MemorySpaceType> &v)
    {
      Assert(false, ExcNotImplemented());
      (void)v;
    }



    template <typename Number, typename MemorySpaceType>
    Vector<Number, MemorySpaceType> &
    Vector<Number, MemorySpaceType>::operator=(const Number s)
    {
      const size_type this_size = local_size();
      if (this_size > 0)
        std::fill_n(data.values.get(), this_size, s);

      if (s == Number())
        zero_out_ghosts();

      return *this;
    }



    template <typename Number, typename MemorySpaceType>
    void
    Vector<Number, MemorySpaceType>::reinit(const VectorSpaceVector<Number> &V,
                                            const bool omit_zeroing_entries)
    {
      Assert(false, ExcNotImplemented());
      (void)V;
      (void)omit_zeroing_entries;
    }



    template <typename Number, typename MemorySpaceType>
    Vector<Number, MemorySpaceType> &
    Vector<Number, MemorySpaceType>::
    operator+=(const VectorSpaceVector<Number> &vv)
    {
      Assert(false, ExcNotImplemented());
      (void)vv;

      return *this;
    }



    template <typename Number, typename MemorySpaceType>
    Vector<Number, MemorySpaceType> &
    Vector<Number, MemorySpaceType>::
    operator-=(const VectorSpaceVector<Number> &vv)
    {
      Assert(false, ExcNotImplemented());
      (void)vv;

      return *this;
    }



    template <typename Number, typename MemorySpaceType>
    void
    Vector<Number, MemorySpaceType>::add(const Number a)
    {
      Assert(false, ExcNotImplemented());
      (void)a;
    }



    template <typename Number, typename MemorySpaceType>
    void
    Vector<Number, MemorySpaceType>::add_local(
      const Number                     a,
      const VectorSpaceVector<Number> &vv)
    {
      Assert(false, ExcNotImplemented());
      (void)a;
      (void)vv;
    }



    template <typename Number, typename MemorySpaceType>
    void
    Vector<Number, MemorySpaceType>::add(const Number                     a,
                                         const VectorSpaceVector<Number> &vv)
    {
      Assert(false, ExcNotImplemented());
      (void)a;
      (void)vv;
    }



    template <typename Number, typename MemorySpaceType>
    void
    Vector<Number, MemorySpaceType>::add(const Number                     a,
                                         const VectorSpaceVector<Number> &vv,
                                         const Number                     b,
                                         const VectorSpaceVector<Number> &ww)
    {
      Assert(false, ExcNotImplemented());
      (void)a;
      (void)vv;
      (void)b;
      (void)ww;
    }



    template <typename Number, typename MemorySpaceType>
    void
    Vector<Number, MemorySpaceType>::add(const std::vector<size_type> &indices,
                                         const std::vector<Number> &   values)
    {
      Assert(false, ExcNotImplemented());
      (void)indices;
      (void)values;
    }



    template <typename Number, typename MemorySpaceType>
    void
    Vector<Number, MemorySpaceType>::sadd(
      const Number                           x,
      const Vector<Number, MemorySpaceType> &v)
    {
      Assert(false, ExcNotImplemented());
      (void)x;
      (void)v;
    }



    template <typename Number, typename MemorySpaceType>
    void
    Vector<Number, MemorySpaceType>::sadd_local(
      const Number                     x,
      const Number                     a,
      const VectorSpaceVector<Number> &vv)
    {
      Assert(false, ExcNotImplemented());
      (void)x;
      (void)a;
      (void)vv;
    }



    template <typename Number, typename MemorySpaceType>
    void
    Vector<Number, MemorySpaceType>::sadd(const Number                     x,
                                          const Number                     a,
                                          const VectorSpaceVector<Number> &vv)
    {
      Assert(false, ExcNotImplemented());
      (void)x;
      (void)a;
      (void)vv;
    }



    template <typename Number, typename MemorySpaceType>
    void
    Vector<Number, MemorySpaceType>::sadd(
      const Number                           x,
      const Number                           a,
      const Vector<Number, MemorySpaceType> &v,
      const Number                           b,
      const Vector<Number, MemorySpaceType> &w)
    {
      Assert(false, ExcNotImplemented());
      (void)x;
      (void)a;
      (void)v;
      (void)b;
      (void)w;
    }



    template <typename Number, typename MemorySpaceType>
    Vector<Number, MemorySpaceType> &
    Vector<Number, MemorySpaceType>::operator*=(const Number factor)
    {
      Assert(false, ExcNotImplemented());
      (void)factor;

      return *this;
    }



    template <typename Number, typename MemorySpaceType>
    Vector<Number, MemorySpaceType> &
    Vector<Number, MemorySpaceType>::operator/=(const Number factor)
    {
      Assert(false, ExcNotImplemented());
      (void)factor;

      return *this;
    }



    template <typename Number, typename MemorySpaceType>
    void
    Vector<Number, MemorySpaceType>::scale(const VectorSpaceVector<Number> &vv)
    {
      Assert(false, ExcNotImplemented());
      (void)vv;
    }



    template <typename Number, typename MemorySpaceType>
    void
    Vector<Number, MemorySpaceType>::equ(const Number                     a,
                                         const VectorSpaceVector<Number> &vv)
    {
      Assert(false, ExcNotImplemented());
      (void)a;
      (void)vv;
    }



    template <typename Number, typename MemorySpaceType>
    void
    Vector<Number, MemorySpaceType>::equ(
      const Number                           a,
      const Vector<Number, MemorySpaceType> &v,
      const Number                           b,
      const Vector<Number, MemorySpaceType> &w)
    {
      Assert(false, ExcNotImplemented());
      (void)a;
      (void)v;
      (void)b;
      (void)w;
    }



    template <typename Number, typename MemorySpaceType>
    bool
    Vector<Number, MemorySpaceType>::all_zero() const
    {
      Assert(false, ExcNotImplemented());
      return false;
    }



    template <typename Number, typename MemorySpaceType>
    template <typename Number2>
    Number
    Vector<Number, MemorySpaceType>::inner_product_local(
      const Vector<Number2, MemorySpaceType> &v) const
    {
      Assert(false, ExcNotImplemented());
      (void)v;
    }



    template <typename Number, typename MemorySpaceType>
    Number Vector<Number, MemorySpaceType>::
           operator*(const VectorSpaceVector<Number> &vv) const
    {
      Assert(false, ExcNotImplemented());
      (void)vv;
      return 0;
    }



    template <typename Number, typename MemorySpaceType>
    typename Vector<Number, MemorySpaceType>::real_type
    Vector<Number, MemorySpaceType>::norm_sqr_local() const
    {
      Assert(false, ExcNotImplemented());
      return 0;
    }



    template <typename Number, typename MemorySpaceType>
    Number
    Vector<Number, MemorySpaceType>::mean_value_local() const
    {
      Assert(false, ExcNotImplemented());
      return 0;
    }



    template <typename Number, typename MemorySpaceType>
    Number
    Vector<Number, MemorySpaceType>::mean_value() const
    {
      Assert(false, ExcNotImplemented());
      return 0;
    }



    template <typename Number, typename MemorySpaceType>
    typename Vector<Number, MemorySpaceType>::real_type
    Vector<Number, MemorySpaceType>::l1_norm_local() const
    {
      Assert(false, ExcNotImplemented());
      return 0;
    }



    template <typename Number, typename MemorySpaceType>
    typename Vector<Number, MemorySpaceType>::real_type
    Vector<Number, MemorySpaceType>::l1_norm() const
    {
      Assert(false, ExcNotImplemented());
      return 0;
    }



    template <typename Number, typename MemorySpaceType>
    typename Vector<Number, MemorySpaceType>::real_type
    Vector<Number, MemorySpaceType>::norm_sqr() const
    {
      Assert(false, ExcNotImplemented());
      return 0;
    }



    template <typename Number, typename MemorySpaceType>
    typename Vector<Number, MemorySpaceType>::real_type
    Vector<Number, MemorySpaceType>::l2_norm() const
    {
      Assert(false, ExcNotImplemented());
      return 0;
    }



    template <typename Number, typename MemorySpaceType>
    typename Vector<Number, MemorySpaceType>::real_type
    Vector<Number, MemorySpaceType>::lp_norm_local(const real_type p) const
    {
      Assert(false, ExcNotImplemented());
      (void)p;
      return 0;
    }



    template <typename Number, typename MemorySpaceType>
    typename Vector<Number, MemorySpaceType>::real_type
    Vector<Number, MemorySpaceType>::lp_norm(const real_type p) const
    {
      Assert(false, ExcNotImplemented());
      (void)p;
      return 0;
    }



    template <typename Number, typename MemorySpaceType>
    typename Vector<Number, MemorySpaceType>::real_type
    Vector<Number, MemorySpaceType>::linfty_norm_local() const
    {
      Assert(false, ExcNotImplemented());
      return 0;
    }



    template <typename Number, typename MemorySpaceType>
    inline typename Vector<Number, MemorySpaceType>::real_type
    Vector<Number, MemorySpaceType>::linfty_norm() const
    {
      Assert(false, ExcNotImplemented());
      return 0;
    }



    template <typename Number, typename MemorySpaceType>
    Number
    Vector<Number, MemorySpaceType>::add_and_dot_local(
      const Number                           a,
      const Vector<Number, MemorySpaceType> &v,
      const Vector<Number, MemorySpaceType> &w)
    {
      Assert(false, ExcNotImplemented());
      (void)a;
      (void)v;
      (void)w;

      return 0;
    }



    template <typename Number, typename MemorySpaceType>
    Number
    Vector<Number, MemorySpaceType>::add_and_dot(
      const Number                     a,
      const VectorSpaceVector<Number> &vv,
      const VectorSpaceVector<Number> &ww)
    {
      Assert(false, ExcNotImplemented());
      (void)a;
      (void)vv;
      (void)ww;
      return 0;
    }



    template <typename Number, typename MemorySpaceType>
    inline bool
    Vector<Number, MemorySpaceType>::partitioners_are_compatible(
      const Utilities::MPI::Partitioner &part) const
    {
      Assert(false, ExcNotImplemented());
      return partitioner->is_compatible(part);
    }



    template <typename Number, typename MemorySpaceType>
    inline bool
    Vector<Number, MemorySpaceType>::partitioners_are_globally_compatible(
      const Utilities::MPI::Partitioner &part) const
    {
      Assert(false, ExcNotImplemented());
      return partitioner->is_globally_compatible(part);
    }



    template <typename Number, typename MemorySpaceType>
    std::size_t
    Vector<Number, MemorySpaceType>::memory_consumption() const
    {
      Assert(false, ExcNotImplemented());
      return 0;
    }



    template <typename Number, typename MemorySpaceType>
    void
    Vector<Number, MemorySpaceType>::print(std::ostream &     out,
                                           const unsigned int precision,
                                           const bool         scientific,
                                           const bool         across) const
    {
      Assert(false, ExcNotImplemented());
      (void)out;
      (void)precision;
      (void)scientific;
      (void)across;
    }

  } // end of namespace SharedMPI
} // end of namespace LinearAlgebra


DEAL_II_NAMESPACE_CLOSE

#endif
