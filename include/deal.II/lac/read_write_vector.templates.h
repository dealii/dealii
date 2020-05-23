// ---------------------------------------------------------------------
//
// Copyright (C) 2015 - 2020 by the deal.II authors
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

#ifndef dealii_parallel_vector_templates_h
#define dealii_parallel_vector_templates_h


#include <deal.II/base/config.h>

#include <deal.II/base/partitioner.h>

#include <deal.II/lac/exceptions.h>
#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/lac/read_write_vector.h>
#include <deal.II/lac/vector_operations_internal.h>

#include <boost/io/ios_state.hpp>

#ifdef DEAL_II_WITH_PETSC
#  include <deal.II/lac/petsc_vector.h>
#endif

#ifdef DEAL_II_WITH_TRILINOS
#  include <deal.II/lac/trilinos_epetra_communication_pattern.h>
#  include <deal.II/lac/trilinos_vector.h>

#  include <Epetra_Import.h>
#endif

#ifdef DEAL_II_WITH_CUDA
#  include <deal.II/lac/cuda_vector.h>
#endif

DEAL_II_NAMESPACE_OPEN


namespace LinearAlgebra
{
  namespace internal
  {
    // In the import_from_ghosted_array_finish we need to calculate the
    // maximal and minimal value for the given number type, which is not
    // straight forward for complex numbers. Therefore, comparison of complex
    // numbers is prohibited and throws an assert.
    template <typename Number>
    Number
    get_min(const Number a, const Number b)
    {
      return std::min(a, b);
    }

    template <typename Number>
    std::complex<Number>
    get_min(const std::complex<Number> a, const std::complex<Number>)
    {
      AssertThrow(false,
                  ExcMessage("VectorOperation::min not "
                             "implemented for complex numbers"));
      return a;
    }

    template <typename Number>
    Number
    get_max(const Number a, const Number b)
    {
      return std::max(a, b);
    }

    template <typename Number>
    std::complex<Number>
    get_max(const std::complex<Number> a, const std::complex<Number>)
    {
      AssertThrow(false,
                  ExcMessage("VectorOperation::max not "
                             "implemented for complex numbers"));
      return a;
    }



    template <typename Number, typename MemorySpace>
    struct read_write_vector_functions
    {
      static void
      import(const std::shared_ptr<const ::dealii::Utilities::MPI::Partitioner>
               &                                     communication_pattern,
             const Number *                          values,
             const ::dealii::VectorOperation::values operation,
             ::dealii::LinearAlgebra::ReadWriteVector<Number> &rw_vector)
      {
        (void)communication_pattern;
        (void)values;
        (void)operation;
        (void)rw_vector;

        static_assert(
          std::is_same<MemorySpace, ::dealii::MemorySpace::Host>::value ||
            std::is_same<MemorySpace, ::dealii::MemorySpace::CUDA>::value,
          "MemorySpace should be Host or CUDA");
      }
    };



    template <typename Number>
    struct read_write_vector_functions<Number, ::dealii::MemorySpace::Host>
    {
      using size_type = types::global_dof_index;


      static void
      import(const std::shared_ptr<const ::dealii::Utilities::MPI::Partitioner>
               &                                     communication_pattern,
             const Number *                          values,
             const ::dealii::VectorOperation::values operation,
             ::dealii::LinearAlgebra::ReadWriteVector<Number> &rw_vector)
      {
        distributed::Vector<Number, ::dealii::MemorySpace::Host> tmp_vector(
          communication_pattern);

        const unsigned int n_elements = communication_pattern->local_size();
        std::copy(values, values + n_elements, tmp_vector.begin());
        tmp_vector.update_ghost_values();

        const IndexSet &stored = rw_vector.get_stored_elements();
        if (operation == VectorOperation::add)
          for (size_type i = 0; i < stored.n_elements(); ++i)
            rw_vector.local_element(i) +=
              tmp_vector(stored.nth_index_in_set(i));
        else if (operation == VectorOperation::min)
          for (size_type i = 0; i < stored.n_elements(); ++i)
            rw_vector.local_element(i) =
              get_min(tmp_vector(stored.nth_index_in_set(i)),
                      rw_vector.local_element(i));
        else if (operation == VectorOperation::max)
          for (size_type i = 0; i < stored.n_elements(); ++i)
            rw_vector.local_element(i) =
              get_max(tmp_vector(stored.nth_index_in_set(i)),
                      rw_vector.local_element(i));
        else
          for (size_type i = 0; i < stored.n_elements(); ++i)
            rw_vector.local_element(i) = tmp_vector(stored.nth_index_in_set(i));
      }
    };



#ifdef DEAL_II_COMPILER_CUDA_AWARE
    template <typename Number>
    struct read_write_vector_functions<Number, ::dealii::MemorySpace::CUDA>
    {
      using size_type = types::global_dof_index;

      static void
      import(const std::shared_ptr<const ::dealii::Utilities::MPI::Partitioner>
               &                                     communication_pattern,
             const Number *                          values,
             const ::dealii::VectorOperation::values operation,
             ::dealii::LinearAlgebra::ReadWriteVector<Number> &rw_vector)
      {
        distributed::Vector<Number, ::dealii::MemorySpace::Host> tmp_vector(
          communication_pattern);

        const unsigned int n_elements = communication_pattern->local_size();
        cudaError_t        cuda_error_code = cudaMemcpy(tmp_vector.begin(),
                                                 values,
                                                 n_elements * sizeof(Number),
                                                 cudaMemcpyDeviceToHost);
        AssertCuda(cuda_error_code);
        tmp_vector.update_ghost_values();

        const IndexSet &stored = rw_vector.get_stored_elements();
        if (operation == VectorOperation::add)
          for (size_type i = 0; i < stored.n_elements(); ++i)
            rw_vector.local_element(i) +=
              tmp_vector(stored.nth_index_in_set(i));
        else if (operation == VectorOperation::min)
          for (size_type i = 0; i < stored.n_elements(); ++i)
            rw_vector.local_element(i) =
              get_min(tmp_vector(stored.nth_index_in_set(i)),
                      rw_vector.local_element(i));
        else if (operation == VectorOperation::max)
          for (size_type i = 0; i < stored.n_elements(); ++i)
            rw_vector.local_element(i) =
              get_max(tmp_vector(stored.nth_index_in_set(i)),
                      rw_vector.local_element(i));
        else
          for (size_type i = 0; i < stored.n_elements(); ++i)
            rw_vector.local_element(i) = tmp_vector(stored.nth_index_in_set(i));
      }
    };
#endif
  } // namespace internal


  template <typename Number>
  void
  ReadWriteVector<Number>::resize_val(const size_type new_alloc_size)
  {
    if (new_alloc_size == 0)
      {
        values.reset();
        thread_loop_partitioner =
          std::make_shared<parallel::internal::TBBPartitioner>();
      }
    else
      {
        Number *new_values;
        Utilities::System::posix_memalign(reinterpret_cast<void **>(
                                            &new_values),
                                          64,
                                          sizeof(Number) * new_alloc_size);
        values.reset(new_values);

        if (new_alloc_size >= 4 * dealii::internal::VectorImplementation::
                                    minimum_parallel_grain_size)
          thread_loop_partitioner =
            std::make_shared<parallel::internal::TBBPartitioner>();
      }
  }



  template <typename Number>
  void
  ReadWriteVector<Number>::reinit(const size_type size,
                                  const bool      omit_zeroing_entries)
  {
    // check whether we need to reallocate
    resize_val(size);

    stored_elements = complete_index_set(size);
    stored_elements.compress();

    // set entries to zero if so requested
    if (omit_zeroing_entries == false)
      this->operator=(Number());

    // reset the communication patter
    source_stored_elements.clear();
    comm_pattern.reset();
  }



  template <typename Number>
  template <typename Number2>
  void
  ReadWriteVector<Number>::reinit(const ReadWriteVector<Number2> &v,
                                  const bool omit_zeroing_entries)
  {
    resize_val(v.n_elements());

    stored_elements = v.get_stored_elements();

    if (omit_zeroing_entries == false)
      this->operator=(Number());

    // reset the communication patter
    source_stored_elements.clear();
    comm_pattern.reset();
  }



  template <typename Number>
  void
  ReadWriteVector<Number>::reinit(const IndexSet &locally_stored_indices,
                                  const bool      omit_zeroing_entries)
  {
    stored_elements = locally_stored_indices;

    // set vector size and allocate memory
    resize_val(stored_elements.n_elements());

    // initialize to zero
    if (omit_zeroing_entries == false)
      this->operator=(Number());

    // reset the communication pattern
    source_stored_elements.clear();
    comm_pattern.reset();
  }



#if defined(DEAL_II_WITH_TRILINOS) && defined(DEAL_II_WITH_MPI)
  template <typename Number>
  void
  ReadWriteVector<Number>::reinit(
    const TrilinosWrappers::MPI::Vector &trilinos_vec)
  {
    // TODO: We could avoid copying the data by just using a view into the
    // trilinos data but only if Number=double. Also update documentation that
    // the argument's lifetime needs to be longer then. If we do this, we need
    // to think about whether the view should be read/write.

    stored_elements = IndexSet(trilinos_vec.trilinos_partitioner());

    resize_val(stored_elements.n_elements());

    TrilinosScalar *start_ptr;
    int             leading_dimension;
    int ierr = trilinos_vec.trilinos_vector().ExtractView(&start_ptr,
                                                          &leading_dimension);
    AssertThrow(ierr == 0, ExcTrilinosError(ierr));

    std::copy(start_ptr, start_ptr + leading_dimension, values.get());

    // reset the communication pattern
    source_stored_elements.clear();
    comm_pattern.reset();
  }
#endif



  template <typename Number>
  template <typename Functor>
  void
  ReadWriteVector<Number>::apply(const Functor &func)
  {
    FunctorTemplate<Functor> functor(*this, func);
    dealii::internal::VectorOperations::parallel_for(functor,
                                                     0,
                                                     n_elements(),
                                                     thread_loop_partitioner);
  }



  template <typename Number>
  ReadWriteVector<Number> &
  ReadWriteVector<Number>::operator=(const ReadWriteVector<Number> &in_vector)
  {
    if (PointerComparison::equal(this, &in_vector))
      return *this;

    thread_loop_partitioner = in_vector.thread_loop_partitioner;
    if (n_elements() != in_vector.n_elements())
      reinit(in_vector, true);

    dealii::internal::VectorOperations::Vector_copy<Number, Number> copier(
      in_vector.values.get(), values.get());
    dealii::internal::VectorOperations::parallel_for(copier,
                                                     0,
                                                     n_elements(),
                                                     thread_loop_partitioner);

    return *this;
  }



  template <typename Number>
  template <typename Number2>
  ReadWriteVector<Number> &
  ReadWriteVector<Number>::operator=(const ReadWriteVector<Number2> &in_vector)
  {
    thread_loop_partitioner = in_vector.thread_loop_partitioner;
    if (n_elements() != in_vector.n_elements())
      reinit(in_vector, true);

    dealii::internal::VectorOperations::Vector_copy<Number, Number2> copier(
      in_vector.values.get(), values.get());
    dealii::internal::VectorOperations::parallel_for(copier,
                                                     0,
                                                     n_elements(),
                                                     thread_loop_partitioner);

    return *this;
  }



  template <typename Number>
  ReadWriteVector<Number> &
  ReadWriteVector<Number>::operator=(const Number s)
  {
    Assert(s == static_cast<Number>(0),
           ExcMessage("Only 0 can be assigned to a vector."));
    (void)s;

    const size_type this_size = n_elements();
    if (this_size > 0)
      {
        dealii::internal::VectorOperations::Vector_set<Number> setter(
          Number(), values.get());
        dealii::internal::VectorOperations::parallel_for(
          setter, 0, this_size, thread_loop_partitioner);
      }

    return *this;
  }



  template <typename Number>
  template <typename MemorySpace>
  void
  ReadWriteVector<Number>::import(
    const distributed::Vector<Number, MemorySpace> &vec,
    VectorOperation::values                         operation,
    const std::shared_ptr<const CommunicationPatternBase>
      &communication_pattern)
  {
    // If no communication pattern is given, create one. Otherwise, use the
    // given one.
    std::shared_ptr<const Utilities::MPI::Partitioner> comm_pattern;
    if (communication_pattern.get() == nullptr)
      {
        comm_pattern = std::make_shared<Utilities::MPI::Partitioner>(
          vec.locally_owned_elements(),
          get_stored_elements(),
          vec.get_mpi_communicator());
      }
    else
      {
        comm_pattern =
          std::dynamic_pointer_cast<const Utilities::MPI::Partitioner>(
            communication_pattern);
        AssertThrow(comm_pattern != nullptr,
                    ExcMessage("The communication pattern is not of type "
                               "Utilities::MPI::Partitioner."));
      }


    internal::read_write_vector_functions<Number, MemorySpace>::import(
      comm_pattern, vec.begin(), operation, *this);
  }



#ifdef DEAL_II_WITH_PETSC
  namespace internal
  {
    template <typename PETSC_Number, typename Number>
    void
    copy_petsc_vector(const PETSC_Number *petsc_start_ptr,
                      const PETSC_Number *petsc_end_ptr,
                      Number *            ptr)
    {
      std::copy(petsc_start_ptr, petsc_end_ptr, ptr);
    }

    template <typename PETSC_Number, typename Number>
    void
    copy_petsc_vector(const std::complex<PETSC_Number> *petsc_start_ptr,
                      const std::complex<PETSC_Number> *petsc_end_ptr,
                      std::complex<Number> *            ptr)
    {
      std::copy(petsc_start_ptr, petsc_end_ptr, ptr);
    }

    template <typename PETSC_Number, typename Number>
    void
    copy_petsc_vector(const std::complex<PETSC_Number> * /*petsc_start_ptr*/,
                      const std::complex<PETSC_Number> * /*petsc_end_ptr*/,
                      Number * /*ptr*/)
    {
      AssertThrow(false, ExcMessage("Tried to copy complex -> real"));
    }
  } // namespace internal



  template <typename Number>
  void
  ReadWriteVector<Number>::import(
    const PETScWrappers::MPI::Vector &petsc_vec,
    VectorOperation::values /*operation*/,
    const std::shared_ptr<const CommunicationPatternBase>
      & /*communication_pattern*/)
  {
    // TODO: this works only if no communication is needed.
    Assert(petsc_vec.locally_owned_elements() == stored_elements,
           StandardExceptions::ExcInvalidState());

    // get a representation of the vector and copy it
    PetscScalar *  start_ptr;
    PetscErrorCode ierr =
      VecGetArray(static_cast<const Vec &>(petsc_vec), &start_ptr);
    AssertThrow(ierr == 0, ExcPETScError(ierr));

    const size_type vec_size = petsc_vec.local_size();
    internal::copy_petsc_vector(start_ptr, start_ptr + vec_size, begin());

    // restore the representation of the vector
    ierr = VecRestoreArray(static_cast<const Vec &>(petsc_vec), &start_ptr);
    AssertThrow(ierr == 0, ExcPETScError(ierr));
  }
#endif



#if defined(DEAL_II_WITH_TRILINOS) && defined(DEAL_II_WITH_MPI)
#  ifdef DEAL_II_TRILINOS_WITH_TPETRA
  template <typename Number>
  void
  ReadWriteVector<Number>::import(
    const Tpetra::Vector<Number, int, types::global_dof_index> &vector,
    const IndexSet &                                            source_elements,
    VectorOperation::values                                     operation,
    const MPI_Comm &                                            mpi_comm,
    const std::shared_ptr<const CommunicationPatternBase>
      &communication_pattern)
  {
    std::shared_ptr<const TpetraWrappers::CommunicationPattern>
      tpetra_comm_pattern;

    // If no communication pattern is given, create one. Otherwise, use the one
    // given.
    if (communication_pattern == nullptr)
      {
        // The first time import is called, we create a communication pattern.
        // Check if the communication pattern already exists and if it can be
        // reused.
        if ((source_elements.size() == source_stored_elements.size()) &&
            (source_elements == source_stored_elements))
          {
            tpetra_comm_pattern = std::dynamic_pointer_cast<
              const TpetraWrappers::CommunicationPattern>(comm_pattern);
            if (tpetra_comm_pattern == nullptr)
              tpetra_comm_pattern =
                std::make_shared<const TpetraWrappers::CommunicationPattern>(
                  create_tpetra_comm_pattern(source_elements, mpi_comm));
          }
        else
          tpetra_comm_pattern =
            std::make_shared<const TpetraWrappers::CommunicationPattern>(
              create_tpetra_comm_pattern(source_elements, mpi_comm));
      }
    else
      {
        tpetra_comm_pattern =
          std::dynamic_pointer_cast<const TpetraWrappers::CommunicationPattern>(
            communication_pattern);
        AssertThrow(tpetra_comm_pattern != nullptr,
                    ExcMessage(
                      std::string("The communication pattern is not of type ") +
                      "LinearAlgebra::TpetraWrappers::CommunicationPattern."));
      }

    Tpetra::Export<int, types::global_dof_index> tpetra_export(
      tpetra_comm_pattern->get_tpetra_export());

    Tpetra::Vector<Number, int, types::global_dof_index> target_vector(
      tpetra_export.getSourceMap());
    target_vector.doImport(vector, tpetra_export, Tpetra::REPLACE);

    const auto *new_values = target_vector.getData().get();
    const auto  size       = target_vector.getLocalLength();

    using size_type = typename std::decay<decltype(size)>::type;

    Assert(size == 0 || values != nullptr, ExcInternalError("Export failed."));
    AssertDimension(size, stored_elements.n_elements());

    switch (operation)
      {
        case VectorOperation::insert:
          for (size_type i = 0; i < size; ++i)
            values[i] = new_values[i];
          break;

        case VectorOperation::add:
          for (size_type i = 0; i < size; ++i)
            values[i] += new_values[i];
          break;

        case VectorOperation::min:
          // To ensure that this code also compiles with complex
          // numbers, we only compare the real part of the
          // variable. Note that min/max do not make sense with complex
          // numbers.
          for (size_type i = 0; i < size; ++i)
            {
              Assert(
                std::imag(new_values[i]) == 0.,
                ExcMessage(
                  "VectorOperation::min is not defined if there is an imaginary part!)"));
              Assert(
                std::imag(values[i]) == 0.,
                ExcMessage(
                  "VectorOperation::min is not defined if there is an imaginary part!)"));
              if (std::real(new_values[i]) - std::real(values[i]) < 0.0)
                values[i] = new_values[i];
            }
          break;

        case VectorOperation::max:
          for (size_type i = 0; i < size; ++i)
            {
              Assert(
                std::imag(new_values[i]) == 0.,
                ExcMessage(
                  "VectorOperation::max is not defined if there is an imaginary part!)"));
              Assert(
                std::imag(values[i]) == 0.,
                ExcMessage(
                  "VectorOperation::max is not defined if there is an imaginary part!)"));
              if (std::real(new_values[i]) - std::real(values[i]) > 0.0)
                values[i] = new_values[i];
            }
          break;

        default:
          AssertThrow(false, ExcNotImplemented());
      }
  }
#  endif



  template <typename Number>
  void
  ReadWriteVector<Number>::import(
    const Epetra_MultiVector &multivector,
    const IndexSet &          source_elements,
    VectorOperation::values   operation,
    const MPI_Comm &          mpi_comm,
    const std::shared_ptr<const CommunicationPatternBase>
      &communication_pattern)
  {
    std::shared_ptr<const EpetraWrappers::CommunicationPattern>
      epetra_comm_pattern;

    // If no communication pattern is given, create one. Otherwise, use the one
    // given.
    if (communication_pattern == nullptr)
      {
        // The first time import is called, we create a communication pattern.
        // Check if the communication pattern already exists and if it can be
        // reused.
        if ((source_elements.size() == source_stored_elements.size()) &&
            (source_elements == source_stored_elements))
          {
            epetra_comm_pattern = std::dynamic_pointer_cast<
              const EpetraWrappers::CommunicationPattern>(comm_pattern);
            if (epetra_comm_pattern == nullptr)
              epetra_comm_pattern =
                std::make_shared<const EpetraWrappers::CommunicationPattern>(
                  create_epetra_comm_pattern(source_elements, mpi_comm));
          }
        else
          epetra_comm_pattern =
            std::make_shared<const EpetraWrappers::CommunicationPattern>(
              create_epetra_comm_pattern(source_elements, mpi_comm));
      }
    else
      {
        epetra_comm_pattern =
          std::dynamic_pointer_cast<const EpetraWrappers::CommunicationPattern>(
            communication_pattern);
        AssertThrow(epetra_comm_pattern != nullptr,
                    ExcMessage(
                      std::string("The communication pattern is not of type ") +
                      "LinearAlgebra::EpetraWrappers::CommunicationPattern."));
      }

    Epetra_Import import(epetra_comm_pattern->get_epetra_import());

    Epetra_FEVector target_vector(import.TargetMap());

    if (operation == VectorOperation::insert)
      {
        const int err = target_vector.Import(multivector, import, Insert);
        AssertThrow(err == 0,
                    ExcMessage("Epetra Import() failed with error code: " +
                               std::to_string(err)));

        const double *new_values = target_vector.Values();
        const int     size       = target_vector.MyLength();
        Assert(size == 0 || values != nullptr,
               ExcInternalError("Import failed."));

        for (int i = 0; i < size; ++i)
          values[i] = new_values[i];
      }
    else if (operation == VectorOperation::add)
      {
        const int err = target_vector.Import(multivector, import, Add);
        AssertThrow(err == 0,
                    ExcMessage("Epetra Import() failed with error code: " +
                               std::to_string(err)));

        const double *new_values = target_vector.Values();
        const int     size       = target_vector.MyLength();
        Assert(size == 0 || values != nullptr,
               ExcInternalError("Import failed."));

        for (int i = 0; i < size; ++i)
          values[i] += new_values[i];
      }
    else if (operation == VectorOperation::min)
      {
        const int err = target_vector.Import(multivector, import, Add);
        AssertThrow(err == 0,
                    ExcMessage("Epetra Import() failed with error code: " +
                               std::to_string(err)));

        const double *new_values = target_vector.Values();
        const int     size       = target_vector.MyLength();
        Assert(size == 0 || values != nullptr,
               ExcInternalError("Import failed."));

        // To ensure that this code also compiles with complex
        // numbers, we only compare the real part of the
        // variable. Note that min/max do not make sense with complex
        // numbers.
        for (int i = 0; i < size; ++i)
          {
            Assert(
              std::imag(new_values[i]) == 0.,
              ExcMessage(
                "VectorOperation::min is not defined if there is an imaginary part!)"));
            Assert(
              std::imag(values[i]) == 0.,
              ExcMessage(
                "VectorOperation::min is not defined if there is an imaginary part!)"));
            if (std::real(new_values[i]) - std::real(values[i]) < 0.0)
              values[i] = new_values[i];
          }
      }
    else if (operation == VectorOperation::max)
      {
        const int err = target_vector.Import(multivector, import, Add);
        AssertThrow(err == 0,
                    ExcMessage("Epetra Import() failed with error code: " +
                               std::to_string(err)));

        const double *new_values = target_vector.Values();
        const int     size       = target_vector.MyLength();
        Assert(size == 0 || values != nullptr,
               ExcInternalError("Import failed."));

        for (int i = 0; i < size; ++i)
          {
            Assert(
              std::imag(new_values[i]) == 0.,
              ExcMessage(
                "VectorOperation::max is not defined if there is an imaginary part!)"));
            Assert(
              std::imag(values[i]) == 0.,
              ExcMessage(
                "VectorOperation::max is not defined if there is an imaginary part!)"));
            if (std::real(new_values[i]) - std::real(values[i]) > 0.0)
              values[i] = new_values[i];
          }
      }
    else
      AssertThrow(false, ExcNotImplemented());
  }


  template <typename Number>
  void
  ReadWriteVector<Number>::import(
    const TrilinosWrappers::MPI::Vector &trilinos_vec,
    VectorOperation::values              operation,
    const std::shared_ptr<const CommunicationPatternBase>
      &communication_pattern)
  {
    // While the import does work with Trilinos 12.8.x, it fails with 12.4.x. To
    // be safe, we disable it here. Note that it would be a useful case, as
    // ReadWriteVector is supposed to replace ghosted vectors anyways.
    AssertThrow(
      !trilinos_vec.has_ghost_elements(),
      ExcMessage(
        "Import() from TrilinosWrappers::MPI::Vector with ghost entries is not supported!"));
    import(trilinos_vec.trilinos_vector(),
           trilinos_vec.locally_owned_elements(),
           operation,
           trilinos_vec.get_mpi_communicator(),
           communication_pattern);
  }



#  ifdef DEAL_II_TRILINOS_WITH_TPETRA
  template <typename Number>
  void
  ReadWriteVector<Number>::import(
    const LinearAlgebra::TpetraWrappers::Vector<Number> &trilinos_vec,
    VectorOperation::values                              operation,
    const std::shared_ptr<const CommunicationPatternBase>
      &communication_pattern)
  {
    import(trilinos_vec.trilinos_vector(),
           trilinos_vec.locally_owned_elements(),
           operation,
           trilinos_vec.get_mpi_communicator(),
           communication_pattern);
  }
#  endif



  template <typename Number>
  void
  ReadWriteVector<Number>::import(
    const LinearAlgebra::EpetraWrappers::Vector &trilinos_vec,
    VectorOperation::values                      operation,
    const std::shared_ptr<const CommunicationPatternBase>
      &communication_pattern)
  {
    import(trilinos_vec.trilinos_vector(),
           trilinos_vec.locally_owned_elements(),
           operation,
           trilinos_vec.get_mpi_communicator(),
           communication_pattern);
  }
#endif



#ifdef DEAL_II_COMPILER_CUDA_AWARE
  template <typename Number>
  void
  ReadWriteVector<Number>::import(
    const LinearAlgebra::CUDAWrappers::Vector<Number> &cuda_vec,
    VectorOperation::values                            operation,
    const std::shared_ptr<const CommunicationPatternBase> &)
  {
    const unsigned int n_elements = stored_elements.n_elements();
    if (operation == VectorOperation::insert)
      {
        cudaError_t error_code = cudaMemcpy(values.get(),
                                            cuda_vec.get_values(),
                                            n_elements * sizeof(Number),
                                            cudaMemcpyDeviceToHost);
        AssertCuda(error_code);
      }
    else if (operation == VectorOperation::add)
      {
        // Copy the vector from the device to a temporary vector on the host
        std::vector<Number> tmp(n_elements);
        cudaError_t         error_code = cudaMemcpy(tmp.data(),
                                            cuda_vec.get_values(),
                                            n_elements * sizeof(Number),
                                            cudaMemcpyDeviceToHost);
        AssertCuda(error_code);

        // Add the two vectors
        for (unsigned int i = 0; i < n_elements; ++i)
          values[i] += tmp[i];
      }
    else if (operation == VectorOperation::min)
      {
        // Copy the vector from the device to a temporary vector on the host
        std::vector<Number> tmp(n_elements);
        cudaError_t         error_code = cudaMemcpy(tmp.data(),
                                            cuda_vec.get_values(),
                                            n_elements * sizeof(Number),
                                            cudaMemcpyDeviceToHost);
        AssertCuda(error_code);

        // To ensure that this code also compiles with complex
        // numbers, we only compare the real part of the
        // variable. Note that min/max do not make sense with complex
        // numbers.
        for (unsigned int i = 0; i < n_elements; ++i)
          {
            Assert(
              std::imag(tmp[i]) == 0.,
              ExcMessage(
                "VectorOperation::min is not defined if there is an imaginary part!)"));
            Assert(
              std::imag(values[i]) == 0.,
              ExcMessage(
                "VectorOperation::min is not defined if there is an imaginary part!)"));
            if (std::real(tmp[i]) - std::real(values[i]) < 0.0)
              values[i] = tmp[i];
          }
      }
    else if (operation == VectorOperation::max)
      {
        // Copy the vector from the device to a temporary vector on the host
        std::vector<Number> tmp(n_elements);
        cudaError_t         error_code = cudaMemcpy(tmp.data(),
                                            cuda_vec.get_values(),
                                            n_elements * sizeof(Number),
                                            cudaMemcpyDeviceToHost);
        AssertCuda(error_code);

        for (unsigned int i = 0; i < n_elements; ++i)
          {
            Assert(
              std::imag(tmp[i]) == 0.,
              ExcMessage(
                "VectorOperation::max is not defined if there is an imaginary part!)"));
            Assert(
              std::imag(values[i]) == 0.,
              ExcMessage(
                "VectorOperation::max is not defined if there is an imaginary part!)"));
            if (std::real(tmp[i]) - std::real(values[i]) > 0.0)
              values[i] = tmp[i];
          }
      }
    else
      AssertThrow(false, ExcNotImplemented());
  }
#endif



  template <typename Number>
  void
  ReadWriteVector<Number>::swap(ReadWriteVector<Number> &v)
  {
    std::swap(stored_elements, v.stored_elements);
    std::swap(values, v.values);
  }



  template <typename Number>
  std::size_t
  ReadWriteVector<Number>::memory_consumption() const
  {
    std::size_t memory = sizeof(*this);
    memory += sizeof(Number) * static_cast<std::size_t>(this->n_elements());

    memory += stored_elements.memory_consumption();

    return memory;
  }



  template <typename Number>
  void
  ReadWriteVector<Number>::print(std::ostream &     out,
                                 const unsigned int precision,
                                 const bool         scientific) const
  {
    AssertThrow(out, ExcIO());
    boost::io::ios_flags_saver restore_flags(out);

    out.precision(precision);
    if (scientific)
      out.setf(std::ios::scientific, std::ios::floatfield);
    else
      out.setf(std::ios::fixed, std::ios::floatfield);

    out << "IndexSet: ";
    stored_elements.print(out);
    out << std::endl;
    unsigned int i = 0;
    for (const auto idx : this->stored_elements)
      out << "[" << idx << "]: " << values[i++] << '\n';
    out << std::flush;

    AssertThrow(out, ExcIO());
  }



#if defined(DEAL_II_WITH_TRILINOS) && defined(DEAL_II_WITH_MPI)
#  ifdef DEAL_II_TRILINOS_WITH_TPETRA
  template <typename Number>
  TpetraWrappers::CommunicationPattern
  ReadWriteVector<Number>::create_tpetra_comm_pattern(
    const IndexSet &source_index_set,
    const MPI_Comm &mpi_comm)
  {
    source_stored_elements = source_index_set;
    TpetraWrappers::CommunicationPattern epetra_comm_pattern(
      source_stored_elements, stored_elements, mpi_comm);
    comm_pattern = std::make_shared<TpetraWrappers::CommunicationPattern>(
      source_stored_elements, stored_elements, mpi_comm);

    return epetra_comm_pattern;
  }
#  endif



  template <typename Number>
  EpetraWrappers::CommunicationPattern
  ReadWriteVector<Number>::create_epetra_comm_pattern(
    const IndexSet &source_index_set,
    const MPI_Comm &mpi_comm)
  {
    source_stored_elements = source_index_set;
    EpetraWrappers::CommunicationPattern epetra_comm_pattern(
      source_stored_elements, stored_elements, mpi_comm);
    comm_pattern = std::make_shared<EpetraWrappers::CommunicationPattern>(
      source_stored_elements, stored_elements, mpi_comm);

    return epetra_comm_pattern;
  }
#endif
} // end of namespace LinearAlgebra


DEAL_II_NAMESPACE_CLOSE

#endif
