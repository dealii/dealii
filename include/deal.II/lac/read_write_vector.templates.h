// ---------------------------------------------------------------------
//
// Copyright (C) 2015-2016 by the deal.II authors
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

#ifndef dealii__parallel_vector_templates_h
#define dealii__parallel_vector_templates_h


#include <deal.II/base/config.h>
#include <deal.II/lac/read_write_vector.h>
#include <deal.II/lac/petsc_parallel_vector.h>
#include <deal.II/lac/trilinos_vector.h>
#include <deal.II/lac/trilinos_epetra_vector.h>

#ifdef DEAL_II_WITH_TRILINOS
#  include <deal.II/lac/trilinos_epetra_communication_pattern.h>
#  include "Epetra_Import.h"
#endif

DEAL_II_NAMESPACE_OPEN


namespace LinearAlgebra
{
  template <typename Number>
  void
  ReadWriteVector<Number>::resize_val (const size_type new_alloc_size)
  {
    if (new_alloc_size == 0)
      {
        if (val != NULL)
          free(val);
        val = NULL;
      }
    else
      {
        if (val != NULL)
          free(val);

        Utilities::System::posix_memalign ((void **)&val, 64, sizeof(Number)*new_alloc_size);
      }
  }



  template <typename Number>
  void
  ReadWriteVector<Number>::reinit (const size_type size,
                                   const bool      omit_zeroing_entries)
  {
    // check whether we need to reallocate
    resize_val(size);

    stored_elements = complete_index_set(size);
    stored_elements.compress();

    // set entries to zero if so requested
    if (omit_zeroing_entries == false)
      this->operator = (Number());

    // reset the communication patter
    source_stored_elements.clear();
    comm_pattern.reset();
  }



  template <typename Number>
  template <typename Number2>
  void
  ReadWriteVector<Number>::reinit (const ReadWriteVector<Number2> &v,
                                   const bool             omit_zeroing_entries)
  {
    resize_val(v.n_elements());

    stored_elements = v.get_stored_elements();

    if (omit_zeroing_entries == false)
      this->operator= (Number());

    // reset the communication patter
    source_stored_elements.clear();
    comm_pattern.reset();
  }



  template <typename Number>
  void
  ReadWriteVector<Number>::reinit (const IndexSet &locally_stored_indices,
                                   const bool      omit_zeroing_entries)
  {
    stored_elements = locally_stored_indices;

    // set vector size and allocate memory
    resize_val(stored_elements.n_elements());

    // initialize to zero
    if (omit_zeroing_entries == false)
      this->operator= (Number());

    // reset the communication patter
    source_stored_elements.clear();
    comm_pattern.reset();
  }



#ifdef DEAL_II_WITH_PETSC
  namespace internal
  {
    template <typename PETSC_Number, typename Number>
    void copy_petsc_vector (const PETSC_Number *petsc_start_ptr,
                            const PETSC_Number *petsc_end_ptr,
                            Number *ptr)
    {
      std::copy(petsc_start_ptr, petsc_end_ptr, ptr);
    }

    template <typename PETSC_Number, typename Number>
    void copy_petsc_vector (const std::complex<PETSC_Number> *petsc_start_ptr,
                            const std::complex<PETSC_Number> *petsc_end_ptr,
                            std::complex<Number> *ptr)
    {
      std::copy(petsc_start_ptr, petsc_end_ptr, ptr);
    }

    template <typename PETSC_Number, typename Number>
    void copy_petsc_vector (const std::complex<PETSC_Number> * /*petsc_start_ptr*/,
                            const std::complex<PETSC_Number> * /*petsc_end_ptr*/,
                            Number * /*ptr*/)
    {
      AssertThrow(false, ExcMessage("Tried to copy complex -> real"));
    }
  }



  template <typename Number>
  void
  ReadWriteVector<Number>::import(const PETScWrappers::MPI::Vector &petsc_vec,
                                  VectorOperation::values /*operation*/,
                                  std_cxx11::shared_ptr<const CommunicationPatternBase> /*communication_pattern*/)
  {
    //TODO: this works only if no communication is needed.
    Assert(petsc_vec.locally_owned_elements() == stored_elements,
           StandardExceptions::ExcInvalidState());

    // get a representation of the vector and copy it
    PetscScalar *start_ptr;
    int ierr = VecGetArray (static_cast<const Vec &>(petsc_vec), &start_ptr);
    AssertThrow (ierr == 0, ExcPETScError(ierr));

    const size_type vec_size = petsc_vec.local_size();
    internal::copy_petsc_vector (start_ptr, start_ptr + vec_size, begin());

    // restore the representation of the vector
    ierr = VecRestoreArray (static_cast<const Vec &>(petsc_vec), &start_ptr);
    AssertThrow (ierr == 0, ExcPETScError(ierr));
  }
#endif



#if defined(DEAL_II_WITH_TRILINOS) && defined(DEAL_II_WITH_MPI)
  template <typename Number>
  void
  ReadWriteVector<Number>::import(const Epetra_MultiVector        &multivector,
                                  const IndexSet                  &source_elements,
                                  VectorOperation::values          operation,
                                  const MPI_Comm                  &mpi_comm,
                                  std_cxx11::shared_ptr<const CommunicationPatternBase> communication_pattern)
  {
    std_cxx11::shared_ptr<const EpetraWrappers::CommunicationPattern> epetra_comm_pattern;

    // If no communication pattern is given, create one. Otherwise, use the one
    // given.
    if (communication_pattern == NULL)
      {
        // The first time import is called, we create a communication pattern.
        // Check if the communication pattern already exists and if it can be
        // reused.
        if ((source_elements.size() == source_stored_elements.size()) &&
            (source_elements == source_stored_elements))
          {
            epetra_comm_pattern =
              std_cxx11::dynamic_pointer_cast<const EpetraWrappers::CommunicationPattern> (comm_pattern);
            if (epetra_comm_pattern == NULL)
              epetra_comm_pattern = std_cxx11::make_shared<const EpetraWrappers::CommunicationPattern>(
                                      create_epetra_comm_pattern(source_elements, mpi_comm));
          }
        else
          epetra_comm_pattern = std_cxx11::make_shared<const EpetraWrappers::CommunicationPattern>(
                                  create_epetra_comm_pattern(source_elements, mpi_comm));
      }
    else
      {
        epetra_comm_pattern = std_cxx11::dynamic_pointer_cast<const EpetraWrappers::CommunicationPattern> (
                                communication_pattern);
        AssertThrow(epetra_comm_pattern != NULL,
                    ExcMessage(std::string("The communication pattern is not of type ") +
                               "LinearAlgebra::EpetraWrappers::CommunicationPattern."));
      }

    Epetra_Import import(epetra_comm_pattern->get_epetra_import());

    Epetra_FEVector target_vector(import.TargetMap());

    if (operation==VectorOperation::insert)
      {
        target_vector.Import(multivector, import, Insert);
        double *values = target_vector.Values();
        for (int i=0; i<target_vector.MyLength(); ++i)
          val[i] = values[i];
      }
    else
      {
        target_vector.Import(multivector, import, Add);
        double *values = target_vector.Values();
        for (int i=0; i<target_vector.MyLength(); ++i)
          val[i] += values[i];
      }
  }


  template <typename Number>
  void
  ReadWriteVector<Number>::import(const TrilinosWrappers::MPI::Vector            &trilinos_vec,
                                  VectorOperation::values                         operation,
                                  std_cxx11::shared_ptr<const CommunicationPatternBase> communication_pattern)
  {
    import(trilinos_vec.trilinos_vector(), trilinos_vec.locally_owned_elements(),
           operation, trilinos_vec.get_mpi_communicator(), communication_pattern);
  }



  template <typename Number>
  void
  ReadWriteVector<Number>::import(const LinearAlgebra::EpetraWrappers::Vector    &trilinos_vec,
                                  VectorOperation::values                         operation,
                                  std_cxx11::shared_ptr<const CommunicationPatternBase> communication_pattern)
  {
    import(trilinos_vec.trilinos_vector(), trilinos_vec.locally_owned_elements(),
           operation, trilinos_vec.get_mpi_communicator(), communication_pattern);
  }
#endif



  template <typename Number>
  void
  ReadWriteVector<Number>::swap (ReadWriteVector<Number> &v)
  {
    std::swap(stored_elements, v.stored_elements);
    std::swap(val, v.val);
  }



  template <typename Number>
  std::size_t
  ReadWriteVector<Number>::memory_consumption () const
  {
    std::size_t memory = sizeof(*this);
    memory += sizeof (Number) * static_cast<std::size_t>(this->n_elements());

    memory += stored_elements.memory_consumption();

    return memory;
  }



  template <typename Number>
  void
  ReadWriteVector<Number>::print (std::ostream      &out,
                                  const unsigned int precision,
                                  const bool         scientific) const
  {
    AssertThrow (out, ExcIO());
    std::ios::fmtflags old_flags = out.flags();
    unsigned int old_precision = out.precision (precision);

    out.precision (precision);
    if (scientific)
      out.setf (std::ios::scientific, std::ios::floatfield);
    else
      out.setf (std::ios::fixed, std::ios::floatfield);

    out << "IndexSet: ";
    stored_elements.print(out);
    out << std::endl;
    for (unsigned int i=0; i<this->n_elements(); ++i)
      out << val[i] << std::endl;
    out << std::flush;


    AssertThrow (out, ExcIO());
    // reset output format
    out.flags (old_flags);
    out.precision(old_precision);
  }



#if defined(DEAL_II_WITH_TRILINOS) && defined(DEAL_II_WITH_MPI)
  template <typename Number>
  EpetraWrappers::CommunicationPattern
  ReadWriteVector<Number>::create_epetra_comm_pattern(const IndexSet &source_index_set,
                                                      const MPI_Comm &mpi_comm)
  {
    source_stored_elements = source_index_set;
    EpetraWrappers::CommunicationPattern epetra_comm_pattern(
      source_stored_elements, stored_elements, mpi_comm);
    comm_pattern.reset(new EpetraWrappers::CommunicationPattern(
                         source_stored_elements, stored_elements, mpi_comm));

    return epetra_comm_pattern;
  }
#endif
} // end of namespace LinearAlgebra


DEAL_II_NAMESPACE_CLOSE

#endif
