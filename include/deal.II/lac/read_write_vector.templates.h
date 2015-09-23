// ---------------------------------------------------------------------
//
// Copyright (C) 2015 by the deal.II authors
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

#ifdef DEAL_II_WITH_TRILINOS
#  include <Epetra_Import.h>
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
                                   const bool      fast)
  {
    // check whether we need to reallocate
    resize_val(size);

    stored_elements = complete_index_set(size);

    // set entries to zero if so requested
    if (fast == false)
      this->operator = (Number());
  }



  template <typename Number>
  template <typename Number2>
  void
  ReadWriteVector<Number>::reinit (const ReadWriteVector<Number2> &v,
                                   const bool             fast)
  {
    resize_val(v.n_elements());

    stored_elements = v.get_stored_elements();

    if (fast == false)
      this->operator= (Number());
  }



  template <typename Number>
  void
  ReadWriteVector<Number>::reinit (const IndexSet &locally_stored_indices,
                                   const bool      fast)
  {
    stored_elements = locally_stored_indices;

    // set vector size and allocate memory
    resize_val(stored_elements.n_elements());

    // initialize to zero
    if (fast == false)
      this->operator= (Number());
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
    void copy_petsc_vector (const std::complex<PETSC_Number> *petsc_start_ptr,
                            const std::complex<PETSC_Number> *petsc_end_ptr,
                            Number *ptr)
    {
      AssertThrow(false, ExcMessage("Tried to copy complex -> real"));
    }
  }



  template <typename Number>
  ReadWriteVector<Number> &
  ReadWriteVector<Number>::operator = (const PETScWrappers::MPI::Vector &petsc_vec)
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

    // return a pointer to this object per normal c++ operator overloading
    // semantics
    return *this;
  }
#endif



#ifdef DEAL_II_WITH_TRILINOS
  template <typename Number>
  ReadWriteVector<Number> &
  ReadWriteVector<Number>::operator = (const TrilinosWrappers::MPI::Vector &trilinos_vec)
  {
    // Copy the local elements
    std::vector<size_type> trilinos_indices, readwrite_indices;
    trilinos_vec.locally_owned_elements().fill_index_vector(trilinos_indices);
    stored_elements.fill_index_vector(readwrite_indices);
    std::vector<size_type> intersection(trilinos_indices.size());
    std::vector<size_type>::iterator end;
    end = std::set_intersection(trilinos_indices.begin(), trilinos_indices.end(),
                                readwrite_indices.begin(), readwrite_indices.end(),
                                intersection.begin());
    const size_t intersection_size = end-intersection.begin();
    intersection.resize(intersection_size);
    std::vector<size_t> trilinos_position(intersection_size);
    unsigned int j=0;
    for (size_t i=0; i<intersection_size; ++i)
      {
        while (intersection[i]!=trilinos_indices[j])
          ++j;
        trilinos_position[i] = j;
      }
    double *values = trilinos_vec.trilinos_vector().Values();
    for (size_t i=0; i<trilinos_position.size(); ++i)
      val[global_to_local(intersection[i])] = values[trilinos_position[i]];


#ifdef DEAL_II_WITH_MPI
    // Copy the off-processor elements if necessary
    if (intersection_size != n_elements())
      {
        Epetra_BlockMap source_map = trilinos_vec.trilinos_vector().Map();
        std::vector<size_type> off_proc_indices(readwrite_indices.size());
        // Subtract the elements of readwrite that are in intersection.
        end = std::set_difference(readwrite_indices.begin(), readwrite_indices.end(),
                                  intersection.begin(), intersection.end(), off_proc_indices.begin());
        off_proc_indices.resize(end-off_proc_indices.begin());
        // Cast size_type to TrilinosWrappers::type::int_type
        TrilinosWrappers::types::int_type *off_proc_array =
          new TrilinosWrappers::types::int_type [off_proc_indices.size()];
        for (size_t i=0; i<off_proc_indices.size(); ++i)
          off_proc_array[i] = static_cast<TrilinosWrappers::types::int_type> (
                                off_proc_indices[i]);
        Epetra_Map target_map(off_proc_indices.size(),off_proc_indices.size(),
                              off_proc_array,0,source_map.Comm());
        delete [] off_proc_array;
        Epetra_Import import(target_map, source_map);
        Epetra_FEVector target_vector(target_map);
        target_vector.Import(trilinos_vec.trilinos_vector(), import, Insert);
        values = target_vector.Values();
        for (unsigned int i=0; i<off_proc_indices.size(); ++i)
          val[global_to_local(off_proc_indices[i])] = values[i];
      }
#endif

    // return a pointer to this object per normal c++ operator overloading
    // semantics
    return *this;
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
} // end of namespace LinearAlgebra


DEAL_II_NAMESPACE_CLOSE

#endif
