// ---------------------------------------------------------------------
//
// Copyright (C) 1999 - 2018 by the deal.II authors
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

#ifndef dealii_parallel_block_vector_templates_h
#define dealii_parallel_block_vector_templates_h


#include <deal.II/base/config.h>

#include <deal.II/lac/la_parallel_block_vector.h>
#include <deal.II/lac/lapack_support.h>
#include <deal.II/lac/petsc_block_vector.h>
#include <deal.II/lac/read_write_vector.h>
#include <deal.II/lac/trilinos_parallel_block_vector.h>
#include <deal.II/lac/vector.h>


DEAL_II_NAMESPACE_OPEN

// Forward type declaration to have special treatment of
// LAPACKFullMatrix<number>
// in multivector_inner_product()
#ifndef DOXYGEN
template <typename Number>
class LAPACKFullMatrix;
#endif

namespace LinearAlgebra
{
  namespace distributed
  {
    template <typename Number>
    BlockVector<Number>::BlockVector(const size_type n_blocks,
                                     const size_type block_size)
    {
      reinit(n_blocks, block_size);
    }



    template <typename Number>
    BlockVector<Number>::BlockVector(const std::vector<size_type> &n)
    {
      reinit(n, false);
    }


    template <typename Number>
    BlockVector<Number>::BlockVector(const std::vector<IndexSet> &local_ranges,
                                     const std::vector<IndexSet> &ghost_indices,
                                     const MPI_Comm               communicator)
    {
      std::vector<size_type> sizes(local_ranges.size());
      for (unsigned int i = 0; i < local_ranges.size(); ++i)
        sizes[i] = local_ranges[i].size();

      this->block_indices.reinit(sizes);
      this->components.resize(this->n_blocks());

      for (unsigned int i = 0; i < this->n_blocks(); ++i)
        this->block(i).reinit(local_ranges[i], ghost_indices[i], communicator);
    }


    template <typename Number>
    BlockVector<Number>::BlockVector(const std::vector<IndexSet> &local_ranges,
                                     const MPI_Comm               communicator)
    {
      std::vector<size_type> sizes(local_ranges.size());
      for (unsigned int i = 0; i < local_ranges.size(); ++i)
        sizes[i] = local_ranges[i].size();

      this->block_indices.reinit(sizes);
      this->components.resize(this->n_blocks());

      for (unsigned int i = 0; i < this->n_blocks(); ++i)
        this->block(i).reinit(local_ranges[i], communicator);
    }



    template <typename Number>
    BlockVector<Number>::BlockVector(const BlockVector<Number> &v)
      : BlockVectorBase<Vector<Number>>()
    {
      this->components.resize(v.n_blocks());
      this->block_indices = v.block_indices;

      for (size_type i = 0; i < this->n_blocks(); ++i)
        this->components[i] = v.components[i];
    }



    template <typename Number>
    template <typename OtherNumber>
    BlockVector<Number>::BlockVector(const BlockVector<OtherNumber> &v)
    {
      reinit(v, true);
      *this = v;
    }



    template <typename Number>
    void
    BlockVector<Number>::reinit(const size_type n_bl,
                                const size_type bl_sz,
                                const bool      omit_zeroing_entries)
    {
      std::vector<size_type> n(n_bl, bl_sz);
      reinit(n, omit_zeroing_entries);
    }


    template <typename Number>
    void
    BlockVector<Number>::reinit(const std::vector<size_type> &n,
                                const bool omit_zeroing_entries)
    {
      this->block_indices.reinit(n);
      if (this->components.size() != this->n_blocks())
        this->components.resize(this->n_blocks());

      for (size_type i = 0; i < this->n_blocks(); ++i)
        this->components[i].reinit(n[i], omit_zeroing_entries);
    }



    template <typename Number>
    template <typename Number2>
    void
    BlockVector<Number>::reinit(const BlockVector<Number2> &v,
                                const bool omit_zeroing_entries)
    {
      this->block_indices = v.get_block_indices();
      if (this->components.size() != this->n_blocks())
        this->components.resize(this->n_blocks());

      for (unsigned int i = 0; i < this->n_blocks(); ++i)
        this->block(i).reinit(v.block(i), omit_zeroing_entries);
    }



    template <typename Number>
    BlockVector<Number> &
    BlockVector<Number>::operator=(const value_type s)
    {
      AssertIsFinite(s);

      BaseClass::operator=(s);
      return *this;
    }



    template <typename Number>
    BlockVector<Number> &
    BlockVector<Number>::operator=(const BlockVector &v)
    {
      // we only allow assignment to vectors with the same number of blocks
      // or to an empty BlockVector
      Assert(this->n_blocks() == 0 || this->n_blocks() == v.n_blocks(),
             ExcDimensionMismatch(this->n_blocks(), v.n_blocks()));

      if (this->n_blocks() != v.n_blocks())
        reinit(v.n_blocks(), true);

      for (size_type i = 0; i < this->n_blocks(); ++i)
        this->components[i] = v.block(i);

      this->collect_sizes();
      return *this;
    }



    template <typename Number>
    BlockVector<Number> &
    BlockVector<Number>::operator=(const Vector<Number> &v)
    {
      BaseClass::operator=(v);
      return *this;
    }



    template <typename Number>
    template <typename Number2>
    BlockVector<Number> &
    BlockVector<Number>::operator=(const BlockVector<Number2> &v)
    {
      reinit(v, true);
      BaseClass::operator=(v);
      return *this;
    }



#ifdef DEAL_II_WITH_PETSC

    namespace petsc_helpers
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
    } // namespace petsc_helpers



    template <typename Number>
    BlockVector<Number> &
    BlockVector<Number>::
    operator=(const PETScWrappers::MPI::BlockVector &petsc_vec)
    {
      AssertDimension(this->n_blocks(), petsc_vec.n_blocks());
      for (unsigned int i = 0; i < this->n_blocks(); ++i)
        {
          // We would like to use the same compact infrastructure as for the
          // Trilinos vector below, but the interface through ReadWriteVector
          // does not support overlapping (ghosted) PETSc vectors, which we need
          // for backward compatibility.

          Assert(petsc_vec.block(i).locally_owned_elements() ==
                   this->block(i).locally_owned_elements(),
                 StandardExceptions::ExcInvalidState());

          // get a representation of the vector and copy it
          PetscScalar *  start_ptr;
          PetscErrorCode ierr =
            VecGetArray(static_cast<const Vec &>(petsc_vec.block(i)),
                        &start_ptr);
          AssertThrow(ierr == 0, ExcPETScError(ierr));

          const size_type vec_size = this->block(i).local_size();
          petsc_helpers::copy_petsc_vector(start_ptr,
                                           start_ptr + vec_size,
                                           this->block(i).begin());

          // restore the representation of the vector
          ierr = VecRestoreArray(static_cast<const Vec &>(petsc_vec.block(i)),
                                 &start_ptr);
          AssertThrow(ierr == 0, ExcPETScError(ierr));

          // spread ghost values between processes?
          if (this->block(i).vector_is_ghosted ||
              petsc_vec.block(i).has_ghost_elements())
            this->block(i).update_ghost_values();
        }

      return *this;
    }

#endif



#ifdef DEAL_II_WITH_TRILINOS

    template <typename Number>
    BlockVector<Number> &
    BlockVector<Number>::
    operator=(const TrilinosWrappers::MPI::BlockVector &trilinos_vec)
    {
      AssertDimension(this->n_blocks(), trilinos_vec.n_blocks());
      for (unsigned int i = 0; i < this->n_blocks(); ++i)
        {
          const auto &partitioner  = this->block(i).get_partitioner();
          IndexSet    combined_set = partitioner->locally_owned_range();
          combined_set.add_indices(partitioner->ghost_indices());
          ReadWriteVector<Number> rw_vector(combined_set);
          rw_vector.import(trilinos_vec.block(i), VectorOperation::insert);
          this->block(i).import(rw_vector, VectorOperation::insert);

          if (this->block(i).has_ghost_elements() ||
              trilinos_vec.block(i).has_ghost_elements())
            this->block(i).update_ghost_values();
        }

      return *this;
    }

#endif



    template <typename Number>
    void
    BlockVector<Number>::compress(::dealii::VectorOperation::values operation)
    {
      const unsigned int n_chunks =
        (this->n_blocks() + communication_block_size - 1) /
        communication_block_size;
      for (unsigned int c = 0; c < n_chunks; ++c)
        {
          const unsigned int start = c * communication_block_size;
          const unsigned int end =
            std::min(start + communication_block_size, this->n_blocks());

          // start all requests for all blocks before finishing the transfers as
          // this saves repeated synchronizations. In order to avoid conflict
          // with possible other ongoing communication requests (from
          // LA::distributed::Vector that supports unfinished requests), add
          // 100 to the communication tag (the first 100 can be used by normal
          // vectors).
          for (unsigned int block = start; block < end; ++block)
            this->block(block).compress_start(block - start + 100, operation);
          for (unsigned int block = start; block < end; ++block)
            this->block(block).compress_finish(operation);
        }
    }



    template <typename Number>
    void
    BlockVector<Number>::update_ghost_values() const
    {
      const unsigned int n_chunks =
        (this->n_blocks() + communication_block_size - 1) /
        communication_block_size;
      for (unsigned int c = 0; c < n_chunks; ++c)
        {
          const unsigned int start = c * communication_block_size;
          const unsigned int end =
            std::min(start + communication_block_size, this->n_blocks());

          // In order to avoid conflict with possible other ongoing
          // communication requests (from LA::distributed::Vector that supports
          // unfinished requests), add 100 to the communication tag (the first
          // 100 can be used by normal vectors)
          for (unsigned int block = start; block < end; ++block)
            this->block(block).update_ghost_values_start(block - start + 100);
          for (unsigned int block = start; block < end; ++block)
            this->block(block).update_ghost_values_finish();
        }
    }



    template <typename Number>
    void
    BlockVector<Number>::zero_out_ghosts() const
    {
      for (unsigned int block = 0; block < this->n_blocks(); ++block)
        this->block(block).zero_out_ghosts();
    }



    template <typename Number>
    bool
    BlockVector<Number>::has_ghost_elements() const
    {
      bool has_ghost_elements = false;
      for (unsigned int block = 0; block < this->n_blocks(); ++block)
        if (this->block(block).has_ghost_elements() == true)
          has_ghost_elements = true;
      return has_ghost_elements;
    }



    template <typename Number>
    void
    BlockVector<Number>::reinit(const VectorSpaceVector<Number> &V,
                                const bool omit_zeroing_entries)
    {
      Assert(dynamic_cast<const BlockVector<Number> *>(&V) != nullptr,
             ExcVectorTypeNotCompatible());
      const BlockVector<Number> &down_V =
        dynamic_cast<const BlockVector<Number> &>(V);
      reinit(down_V, omit_zeroing_entries);
    }


    template <typename Number>
    BlockVector<Number> &
    BlockVector<Number>::operator*=(const Number factor)
    {
      for (unsigned int block = 0; block < this->n_blocks(); ++block)
        this->block(block) *= factor;
      return *this;
    }



    template <typename Number>
    BlockVector<Number> &
    BlockVector<Number>::operator/=(const Number factor)
    {
      operator*=(static_cast<Number>(1.) / factor);
      return *this;
    }



    template <typename Number>
    void
    BlockVector<Number>::scale(const VectorSpaceVector<Number> &vv)
    {
      // Downcast. Throws an exception if invalid.
      Assert(dynamic_cast<const BlockVector<Number> *>(&vv) != nullptr,
             ExcVectorTypeNotCompatible());
      const BlockVector<Number> &v =
        dynamic_cast<const BlockVector<Number> &>(vv);
      AssertDimension(this->n_blocks(), v.n_blocks());
      for (unsigned int block = 0; block < this->n_blocks(); ++block)
        this->block(block).scale(v.block(block));
    }



    template <typename Number>
    void
    BlockVector<Number>::equ(const Number                     a,
                             const VectorSpaceVector<Number> &vv)
    {
      // Downcast. Throws an exception if invalid.
      Assert(dynamic_cast<const BlockVector<Number> *>(&vv) != nullptr,
             ExcVectorTypeNotCompatible());
      const BlockVector<Number> &v =
        dynamic_cast<const BlockVector<Number> &>(vv);
      AssertDimension(this->n_blocks(), v.n_blocks());
      for (unsigned int block = 0; block < this->n_blocks(); ++block)
        this->block(block).equ(a, v.block(block));
    }



    template <typename Number>
    BlockVector<Number> &
    BlockVector<Number>::operator+=(const VectorSpaceVector<Number> &vv)
    {
      // Downcast. Throws an exception if invalid.
      Assert(dynamic_cast<const BlockVector<Number> *>(&vv) != nullptr,
             ExcVectorTypeNotCompatible());
      const BlockVector<Number> &v =
        dynamic_cast<const BlockVector<Number> &>(vv);
      AssertDimension(this->n_blocks(), v.n_blocks());
      for (unsigned int block = 0; block < this->n_blocks(); ++block)
        this->block(block) += v.block(block);

      return *this;
    }



    template <typename Number>
    BlockVector<Number> &
    BlockVector<Number>::operator-=(const VectorSpaceVector<Number> &vv)
    {
      // Downcast. Throws an exception if invalid.
      Assert(dynamic_cast<const BlockVector<Number> *>(&vv) != nullptr,
             ExcVectorTypeNotCompatible());
      const BlockVector<Number> &v =
        dynamic_cast<const BlockVector<Number> &>(vv);
      AssertDimension(this->n_blocks(), v.n_blocks());
      for (unsigned int block = 0; block < this->n_blocks(); ++block)
        this->block(block) -= v.block(block);

      return *this;
    }



    template <typename Number>
    void
    BlockVector<Number>::add(const Number a)
    {
      for (unsigned int block = 0; block < this->n_blocks(); ++block)
        this->block(block).add(a);
    }



    template <typename Number>
    void
    BlockVector<Number>::add(const Number                     a,
                             const VectorSpaceVector<Number> &vv)
    {
      // Downcast. Throws an exception if invalid.
      Assert(dynamic_cast<const BlockVector<Number> *>(&vv) != nullptr,
             ExcVectorTypeNotCompatible());
      const BlockVector<Number> &v =
        dynamic_cast<const BlockVector<Number> &>(vv);
      AssertDimension(this->n_blocks(), v.n_blocks());
      for (unsigned int block = 0; block < this->n_blocks(); ++block)
        this->block(block).add(a, v.block(block));
    }



    template <typename Number>
    void
    BlockVector<Number>::add(const Number                     a,
                             const VectorSpaceVector<Number> &vv,
                             const Number                     b,
                             const VectorSpaceVector<Number> &ww)
    {
      // Downcast. Throws an exception if invalid.
      Assert(dynamic_cast<const BlockVector<Number> *>(&vv) != nullptr,
             ExcVectorTypeNotCompatible());
      const BlockVector<Number> &v =
        dynamic_cast<const BlockVector<Number> &>(vv);
      AssertDimension(this->n_blocks(), v.n_blocks());
      Assert(dynamic_cast<const BlockVector<Number> *>(&ww) != nullptr,
             ExcVectorTypeNotCompatible());
      const BlockVector<Number> &w =
        dynamic_cast<const BlockVector<Number> &>(ww);
      AssertDimension(this->n_blocks(), v.n_blocks());

      for (unsigned int block = 0; block < this->n_blocks(); ++block)
        this->block(block).add(a, v.block(block), b, w.block(block));
    }



    template <typename Number>
    void
    BlockVector<Number>::sadd(const Number                     x,
                              const Number                     a,
                              const VectorSpaceVector<Number> &vv)
    {
      // Downcast. Throws an exception if invalid.
      Assert(dynamic_cast<const BlockVector<Number> *>(&vv) != nullptr,
             ExcVectorTypeNotCompatible());
      const BlockVector<Number> &v =
        dynamic_cast<const BlockVector<Number> &>(vv);
      AssertDimension(this->n_blocks(), v.n_blocks());
      for (unsigned int block = 0; block < this->n_blocks(); ++block)
        this->block(block).sadd(x, a, v.block(block));
    }



    template <typename Number>
    void
    BlockVector<Number>::sadd(const Number x, const BlockVector<Number> &v)
    {
      AssertDimension(this->n_blocks(), v.n_blocks());
      for (unsigned int block = 0; block < this->n_blocks(); ++block)
        this->block(block).sadd(x, v.block(block));
    }



    template <typename Number>
    template <typename OtherNumber>
    void
    BlockVector<Number>::add(const std::vector<size_type> &       indices,
                             const ::dealii::Vector<OtherNumber> &values)
    {
      for (size_type i = 0; i < indices.size(); ++i)
        (*this)(indices[i]) += values[i];
    }



    template <typename Number>
    void
    BlockVector<Number>::add(const std::vector<size_type> &indices,
                             const std::vector<Number> &   values)
    {
      for (size_type i = 0; i < indices.size(); ++i)
        (*this)(indices[i]) += values[i];
    }



    template <typename Number>
    bool
    BlockVector<Number>::all_zero() const
    {
      Assert(this->n_blocks() > 0, ExcEmptyObject());

      // use int instead of bool. in order to make global reduction operations
      // work also when MPI_Init was not called, only call MPI_Allreduce
      // commands when there is more than one processor (note that reinit()
      // functions handle this case correctly through the job_supports_mpi()
      // query). this is the same in all the functions below
      int local_result = -1;
      for (unsigned int i = 0; i < this->n_blocks(); ++i)
        local_result =
          std::max(local_result,
                   (this->block(i).linfty_norm_local() == 0) ? -1 : 0);

      if (this->block(0).partitioner->n_mpi_processes() > 1)
        return -Utilities::MPI::max(
          local_result, this->block(0).partitioner->get_mpi_communicator());
      else
        return local_result;
    }



    template <typename Number>
    Number BlockVector<Number>::
           operator*(const VectorSpaceVector<Number> &vv) const
    {
      Assert(this->n_blocks() > 0, ExcEmptyObject());

      // Downcast. Throws an exception if invalid.
      Assert(dynamic_cast<const BlockVector<Number> *>(&vv) != nullptr,
             ExcVectorTypeNotCompatible());
      const BlockVector<Number> &v =
        dynamic_cast<const BlockVector<Number> &>(vv);
      AssertDimension(this->n_blocks(), v.n_blocks());

      Number local_result = Number();
      for (unsigned int i = 0; i < this->n_blocks(); ++i)
        local_result += this->block(i).inner_product_local(v.block(i));

      if (this->block(0).partitioner->n_mpi_processes() > 1)
        return Utilities::MPI::sum(
          local_result, this->block(0).partitioner->get_mpi_communicator());
      else
        return local_result;
    }



    template <typename Number>
    inline Number
    BlockVector<Number>::mean_value() const
    {
      Assert(this->n_blocks() > 0, ExcEmptyObject());

      Number local_result = Number();
      for (unsigned int i = 0; i < this->n_blocks(); ++i)
        local_result +=
          this->block(i).mean_value_local() *
          static_cast<real_type>(this->block(i).partitioner->local_size());

      if (this->block(0).partitioner->n_mpi_processes() > 1)
        return Utilities::MPI::sum(
                 local_result,
                 this->block(0).partitioner->get_mpi_communicator()) /
               static_cast<real_type>(this->size());
      else
        return local_result / static_cast<real_type>(this->size());
    }



    template <typename Number>
    inline typename BlockVector<Number>::real_type
    BlockVector<Number>::l1_norm() const
    {
      Assert(this->n_blocks() > 0, ExcEmptyObject());

      real_type local_result = real_type();
      for (unsigned int i = 0; i < this->n_blocks(); ++i)
        local_result += this->block(i).l1_norm_local();

      if (this->block(0).partitioner->n_mpi_processes() > 1)
        return Utilities::MPI::sum(
          local_result, this->block(0).partitioner->get_mpi_communicator());
      else
        return local_result;
    }



    template <typename Number>
    inline typename BlockVector<Number>::real_type
    BlockVector<Number>::norm_sqr() const
    {
      Assert(this->n_blocks() > 0, ExcEmptyObject());

      real_type local_result = real_type();
      for (unsigned int i = 0; i < this->n_blocks(); ++i)
        local_result += this->block(i).norm_sqr_local();

      if (this->block(0).partitioner->n_mpi_processes() > 1)
        return Utilities::MPI::sum(
          local_result, this->block(0).partitioner->get_mpi_communicator());
      else
        return local_result;
    }



    template <typename Number>
    inline typename BlockVector<Number>::real_type
    BlockVector<Number>::l2_norm() const
    {
      return std::sqrt(norm_sqr());
    }



    template <typename Number>
    inline typename BlockVector<Number>::real_type
    BlockVector<Number>::lp_norm(const real_type p) const
    {
      Assert(this->n_blocks() > 0, ExcEmptyObject());

      real_type local_result = real_type();
      for (unsigned int i = 0; i < this->n_blocks(); ++i)
        local_result += std::pow(this->block(i).lp_norm_local(p), p);

      if (this->block(0).partitioner->n_mpi_processes() > 1)
        return std::pow(Utilities::MPI::sum(
                          local_result,
                          this->block(0).partitioner->get_mpi_communicator()),
                        static_cast<real_type>(1.0 / p));
      else
        return std::pow(local_result, static_cast<real_type>(1.0 / p));
    }



    template <typename Number>
    inline typename BlockVector<Number>::real_type
    BlockVector<Number>::linfty_norm() const
    {
      Assert(this->n_blocks() > 0, ExcEmptyObject());

      real_type local_result = real_type();
      for (unsigned int i = 0; i < this->n_blocks(); ++i)
        local_result =
          std::max(local_result, this->block(i).linfty_norm_local());

      if (this->block(0).partitioner->n_mpi_processes() > 1)
        return Utilities::MPI::max(
          local_result, this->block(0).partitioner->get_mpi_communicator());
      else
        return local_result;
    }



    template <typename Number>
    inline Number
    BlockVector<Number>::add_and_dot(const Number                     a,
                                     const VectorSpaceVector<Number> &vv,
                                     const VectorSpaceVector<Number> &ww)
    {
      // Downcast. Throws an exception if invalid.
      Assert(dynamic_cast<const BlockVector<Number> *>(&vv) != nullptr,
             ExcVectorTypeNotCompatible());
      const BlockVector<Number> &v =
        dynamic_cast<const BlockVector<Number> &>(vv);
      AssertDimension(this->n_blocks(), v.n_blocks());

      // Downcast. Throws an exception if invalid.
      Assert(dynamic_cast<const BlockVector<Number> *>(&ww) != nullptr,
             ExcVectorTypeNotCompatible());
      const BlockVector<Number> &w =
        dynamic_cast<const BlockVector<Number> &>(ww);
      AssertDimension(this->n_blocks(), w.n_blocks());

      Number local_result = Number();
      for (unsigned int i = 0; i < this->n_blocks(); ++i)
        local_result +=
          this->block(i).add_and_dot_local(a, v.block(i), w.block(i));

      if (this->block(0).partitioner->n_mpi_processes() > 1)
        return Utilities::MPI::sum(
          local_result, this->block(0).partitioner->get_mpi_communicator());
      else
        return local_result;
    }



    template <typename Number>
    inline void
    BlockVector<Number>::swap(BlockVector<Number> &v)
    {
      Assert(this->n_blocks() == v.n_blocks(),
             ExcDimensionMismatch(this->n_blocks(), v.n_blocks()));

      for (size_type i = 0; i < this->n_blocks(); ++i)
        dealii::swap(this->components[i], v.components[i]);
      dealii::swap(this->block_indices, v.block_indices);
    }



    template <typename Number>
    typename BlockVector<Number>::size_type
    BlockVector<Number>::size() const
    {
      return this->block_indices.total_size();
    }



    template <typename Number>
    inline void
    BlockVector<Number>::import(const LinearAlgebra::ReadWriteVector<Number> &,
                                VectorOperation::values,
                                std::shared_ptr<const CommunicationPatternBase>)
    {
      AssertThrow(false, ExcNotImplemented());
    }



    template <typename Number>
    IndexSet
    BlockVector<Number>::locally_owned_elements() const
    {
      IndexSet is(size());

      // copy index sets from blocks into the global one, shifted by the
      // appropriate amount for each block
      for (unsigned int b = 0; b < this->n_blocks(); ++b)
        {
          IndexSet x = this->block(b).locally_owned_elements();
          is.add_indices(x, this->block_indices.block_start(b));
        }

      is.compress();

      return is;
    }

    template <typename Number>
    void
    BlockVector<Number>::print(std::ostream &     out,
                               const unsigned int precision,
                               const bool         scientific,
                               const bool         across) const
    {
      for (unsigned int b = 0; b < this->n_blocks(); ++b)
        this->block(b).print(out, precision, scientific, across);
    }



    template <typename Number>
    std::size_t
    BlockVector<Number>::memory_consumption() const
    {
      return (MemoryConsumption::memory_consumption(this->block_indices) +
              MemoryConsumption::memory_consumption(this->components));
    }



    namespace internal
    {
      template <typename FullMatrixType>
      inline void
      set_symmetric(FullMatrixType &, const bool)
      {}

      template <typename NumberType>
      inline void
      set_symmetric(LAPACKFullMatrix<NumberType> &matrix, const bool symmetric)
      {
        if (symmetric)
          matrix.set_property(LAPACKSupport::symmetric);
        else
          matrix.set_property(LAPACKSupport::general);
      }
    } // namespace internal



    template <typename Number>
    template <typename FullMatrixType>
    void
    BlockVector<Number>::multivector_inner_product(FullMatrixType &matrix,
                                                   const BlockVector<Number> &V,
                                                   const bool symmetric) const
    {
      const unsigned int m = this->n_blocks();
      const unsigned int n = V.n_blocks();

      // in case one vector is empty and the second one is not, the
      // FullMatrix resized to (m,n) will have 0 both in m() and n()
      // which is how TableBase<N,T>::reinit() works as of deal.ii@8.5.0.
      // Since in this case there is nothing to do anyway -- return immediately.
      if (n == 0 || m == 0)
        return;

      Assert(matrix.m() == m, ExcDimensionMismatch(matrix.m(), m));
      Assert(matrix.n() == n, ExcDimensionMismatch(matrix.n(), n));

      // reset the matrix
      matrix = typename FullMatrixType::value_type(0.0);

      internal::set_symmetric(matrix, symmetric);
      if (symmetric)
        {
          Assert(m == n, ExcDimensionMismatch(m, n));

          for (unsigned int i = 0; i < m; i++)
            for (unsigned int j = i; j < n; j++)
              matrix(i, j) = this->block(i).inner_product_local(V.block(j));

          for (unsigned int i = 0; i < m; i++)
            for (unsigned int j = i + 1; j < n; j++)
              matrix(j, i) = matrix(i, j);
        }
      else
        {
          for (unsigned int i = 0; i < m; ++i)
            for (unsigned int j = 0; j < n; ++j)
              matrix(i, j) = this->block(i).inner_product_local(V.block(j));
        }

      Utilities::MPI::sum(matrix,
                          this->block(0).get_mpi_communicator(),
                          matrix);
    }



    template <typename Number>
    template <typename FullMatrixType>
    Number
    BlockVector<Number>::multivector_inner_product_with_metric(
      const FullMatrixType &     matrix,
      const BlockVector<Number> &V,
      const bool                 symmetric) const
    {
      Number res = Number(0.);

      const unsigned int m = this->n_blocks();
      const unsigned int n = V.n_blocks();

      // in case one vector is empty and the second one is not, the
      // FullMatrix resized to (m,n) will have 0 both in m() and n()
      // which is how TableBase<N,T>::reinit() works.
      // Since in this case there is nothing to do anyway -- return immediately.
      if (n == 0 || m == 0)
        return res;

      Assert(matrix.m() == m, ExcDimensionMismatch(matrix.m(), m));
      Assert(matrix.n() == n, ExcDimensionMismatch(matrix.n(), n));

      if (symmetric)
        {
          Assert(m == n, ExcDimensionMismatch(m, n));

          for (unsigned int i = 0; i < m; i++)
            {
              res +=
                matrix(i, i) * this->block(i).inner_product_local(V.block(i));
              for (unsigned int j = i + 1; j < n; j++)
                res += 2. * matrix(i, j) *
                       this->block(i).inner_product_local(V.block(j));
            }
        }
      else
        {
          for (unsigned int i = 0; i < m; i++)
            for (unsigned int j = 0; j < n; j++)
              res +=
                matrix(i, j) * this->block(i).inner_product_local(V.block(j));
        }

      return Utilities::MPI::sum(res, this->block(0).get_mpi_communicator());
    }



    template <typename Number>
    template <typename FullMatrixType>
    void
    BlockVector<Number>::mmult(BlockVector<Number> & V,
                               const FullMatrixType &matrix,
                               const Number          s,
                               const Number          b) const
    {
      const unsigned int m = this->n_blocks();
      const unsigned int n = V.n_blocks();

      // in case one vector is empty and the second one is not, the
      // FullMatrix resized to (m,n) will have 0 both in m() and n()
      // which is how TableBase<N,T>::reinit() works.
      // Since in this case there is nothing to do anyway -- return immediately.
      if (n == 0 || m == 0)
        return;

      Assert(matrix.m() == m, ExcDimensionMismatch(matrix.m(), m));
      Assert(matrix.n() == n, ExcDimensionMismatch(matrix.n(), n));

      for (unsigned int i = 0; i < n; i++)
        {
          // below we make this work gracefully for identity-like matrices in
          // which case the two loops over j won't do any work as A(j,i)==0
          const unsigned int k = std::min(i, m - 1);
          V.block(i).sadd_local(s, matrix(k, i) * b, this->block(k));
          for (unsigned int j = 0; j < k; j++)
            V.block(i).add_local(matrix(j, i) * b, this->block(j));
          for (unsigned int j = k + 1; j < m; j++)
            V.block(i).add_local(matrix(j, i) * b, this->block(j));
        }

      if (V.block(0).vector_is_ghosted)
        {
          for (unsigned int i = 0; i < n; i++)
            Assert(V.block(i).vector_is_ghosted,
                   ExcMessage(
                     "All blocks should be either in ghosted state or not."));

          V.update_ghost_values();
        }
    }



  } // end of namespace distributed

} // namespace LinearAlgebra

DEAL_II_NAMESPACE_CLOSE

#endif
