// ---------------------------------------------------------------------
//
// Copyright (C) 1999 - 2016 by the deal.II authors
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

#ifndef dealii__parallel_block_vector_templates_h
#define dealii__parallel_block_vector_templates_h


#include <deal.II/base/config.h>
#include <deal.II/lac/la_parallel_block_vector.h>
#include <deal.II/lac/petsc_parallel_block_vector.h>
#include <deal.II/lac/trilinos_parallel_block_vector.h>


DEAL_II_NAMESPACE_OPEN


namespace LinearAlgebra
{
  namespace distributed
  {
    template <typename Number>
    BlockVector<Number>::BlockVector (const size_type n_blocks,
                                      const size_type block_size)
    {
      reinit (n_blocks, block_size);
    }



    template <typename Number>
    BlockVector<Number>::BlockVector (const std::vector<size_type> &n)
    {
      reinit (n, false);
    }


    template <typename Number>
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
    BlockVector<Number>::BlockVector (const BlockVector<OtherNumber> &v)
    {
      reinit (v, true);
      *this = v;
    }

#endif



    template <typename Number>
    void BlockVector<Number>::reinit (const size_type n_bl,
                                      const size_type bl_sz,
                                      const bool         omit_zeroing_entries)
    {
      std::vector<size_type> n(n_bl, bl_sz);
      reinit(n, omit_zeroing_entries);
    }


    template <typename Number>
    void BlockVector<Number>::reinit (const std::vector<size_type> &n,
                                      const bool                    omit_zeroing_entries)
    {
      this->block_indices.reinit (n);
      if (this->components.size() != this->n_blocks())
        this->components.resize(this->n_blocks());

      for (size_type i=0; i<this->n_blocks(); ++i)
        this->components[i].reinit(n[i], omit_zeroing_entries);
    }



    template <typename Number>
    template <typename Number2>
    void BlockVector<Number>::reinit (const BlockVector<Number2> &v,
                                      const bool omit_zeroing_entries)
    {
      this->block_indices = v.get_block_indices();
      if (this->components.size() != this->n_blocks())
        this->components.resize(this->n_blocks());

      for (unsigned int i=0; i<this->n_blocks(); ++i)
        this->block(i).reinit(v.block(i), omit_zeroing_entries);
    }



    template <typename Number>
    BlockVector<Number>::~BlockVector ()
    {}



    template <typename Number>
    BlockVector<Number> &
    BlockVector<Number>::operator = (const value_type s)
    {
      AssertIsFinite(s);

      BaseClass::operator = (s);
      return *this;
    }



    template <typename Number>
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
    BlockVector<Number> &
    BlockVector<Number>::operator = (const Vector<Number> &v)
    {
      BaseClass::operator = (v);
      return *this;
    }



    template <typename Number>
    template <typename Number2>
    BlockVector<Number> &
    BlockVector<Number>::operator = (const BlockVector<Number2> &v)
    {
      reinit (v, true);
      BaseClass::operator = (v);
      return *this;
    }



#ifdef DEAL_II_WITH_PETSC

    template <typename Number>
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
    void
    BlockVector<Number>::compress (::dealii::VectorOperation::values operation)
    {
      // start all requests for all blocks before finishing the transfers as
      // this saves repeated synchronizations. In order to avoid conflict with
      // possible other ongoing communication requests (from
      // LA::distributed::Vector that supports unfinished requests), add an
      // arbitrary number 8273 to the communication tag
      for (unsigned int block=0; block<this->n_blocks(); ++block)
        this->block(block).compress_start(block + 8273, operation);
      for (unsigned int block=0; block<this->n_blocks(); ++block)
        this->block(block).compress_finish(operation);
    }



    template <typename Number>
    void
    BlockVector<Number>::update_ghost_values () const
    {
      // In order to avoid conflict with possible other ongoing communication
      // requests (from LA::distributed::Vector that supports unfinished
      // requests), add an arbitrary number 9923 to the communication tag
      for (unsigned int block=0; block<this->n_blocks(); ++block)
        this->block(block).update_ghost_values_start(block + 9923);
      for (unsigned int block=0; block<this->n_blocks(); ++block)
        this->block(block).update_ghost_values_finish();
    }



    template <typename Number>
    void
    BlockVector<Number>::zero_out_ghosts ()
    {
      for (unsigned int block=0; block<this->n_blocks(); ++block)
        this->block(block).zero_out_ghosts();
    }



    template <typename Number>
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
    BlockVector<Number> &
    BlockVector<Number>::operator *= (const Number factor)
    {
      for (unsigned int block=0; block<this->n_blocks(); ++block)
        this->block(block) *= factor;
      return *this;
    }



    template <typename Number>
    BlockVector<Number> &
    BlockVector<Number>::operator /= (const Number factor)
    {
      operator *= (static_cast<Number>(1.)/factor);
      return *this;
    }



    template <typename Number>
    void
    BlockVector<Number>::scale (const VectorSpaceVector<Number> &vv)
    {
      // Downcast. Throws an exception if invalid.
      Assert(dynamic_cast<const BlockVector<Number> *>(&vv)!=NULL,
             ExcVectorTypeNotCompatible());
      const BlockVector<Number> &v = dynamic_cast<const BlockVector<Number> &>(vv);
      AssertDimension(this->n_blocks(), v.n_blocks());
      for (unsigned int block=0; block<this->n_blocks(); ++block)
        this->block(block).scale(v.block(block));
    }



    template <typename Number>
    void
    BlockVector<Number>::equ (const Number a,
                              const VectorSpaceVector<Number> &vv)
    {
      // Downcast. Throws an exception if invalid.
      Assert(dynamic_cast<const BlockVector<Number> *>(&vv)!=NULL,
             ExcVectorTypeNotCompatible());
      const BlockVector<Number> &v = dynamic_cast<const BlockVector<Number> &>(vv);
      AssertDimension(this->n_blocks(), v.n_blocks());
      for (unsigned int block=0; block<this->n_blocks(); ++block)
        this->block(block).equ(a, v.block(block));
    }



    template <typename Number>
    void
    BlockVector<Number>::equ (const Number a,
                              const BlockVector<Number> &v,
                              const Number b,
                              const BlockVector<Number> &w)
    {
      AssertDimension(this->n_blocks(), v.n_blocks());
      AssertDimension(this->n_blocks(), w.n_blocks());
      for (unsigned int block=0; block<this->n_blocks(); ++block)
        this->block(block).equ(a, v.block(block), b, w.block(block));
    }



    template <typename Number>
    BlockVector<Number> &
    BlockVector<Number>::operator += (const VectorSpaceVector<Number> &vv)
    {
      // Downcast. Throws an exception if invalid.
      Assert(dynamic_cast<const BlockVector<Number> *>(&vv)!=NULL,
             ExcVectorTypeNotCompatible());
      const BlockVector<Number> &v = dynamic_cast<const BlockVector<Number> &>(vv);
      AssertDimension(this->n_blocks(), v.n_blocks());
      for (unsigned int block=0; block<this->n_blocks(); ++block)
        this->block(block) += v.block(block);

      return *this;
    }



    template <typename Number>
    BlockVector<Number> &
    BlockVector<Number>::operator -= (const VectorSpaceVector<Number> &vv)
    {
      // Downcast. Throws an exception if invalid.
      Assert(dynamic_cast<const BlockVector<Number> *>(&vv)!=NULL,
             ExcVectorTypeNotCompatible());
      const BlockVector<Number> &v = dynamic_cast<const BlockVector<Number> &>(vv);
      AssertDimension(this->n_blocks(), v.n_blocks());
      for (unsigned int block=0; block<this->n_blocks(); ++block)
        this->block(block) -= v.block(block);

      return *this;
    }



    template <typename Number>
    void
    BlockVector<Number>::add (const Number a)
    {
      for (unsigned int block=0; block<this->n_blocks(); ++block)
        this->block(block).add(a);
    }



    template <typename Number>
    void
    BlockVector<Number>::add (const Number a,
                              const VectorSpaceVector<Number> &vv)
    {
      // Downcast. Throws an exception if invalid.
      Assert(dynamic_cast<const BlockVector<Number> *>(&vv)!=NULL,
             ExcVectorTypeNotCompatible());
      const BlockVector<Number> &v = dynamic_cast<const BlockVector<Number> &>(vv);
      AssertDimension(this->n_blocks(), v.n_blocks());
      for (unsigned int block=0; block<this->n_blocks(); ++block)
        this->block(block).add(a, v.block(block));
    }



    template <typename Number>
    void
    BlockVector<Number>::add (const Number a,
                              const VectorSpaceVector<Number> &vv,
                              const Number b,
                              const VectorSpaceVector<Number> &ww)
    {
      // Downcast. Throws an exception if invalid.
      Assert(dynamic_cast<const BlockVector<Number> *>(&vv)!=NULL,
             ExcVectorTypeNotCompatible());
      const BlockVector<Number> &v = dynamic_cast<const BlockVector<Number> &>(vv);
      AssertDimension(this->n_blocks(), v.n_blocks());
      Assert(dynamic_cast<const BlockVector<Number> *>(&ww)!=NULL,
             ExcVectorTypeNotCompatible());
      const BlockVector<Number> &w = dynamic_cast<const BlockVector<Number> &>(ww);
      AssertDimension(this->n_blocks(), v.n_blocks());

      for (unsigned int block=0; block<this->n_blocks(); ++block)
        this->block(block).add(a, v.block(block), b, w.block(block));
    }



    template <typename Number>
    void
    BlockVector<Number>::sadd (const Number x,
                               const Number a,
                               const VectorSpaceVector<Number> &vv)
    {
      // Downcast. Throws an exception if invalid.
      Assert(dynamic_cast<const BlockVector<Number> *>(&vv)!=NULL,
             ExcVectorTypeNotCompatible());
      const BlockVector<Number> &v = dynamic_cast<const BlockVector<Number> &>(vv);
      AssertDimension(this->n_blocks(), v.n_blocks());
      for (unsigned int block=0; block<this->n_blocks(); ++block)
        this->block(block).sadd(x, a, v.block(block));
    }



    template <typename Number>
    void
    BlockVector<Number>::sadd (const Number x,
                               const BlockVector<Number> &v)
    {
      AssertDimension(this->n_blocks(), v.n_blocks());
      for (unsigned int block=0; block<this->n_blocks(); ++block)
        this->block(block).sadd(x, v.block(block));
    }



    template <typename Number>
    void
    BlockVector<Number>::sadd (const Number x,
                               const Number a,
                               const BlockVector<Number> &v,
                               const Number b,
                               const BlockVector<Number> &w)
    {
      AssertDimension(this->n_blocks(), v.n_blocks());
      AssertDimension(this->n_blocks(), w.n_blocks());
      for (unsigned int block=0; block<this->n_blocks(); ++block)
        this->block(block).sadd(x, a, v.block(block), b, w.block(block));
    }



    template <typename Number>
    template <typename OtherNumber>
    void
    BlockVector<Number>::add (const std::vector<size_type>        &indices,
                              const ::dealii::Vector<OtherNumber> &values)
    {
      for (size_type i=0; i<indices.size(); ++i)
        (*this)(indices[i]) += values[i];
    }



    template <typename Number>
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
    Number
    BlockVector<Number>::operator * (const VectorSpaceVector<Number> &vv) const
    {
      Assert (this->n_blocks() > 0, ExcEmptyObject());

      // Downcast. Throws an exception if invalid.
      Assert(dynamic_cast<const BlockVector<Number> *>(&vv)!=NULL,
             ExcVectorTypeNotCompatible());
      const BlockVector<Number> &v = dynamic_cast<const BlockVector<Number> &>(vv);
      AssertDimension(this->n_blocks(), v.n_blocks());

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
      Assert (this->n_blocks() > 0, ExcEmptyObject());

      real_type local_result = real_type();
      for (unsigned int i=0; i<this->n_blocks(); ++i)
        local_result += this->block(i).norm_sqr_local();

      if (this->block(0).partitioner->n_mpi_processes() > 1)
        return std::sqrt(Utilities::MPI::sum (local_result,
                                              this->block(0).partitioner->get_communicator()));
      else
        return std::sqrt(local_result);
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
    Number
    BlockVector<Number>::add_and_dot (const Number                     a,
                                      const VectorSpaceVector<Number> &vv,
                                      const VectorSpaceVector<Number> &ww)
    {
      // Downcast. Throws an exception if invalid.
      Assert(dynamic_cast<const BlockVector<Number> *>(&vv)!=NULL,
             ExcVectorTypeNotCompatible());
      const BlockVector<Number> &v = dynamic_cast<const BlockVector<Number> &>(vv);
      AssertDimension(this->n_blocks(), v.n_blocks());

      // Downcast. Throws an exception if invalid.
      Assert(dynamic_cast<const BlockVector<Number> *>(&ww)!=NULL,
             ExcVectorTypeNotCompatible());
      const BlockVector<Number> &w = dynamic_cast<const BlockVector<Number> &>(ww);
      AssertDimension(this->n_blocks(), w.n_blocks());

      Number local_result = Number();
      for (unsigned int i=0; i<this->n_blocks(); ++i)
        local_result += this->block(i).add_and_dot_local(a, v.block(i), w.block(i));

      if (this->block(0).partitioner->n_mpi_processes() > 1)
        return Utilities::MPI::sum (local_result,
                                    this->block(0).partitioner->get_communicator());
      else
        return local_result;
    }



    template <typename Number>
    inline
    void
    BlockVector<Number>::swap (BlockVector<Number> &v)
    {
      Assert (this->n_blocks() == v.n_blocks(),
              ExcDimensionMismatch(this->n_blocks(), v.n_blocks()));

      for (size_type i=0; i<this->n_blocks(); ++i)
        dealii::swap (this->components[i], v.components[i]);
      dealii::swap (this->block_indices, v.block_indices);
    }



    template <typename Number>
    typename BlockVector<Number>::size_type
    BlockVector<Number>::size () const
    {
      return this->block_indices.total_size();
    }



    template <typename Number>
    inline
    void
    BlockVector<Number>::import(const LinearAlgebra::ReadWriteVector<Number> &,
                                VectorOperation::values,
                                std_cxx11::shared_ptr<const CommunicationPatternBase>)
    {
      AssertThrow(false, ExcNotImplemented());
    }



    template <typename Number>
    IndexSet
    BlockVector<Number>::locally_owned_elements () const
    {
      IndexSet is (size());

      // copy index sets from blocks into the global one, shifted by the
      // appropriate amount for each block
      for (unsigned int b=0; b<this->n_blocks(); ++b)
        {
          IndexSet x = this->block(b).locally_owned_elements();
          is.add_indices(x, this->block_indices.block_start(b));
        }

      is.compress();

      return is;
    }

    template <typename Number>
    void
    BlockVector<Number>::print(std::ostream &out,
                               const unsigned int precision,
                               const bool scientific,
                               const bool across) const
    {
      for (unsigned int b=0; b<this->n_blocks(); ++b)
        this->block(b).print(out, precision, scientific, across);
    }



    template <typename Number>
    std::size_t
    BlockVector<Number>::memory_consumption () const
    {
      std::size_t mem = sizeof(this->n_blocks());
      for (size_type i=0; i<this->components.size(); ++i)
        mem += MemoryConsumption::memory_consumption (this->components[i]);
      mem += MemoryConsumption::memory_consumption (this->block_indices);
      return mem;
    }

  } // end of namespace distributed

} // end of namespace parallel

DEAL_II_NAMESPACE_CLOSE

#endif
