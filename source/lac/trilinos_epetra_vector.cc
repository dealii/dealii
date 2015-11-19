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

#include <deal.II/lac/trilinos_epetra_vector.h>

#ifdef DEAL_II_WITH_TRILINOS
#ifdef DEAL_II_WITH_MPI

#include <deal.II/base/index_set.h>
#include "Epetra_Import.h"
#include "Epetra_Map.h"
#include "Epetra_MpiComm.h"


DEAL_II_NAMESPACE_OPEN

namespace LinearAlgebra
{
  namespace EpetraWrappers
  {
    Vector::Vector()
      :
      vector(new Epetra_FEVector(Epetra_Map(0,1,0,Epetra_MpiComm(MPI_COMM_SELF))))
    {}



    Vector::Vector(const Vector &V)
      :
      vector(new Epetra_FEVector(V.trilinos_vector()))
    {}



    Vector::Vector(const IndexSet &parallel_partitioner,
                   const MPI_Comm &communicator)
      :
      vector(new Epetra_FEVector(parallel_partitioner.make_trilinos_map(communicator,false)))
    {}



    void Vector::reinit(const IndexSet &parallel_partitioner,
                        const MPI_Comm &communicator,
                        const bool      omit_zeroing_entries)
    {
      Epetra_Map input_map = parallel_partitioner.make_trilinos_map(communicator,false);
      if (vector->Map().SameAs(input_map)==false)
        vector.reset(new Epetra_FEVector(input_map));
      else if (omit_zeroing_entries==false)
        {
          const int ierr = vector->PutScalar(0.);
          (void) ierr;
          Assert(ierr==0, ExcTrilinosError(ierr));
        }
    }



    Vector &Vector::operator= (const Vector &V)
    {
      // Distinguish three cases:
      //  - First case: both vectors have the same layout.
      //  - Second case: both vectors have the same size but different layout.
      //  - Third case: the vectors have different size.
      if (vector->Map().SameAs(V.trilinos_vector().Map()))
        *vector = V.trilinos_vector();
      else
        {
          if (size()==V.size())
            {
              Epetra_Import data_exchange(vector->Map(), V.trilinos_vector().Map());

              const int ierr = vector->Import(V.trilinos_vector(), data_exchange, Insert);
              Assert(ierr==0, ExcTrilinosError(ierr));
            }
          else
            vector.reset(new Epetra_FEVector(V.trilinos_vector()));
        }

      return *this;
    }



    Vector &Vector::operator= (const ::dealii::LinearAlgebra::ReadWriteVector<double> &V)
    {
      IndexSet elements_to_write(V.get_stored_elements());
      IndexSet locally_owned_elements(this->locally_owned_elements());
      IndexSet off_proc_elements(elements_to_write);
      off_proc_elements.subtract_set(locally_owned_elements);
      IndexSet on_proc_elements(elements_to_write);
      on_proc_elements.subtract_set(off_proc_elements);

      // Write the elements that locally owned.
      IndexSet::ElementIterator elem = on_proc_elements.begin(),
                                elem_end = on_proc_elements.end();
      for (; elem!=elem_end; ++elem)
        {
          const TrilinosWrappers::types::int_type local_index =
            vector->Map().LID(static_cast<TrilinosWrappers::types::int_type>(*elem));
          (*vector)[0][local_index] = V[*elem];
        }


      // Write the elements off-processor. Cannot use Import function of Trilinos
      // because V is not a Trilinos object. The most straightforward to send
      // the off-processor to the other processors is to AllGather. This way
      // every processor has access to all the off-processor elements. The
      // problem with this method is that if every element in V is
      // off-processor, it is possible that the buffers used to receive the data
      // has twice the global size of the current Vector.
      // TODO This won't scale past a few 10,000's of processors.

      // Get recv_count and the size of the buffer needed to receive the data.
      const Epetra_MpiComm *epetra_comm
        = dynamic_cast<const Epetra_MpiComm *>(&(vector->Comm()));
      MPI_Comm comm = epetra_comm->GetMpiComm();
      int send_buffer_size(off_proc_elements.n_elements());
      int n_procs(0);
      MPI_Comm_size(comm, &n_procs);
      std::vector<int> recv_count(n_procs);
      MPI_Allgather(&send_buffer_size, 1, MPI_INT, &recv_count, 1,
                    MPI_INT, comm);
      std::vector<int> displacement(n_procs,0);
      for (int i=1; i<n_procs; ++i)
        displacement[i] = displacement[i-1] + recv_count[i-1];
      const unsigned int recv_buffer_size = displacement[n_procs-1] +
                                            recv_count[n_procs-1];

      // Write the indices in the send buffer
      std::vector<types::global_dof_index> send_index_buffer(send_buffer_size);
      elem = off_proc_elements.begin();
      for (int i=0; i<send_buffer_size; ++i)
        {
          send_index_buffer[i] = *elem;
          ++elem;
        }
      // Send the index of the off-processor elements
      std::vector<types::global_dof_index> recv_index_buffer(recv_buffer_size);
      MPI_Allgatherv(&send_index_buffer[0], send_buffer_size, DEAL_II_DOF_INDEX_MPI_TYPE,
                     &recv_index_buffer[0], &recv_count[0], &displacement[0],
                     DEAL_II_DOF_INDEX_MPI_TYPE, comm);

      // Write the values in the send buffer
      std::vector<double> send_value_buffer(send_buffer_size);
      elem = off_proc_elements.begin();
      for (int i=0; i<send_buffer_size; ++i)
        {
          send_value_buffer[i] = V[*elem];
          ++elem;
        }
      // Send the value of the off-processor elements
      std::vector<double> recv_value_buffer(recv_buffer_size);
      MPI_Allgatherv(&send_value_buffer[0], send_buffer_size, MPI_DOUBLE,
                     &recv_value_buffer[0], &recv_count[0], &displacement[0], MPI_DOUBLE, comm);

      // Write the data in the vector
      for (unsigned int i=0; i<recv_buffer_size; ++i)
        {
          if (locally_owned_elements.is_element(recv_index_buffer[i]))
            {
              const TrilinosWrappers::types::int_type local_index =
                vector->Map().LID(static_cast<TrilinosWrappers::types::int_type>(
                                    recv_index_buffer[i]));
              (*vector)[0][local_index] = recv_value_buffer[i];
            }
        }

      return *this;
    }



    VectorSpaceVector<double> &Vector::operator*= (const double factor)
    {
      AssertIsFinite(factor);
      vector->Scale(factor);

      return *this;
    }



    VectorSpaceVector<double> &Vector::operator/= (const double factor)
    {
      AssertIsFinite(factor);
      Assert(factor!=0., ExcZero());
      *this *= 1./factor;

      return *this;
    }



    VectorSpaceVector<double> &Vector::operator+= (const VectorSpaceVector<double> &V)
    {
      // Check that casting will work.
      Assert(dynamic_cast<const Vector *>(&V)!=NULL, ExcVectorTypeNotCompatible());

      // Downcast V. If fails, throws an exception.
      const Vector &down_V = dynamic_cast<const Vector &>(V);
      // If the maps are the same we can Update right away.
      if (vector->Map().SameAs(down_V.trilinos_vector().Map()))
        vector->Update(1., down_V.trilinos_vector(), 1.);
      else
        {
          Assert(this->size()==down_V.size(),
                 ExcDimensionMismatch(this->size(), down_V.size()));

#if DEAL_II_TRILINOS_VERSION_GTE(11,11,0)
          Epetra_Import data_exchange (vector->Map(), down_V.trilinos_vector().Map());
          int ierr = vector->Import(down_V.trilinos_vector(), data_exchange, Epetra_AddLocalAlso);
          (void) ierr;
          Assert(ierr==0, ExcTrilinosError(ierr));
#else
          // In versions older than 11.11 the Import function is broken for adding
          // Hence, we provide a workaround in this case

          Epetra_MultiVector dummy(vector->Map(), 1, false);
          Epetra_Import data_exchange(dummy.Map(), down_V.trilinos_vector().Map());

          int ierr = dummy.Import(down_V.trilinos_vector(), data_exchange, Insert);
          (void) ierr;
          Assert(ierr==0, ExcTrilinosError(ierr));

          ierr = vector->Update(1.0, dummy, 1.0);
          (void) ierr;
          Assert(ierr==0, ExcTrilinosError(ierr));
#endif
        }

      return *this;
    }



    VectorSpaceVector<double> &Vector::operator-= (const VectorSpaceVector<double> &V)
    {
      this->add(-1.,V);

      return *this;
    }



    double Vector::operator* (const VectorSpaceVector<double> &V)
    {
      // Check that casting will work.
      Assert(dynamic_cast<const Vector *>(&V)!=NULL,
             ExcVectorTypeNotCompatible());

      // Downcast V. If fails, throws an exception.
      const Vector &down_V = dynamic_cast<const Vector &>(V);
      Assert(this->size()==down_V.size(),
             ExcDimensionMismatch(this->size(), down_V.size()));
      Assert(vector->Map().SameAs(down_V.trilinos_vector().Map()),
             ExcDifferentParallelPartitioning());

      double result(0.);
      const int ierr = vector->Dot(down_V.trilinos_vector(), &result);
      (void) ierr;
      Assert(ierr==0, ExcTrilinosError(ierr));

      return result;
    }



    void Vector::add(const double a)
    {
      AssertIsFinite(a);
      const unsigned local_size(vector->MyLength());
      for (unsigned int i=0; i<local_size; ++i)
        (*vector)[0][i] += a;
    }



    void Vector::add(const double a, const VectorSpaceVector<double> &V)
    {
      // Check that casting will work.
      Assert(dynamic_cast<const Vector *>(&V)!=NULL,
             ExcVectorTypeNotCompatible());

      // Downcast V. If fails, throws an exception.
      const Vector &down_V = dynamic_cast<const Vector &>(V);
      AssertIsFinite(a);
      Assert(vector->Map().SameAs(down_V.trilinos_vector().Map()),
             ExcDifferentParallelPartitioning());

      const int ierr = vector->Update(a, down_V.trilinos_vector(), 1.);
      Assert(ierr==0, ExcTrilinosError(ierr));
    }



    void Vector::add(const double a, const VectorSpaceVector<double> &V,
                     const double b, const VectorSpaceVector<double> &W)
    {
      // Check that casting will work.
      Assert(dynamic_cast<const Vector *>(&V)!=NULL,
             ExcVectorTypeNotCompatible());
      // Check that casting will work.
      Assert(dynamic_cast<const Vector *>(&W)!=NULL,
             ExcVectorTypeNotCompatible());

      // Downcast V. If fails, throws an exception.
      const Vector &down_V = dynamic_cast<const Vector &>(V);
      // Downcast W. If fails, throws an exception.
      const Vector &down_W = dynamic_cast<const Vector &>(W);
      Assert(vector->Map().SameAs(down_V.trilinos_vector().Map()),
             ExcDifferentParallelPartitioning());
      Assert(vector->Map().SameAs(down_W.trilinos_vector().Map()),
             ExcDifferentParallelPartitioning());
      AssertIsFinite(a);
      AssertIsFinite(b);

      const int ierr = vector->Update(a, down_V.trilinos_vector(), b,
                                      down_W.trilinos_vector(), 1.);
      Assert(ierr==0, ExcTrilinosError(ierr));
    }



    void Vector::sadd(const double s, const double a,
                      const VectorSpaceVector<double> &V)
    {
      // Check that casting will work.
      Assert(dynamic_cast<const Vector *>(&V)!=NULL,
             ExcVectorTypeNotCompatible());

      *this *= s;
      // Downcast V. It fails, throws an exception.
      const Vector &down_V = dynamic_cast<const Vector &>(V);
      Vector tmp(down_V);
      tmp *= a;
      *this += tmp;
    }



    void Vector::scale(const VectorSpaceVector<double> &scaling_factors)
    {
      // Check that casting will work.
      Assert(dynamic_cast<const Vector *>(&scaling_factors)!=NULL,
             ExcVectorTypeNotCompatible());

      // Downcast scaling_factors. If fails, throws an exception.
      const Vector &down_scaling_factors =
        dynamic_cast<const Vector &>(scaling_factors);
      Assert(vector->Map().SameAs(down_scaling_factors.trilinos_vector().Map()),
             ExcDifferentParallelPartitioning());

      const int ierr = vector->Multiply(1.0, down_scaling_factors.trilinos_vector(),
                                        *vector, 0.0);
      Assert(ierr==0, ExcTrilinosError(ierr));
    }



    void Vector::equ(const double a, const VectorSpaceVector<double> &V)
    {
      // Check that casting will work.
      Assert(dynamic_cast<const Vector *>(&V)!=NULL,
             ExcVectorTypeNotCompatible());

      // Downcast V. If fails, throws an exception.
      const Vector &down_V = dynamic_cast<const Vector &>(V);
      // If we don't have the same map, copy.
      if (vector->Map().SameAs(down_V.trilinos_vector().Map())==false)
        this->sadd(0., a, V);
      else
        {
          // Otherwise, just update
          int ierr = vector->Update(a, down_V.trilinos_vector(), 0.);
          Assert(ierr==0, ExcTrilinosError(ierr));
        }
    }



    double Vector::l1_norm()
    {
      double norm(0.);
      int ierr = vector->Norm1(&norm);
      Assert(ierr==0, ExcTrilinosError(ierr));

      return norm;
    }



    double Vector::l2_norm()
    {
      double norm(0.);
      int ierr = vector->Norm2(&norm);
      Assert(ierr==0, ExcTrilinosError(ierr));

      return norm;
    }



    double Vector::linfty_norm()
    {
      double norm(0.);
      int ierr = vector->NormInf(&norm);
      Assert(ierr==0, ExcTrilinosError(ierr));

      return norm;
    }



    double Vector::add_and_dot(const double a,
                               const VectorSpaceVector<double> &V,
                               const VectorSpaceVector<double> &W)
    {
      this->add(a, V);

      return *this * W;
    }



    Vector::size_type Vector::size() const
    {
#ifndef DEAL_II_WITH_64BIT_INDICES
      return vector->GlobalLength();
#else
      return vector->GlobalLength64();
#endif
    }



    MPI_Comm Vector::get_mpi_communicator() const
    {
      const Epetra_MpiComm *epetra_comm
        = dynamic_cast<const Epetra_MpiComm *>(&(vector->Comm()));
      return epetra_comm->GetMpiComm();
    }



    const ::dealii::IndexSet Vector::locally_owned_elements() const
    {
      const Epetra_Map *map = dynamic_cast<const Epetra_Map *>(&(vector->Map()));
      return IndexSet(*map);
    }



    const Epetra_FEVector &Vector::trilinos_vector() const
    {
      return *vector;
    }



    Epetra_FEVector &Vector::trilinos_vector()
    {
      return *vector;
    }



    void Vector::print(std::ostream &out,
                       const unsigned int precision,
                       const bool scientific,
                       const bool across) const
    {
      AssertThrow(out, ExcIO());

      // Get a representation of the vector and loop over all
      // the elements
      double *val;
      int leading_dimension;
      int ierr = vector->ExtractView(&val, &leading_dimension);

      Assert(ierr==0, ExcTrilinosError(ierr));
      out.precision (precision);
      if (scientific)
        out.setf(std::ios::scientific, std::ios::floatfield);
      else
        out.setf(std::ios::fixed, std::ios::floatfield);

      if (across)
        for (size_type i=0; i<size(); ++i)
          out << val[i] << ' ';
      else
        for (size_type i=0; i<size(); ++i)
          out << val[i] << std::endl;
      out << std::endl;

      // restore the representation
      // of the vector
      AssertThrow(out, ExcIO());
    }



    std::size_t Vector::memory_consumption() const
    {
      return sizeof(*this)
             + vector->MyLength()*(sizeof(double)+
                                   sizeof(TrilinosWrappers::types::int_type));
    }
  }
}

DEAL_II_NAMESPACE_CLOSE

#endif

#endif
