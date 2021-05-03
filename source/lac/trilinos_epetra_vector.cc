// ---------------------------------------------------------------------
//
// Copyright (C) 2015 - 2019 by the deal.II authors
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

#include <deal.II/lac/trilinos_epetra_vector.h>

#ifdef DEAL_II_WITH_TRILINOS

#  ifdef DEAL_II_WITH_MPI

#    include <deal.II/base/index_set.h>

#    include <deal.II/lac/read_write_vector.h>

#    include <boost/io/ios_state.hpp>

#    include <Epetra_Import.h>
#    include <Epetra_Map.h>
#    include <Epetra_MpiComm.h>

#    include <memory>


DEAL_II_NAMESPACE_OPEN

namespace LinearAlgebra
{
  namespace EpetraWrappers
  {
    Vector::Vector()
      : Subscriptor()
      , vector(new Epetra_FEVector(
          Epetra_Map(0, 0, 0, Utilities::Trilinos::comm_self())))
    {}



    Vector::Vector(const Vector &V)
      : Subscriptor()
      , vector(new Epetra_FEVector(V.trilinos_vector()))
    {}



    Vector::Vector(const IndexSet &parallel_partitioner,
                   const MPI_Comm &communicator)
      : Subscriptor()
      , vector(new Epetra_FEVector(
          parallel_partitioner.make_trilinos_map(communicator, false)))
    {}



    void
    Vector::reinit(const IndexSet &parallel_partitioner,
                   const MPI_Comm &communicator,
                   const bool      omit_zeroing_entries)
    {
      Epetra_Map input_map =
        parallel_partitioner.make_trilinos_map(communicator, false);
      if (vector->Map().SameAs(input_map) == false)
        vector = std::make_unique<Epetra_FEVector>(input_map);
      else if (omit_zeroing_entries == false)
        {
          const int ierr = vector->PutScalar(0.);
          Assert(ierr == 0, ExcTrilinosError(ierr));
          (void)ierr;
        }
    }



    void
    Vector::reinit(const VectorSpaceVector<double> &V,
                   const bool                       omit_zeroing_entries)
    {
      // Check that casting will work.
      Assert(dynamic_cast<const Vector *>(&V) != nullptr,
             ExcVectorTypeNotCompatible());

      // Downcast V. If fails, throws an exception.
      const Vector &down_V = dynamic_cast<const Vector &>(V);

      reinit(down_V.locally_owned_elements(),
             down_V.get_mpi_communicator(),
             omit_zeroing_entries);
    }



    Vector &
    Vector::operator=(const Vector &V)
    {
      // Distinguish three cases:
      //  - First case: both vectors have the same layout.
      //  - Second case: both vectors have the same size but different layout.
      //  - Third case: the vectors have different size.
      if (vector->Map().SameAs(V.trilinos_vector().Map()))
        *vector = V.trilinos_vector();
      else
        {
          if (size() == V.size())
            {
              Epetra_Import data_exchange(vector->Map(),
                                          V.trilinos_vector().Map());

              const int ierr =
                vector->Import(V.trilinos_vector(), data_exchange, Insert);
              Assert(ierr == 0, ExcTrilinosError(ierr));
              (void)ierr;
            }
          else
            vector = std::make_unique<Epetra_FEVector>(V.trilinos_vector());
        }

      return *this;
    }



    Vector &
    Vector::operator=(const double s)
    {
      Assert(s == 0., ExcMessage("Only 0 can be assigned to a vector."));

      const int ierr = vector->PutScalar(s);
      Assert(ierr == 0, ExcTrilinosError(ierr));
      (void)ierr;

      return *this;
    }



    void
    Vector::import(
      const ReadWriteVector<double> &V,
      VectorOperation::values        operation,
      std::shared_ptr<const Utilities::MPI::CommunicationPatternBase>
        communication_pattern)
    {
      // If no communication pattern is given, create one. Otherwise, use the
      // one given.
      if (communication_pattern == nullptr)
        {
          // The first time import is called, a communication pattern is
          // created. Check if the communication pattern already exists and if
          // it can be reused.
          if ((source_stored_elements.size() !=
               V.get_stored_elements().size()) ||
              (source_stored_elements != V.get_stored_elements()))
            {
              create_epetra_comm_pattern(
                V.get_stored_elements(),
                dynamic_cast<const Epetra_MpiComm &>(vector->Comm()).Comm());
            }
        }
      else
        {
          epetra_comm_pattern =
            std::dynamic_pointer_cast<const CommunicationPattern>(
              communication_pattern);
          AssertThrow(
            epetra_comm_pattern != nullptr,
            ExcMessage(
              std::string("The communication pattern is not of type ") +
              "LinearAlgebra::EpetraWrappers::CommunicationPattern."));
        }

      Epetra_Import import(epetra_comm_pattern->get_epetra_import());

      // The TargetMap and the SourceMap have their roles inverted.
      Epetra_FEVector source_vector(import.TargetMap());
      double *        values = source_vector.Values();
      std::copy(V.begin(), V.end(), values);

      if (operation == VectorOperation::insert)
        vector->Export(source_vector, import, Insert);
      else if (operation == VectorOperation::add)
        vector->Export(source_vector, import, Add);
      else if (operation == VectorOperation::max)
        vector->Export(source_vector, import, Epetra_Max);
      else if (operation == VectorOperation::min)
        vector->Export(source_vector, import, Epetra_Min);
      else
        AssertThrow(false, ExcNotImplemented());
    }



    Vector &
    Vector::operator*=(const double factor)
    {
      AssertIsFinite(factor);
      vector->Scale(factor);

      return *this;
    }



    Vector &
    Vector::operator/=(const double factor)
    {
      AssertIsFinite(factor);
      Assert(factor != 0., ExcZero());
      *this *= 1. / factor;

      return *this;
    }



    Vector &
    Vector::operator+=(const VectorSpaceVector<double> &V)
    {
      // Check that casting will work.
      Assert(dynamic_cast<const Vector *>(&V) != nullptr,
             ExcVectorTypeNotCompatible());

      // Downcast V. If fails, throws an exception.
      const Vector &down_V = dynamic_cast<const Vector &>(V);
      // If the maps are the same we can Update right away.
      if (vector->Map().SameAs(down_V.trilinos_vector().Map()))
        {
          const int ierr = vector->Update(1., down_V.trilinos_vector(), 1.);
          Assert(ierr == 0, ExcTrilinosError(ierr));
          (void)ierr;
        }
      else
        {
          Assert(this->size() == down_V.size(),
                 ExcDimensionMismatch(this->size(), down_V.size()));

#    if DEAL_II_TRILINOS_VERSION_GTE(11, 11, 0)
          Epetra_Import data_exchange(vector->Map(),
                                      down_V.trilinos_vector().Map());
          const int     ierr = vector->Import(down_V.trilinos_vector(),
                                          data_exchange,
                                          Epetra_AddLocalAlso);
          Assert(ierr == 0, ExcTrilinosError(ierr));
          (void)ierr;
#    else
          // In versions older than 11.11 the Import function is broken for
          // adding Hence, we provide a workaround in this case

          Epetra_MultiVector dummy(vector->Map(), 1, false);
          Epetra_Import      data_exchange(dummy.Map(),
                                      down_V.trilinos_vector().Map());

          int ierr =
            dummy.Import(down_V.trilinos_vector(), data_exchange, Insert);
          Assert(ierr == 0, ExcTrilinosError(ierr));

          ierr = vector->Update(1.0, dummy, 1.0);
          Assert(ierr == 0, ExcTrilinosError(ierr));
          (void)ierr;
#    endif
        }

      return *this;
    }



    Vector &
    Vector::operator-=(const VectorSpaceVector<double> &V)
    {
      this->add(-1., V);

      return *this;
    }



    double Vector::operator*(const VectorSpaceVector<double> &V) const
    {
      // Check that casting will work.
      Assert(dynamic_cast<const Vector *>(&V) != nullptr,
             ExcVectorTypeNotCompatible());

      // Downcast V. If fails, throws an exception.
      const Vector &down_V = dynamic_cast<const Vector &>(V);
      Assert(this->size() == down_V.size(),
             ExcDimensionMismatch(this->size(), down_V.size()));
      Assert(vector->Map().SameAs(down_V.trilinos_vector().Map()),
             ExcDifferentParallelPartitioning());

      double    result(0.);
      const int ierr = vector->Dot(down_V.trilinos_vector(), &result);
      Assert(ierr == 0, ExcTrilinosError(ierr));
      (void)ierr;

      return result;
    }



    void
    Vector::add(const double a)
    {
      AssertIsFinite(a);
      const unsigned local_size(vector->MyLength());
      for (unsigned int i = 0; i < local_size; ++i)
        (*vector)[0][i] += a;
    }



    void
    Vector::add(const double a, const VectorSpaceVector<double> &V)
    {
      // Check that casting will work.
      Assert(dynamic_cast<const Vector *>(&V) != nullptr,
             ExcVectorTypeNotCompatible());

      // Downcast V. If fails, throws an exception.
      const Vector &down_V = dynamic_cast<const Vector &>(V);
      AssertIsFinite(a);
      Assert(vector->Map().SameAs(down_V.trilinos_vector().Map()),
             ExcDifferentParallelPartitioning());

      const int ierr = vector->Update(a, down_V.trilinos_vector(), 1.);
      Assert(ierr == 0, ExcTrilinosError(ierr));
      (void)ierr;
    }



    void
    Vector::add(const double                     a,
                const VectorSpaceVector<double> &V,
                const double                     b,
                const VectorSpaceVector<double> &W)
    {
      // Check that casting will work.
      Assert(dynamic_cast<const Vector *>(&V) != nullptr,
             ExcVectorTypeNotCompatible());
      // Check that casting will work.
      Assert(dynamic_cast<const Vector *>(&W) != nullptr,
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

      const int ierr = vector->Update(
        a, down_V.trilinos_vector(), b, down_W.trilinos_vector(), 1.);
      Assert(ierr == 0, ExcTrilinosError(ierr));
      (void)ierr;
    }



    void
    Vector::sadd(const double                     s,
                 const double                     a,
                 const VectorSpaceVector<double> &V)
    {
      // Check that casting will work.
      Assert(dynamic_cast<const Vector *>(&V) != nullptr,
             ExcVectorTypeNotCompatible());

      *this *= s;
      // Downcast V. It fails, throws an exception.
      const Vector &down_V = dynamic_cast<const Vector &>(V);
      Vector        tmp(down_V);
      tmp *= a;
      *this += tmp;
    }



    void
    Vector::scale(const VectorSpaceVector<double> &scaling_factors)
    {
      // Check that casting will work.
      Assert(dynamic_cast<const Vector *>(&scaling_factors) != nullptr,
             ExcVectorTypeNotCompatible());

      // Downcast scaling_factors. If fails, throws an exception.
      const Vector &down_scaling_factors =
        dynamic_cast<const Vector &>(scaling_factors);
      Assert(vector->Map().SameAs(down_scaling_factors.trilinos_vector().Map()),
             ExcDifferentParallelPartitioning());

      const int ierr = vector->Multiply(1.0,
                                        down_scaling_factors.trilinos_vector(),
                                        *vector,
                                        0.0);
      Assert(ierr == 0, ExcTrilinosError(ierr));
      (void)ierr;
    }



    void
    Vector::equ(const double a, const VectorSpaceVector<double> &V)
    {
      // Check that casting will work.
      Assert(dynamic_cast<const Vector *>(&V) != nullptr,
             ExcVectorTypeNotCompatible());

      // Downcast V. If fails, throws an exception.
      const Vector &down_V = dynamic_cast<const Vector &>(V);
      // If we don't have the same map, copy.
      if (vector->Map().SameAs(down_V.trilinos_vector().Map()) == false)
        this->sadd(0., a, V);
      else
        {
          // Otherwise, just update
          int ierr = vector->Update(a, down_V.trilinos_vector(), 0.);
          Assert(ierr == 0, ExcTrilinosError(ierr));
          (void)ierr;
        }
    }



    bool
    Vector::all_zero() const
    {
      // get a representation of the vector and
      // loop over all the elements
      double *      start_ptr = (*vector)[0];
      const double *ptr = start_ptr, *eptr = start_ptr + vector->MyLength();
      unsigned int  flag = 0;
      while (ptr != eptr)
        {
          if (*ptr != 0)
            {
              flag = 1;
              break;
            }
          ++ptr;
        }

      // Check that the vector is zero on _all_ processors.
      const Epetra_MpiComm *mpi_comm =
        dynamic_cast<const Epetra_MpiComm *>(&vector->Map().Comm());
      Assert(mpi_comm != nullptr, ExcInternalError());
      unsigned int num_nonzero = Utilities::MPI::sum(flag, mpi_comm->Comm());

      return num_nonzero == 0;
    }



    double
    Vector::mean_value() const
    {
      double mean_value(0.);

      int ierr = vector->MeanValue(&mean_value);
      Assert(ierr == 0, ExcTrilinosError(ierr));
      (void)ierr;

      return mean_value;
    }



    double
    Vector::l1_norm() const
    {
      double norm(0.);
      int    ierr = vector->Norm1(&norm);
      Assert(ierr == 0, ExcTrilinosError(ierr));
      (void)ierr;

      return norm;
    }



    double
    Vector::l2_norm() const
    {
      double norm(0.);
      int    ierr = vector->Norm2(&norm);
      Assert(ierr == 0, ExcTrilinosError(ierr));
      (void)ierr;

      return norm;
    }



    double
    Vector::linfty_norm() const
    {
      double norm(0.);
      int    ierr = vector->NormInf(&norm);
      Assert(ierr == 0, ExcTrilinosError(ierr));
      (void)ierr;

      return norm;
    }



    double
    Vector::add_and_dot(const double                     a,
                        const VectorSpaceVector<double> &V,
                        const VectorSpaceVector<double> &W)
    {
      this->add(a, V);

      return *this * W;
    }



    Vector::size_type
    Vector::size() const
    {
#    ifndef DEAL_II_WITH_64BIT_INDICES
      return vector->GlobalLength();
#    else
      return vector->GlobalLength64();
#    endif
    }



    Vector::size_type
    Vector::locally_owned_size() const
    {
      return vector->MyLength();
    }



    MPI_Comm
    Vector::get_mpi_communicator() const
    {
      const Epetra_MpiComm *epetra_comm =
        dynamic_cast<const Epetra_MpiComm *>(&(vector->Comm()));
      Assert(epetra_comm != nullptr, ExcInternalError());
      return epetra_comm->GetMpiComm();
    }



    ::dealii::IndexSet
    Vector::locally_owned_elements() const
    {
      IndexSet is(size());

      // easy case: local range is contiguous
      if (vector->Map().LinearMap())
        {
#    ifndef DEAL_II_WITH_64BIT_INDICES
          is.add_range(vector->Map().MinMyGID(), vector->Map().MaxMyGID() + 1);
#    else
          is.add_range(vector->Map().MinMyGID64(),
                       vector->Map().MaxMyGID64() + 1);
#    endif
        }
      else if (vector->Map().NumMyElements() > 0)
        {
          const size_type n_indices = vector->Map().NumMyElements();
#    ifndef DEAL_II_WITH_64BIT_INDICES
          unsigned int *vector_indices =
            reinterpret_cast<unsigned int *>(vector->Map().MyGlobalElements());
#    else
          size_type *vector_indices =
            reinterpret_cast<size_type *>(vector->Map().MyGlobalElements64());
#    endif
          is.add_indices(vector_indices, vector_indices + n_indices);
        }
      is.compress();

      return is;
    }



    const Epetra_FEVector &
    Vector::trilinos_vector() const
    {
      return *vector;
    }



    Epetra_FEVector &
    Vector::trilinos_vector()
    {
      return *vector;
    }



    void
    Vector::print(std::ostream &     out,
                  const unsigned int precision,
                  const bool         scientific,
                  const bool         across) const
    {
      AssertThrow(out, ExcIO());
      boost::io::ios_flags_saver restore_flags(out);

      // Get a representation of the vector and loop over all
      // the elements
      double *val;
      int     leading_dimension;
      int     ierr = vector->ExtractView(&val, &leading_dimension);

      Assert(ierr == 0, ExcTrilinosError(ierr));
      (void)ierr;
      out.precision(precision);
      if (scientific)
        out.setf(std::ios::scientific, std::ios::floatfield);
      else
        out.setf(std::ios::fixed, std::ios::floatfield);

      if (across)
        for (int i = 0; i < vector->MyLength(); ++i)
          out << val[i] << ' ';
      else
        for (int i = 0; i < vector->MyLength(); ++i)
          out << val[i] << std::endl;
      out << std::endl;

      // restore the representation
      // of the vector
      AssertThrow(out, ExcIO());
    }



    std::size_t
    Vector::memory_consumption() const
    {
      return sizeof(*this) +
             vector->MyLength() *
               (sizeof(double) + sizeof(TrilinosWrappers::types::int_type));
    }



    void
    Vector::create_epetra_comm_pattern(const IndexSet &source_index_set,
                                       const MPI_Comm &mpi_comm)
    {
      source_stored_elements = source_index_set;
      epetra_comm_pattern =
        std::make_shared<CommunicationPattern>(locally_owned_elements(),
                                               source_index_set,
                                               mpi_comm);
    }
  } // namespace EpetraWrappers
} // namespace LinearAlgebra

DEAL_II_NAMESPACE_CLOSE

#  endif

#endif
