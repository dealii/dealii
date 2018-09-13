// ---------------------------------------------------------------------
//
// Copyright (C) 2015 - 2018 by the deal.II authors
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

#include <deal.II/base/std_cxx14/memory.h>

#include <deal.II/lac/trilinos_tpetra_vector.h>

#ifdef DEAL_II_WITH_TRILINOS

#  ifdef DEAL_II_WITH_MPI

#    include <deal.II/base/index_set.h>

#    include <deal.II/lac/read_write_vector.h>

#    include <boost/io/ios_state.hpp>

#    include <Teuchos_DefaultMpiComm.hpp>
#    include <Tpetra_Import_decl.hpp>
#    include <Tpetra_Map_decl.hpp>

#    include <memory>


DEAL_II_NAMESPACE_OPEN

namespace LinearAlgebra
{
  namespace TpetraWrappers
  {
    Vector::Vector()
      : vector(new Tpetra::Vector<>(Teuchos::RCP<Tpetra::Map<>>(
          new Tpetra::Map<>(0, 0, Utilities::Trilinos::tpetra_comm_self()))))
    {}



    Vector::Vector(const Vector &V)
      : Subscriptor()
      , vector(new Tpetra::Vector<>(V.trilinos_vector(), Teuchos::Copy))
    {}



    Vector::Vector(const IndexSet &parallel_partitioner,
                   const MPI_Comm &communicator)
      : vector(new Tpetra::Vector<>(Teuchos::rcp(new Tpetra::Map<>(
          parallel_partitioner.make_tpetra_map(communicator, false)))))
    {}



    void
    Vector::reinit(const IndexSet &parallel_partitioner,
                   const MPI_Comm &communicator,
                   const bool      omit_zeroing_entries)
    {
      Tpetra::Map<> input_map =
        parallel_partitioner.make_tpetra_map(communicator, false);
      if (vector->getMap()->isSameAs(input_map) == false)
        vector = std_cxx14::make_unique<Tpetra::Vector<>>(
          Teuchos::rcp(new Tpetra::Map<>(input_map)));
      else if (omit_zeroing_entries == false)
        {
          vector->putScalar(0.);
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
      if (vector->getMap()->isSameAs(*(V.trilinos_vector().getMap())))
        *vector = V.trilinos_vector();
      else
        {
          if (size() == V.size())
            {
              Tpetra::Import<> data_exchange(vector->getMap(),
                                             V.trilinos_vector().getMap());

              vector->doImport(V.trilinos_vector(),
                               data_exchange,
                               Tpetra::REPLACE);
            }
          else
            vector =
              std_cxx14::make_unique<Tpetra::Vector<>>(V.trilinos_vector());
        }

      return *this;
    }



    Vector &
    Vector::operator=(const double s)
    {
      Assert(s == 0., ExcMessage("Only 0 can be assigned to a vector."));

      vector->putScalar(s);

      return *this;
    }



    void
    Vector::import(
      const ReadWriteVector<double> &                 V,
      VectorOperation::values                         operation,
      std::shared_ptr<const CommunicationPatternBase> communication_pattern)
    {
      // If no communication pattern is given, create one. Otherwsie, use the
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
              const Teuchos::MpiComm<int> *mpi_comm =
                dynamic_cast<const Teuchos::MpiComm<int> *>(
                  vector->getMap()->getComm().get());
              Assert(mpi_comm != nullptr, ExcInternalError());
              create_tpetra_comm_pattern(V.get_stored_elements(),
                                         *(mpi_comm->getRawMpiComm())());
            }
        }
      else
        {
          tpetra_comm_pattern = std::dynamic_pointer_cast<
            const TpetraWrappers::CommunicationPattern>(communication_pattern);
          AssertThrow(
            tpetra_comm_pattern != nullptr,
            ExcMessage(
              std::string("The communication pattern is not of type ") +
              "LinearAlgebra::TpetraWrappers::CommunicationPattern."));
        }

      Tpetra::Export<> tpetra_export(tpetra_comm_pattern->get_tpetra_export());
      Tpetra::Vector<> source_vector(tpetra_export.getSourceMap());

      source_vector.sync<Kokkos::HostSpace>();
      auto x_2d = source_vector.getLocalView<Kokkos::HostSpace>();
      auto x_1d = Kokkos::subview(x_2d, Kokkos::ALL(), 0);
      source_vector.modify<Kokkos::HostSpace>();
      const size_t localLength = source_vector.getLocalLength();
      auto         values_it   = V.begin();
      for (size_t k = 0; k < localLength; ++k)
        x_1d(k) = *values_it++;
      source_vector.sync<Tpetra::Vector<double>::device_type::memory_space>();
      if (operation == VectorOperation::insert)
        vector->doExport(source_vector, tpetra_export, Tpetra::REPLACE);
      else if (operation == VectorOperation::add)
        vector->doExport(source_vector, tpetra_export, Tpetra::ADD);
      else
        AssertThrow(false, ExcNotImplemented());
    }



    Vector &
    Vector::operator*=(const double factor)
    {
      AssertIsFinite(factor);
      vector->scale(factor);

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
      if (vector->getMap()->isSameAs(*(down_V.trilinos_vector().getMap())))
        {
          vector->update(1., down_V.trilinos_vector(), 1.);
        }
      else
        {
          Assert(this->size() == down_V.size(),
                 ExcDimensionMismatch(this->size(), down_V.size()));

          // TODO: The code doesn't work as expected so we use a workaround.
          /*Tpetra::Export<> data_exchange(vector->getMap(),
                                         down_V.trilinos_vector().getMap());
          vector->doExport(down_V.trilinos_vector(),
                           data_exchange,
                           Tpetra::ADD);*/

          Tpetra::Vector<> dummy(vector->getMap(), false);
          Tpetra::Import<> data_exchange(dummy.getMap(),
                                         down_V.trilinos_vector().getMap());

          dummy.doExport(down_V.trilinos_vector(),
                         data_exchange,
                         Tpetra::REPLACE);

          vector->update(1.0, dummy, 1.0);
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
      Assert(vector->getMap()->isSameAs(*down_V.trilinos_vector().getMap()),
             ExcDifferentParallelPartitioning());

      return vector->dot(down_V.trilinos_vector());
    }



    void
    Vector::add(const double a)
    {
      AssertIsFinite(a);

      vector->sync<Kokkos::HostSpace>();
      auto vector_2d = vector->getLocalView<Kokkos::HostSpace>();
      auto vector_1d = Kokkos::subview(vector_2d, Kokkos::ALL(), 0);
      vector->modify<Kokkos::HostSpace>();
      const size_t localLength = vector->getLocalLength();
      for (size_t k = 0; k < localLength; ++k)
        {
          vector_1d(k) += a;
        }
      vector->sync<Tpetra::Vector<double>::device_type::memory_space>();
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
      Assert(vector->getMap()->isSameAs(*(down_V.trilinos_vector().getMap())),
             ExcDifferentParallelPartitioning());

      vector->update(a, down_V.trilinos_vector(), 1.);
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
      Assert(vector->getMap()->isSameAs(*(down_V.trilinos_vector().getMap())),
             ExcDifferentParallelPartitioning());
      Assert(vector->getMap()->isSameAs(*(down_W.trilinos_vector().getMap())),
             ExcDifferentParallelPartitioning());
      AssertIsFinite(a);
      AssertIsFinite(b);

      vector->update(
        a, down_V.trilinos_vector(), b, down_W.trilinos_vector(), 1.);
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
      Assert(vector->getMap()->isSameAs(
               *(down_scaling_factors.trilinos_vector().getMap())),
             ExcDifferentParallelPartitioning());

      vector->elementWiseMultiply(1.,
                                  *down_scaling_factors.vector,
                                  *vector,
                                  0.);
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
      if (vector->getMap()->isSameAs(*down_V.trilinos_vector().getMap()) ==
          false)
        this->sadd(0., a, V);
      else
        {
          // Otherwise, just update
          vector->update(a, down_V.trilinos_vector(), 0.);
        }
    }



    bool
    Vector::all_zero() const
    {
      // get a representation of the vector and
      // loop over all the elements
      double *      start_ptr = vector->getDataNonConst().get();
      const double *ptr       = start_ptr,
                   *eptr      = start_ptr + vector->getLocalLength();
      unsigned int flag       = 0;
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
      const Teuchos::MpiComm<int> *mpi_comm =
        dynamic_cast<const Teuchos::MpiComm<int> *>(
          vector->getMap()->getComm().get());
      Assert(mpi_comm != nullptr, ExcInternalError());
      unsigned int num_nonzero =
        Utilities::MPI::sum(flag, *(mpi_comm->getRawMpiComm())());

      return num_nonzero == 0;
    }



    double
    Vector::mean_value() const
    {
      return vector->meanValue();
    }



    double
    Vector::l1_norm() const
    {
      return vector->norm1();
    }



    double
    Vector::l2_norm() const
    {
      return vector->norm2();
    }



    double
    Vector::linfty_norm() const
    {
      return vector->normInf();
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
      return vector->getGlobalLength();
    }



    MPI_Comm
    Vector::get_mpi_communicator() const
    {
      const auto tpetra_comm = dynamic_cast<const Teuchos::MpiComm<int> *>(
        vector->getMap()->getComm().get());
      Assert(tpetra_comm != nullptr, ExcInternalError());
      return *(tpetra_comm->getRawMpiComm())();
    }



    ::dealii::IndexSet
    Vector::locally_owned_elements() const
    {
      IndexSet is(size());

      // easy case: local range is contiguous
      if (vector->getMap()->isContiguous())
        {
#    ifndef DEAL_II_WITH_64BIT_INDICES
          is.add_range(vector->getMap()->getMinGlobalIndex(),
                       vector->getMap()->getMaxGlobalIndex() + 1);
#    else
          is.add_range(vector->getMap()->getMinGlobalIndex(),
                       vector->getMap()->getMaxGlobalIndex() + 1);
#    endif
        }
      else if (vector->getLocalLength() > 0)
        {
          const size_type n_indices = vector->getLocalLength();
          auto vector_indices       = vector->getMap()->getMyGlobalIndices();
          is.add_indices((unsigned int *)&vector_indices[0],
                         (unsigned int *)&vector_indices[0] + n_indices);
        }
      is.compress();

      return is;
    }



    const Tpetra::Vector<> &
    Vector::trilinos_vector() const
    {
      return *vector;
    }



    Tpetra::Vector<> &
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
      const auto val = vector->get1dView();

      out.precision(precision);
      if (scientific)
        out.setf(std::ios::scientific, std::ios::floatfield);
      else
        out.setf(std::ios::fixed, std::ios::floatfield);

      vector->sync<Kokkos::HostSpace>();
      auto         vector_2d    = vector->getLocalView<Kokkos::HostSpace>();
      auto         vector_1d    = Kokkos::subview(vector_2d, Kokkos::ALL(), 0);
      const size_t local_length = vector->getLocalLength();

      if (across)
        for (unsigned int i = 0; i < local_length; ++i)
          out << vector_1d(i) << ' ';
      else
        for (unsigned int i = 0; i < local_length; ++i)
          out << vector_1d(i) << std::endl;
      out << std::endl;

      // restore the representation
      // of the vector
      AssertThrow(out, ExcIO());
    }



    std::size_t
    Vector::memory_consumption() const
    {
      return sizeof(*this) +
             vector->getLocalLength() *
               (sizeof(double) + sizeof(TrilinosWrappers::types::int_type));
    }



    void
    Vector::create_tpetra_comm_pattern(const IndexSet &source_index_set,
                                       const MPI_Comm &mpi_comm)
    {
      source_stored_elements = source_index_set;
      tpetra_comm_pattern =
        std::make_shared<TpetraWrappers::CommunicationPattern>(
          locally_owned_elements(), source_index_set, mpi_comm);
    }
  } // namespace TpetraWrappers
} // namespace LinearAlgebra

DEAL_II_NAMESPACE_CLOSE

#  endif

#endif
