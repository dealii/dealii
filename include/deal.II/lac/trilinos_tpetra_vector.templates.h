// ---------------------------------------------------------------------
//
// Copyright (C) 2018 - 2019 by the deal.II authors
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

#ifndef dealii_trilinos_tpetra_vector_templates_h
#define dealii_trilinos_tpetra_vector_templates_h

#include <deal.II/base/config.h>

#include <deal.II/lac/trilinos_tpetra_vector.h>

#ifdef DEAL_II_TRILINOS_WITH_TPETRA

#  ifdef DEAL_II_WITH_MPI

#    include <deal.II/base/index_set.h>

#    include <deal.II/lac/read_write_vector.h>

#    include <boost/io/ios_state.hpp>

#    include <Teuchos_DefaultMpiComm.hpp>
#    include <Tpetra_Import_def.hpp>
#    include <Tpetra_Map_def.hpp>

#    include <memory>


DEAL_II_NAMESPACE_OPEN

namespace LinearAlgebra
{
  namespace TpetraWrappers
  {
    template <typename Number>
    Vector<Number>::Vector()
      : Subscriptor()
      , vector(new Tpetra::Vector<Number, int, types::global_dof_index>(
          Teuchos::RCP<Tpetra::Map<int, types::global_dof_index>>(
            new Tpetra::Map<int, types::global_dof_index>(
              0,
              0,
              Utilities::Trilinos::tpetra_comm_self()))))
    {}



    template <typename Number>
    Vector<Number>::Vector(const Vector<Number> &V)
      : Subscriptor()
      , vector(new Tpetra::Vector<Number, int, types::global_dof_index>(
          V.trilinos_vector(),
          Teuchos::Copy))
    {}



    template <typename Number>
    Vector<Number>::Vector(const IndexSet &parallel_partitioner,
                           const MPI_Comm &communicator)
      : Subscriptor()
      , vector(new Tpetra::Vector<Number, int, types::global_dof_index>(
          Teuchos::rcp(new Tpetra::Map<int, types::global_dof_index>(
            parallel_partitioner.make_tpetra_map(communicator, false)))))
    {}



    template <typename Number>
    void
    Vector<Number>::reinit(const IndexSet &parallel_partitioner,
                           const MPI_Comm &communicator,
                           const bool      omit_zeroing_entries)
    {
      Tpetra::Map<int, types::global_dof_index> input_map =
        parallel_partitioner.make_tpetra_map(communicator, false);
      if (vector->getMap()->isSameAs(input_map) == false)
        vector = std::make_unique<
          Tpetra::Vector<Number, int, types::global_dof_index>>(Teuchos::rcp(
          new Tpetra::Map<int, types::global_dof_index>(input_map)));
      else if (omit_zeroing_entries == false)
        {
          vector->putScalar(0.);
        }
    }



    template <typename Number>
    void
    Vector<Number>::reinit(const VectorSpaceVector<Number> &V,
                           const bool omit_zeroing_entries)
    {
      // Check that casting will work.
      Assert(dynamic_cast<const Vector<Number> *>(&V) != nullptr,
             ExcVectorTypeNotCompatible());

      // Downcast V. If fails, throws an exception.
      const Vector<Number> &down_V = dynamic_cast<const Vector<Number> &>(V);

      reinit(down_V.locally_owned_elements(),
             down_V.get_mpi_communicator(),
             omit_zeroing_entries);
    }



    template <typename Number>
    Vector<Number> &
    Vector<Number>::operator=(const Vector<Number> &V)
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
              Tpetra::Import<int, types::global_dof_index> data_exchange(
                vector->getMap(), V.trilinos_vector().getMap());

              vector->doImport(V.trilinos_vector(),
                               data_exchange,
                               Tpetra::REPLACE);
            }
          else
            vector = std::make_unique<
              Tpetra::Vector<Number, int, types::global_dof_index>>(
              V.trilinos_vector());
        }

      return *this;
    }



    template <typename Number>
    Vector<Number> &
    Vector<Number>::operator=(const Number s)
    {
      Assert(s == Number(0.),
             ExcMessage("Only 0 can be assigned to a vector."));

      vector->putScalar(s);

      return *this;
    }



    template <typename Number>
    void
    Vector<Number>::import(
      const ReadWriteVector<Number> &                 V,
      VectorOperation::values                         operation,
      std::shared_ptr<const CommunicationPatternBase> communication_pattern)
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

      Tpetra::Export<int, types::global_dof_index> tpetra_export(
        tpetra_comm_pattern->get_tpetra_export());
      Tpetra::Vector<Number, int, types::global_dof_index> source_vector(
        tpetra_export.getSourceMap());

      source_vector.template sync<Kokkos::HostSpace>();
      auto x_2d = source_vector.template getLocalView<Kokkos::HostSpace>();
      auto x_1d = Kokkos::subview(x_2d, Kokkos::ALL(), 0);
      source_vector.template modify<Kokkos::HostSpace>();
      const size_t localLength = source_vector.getLocalLength();
      auto         values_it   = V.begin();
      for (size_t k = 0; k < localLength; ++k)
        x_1d(k) = *values_it++;
      source_vector.template sync<
        typename Tpetra::Vector<Number, int, types::global_dof_index>::
          device_type::memory_space>();
      if (operation == VectorOperation::insert)
        vector->doExport(source_vector, tpetra_export, Tpetra::REPLACE);
      else if (operation == VectorOperation::add)
        vector->doExport(source_vector, tpetra_export, Tpetra::ADD);
      else
        AssertThrow(false, ExcNotImplemented());
    }



    template <typename Number>
    Vector<Number> &
    Vector<Number>::operator*=(const Number factor)
    {
      AssertIsFinite(factor);
      vector->scale(factor);

      return *this;
    }



    template <typename Number>
    Vector<Number> &
    Vector<Number>::operator/=(const Number factor)
    {
      AssertIsFinite(factor);
      Assert(factor != Number(0.), ExcZero());
      *this *= Number(1.) / factor;

      return *this;
    }



    template <typename Number>
    Vector<Number> &
    Vector<Number>::operator+=(const VectorSpaceVector<Number> &V)
    {
      // Check that casting will work.
      Assert(dynamic_cast<const Vector<Number> *>(&V) != nullptr,
             ExcVectorTypeNotCompatible());

      // Downcast V. If fails, throws an exception.
      const Vector<Number> &down_V = dynamic_cast<const Vector<Number> &>(V);
      // If the maps are the same we can update right away.
      if (vector->getMap()->isSameAs(*(down_V.trilinos_vector().getMap())))
        {
          vector->update(1., down_V.trilinos_vector(), 1.);
        }
      else
        {
          Assert(this->size() == down_V.size(),
                 ExcDimensionMismatch(this->size(), down_V.size()));

          // TODO: Tpetra doesn't have a combine mode that also updates local
          // elements, maybe there is a better workaround.
          Tpetra::Vector<Number, int, types::global_dof_index> dummy(
            vector->getMap(), false);
          Tpetra::Import<int, types::global_dof_index> data_exchange(
            down_V.trilinos_vector().getMap(), dummy.getMap());

          dummy.doImport(down_V.trilinos_vector(),
                         data_exchange,
                         Tpetra::INSERT);

          vector->update(1.0, dummy, 1.0);
        }

      return *this;
    }



    template <typename Number>
    Vector<Number> &
    Vector<Number>::operator-=(const VectorSpaceVector<Number> &V)
    {
      this->add(-1., V);

      return *this;
    }



    template <typename Number>
    Number Vector<Number>::operator*(const VectorSpaceVector<Number> &V) const
    {
      // Check that casting will work.
      Assert(dynamic_cast<const Vector<Number> *>(&V) != nullptr,
             ExcVectorTypeNotCompatible());

      // Downcast V. If fails, throws an exception.
      const Vector<Number> &down_V = dynamic_cast<const Vector<Number> &>(V);
      Assert(this->size() == down_V.size(),
             ExcDimensionMismatch(this->size(), down_V.size()));
      Assert(vector->getMap()->isSameAs(*down_V.trilinos_vector().getMap()),
             ExcDifferentParallelPartitioning());

      return vector->dot(down_V.trilinos_vector());
    }



    template <typename Number>
    void
    Vector<Number>::add(const Number a)
    {
      AssertIsFinite(a);

      vector->template sync<Kokkos::HostSpace>();
      auto vector_2d = vector->template getLocalView<Kokkos::HostSpace>();
      auto vector_1d = Kokkos::subview(vector_2d, Kokkos::ALL(), 0);
      vector->template modify<Kokkos::HostSpace>();
      const size_t localLength = vector->getLocalLength();
      for (size_t k = 0; k < localLength; ++k)
        {
          vector_1d(k) += a;
        }
      vector->template sync<
        typename Tpetra::Vector<Number, int, types::global_dof_index>::
          device_type::memory_space>();
    }



    template <typename Number>
    void
    Vector<Number>::add(const Number a, const VectorSpaceVector<Number> &V)
    {
      // Check that casting will work.
      Assert(dynamic_cast<const Vector<Number> *>(&V) != nullptr,
             ExcVectorTypeNotCompatible());

      // Downcast V. If fails, throws an exception.
      const Vector<Number> &down_V = dynamic_cast<const Vector<Number> &>(V);
      AssertIsFinite(a);
      Assert(vector->getMap()->isSameAs(*(down_V.trilinos_vector().getMap())),
             ExcDifferentParallelPartitioning());

      vector->update(a, down_V.trilinos_vector(), 1.);
    }



    template <typename Number>
    void
    Vector<Number>::add(const Number                     a,
                        const VectorSpaceVector<Number> &V,
                        const Number                     b,
                        const VectorSpaceVector<Number> &W)
    {
      // Check that casting will work.
      Assert(dynamic_cast<const Vector<Number> *>(&V) != nullptr,
             ExcVectorTypeNotCompatible());
      // Check that casting will work.
      Assert(dynamic_cast<const Vector<Number> *>(&W) != nullptr,
             ExcVectorTypeNotCompatible());

      // Downcast V. If fails, throws an exception.
      const Vector<Number> &down_V = dynamic_cast<const Vector<Number> &>(V);
      // Downcast W. If fails, throws an exception.
      const Vector<Number> &down_W = dynamic_cast<const Vector<Number> &>(W);
      Assert(vector->getMap()->isSameAs(*(down_V.trilinos_vector().getMap())),
             ExcDifferentParallelPartitioning());
      Assert(vector->getMap()->isSameAs(*(down_W.trilinos_vector().getMap())),
             ExcDifferentParallelPartitioning());
      AssertIsFinite(a);
      AssertIsFinite(b);

      vector->update(
        a, down_V.trilinos_vector(), b, down_W.trilinos_vector(), 1.);
    }



    template <typename Number>
    void
    Vector<Number>::sadd(const Number                     s,
                         const Number                     a,
                         const VectorSpaceVector<Number> &V)
    {
      // Check that casting will work.
      Assert(dynamic_cast<const Vector<Number> *>(&V) != nullptr,
             ExcVectorTypeNotCompatible());

      *this *= s;
      // Downcast V. It fails, throws an exception.
      const Vector<Number> &down_V = dynamic_cast<const Vector<Number> &>(V);
      Vector<Number>        tmp(down_V);
      tmp *= a;
      *this += tmp;
    }



    template <typename Number>
    void
    Vector<Number>::scale(const VectorSpaceVector<Number> &scaling_factors)
    {
      // Check that casting will work.
      Assert(dynamic_cast<const Vector<Number> *>(&scaling_factors) != nullptr,
             ExcVectorTypeNotCompatible());

      // Downcast scaling_factors. If fails, throws an exception.
      const Vector<Number> &down_scaling_factors =
        dynamic_cast<const Vector<Number> &>(scaling_factors);
      Assert(vector->getMap()->isSameAs(
               *(down_scaling_factors.trilinos_vector().getMap())),
             ExcDifferentParallelPartitioning());

      vector->elementWiseMultiply(1.,
                                  *down_scaling_factors.vector,
                                  *vector,
                                  0.);
    }



    template <typename Number>
    void
    Vector<Number>::equ(const Number a, const VectorSpaceVector<Number> &V)
    {
      // Check that casting will work.
      Assert(dynamic_cast<const Vector<Number> *>(&V) != nullptr,
             ExcVectorTypeNotCompatible());

      // Downcast V. If fails, throws an exception.
      const Vector<Number> &down_V = dynamic_cast<const Vector<Number> &>(V);
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



    template <typename Number>
    bool
    Vector<Number>::all_zero() const
    {
      // get a representation of the vector and
      // loop over all the elements
      Number *      start_ptr = vector->getDataNonConst().get();
      const Number *ptr       = start_ptr,
                   *eptr      = start_ptr + vector->getLocalLength();
      unsigned int flag       = 0;
      while (ptr != eptr)
        {
          if (*ptr != Number(0))
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



    template <typename Number>
    Number
    Vector<Number>::mean_value() const
    {
      return vector->meanValue();
    }



    template <typename Number>
    typename LinearAlgebra::VectorSpaceVector<Number>::real_type
    Vector<Number>::l1_norm() const
    {
      return vector->norm1();
    }



    template <typename Number>
    typename LinearAlgebra::VectorSpaceVector<Number>::real_type
    Vector<Number>::l2_norm() const
    {
      return vector->norm2();
    }



    template <typename Number>
    typename LinearAlgebra::VectorSpaceVector<Number>::real_type
    Vector<Number>::linfty_norm() const
    {
      return vector->normInf();
    }



    template <typename Number>
    Number
    Vector<Number>::add_and_dot(const Number                     a,
                                const VectorSpaceVector<Number> &V,
                                const VectorSpaceVector<Number> &W)
    {
      this->add(a, V);

      return *this * W;
    }



    template <typename Number>
    typename Vector<Number>::size_type
    Vector<Number>::size() const
    {
      return vector->getGlobalLength();
    }



    template <typename Number>
    MPI_Comm
    Vector<Number>::get_mpi_communicator() const
    {
      const auto tpetra_comm = dynamic_cast<const Teuchos::MpiComm<int> *>(
        vector->getMap()->getComm().get());
      Assert(tpetra_comm != nullptr, ExcInternalError());
      return *(tpetra_comm->getRawMpiComm())();
    }



    template <typename Number>
    ::dealii::IndexSet
    Vector<Number>::locally_owned_elements() const
    {
      IndexSet is(size());

      // easy case: local range is contiguous
      if (vector->getMap()->isContiguous())
        {
          is.add_range(vector->getMap()->getMinGlobalIndex(),
                       vector->getMap()->getMaxGlobalIndex() + 1);
        }
      else if (vector->getLocalLength() > 0)
        {
          const size_type n_indices = vector->getLocalLength();
          std::vector<types::global_dof_index> vector_indices;
          vector_indices.reserve(n_indices);
          for (unsigned int i = 0; i < n_indices; ++i)
            vector_indices.push_back(vector->getMap()->getGlobalElement(i));

          is.add_indices(vector_indices.data(),
                         vector_indices.data() + n_indices);
        }
      is.compress();

      return is;
    }



    template <typename Number>
    const Tpetra::Vector<Number, int, types::global_dof_index> &
    Vector<Number>::trilinos_vector() const
    {
      return *vector;
    }



    template <typename Number>
    Tpetra::Vector<Number, int, types::global_dof_index> &
    Vector<Number>::trilinos_vector()
    {
      return *vector;
    }



    template <typename Number>
    void
    Vector<Number>::print(std::ostream &     out,
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

      vector->template sync<Kokkos::HostSpace>();
      auto vector_2d = vector->template getLocalView<Kokkos::HostSpace>();
      auto vector_1d = Kokkos::subview(vector_2d, Kokkos::ALL(), 0);
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



    template <typename Number>
    std::size_t
    Vector<Number>::memory_consumption() const
    {
      return sizeof(*this) +
             vector->getLocalLength() *
               (sizeof(Number) + sizeof(TrilinosWrappers::types::int_type));
    }



    template <typename Number>
    void
    Vector<Number>::create_tpetra_comm_pattern(const IndexSet &source_index_set,
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

#endif
