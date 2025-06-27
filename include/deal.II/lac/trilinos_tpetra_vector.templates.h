// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2018 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_trilinos_tpetra_vector_templates_h
#define dealii_trilinos_tpetra_vector_templates_h

#include <deal.II/base/config.h>

#include <deal.II/base/mpi.h>

#include <deal.II/lac/trilinos_tpetra_vector.h>

#ifdef DEAL_II_TRILINOS_WITH_TPETRA

#  include <deal.II/base/index_set.h>
#  include <deal.II/base/trilinos_utilities.h>

#  include <deal.II/lac/read_write_vector.h>

#  include <boost/io/ios_state.hpp>

#  include <Teuchos_DefaultMpiComm.hpp>
#  include <Tpetra_Import_def.hpp>
#  include <Tpetra_Map_def.hpp>

#  include <memory>


DEAL_II_NAMESPACE_OPEN

namespace LinearAlgebra
{
  namespace TpetraWrappers
  {
    namespace internal
    {
      template <typename Number, typename MemorySpace>
      VectorReference<Number, MemorySpace>::operator Number() const
      {
        // Get the local index
        const TrilinosWrappers::types::int_type local_index =
          vector.vector->getMap()->getLocalElement(
            static_cast<TrilinosWrappers::types::int_type>(index));

#  if DEAL_II_TRILINOS_VERSION_GTE(14, 0, 0)
        Assert(
          (local_index != Teuchos::OrdinalTraits<int>::invalid()),
          (typename Vector<Number, MemorySpace>::ExcAccessToNonLocalElement(
            index,
            vector.vector->getMap()->getLocalNumElements(),
            vector.vector->getMap()->getMinLocalIndex(),
            vector.vector->getMap()->getMaxLocalIndex())));
#  else
        Assert(
          (local_index != Teuchos::OrdinalTraits<int>::invalid()),
          (typename Vector<Number, MemorySpace>::ExcAccessToNonLocalElement(
            index,
            vector.vector->getMap()->getNodeNumElements(),
            vector.vector->getMap()->getMinLocalIndex(),
            vector.vector->getMap()->getMaxLocalIndex())));
#  endif
        return vector.vector->getData()[local_index];
      }
    } // namespace internal



    template <typename Number, typename MemorySpace>
    Vector<Number, MemorySpace>::Vector()
      : compressed(true)
      , has_ghost(false)
      , vector(Utilities::Trilinos::internal::make_rcp<
               TpetraTypes::VectorType<Number, MemorySpace>>(
          Utilities::Trilinos::internal::make_rcp<TpetraTypes::MapType<
            MemorySpace>>(0, 0, Utilities::Trilinos::tpetra_comm_self())))
    {}



    template <typename Number, typename MemorySpace>
    Vector<Number, MemorySpace>::Vector(const Vector<Number, MemorySpace> &V)
      : compressed(V.compressed)
      , has_ghost(V.has_ghost)
      , vector(Utilities::Trilinos::internal::make_rcp<
               TpetraTypes::VectorType<Number, MemorySpace>>(*V.vector,
                                                             Teuchos::Copy))
    {
      if (!V.nonlocal_vector.is_null())
        nonlocal_vector = Utilities::Trilinos::internal::make_rcp<
          TpetraTypes::VectorType<Number, MemorySpace>>(*V.nonlocal_vector,
                                                        Teuchos::Copy);
    }



    template <typename Number, typename MemorySpace>
    Vector<Number, MemorySpace>::Vector(
      const Teuchos::RCP<TpetraTypes::VectorType<Number, MemorySpace>> V)
      : compressed(true)
      , has_ghost(V->getMap()->isOneToOne() == false)
      , vector(V)
    {}



    template <typename Number, typename MemorySpace>
    Vector<Number, MemorySpace>::Vector(const IndexSet &parallel_partitioner,
                                        const MPI_Comm  communicator)
      : compressed(true)
      , has_ghost(false)
      , vector(Utilities::Trilinos::internal::make_rcp<
               TpetraTypes::VectorType<Number, MemorySpace>>(
          parallel_partitioner.make_tpetra_map_rcp<
            TpetraTypes::NodeType<MemorySpace>>(communicator, true)))
    {}



    template <typename Number, typename MemorySpace>
    Vector<Number, MemorySpace>::Vector(const IndexSet &locally_owned_entries,
                                        const IndexSet &ghost_entries,
                                        const MPI_Comm  communicator,
                                        const bool      vector_writable)
    {
      if (!vector_writable)
        {
          IndexSet parallel_partitioner = locally_owned_entries;
          parallel_partitioner.add_indices(ghost_entries);

          vector = Utilities::Trilinos::internal::make_rcp<
            TpetraTypes::VectorType<Number, MemorySpace>>(
            parallel_partitioner
              .template make_tpetra_map_rcp<TpetraTypes::NodeType<MemorySpace>>(
                communicator, true));

          compressed = true;
        }
      else
        {
          Teuchos::RCP<TpetraTypes::MapType<MemorySpace>> map =
            locally_owned_entries
              .template make_tpetra_map_rcp<TpetraTypes::NodeType<MemorySpace>>(
                communicator, false);
          vector = Utilities::Trilinos::internal::make_rcp<
            TpetraTypes::VectorType<Number, MemorySpace>>(map);

          IndexSet nonlocal_entries(ghost_entries);
          nonlocal_entries.subtract_set(locally_owned_entries);
          nonlocal_vector = Utilities::Trilinos::internal::make_rcp<
            TpetraTypes::VectorType<Number, MemorySpace>>(
            nonlocal_entries
              .template make_tpetra_map_rcp<TpetraTypes::NodeType<MemorySpace>>(
                communicator, true));

          compressed = false;
        }

      has_ghost = (vector->getMap()->isOneToOne() == false);

      if constexpr (running_in_debug_mode())
        {
          MPI_Comm comm = Utilities::Trilinos::teuchos_comm_to_mpi_comm(
            vector->getMap()->getComm());
          const size_type n_elements_global =
            Utilities::MPI::sum(vector->getLocalLength(), comm);
          Assert(has_ghost || n_elements_global == size(), ExcInternalError());
        }
    }



    template <typename Number, typename MemorySpace>
    void
    Vector<Number, MemorySpace>::clear()
    {
      vector = Utilities::Trilinos::internal::make_rcp<
        TpetraTypes::VectorType<Number, MemorySpace>>(
        Utilities::Trilinos::internal::make_rcp<
          TpetraTypes::MapType<MemorySpace>>(
          0, 0, Utilities::Trilinos::tpetra_comm_self()));
      has_ghost  = false;
      compressed = true;
      nonlocal_vector.reset();
    }



    template <typename Number, typename MemorySpace>
    void
    Vector<Number, MemorySpace>::reinit(const IndexSet &parallel_partitioner,
                                        const MPI_Comm  communicator,
                                        const bool /*omit_zeroing_entries*/)
    {
      vector.reset();
      nonlocal_vector.reset();

      compressed = true;
      has_ghost  = false;
      vector     = Utilities::Trilinos::internal::make_rcp<
        TpetraTypes::VectorType<Number, MemorySpace>>(
        parallel_partitioner
          .template make_tpetra_map_rcp<TpetraTypes::NodeType<MemorySpace>>(
            communicator, true));
    }



    template <typename Number, typename MemorySpace>
    void
    Vector<Number, MemorySpace>::reinit(const IndexSet &locally_owned_entries,
                                        const IndexSet &ghost_entries,
                                        const MPI_Comm  communicator,
                                        const bool      vector_writable)
    {
      // release memory before reallocation
      nonlocal_vector.reset();

      if (!vector_writable)
        {
          vector.reset();

          IndexSet parallel_partitioner = locally_owned_entries;
          parallel_partitioner.add_indices(ghost_entries);

          vector = Utilities::Trilinos::internal::make_rcp<
            TpetraTypes::VectorType<Number, MemorySpace>>(
            parallel_partitioner
              .template make_tpetra_map_rcp<TpetraTypes::NodeType<MemorySpace>>(
                communicator, true));

          compressed = true;
        }
      else
        {
          Teuchos::RCP<TpetraTypes::MapType<MemorySpace>> map =
            locally_owned_entries
              .template make_tpetra_map_rcp<TpetraTypes::NodeType<MemorySpace>>(
                communicator, false);

          if (!vector->getMap()->isSameAs(*map))
            {
              vector.reset();
              vector = Utilities::Trilinos::internal::make_rcp<
                TpetraTypes::VectorType<Number, MemorySpace>>(map);
            }
          else
            vector->putScalar(0);

          IndexSet nonlocal_entries(ghost_entries);
          nonlocal_entries.subtract_set(locally_owned_entries);

          nonlocal_vector = Utilities::Trilinos::internal::make_rcp<
            TpetraTypes::VectorType<Number, MemorySpace>>(
            nonlocal_entries
              .template make_tpetra_map_rcp<TpetraTypes::NodeType<MemorySpace>>(
                communicator, true));

          compressed = false;
        }

      has_ghost = (vector->getMap()->isOneToOne() == false);
    }



    template <typename Number, typename MemorySpace>
    void
    Vector<Number, MemorySpace>::reinit(const Vector<Number, MemorySpace> &V,
                                        const bool omit_zeroing_entries)
    {
      reinit(V.locally_owned_elements(),
             V.get_mpi_communicator(),
             omit_zeroing_entries);
    }



    template <typename Number, typename MemorySpace>
    void
    Vector<Number, MemorySpace>::extract_subvector_to(
      const ArrayView<const types::global_dof_index> &indices,
      const ArrayView<Number>                        &elements) const
    {
      AssertDimension(indices.size(), elements.size());

      auto vector_2d = vector->template getLocalView<Kokkos::HostSpace>(
        Tpetra::Access::ReadOnly);
      auto vector_1d = Kokkos::subview(vector_2d, Kokkos::ALL(), 0);

      for (unsigned int i = 0; i < indices.size(); ++i)
        {
          AssertIndexRange(indices[i], size());
          const size_type                   row = indices[i];
          TrilinosWrappers::types::int_type local_row =
            vector->getMap()->getLocalElement(row);


#  if DEAL_II_TRILINOS_VERSION_GTE(14, 0, 0)
          Assert(
            local_row != Teuchos::OrdinalTraits<int>::invalid(),
            ExcAccessToNonLocalElement(row,
                                       vector->getMap()->getLocalNumElements(),
                                       vector->getMap()->getMinLocalIndex(),
                                       vector->getMap()->getMaxLocalIndex()));
#  else
          Assert(
            local_row != Teuchos::OrdinalTraits<int>::invalid(),
            ExcAccessToNonLocalElement(row,
                                       vector->getMap()->getNodeNumElements(),
                                       vector->getMap()->getMinLocalIndex(),
                                       vector->getMap()->getMaxLocalIndex()));

#  endif

          if (local_row != Teuchos::OrdinalTraits<int>::invalid())
            elements[i] = vector_1d(local_row);
        }
    }



    template <typename Number, typename MemorySpace>
    Vector<Number, MemorySpace> &
    Vector<Number, MemorySpace>::operator=(const Vector<Number, MemorySpace> &V)
    {
      // Distinguish three cases:
      //  - First case: both vectors have the same layout.
      //  - Second case: both vectors have the same size but different layout.
      //  - Third case: the vectors have different size.

      AssertThrow(V.compressed,
                  ExcMessage("Cannot copy-assign from a vector that has not "
                             "been compressed. Please call compress() on the "
                             "source vector before using operator=()."));

      if (vector->getMap()->isSameAs(*V.vector->getMap()))
        {
          // Create a read-only Kokkos view from the source vector
          auto source_vector_2d =
            V.vector->template getLocalView<Kokkos::HostSpace>(
              Tpetra::Access::ReadOnly);

          auto source_vector_1d =
            Kokkos::subview(source_vector_2d, Kokkos::ALL(), 0);

          // Create a read/write Kokkos view from the target vector
          auto target_vector_2d =
            vector->template getLocalView<Kokkos::HostSpace>(
              Tpetra::Access::ReadWrite);
          auto target_vector_1d =
            Kokkos::subview(target_vector_2d, Kokkos::ALL(), 0);

          // Copy the data
          Kokkos::deep_copy(target_vector_1d, source_vector_1d);
        }
      else if (size() == V.size())
        {
          // We expect that at least one vector has a one-to-one map, otherwise
          // we can neither call Import nor Export.
          if (V.vector->getMap()->isOneToOne())
            {
              Teuchos::RCP<const TpetraTypes::ImportType<MemorySpace>>
                importer =
                  Tpetra::createImport(V.vector->getMap(), vector->getMap());

              // Since we are distributing the vector from a one-to-one map
              // we can always use the VectorOperation::insert / Tpetra::INSERT
              // here.
              vector->doImport(*V.vector, *importer, Tpetra::INSERT);
            }
          else if (vector->getMap()->isOneToOne())
            {
              Teuchos::RCP<const TpetraTypes::ExportType<MemorySpace>>
                exporter =
                  Tpetra::createExport(V.vector->getMap(), vector->getMap());

              vector->doExport(*V.vector, *exporter, Tpetra::INSERT);
            }
          else
            {
              Assert(false,
                     ExcMessage(
                       "You are trying to map one vector distributed "
                       "between processors, where some elements belong "
                       "to multiple processors, onto another distribution "
                       "pattern, where some elements belong to multiple "
                       "processors. It is unclear how to deal with elements "
                       "in the vector belonging to multiple processors. "
                       "Therefore, compress() must be called on this "
                       "vector first."));
            }
        }
      else
        {
          vector.reset();
          vector = Utilities::Trilinos::internal::make_rcp<
            TpetraTypes::VectorType<Number, MemorySpace>>(*V.vector,
                                                          Teuchos::Copy);

          compressed             = V.compressed;
          has_ghost              = V.has_ghost;
          source_stored_elements = V.source_stored_elements;
          tpetra_comm_pattern    = V.tpetra_comm_pattern;
        }

      return *this;
    }



    template <typename Number, typename MemorySpace>
    template <typename OtherNumber>
    Vector<Number, MemorySpace> &
    Vector<Number, MemorySpace>::operator=(const dealii::Vector<OtherNumber> &V)
    {
      static_assert(
        std::is_same<Number, OtherNumber>::value,
        "TpetraWrappers::Vector and dealii::Vector must use the same number type here.");

      vector.reset();
      nonlocal_vector.reset();

      Teuchos::Array<OtherNumber> vector_data(V.begin(), V.end());
      vector = Utilities::Trilinos::internal::make_rcp<
        TpetraTypes::VectorType<Number, MemorySpace>>(
        V.locally_owned_elements()
          .template make_tpetra_map_rcp<TpetraTypes::NodeType<MemorySpace>>(),
        vector_data);

      has_ghost  = false;
      compressed = true;

      return *this;
    }



    template <typename Number, typename MemorySpace>
    Vector<Number, MemorySpace> &
    Vector<Number, MemorySpace>::operator=(const Number s)
    {
      Assert(s == Number(0.0),
             ExcMessage("Only 0 can be assigned to a vector."));

      // As checked above, we are only allowed to use d==0.0, so pass
      // a constant zero (instead of a run-time value 'd' that *happens* to
      // have a zero value) to the underlying class in hopes that the compiler
      // can optimize this somehow.
      vector->putScalar(/*s=*/0.0);

      if (!nonlocal_vector.is_null())
        nonlocal_vector->putScalar(/*s=*/0.0);

      return *this;
    }



    template <typename Number, typename MemorySpace>
    void
    Vector<Number, MemorySpace>::import_elements(
      const ReadWriteVector<Number> &V,
      VectorOperation::values        operation,
      const Teuchos::RCP<const Utilities::MPI::CommunicationPatternBase>
        &communication_pattern)
    {
      // If no communication pattern is given, create one. Otherwise, use the
      // one given.
      if (communication_pattern.is_null())
        {
          // The first time import is called, a communication pattern is
          // created. Check if the communication pattern already exists and if
          // it can be reused.
          if ((source_stored_elements.size() !=
               V.get_stored_elements().size()) ||
              (source_stored_elements != V.get_stored_elements()))
            {
              create_tpetra_comm_pattern(
                V.get_stored_elements(),
                Utilities::Trilinos::teuchos_comm_to_mpi_comm(
                  vector->getMap()->getComm()));
            }
        }
      else
        {
          tpetra_comm_pattern = Teuchos::rcp_dynamic_cast<
            const TpetraWrappers::CommunicationPattern<MemorySpace>>(
            communication_pattern);

          AssertThrow(
            !tpetra_comm_pattern.is_null(),
            ExcMessage("The communication pattern is not of type "
                       "LinearAlgebra::TpetraWrappers::CommunicationPattern."));
        }

      Teuchos::RCP<const TpetraTypes::ExportType<MemorySpace>> tpetra_export =
        tpetra_comm_pattern->get_tpetra_export_rcp();

      TpetraTypes::VectorType<Number, MemorySpace> source_vector(
        tpetra_export->getSourceMap());

      {
        auto x_2d = source_vector.template getLocalView<Kokkos::HostSpace>(
          Tpetra::Access::ReadWrite);
        auto x_1d = Kokkos::subview(x_2d, Kokkos::ALL(), 0);

        const size_t localLength = source_vector.getLocalLength();
        auto         values_it   = V.begin();
        for (size_t k = 0; k < localLength; ++k)
          x_1d(k) = *values_it++;
      }
      Tpetra::CombineMode tpetra_operation = Tpetra::ZERO;
      if (operation == VectorOperation::insert)
        tpetra_operation = Tpetra::INSERT;
      else if (operation == VectorOperation::add)
        tpetra_operation = Tpetra::ADD;
      else
        DEAL_II_NOT_IMPLEMENTED();

      vector->doExport(source_vector, *tpetra_export, tpetra_operation);
    }



    template <typename Number, typename MemorySpace>
    void
    Vector<Number, MemorySpace>::import_elements(
      const ReadWriteVector<Number> &V,
      VectorOperation::values        operation,
      const std::shared_ptr<const Utilities::MPI::CommunicationPatternBase> &)
    {
      import_elements(V, operation);
    }



    template <typename Number, typename MemorySpace>
    void
    Vector<Number, MemorySpace>::import_elements(
      const ReadWriteVector<Number> &V,
      VectorOperation::values        operation)
    {
      // Create an empty CommunicationPattern
      const Teuchos::RCP<const Utilities::MPI::CommunicationPatternBase>
        communication_pattern_empty;

      import_elements(V, operation, communication_pattern_empty);
    }



    template <typename Number, typename MemorySpace>
    Vector<Number, MemorySpace> &
    Vector<Number, MemorySpace>::operator*=(const Number factor)
    {
      AssertIsFinite(factor);
      vector->scale(factor);

      return *this;
    }



    template <typename Number, typename MemorySpace>
    Vector<Number, MemorySpace> &
    Vector<Number, MemorySpace>::operator/=(const Number factor)
    {
      AssertIsFinite(factor);
      Assert(factor != Number(0.), ExcZero());
      *this *= Number(1.) / factor;

      return *this;
    }



    template <typename Number, typename MemorySpace>
    Vector<Number, MemorySpace> &
    Vector<Number, MemorySpace>::operator+=(
      const Vector<Number, MemorySpace> &V)
    {
      // If the maps are the same we can update right away.
      if (vector->getMap()->isSameAs(*(V.trilinos_vector().getMap())))
        {
          vector->update(1., V.trilinos_vector(), 1.);
        }
      else
        {
          Assert(this->size() == V.size(),
                 ExcDimensionMismatch(this->size(), V.size()));

          // TODO: Tpetra doesn't have a combine mode that also updates local
          // elements, maybe there is a better workaround.
          Tpetra::Vector<Number,
                         int,
                         types::signed_global_dof_index,
                         TpetraTypes::NodeType<MemorySpace>>
                                               dummy(vector->getMap(), false);
          TpetraTypes::ImportType<MemorySpace> data_exchange(
            V.trilinos_vector().getMap(), dummy.getMap());
          dummy.doImport(V.trilinos_vector(), data_exchange, Tpetra::INSERT);

          vector->update(1.0, dummy, 1.0);
        }

      return *this;
    }



    template <typename Number, typename MemorySpace>
    Vector<Number, MemorySpace> &
    Vector<Number, MemorySpace>::operator-=(
      const Vector<Number, MemorySpace> &V)
    {
      this->add(-1., V);

      return *this;
    }



    template <typename Number, typename MemorySpace>
    Number
    Vector<Number, MemorySpace>::operator*(
      const Vector<Number, MemorySpace> &V) const
    {
      Assert(this->size() == V.size(),
             ExcDimensionMismatch(this->size(), V.size()));
      Assert(vector->getMap()->isSameAs(*V.trilinos_vector().getMap()),
             ExcDifferentParallelPartitioning());
      Assert(!has_ghost_elements(), ExcGhostsPresent());

      return vector->dot(V.trilinos_vector());
    }



    template <typename Number, typename MemorySpace>
    Number
    Vector<Number, MemorySpace>::operator()(const size_type index) const
    {
      // Get the local index
      const TrilinosWrappers::types::int_type local_index =
        vector->getMap()->getLocalElement(
          static_cast<TrilinosWrappers::types::int_type>(index));

      Number value = 0.0;

      // If the element is not present on the current processor, we can't
      // continue. This is the main difference to the el() function.
      if (local_index == Teuchos::OrdinalTraits<int>::invalid())
        {
#  if DEAL_II_TRILINOS_VERSION_GTE(14, 0, 0)
          Assert(
            false,
            ExcAccessToNonLocalElement(index,
                                       vector->getMap()->getLocalNumElements(),
                                       vector->getMap()->getMinLocalIndex(),
                                       vector->getMap()->getMaxLocalIndex()));
#  else
          Assert(
            false,
            ExcAccessToNonLocalElement(index,
                                       vector->getMap()->getNodeNumElements(),
                                       vector->getMap()->getMinLocalIndex(),
                                       vector->getMap()->getMaxLocalIndex()));
#  endif
        }
      else
        value = vector->getData()[local_index];

      return value;
    }



    template <typename Number, typename MemorySpace>
    void
    Vector<Number, MemorySpace>::add(const Number a)
    {
      AssertIsFinite(a);

      // if we have ghost values, do not allow
      // writing to this vector at all.
      Assert(!has_ghost_elements(), ExcGhostsPresent());

      auto vector_2d = vector->template getLocalView<Kokkos::HostSpace>(
        Tpetra::Access::ReadWrite);
      auto vector_1d = Kokkos::subview(vector_2d, Kokkos::ALL(), 0);

      const size_t localLength = vector->getLocalLength();
      for (size_t k = 0; k < localLength; ++k)
        {
          vector_1d(k) += a;
        }
    }



    template <typename Number, typename MemorySpace>
    void
    Vector<Number, MemorySpace>::add(const Number                       a,
                                     const Vector<Number, MemorySpace> &V)
    {
      AssertIsFinite(a);

      Assert(vector->getMap()->isSameAs(*(V.trilinos_vector().getMap())),
             ExcDifferentParallelPartitioning());

      // if we have ghost values, do not allow
      // writing to this vector at all.
      Assert(!has_ghost_elements(), ExcGhostsPresent());

      vector->update(a, V.trilinos_vector(), 1.);
    }



    template <typename Number, typename MemorySpace>
    void
    Vector<Number, MemorySpace>::add(const Number                       a,
                                     const Vector<Number, MemorySpace> &V,
                                     const Number                       b,
                                     const Vector<Number, MemorySpace> &W)
    {
      AssertIsFinite(a);
      AssertIsFinite(b);

      Assert(vector->getMap()->isSameAs(*(V.trilinos_vector().getMap())),
             ExcDifferentParallelPartitioning());
      Assert(vector->getMap()->isSameAs(*(W.trilinos_vector().getMap())),
             ExcDifferentParallelPartitioning());

      // if we have ghost values, do not allow
      // writing to this vector at all.
      Assert(!has_ghost_elements(), ExcGhostsPresent());

      vector->update(a, V.trilinos_vector(), b, W.trilinos_vector(), 1.);
    }



    template <typename Number, typename MemorySpace>
    void
    Vector<Number, MemorySpace>::sadd(const Number                       s,
                                      const Number                       a,
                                      const Vector<Number, MemorySpace> &V)
    {
      AssertIsFinite(s);
      AssertIsFinite(a);

      // if we have ghost values, do not allow
      // writing to this vector at all.
      Assert(!has_ghost_elements(), ExcGhostsPresent());

      *this *= s;

      Vector<Number, MemorySpace> tmp(V);
      tmp *= a;
      *this += tmp;
    }



    template <typename Number, typename MemorySpace>
    void
    Vector<Number, MemorySpace>::scale(
      const Vector<Number, MemorySpace> &scaling_factors)
    {
      Assert(vector->getMap()->isSameAs(
               *(scaling_factors.trilinos_vector().getMap())),
             ExcDifferentParallelPartitioning());

      // if we have ghost values, do not allow
      // writing to this vector at all.
      Assert(!has_ghost_elements(), ExcGhostsPresent());

      vector->elementWiseMultiply(1., *scaling_factors.vector, *vector, 0.);
    }



    template <typename Number, typename MemorySpace>
    void
    Vector<Number, MemorySpace>::equ(const Number                       a,
                                     const Vector<Number, MemorySpace> &V)
    {
      AssertIsFinite(a);

      // if we have ghost values, do not allow
      // writing to this vector at all.
      Assert(!has_ghost_elements(), ExcGhostsPresent());

      // If we don't have the same map, copy.
      if (vector->getMap()->isSameAs(*V.trilinos_vector().getMap()) == false)
        this->sadd(0., a, V);
      else
        {
          // Otherwise, just update
          vector->update(a, V.trilinos_vector(), 0.);
        }
    }



    template <typename Number, typename MemorySpace>
    bool
    Vector<Number, MemorySpace>::all_zero() const
    {
      // get a representation of the vector and
      // loop over all the elements
      Teuchos::ArrayRCP<const Number> data       = vector->getData();
      const size_type                 n_elements = vector->getLocalLength();
      unsigned int                    flag       = 0;
      for (size_type i = 0; i < n_elements; ++i)
        {
          if (data[i] != Number(0))
            {
              flag = 1;
              break;
            }
        }

      // Check that the vector is zero on _all_ processors.
      unsigned int num_nonzero =
        Utilities::MPI::sum(flag,
                            Utilities::Trilinos::teuchos_comm_to_mpi_comm(
                              vector->getMap()->getComm()));

      return num_nonzero == 0;
    }



    template <typename Number, typename MemorySpace>
    bool
    Vector<Number, MemorySpace>::is_non_negative() const
    {
      if constexpr (!std::is_same_v<Number, std::complex<double>> &&
                    !std::is_same_v<Number, std::complex<float>>)
        {
          // get a representation of the vector and
          // loop over all the elements
          Teuchos::ArrayRCP<const Number> data       = vector->getData();
          const size_type                 n_elements = vector->getLocalLength();
          unsigned int                    flag       = 0;
          for (size_type i = 0; i < n_elements; ++i)
            {
              if (data[i] < Number(0))
                {
                  flag = 1;
                  break;
                }
            }

          // Check that the vector is non-negative on _all_ processors.
          unsigned int num_negative =
            Utilities::MPI::sum(flag,
                                Utilities::Trilinos::teuchos_comm_to_mpi_comm(
                                  vector->getMap()->getComm()));
          return num_negative == 0;
        }
      Assert(false,
             ExcMessage("You can't ask a complex value "
                        "whether it is non-negative."));
      return true;
    }



    template <typename Number, typename MemorySpace>
    Number
    Vector<Number, MemorySpace>::mean_value() const
    {
      Assert(!has_ghost_elements(), ExcGhostsPresent());

      return vector->meanValue();
    }



    template <typename Number, typename MemorySpace>
    typename Vector<Number, MemorySpace>::real_type
    Vector<Number, MemorySpace>::l1_norm() const
    {
      Assert(!has_ghost_elements(), ExcGhostsPresent());

      return vector->norm1();
    }



    template <typename Number, typename MemorySpace>
    typename Vector<Number, MemorySpace>::real_type
    Vector<Number, MemorySpace>::l2_norm() const
    {
      Assert(!has_ghost_elements(), ExcGhostsPresent());

      return vector->norm2();
    }



    template <typename Number, typename MemorySpace>
    typename Vector<Number, MemorySpace>::real_type
    Vector<Number, MemorySpace>::linfty_norm() const
    {
      Assert(!has_ghost_elements(), ExcGhostsPresent());

      return vector->normInf();
    }



    template <typename Number, typename MemorySpace>
    typename Vector<Number, MemorySpace>::real_type
    Vector<Number, MemorySpace>::norm_sqr() const
    {
      Vector<Number, MemorySpace>::real_type d = l2_norm();
      return d * d;
    }


    template <typename Number, typename MemorySpace>
    Number
    Vector<Number, MemorySpace>::add_and_dot(
      const Number                       a,
      const Vector<Number, MemorySpace> &V,
      const Vector<Number, MemorySpace> &W)
    {
      AssertIsFinite(a);

      this->add(a, V);

      return *this * W;
    }



    template <typename Number, typename MemorySpace>
    bool
    Vector<Number, MemorySpace>::operator==(
      const Vector<Number, MemorySpace> &v) const
    {
      Assert(size() == v.size(), ExcDimensionMismatch(size(), v.size()));

      const size_t this_local_length  = vector->getLocalLength();
      const size_t other_local_length = v.vector->getLocalLength();
      if (this_local_length != other_local_length)
        return false;

      auto this_vector_2d = vector->template getLocalView<Kokkos::HostSpace>(
        Tpetra::Access::ReadOnly);
      auto other_vector_2d = v.vector->template getLocalView<Kokkos::HostSpace>(
        Tpetra::Access::ReadOnly);

      auto this_vector_1d  = Kokkos::subview(this_vector_2d, Kokkos::ALL(), 0);
      auto other_vector_1d = Kokkos::subview(other_vector_2d, Kokkos::ALL(), 0);

      for (size_type i = 0; i < this_local_length; ++i)
        if (this_vector_1d(i) != other_vector_1d(i))
          return false;

      return true;
    }



    template <typename Number, typename MemorySpace>
    bool
    Vector<Number, MemorySpace>::operator!=(
      const Vector<Number, MemorySpace> &v) const
    {
      return (!(*this == v));
    }



    template <typename Number, typename MemorySpace>
    typename Vector<Number, MemorySpace>::size_type
    Vector<Number, MemorySpace>::size() const
    {
      return vector->getGlobalLength();
    }



    template <typename Number, typename MemorySpace>
    typename Vector<Number, MemorySpace>::size_type
    Vector<Number, MemorySpace>::locally_owned_size() const
    {
      return vector->getLocalLength();
    }



    template <typename Number, typename MemorySpace>
    std::pair<typename Vector<Number, MemorySpace>::size_type,
              typename Vector<Number, MemorySpace>::size_type>
    Vector<Number, MemorySpace>::local_range() const
    {
      const size_type begin = vector->getMap()->getMinGlobalIndex();
      const size_type end   = vector->getMap()->getMaxGlobalIndex() + 1;

      if constexpr (running_in_debug_mode())
        {
          const size_type n_local_elements =
#  if DEAL_II_TRILINOS_VERSION_GTE(14, 0, 0)
            vector->getMap()->getLocalNumElements();
#  else
            vector->getMap()->getNodeNumElements();
#  endif
          Assert(
            end - begin == n_local_elements,
            ExcMessage(
              "This function only makes sense if the elements that this "
              "vector stores on the current processor form a contiguous range. "
              "This does not appear to be the case for the current vector."));
        }

      return std::make_pair(begin, end);
    }



    template <typename Number, typename MemorySpace>
    bool
    Vector<Number, MemorySpace>::in_local_range(const size_type index) const
    {
      std::pair<size_type, size_type> range = local_range();

      return ((index >= range.first) && (index < range.second));
    }



    template <typename Number, typename MemorySpace>
    MPI_Comm
    Vector<Number, MemorySpace>::get_mpi_communicator() const
    {
      return Utilities::Trilinos::teuchos_comm_to_mpi_comm(
        vector->getMap()->getComm());
    }



    template <typename Number, typename MemorySpace>
    IndexSet
    Vector<Number, MemorySpace>::locally_owned_elements() const
    {
      return IndexSet(vector->getMap());
    }



    template <typename Number, typename MemorySpace>
    void
    Vector<Number, MemorySpace>::compress(
      const VectorOperation::values operation)
    {
      Assert(has_ghost == false,
             ExcMessage(
               "Calling compress() is only useful if a vector "
               "has been written into, but this is a vector with ghost "
               "elements and consequently is read-only. It does "
               "not make sense to call compress() for such "
               "vectors."));

      if (!compressed)
        {
          Tpetra::CombineMode tpetra_operation = Tpetra::ZERO;
          if (operation == VectorOperation::insert)
            tpetra_operation = Tpetra::INSERT;
          else if (operation == VectorOperation::add)
            tpetra_operation = Tpetra::ADD;
          else
            DEAL_II_NOT_IMPLEMENTED();

          Teuchos::RCP<const TpetraTypes::ExportType<MemorySpace>> exporter =
            Tpetra::createExport(nonlocal_vector->getMap(), vector->getMap());
          vector->doExport(*nonlocal_vector, *exporter, tpetra_operation);

          compressed = true;
        }
    }



    template <typename Number, typename MemorySpace>
    const TpetraTypes::VectorType<Number, MemorySpace> &
    Vector<Number, MemorySpace>::trilinos_vector() const
    {
      return *vector;
    }



    template <typename Number, typename MemorySpace>
    TpetraTypes::VectorType<Number, MemorySpace> &
    Vector<Number, MemorySpace>::trilinos_vector()
    {
      return *vector;
    }



    template <typename Number, typename MemorySpace>
    Teuchos::RCP<TpetraTypes::VectorType<Number, MemorySpace>>
    Vector<Number, MemorySpace>::trilinos_rcp()
    {
      return vector;
    }



    template <typename Number, typename MemorySpace>
    Teuchos::RCP<const TpetraTypes::VectorType<Number, MemorySpace>>
    Vector<Number, MemorySpace>::trilinos_rcp() const
    {
      return vector.getConst();
    }



    template <typename Number, typename MemorySpace>
    void
    Vector<Number, MemorySpace>::print(std::ostream      &out,
                                       const unsigned int precision,
                                       const bool         scientific,
                                       const bool         across) const
    {
      AssertThrow(out.fail() == false, ExcIO());
      boost::io::ios_flags_saver restore_flags(out);

      out.precision(precision);
      if (scientific)
        out.setf(std::ios::scientific, std::ios::floatfield);
      else
        out.setf(std::ios::fixed, std::ios::floatfield);

      auto vector_2d = vector->template getLocalView<Kokkos::HostSpace>(
        Tpetra::Access::ReadOnly);

      auto         vector_1d    = Kokkos::subview(vector_2d, Kokkos::ALL(), 0);
      const size_t local_length = vector->getLocalLength();

      if (size() != local_length)
        {
          out << "size:" << size() << " locally_owned_size:" << local_length
              << " :" << std::endl;
          for (size_type i = 0; i < local_length; ++i)
            {
              const TrilinosWrappers::types::int_type global_row =
                vector->getMap()->getGlobalElement(i);
              out << "[" << global_row << "]: " << vector_1d(i) << std::endl;
            }
        }
      else
        {
          if (across)
            for (unsigned int i = 0; i < local_length; ++i)
              out << vector_1d(i) << ' ';
          else
            for (unsigned int i = 0; i < local_length; ++i)
              out << vector_1d(i) << std::endl;
          out << std::endl;
        }

      // restore the representation
      // of the vector
      AssertThrow(out.fail() == false, ExcIO());
    }



    template <typename Number, typename MemorySpace>
    MPI_Comm
    Vector<Number, MemorySpace>::mpi_comm() const
    {
      return Utilities::Trilinos::teuchos_comm_to_mpi_comm(
        vector->getMap()->getComm());
    }


    template <typename Number, typename MemorySpace>
    std::size_t
    Vector<Number, MemorySpace>::memory_consumption() const
    {
      return sizeof(*this) +
             vector->getLocalLength() *
               (sizeof(Number) + sizeof(TrilinosWrappers::types::int_type));
    }



    template <typename Number, typename MemorySpace>
    void
    Vector<Number, MemorySpace>::create_tpetra_comm_pattern(
      const IndexSet &source_index_set,
      const MPI_Comm  mpi_comm)
    {
      source_stored_elements = source_index_set;

      tpetra_comm_pattern =
        Teuchos::rcp(new TpetraWrappers::CommunicationPattern<MemorySpace>(
          locally_owned_elements(), source_index_set, mpi_comm));
    }
  } // namespace TpetraWrappers
} // namespace LinearAlgebra

DEAL_II_NAMESPACE_CLOSE

#else

// Make sure the scripts that create the C++20 module input files have
// something to latch on if the preprocessor #ifdef above would
// otherwise lead to an empty content of the file.
DEAL_II_NAMESPACE_OPEN
DEAL_II_NAMESPACE_CLOSE

#endif

#endif
