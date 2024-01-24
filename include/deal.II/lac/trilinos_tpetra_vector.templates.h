// ---------------------------------------------------------------------
//
// Copyright (C) 2018 - 2023 by the deal.II authors
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
    template <typename Number>
    Vector<Number>::Vector()
      : Subscriptor()
      , compressed(true)
      , has_ghost(false)
      , vector(Utilities::Trilinos::internal::make_rcp<VectorType>(
          Utilities::Trilinos::internal::make_rcp<MapType>(
            0,
            0,
            Utilities::Trilinos::tpetra_comm_self())))
    {}



    template <typename Number>
    Vector<Number>::Vector(const Vector<Number> &V)
      : Subscriptor()
      , compressed(V.compressed)
      , has_ghost(V.has_ghost)
      , vector(V.vector)
      , nonlocal_vector(V.nonlocal_vector)
    {}



    template <typename Number>
    Vector<Number>::Vector(const Teuchos::RCP<VectorType> V)
      : Subscriptor()
      , compressed(true)
      , has_ghost(V->getMap()->isOneToOne() == false)
      , vector(V)
    {}



    template <typename Number>
    Vector<Number>::Vector(const IndexSet &parallel_partitioner,
                           const MPI_Comm  communicator)
      : Subscriptor()
      , compressed(true)
      , has_ghost(false)
      , vector(Utilities::Trilinos::internal::make_rcp<VectorType>(
          parallel_partitioner.make_tpetra_map_rcp(communicator, true)))
    {}



    template <typename Number>
    Vector<Number>::Vector(const IndexSet &locally_owned_entries,
                           const IndexSet &ghost_entries,
                           const MPI_Comm  communicator,
                           const bool      vector_writable)
      : Subscriptor()
    {
      if (!vector_writable)
        {
          IndexSet parallel_partitioner = locally_owned_entries;
          parallel_partitioner.add_indices(ghost_entries);

          vector = Utilities::Trilinos::internal::make_rcp<VectorType>(
            parallel_partitioner.make_tpetra_map_rcp(communicator, true));

          compressed = true;
        }
      else
        {
          Teuchos::RCP<MapType> map =
            locally_owned_entries.make_tpetra_map_rcp(communicator, false);
          vector = Utilities::Trilinos::internal::make_rcp<VectorType>(map);

          IndexSet nonlocal_entries(ghost_entries);
          nonlocal_entries.subtract_set(locally_owned_entries);
          nonlocal_vector = Utilities::Trilinos::internal::make_rcp<VectorType>(
            nonlocal_entries.make_tpetra_map_rcp(communicator, true));

          compressed = false;
        }

      has_ghost = (vector->getMap()->isOneToOne() == false);

#  ifdef DEBUG
      MPI_Comm comm = Utilities::Trilinos::teuchos_comm_to_mpi_comm(
        vector->getMap()->getComm());
      const size_type n_elements_global =
        Utilities::MPI::sum(vector->getLocalLength(), comm);
      Assert(has_ghost || n_elements_global == size(), ExcInternalError());
#  endif
    }



    template <typename Number>
    void
    Vector<Number>::reinit(const IndexSet &parallel_partitioner,
                           const MPI_Comm  communicator,
                           const bool /*omit_zeroing_entries*/)
    {
      vector.reset();
      nonlocal_vector.reset();

      compressed = true;
      has_ghost  = false;
      vector     = Utilities::Trilinos::internal::make_rcp<VectorType>(
        parallel_partitioner.make_tpetra_map_rcp(communicator, true));
    }



    template <typename Number>
    void
    Vector<Number>::reinit(const IndexSet &locally_owned_entries,
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

          vector = Utilities::Trilinos::internal::make_rcp<VectorType>(
            parallel_partitioner.make_tpetra_map_rcp(communicator, true));

          compressed = true;
        }
      else
        {
          Teuchos::RCP<MapType> map =
            locally_owned_entries.make_tpetra_map_rcp(communicator, false);

          if (!vector->getMap()->isSameAs(*map))
            {
              vector.reset();
              vector = Utilities::Trilinos::internal::make_rcp<VectorType>(map);
            }
          else
            vector->putScalar(0);

          IndexSet nonlocal_entries(ghost_entries);
          nonlocal_entries.subtract_set(locally_owned_entries);

          nonlocal_vector = Utilities::Trilinos::internal::make_rcp<VectorType>(
            nonlocal_entries.make_tpetra_map_rcp(communicator, true));

          compressed = false;
        }

      has_ghost = (vector->getMap()->isOneToOne() == false);
    }



    template <typename Number>
    void
    Vector<Number>::reinit(const Vector<Number> &V,
                           const bool            omit_zeroing_entries)
    {
      reinit(V.locally_owned_elements(),
             V.get_mpi_communicator(),
             omit_zeroing_entries);
    }



    template <typename Number>
    void
    Vector<Number>::extract_subvector_to(
      const ArrayView<const types::global_dof_index> &indices,
      ArrayView<Number>                              &elements) const
    {
      AssertDimension(indices.size(), elements.size());

#  if DEAL_II_TRILINOS_VERSION_GTE(13, 2, 0)
      auto vector_2d = vector->template getLocalView<Kokkos::HostSpace>(
        Tpetra::Access::ReadOnly);
#  else
      /*
       * For Trilinos older than 13.2 we would normally have to call
       * vector.template sync<Kokkos::HostSpace>() at this place in order
       * to sync between memory spaces. This is necessary for GPU support.
       * Unfortunately, we are in a const context here and cannot call to
       * sync() (which is a non-const member function).
       *
       * Let us choose to simply ignore this problem for such an old
       * Trilinos version.
       */
      auto vector_2d = vector->template getLocalView<Kokkos::HostSpace>();
#  endif
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



    template <typename Number>
    Vector<Number> &
    Vector<Number>::operator=(const Vector<Number> &V)
    {
      // Distinguish three cases:
      //  - First case: both vectors have the same layout.
      //  - Second case: both vectors have the same size but different layout.
      //  - Third case: the vectors have different size.
      if (vector->getMap()->isSameAs(*V.vector->getMap()))
        {
          *vector = *V.vector;
        }
      else if (size() == V.size())
        {
          // We expect the origin vector to have a one-to-one map, otherwise
          // we can not call Import
          Assert(V.vector->getMap()->isOneToOne(),
                 ExcMessage(
                   "You are trying to map one vector distributed "
                   "between processors, where some elements belong "
                   "to multiple processors, onto another distribution "
                   "pattern, where some elements belong to multiple "
                   "processors. It is unclear how to deal with elements "
                   "in the vector belonging to multiple processors. "
                   "Therefore, compress() must be called on this "
                   "vector first."));

          Teuchos::RCP<const ImportType> importer =
            Tpetra::createImport(V.vector->getMap(), vector->getMap());

          // Since we are distributing the vector from a one-to-one map
          // we can always use the VectorOperation::insert / Tpetra::INSERT
          // here.
          vector->doImport(*V.vector, *importer, Tpetra::INSERT);
        }
      else
        {
          vector.reset();
          vector = Utilities::Trilinos::internal::make_rcp<VectorType>(
            V.vector->getMap());
          Tpetra::deep_copy(*vector, *V.vector);

          compressed             = V.compressed;
          has_ghost              = V.has_ghost;
          source_stored_elements = V.source_stored_elements;
          tpetra_comm_pattern    = V.tpetra_comm_pattern;
        }

      return *this;
    }



    template <typename Number>
    Vector<Number> &
    Vector<Number>::operator=(const Number s)
    {
      (void)s;
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



    template <typename Number>
    void
    Vector<Number>::import_elements(
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
            const TpetraWrappers::CommunicationPattern>(communication_pattern);

          AssertThrow(
            !tpetra_comm_pattern.is_null(),
            ExcMessage(
              std::string("The communication pattern is not of type ") +
              "LinearAlgebra::TpetraWrappers::CommunicationPattern."));
        }

      Teuchos::RCP<const ExportType> tpetra_export =
        tpetra_comm_pattern->get_tpetra_export_rcp();

      VectorType source_vector(tpetra_export->getSourceMap());

      {
#  if DEAL_II_TRILINOS_VERSION_GTE(13, 2, 0)
        auto x_2d = source_vector.template getLocalView<Kokkos::HostSpace>(
          Tpetra::Access::ReadWrite);
#  else
        source_vector.template sync<Kokkos::HostSpace>();
        auto x_2d = source_vector.template getLocalView<Kokkos::HostSpace>();
#  endif
        auto x_1d = Kokkos::subview(x_2d, Kokkos::ALL(), 0);
#  if !DEAL_II_TRILINOS_VERSION_GTE(13, 2, 0)
        source_vector.template modify<Kokkos::HostSpace>();
#  endif
        const size_t localLength = source_vector.getLocalLength();
        auto         values_it   = V.begin();
        for (size_t k = 0; k < localLength; ++k)
          x_1d(k) = *values_it++;
#  if !DEAL_II_TRILINOS_VERSION_GTE(13, 2, 0)
        source_vector.template sync<
          typename Tpetra::Vector<Number, int, types::signed_global_dof_index>::
            device_type::memory_space>();
#  endif
      }
      Tpetra::CombineMode tpetra_operation = Tpetra::ZERO;
      if (operation == VectorOperation::insert)
        tpetra_operation = Tpetra::INSERT;
      else if (operation == VectorOperation::add)
        tpetra_operation = Tpetra::ADD;
      else
        Assert(false, ExcNotImplemented());

      vector->doExport(source_vector, *tpetra_export, tpetra_operation);
    }



    template <typename Number>
    void
    Vector<Number>::import_elements(
      const ReadWriteVector<Number> &V,
      VectorOperation::values        operation,
      const std::shared_ptr<const Utilities::MPI::CommunicationPatternBase> &)
    {
      import_elements(V, operation);
    }



    template <typename Number>
    void
    Vector<Number>::import_elements(const ReadWriteVector<Number> &V,
                                    VectorOperation::values        operation)
    {
      // Create an empty CommunicationPattern
      const Teuchos::RCP<const Utilities::MPI::CommunicationPatternBase>
        communication_pattern_empty;

      import_elements(V, operation, communication_pattern_empty);
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
    Vector<Number>::operator+=(const Vector<Number> &V)
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
          Tpetra::Vector<Number, int, types::signed_global_dof_index> dummy(
            vector->getMap(), false);
          ImportType data_exchange(V.trilinos_vector().getMap(),
                                   dummy.getMap());
          dummy.doImport(V.trilinos_vector(), data_exchange, Tpetra::INSERT);

          vector->update(1.0, dummy, 1.0);
        }

      return *this;
    }



    template <typename Number>
    Vector<Number> &
    Vector<Number>::operator-=(const Vector<Number> &V)
    {
      this->add(-1., V);

      return *this;
    }



    template <typename Number>
    Number
    Vector<Number>::operator*(const Vector<Number> &V) const
    {
      Assert(this->size() == V.size(),
             ExcDimensionMismatch(this->size(), V.size()));
      Assert(vector->getMap()->isSameAs(*V.trilinos_vector().getMap()),
             ExcDifferentParallelPartitioning());
      Assert(!has_ghost_elements(), ExcGhostsPresent());

      return vector->dot(V.trilinos_vector());
    }



    template <typename Number>
    Number
    Vector<Number>::operator()(const size_type index) const
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



    template <typename Number>
    void
    Vector<Number>::add(const Number a)
    {
      AssertIsFinite(a);

      // if we have ghost values, do not allow
      // writing to this vector at all.
      Assert(!has_ghost_elements(), ExcGhostsPresent());

#  if DEAL_II_TRILINOS_VERSION_GTE(13, 2, 0)
      auto vector_2d = vector->template getLocalView<Kokkos::HostSpace>(
        Tpetra::Access::ReadWrite);
#  else
      vector->template sync<Kokkos::HostSpace>();
      auto vector_2d = vector->template getLocalView<Kokkos::HostSpace>();
#  endif
      auto vector_1d = Kokkos::subview(vector_2d, Kokkos::ALL(), 0);
#  if !DEAL_II_TRILINOS_VERSION_GTE(13, 2, 0)
      vector->template modify<Kokkos::HostSpace>();
#  endif
      const size_t localLength = vector->getLocalLength();
      for (size_t k = 0; k < localLength; ++k)
        {
          vector_1d(k) += a;
        }
#  if !DEAL_II_TRILINOS_VERSION_GTE(13, 2, 0)
      vector->template sync<
        typename Tpetra::Vector<Number, int, types::signed_global_dof_index>::
          device_type::memory_space>();
#  endif
    }



    template <typename Number>
    void
    Vector<Number>::add(const Number a, const Vector<Number> &V)
    {
      AssertIsFinite(a);

      Assert(vector->getMap()->isSameAs(*(V.trilinos_vector().getMap())),
             ExcDifferentParallelPartitioning());

      // if we have ghost values, do not allow
      // writing to this vector at all.
      Assert(!has_ghost_elements(), ExcGhostsPresent());

      vector->update(a, V.trilinos_vector(), 1.);
    }



    template <typename Number>
    void
    Vector<Number>::add(const Number          a,
                        const Vector<Number> &V,
                        const Number          b,
                        const Vector<Number> &W)
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



    template <typename Number>
    void
    Vector<Number>::sadd(const Number          s,
                         const Number          a,
                         const Vector<Number> &V)
    {
      AssertIsFinite(s);
      AssertIsFinite(a);

      // if we have ghost values, do not allow
      // writing to this vector at all.
      Assert(!has_ghost_elements(), ExcGhostsPresent());

      *this *= s;

      Vector<Number> tmp(V);
      tmp *= a;
      *this += tmp;
    }



    template <typename Number>
    void
    Vector<Number>::scale(const Vector<Number> &scaling_factors)
    {
      Assert(vector->getMap()->isSameAs(
               *(scaling_factors.trilinos_vector().getMap())),
             ExcDifferentParallelPartitioning());

      // if we have ghost values, do not allow
      // writing to this vector at all.
      Assert(!has_ghost_elements(), ExcGhostsPresent());

      vector->elementWiseMultiply(1., *scaling_factors.vector, *vector, 0.);
    }



    template <typename Number>
    void
    Vector<Number>::equ(const Number a, const Vector<Number> &V)
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



    template <typename Number>
    bool
    Vector<Number>::all_zero() const
    {
      // get a representation of the vector and
      // loop over all the elements
      Number       *start_ptr = vector->getDataNonConst().get();
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
      unsigned int num_nonzero =
        Utilities::MPI::sum(flag,
                            Utilities::Trilinos::teuchos_comm_to_mpi_comm(
                              vector->getMap()->getComm()));

      return num_nonzero == 0;
    }



    template <typename Number>
    Number
    Vector<Number>::mean_value() const
    {
      Assert(!has_ghost_elements(), ExcGhostsPresent());

      return vector->meanValue();
    }



    template <typename Number>
    typename Vector<Number>::real_type
    Vector<Number>::l1_norm() const
    {
      Assert(!has_ghost_elements(), ExcGhostsPresent());

      return vector->norm1();
    }



    template <typename Number>
    typename Vector<Number>::real_type
    Vector<Number>::l2_norm() const
    {
      Assert(!has_ghost_elements(), ExcGhostsPresent());

      return vector->norm2();
    }



    template <typename Number>
    typename Vector<Number>::real_type
    Vector<Number>::linfty_norm() const
    {
      Assert(!has_ghost_elements(), ExcGhostsPresent());

      return vector->normInf();
    }



    template <typename Number>
    Number
    Vector<Number>::add_and_dot(const Number          a,
                                const Vector<Number> &V,
                                const Vector<Number> &W)
    {
      AssertIsFinite(a);

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
    typename Vector<Number>::size_type
    Vector<Number>::locally_owned_size() const
    {
      return vector->getLocalLength();
    }



    template <typename Number>
    MPI_Comm
    Vector<Number>::get_mpi_communicator() const
    {
      return Utilities::Trilinos::teuchos_comm_to_mpi_comm(
        vector->getMap()->getComm());
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
    void
    Vector<Number>::compress(const VectorOperation::values operation)
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
            Assert(false, ExcNotImplemented());

          Teuchos::RCP<const ExportType> exporter =
            Tpetra::createExport(nonlocal_vector->getMap(), vector->getMap());
          vector->doExport(*nonlocal_vector, *exporter, tpetra_operation);

          compressed = true;
        }
    }



    template <typename Number>
    const Tpetra::Vector<Number, int, types::signed_global_dof_index> &
    Vector<Number>::trilinos_vector() const
    {
      return *vector;
    }



    template <typename Number>
    Tpetra::Vector<Number, int, types::signed_global_dof_index> &
    Vector<Number>::trilinos_vector()
    {
      return *vector;
    }



    template <typename Number>
    Teuchos::RCP<Tpetra::Vector<Number, int, types::signed_global_dof_index>>
    Vector<Number>::trilinos_rcp()
    {
      return vector;
    }



    template <typename Number>
    Teuchos::RCP<
      const Tpetra::Vector<Number, int, types::signed_global_dof_index>>
    Vector<Number>::trilinos_rcp() const
    {
      return vector.getConst();
    }



    template <typename Number>
    void
    Vector<Number>::print(std::ostream      &out,
                          const unsigned int precision,
                          const bool         scientific,
                          const bool         across) const
    {
      AssertThrow(out.fail() == false, ExcIO());
      boost::io::ios_flags_saver restore_flags(out);

      // Get a representation of the vector and loop over all
      // the elements
      const auto val = vector->get1dView();

      out.precision(precision);
      if (scientific)
        out.setf(std::ios::scientific, std::ios::floatfield);
      else
        out.setf(std::ios::fixed, std::ios::floatfield);

#  if DEAL_II_TRILINOS_VERSION_GTE(13, 2, 0)
      auto vector_2d = vector->template getLocalView<Kokkos::HostSpace>(
        Tpetra::Access::ReadOnly);
#  else
      vector->template sync<Kokkos::HostSpace>();
      auto vector_2d = vector->template getLocalView<Kokkos::HostSpace>();
#  endif
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


    template <typename Number>
    MPI_Comm
    Vector<Number>::mpi_comm() const
    {
      return Utilities::Trilinos::teuchos_comm_to_mpi_comm(
        vector->getMap()->getComm());
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
                                               const MPI_Comm  mpi_comm)
    {
      source_stored_elements = source_index_set;

      tpetra_comm_pattern =
        Teuchos::rcp(new TpetraWrappers::CommunicationPattern(
          locally_owned_elements(), source_index_set, mpi_comm));
    }
  } // namespace TpetraWrappers
} // namespace LinearAlgebra

DEAL_II_NAMESPACE_CLOSE

#endif

#endif
