// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2023 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#include <deal.II/lac/petsc_communication_pattern.h>

#ifdef DEAL_II_WITH_PETSC

#  include <deal.II/base/mpi.h>

#  include <deal.II/lac/exceptions.h>

DEAL_II_NAMESPACE_OPEN
// Shorthand notation for PETSc error codes.
#  define AssertPETSc(code)                          \
    do                                               \
      {                                              \
        PetscErrorCode ierr = (code);                \
        AssertThrow(ierr == 0, ExcPETScError(ierr)); \
      }                                              \
    while (false)

namespace PETScWrappers
{
  CommunicationPattern::CommunicationPattern()
    : sf(nullptr)
  {}



  CommunicationPattern::~CommunicationPattern()
  {
    clear();
  }



  void
  CommunicationPattern::reinit(const types::global_dof_index local_size,
                               const IndexSet               &ghost_indices,
                               const MPI_Comm                communicator)
  {
    clear();
    AssertIndexRange(local_size, ghost_indices.size() + 1);
    // If the size of the index set can be converted to a PetscInt then every
    // index can also be converted
    AssertThrowIntegerConversion(ghost_indices.size(),
                                 static_cast<PetscInt>(ghost_indices.size()));

    PetscLayout layout;
    AssertPETSc(PetscLayoutCreate(communicator, &layout));
    const auto petsc_local_size = static_cast<PetscInt>(local_size);
    AssertThrowIntegerConversion(local_size, petsc_local_size);
    AssertPETSc(PetscLayoutSetLocalSize(layout, petsc_local_size));
    AssertPETSc(PetscLayoutSetUp(layout));

    PetscInt start, end;
    AssertPETSc(PetscLayoutGetRange(layout, &start, &end));

    IndexSet want;
    want.add_range(start, end);
    want.add_indices(ghost_indices);
    want.compress();

    const PetscInt *idxs;
    PetscInt        n;
    IS              is = want.make_petsc_is(communicator);
    AssertPETSc(ISGetLocalSize(is, &n));
    AssertPETSc(ISGetIndices(is, &idxs));

    AssertPETSc(PetscSFCreate(communicator, &sf));
    AssertPETSc(
      PetscSFSetGraphLayout(sf, layout, n, nullptr, PETSC_OWN_POINTER, idxs));
    AssertPETSc(PetscSFSetUp(sf));

    AssertPETSc(ISRestoreIndices(is, &idxs));
    AssertPETSc(ISDestroy(&is));
    AssertPETSc(PetscLayoutDestroy(&layout));
  }



  void
  CommunicationPattern::reinit(const IndexSet &locally_owned_indices,
                               const IndexSet &ghost_indices,
                               const MPI_Comm  communicator)
  {
    // If the sizes of the index sets can be converted to PetscInts then every
    // index can also be converted
    AssertThrowIntegerConversion(static_cast<PetscInt>(
                                   locally_owned_indices.size()),
                                 locally_owned_indices.size());
    AssertThrowIntegerConversion(static_cast<PetscInt>(ghost_indices.size()),
                                 ghost_indices.size());

    const auto            in_deal = locally_owned_indices.get_index_vector();
    std::vector<PetscInt> in_petsc(in_deal.begin(), in_deal.end());

    const auto            out_deal = ghost_indices.get_index_vector();
    std::vector<PetscInt> out_petsc(out_deal.begin(), out_deal.end());

    std::vector<PetscInt> dummy;

    this->do_reinit(in_petsc, dummy, out_petsc, dummy, communicator);
  }



  void
  CommunicationPattern::reinit(
    const std::vector<types::global_dof_index> &indices_has,
    const std::vector<types::global_dof_index> &indices_want,
    const MPI_Comm                              communicator)
  {
    // Clean vectors from numbers::invalid_dof_index (indicating padding)
    std::vector<PetscInt> indices_has_clean, indices_has_loc;
    std::vector<PetscInt> indices_want_clean, indices_want_loc;
    indices_want_clean.reserve(indices_want.size());
    indices_want_loc.reserve(indices_want.size());
    indices_has_clean.reserve(indices_has.size());
    indices_has_loc.reserve(indices_has.size());

    PetscInt loc         = 0;
    bool     has_invalid = false;
    for (const auto i : indices_has)
      {
        if (i != numbers::invalid_dof_index)
          {
            const auto petsc_i = static_cast<PetscInt>(i);
            AssertThrowIntegerConversion(i, petsc_i);
            indices_has_clean.push_back(petsc_i);
            indices_has_loc.push_back(loc);
          }
        else
          has_invalid = true;
        ++loc;
      }
    if (!has_invalid)
      indices_has_loc.clear();

    loc         = 0;
    has_invalid = false;
    for (const auto i : indices_want)
      {
        if (i != numbers::invalid_dof_index)
          {
            const auto petsc_i = static_cast<PetscInt>(i);
            AssertThrowIntegerConversion(i, petsc_i);
            indices_want_clean.push_back(petsc_i);
            indices_want_loc.push_back(loc);
          }
        else
          has_invalid = true;
        ++loc;
      }
    if (!has_invalid)
      indices_want_loc.clear();

    this->do_reinit(indices_has_clean,
                    indices_has_loc,
                    indices_want_clean,
                    indices_want_loc,
                    communicator);
  }



  void
  CommunicationPattern::do_reinit(const std::vector<PetscInt> &inidx,
                                  const std::vector<PetscInt> &inloc,
                                  const std::vector<PetscInt> &outidx,
                                  const std::vector<PetscInt> &outloc,
                                  const MPI_Comm               communicator)
  {
    clear();

    // inidx is assumed to be unstructured and non-overlapping.
    // However, it may have holes in it and not be a full cover.
    //
    // We create two PETSc SFs and compose them to get
    // the final communication pattern
    //
    //  sf1 : local distributed to tmp
    //  sf2 : tmp to local with ghosts
    //  sf(x) = sf2(sf1(x))
    PetscSF sf1, sf2;

    // First create an SF where leaves are inidx (at location inloc)
    // and roots are unique indices in contiguous way
    // Code adapted from MatZeroRowsMapLocal_Private in PETSc
    PetscInt n  = static_cast<PetscInt>(inidx.size());
    PetscInt lN = n > 0 ? *std::max_element(inidx.begin(), inidx.end()) : -1;
    PetscInt N, nl;

    Utilities::MPI::internal::all_reduce<PetscInt>(
      MPI_MAX,
      ArrayView<const PetscInt>(&lN, 1),
      communicator,
      ArrayView<PetscInt>(&N, 1));

    PetscSFNode *remotes;
    AssertPETSc(PetscMalloc1(n, &remotes));

    PetscLayout layout;
    AssertPETSc(PetscLayoutCreate(communicator, &layout));
    AssertPETSc(PetscLayoutSetSize(layout, N + 1));
    AssertPETSc(PetscLayoutSetUp(layout));
    AssertPETSc(PetscLayoutGetLocalSize(layout, &nl));

    const PetscInt *ranges;
    AssertPETSc(PetscLayoutGetRanges(layout, &ranges));

    PetscInt cnt = 0;
#  if DEAL_II_PETSC_VERSION_GTE(3, 13, 0)
    PetscMPIInt owner = 0;
#  else
    PetscInt owner = 0;
#  endif
    for (const auto idx : inidx)
      {
        // short-circuit the search if the last owner owns this index too
        if (idx < ranges[owner] || ranges[owner + 1] <= idx)
          {
            AssertPETSc(PetscLayoutFindOwner(layout, idx, &owner));
          }
        remotes[cnt].rank  = owner;
        remotes[cnt].index = idx - ranges[owner];
        ++cnt;
      }

    AssertPETSc(PetscSFCreate(communicator, &sf2));
    AssertPETSc(PetscSFSetGraph(sf2,
                                nl,
                                n,
                                const_cast<PetscInt *>(
                                  inloc.size() > 0 ? inloc.data() : nullptr),
                                PETSC_COPY_VALUES,
                                remotes,
                                PETSC_OWN_POINTER));
    AssertPETSc(PetscSFSetUp(sf2));
    // We need to invert root and leaf space to create the first SF
    AssertPETSc(PetscSFCreateInverseSF(sf2, &sf1));
    AssertPETSc(PetscSFDestroy(&sf2));

    // Now create the SF from the contiguous space to the local output space
    n = static_cast<PetscInt>(outidx.size());
    AssertPETSc(PetscSFCreate(communicator, &sf2));
    AssertPETSc(PetscSFSetGraphLayout(
      sf2,
      layout,
      n,
      const_cast<PetscInt *>(outloc.size() > 0 ? outloc.data() : nullptr),
      PETSC_COPY_VALUES,
      const_cast<PetscInt *>(n > 0 ? outidx.data() : nullptr)));
    AssertPETSc(PetscSFSetUp(sf2));

    // The final SF is the composition of the two
    AssertPETSc(PetscSFCompose(sf1, sf2, &sf));

    // Cleanup
    AssertPETSc(PetscLayoutDestroy(&layout));
    AssertPETSc(PetscSFDestroy(&sf1));
    AssertPETSc(PetscSFDestroy(&sf2));
  }



  void
  CommunicationPattern::clear()
  {
    AssertPETSc(PetscSFDestroy(&sf));
  }



  MPI_Comm
  CommunicationPattern::get_mpi_communicator() const
  {
    return PetscObjectComm(reinterpret_cast<PetscObject>(sf));
  }



  template <typename Number>
  void
  CommunicationPattern::export_to_ghosted_array_start(
    const ArrayView<const Number> &src,
    const ArrayView<Number>       &dst) const
  {
    auto datatype = Utilities::MPI::mpi_type_id_for_type<Number>;

#  if DEAL_II_PETSC_VERSION_LT(3, 15, 0)
    AssertPETSc(PetscSFBcastBegin(sf, datatype, src.data(), dst.data()));
#  else
    AssertPETSc(
      PetscSFBcastBegin(sf, datatype, src.data(), dst.data(), MPI_REPLACE));
#  endif
  }



  template <typename Number>
  void
  CommunicationPattern::export_to_ghosted_array_finish(
    const ArrayView<const Number> &src,
    const ArrayView<Number>       &dst) const
  {
    auto datatype = Utilities::MPI::mpi_type_id_for_type<Number>;

#  if DEAL_II_PETSC_VERSION_LT(3, 15, 0)
    AssertPETSc(PetscSFBcastEnd(sf, datatype, src.data(), dst.data()));
#  else
    AssertPETSc(
      PetscSFBcastEnd(sf, datatype, src.data(), dst.data(), MPI_REPLACE));
#  endif
  }



  template <typename Number>
  void
  CommunicationPattern::export_to_ghosted_array(
    const ArrayView<const Number> &src,
    const ArrayView<Number>       &dst) const
  {
    export_to_ghosted_array_start(src, dst);
    export_to_ghosted_array_finish(src, dst);
  }



  template <typename Number>
  void
  CommunicationPattern::import_from_ghosted_array_start(
    const VectorOperation::values  op,
    const ArrayView<const Number> &src,
    const ArrayView<Number>       &dst) const
  {
    MPI_Op mpiop    = (op == VectorOperation::insert) ? MPI_REPLACE : MPI_SUM;
    auto   datatype = Utilities::MPI::mpi_type_id_for_type<Number>;

    AssertPETSc(
      PetscSFReduceBegin(sf, datatype, src.data(), dst.data(), mpiop));
  }



  template <typename Number>
  void
  CommunicationPattern::import_from_ghosted_array_finish(
    const VectorOperation::values  op,
    const ArrayView<const Number> &src,
    const ArrayView<Number>       &dst) const
  {
    MPI_Op mpiop    = (op == VectorOperation::insert) ? MPI_REPLACE : MPI_SUM;
    auto   datatype = Utilities::MPI::mpi_type_id_for_type<Number>;

    AssertPETSc(PetscSFReduceEnd(sf, datatype, src.data(), dst.data(), mpiop));
  }



  template <typename Number>
  void
  CommunicationPattern::import_from_ghosted_array(
    const VectorOperation::values  op,
    const ArrayView<const Number> &src,
    const ArrayView<Number>       &dst) const
  {
    import_from_ghosted_array_start(op, src, dst);
    import_from_ghosted_array_finish(op, src, dst);
  }

#  ifndef DOXYGEN

  // Partitioner

  Partitioner::Partitioner()
    : ghost()
    , larger_ghost()
    , ghost_indices_data()
    , n_ghost_indices_data(numbers::invalid_dof_index)
    , n_ghost_indices_larger(numbers::invalid_dof_index)
  {}



  void
  Partitioner::reinit(const IndexSet &locally_owned_indices,
                      const IndexSet &ghost_indices,
                      const MPI_Comm  communicator)
  {
    ghost_indices_data = ghost_indices;
    ghost_indices_data.subtract_set(locally_owned_indices);
    ghost_indices_data.compress();

    ghost.reinit(locally_owned_indices, ghost_indices_data, communicator);
    larger_ghost.clear();

    n_ghost_indices_data   = ghost_indices_data.n_elements();
    n_ghost_indices_larger = numbers::invalid_dof_index;
  }

  void
  Partitioner::reinit(const IndexSet &locally_owned_indices,
                      const IndexSet &ghost_indices,
                      const IndexSet &larger_ghost_indices,
                      const MPI_Comm  communicator)
  {
    ghost_indices_data = ghost_indices;
    ghost_indices_data.subtract_set(locally_owned_indices);
    ghost_indices_data.compress();

    std::vector<types::global_dof_index> expanded_ghost_indices(
      larger_ghost_indices.n_elements(), numbers::invalid_dof_index);
    for (auto index : ghost_indices_data)
      {
        Assert(larger_ghost_indices.is_element(index),
               ExcMessage("The given larger ghost index set must contain "
                          "all indices in the actual index set."));
        auto tmp_index = larger_ghost_indices.index_within_set(index);
        expanded_ghost_indices[tmp_index] = index;
      }

    ghost.reinit(locally_owned_indices, ghost_indices_data, communicator);
    larger_ghost.reinit(locally_owned_indices.get_index_vector(),
                        expanded_ghost_indices,
                        communicator);
    n_ghost_indices_data   = ghost_indices_data.n_elements();
    n_ghost_indices_larger = larger_ghost_indices.n_elements();
  }

  MPI_Comm
  Partitioner::get_mpi_communicator() const
  {
    return ghost.get_mpi_communicator();
  }

  template <typename Number>
  void
  Partitioner::export_to_ghosted_array_start(const ArrayView<const Number> &src,
                                             const ArrayView<Number> &dst) const
  {
    if (dst.size() == n_ghost_indices_larger)
      {
        larger_ghost.export_to_ghosted_array_start(src, dst);
      }
    else
      {
        ghost.export_to_ghosted_array_start(src, dst);
      }
  }

  template <typename Number>
  void
  Partitioner::export_to_ghosted_array_finish(
    const ArrayView<const Number> &src,
    const ArrayView<Number>       &dst) const
  {
    if (dst.size() == n_ghost_indices_larger)
      {
        larger_ghost.export_to_ghosted_array_finish(src, dst);
      }
    else
      {
        ghost.export_to_ghosted_array_finish(src, dst);
      }
  }

  template <typename Number>
  void
  Partitioner::export_to_ghosted_array(const ArrayView<const Number> &src,
                                       const ArrayView<Number>       &dst) const
  {
    export_to_ghosted_array_start(src, dst);
    export_to_ghosted_array_finish(src, dst);
  }

  template <typename Number>
  void
  Partitioner::import_from_ghosted_array_start(
    const VectorOperation::values  op,
    const ArrayView<const Number> &src,
    const ArrayView<Number>       &dst) const
  {
    if (src.size() == n_ghost_indices_larger)
      {
        larger_ghost.import_from_ghosted_array_start(op, src, dst);
      }
    else
      {
        ghost.import_from_ghosted_array_start(op, src, dst);
      }
  }

  template <typename Number>
  void
  Partitioner::import_from_ghosted_array_finish(
    const VectorOperation::values  op,
    const ArrayView<const Number> &src,
    const ArrayView<Number>       &dst) const
  {
    if (src.size() == n_ghost_indices_larger)
      {
        larger_ghost.import_from_ghosted_array_finish(op, src, dst);
      }
    else
      {
        ghost.import_from_ghosted_array_finish(op, src, dst);
      }
  }

  template <typename Number>
  void
  Partitioner::import_from_ghosted_array(const VectorOperation::values  op,
                                         const ArrayView<const Number> &src,
                                         const ArrayView<Number> &dst) const
  {
    import_from_ghosted_array_start(op, src, dst);
    import_from_ghosted_array_finish(op, src, dst);
  }
#  endif
} // namespace PETScWrappers

// Explicit instantiations
#  include "lac/petsc_communication_pattern.inst"

DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_WITH_PETSC
