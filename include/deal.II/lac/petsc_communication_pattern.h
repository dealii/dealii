// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2023 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_petsc_communication_pattern_h
#define dealii_petsc_communication_pattern_h

#include <deal.II/base/config.h>

#ifdef DEAL_II_WITH_PETSC

#  include <deal.II/base/array_view.h>
#  include <deal.II/base/communication_pattern_base.h>
#  include <deal.II/base/index_set.h>

#  include <deal.II/lac/vector_operation.h>

#  include <petscsf.h>

#  include <vector>

DEAL_II_NAMESPACE_OPEN

/**
 * @addtogroup PETScWrappers
 * @{
 */
namespace PETScWrappers
{
  /**
   * CommunicationPattern implementation based on the PetscSF object.
   * This class implements the same communication patterns of
   * Utilities::MPI::NoncontiguousPartitioner, internally using PetscSF
   * API calls.
   *
   * For additional information, see the paper @cite zhang2022petscsf
   */
  class CommunicationPattern : public Utilities::MPI::CommunicationPatternBase
  {
  public:
    /**
     * Default constructor.
     */
    CommunicationPattern();

    /**
     * Destructor.
     */
    virtual ~CommunicationPattern() override;

    virtual void
    reinit(const IndexSet &locally_owned_indices,
           const IndexSet &ghost_indices,
           const MPI_Comm  communicator) override;

    /**
     * Reinitialize the communication pattern. The argument @p indices_locally_owned
     * and @p indices_want indicates the owned and required dofs, respectively.
     * This allows the indices to not be sorted and to include entries with
     * the value numbers::invalid_dof_index which do not take part of the index
     * exchange but are present in the data vectors as padding.
     *
     * The export_to_ghost_array will populate an array with the values
     * associated to the @p indices_want only.
     *
     * This emulates the corresponding constructor in
     * Utilities::MPI::NoncontiguousPartitioner.
     */
    void
    reinit(const std::vector<types::global_dof_index> &indices_locally_owned,
           const std::vector<types::global_dof_index> &indices_want,
           const MPI_Comm                              communicator);

    /**
     * Reinitialization that takes the number of locally-owned degrees of
     * freedom @p local_size and an index set for the required ghost indices
     * @p ghost_indices.
     *
     * The local index range is translated to global indices in an ascending
     * and one-to-one fashion, i.e., the indices of process $p$ sit exactly
     * between the indices of the processes $p-1$ and $p+1$, respectively.
     *
     * The export_to_ghost_array will populate an array containing
     * values from locally-owned AND ghost indices, as for the relevant
     * set of dofs of a usual FEM simulation.
     */
    void
    reinit(const types::global_dof_index local_size,
           const IndexSet               &ghost_indices,
           const MPI_Comm                communicator);

    /**
     * Fill the vector @p ghost_array according to the precomputed communication
     * pattern with values from @p locally_owned_array.
     */
    template <typename Number>
    void
    export_to_ghosted_array(const ArrayView<const Number> &locally_owned_array,
                            const ArrayView<Number>       &ghost_array) const;

    /**
     * Start the communication round to fill the vector @p ghost_array according
     * to the precomputed communication pattern with values from
     * @p locally_owned_array. It can be overlapped with other communications.
     */
    template <typename Number>
    void
    export_to_ghosted_array_start(
      const ArrayView<const Number> &locally_owned_array,
      const ArrayView<Number>       &ghost_array) const;

    /**
     * Finish the communication round to fill the vector @p ghost_array according
     * to the precomputed communication pattern with values from
     * @p locally_owned_array. It can be overlapped with other communications.
     */
    template <typename Number>
    void
    export_to_ghosted_array_finish(
      const ArrayView<const Number> &locally_owned_array,
      const ArrayView<Number>       &ghost_array) const;

    /**
     * Modify the vector @p locally_owned_array according to the precomputed communication
     * pattern and the operation @p op with values from @p ghost_array.
     */
    template <typename Number>
    void
    import_from_ghosted_array(
      const VectorOperation::values  op,
      const ArrayView<const Number> &ghost_array,
      const ArrayView<Number>       &locally_owned_array) const;

    /**
     * Start the communication round to modify the vector @p locally_owned_array according
     * to the precomputed communication pattern and the operation @p op with
     * values from
     * @p ghost_array. It can be overlapped with other communications.
     */
    template <typename Number>
    void
    import_from_ghosted_array_start(
      const VectorOperation::values  op,
      const ArrayView<const Number> &ghost_array,
      const ArrayView<Number>       &locally_owned_array) const;

    /**
     * Finish the communication round to modify the vector @p locally_owned_array according
     * to the precomputed communication pattern and the operation @p op with
     * values from
     * @p ghost_array. It can be overlapped with other communications.
     */
    template <typename Number>
    void
    import_from_ghosted_array_finish(
      const VectorOperation::values  op,
      const ArrayView<const Number> &ghost_array,
      const ArrayView<Number>       &locally_owned_array) const;

    /**
     * Return the underlying MPI communicator.
     */
    MPI_Comm
    get_mpi_communicator() const override;

    /**
     * Conversion operator to gain access to the underlying PETSc object.
     */
    operator const PetscSF &() const
    {
      return sf;
    }

    /**
     * Reset the object.
     */
    void
    clear();

  protected:
    /**
     * A generic PetscSF object that will perform the communication.
     */
    PetscSF sf;

    /**
     * General setup
     */
    void
    do_reinit(const std::vector<PetscInt> &inidx,
              const std::vector<PetscInt> &inloc,
              const std::vector<PetscInt> &outidx,
              const std::vector<PetscInt> &outloc,
              const MPI_Comm               communicator);
  };

  /**
   * Partitioner implementation based on the PetscSF object.
   * This class implements the same communication patterns of
   * Utilities::MPI::Partitioner, internally using PetscSF
   * API calls.
   * Differently from the Utilities::MPI::Partitioner implementation, here we
   * don't need to specify the communication channel, the temporary storage,
   * and the MPI requests within export and import functions. Moreover,
   * the import API does not zero the input ghost array.
   *
   * For additional information, see the paper @cite zhang2022petscsf
   */
  class Partitioner : public Utilities::MPI::CommunicationPatternBase
  {
  public:
    /**
     * Default constructor.
     */
    Partitioner();

    /**
     * Destructor.
     */
    virtual ~Partitioner() override = default;

    /**
     * Reinitialize the partitioner. As for the Utilities::MPI::Partitioner,
     * any entry of @p ghost_indices that is also present in
     * @p locally_owned_indices is discarded.
     */
    virtual void
    reinit(const IndexSet &locally_owned_indices,
           const IndexSet &ghost_indices,
           const MPI_Comm  communicator) override;

    /**
     * Reinitialize the partitioner. As for the Utilities::MPI::Partitioner,
     * any entry of @p ghost_indices that is also present in
     * @p locally_owned_indices is discarded. This reinitialization will allow
     * to perform communications either using a ghost data array of the size
     * of @p ghost_indices or of @p larger_ghost_indices.
     */
    void
    reinit(const IndexSet &locally_owned_indices,
           const IndexSet &ghost_indices,
           const IndexSet &larger_ghost_indices,
           const MPI_Comm  communicator);

    /**
     * Return the actual number of ghost indices.
     */
    unsigned int
    n_ghost_indices() const
    {
      return n_ghost_indices_data;
    }

    /**
     * Return an IndexSet representation of the actual ghost indices.
     */
    const IndexSet &
    ghost_indices() const
    {
      return ghost_indices_data;
    }

    /**
     * Fill the vector @p ghost_array according to the precomputed communication
     * pattern with values from @p locally_owned_array.
     */
    template <typename Number>
    void
    export_to_ghosted_array(const ArrayView<const Number> &locally_owned_array,
                            const ArrayView<Number>       &ghost_array) const;

    /**
     * Start the communication round to fill the vector @p ghost_array according
     * to the precomputed communication pattern with values from
     * @p locally_owned_array. It can be overlapped with other communications.
     * Differently from the Utilities::MPI::Partitioner implementation, here we
     * don't need to specify the communication channel, the temporary storage,
     * and the MPI requests.
     */
    template <typename Number>
    void
    export_to_ghosted_array_start(
      const ArrayView<const Number> &locally_owned_array,
      const ArrayView<Number>       &ghost_array) const;

    /**
     * Finish the communication round to fill the vector @p ghost_array according
     * to the precomputed communication pattern with values from
     * @p locally_owned_array. It can be overlapped with other communications.
     * Differently from the Utilities::MPI::Partitioner implementation, here we
     * don't need to specify the communication channel, the temporary storage,
     * and the MPI requests.
     */
    template <typename Number>
    void
    export_to_ghosted_array_finish(
      const ArrayView<const Number> &locally_owned_array,
      const ArrayView<Number>       &ghost_array) const;

    /**
     * Modify the vector @p locally_owned_array according to the precomputed communication
     * pattern and the operation @p op with values from @p ghost_array.
     */
    template <typename Number>
    void
    import_from_ghosted_array(
      const VectorOperation::values  op,
      const ArrayView<const Number> &ghost_array,
      const ArrayView<Number>       &locally_owned_array) const;

    /**
     * Start the communication round to modify the vector @p locally_owned_array according
     * to the precomputed communication pattern and the operation @p op with
     * values from
     * @p ghost_array. It can be overlapped with other communications.
     * Differently from the Utilities::MPI::Partitioner implementation, here we
     * don't need to specify the communication channel, the temporary storage,
     * and the MPI requests.
     */
    template <typename Number>
    void
    import_from_ghosted_array_start(
      const VectorOperation::values  op,
      const ArrayView<const Number> &ghost_array,
      const ArrayView<Number>       &locally_owned_array) const;

    /**
     * Finish the communication round to modify the vector @p locally_owned_array according
     * to the precomputed communication pattern and the operation @p op with
     * values from
     * @p ghost_array. It can be overlapped with other communications.
     * Differently from the Utilities::MPI::Partitioner implementation, here we
     * don't need to specify the communication channel, the temporary storage,
     * and the MPI requests.
     */
    template <typename Number>
    void
    import_from_ghosted_array_finish(
      const VectorOperation::values  op,
      const ArrayView<const Number> &ghost_array,
      const ArrayView<Number>       &locally_owned_array) const;

    /**
     * Return the underlying MPI communicator.
     */
    MPI_Comm
    get_mpi_communicator() const override;

  protected:
    CommunicationPattern    ghost, larger_ghost;
    IndexSet                ghost_indices_data;
    types::global_dof_index n_ghost_indices_data;
    types::global_dof_index n_ghost_indices_larger;
  };

} // namespace PETScWrappers

/** @} */

DEAL_II_NAMESPACE_CLOSE

#else

// Make sure the scripts that create the C++20 module input files have
// something to latch on if the preprocessor #ifdef above would
// otherwise lead to an empty content of the file.
DEAL_II_NAMESPACE_OPEN
DEAL_II_NAMESPACE_CLOSE

#endif

#endif
