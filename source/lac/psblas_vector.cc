// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2019 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


#include "deal.II/base/exception_macros.h"
#include "deal.II/base/exceptions.h"
#include <deal.II/base/memory_consumption.h>
#include <deal.II/base/mpi.h>

#include <deal.II/lac/psblas_vector.h>

#ifdef DEAL_II_WITH_PSBLAS
#  include <psb_base_cbind.h>
#  include <psb_c_base.h>
#  include <psb_c_dbase.h>

DEAL_II_NAMESPACE_OPEN

namespace PSCToolkit
{

  Vector::Vector()
    : psblas_vector(nullptr)
    , psblas_context(InitFinalize::get_psblas_context())
    , psblas_descriptor(nullptr)
    , ghosted(false)
    , state(internal::State::Default)
    , last_action(VectorOperation::unknown)
  {}


  Vector::Vector(const Vector &v)
    : Vector::Vector()
  {
    if (v.has_ghost_elements())
      reinit(v.owned_elements, v.ghost_indices, v.communicator);
    else
      reinit(v.owned_elements, v.communicator);

    this->operator=(v);
  }


  Vector::~Vector()
  {
    if (state != internal::State::Default)
      {
        int err;
        err = psb_c_dgefree(psblas_vector, psblas_descriptor.get());
        Assert(err == 0, internal::ExcFreePSBLASVector(err));
      }
  }



  Vector::Vector(const IndexSet &local_partitioning, const MPI_Comm comm)
  {
    communicator  = comm;
    psblas_vector = nullptr;
    ghosted       = false;
    // forward the call to reinit function below.
    reinit(local_partitioning, communicator);
  }



  Vector::Vector(const IndexSet &local_partitioning,
                 const IndexSet &ghost_indices,
                 const MPI_Comm  comm)
  {
    communicator  = comm;
    psblas_vector = nullptr;
    ghosted       = true;
    // forward the call to reinit function taking ghost indices.
    reinit(local_partitioning, ghost_indices, communicator);
  }



  void
  Vector::reinit(const IndexSet &local_partitioning,
                 const MPI_Comm  comm,
                 const bool      omit_zeroing_entries)
  {
    Assert(communicator != MPI_COMM_NULL,
           ExcMessage("MPI_COMM_NULL passed to Vector::reinit()."));
    communicator   = comm;
    ghosted        = false;
    owned_elements = local_partitioning;

    int err;

    // we do not need to create a new descriptor if the existing one has already
    // been assembled
    if (psblas_descriptor.get() == nullptr)
      {
        psblas_descriptor = std::shared_ptr<psb_c_descriptor>(
          psb_c_new_descriptor(), internal::PSBLASDescriptorDeleter{});

        // Use get_index_vector() from IndexSet to get the indexes
        const std::vector<types::global_dof_index> &indexes =
          local_partitioning.get_index_vector();

        psb_i_t number_of_local_indexes =
          indexes.size(); // Number of local indexes
        // Copy the indexes into a psb_l_t array called vl
        psb_l_t *vl =
          (psb_l_t *)malloc(number_of_local_indexes * sizeof(psb_l_t));
        for (psb_i_t i = 0; i < number_of_local_indexes; ++i)
          {
            vl[i] = static_cast<psb_l_t>(indexes[i]);
          }

        // Insert the indexes into the descriptor
        psblas_context = InitFinalize::get_psblas_context();
        err            = psb_c_cdall_vl(number_of_local_indexes,
                             vl,
                             *psblas_context,
                             psblas_descriptor.get());

        // Free the vl array
        free(vl);
        Assert(err == 0, internal::ExcInitializePSBLASDescriptor(err));
      }

    // Create a new PSBLAS vector and allocate mem space for vector
    psblas_vector = psb_c_new_dvector();

    err = psb_c_dgeall_remote(psblas_vector, psblas_descriptor.get());
    Assert(err == 0, internal::ExcInitializePSBLASVector(err));

    if (omit_zeroing_entries == false)
      {
        err = psb_c_dvect_set_scal(psblas_vector, 0.0);
        Assert(err == 0,
               internal::ExcCallingPSBLASFunction(err, "psb_c_dvect_set_scal"));
      }
    state       = internal::State::Build;
    last_action = VectorOperation::unknown;
  }



  void
  Vector::reinit(const IndexSet &local_partitioning,
                 const IndexSet &ghosts,
                 const MPI_Comm  comm)
  {
    Assert(comm != MPI_COMM_NULL,
           ExcMessage("MPI_COMM_NULL passed to Vector::reinit()."));
    communicator   = comm;
    ghosted        = true;
    owned_elements = local_partitioning;

    ghost_indices = ghosts;
    ghost_indices.subtract_set(local_partitioning);

    int err;
    if (psblas_descriptor.get() == nullptr)
      {
        psblas_descriptor = std::shared_ptr<psb_c_descriptor>(
          psb_c_new_descriptor(), internal::PSBLASDescriptorDeleter{});

        // Use get_index_vector() from IndexSet to get the indexes
        const std::vector<types::global_dof_index> &indexes =
          local_partitioning.get_index_vector();

        psb_i_t number_of_local_indexes =
          indexes.size(); // Number of local indexes
        // Copy the indexes into a psb_l_t array called vl
        psb_l_t *vl =
          (psb_l_t *)malloc(number_of_local_indexes * sizeof(psb_l_t));
        psb_i_t *lidx =
          (psb_i_t *)malloc(number_of_local_indexes * sizeof(psb_i_t));


        for (psb_i_t i = 0; i < number_of_local_indexes; ++i)
          {
            vl[i]   = static_cast<psb_l_t>(indexes[i]);
            lidx[i] = i;
          }

        // Ghost case. From the manual:
        // 1) psb_cdall() with vl (global) and lidx (local)
        // 2) psb_cdins(nz,ja,desc,info,lidx=lidx), where:
        // - ja contains halo indices
        // - lidx: corresponding local indices

        // Insert the indexes into the descriptor
        psblas_context = InitFinalize::get_psblas_context();
        err            = psb_c_cdall_vl_lidx(number_of_local_indexes,
                                  vl,
                                  lidx,
                                  *psblas_context,
                                  psblas_descriptor.get());

        // Free vl array
        free(vl);
        free(lidx);
        Assert(err == 0, internal::ExcInitializePSBLASDescriptor(err));

        // ... insert the ghost indices ...
        const std::vector<types::global_dof_index> &ghost_indexes =
          ghost_indices.get_index_vector();
        const int number_of_ghost_indices = ghost_indexes.size();

        psb_l_t *global_ghost_indices =
          (psb_l_t *)malloc(number_of_ghost_indices * sizeof(psb_l_t));
        psb_i_t *local_ghost_indices =
          (psb_i_t *)malloc(number_of_ghost_indices * sizeof(psb_i_t));

        psb_i_t extended_idx_counter = number_of_local_indexes;
        for (psb_i_t i = 0; i < number_of_ghost_indices; ++i)
          {
            global_ghost_indices[i] = static_cast<psb_l_t>(ghost_indexes[i]);
            local_ghost_indices[i]  = extended_idx_counter++;
          }

        err = psb_c_cdins_lidx(number_of_ghost_indices,
                               global_ghost_indices,
                               local_ghost_indices,
                               psblas_descriptor.get());
        free(global_ghost_indices);
        free(local_ghost_indices);

        Assert(err == 0,
               internal::ExcCallingPSBLASFunction(err, "psb_c_cdins_lidx"));
      }
    // ... create and finalize vector
    psblas_vector = psb_c_new_dvector();
    err           = psb_c_dgeall_remote(psblas_vector, psblas_descriptor.get());

    // ...and descriptor
    err = psb_c_cdasb(psblas_descriptor.get());

    Assert(err == 0, internal::ExcInitializePSBLASVector(err));
    state       = internal::State::Build;
    last_action = VectorOperation::unknown;
  }



  void
  Vector::reinit(const Vector &v, const bool omit_zeroing_entries)
  {
    psblas_descriptor = v.psblas_descriptor;
    if (v.has_ghost_elements())
      {
        reinit(v.locally_owned_elements(),
               v.ghost_indices,
               v.get_mpi_communicator());
        if (!omit_zeroing_entries)
          {
            int err = psb_c_dvect_set_scal(psblas_vector, 0.0);
            Assert(err == 0,
                   internal::ExcCallingPSBLASFunction(err,
                                                      "psb_c_dvect_set_scal"));
          }
      }
    else
      {
        reinit(v.owned_elements,
               v.get_mpi_communicator(),
               omit_zeroing_entries);
      }
  }



  Vector &
  Vector::operator=(const Vector &v)
  {
    Assert(v.last_action == VectorOperation::unknown,
           internal::ExcWrongMode(VectorOperation::unknown, v.last_action));
    Assert(last_action == VectorOperation::unknown,
           internal::ExcWrongMode(VectorOperation::unknown, last_action));
    Assert(v.state == internal::State::Assembled,
           internal::ExcInvalidStateAssembled(v.state));

    // Vectors can have different sizes. If so, we first resize the vector
    if (size() != v.size())
      {
        //  since v has been assembled, we can recycle its descriptor
        psblas_descriptor = v.psblas_descriptor;
        if (v.has_ghost_elements())
          reinit(v.locally_owned_elements(),
                 v.ghost_indices,
                 v.get_mpi_communicator());
        else
          reinit(v.owned_elements, v.get_mpi_communicator(), true);
      }

    // copy local values using axpby routine
    int err = psb_c_dgeaxpby(
      1.0, v.psblas_vector, 0.0, psblas_vector, psblas_descriptor.get());
    AssertThrow(err == 0, internal::ExcAXPBY(err));

    // If current vector has ghost elements, update them
    if (has_ghost_elements())
      update_ghost_values();

    state = internal::State::Assembled;

    return *this;
  }



  Vector &
  Vector::operator=(const value_type s)
  {
    AssertIsFinite(s);
    Assert(state != internal::State::Default, internal::ExcInvalidState(state));
    int err = psb_c_dvect_set_scal(psblas_vector, s);
    Assert(err == 0,
           internal::ExcCallingPSBLASFunction(err, "psb_c_dvect_set_scal"));
    return *this;
  }



  MPI_Comm
  Vector::get_mpi_communicator() const
  {
    return communicator;
  }



  psb_c_descriptor *
  Vector::get_psblas_descriptor() const
  {
    return psblas_descriptor.get();
  }



  psb_c_dvector *
  Vector::get_psblas_vector() const
  {
    return psblas_vector;
  }



  void
  Vector::clear()
  {
    if (state != internal::State::Default)
      {
        int err = psb_c_dgefree(psblas_vector, psblas_descriptor.get());
        Assert(err == 0, internal::ExcFreePSBLASVector(err));

        // Reset the vector
        psblas_vector = nullptr;
        owned_elements.clear();
        owned_elements.set_size(0);
        ghost_indices.clear();
        owned_elements.set_size(0);
        state       = internal::State::Default;
        last_action = VectorOperation::unknown;
      }
  }



  Vector::value_type
  Vector::linfty_norm() const
  {
    Assert(state == internal::State::Assembled,
           internal::ExcInvalidStateAssembled(state));
    return psb_c_dgenrmi(psblas_vector, psblas_descriptor.get());
  }



  Vector::value_type
  Vector::l1_norm() const
  {
    Assert(state == internal::State::Assembled,
           internal::ExcInvalidStateAssembled(state));
    return psb_c_dgeasum(psblas_vector, psblas_descriptor.get());
  }



  Vector::value_type
  Vector::l2_norm() const
  {
    Assert(state == internal::State::Assembled,
           internal::ExcInvalidStateAssembled(state));
    return psb_c_dgenrm2(psblas_vector, psblas_descriptor.get());
  }



  bool
  Vector::all_zero() const
  {
    // we get a pointer to the underlying vector and check if all
    // entries are zero.
    const value_type *start_ptr = psb_c_dvect_f_get_pnt(psblas_vector);
    Assert(start_ptr != nullptr,
           ExcMessage("Error getting underlying PSBLAS vector."));

    const value_type *ptr = start_ptr, *eptr = start_ptr + locally_owned_size();
    bool              flag = true;
    while (ptr != eptr)
      {
        if (*ptr != value_type())
          {
            flag = false;
            break;
          }
        ++ptr;
      }

    unsigned int has_nonzero = flag ? 0 : 1;

    // check that the vector is zero on all processors
    unsigned int num_nonzero = Utilities::MPI::sum(has_nonzero, communicator);
    return num_nonzero == 0;
  }



  void
  Vector::equ(const value_type a, const Vector &v)
  {
    AssertDimension(size(), v.size());
    Assert(!has_ghost_elements(), ExcGhostsPresent());
    AssertIsFinite(a);

    // do the analogous to VecAXPBY(vector, a, 0.0, v.vector) in PETSc
    int err = psb_c_dgeaxpby(
      a, v.psblas_vector, 0.0, psblas_vector, psblas_descriptor.get());
    AssertThrow(err == 0, internal::ExcAXPBY(err));
  }



  Vector::size_type
  Vector::locally_owned_size() const
  {
    return owned_elements.n_elements();
  }



  Vector::value_type
  Vector::operator*(const Vector &v) const
  {
    AssertDimension(size(), v.size());
    return psb_c_dgedot(psblas_vector,
                        v.psblas_vector,
                        psblas_descriptor.get());
  }



  Vector &
  Vector::operator-=(const Vector &v)
  {
    AssertDimension(size(), v.size());
    Assert(!has_ghost_elements(), ExcGhostsPresent());
    int err = psb_c_dgeaxpby(value_type(-1.0),
                             v.psblas_vector,
                             1.0,
                             psblas_vector,
                             psblas_descriptor.get());
    AssertThrow(err == 0, internal::ExcAXPBY(err));
    return *this;
  }


  void
  Vector::set(const std::vector<Vector::size_type>  &indices,
              const std::vector<Vector::value_type> &values)
  {
    Assert(!has_ghost_elements(), ExcGhostsPresent());
    Assert(state != internal::State::Default, internal::ExcInvalidState(state));
    AssertDimension(indices.size(), values.size());

    psb_i_t nz = indices.size(); // Number of non-zero entries

    // Allocate memory for row indices and values. We need to subtract the
    // current value in order to set the value.
    psb_l_t *irw = (psb_l_t *)malloc(nz * sizeof(psb_l_t));
    psb_d_t *val = (psb_d_t *)malloc(nz * sizeof(psb_d_t));
    for (psb_i_t i = 0; i < nz; ++i)
      {
        irw[i] = static_cast<psb_l_t>(indices[i]);
        val[i] = values[i] -
                 psb_c_dgetelem(psblas_vector, irw[i], psblas_descriptor.get());
      }

    int err =
      psb_c_dgeins(nz /*nz*/, irw, val, psblas_vector, psblas_descriptor.get());

    // Free allocated memory
    free(irw);
    free(val);
    Assert(err == 0, internal::ExcInsertionInPSBLASVector(err));
  }



  void
  Vector::add(const std::vector<Vector::size_type>  &indices,
              const std::vector<Vector::value_type> &values)
  {
    Assert(!has_ghost_elements(), ExcGhostsPresent());
    Assert(state != internal::State::Default, internal::ExcInvalidState(state));
    Assert(indices.size() == values.size(),
           ExcMessage("Indices and values size mismatch."));

    psb_i_t nz = indices.size(); // Number of non-zero entries

    // Allocate memory for row indices and values
    psb_l_t *irw = (psb_l_t *)malloc(nz * sizeof(psb_l_t));
    psb_d_t *val = (psb_d_t *)malloc(nz * sizeof(psb_d_t));
    for (psb_i_t i = 0; i < nz; ++i)
      {
        irw[i] = static_cast<psb_l_t>(indices[i]);
        val[i] = values[i];
      }

    int err =
      psb_c_dgeins(nz /*nz*/, irw, val, psblas_vector, psblas_descriptor.get());
    // Free allocated memory
    free(irw);
    free(val);
    Assert(err == 0, internal::ExcInsertionInPSBLASVector(err));
  }



  void
  Vector::add(const value_type s, const Vector &V)
  {
    AssertIsFinite(s);
    Assert(!has_ghost_elements(), ExcGhostsPresent());
    int err = psb_c_dgeaxpby(s,
                             V.psblas_vector,
                             value_type(1.0),
                             psblas_vector,
                             psblas_descriptor.get());
    Assert(err == 0, internal::ExcAXPBY(err));
  }



  void
  Vector::sadd(const value_type s, const Vector &V)
  {
    Assert(!has_ghost_elements(), ExcGhostsPresent());
    AssertIsFinite(s);
    int err = psb_c_dgeaxpby(value_type(1.0),
                             V.psblas_vector,
                             s,
                             psblas_vector,
                             psblas_descriptor.get());
    Assert(err == 0, internal::ExcAXPBY(err));
  }


  void
  Vector::sadd(const value_type s, const value_type a, const Vector &V)
  {
    Assert(!has_ghost_elements(), ExcGhostsPresent());
    AssertIsFinite(s);
    int err = psb_c_dgeaxpby(
      a, V.psblas_vector, s, psblas_vector, psblas_descriptor.get());
    Assert(err == 0, internal::ExcAXPBY(err));
  }



  Vector::value_type
  Vector::add_and_dot(const value_type a, const Vector &V, const Vector &W)
  {
    // forward the call to add()...
    add(a, V);
    // ... and then do the dot product
    return psb_c_dgedot(psblas_vector,
                        W.psblas_vector,
                        psblas_descriptor.get());
  }



  void
  Vector::swap(Vector &v)
  {
    // simply swap pointers
    std::swap(psblas_vector, v.psblas_vector);
    std::swap(psblas_descriptor, v.psblas_descriptor);
    std::swap(ghosted, v.ghosted);
    std::swap(this->last_action, v.last_action);
    // We use a temp variable to swap the IndexSets
    IndexSet temp(ghost_indices);
    ghost_indices   = v.ghost_indices;
    v.ghost_indices = temp;
  }



  void
  Vector::compress(const VectorOperation::values operation)
  {
    Assert(has_ghost_elements() == false,
           ExcMessage("Calling compress() is only useful if a vector "
                      "has been written into, but this is a vector with ghost "
                      "elements and consequently is read-only. It does "
                      "not make sense to call compress() for such "
                      "vectors."));
    AssertThrow(operation == VectorOperation::values::add, ExcNotImplemented());
    // We check the state of the vector...
    Assert(state != internal::State::Default, internal::ExcInvalidDefault());
    // ... and if the last action was compatible with what we want to do
    AssertThrow(
      last_action == VectorOperation::unknown || last_action == operation,
      ExcMessage(
        "Missing compress() or calling with wrong VectorOperation argument."));

    // We check if the descriptor has already been assembled somewhere
    int err;
    if (!psb_c_cd_is_asb(psblas_descriptor.get()))
      {
        err = psb_c_cdasb(psblas_descriptor.get());
        Assert(err == 0, internal::ExcAssemblePSBLASDescriptor(err));
      }

    // finally, we perform the assemble operation
    err = psb_c_dgeasb(psblas_vector, psblas_descriptor.get());

    Assert(err == 0, internal::ExcAssemblePSBLASVector(err));
    state       = internal::State::Assembled;
    last_action = VectorOperation::unknown;
  }



  void
  Vector::update_ghost_values() const
  {
    if (ghosted)
      {
        int err = psb_c_dhalo(psblas_vector, psblas_descriptor.get());
        Assert(err == 0, ExcMessage("Error while updating ghost values."));
      }
  }


  std::size_t
  Vector::memory_consumption() const
  {
    std::size_t mem = MemoryConsumption::memory_consumption(ghosted) +
                      MemoryConsumption::memory_consumption(ghost_indices);

    // mem += psb_c_sizeof(psblas_vector); TODO[MF]: missing from PSBLAS's
    // interface
    return mem;
  }


} // namespace PSCToolkit

DEAL_II_NAMESPACE_CLOSE
#endif
