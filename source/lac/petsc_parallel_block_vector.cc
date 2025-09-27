// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2004 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#include <deal.II/lac/petsc_block_vector.h>

#ifdef DEAL_II_WITH_PETSC

#  include <deal.II/lac/petsc_compatibility.h>

DEAL_II_NAMESPACE_OPEN

namespace PETScWrappers
{
  namespace MPI
  {
    // Check that the class we declare here satisfies the
    // vector-space-vector concept. If we catch it here,
    // any mistake in the vector class declaration would
    // show up in uses of this class later on as well.
#  ifdef DEAL_II_HAVE_CXX20
    static_assert(concepts::is_vector_space_vector<BlockVector>);
#  endif

    using size_type = types::global_dof_index;

    BlockVector::~BlockVector()
    {
      PetscErrorCode ierr = VecDestroy(&petsc_nest_vector);
      AssertThrow(ierr == 0, ExcPETScError(ierr));
    }

    void
    BlockVector::reinit(const unsigned int num_blocks)
    {
      std::vector<size_type> block_sizes(num_blocks, 0);
      this->block_indices.reinit(block_sizes);
      if (this->components.size() != this->n_blocks())
        this->components.resize(this->n_blocks());

      for (unsigned int i = 0; i < this->n_blocks(); ++i)
        components[i].reinit(MPI_COMM_SELF, 0, 0);

      collect_sizes();
    }

    void
    BlockVector::reinit(Vec v)
    {
      PetscBool isnest;

      PetscErrorCode ierr =
        PetscObjectTypeCompare(reinterpret_cast<PetscObject>(v),
                               VECNEST,
                               &isnest);
      AssertThrow(ierr == 0, ExcPETScError(ierr));
      std::vector<Vec> sv;
      if (isnest)
        {
          PetscInt nb;
          ierr = VecNestGetSize(v, &nb);
          AssertThrow(ierr == 0, ExcPETScError(ierr));
          for (PetscInt i = 0; i < nb; ++i)
            {
              Vec vv;
              ierr = VecNestGetSubVec(v, i, &vv);
              sv.push_back(vv);
            }
        }
      else
        {
          sv.push_back(v);
        }

      auto nb = sv.size();

      std::vector<size_type> block_sizes(nb, 0);
      this->block_indices.reinit(block_sizes);

      this->components.resize(nb);
      for (unsigned int i = 0; i < nb; ++i)
        {
          this->components[i].reinit(sv[i]);
        }

      BlockVectorBase::collect_sizes();
      if (!isnest)
        setup_nest_vec();
      else
        {
          ierr = PetscObjectReference(reinterpret_cast<PetscObject>(v));
          AssertThrow(ierr == 0, ExcPETScError(ierr));
          PetscErrorCode ierr = VecDestroy(&petsc_nest_vector);
          AssertThrow(ierr == 0, ExcPETScError(ierr));
          petsc_nest_vector = v;
        }
    }

    Vec &
    BlockVector::petsc_vector()
    {
      return petsc_nest_vector;
    }

    BlockVector::operator const Vec &() const
    {
      return petsc_nest_vector;
    }

    void
    BlockVector::collect_sizes()
    {
      BlockVectorBase::collect_sizes();
      setup_nest_vec();
    }

    void
    BlockVector::compress(VectorOperation::values operation)
    {
      BlockVectorBase::compress(operation);
      petsc_increment_state_counter(petsc_nest_vector);
    }

    void
    BlockVector::setup_nest_vec()
    {
      PetscErrorCode ierr = VecDestroy(&petsc_nest_vector);
      AssertThrow(ierr == 0, ExcPETScError(ierr));

      auto n = this->n_blocks();

      std::vector<Vec> pcomponents(n);
      for (unsigned int i = 0; i < n; i++)
        pcomponents[i] = this->components[i].petsc_vector();

      MPI_Comm comm =
        pcomponents.size() > 0 ?
          PetscObjectComm(reinterpret_cast<PetscObject>(pcomponents[0])) :
          PETSC_COMM_SELF;

      ierr =
        VecCreateNest(comm, n, nullptr, pcomponents.data(), &petsc_nest_vector);
      AssertThrow(ierr == 0, ExcPETScError(ierr));
    }
  } // namespace MPI

} // namespace PETScWrappers

DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_WITH_PETSC
