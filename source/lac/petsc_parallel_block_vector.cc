// ---------------------------------------------------------------------
//
// Copyright (C) 2004 - 2018 by the deal.II authors
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

#include <deal.II/lac/petsc_block_vector.h>

#ifdef DEAL_II_WITH_PETSC

DEAL_II_NAMESPACE_OPEN

namespace PETScWrappers
{
  namespace MPI
  {
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
    BlockVector::assign_petsc_vector(Vec v)
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
          this->components[i].assign_petsc_vector(sv[i]);
        }

      this->collect_sizes();
    }

    Vec &
    BlockVector::petsc_vector()
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
        VecCreateNest(comm, n, NULL, pcomponents.data(), &petsc_nest_vector);
      AssertThrow(ierr == 0, ExcPETScError(ierr));
    }
  } // namespace MPI

} // namespace PETScWrappers

DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_WITH_PETSC
