// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_trilinos_tpetra_precondition_templates_h
#define dealii_trilinos_tpetra_precondition_templates_h


#include <deal.II/base/config.h>

#include <deal.II/base/index_set.h>

#include <string>

#ifdef DEAL_II_TRILINOS_WITH_TPETRA
#  ifdef DEAL_II_TRILINOS_WITH_IFPACK2

#    include <deal.II/lac/trilinos_tpetra_precondition.h>
#    include <deal.II/lac/trilinos_tpetra_sparse_matrix.h>

#    include <Ifpack2_Factory.hpp>
#    include <Ifpack2_Preconditioner.hpp>

DEAL_II_NAMESPACE_OPEN


namespace LinearAlgebra
{
  namespace TpetraWrappers
  {
    /* ---------------------- PreconditionBase ------------------------ */
    template <typename Number, typename MemorySpace>
    void
    PreconditionBase<Number, MemorySpace>::clear()
    {
      preconditioner.reset();
    }



    template <typename Number, typename MemorySpace>
    Teuchos::RCP<Tpetra::Operator<
      Number,
      int,
      types::signed_global_dof_index,
      typename PreconditionBase<Number, MemorySpace>::NodeType>>
    PreconditionBase<Number, MemorySpace>::trilinos_operator() const
    {
      return preconditioner;
    }



    template <typename Number, typename MemorySpace>
    IndexSet
    PreconditionBase<Number, MemorySpace>::locally_owned_domain_indices() const
    {
      return IndexSet(preconditioner->getDomainMap());
    }



    template <typename Number, typename MemorySpace>
    IndexSet
    PreconditionBase<Number, MemorySpace>::locally_owned_range_indices() const
    {
      return IndexSet(preconditioner->getRangeMap());
    }



    template <typename Number, typename MemorySpace>
    inline void
    PreconditionBase<Number, MemorySpace>::vmult(
      Vector<Number, MemorySpace>       &dst,
      const Vector<Number, MemorySpace> &src) const
    {
      Assert(dst.trilinos_vector().getMap()->isSameAs(
               *(preconditioner->getRangeMap())),
             ExcNonMatchingMaps("dst"));
      Assert(src.trilinos_rcp()->getMap()->isSameAs(
               *(preconditioner->getDomainMap())),
             ExcNonMatchingMaps("src"));

      preconditioner->apply(src.trilinos_vector(), dst.trilinos_vector());
    }



    template <typename Number, typename MemorySpace>
    inline void
    PreconditionBase<Number, MemorySpace>::Tvmult(
      Vector<Number, MemorySpace>       &dst,
      const Vector<Number, MemorySpace> &src) const
    {
      Assert(preconditioner->hasTransposeApply(), ExcTransposeNotSupported());

      Assert(dst.trilinos_rcp()->getMap()->isSameAs(
               *(preconditioner->getDomainMap())),
             ExcNonMatchingMaps("dst"));

      Assert(src.trilinos_rcp()->getMap()->isSameAs(
               *(preconditioner->getRangeMap())),
             ExcNonMatchingMaps("src"));

      preconditioner->apply(src.trilinos_vector(),
                            dst.trilinos_vector(),
                            Teuchos::TRANS);
    }



    /* ---------------------- PreconditionIfpackBase ------------------------ */
    template <typename Number, typename MemorySpace>
    PreconditionIfpackBase<Number, MemorySpace>::PreconditionIfpackBase(
      const std::string &preconditioner_type)
      : preconditioner_type(preconditioner_type)
    {}



    template <typename Number, typename MemorySpace>
    void
    PreconditionIfpackBase<Number, MemorySpace>::initialize(
      const SparseMatrix<Number, MemorySpace> &A)
    {
      Ifpack2::Factory factory;
      Teuchos::RCP<Ifpack2::Preconditioner<
        Number,
        int,
        dealii::types::signed_global_dof_index,
        typename PreconditionBase<Number, MemorySpace>::NodeType>>
        ifpack2Preconditioner;
      ifpack2Preconditioner =
        factory.create(preconditioner_type, A.trilinos_rcp());
      ifpack2Preconditioner->setParameters(this->parameter_list);
      ifpack2Preconditioner->initialize();
      ifpack2Preconditioner->compute();
      this->preconditioner = ifpack2Preconditioner;
    }



    /* ---------------------- PreconditionJacobi ------------------------ */

    template <typename Number, typename MemorySpace>
    PreconditionJacobi<Number, MemorySpace>::AdditionalData::AdditionalData(
      const double omega,
      const int    n_sweeps,
      const bool   fix_diagonal,
      const double min_diagonal)
      : omega(omega)
      , n_sweeps(n_sweeps)
      , fix_diagonal(fix_diagonal)
      , min_diagonal(min_diagonal)
    {}



    template <typename Number, typename MemorySpace>
    PreconditionJacobi<Number, MemorySpace>::PreconditionJacobi()
      : PreconditionIfpackBase<Number, MemorySpace>("RELAXATION")
    {}



    template <typename Number, typename MemorySpace>
    void
    PreconditionJacobi<Number, MemorySpace>::initialize(
      const SparseMatrix<Number, MemorySpace> &A,
      const AdditionalData                    &ad)
    {
      this->parameter_list.set("relaxation: type", "Jacobi");
      this->parameter_list.set("relaxation: damping factor", ad.omega);
      this->parameter_list.set("relaxation: sweeps", ad.n_sweeps);
      this->parameter_list.set("relaxation: fix tiny diagonal entries",
                               ad.fix_diagonal);
      this->parameter_list.set("relaxation: min diagonal value",
                               ad.min_diagonal);
      PreconditionIfpackBase<Number, MemorySpace>::initialize(A);
    }



    /* ---------------------- PreconditionL1Jacobi ------------------------ */
    template <typename Number, typename MemorySpace>
    PreconditionL1Jacobi<Number, MemorySpace>::AdditionalData::AdditionalData(
      const double omega,
      const int    n_sweeps,
      const bool   fix_diagonal,
      const double min_diagonal)
      : omega(omega)
      , n_sweeps(n_sweeps)
      , fix_diagonal(fix_diagonal)
      , min_diagonal(min_diagonal)
    {}



    template <typename Number, typename MemorySpace>
    PreconditionL1Jacobi<Number, MemorySpace>::PreconditionL1Jacobi()
      : PreconditionIfpackBase<Number, MemorySpace>("RELAXATION")
    {}



    template <typename Number, typename MemorySpace>
    void
    PreconditionL1Jacobi<Number, MemorySpace>::initialize(
      const SparseMatrix<Number, MemorySpace> &A,
      const AdditionalData                    &ad)
    {
      this->parameter_list.set("relaxation: type", "Jacobi");
      this->parameter_list.set("relaxation: use l1", true);
      this->parameter_list.set("relaxation: damping factor", ad.omega);
      this->parameter_list.set("relaxation: sweeps", ad.n_sweeps);
      this->parameter_list.set("relaxation: fix tiny diagonal entries",
                               ad.fix_diagonal);
      this->parameter_list.set("relaxation: min diagonal value",
                               ad.min_diagonal);
      PreconditionIfpackBase<Number, MemorySpace>::initialize(A);
    }



  } // namespace TpetraWrappers
} // namespace LinearAlgebra

DEAL_II_NAMESPACE_CLOSE

#  endif // DEAL_II_TRILINOS_WITH_IFPACK2
#endif   // DEAL_II_TRILINOS_WITH_TPETRA

#endif
