// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2024 - 2025 by the deal.II authors
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

#include <deal.II/base/exceptions.h>
#include <deal.II/base/index_set.h>

#include <deal.II/lac/trilinos_tpetra_types.h>

#include <string>

#ifdef DEAL_II_TRILINOS_WITH_TPETRA
#  ifdef DEAL_II_TRILINOS_WITH_IFPACK2
#    include <deal.II/lac/trilinos_tpetra_precondition.h>
#    include <deal.II/lac/trilinos_tpetra_sparse_matrix.h>

#    include <Ifpack2_Factory.hpp>
#    include <Ifpack2_IdentitySolver.hpp>
#    include <Ifpack2_IdentitySolver_decl.hpp>
#    include <Ifpack2_Preconditioner.hpp>
#    include <Kokkos_Core.hpp>
#    include <Kokkos_DualView.hpp>
#    include <Teuchos_RCPDecl.hpp>
#    include <Tpetra_CrsMatrix_decl.hpp>
#    include <Tpetra_Vector_decl.hpp>

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
    const TpetraTypes::LinearOperator<Number, MemorySpace> &
    PreconditionBase<Number, MemorySpace>::trilinos_operator() const
    {
      return *preconditioner;
    }



    template <typename Number, typename MemorySpace>
    Teuchos::RCP<TpetraTypes::LinearOperator<Number, MemorySpace>>
    PreconditionBase<Number, MemorySpace>::trilinos_rcp() const
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



    template <typename Number, typename MemorySpace>
    inline void
    PreconditionBase<Number, MemorySpace>::vmult(
      dealii::Vector<Number> &dst,
      dealii::Vector<Number> &src) const
    {
      AssertDimension(dst.size(),
                      preconditioner->getRangeMap()->getGlobalNumElements());
      AssertDimension(src.size(),
                      preconditioner->getDomainMap()->getGlobalNumElements());

      // Construct Tpetra Vectors as views to the actual dealii::Vector data
      // Since all of Tpetra is based on Kokkos::DualView this is a bit of a
      // process

      // 1. Make host view with underlying dealii::Vector data.
      //    (This is always a n x 1 matrix since multivector expects that.)

      // Tpetra uses Kokkos::complex<> internally for std::complex<>
      using value_type = typename TpetraTypes::HostViewType<Number>::value_type;
      static_assert((numbers::NumberTraits<Number>::is_complex &&
                     sizeof(value_type) == sizeof(Number)) ||
                    (!numbers::NumberTraits<Number>::is_complex &&
                     std::is_same_v<value_type, Number>));

      TpetraTypes::HostViewType<Number> host_src(
        reinterpret_cast<value_type *>(src.data()), src.size(), 1);
      TpetraTypes::HostViewType<Number> host_dst(
        reinterpret_cast<value_type *>(dst.data()), dst.size(), 1);

      // 2. Create device mirror (if Device == Host, no alloc )
      auto dev_src = Kokkos::create_mirror_view_and_copy(
        typename MemorySpace::kokkos_space::execution_space(), host_src);

      auto dev_dst = Kokkos::create_mirror_view_and_copy(
        typename MemorySpace::kokkos_space::execution_space(), host_dst);

      // 3. Create dual views from the previous ones
      TpetraTypes::DualViewType<Number, MemorySpace> dual_src(dev_src,
                                                              host_src);
      TpetraTypes::DualViewType<Number, MemorySpace> dual_dst(dev_dst,
                                                              host_dst);

      // 4. Create Tpetra Vector from
      TpetraTypes::VectorType<Number, MemorySpace> tpetra_dst(
        preconditioner->getRangeMap(), dual_dst);
      TpetraTypes::VectorType<Number, MemorySpace> tpetra_src(
        preconditioner->getDomainMap(), dual_src);

      // 5. Finally, apply our preconditioner
      preconditioner->apply(tpetra_src, tpetra_dst);
    }



    template <typename Number, typename MemorySpace>
    inline void
    PreconditionBase<Number, MemorySpace>::Tvmult(
      dealii::Vector<Number> &dst,
      dealii::Vector<Number> &src) const
    {
      AssertDimension(dst.size(),
                      preconditioner->getDomainMap()->getGlobalNumElements());
      AssertDimension(src.size(),
                      preconditioner->getRangeMap()->getGlobalNumElements());

      Assert(preconditioner->hasTransposeApply(), ExcTransposeNotSupported());
      // Construct Tpetra Vectors as views to the actual dealii::Vector data
      // Since all of Tpetra is based on Kokkos::DualView this is a bit of a
      // process

      // 1. Make host view with underlying dealii::Vector data.
      //    (This is always a n x 1 matrix since multivector expects that.)

      // Tpetra uses Kokkos::complex<> internally for std::complex<>
      using value_type = typename TpetraTypes::HostViewType<Number>::value_type;
      static_assert((numbers::NumberTraits<Number>::is_complex &&
                     sizeof(value_type) == sizeof(Number)) ||
                    (!numbers::NumberTraits<Number>::is_complex &&
                     std::is_same_v<value_type, Number>));

      TpetraTypes::HostViewType<Number> host_src(
        reinterpret_cast<value_type *>(src.data()), src.size(), 1);
      TpetraTypes::HostViewType<Number> host_dst(
        reinterpret_cast<value_type *>(dst.data()), dst.size(), 1);

      // 2. Create device mirror (if Device == Host, no alloc )
      auto dev_src = Kokkos::create_mirror_view_and_copy(
        typename MemorySpace::kokkos_space::execution_space(), host_src);

      auto dev_dst = Kokkos::create_mirror_view_and_copy(
        typename MemorySpace::kokkos_space::execution_space(), host_dst);

      // 3. Create dual views from the previous ones
      TpetraTypes::DualViewType<Number, MemorySpace> dual_src(dev_src,
                                                              host_src);
      TpetraTypes::DualViewType<Number, MemorySpace> dual_dst(dev_dst,
                                                              host_dst);

      // 4. Create Tpetra Vector from
      TpetraTypes::VectorType<Number, MemorySpace> tpetra_dst(
        preconditioner->getRangeMap(), dual_dst);
      TpetraTypes::VectorType<Number, MemorySpace> tpetra_src(
        preconditioner->getDomainMap(), dual_src);

      // 5. Finally, apply our preconditioner
      preconditioner->apply(tpetra_src, tpetra_dst, Teuchos::TRANS);
    }



    /* ---------------------- PreconditionIdentity ------------------------ */
    template <typename Number, typename MemorySpace>
    void
    PreconditionIdentity<Number, MemorySpace>::initialize(
      const SparseMatrix<Number, MemorySpace> &A)
    {
      Teuchos::RCP<Ifpack2::IdentitySolver<
        TpetraTypes::RowMatrixType<Number, MemorySpace>>>
        ifpack2Preconditioner;
      ifpack2Preconditioner =
        rcp(new Ifpack2::IdentitySolver<
            TpetraTypes::RowMatrixType<Number, MemorySpace>>(A.trilinos_rcp()));

      ifpack2Preconditioner->initialize();
      ifpack2Preconditioner->compute();
      this->preconditioner = ifpack2Preconditioner;
    }



    /* ---------------------- PreconditionIfpackBase ------------------------ */
    template <typename Number, typename MemorySpace>
    PreconditionIfpackBase<Number, MemorySpace>::PreconditionIfpackBase(
      const std::string &preconditioner_type)
      : preconditioner_type(preconditioner_type)
    {
#    if DEAL_II_TRILINOS_VERSION_GTE(13, 0, 0)
      Ifpack2::Factory factory;
      AssertThrow(
        (factory.isSupported<TpetraTypes::MatrixType<Number, MemorySpace>>(
          preconditioner_type)),
        ExcTrilinosIpack2PreconditionerUnsupported(preconditioner_type));
#    endif
    }



    template <typename Number, typename MemorySpace>
    void
    PreconditionIfpackBase<Number, MemorySpace>::initialize(
      const SparseMatrix<Number, MemorySpace> &A)
    {
      Ifpack2::Factory factory;
      Teuchos::RCP<TpetraTypes::Ifpack2PreconType<Number, MemorySpace>>
        ifpack2Preconditioner;
      ifpack2Preconditioner =
        factory.create(preconditioner_type, A.trilinos_rcp());
      ifpack2Preconditioner->setParameters(this->parameter_list);
      ifpack2Preconditioner->initialize();
      ifpack2Preconditioner->compute();
      this->preconditioner = ifpack2Preconditioner;
    }



    /* ---------------------- PreconditionIfpack ------------------------ */
    template <typename Number, typename MemorySpace>
    PreconditionIfpack<Number, MemorySpace>::PreconditionIfpack(
      const std::string &preconditioner_type)
      : PreconditionIfpackBase<Number, MemorySpace>(preconditioner_type)
    {}



    template <typename Number, typename MemorySpace>
    void
    PreconditionIfpack<Number, MemorySpace>::set_parameter_list(
      Teuchos::ParameterList &parameter_list)
    {
      this->parameter_list.setParameters(parameter_list);
    }



    /* ---------------------- PreconditionJacobi ------------------------ */
    template <typename Number, typename MemorySpace>
    PreconditionJacobi<Number, MemorySpace>::AdditionalData::AdditionalData(
      const double omega,
      const bool   fix_diagonal,
      const double min_diagonal,
      const int    n_sweeps)
      : omega(omega)
      , fix_diagonal(fix_diagonal)
      , min_diagonal(min_diagonal)
      , n_sweeps(n_sweeps)
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
      const bool   fix_diagonal,
      const double min_diagonal,
      const int    n_sweeps)
      : omega(omega)
      , fix_diagonal(fix_diagonal)
      , min_diagonal(min_diagonal)
      , n_sweeps(n_sweeps)
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



    /* ---------------------- PreconditionL1GaussSeidel -------------------- */
    template <typename Number, typename MemorySpace>
    PreconditionL1GaussSeidel<Number, MemorySpace>::AdditionalData::
      AdditionalData(const double omega,
                     const double eta,
                     const bool   fix_diagonal,
                     const double min_diagonal,
                     const int    n_sweeps)
      : omega(omega)
      , eta(eta)
      , fix_diagonal(fix_diagonal)
      , min_diagonal(min_diagonal)
      , n_sweeps(n_sweeps)
    {}



    template <typename Number, typename MemorySpace>
    PreconditionL1GaussSeidel<Number, MemorySpace>::PreconditionL1GaussSeidel()
      : PreconditionIfpackBase<Number, MemorySpace>("RELAXATION")
    {}



    template <typename Number, typename MemorySpace>
    void
    PreconditionL1GaussSeidel<Number, MemorySpace>::initialize(
      const SparseMatrix<Number, MemorySpace> &A,
      const AdditionalData                    &ad)
    {
      this->parameter_list.set("relaxation: type", "Gauss-Seidel");
      this->parameter_list.set("relaxation: use l1", true);
      this->parameter_list.set("relaxation: damping factor", ad.omega);
      this->parameter_list.set("relaxation: l1 eta", ad.eta);
      this->parameter_list.set("relaxation: sweeps", ad.n_sweeps);
      this->parameter_list.set("relaxation: fix tiny diagonal entries",
                               ad.fix_diagonal);
      this->parameter_list.set("relaxation: min diagonal value",
                               ad.min_diagonal);
      PreconditionIfpackBase<Number, MemorySpace>::initialize(A);
    }



    /* ---------------------- PreconditionSOR ------------------------ */
    template <typename Number, typename MemorySpace>
    PreconditionSOR<Number, MemorySpace>::AdditionalData::AdditionalData(
      const double omega,
      const int    overlap,
      const bool   fix_diagonal,
      const double min_diagonal,
      const int    n_sweeps)
      : omega(omega)
      , overlap(overlap)
      , fix_diagonal(fix_diagonal)
      , min_diagonal(min_diagonal)
      , n_sweeps(n_sweeps)
    {}



    template <typename Number, typename MemorySpace>
    PreconditionSOR<Number, MemorySpace>::PreconditionSOR()
      : PreconditionIfpackBase<Number, MemorySpace>("SCHWARZ")
    {}



    template <typename Number, typename MemorySpace>
    void
    PreconditionSOR<Number, MemorySpace>::initialize(
      const SparseMatrix<Number, MemorySpace> &A,
      const AdditionalData                    &ad)
    {
      this->parameter_list.set("schwarz: subdomain solver name", "RELAXATION");
      {
        Teuchos::ParameterList &relaxParams =
          this->parameter_list.sublist("schwarz: subdomain solver parameters",
                                       false);
        relaxParams.set("relaxation: type", "Gauss-Seidel");
        relaxParams.set("relaxation: damping factor", ad.omega);
        relaxParams.set("relaxation: sweeps", ad.n_sweeps);
        relaxParams.set("relaxation: fix tiny diagonal entries",
                        ad.fix_diagonal);
        relaxParams.set("relaxation: min diagonal value", ad.min_diagonal);
      }
      this->parameter_list.set("schwarz: combine mode", "ADD");
      this->parameter_list.set("schwarz: overlap level", ad.overlap);
      PreconditionIfpackBase<Number, MemorySpace>::initialize(A);
    }



    /* ---------------------- PreconditionSSOR ----------------------*/
    template <typename Number, typename MemorySpace>
    PreconditionSSOR<Number, MemorySpace>::AdditionalData::AdditionalData(
      const double omega,
      const int    overlap,
      const bool   fix_diagonal,
      const double min_diagonal,
      const int    n_sweeps)
      : omega(omega)
      , overlap(overlap)
      , fix_diagonal(fix_diagonal)
      , min_diagonal(min_diagonal)
      , n_sweeps(n_sweeps)
    {}



    template <typename Number, typename MemorySpace>
    PreconditionSSOR<Number, MemorySpace>::PreconditionSSOR()
      : PreconditionIfpackBase<Number, MemorySpace>("SCHWARZ")
    {}



    template <typename Number, typename MemorySpace>
    void
    PreconditionSSOR<Number, MemorySpace>::initialize(
      const SparseMatrix<Number, MemorySpace> &A,
      const AdditionalData                    &ad)
    {
      this->parameter_list.set("schwarz: subdomain solver name", "RELAXATION");
      {
        Teuchos::ParameterList &relaxParams =
          this->parameter_list.sublist("schwarz: subdomain solver parameters",
                                       false);
        relaxParams.set("relaxation: type", "Symmetric Gauss-Seidel");
        relaxParams.set("relaxation: damping factor", ad.omega);
        relaxParams.set("relaxation: sweeps", ad.n_sweeps);
        relaxParams.set("relaxation: fix tiny diagonal entries",
                        ad.fix_diagonal);
        relaxParams.set("relaxation: min diagonal value", ad.min_diagonal);
      }
      this->parameter_list.set("schwarz: combine mode", "ADD");
      this->parameter_list.set("schwarz: overlap level", ad.overlap);
      PreconditionIfpackBase<Number, MemorySpace>::initialize(A);
    }



    /* ---------------------- PreconditionBlockJacobi ----------------------*/
    template <typename Number, typename MemorySpace>
    PreconditionBlockJacobi<Number, MemorySpace>::AdditionalData::
      AdditionalData(const int    n_local_parts,
                     const double omega,
                     const int    block_overlap,
                     const int    n_sweeps)
      : n_local_parts(n_local_parts)
      , omega(omega)
      , block_overlap(block_overlap)
      , n_sweeps(n_sweeps)
    {}



    template <typename Number, typename MemorySpace>
    PreconditionBlockJacobi<Number, MemorySpace>::PreconditionBlockJacobi()
      : PreconditionIfpackBase<Number, MemorySpace>("BLOCK_RELAXATION")
    {}



    template <typename Number, typename MemorySpace>
    void
    PreconditionBlockJacobi<Number, MemorySpace>::initialize(
      const SparseMatrix<Number, MemorySpace> &A,
      const AdditionalData                    &ad)
    {
      this->parameter_list.set("relaxation: type", "Jacobi");
      this->parameter_list.set("relaxation: damping factor", ad.omega);
      this->parameter_list.set("relaxation: sweeps", ad.n_sweeps);
      this->parameter_list.set("partitioner: local parts", ad.n_local_parts);
      this->parameter_list.set("partitioner: overlap", ad.block_overlap);
      PreconditionIfpackBase<Number, MemorySpace>::initialize(A);
    }



    /* ---------------------- PreconditionBlockSOR ----------------------*/
    template <typename Number, typename MemorySpace>
    PreconditionBlockSOR<Number, MemorySpace>::AdditionalData::AdditionalData(
      const int    n_local_parts,
      const double omega,
      const int    overlap,
      const int    n_sweeps)
      : n_local_parts(n_local_parts)
      , omega(omega)
      , overlap(overlap)
      , n_sweeps(n_sweeps)
    {}



    template <typename Number, typename MemorySpace>
    PreconditionBlockSOR<Number, MemorySpace>::PreconditionBlockSOR()
      : PreconditionIfpackBase<Number, MemorySpace>("SCHWARZ")
    {}



    template <typename Number, typename MemorySpace>
    void
    PreconditionBlockSOR<Number, MemorySpace>::initialize(
      const SparseMatrix<Number, MemorySpace> &A,
      const AdditionalData                    &ad)
    {
      this->parameter_list.set("schwarz: subdomain solver name",
                               "BLOCK_RELAXATION");
      {
        Teuchos::ParameterList &relaxParams =
          this->parameter_list.sublist("schwarz: subdomain solver parameters",
                                       false);
        relaxParams.set("relaxation: type", "Gauss-Seidel");
        relaxParams.set("relaxation: damping factor", ad.omega);
        relaxParams.set("relaxation: sweeps", ad.n_sweeps);
        relaxParams.set("partitioner: local parts", ad.n_local_parts);
      }
      this->parameter_list.set("schwarz: combine mode", "ADD");
      this->parameter_list.set("schwarz: overlap level", ad.overlap);
      PreconditionIfpackBase<Number, MemorySpace>::initialize(A);
    }



    /* ---------------------- PreconditionBlockSSOR ----------------------*/
    template <typename Number, typename MemorySpace>
    PreconditionBlockSSOR<Number, MemorySpace>::AdditionalData::AdditionalData(
      const int    n_local_parts,
      const double omega,
      const int    overlap,
      const int    n_sweeps)
      : n_local_parts(n_local_parts)
      , omega(omega)
      , overlap(overlap)
      , n_sweeps(n_sweeps)
    {}



    template <typename Number, typename MemorySpace>
    PreconditionBlockSSOR<Number, MemorySpace>::PreconditionBlockSSOR()
      : PreconditionIfpackBase<Number, MemorySpace>("SCHWARZ")
    {}



    template <typename Number, typename MemorySpace>
    void
    PreconditionBlockSSOR<Number, MemorySpace>::initialize(
      const SparseMatrix<Number, MemorySpace> &A,
      const AdditionalData                    &ad)
    {
      this->parameter_list.set("schwarz: subdomain solver name",
                               "BLOCK_RELAXATION");
      {
        Teuchos::ParameterList &relaxParams =
          this->parameter_list.sublist("schwarz: subdomain solver parameters",
                                       false);
        relaxParams.set("relaxation: type", "Symmetric Gauss-Seidel");
        relaxParams.set("relaxation: damping factor", ad.omega);
        relaxParams.set("relaxation: sweeps", ad.n_sweeps);
        relaxParams.set("partitioner: local parts", ad.n_local_parts);
      }
      this->parameter_list.set("schwarz: combine mode", "ADD");
      this->parameter_list.set("schwarz: overlap level", ad.overlap);
      PreconditionIfpackBase<Number, MemorySpace>::initialize(A);
    }



    /* ---------------------- PreconditionChebyshev ----------------------*/
    template <typename Number, typename MemorySpace>
    PreconditionChebyshev<Number, MemorySpace>::AdditionalData::AdditionalData(
      const int    degree,
      const double max_eigenvalue,
      const double eigenvalue_ratio,
      const double min_eigenvalue,
      const double min_diagonal,
      const bool   nonzero_starting)
      : degree(degree)
      , max_eigenvalue(max_eigenvalue)
      , min_eigenvalue(min_eigenvalue)
      , eigenvalue_ratio(eigenvalue_ratio)
      , min_diagonal(min_diagonal)
      , nonzero_starting(nonzero_starting)
    {}



    template <typename Number, typename MemorySpace>
    PreconditionChebyshev<Number, MemorySpace>::PreconditionChebyshev()
      : PreconditionIfpackBase<Number, MemorySpace>("CHEBYSHEV")
    {}



    template <typename Number, typename MemorySpace>
    void
    PreconditionChebyshev<Number, MemorySpace>::initialize(
      const SparseMatrix<Number, MemorySpace> &A,
      const AdditionalData                    &ad)
    {
      this->parameter_list.set("chebyshev: degree", ad.degree);
      this->parameter_list.set("chebyshev: max eigenvalue", ad.max_eigenvalue);
      this->parameter_list.set("chebyshev: min eigenvalue", ad.min_eigenvalue);
      this->parameter_list.set("chebyshev: ratio eigenvalue",
                               ad.eigenvalue_ratio);
      this->parameter_list.set("chebyshev: min diagonal value",
                               ad.min_diagonal);
      this->parameter_list.set("chebyshev: zero starting solution",
                               !ad.nonzero_starting);
      PreconditionIfpackBase<Number, MemorySpace>::initialize(A);
    }



    /* ---------------------- PreconditionILU ----------------------*/
    template <typename Number, typename MemorySpace>
    PreconditionILU<Number, MemorySpace>::AdditionalData::AdditionalData(
      const int    ilu_fill,
      const double ilu_atol,
      const double ilu_rtol,
      const int    overlap)
      : ilu_fill(ilu_fill)
      , ilu_atol(ilu_atol)
      , ilu_rtol(ilu_rtol)
      , overlap(overlap)
    {}



    template <typename Number, typename MemorySpace>
    PreconditionILU<Number, MemorySpace>::PreconditionILU()
      : PreconditionIfpackBase<Number, MemorySpace>("SCHWARZ")
    {}



    template <typename Number, typename MemorySpace>
    void
    PreconditionILU<Number, MemorySpace>::initialize(
      const SparseMatrix<Number, MemorySpace> &A,
      const AdditionalData                    &ad)
    {
      this->parameter_list.set("schwarz: subdomain solver name", "RILUK");
      {
        Teuchos::ParameterList &rilukParams =
          this->parameter_list.sublist("schwarz: subdomain solver parameters",
                                       false);
        rilukParams.set("fact: iluk level-of-fill", ad.ilu_fill);
        rilukParams.set("fact: absolute threshold", ad.ilu_atol);
        rilukParams.set("fact: relative threshold", ad.ilu_rtol);
      }
      this->parameter_list.set("schwarz: combine mode", "ADD");
      this->parameter_list.set("schwarz: overlap level", ad.overlap);
      PreconditionIfpackBase<Number, MemorySpace>::initialize(A);
    }



    /* ---------------------- PreconditionILUT ----------------------*/
    template <typename Number, typename MemorySpace>
    PreconditionILUT<Number, MemorySpace>::AdditionalData::AdditionalData(
      const double ilut_drop,
      const double ilut_fill,
      const double ilut_atol,
      const double ilut_rtol,
      const int    overlap)
      : ilut_drop(ilut_drop)
      , ilut_fill(ilut_fill)
      , ilut_atol(ilut_atol)
      , ilut_rtol(ilut_rtol)
      , overlap(overlap)
    {}



    template <typename Number, typename MemorySpace>
    PreconditionILUT<Number, MemorySpace>::PreconditionILUT()
      : PreconditionIfpackBase<Number, MemorySpace>("SCHWARZ")
    {}



    template <typename Number, typename MemorySpace>
    void
    PreconditionILUT<Number, MemorySpace>::initialize(
      const SparseMatrix<Number, MemorySpace> &A,
      const AdditionalData                    &ad)
    {
      this->parameter_list.set("schwarz: subdomain solver name", "ILUT");
      {
        Teuchos::ParameterList &ilutParams =
          this->parameter_list.sublist("schwarz: subdomain solver parameters",
                                       false);
        ilutParams.set("fact: drop tolerance", ad.ilut_drop);
        ilutParams.set("fact: ilut level-of-fill", ad.ilut_fill);
        ilutParams.set("fact: absolute threshold", ad.ilut_atol);
        ilutParams.set("fact: relative threshold", ad.ilut_rtol);
      }
      this->parameter_list.set("schwarz: combine mode", "ADD");
      this->parameter_list.set("schwarz: overlap level", ad.overlap);
      PreconditionIfpackBase<Number, MemorySpace>::initialize(A);
    }
  } // namespace TpetraWrappers
} // namespace LinearAlgebra

DEAL_II_NAMESPACE_CLOSE

#  endif // DEAL_II_TRILINOS_WITH_IFPACK2
#else

// Make sure the scripts that create the C++20 module input files have
// something to latch on if the preprocessor #ifdef above would
// otherwise lead to an empty content of the file.
DEAL_II_NAMESPACE_OPEN
DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_TRILINOS_WITH_TPETRA

#endif
