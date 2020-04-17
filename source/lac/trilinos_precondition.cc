// ---------------------------------------------------------------------
//
// Copyright (C) 2008 - 2019 by the deal.II authors
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

#include <deal.II/lac/trilinos_precondition.h>

#ifdef DEAL_II_WITH_TRILINOS

#  include <deal.II/lac/sparse_matrix.h>
#  include <deal.II/lac/trilinos_sparse_matrix.h>
#  include <deal.II/lac/vector.h>

#  include <Epetra_MultiVector.h>
#  include <Ifpack.h>
#  include <Ifpack_Chebyshev.h>
#  include <Teuchos_ParameterList.hpp>
#  include <Teuchos_RCP.hpp>

DEAL_II_NAMESPACE_OPEN

namespace TrilinosWrappers
{
  PreconditionBase::PreconditionBase()
#  ifdef DEAL_II_WITH_MPI
    : communicator(MPI_COMM_SELF)
#  endif
  {}



  PreconditionBase::PreconditionBase(const PreconditionBase &base)
    : Subscriptor()
    , preconditioner(base.preconditioner)
    ,
#  ifdef DEAL_II_WITH_MPI
    communicator(base.communicator)
    ,
#  endif
    vector_distributor(new Epetra_Map(*base.vector_distributor))
  {}



  void
  PreconditionBase::clear()
  {
    preconditioner.reset();
#  ifdef DEAL_II_WITH_MPI
    communicator = MPI_COMM_SELF;
#  endif
    vector_distributor.reset();
  }


  MPI_Comm
  PreconditionBase::get_mpi_communicator() const
  {
#  ifdef DEAL_II_WITH_MPI
    return communicator.Comm();
#  else
    return MPI_COMM_SELF;
#  endif
  }


  Epetra_Operator &
  PreconditionBase::trilinos_operator() const
  {
    AssertThrow(!preconditioner.is_null(),
                ExcMessage("Trying to dereference a null pointer."));
    return (*preconditioner);
  }


  IndexSet
  PreconditionBase::locally_owned_domain_indices() const
  {
    return IndexSet(preconditioner->OperatorDomainMap());
  }


  IndexSet
  PreconditionBase::locally_owned_range_indices() const
  {
    return IndexSet(preconditioner->OperatorRangeMap());
  }

  /* -------------------------- PreconditionJacobi -------------------------- */

  PreconditionJacobi::AdditionalData::AdditionalData(
    const double       omega,
    const double       min_diagonal,
    const unsigned int n_sweeps)
    : omega(omega)
    , min_diagonal(min_diagonal)
    , n_sweeps(n_sweeps)
  {}



  void
  PreconditionJacobi::initialize(const SparseMatrix &  matrix,
                                 const AdditionalData &additional_data)
  {
    // release memory before reallocation
    preconditioner.reset();
    preconditioner.reset(
      Ifpack().Create("point relaxation",
                      const_cast<Epetra_CrsMatrix *>(&matrix.trilinos_matrix()),
                      0));

    Ifpack_Preconditioner *ifpack =
      dynamic_cast<Ifpack_Preconditioner *>(preconditioner.get());
    Assert(ifpack != nullptr,
           ExcMessage("Trilinos could not create this "
                      "preconditioner"));

    int ierr;

    Teuchos::ParameterList parameter_list;
    parameter_list.set("relaxation: sweeps",
                       static_cast<int>(additional_data.n_sweeps));
    parameter_list.set("relaxation: type", "Jacobi");
    parameter_list.set("relaxation: damping factor", additional_data.omega);
    parameter_list.set("relaxation: min diagonal value",
                       additional_data.min_diagonal);

    ierr = ifpack->SetParameters(parameter_list);
    AssertThrow(ierr == 0, ExcTrilinosError(ierr));

    ierr = ifpack->Initialize();
    AssertThrow(ierr == 0, ExcTrilinosError(ierr));

    ierr = ifpack->Compute();
    AssertThrow(ierr == 0, ExcTrilinosError(ierr));
  }



  /* -------------------------- PreconditionSSOR -------------------------- */

  PreconditionSSOR::AdditionalData::AdditionalData(const double omega,
                                                   const double min_diagonal,
                                                   const unsigned int overlap,
                                                   const unsigned int n_sweeps)
    : omega(omega)
    , min_diagonal(min_diagonal)
    , overlap(overlap)
    , n_sweeps(n_sweeps)
  {}



  void
  PreconditionSSOR::initialize(const SparseMatrix &  matrix,
                               const AdditionalData &additional_data)
  {
    preconditioner.reset();
    preconditioner.reset(
      Ifpack().Create("point relaxation",
                      const_cast<Epetra_CrsMatrix *>(&matrix.trilinos_matrix()),
                      additional_data.overlap));

    Ifpack_Preconditioner *ifpack =
      dynamic_cast<Ifpack_Preconditioner *>(preconditioner.get());
    Assert(ifpack != nullptr,
           ExcMessage("Trilinos could not create this "
                      "preconditioner"));

    int ierr;

    Teuchos::ParameterList parameter_list;
    parameter_list.set("relaxation: sweeps",
                       static_cast<int>(additional_data.n_sweeps));
    parameter_list.set("relaxation: type", "symmetric Gauss-Seidel");
    parameter_list.set("relaxation: damping factor", additional_data.omega);
    parameter_list.set("relaxation: min diagonal value",
                       additional_data.min_diagonal);
    parameter_list.set("schwarz: combine mode", "Add");

    ierr = ifpack->SetParameters(parameter_list);
    AssertThrow(ierr == 0, ExcTrilinosError(ierr));

    ierr = ifpack->Initialize();
    AssertThrow(ierr == 0, ExcTrilinosError(ierr));

    ierr = ifpack->Compute();
    AssertThrow(ierr == 0, ExcTrilinosError(ierr));
  }



  /* -------------------------- PreconditionSOR -------------------------- */

  PreconditionSOR::AdditionalData::AdditionalData(const double omega,
                                                  const double min_diagonal,
                                                  const unsigned int overlap,
                                                  const unsigned int n_sweeps)
    : omega(omega)
    , min_diagonal(min_diagonal)
    , overlap(overlap)
    , n_sweeps(n_sweeps)
  {}



  void
  PreconditionSOR::initialize(const SparseMatrix &  matrix,
                              const AdditionalData &additional_data)
  {
    preconditioner.reset();
    preconditioner.reset(
      Ifpack().Create("point relaxation",
                      const_cast<Epetra_CrsMatrix *>(&matrix.trilinos_matrix()),
                      additional_data.overlap));

    Ifpack_Preconditioner *ifpack =
      dynamic_cast<Ifpack_Preconditioner *>(preconditioner.get());
    Assert(ifpack != nullptr,
           ExcMessage("Trilinos could not create this "
                      "preconditioner"));

    int ierr;

    Teuchos::ParameterList parameter_list;
    parameter_list.set("relaxation: sweeps",
                       static_cast<int>(additional_data.n_sweeps));
    parameter_list.set("relaxation: type", "Gauss-Seidel");
    parameter_list.set("relaxation: damping factor", additional_data.omega);
    parameter_list.set("relaxation: min diagonal value",
                       additional_data.min_diagonal);
    parameter_list.set("schwarz: combine mode", "Add");

    ierr = ifpack->SetParameters(parameter_list);
    AssertThrow(ierr == 0, ExcTrilinosError(ierr));

    ierr = ifpack->Initialize();
    AssertThrow(ierr == 0, ExcTrilinosError(ierr));

    ierr = ifpack->Compute();
    AssertThrow(ierr == 0, ExcTrilinosError(ierr));
  }



  /* ----------------------- PreconditionBlockJacobi ---------------------- */

  PreconditionBlockJacobi::AdditionalData::AdditionalData(
    const unsigned int block_size,
    const std::string &block_creation_type,
    const double       omega,
    const double       min_diagonal,
    const unsigned int n_sweeps)
    : block_size(block_size)
    , block_creation_type(block_creation_type)
    , omega(omega)
    , min_diagonal(min_diagonal)
    , n_sweeps(n_sweeps)
  {}



  void
  PreconditionBlockJacobi::initialize(const SparseMatrix &  matrix,
                                      const AdditionalData &additional_data)
  {
    // release memory before reallocation
    preconditioner.reset();

    // Block relaxation setup fails if we have no locally owned rows. As a
    // work-around we just pretend to use point relaxation on those processors:
    preconditioner.reset(
      Ifpack().Create((matrix.trilinos_matrix().NumMyRows() == 0) ?
                        "point relaxation" :
                        "block relaxation",
                      const_cast<Epetra_CrsMatrix *>(&matrix.trilinos_matrix()),
                      0));

    Ifpack_Preconditioner *ifpack =
      dynamic_cast<Ifpack_Preconditioner *>(preconditioner.get());
    Assert(ifpack != nullptr,
           ExcMessage("Trilinos could not create this "
                      "preconditioner"));

    int ierr;

    Teuchos::ParameterList parameter_list;
    parameter_list.set("relaxation: sweeps",
                       static_cast<int>(additional_data.n_sweeps));
    parameter_list.set("relaxation: type", "Jacobi");
    parameter_list.set("relaxation: damping factor", additional_data.omega);
    parameter_list.set("relaxation: min diagonal value",
                       additional_data.min_diagonal);
    parameter_list.set("partitioner: type",
                       additional_data.block_creation_type);
    int n_local_parts =
      (matrix.trilinos_matrix().NumMyRows() + additional_data.block_size - 1) /
      additional_data.block_size;
    parameter_list.set("partitioner: local parts", n_local_parts);

    ierr = ifpack->SetParameters(parameter_list);
    AssertThrow(ierr == 0, ExcTrilinosError(ierr));

    ierr = ifpack->Initialize();
    AssertThrow(ierr == 0, ExcTrilinosError(ierr));

    ierr = ifpack->Compute();
    AssertThrow(ierr == 0, ExcTrilinosError(ierr));
  }



  /* ----------------------- PreconditionBlockSSOR ------------------------ */

  PreconditionBlockSSOR::AdditionalData::AdditionalData(
    const unsigned int block_size,
    const std::string &block_creation_type,
    const double       omega,
    const double       min_diagonal,
    const unsigned int overlap,
    const unsigned int n_sweeps)
    : block_size(block_size)
    , block_creation_type(block_creation_type)
    , omega(omega)
    , min_diagonal(min_diagonal)
    , overlap(overlap)
    , n_sweeps(n_sweeps)
  {}



  void
  PreconditionBlockSSOR::initialize(const SparseMatrix &  matrix,
                                    const AdditionalData &additional_data)
  {
    preconditioner.reset();

    // Block relaxation setup fails if we have no locally owned rows. As a
    // work-around we just pretend to use point relaxation on those processors:
    preconditioner.reset(
      Ifpack().Create((matrix.trilinos_matrix().NumMyRows() == 0) ?
                        "point relaxation" :
                        "block relaxation",
                      const_cast<Epetra_CrsMatrix *>(&matrix.trilinos_matrix()),
                      additional_data.overlap));

    Ifpack_Preconditioner *ifpack =
      dynamic_cast<Ifpack_Preconditioner *>(preconditioner.get());
    Assert(ifpack != nullptr,
           ExcMessage("Trilinos could not create this "
                      "preconditioner"));

    int ierr;

    Teuchos::ParameterList parameter_list;
    parameter_list.set("relaxation: sweeps",
                       static_cast<int>(additional_data.n_sweeps));
    parameter_list.set("relaxation: type", "symmetric Gauss-Seidel");
    parameter_list.set("relaxation: damping factor", additional_data.omega);
    parameter_list.set("relaxation: min diagonal value",
                       additional_data.min_diagonal);
    parameter_list.set("schwarz: combine mode", "Add");
    parameter_list.set("partitioner: type",
                       additional_data.block_creation_type);
    int n_local_parts =
      (matrix.trilinos_matrix().NumMyRows() + additional_data.block_size - 1) /
      additional_data.block_size;
    parameter_list.set("partitioner: local parts", n_local_parts);

    ierr = ifpack->SetParameters(parameter_list);
    AssertThrow(ierr == 0, ExcTrilinosError(ierr));

    ierr = ifpack->Initialize();
    AssertThrow(ierr == 0, ExcTrilinosError(ierr));

    ierr = ifpack->Compute();
    AssertThrow(ierr == 0, ExcTrilinosError(ierr));
  }



  /* ------------------------ PreconditionBlockSOR ------------------------ */

  PreconditionBlockSOR::AdditionalData::AdditionalData(
    const unsigned int block_size,
    const std::string &block_creation_type,
    const double       omega,
    const double       min_diagonal,
    const unsigned int overlap,
    const unsigned int n_sweeps)
    : block_size(block_size)
    , block_creation_type(block_creation_type)
    , omega(omega)
    , min_diagonal(min_diagonal)
    , overlap(overlap)
    , n_sweeps(n_sweeps)
  {}



  void
  PreconditionBlockSOR::initialize(const SparseMatrix &  matrix,
                                   const AdditionalData &additional_data)
  {
    preconditioner.reset();

    // Block relaxation setup fails if we have no locally owned rows. As a
    // work-around we just pretend to use point relaxation on those processors:
    preconditioner.reset(
      Ifpack().Create((matrix.trilinos_matrix().NumMyRows() == 0) ?
                        "point relaxation" :
                        "block relaxation",
                      const_cast<Epetra_CrsMatrix *>(&matrix.trilinos_matrix()),
                      additional_data.overlap));

    Ifpack_Preconditioner *ifpack =
      dynamic_cast<Ifpack_Preconditioner *>(preconditioner.get());
    Assert(ifpack != nullptr,
           ExcMessage("Trilinos could not create this "
                      "preconditioner"));

    int ierr;

    Teuchos::ParameterList parameter_list;
    parameter_list.set("relaxation: sweeps",
                       static_cast<int>(additional_data.n_sweeps));
    parameter_list.set("relaxation: type", "Gauss-Seidel");
    parameter_list.set("relaxation: damping factor", additional_data.omega);
    parameter_list.set("relaxation: min diagonal value",
                       additional_data.min_diagonal);
    parameter_list.set("schwarz: combine mode", "Add");
    parameter_list.set("partitioner: type",
                       additional_data.block_creation_type);
    int n_local_parts =
      (matrix.trilinos_matrix().NumMyRows() + additional_data.block_size - 1) /
      additional_data.block_size;
    parameter_list.set("partitioner: local parts", n_local_parts);

    ierr = ifpack->SetParameters(parameter_list);
    AssertThrow(ierr == 0, ExcTrilinosError(ierr));

    ierr = ifpack->Initialize();
    AssertThrow(ierr == 0, ExcTrilinosError(ierr));

    ierr = ifpack->Compute();
    AssertThrow(ierr == 0, ExcTrilinosError(ierr));
  }



  /* -------------------------- PreconditionIC -------------------------- */

  PreconditionIC::AdditionalData::AdditionalData(const unsigned int ic_fill,
                                                 const double       ic_atol,
                                                 const double       ic_rtol,
                                                 const unsigned int overlap)
    : ic_fill(ic_fill)
    , ic_atol(ic_atol)
    , ic_rtol(ic_rtol)
    , overlap(overlap)
  {}



  void
  PreconditionIC::initialize(const SparseMatrix &  matrix,
                             const AdditionalData &additional_data)
  {
    preconditioner.reset();
    preconditioner.reset(
      Ifpack().Create("IC",
                      const_cast<Epetra_CrsMatrix *>(&matrix.trilinos_matrix()),
                      additional_data.overlap));

    Ifpack_Preconditioner *ifpack =
      dynamic_cast<Ifpack_Preconditioner *>(preconditioner.get());
    Assert(ifpack != nullptr,
           ExcMessage("Trilinos could not create this "
                      "preconditioner"));

    int ierr;

    Teuchos::ParameterList parameter_list;
    parameter_list.set("fact: level-of-fill", additional_data.ic_fill);
    parameter_list.set("fact: absolute threshold", additional_data.ic_atol);
    parameter_list.set("fact: relative threshold", additional_data.ic_rtol);
    parameter_list.set("schwarz: combine mode", "Add");

    ierr = ifpack->SetParameters(parameter_list);
    AssertThrow(ierr == 0, ExcTrilinosError(ierr));

    ierr = ifpack->Initialize();
    AssertThrow(ierr == 0, ExcTrilinosError(ierr));

    ierr = ifpack->Compute();
    AssertThrow(ierr == 0, ExcTrilinosError(ierr));
  }



  /* -------------------------- PreconditionILU -------------------------- */

  PreconditionILU::AdditionalData::AdditionalData(const unsigned int ilu_fill,
                                                  const double       ilu_atol,
                                                  const double       ilu_rtol,
                                                  const unsigned int overlap)
    : ilu_fill(ilu_fill)
    , ilu_atol(ilu_atol)
    , ilu_rtol(ilu_rtol)
    , overlap(overlap)
  {}



  void
  PreconditionILU::initialize(const SparseMatrix &  matrix,
                              const AdditionalData &additional_data)
  {
    preconditioner.reset();
    preconditioner.reset(
      Ifpack().Create("ILU",
                      const_cast<Epetra_CrsMatrix *>(&matrix.trilinos_matrix()),
                      additional_data.overlap));

    Ifpack_Preconditioner *ifpack =
      dynamic_cast<Ifpack_Preconditioner *>(preconditioner.get());
    Assert(ifpack != nullptr,
           ExcMessage("Trilinos could not create this "
                      "preconditioner"));

    int ierr;

    Teuchos::ParameterList parameter_list;
    parameter_list.set("fact: level-of-fill",
                       static_cast<int>(additional_data.ilu_fill));
    parameter_list.set("fact: absolute threshold", additional_data.ilu_atol);
    parameter_list.set("fact: relative threshold", additional_data.ilu_rtol);
    parameter_list.set("schwarz: combine mode", "Add");

    ierr = ifpack->SetParameters(parameter_list);
    AssertThrow(ierr == 0, ExcTrilinosError(ierr));

    ierr = ifpack->Initialize();
    AssertThrow(ierr == 0, ExcTrilinosError(ierr));

    ierr = ifpack->Compute();
    AssertThrow(ierr == 0, ExcTrilinosError(ierr));
  }



  /* -------------------------- PreconditionILUT -------------------------- */

  PreconditionILUT::AdditionalData::AdditionalData(const double       ilut_drop,
                                                   const unsigned int ilut_fill,
                                                   const double       ilut_atol,
                                                   const double       ilut_rtol,
                                                   const unsigned int overlap)
    : ilut_drop(ilut_drop)
    , ilut_fill(ilut_fill)
    , ilut_atol(ilut_atol)
    , ilut_rtol(ilut_rtol)
    , overlap(overlap)
  {}



  void
  PreconditionILUT::initialize(const SparseMatrix &  matrix,
                               const AdditionalData &additional_data)
  {
    preconditioner.reset();
    preconditioner.reset(
      Ifpack().Create("ILUT",
                      const_cast<Epetra_CrsMatrix *>(&matrix.trilinos_matrix()),
                      additional_data.overlap));

    Ifpack_Preconditioner *ifpack =
      dynamic_cast<Ifpack_Preconditioner *>(preconditioner.get());
    Assert(ifpack != nullptr,
           ExcMessage("Trilinos could not create this "
                      "preconditioner"));

    int ierr;

    Teuchos::ParameterList parameter_list;
    parameter_list.set("fact: drop value", additional_data.ilut_drop);
    parameter_list.set("fact: level-of-fill",
                       static_cast<int>(additional_data.ilut_fill));
    parameter_list.set("fact: absolute threshold", additional_data.ilut_atol);
    parameter_list.set("fact: relative threshold", additional_data.ilut_rtol);
    parameter_list.set("schwarz: combine mode", "Add");

    ierr = ifpack->SetParameters(parameter_list);
    AssertThrow(ierr == 0, ExcTrilinosError(ierr));

    ierr = ifpack->Initialize();
    AssertThrow(ierr == 0, ExcTrilinosError(ierr));

    ierr = ifpack->Compute();
    AssertThrow(ierr == 0, ExcTrilinosError(ierr));
  }



  /* ---------------------- PreconditionBlockDirect --------------------- */

  PreconditionBlockwiseDirect::AdditionalData::AdditionalData(
    const unsigned int overlap)
    : overlap(overlap)
  {}



  void
  PreconditionBlockwiseDirect::initialize(const SparseMatrix &  matrix,
                                          const AdditionalData &additional_data)
  {
    preconditioner.reset();
    preconditioner.reset(
      Ifpack().Create("Amesos",
                      const_cast<Epetra_CrsMatrix *>(&matrix.trilinos_matrix()),
                      additional_data.overlap));

    Ifpack_Preconditioner *ifpack =
      dynamic_cast<Ifpack_Preconditioner *>(preconditioner.get());
    Assert(ifpack != nullptr,
           ExcMessage("Trilinos could not create this "
                      "preconditioner"));

    int ierr;

    Teuchos::ParameterList parameter_list;
    parameter_list.set("schwarz: combine mode", "Add");

    ierr = ifpack->SetParameters(parameter_list);
    AssertThrow(ierr == 0, ExcTrilinosError(ierr));

    ierr = ifpack->Initialize();
    AssertThrow(ierr == 0, ExcTrilinosError(ierr));

    ierr = ifpack->Compute();
    AssertThrow(ierr == 0, ExcTrilinosError(ierr));
  }



  /* ---------------------- PreconditionBlockDirect --------------------- */

  PreconditionChebyshev::AdditionalData::AdditionalData(
    const unsigned int degree,
    const double       max_eigenvalue,
    const double       eigenvalue_ratio,
    const double       min_eigenvalue,
    const double       min_diagonal,
    const bool         nonzero_starting)
    : degree(degree)
    , max_eigenvalue(max_eigenvalue)
    , eigenvalue_ratio(eigenvalue_ratio)
    , min_eigenvalue(min_eigenvalue)
    , min_diagonal(min_diagonal)
    , nonzero_starting(nonzero_starting)
  {}



  void
  PreconditionChebyshev::initialize(const SparseMatrix &  matrix,
                                    const AdditionalData &additional_data)
  {
    preconditioner =
      Teuchos::rcp(new Ifpack_Chebyshev(&matrix.trilinos_matrix()));

    Ifpack_Chebyshev *ifpack =
      dynamic_cast<Ifpack_Chebyshev *>(preconditioner.get());
    Assert(ifpack != nullptr,
           ExcMessage("Trilinos could not create this "
                      "preconditioner"));

    int ierr;

    Teuchos::ParameterList parameter_list;
    parameter_list.set("chebyshev: ratio eigenvalue",
                       additional_data.eigenvalue_ratio);
    parameter_list.set("chebyshev: min eigenvalue",
                       additional_data.min_eigenvalue);
    parameter_list.set("chebyshev: max eigenvalue",
                       additional_data.max_eigenvalue);
    parameter_list.set("chebyshev: degree",
                       static_cast<int>(additional_data.degree));
    parameter_list.set("chebyshev: min diagonal value",
                       additional_data.min_diagonal);
    parameter_list.set("chebyshev: zero starting solution",
                       !additional_data.nonzero_starting);

    ierr = ifpack->SetParameters(parameter_list);
    AssertThrow(ierr == 0, ExcTrilinosError(ierr));

    ierr = ifpack->Initialize();
    AssertThrow(ierr == 0, ExcTrilinosError(ierr));

    ierr = ifpack->Compute();
    AssertThrow(ierr == 0, ExcTrilinosError(ierr));
  }



  /* -------------------------- PreconditionIdentity --------------------- */

  void
  PreconditionIdentity::initialize(const SparseMatrix &matrix,
                                   const AdditionalData &)
  {
    // What follows just configures a dummy preconditioner that
    // sets up the domain and range maps, as well as the communicator.
    // It is never used as the vmult, Tvmult operations are
    // given a custom definition.
    // Note: This is only required in order to wrap this
    // preconditioner in a LinearOperator without an exemplar
    // matrix.

    // From PreconditionJacobi:
    // release memory before reallocation
    preconditioner.reset();
    preconditioner.reset(
      Ifpack().Create("point relaxation",
                      const_cast<Epetra_CrsMatrix *>(&matrix.trilinos_matrix()),
                      0));

    Ifpack_Preconditioner *ifpack =
      dynamic_cast<Ifpack_Preconditioner *>(preconditioner.get());
    Assert(ifpack != nullptr,
           ExcMessage("Trilinos could not create this "
                      "preconditioner"));

    int ierr;

    Teuchos::ParameterList parameter_list;
    parameter_list.set("relaxation: sweeps", 1);
    parameter_list.set("relaxation: type", "Jacobi");
    parameter_list.set("relaxation: damping factor", 1.0);
    parameter_list.set("relaxation: min diagonal value", 0.0);

    ierr = ifpack->SetParameters(parameter_list);
    AssertThrow(ierr == 0, ExcTrilinosError(ierr));

    ierr = ifpack->Initialize();
    AssertThrow(ierr == 0, ExcTrilinosError(ierr));

    ierr = ifpack->Compute();
    AssertThrow(ierr == 0, ExcTrilinosError(ierr));
  }

  void
  PreconditionIdentity::vmult(MPI::Vector &dst, const MPI::Vector &src) const
  {
    dst = src;
  }

  void
  PreconditionIdentity::Tvmult(MPI::Vector &dst, const MPI::Vector &src) const
  {
    dst = src;
  }

  void
  PreconditionIdentity::vmult(dealii::Vector<double> &      dst,
                              const dealii::Vector<double> &src) const
  {
    dst = src;
  }

  void
  PreconditionIdentity::Tvmult(dealii::Vector<double> &      dst,
                               const dealii::Vector<double> &src) const
  {
    dst = src;
  }

#  ifndef DOXYGEN
  void
  PreconditionIdentity::vmult(
    LinearAlgebra::distributed::Vector<double> &      dst,
    const LinearAlgebra::distributed::Vector<double> &src) const
  {
    dst = src;
  }

  void
  PreconditionIdentity::Tvmult(
    LinearAlgebra::distributed::Vector<double> &      dst,
    const LinearAlgebra::distributed::Vector<double> &src) const
  {
    dst = src;
  }
#  endif // DOXYGEN
} // namespace TrilinosWrappers

DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_WITH_TRILINOS
