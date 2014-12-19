// ---------------------------------------------------------------------
//
// Copyright (C) 2008 - 2014 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------

#include <deal.II/lac/trilinos_precondition.h>

#ifdef DEAL_II_WITH_TRILINOS

#  include <deal.II/lac/vector.h>
#  include <deal.II/lac/sparse_matrix.h>
#  include <deal.II/lac/trilinos_sparse_matrix.h>

#  include <Ifpack.h>
#  include <Ifpack_Chebyshev.h>
#  include <Teuchos_ParameterList.hpp>
#  include <Epetra_MultiVector.h>
#  include <ml_include.h>
#  include <ml_MultiLevelPreconditioner.h>


DEAL_II_NAMESPACE_OPEN

namespace TrilinosWrappers
{
  namespace
  {
#ifndef DEAL_II_WITH_64BIT_INDICES
    int n_global_rows (const Epetra_RowMatrix &matrix)
    {
      return matrix.NumGlobalRows();
    }

    int global_length (const Epetra_MultiVector &vector)
    {
      return vector.GlobalLength();
    }

    int gid(const Epetra_Map &map, unsigned int i)
    {
      return map.GID(i);
    }
#else
    long long int n_global_rows (const Epetra_RowMatrix &matrix)
    {
      return matrix.NumGlobalRows64();
    }

    long long int global_length (const Epetra_MultiVector &vector)
    {
      return vector.GlobalLength64();
    }

    long long int gid(const Epetra_Map &map, dealii::types::global_dof_index i)
    {
      return map.GID64(i);
    }
#endif
  }

  PreconditionBase::PreconditionBase()
#ifdef DEAL_II_WITH_MPI
    :
    communicator (MPI_COMM_SELF)
#endif
  {}



  PreconditionBase::PreconditionBase(const PreconditionBase &base)
    :
    Subscriptor (),
    preconditioner (base.preconditioner),
#ifdef DEAL_II_WITH_MPI
    communicator (base.communicator),
#endif
    vector_distributor (new Epetra_Map(*base.vector_distributor))
  {}



  PreconditionBase::~PreconditionBase()
  {}



  void PreconditionBase::clear ()
  {
    preconditioner.reset();
#ifdef DEAL_II_WITH_MPI
    communicator = MPI_COMM_SELF;
#endif
    vector_distributor.reset();
  }



  /* -------------------------- PreconditionJacobi -------------------------- */

  PreconditionJacobi::AdditionalData::
  AdditionalData (const double omega,
                  const double min_diagonal,
                  const unsigned int n_sweeps)
    :
    omega (omega),
    min_diagonal (min_diagonal),
    n_sweeps     (n_sweeps)
  {}



  void
  PreconditionJacobi::initialize (const SparseMatrix   &matrix,
                                  const AdditionalData &additional_data)
  {
    // release memory before reallocation
    preconditioner.reset ();
    preconditioner.reset (Ifpack().Create
                          ("point relaxation",
                           const_cast<Epetra_CrsMatrix *>(&matrix.trilinos_matrix()),
                           0));

    Ifpack_Preconditioner *ifpack = static_cast<Ifpack_Preconditioner *>
                                    (preconditioner.get());
    Assert (ifpack != 0, ExcMessage ("Trilinos could not create this "
                                     "preconditioner"));

    int ierr;

    Teuchos::ParameterList parameter_list;
    parameter_list.set ("relaxation: sweeps", static_cast<int>(additional_data.n_sweeps));
    parameter_list.set ("relaxation: type", "Jacobi");
    parameter_list.set ("relaxation: damping factor", additional_data.omega);
    parameter_list.set ("relaxation: min diagonal value",
                        additional_data.min_diagonal);

    ierr = ifpack->SetParameters(parameter_list);
    AssertThrow (ierr == 0, ExcTrilinosError(ierr));

    ierr = ifpack->Initialize();
    AssertThrow (ierr == 0, ExcTrilinosError(ierr));

    ierr = ifpack->Compute();
    AssertThrow (ierr == 0, ExcTrilinosError(ierr));
  }



  /* -------------------------- PreconditionSSOR -------------------------- */

  PreconditionSSOR::AdditionalData::
  AdditionalData (const double       omega,
                  const double       min_diagonal,
                  const unsigned int overlap,
                  const unsigned int n_sweeps)
    :
    omega        (omega),
    min_diagonal (min_diagonal),
    overlap      (overlap),
    n_sweeps     (n_sweeps)
  {}



  void
  PreconditionSSOR::initialize (const SparseMatrix   &matrix,
                                const AdditionalData &additional_data)
  {
    preconditioner.reset ();
    preconditioner.reset (Ifpack().Create
                          ("point relaxation",
                           const_cast<Epetra_CrsMatrix *>(&matrix.trilinos_matrix()),
                           additional_data.overlap));

    Ifpack_Preconditioner *ifpack = static_cast<Ifpack_Preconditioner *>
                                    (preconditioner.get());
    Assert (ifpack != 0, ExcMessage ("Trilinos could not create this "
                                     "preconditioner"));

    int ierr;

    Teuchos::ParameterList parameter_list;
    parameter_list.set ("relaxation: sweeps", static_cast<int>(additional_data.n_sweeps));
    parameter_list.set ("relaxation: type", "symmetric Gauss-Seidel");
    parameter_list.set ("relaxation: damping factor", additional_data.omega);
    parameter_list.set ("relaxation: min diagonal value",
                        additional_data.min_diagonal);
    parameter_list.set ("schwarz: combine mode", "Add");

    ierr = ifpack->SetParameters(parameter_list);
    AssertThrow (ierr == 0, ExcTrilinosError(ierr));

    ierr = ifpack->Initialize();
    AssertThrow (ierr == 0, ExcTrilinosError(ierr));

    ierr = ifpack->Compute();
    AssertThrow (ierr == 0, ExcTrilinosError(ierr));
  }



  /* -------------------------- PreconditionSOR -------------------------- */

  PreconditionSOR::AdditionalData::
  AdditionalData (const double       omega,
                  const double       min_diagonal,
                  const unsigned int overlap,
                  const unsigned int n_sweeps)
    :
    omega        (omega),
    min_diagonal (min_diagonal),
    overlap      (overlap),
    n_sweeps     (n_sweeps)
  {}



  void
  PreconditionSOR::initialize (const SparseMatrix   &matrix,
                               const AdditionalData &additional_data)
  {
    preconditioner.reset ();
    preconditioner.reset (Ifpack().Create
                          ("point relaxation",
                           const_cast<Epetra_CrsMatrix *>(&matrix.trilinos_matrix()),
                           additional_data.overlap));

    Ifpack_Preconditioner *ifpack = static_cast<Ifpack_Preconditioner *>
                                    (preconditioner.get());
    Assert (ifpack != 0, ExcMessage ("Trilinos could not create this "
                                     "preconditioner"));

    int ierr;

    Teuchos::ParameterList parameter_list;
    parameter_list.set ("relaxation: sweeps", static_cast<int>(additional_data.n_sweeps));
    parameter_list.set ("relaxation: type", "Gauss-Seidel");
    parameter_list.set ("relaxation: damping factor", additional_data.omega);
    parameter_list.set ("relaxation: min diagonal value",
                        additional_data.min_diagonal);
    parameter_list.set ("schwarz: combine mode", "Add");

    ierr = ifpack->SetParameters(parameter_list);
    AssertThrow (ierr == 0, ExcTrilinosError(ierr));

    ierr = ifpack->Initialize();
    AssertThrow (ierr == 0, ExcTrilinosError(ierr));

    ierr = ifpack->Compute();
    AssertThrow (ierr == 0, ExcTrilinosError(ierr));
  }



  /* ----------------------- PreconditionBlockJacobi ---------------------- */

  PreconditionBlockJacobi::AdditionalData::
  AdditionalData (const unsigned int block_size,
                  const std::string  block_creation_type,
                  const double omega,
                  const double min_diagonal,
                  const unsigned int n_sweeps)
    :
    block_size(block_size),
    block_creation_type(block_creation_type),
    omega (omega),
    min_diagonal (min_diagonal),
    n_sweeps     (n_sweeps)
  {}



  void
  PreconditionBlockJacobi::initialize (const SparseMatrix   &matrix,
                                       const AdditionalData &additional_data)
  {
    // release memory before reallocation
    preconditioner.reset ();
    preconditioner.reset (Ifpack().Create
                          ("block relaxation",
                           const_cast<Epetra_CrsMatrix *>(&matrix.trilinos_matrix()),
                           0));

    Ifpack_Preconditioner *ifpack = static_cast<Ifpack_Preconditioner *>
                                    (preconditioner.get());
    Assert (ifpack != 0, ExcMessage ("Trilinos could not create this "
                                     "preconditioner"));

    int ierr;

    Teuchos::ParameterList parameter_list;
    parameter_list.set ("relaxation: sweeps", static_cast<int>(additional_data.n_sweeps));
    parameter_list.set ("relaxation: type", "Jacobi");
    parameter_list.set ("relaxation: damping factor", additional_data.omega);
    parameter_list.set ("relaxation: min diagonal value",
                        additional_data.min_diagonal);
    parameter_list.set ("partitioner: type", additional_data.block_creation_type);
    int n_local_parts = (matrix.trilinos_matrix().NumMyRows()+additional_data.
                         block_size-1)/additional_data.block_size;
    parameter_list.set ("partitioner: local parts", n_local_parts);

    ierr = ifpack->SetParameters(parameter_list);
    AssertThrow (ierr == 0, ExcTrilinosError(ierr));

    ierr = ifpack->Initialize();
    AssertThrow (ierr == 0, ExcTrilinosError(ierr));

    ierr = ifpack->Compute();
    AssertThrow (ierr == 0, ExcTrilinosError(ierr));
  }



  /* ----------------------- PreconditionBlockSSOR ------------------------ */

  PreconditionBlockSSOR::AdditionalData::
  AdditionalData (const unsigned int block_size,
                  const std::string  block_creation_type,
                  const double       omega,
                  const double       min_diagonal,
                  const unsigned int overlap,
                  const unsigned int n_sweeps)
    :
    block_size(block_size),
    block_creation_type(block_creation_type),
    omega        (omega),
    min_diagonal (min_diagonal),
    overlap      (overlap),
    n_sweeps     (n_sweeps)
  {}



  void
  PreconditionBlockSSOR::initialize (const SparseMatrix   &matrix,
                                     const AdditionalData &additional_data)
  {
    preconditioner.reset ();
    preconditioner.reset (Ifpack().Create
                          ("block relaxation",
                           const_cast<Epetra_CrsMatrix *>(&matrix.trilinos_matrix()),
                           additional_data.overlap));

    Ifpack_Preconditioner *ifpack = static_cast<Ifpack_Preconditioner *>
                                    (preconditioner.get());
    Assert (ifpack != 0, ExcMessage ("Trilinos could not create this "
                                     "preconditioner"));

    int ierr;

    Teuchos::ParameterList parameter_list;
    parameter_list.set ("relaxation: sweeps", static_cast<int>(additional_data.n_sweeps));
    parameter_list.set ("relaxation: type", "symmetric Gauss-Seidel");
    parameter_list.set ("relaxation: damping factor", additional_data.omega);
    parameter_list.set ("relaxation: min diagonal value",
                        additional_data.min_diagonal);
    parameter_list.set ("schwarz: combine mode", "Add");
    parameter_list.set ("partitioner: type", additional_data.block_creation_type);
    int n_local_parts = (matrix.trilinos_matrix().NumMyRows()+additional_data.
                         block_size-1)/additional_data.block_size;
    parameter_list.set ("partitioner: local parts", n_local_parts);

    ierr = ifpack->SetParameters(parameter_list);
    AssertThrow (ierr == 0, ExcTrilinosError(ierr));

    ierr = ifpack->Initialize();
    AssertThrow (ierr == 0, ExcTrilinosError(ierr));

    ierr = ifpack->Compute();
    AssertThrow (ierr == 0, ExcTrilinosError(ierr));
  }



  /* ------------------------ PreconditionBlockSOR ------------------------ */

  PreconditionBlockSOR::AdditionalData::
  AdditionalData (const unsigned int block_size,
                  const std::string  block_creation_type,
                  const double       omega,
                  const double       min_diagonal,
                  const unsigned int overlap,
                  const unsigned int n_sweeps)
    :
    block_size(block_size),
    block_creation_type(block_creation_type),
    omega        (omega),
    min_diagonal (min_diagonal),
    overlap      (overlap),
    n_sweeps     (n_sweeps)
  {}



  void
  PreconditionBlockSOR::initialize (const SparseMatrix   &matrix,
                                    const AdditionalData &additional_data)
  {
    preconditioner.reset ();
    preconditioner.reset (Ifpack().Create
                          ("block relaxation",
                           const_cast<Epetra_CrsMatrix *>(&matrix.trilinos_matrix()),
                           additional_data.overlap));

    Ifpack_Preconditioner *ifpack = static_cast<Ifpack_Preconditioner *>
                                    (preconditioner.get());
    Assert (ifpack != 0, ExcMessage ("Trilinos could not create this "
                                     "preconditioner"));

    int ierr;

    Teuchos::ParameterList parameter_list;
    parameter_list.set ("relaxation: sweeps", static_cast<int>(additional_data.n_sweeps));
    parameter_list.set ("relaxation: type", "Gauss-Seidel");
    parameter_list.set ("relaxation: damping factor", additional_data.omega);
    parameter_list.set ("relaxation: min diagonal value",
                        additional_data.min_diagonal);
    parameter_list.set ("schwarz: combine mode", "Add");
    parameter_list.set ("partitioner: type", additional_data.block_creation_type);
    int n_local_parts = (matrix.trilinos_matrix().NumMyRows()+additional_data.
                         block_size-1)/additional_data.block_size;
    parameter_list.set ("partitioner: local parts", n_local_parts);

    ierr = ifpack->SetParameters(parameter_list);
    AssertThrow (ierr == 0, ExcTrilinosError(ierr));

    ierr = ifpack->Initialize();
    AssertThrow (ierr == 0, ExcTrilinosError(ierr));

    ierr = ifpack->Compute();
    AssertThrow (ierr == 0, ExcTrilinosError(ierr));
  }



  /* -------------------------- PreconditionIC -------------------------- */

  PreconditionIC::AdditionalData::
  AdditionalData (const unsigned int ic_fill,
                  const double       ic_atol,
                  const double       ic_rtol,
                  const unsigned int overlap)
    :
    ic_fill (ic_fill),
    ic_atol (ic_atol),
    ic_rtol (ic_rtol),
    overlap (overlap)
  {}



  void
  PreconditionIC::initialize (const SparseMatrix   &matrix,
                              const AdditionalData &additional_data)
  {
    preconditioner.reset ();
    preconditioner.reset (Ifpack().Create
                          ("IC",
                           const_cast<Epetra_CrsMatrix *>(&matrix.trilinos_matrix()),
                           additional_data.overlap));

    Ifpack_Preconditioner *ifpack = static_cast<Ifpack_Preconditioner *>
                                    (preconditioner.get());
    Assert (ifpack != 0, ExcMessage ("Trilinos could not create this "
                                     "preconditioner"));

    int ierr;

    Teuchos::ParameterList parameter_list;
    parameter_list.set ("fact: level-of-fill",(int)additional_data.ic_fill);
    parameter_list.set ("fact: absolute threshold",additional_data.ic_atol);
    parameter_list.set ("fact: relative threshold",additional_data.ic_rtol);
    parameter_list.set ("schwarz: combine mode", "Add");

    ierr = ifpack->SetParameters(parameter_list);
    AssertThrow (ierr == 0, ExcTrilinosError(ierr));

    ierr = ifpack->Initialize();
    AssertThrow (ierr == 0, ExcTrilinosError(ierr));

    ierr = ifpack->Compute();
    AssertThrow (ierr == 0, ExcTrilinosError(ierr));
  }



  /* -------------------------- PreconditionILU -------------------------- */

  PreconditionILU::AdditionalData::
  AdditionalData (const unsigned int ilu_fill,
                  const double       ilu_atol,
                  const double       ilu_rtol,
                  const unsigned int overlap)
    :
    ilu_fill (ilu_fill),
    ilu_atol (ilu_atol),
    ilu_rtol (ilu_rtol),
    overlap  (overlap)
  {}



  void
  PreconditionILU::initialize (const SparseMatrix   &matrix,
                               const AdditionalData &additional_data)
  {
    preconditioner.reset ();
    preconditioner.reset (Ifpack().Create
                          ("ILU",
                           const_cast<Epetra_CrsMatrix *>(&matrix.trilinos_matrix()),
                           additional_data.overlap));

    Ifpack_Preconditioner *ifpack = static_cast<Ifpack_Preconditioner *>
                                    (preconditioner.get());
    Assert (ifpack != 0, ExcMessage ("Trilinos could not create this "
                                     "preconditioner"));

    int ierr;

    Teuchos::ParameterList parameter_list;
    parameter_list.set ("fact: level-of-fill",(int)additional_data.ilu_fill);
    parameter_list.set ("fact: absolute threshold",additional_data.ilu_atol);
    parameter_list.set ("fact: relative threshold",additional_data.ilu_rtol);
    parameter_list.set ("schwarz: combine mode", "Add");

    ierr = ifpack->SetParameters(parameter_list);
    AssertThrow (ierr == 0, ExcTrilinosError(ierr));

    ierr = ifpack->Initialize();
    AssertThrow (ierr == 0, ExcTrilinosError(ierr));

    ierr = ifpack->Compute();
    AssertThrow (ierr == 0, ExcTrilinosError(ierr));
  }



  /* -------------------------- PreconditionILUT -------------------------- */

  PreconditionILUT::AdditionalData::
  AdditionalData (const double       ilut_drop,
                  const unsigned int ilut_fill,
                  const double       ilut_atol,
                  const double       ilut_rtol,
                  const unsigned int overlap)
    :
    ilut_drop (ilut_drop),
    ilut_fill (ilut_fill),
    ilut_atol (ilut_atol),
    ilut_rtol (ilut_rtol),
    overlap  (overlap)
  {}



  void
  PreconditionILUT::initialize (const SparseMatrix   &matrix,
                                const AdditionalData &additional_data)
  {
    preconditioner.reset ();
    preconditioner.reset (Ifpack().Create
                          ("ILUT",
                           const_cast<Epetra_CrsMatrix *>(&matrix.trilinos_matrix()),
                           additional_data.overlap));

    Ifpack_Preconditioner *ifpack = static_cast<Ifpack_Preconditioner *>
                                    (preconditioner.get());
    Assert (ifpack != 0, ExcMessage ("Trilinos could not create this "
                                     "preconditioner"));

    int ierr;

    Teuchos::ParameterList parameter_list;
    parameter_list.set ("fact: drop value",additional_data.ilut_drop);
    parameter_list.set ("fact: level-of-fill",(int)additional_data.ilut_fill);
    parameter_list.set ("fact: absolute threshold",additional_data.ilut_atol);
    parameter_list.set ("fact: relative threshold",additional_data.ilut_rtol);
    parameter_list.set ("schwarz: combine mode", "Add");

    ierr = ifpack->SetParameters(parameter_list);
    AssertThrow (ierr == 0, ExcTrilinosError(ierr));

    ierr = ifpack->Initialize();
    AssertThrow (ierr == 0, ExcTrilinosError(ierr));

    ierr = ifpack->Compute();
    AssertThrow (ierr == 0, ExcTrilinosError(ierr));
  }



  /* ---------------------- PreconditionBlockDirect --------------------- */

  PreconditionBlockwiseDirect::AdditionalData::
  AdditionalData (const unsigned int overlap)
    :
    overlap  (overlap)
  {}



  void
  PreconditionBlockwiseDirect::initialize (const SparseMatrix   &matrix,
                                           const AdditionalData &additional_data)
  {
    preconditioner.reset ();
    preconditioner.reset (Ifpack().Create
                          ("Amesos",
                           const_cast<Epetra_CrsMatrix *>(&matrix.trilinos_matrix()),
                           additional_data.overlap));

    Ifpack_Preconditioner *ifpack = static_cast<Ifpack_Preconditioner *>
                                    (preconditioner.get());
    Assert (ifpack != 0, ExcMessage ("Trilinos could not create this "
                                     "preconditioner"));

    int ierr;

    Teuchos::ParameterList parameter_list;
    parameter_list.set ("schwarz: combine mode", "Add");

    ierr = ifpack->SetParameters(parameter_list);
    AssertThrow (ierr == 0, ExcTrilinosError(ierr));

    ierr = ifpack->Initialize();
    AssertThrow (ierr == 0, ExcTrilinosError(ierr));

    ierr = ifpack->Compute();
    AssertThrow (ierr == 0, ExcTrilinosError(ierr));
  }



  /* ---------------------- PreconditionBlockDirect --------------------- */

  PreconditionChebyshev::AdditionalData::
  AdditionalData (const unsigned int degree,
                  const double       max_eigenvalue,
                  const double       eigenvalue_ratio,
                  const double       min_eigenvalue,
                  const double       min_diagonal,
                  const bool         nonzero_starting)
    :
    degree  (degree),
    max_eigenvalue (max_eigenvalue),
    eigenvalue_ratio (eigenvalue_ratio),
    min_eigenvalue (min_eigenvalue),
    min_diagonal (min_diagonal),
    nonzero_starting (nonzero_starting)
  {}



  void
  PreconditionChebyshev::initialize (const SparseMatrix   &matrix,
                                     const AdditionalData &additional_data)
  {
    preconditioner.reset ();
    preconditioner.reset (new Ifpack_Chebyshev (&matrix.trilinos_matrix()));

    Ifpack_Chebyshev *ifpack = static_cast<Ifpack_Chebyshev *>
                               (preconditioner.get());
    Assert (ifpack != 0, ExcMessage ("Trilinos could not create this "
                                     "preconditioner"));

    int ierr;

    Teuchos::ParameterList parameter_list;
    parameter_list.set ("chebyshev: ratio eigenvalue",
                        additional_data.eigenvalue_ratio);
    parameter_list.set ("chebyshev: min eigenvalue",
                        additional_data.min_eigenvalue);
    parameter_list.set ("chebyshev: max eigenvalue",
                        additional_data.max_eigenvalue);
    parameter_list.set ("chebyshev: degree",
                        (int)additional_data.degree);
    parameter_list.set ("chebyshev: min diagonal value",
                        additional_data.min_diagonal);
    parameter_list.set ("chebyshev: zero starting solution",
                        !additional_data.nonzero_starting);

    ierr = ifpack->SetParameters(parameter_list);
    AssertThrow (ierr == 0, ExcTrilinosError(ierr));

    ierr = ifpack->Initialize();
    AssertThrow (ierr == 0, ExcTrilinosError(ierr));

    ierr = ifpack->Compute();
    AssertThrow (ierr == 0, ExcTrilinosError(ierr));
  }



  /* -------------------------- PreconditionAMG -------------------------- */

  PreconditionAMG::AdditionalData::
  AdditionalData (const bool                             elliptic,
                  const bool                             higher_order_elements,
                  const unsigned int                     n_cycles,
                  const bool                             w_cycle,
                  const double                           aggregation_threshold,
                  const std::vector<std::vector<bool> > &constant_modes,
                  const unsigned int                     smoother_sweeps,
                  const unsigned int                     smoother_overlap,
                  const bool                             output_details,
                  const char                            *smoother_type,
                  const char                            *coarse_type)
    :
    elliptic (elliptic),
    higher_order_elements (higher_order_elements),
    n_cycles (n_cycles),
    w_cycle (w_cycle),
    aggregation_threshold (aggregation_threshold),
    constant_modes (constant_modes),
    smoother_sweeps (smoother_sweeps),
    smoother_overlap (smoother_overlap),
    output_details (output_details),
    smoother_type (smoother_type),
    coarse_type (coarse_type)
  {}


  PreconditionAMG::~PreconditionAMG()
  {
    preconditioner.reset();
    trilinos_matrix.reset();
  }



  void
  PreconditionAMG:: initialize (const SparseMatrix   &matrix,
                                const AdditionalData &additional_data)
  {
    initialize(matrix.trilinos_matrix(), additional_data);
  }



  void
  PreconditionAMG:: initialize (const Epetra_RowMatrix &matrix,
                                const AdditionalData   &additional_data)
  {
    // Build the AMG preconditioner.
    Teuchos::ParameterList parameter_list;

    if (additional_data.elliptic == true)
      {
        ML_Epetra::SetDefaults("SA",parameter_list);

        // uncoupled mode can give a lot of warnings or even fail when there
        // are too many entries per row and aggreggation gets complicated, but
        // MIS does not work if too few elements are located on one
        // processor. work around these warnings by choosing the different
        // strategies in different situations: for low order, always use the
        // standard choice uncoupled. if higher order, right now we also just
        // use Uncoupled, but we should be aware that maybe MIS might be
        // needed
        if (additional_data.higher_order_elements)
          parameter_list.set("aggregation: type", "Uncoupled");
      }
    else
      {
        ML_Epetra::SetDefaults("NSSA",parameter_list);
        parameter_list.set("aggregation: type", "Uncoupled");
        parameter_list.set("aggregation: block scaling", true);
      }

    parameter_list.set("smoother: type", additional_data.smoother_type);
    parameter_list.set("coarse: type", additional_data.coarse_type);

    parameter_list.set("smoother: sweeps",
                       static_cast<int>(additional_data.smoother_sweeps));
    parameter_list.set("cycle applications",
                       static_cast<int>(additional_data.n_cycles));
    if (additional_data.w_cycle == true)
      parameter_list.set("prec type", "MGW");
    else
      parameter_list.set("prec type", "MGV");

    parameter_list.set("smoother: Chebyshev alpha",10.);
    parameter_list.set("smoother: ifpack overlap",
                       static_cast<int>(additional_data.smoother_overlap));
    parameter_list.set("aggregation: threshold",
                       additional_data.aggregation_threshold);
    parameter_list.set("coarse: max size", 2000);

    if (additional_data.output_details)
      parameter_list.set("ML output", 10);
    else
      parameter_list.set("ML output", 0);

    const Epetra_Map &domain_map = matrix.OperatorDomainMap();

    const size_type constant_modes_dimension =
      additional_data.constant_modes.size();
    Epetra_MultiVector distributed_constant_modes (domain_map,
                                                   constant_modes_dimension > 0 ?
                                                   constant_modes_dimension : 1);
    std::vector<double> dummy (constant_modes_dimension);

    if (constant_modes_dimension > 0)
      {
        const size_type n_rows = n_global_rows(matrix);
        const bool constant_modes_are_global =
          additional_data.constant_modes[0].size() == n_rows;
        const size_type n_relevant_rows =
          constant_modes_are_global ? n_rows : additional_data.constant_modes[0].size();
        const size_type my_size = domain_map.NumMyElements();
        if (constant_modes_are_global == false)
          Assert (n_relevant_rows == my_size,
                  ExcDimensionMismatch(n_relevant_rows, my_size));
        Assert (n_rows ==
                static_cast<size_type>(global_length(distributed_constant_modes)),
                ExcDimensionMismatch(n_rows,
                                     global_length(distributed_constant_modes)));

        // Reshape null space as a contiguous vector of doubles so that
        // Trilinos can read from it.
        for (size_type d=0; d<constant_modes_dimension; ++d)
          for (size_type row=0; row<my_size; ++row)
            {
              TrilinosWrappers::types::int_type global_row_id =
                constant_modes_are_global ? gid(domain_map,row) : row;
              distributed_constant_modes[d][row] =
                additional_data.constant_modes[d][global_row_id];
            }

        parameter_list.set("null space: type", "pre-computed");
        parameter_list.set("null space: dimension",
                           distributed_constant_modes.NumVectors());
        if (my_size > 0)
          parameter_list.set("null space: vectors",
                             distributed_constant_modes.Values());
        // We need to set a valid pointer to data even if there is no data on
        // the current processor. Therefore, pass a dummy in that case
        else
          parameter_list.set("null space: vectors",
                             &dummy[0]);
      }

    initialize (matrix, parameter_list);

    if (additional_data.output_details)
      {
        ML_Epetra::MultiLevelPreconditioner *multilevel_operator =
          dynamic_cast<ML_Epetra::MultiLevelPreconditioner *> (preconditioner.get());
        Assert (multilevel_operator != 0,
                ExcMessage ("Preconditioner setup failed."));
        multilevel_operator->PrintUnused(0);
      }
  }



  void
  PreconditionAMG::initialize (const SparseMatrix           &matrix,
                               const Teuchos::ParameterList &ml_parameters)
  {
    initialize(matrix.trilinos_matrix(), ml_parameters);
  }



  void
  PreconditionAMG::initialize (const Epetra_RowMatrix       &matrix,
                               const Teuchos::ParameterList &ml_parameters)
  {
    preconditioner.reset ();
    preconditioner.reset (new ML_Epetra::MultiLevelPreconditioner
                          (matrix, ml_parameters));
  }



  template <typename number>
  void
  PreconditionAMG::
  initialize (const ::dealii::SparseMatrix<number> &deal_ii_sparse_matrix,
              const AdditionalData                 &additional_data,
              const double                          drop_tolerance,
              const ::dealii::SparsityPattern      *use_this_sparsity)
  {
    preconditioner.reset();
    const size_type n_rows = deal_ii_sparse_matrix.m();

    // Init Epetra Matrix using an
    // equidistributed map; avoid
    // storing the nonzero
    // elements.
    vector_distributor.reset (new Epetra_Map(static_cast<TrilinosWrappers::types::int_type>(n_rows),
                                             0, communicator));

    if (trilinos_matrix.get() == 0)
      trilinos_matrix.reset (new SparseMatrix());

    trilinos_matrix->reinit (*vector_distributor, *vector_distributor,
                             deal_ii_sparse_matrix, drop_tolerance, true,
                             use_this_sparsity);

    initialize (*trilinos_matrix, additional_data);
  }



  void PreconditionAMG::reinit ()
  {
    ML_Epetra::MultiLevelPreconditioner *multilevel_operator =
      dynamic_cast<ML_Epetra::MultiLevelPreconditioner *> (preconditioner.get());
    multilevel_operator->ReComputePreconditioner();
  }



  void PreconditionAMG::clear ()
  {
    PreconditionBase::clear();
    trilinos_matrix.reset();
  }



  PreconditionAMG::size_type
  PreconditionAMG::memory_consumption() const
  {
    unsigned int memory = sizeof(this);

    // todo: find a way to read out ML's data
    // sizes
    if (trilinos_matrix.get() != 0)
      memory += trilinos_matrix->memory_consumption();
    return memory;
  }




  // explicit instantiations
  template void PreconditionAMG::initialize (const ::dealii::SparseMatrix<double> &,
                                             const AdditionalData &, const double,
                                             const ::dealii::SparsityPattern *);
  template void PreconditionAMG::initialize (const ::dealii::SparseMatrix<float> &,
                                             const AdditionalData &, const double,
                                             const ::dealii::SparsityPattern *);




  /* -------------------------- PreconditionAMG -------------------------- */

  void
  PreconditionIdentity::vmult(VectorBase       &dst,
                              const VectorBase &src) const
  {
    dst = src;
  }

  void
  PreconditionIdentity::Tvmult(VectorBase       &dst,
                               const VectorBase &src) const
  {
    dst = src;
  }

  void
  PreconditionIdentity::vmult(dealii::Vector<double>       &dst,
                              const dealii::Vector<double> &src) const
  {
    dst = src;
  }

  void
  PreconditionIdentity::Tvmult(dealii::Vector<double>       &dst,
                               const dealii::Vector<double> &src) const
  {
    dst = src;
  }

  void
  PreconditionIdentity::vmult(parallel::distributed::Vector<double>       &dst,
                              const parallel::distributed::Vector<double> &src) const
  {
    dst = src;
  }

  void
  PreconditionIdentity::Tvmult(parallel::distributed::Vector<double>       &dst,
                               const parallel::distributed::Vector<double> &src) const
  {
    dst = src;
  }
}

DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_WITH_TRILINOS
