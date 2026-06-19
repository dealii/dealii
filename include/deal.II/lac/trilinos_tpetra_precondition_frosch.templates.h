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

#ifndef dealii_trilinos_tpetra_precondition_frosch_templates_h
#define dealii_trilinos_tpetra_precondition_frosch_templates_h

#include <deal.II/base/config.h>

#if defined(DEAL_II_TRILINOS_WITH_TPETRA) && \
  defined(DEAL_II_TRILINOS_WITH_SHYLU_DDFROSCH)

#  include <deal.II/base/geometry_info.h>
#  include <deal.II/base/index_set.h>

#  include <deal.II/distributed/shared_tria.h>
#  include <deal.II/distributed/tria.h>

#  include <deal.II/dofs/dof_handler.h>
#  include <deal.II/dofs/dof_tools.h>

#  include <deal.II/grid/tria.h>

#  include <deal.II/lac/trilinos_tpetra_precondition.h>
#  include <deal.II/lac/trilinos_tpetra_precondition.templates.h>
#  include <deal.II/lac/trilinos_tpetra_types.h>
#  include <deal.II/lac/trilinos_xpetra_types.h>

#  include <FROSch_OneLevelPreconditioner_def.hpp>
#  include <FROSch_SchwarzPreconditioners_fwd.hpp>
#  include <FROSch_Tools_def.hpp>
#  include <ShyLU_DDFROSch_config.h>
#  include <Teuchos_ParameterList.hpp>
#  include <Teuchos_RCP.hpp>

#  include <string>
#endif

DEAL_II_NAMESPACE_OPEN

#if defined(DEAL_II_TRILINOS_WITH_TPETRA) && \
  defined(DEAL_II_TRILINOS_WITH_SHYLU_DDFROSCH)

namespace LinearAlgebra
{
  namespace TpetraWrappers
  {
    /* ---------------------- PreconditionFROSch ------------------------ */
    // Function to convert a Xpetra::Operator into a Tpetra::Operator
    namespace internal
    {
      template <typename Number, typename MemorySpace>
      Teuchos::RCP<TpetraTypes::LinearOperator<Number, MemorySpace>>
      XpetraToTpetra(
        Teuchos::RCP<XpetraTypes::LinearOperator<Number, MemorySpace>>
          xpetra_prec)
      {
        Teuchos::RCP<
          XpetraTypes::XpetraToTpetraLinearOperator<Number, MemorySpace>>
          frosch_tpetra_prec = Teuchos::rcp(
            new XpetraTypes::XpetraToTpetraLinearOperator<Number, MemorySpace>(
              xpetra_prec));

        Teuchos::RCP<TpetraTypes::LinearOperator<Number, MemorySpace>>
          tpetra_prec = Teuchos::rcp_dynamic_cast<
            TpetraTypes::LinearOperator<Number, MemorySpace>>(
            frosch_tpetra_prec);

        return tpetra_prec;
      }


      /**
       * Create a Teuchos::ParameterList and fill it with the information from
       * the AdditionalData struct.
       */
      template <typename Number, typename MemorySpace>
      Teuchos::RCP<Teuchos::ParameterList>
      AdditionalData_to_ParameterList(
        const typename PreconditionFROSch<Number, MemorySpace>::AdditionalData
          &additional_data)
      {
        Teuchos::RCP<Teuchos::ParameterList> parameter_list =
          Teuchos::rcp(new Teuchos::ParameterList("Preconditioner List"));

        parameter_list->set("OverlappingOperator Type",
                            "AlgebraicOverlappingOperator");
        // Note: Teuchos::ParameterList expects an int (and not an unsigned int)
        // here, otherwise this will fail silently.
        parameter_list->set("Dimension", (int)additional_data.dimension);
        parameter_list->set("Overlap", (int)additional_data.overlap);

        std::string coarse_operator_string;
        switch (additional_data.coarse_operator_type)
          {
            case PreconditionFROSch<Number,
                                    MemorySpace>::IPOUHarmonicCoarseOperator:
              coarse_operator_string = "IPOUHarmonicCoarseOperator";
              break;
            case PreconditionFROSch<Number, MemorySpace>::GDSWCoarseOperator:
              coarse_operator_string = "GDSWCoarseOperator";
              break;
            case PreconditionFROSch<Number, MemorySpace>::RGDSWCoarseOperator:
              coarse_operator_string = "RGDSWCoarseOperator";
              break;
            default:
              DEAL_II_NOT_IMPLEMENTED();
          }
        parameter_list->set("CoarseOperator Type", coarse_operator_string);

        { /* sublist: AlgebraicOverlappingOperator */
          Teuchos::RCP<Teuchos::ParameterList> operator_parameter_list =
            Teuchos::rcp(
              new Teuchos::ParameterList("AlgebraicOverlappingOperator"));
          switch (additional_data.combine_values_in_overlap)
            {
              case PreconditionFROSch<Number, MemorySpace>::Restricted:
                operator_parameter_list->set("Combine Values in Overlap",
                                             "Restricted");
                break;
              case PreconditionFROSch<Number, MemorySpace>::Averaging:
                operator_parameter_list->set("Combine Values in Overlap",
                                             "Averaging");
                break;
              case PreconditionFROSch<Number, MemorySpace>::Full:
                operator_parameter_list->set("Combine Values in Overlap",
                                             "Full");
                break;
              default:
                DEAL_II_NOT_IMPLEMENTED();
            }
          operator_parameter_list->set("Adding Layers Strategy", "CrsMatrix");

          { /* sub-sublist: Solver */
            Teuchos::RCP<Teuchos::ParameterList> solver_parameter_list =
              Teuchos::rcp(new Teuchos::ParameterList("Solver"));
            solver_parameter_list->set("SolverType", "Amesos2");
            switch (additional_data.subdomain_solver)
              {
                case PreconditionFROSch<Number, MemorySpace>::KLU:
                  solver_parameter_list->set("Solver", "KLU");
                  break;
                case PreconditionFROSch<Number, MemorySpace>::UMFPACK:
                  solver_parameter_list->set("Solver", "UMFPACK");
                  break;
                case PreconditionFROSch<Number, MemorySpace>::MUMPS:
                  solver_parameter_list->set("Solver", "MUMPS");
                  break;
                case PreconditionFROSch<Number, MemorySpace>::SuperLU_dist:
                  solver_parameter_list->set("Solver", "SuperLU_dist");
                  break;
                default:
                  DEAL_II_NOT_IMPLEMENTED();
              }

            operator_parameter_list->set("Solver", *solver_parameter_list);
          }

          parameter_list->set("AlgebraicOverlappingOperator",
                              *operator_parameter_list);
        }

        { /* sublist: <coarse_operator_type> */
          Teuchos::RCP<Teuchos::ParameterList>
            coarse_operator_type_parameter_list =
              Teuchos::rcp(new Teuchos::ParameterList(coarse_operator_string));

          { /* sub-sublist: ExtensionSolver */
            Teuchos::RCP<Teuchos::ParameterList>
              extension_solver_parameter_list =
                Teuchos::rcp(new Teuchos::ParameterList("ExtensionSolver"));
            extension_solver_parameter_list->set("SolverType", "Amesos2");
            switch (additional_data.extension_solver)
              {
                case PreconditionFROSch<Number, MemorySpace>::KLU:
                  extension_solver_parameter_list->set("Solver", "KLU");
                  break;
                case PreconditionFROSch<Number, MemorySpace>::UMFPACK:
                  extension_solver_parameter_list->set("Solver", "UMFPACK");
                  break;
                case PreconditionFROSch<Number, MemorySpace>::MUMPS:
                  extension_solver_parameter_list->set("Solver", "MUMPS");
                  break;
                case PreconditionFROSch<Number, MemorySpace>::SuperLU_dist:
                  extension_solver_parameter_list->set("Solver",
                                                       "SuperLU_dist");
                  break;
                default:
                  DEAL_II_NOT_IMPLEMENTED();
              }

            coarse_operator_type_parameter_list->set(
              "ExtensionSolver", *extension_solver_parameter_list);
          }

          { /* sub-sublist: CoarseSolver */
            Teuchos::RCP<Teuchos::ParameterList> coarse_solver_parameter_list =
              Teuchos::rcp(new Teuchos::ParameterList("CoarseSolver"));
            coarse_solver_parameter_list->set("SolverType", "Amesos2");
            switch (additional_data.coarse_solver)
              {
                case PreconditionFROSch<Number, MemorySpace>::KLU:
                  coarse_solver_parameter_list->set("Solver", "KLU");
                  break;
                case PreconditionFROSch<Number, MemorySpace>::UMFPACK:
                  coarse_solver_parameter_list->set("Solver", "UMFPACK");
                  break;
                case PreconditionFROSch<Number, MemorySpace>::MUMPS:
                  coarse_solver_parameter_list->set("Solver", "MUMPS");
                  break;
                case PreconditionFROSch<Number, MemorySpace>::SuperLU_dist:
                  coarse_solver_parameter_list->set("Solver", "SuperLU_dist");
                  break;
                default:
                  DEAL_II_NOT_IMPLEMENTED();
              }

            coarse_operator_type_parameter_list->set(
              "CoarseSolver", *coarse_solver_parameter_list);
          }

          parameter_list->set(coarse_operator_string,
                              *coarse_operator_type_parameter_list);
        }

        return parameter_list;
      }

      // Create a Xpetra::CrsMatrix from a dealii::SparseMatrix
      template <typename Number, typename MemorySpace>
      Teuchos::RCP<XpetraTypes::MatrixType<Number, MemorySpace>>
      SparseMatrix_to_CRSMatrix(const SparseMatrix<Number, MemorySpace> &A)
      {
        // create a Xpetra::CrsMatrix Object,
        // which is used to create the Xpetra::Matrix
        Teuchos::RCP<XpetraTypes::CrsMatrixType<Number, MemorySpace>>
          xpetra_crsmatrix = Teuchos::rcp(
            new XpetraTypes::TpetraCrsMatrixType<Number, MemorySpace>(
              Teuchos::rcp_const_cast<
                TpetraTypes::MatrixType<Number, MemorySpace>>(
                A.trilinos_rcp())));

        // Create from the above defined Xpetra::CrsMatrix
        // an Xpetra::Matrix
        Teuchos::RCP<XpetraTypes::MatrixType<Number, MemorySpace>>
          xpetra_matrix = Teuchos::rcp(
            new XpetraTypes::CrsMatrixWrapType<Number, MemorySpace>(
              xpetra_crsmatrix));

        return xpetra_matrix;
      }

    } // namespace internal



    template <typename Number, typename MemorySpace>
    PreconditionFROSch<Number, MemorySpace>::PreconditionFROSch(
      const enum PreconditionerType precondition_type)
      : precondition_type(precondition_type)
      , user_provided_parameter_list(false)
    {}



    template <typename Number, typename MemorySpace>
    void
    PreconditionFROSch<Number, MemorySpace>::set_parameter_list(
      const Teuchos::ParameterList &parameter_list)
    {
      this->parameter_list.setParameters(parameter_list);
      user_provided_parameter_list = true;
    }



    template <typename Number, typename MemorySpace>
    PreconditionFROSch<Number, MemorySpace>::AdditionalData::AdditionalData(
      const unsigned int       overlap,
      const unsigned int       dimension,
      const enum CombineMethod combine_values_in_overlap,
      const enum SolverName    subdomain_solver,
      const enum CoarseType    coarse_operator_type,
      const enum SolverName    extension_solver,
      const enum SolverName    coarse_solver)
      : overlap(overlap)
      , dimension(dimension)
      , combine_values_in_overlap(combine_values_in_overlap)
      , subdomain_solver(subdomain_solver)
      , coarse_operator_type(coarse_operator_type)
      , extension_solver(extension_solver)
      , coarse_solver(coarse_solver)
    {}



    template <typename Number, typename MemorySpace>
    void
    PreconditionFROSch<Number, MemorySpace>::initialize(
      const SparseMatrix<Number, MemorySpace> &A,
      const AdditionalData                    &additional_data)
    {
      // If the user did not provide a parameter list, use the AdditionalData
      if (!user_provided_parameter_list)
        {
          const Teuchos::RCP<Teuchos::ParameterList> parameters =
            internal::AdditionalData_to_ParameterList<Number, MemorySpace>(
              additional_data);

          // store the parameter list
          this->parameter_list.setParameters(*parameters);
        }

      const Teuchos::RCP<XpetraTypes::MatrixType<Number, MemorySpace>>
        xpetra_matrix = internal::SparseMatrix_to_CRSMatrix(A);

      switch (precondition_type)
        {
          case OneLevel:
            {
              // The one-level Schwarz preconditioner object
              Teuchos::RCP<XpetraTypes::FROSchOneLevelType<Number, MemorySpace>>
                prec(new XpetraTypes::FROSchOneLevelType<Number, MemorySpace>(
                  xpetra_matrix,
                  Teuchos::rcpFromRef<Teuchos::ParameterList>(
                    this->parameter_list)));

              // Initialize
              prec->initialize(false);

              // THIS ACTUALLY COMPUTES THE PRECONDITIONER
              prec->compute();

              // convert the FROSch preconditioner into a Xpetra::Operator
              Teuchos::RCP<XpetraTypes::LinearOperator<Number, MemorySpace>>
                xpetra_prec = Teuchos::rcp_dynamic_cast<
                  XpetraTypes::LinearOperator<Number, MemorySpace>>(prec);

              // convert the FROSch preconditioner into a Tpetra::Operator
              // (The OneLevelOperator is derived from the Xpetra::Operator)
              this->preconditioner =
                internal::XpetraToTpetra<Number, MemorySpace>(xpetra_prec);

              break;
            }
          case TwoLevel:
            {
              // The one-level Schwarz preconditioner object
              Teuchos::RCP<XpetraTypes::FROSchTwoLevelType<Number, MemorySpace>>
                prec(new XpetraTypes::FROSchTwoLevelType<Number, MemorySpace>(
                  xpetra_matrix,
                  Teuchos::rcpFromRef<Teuchos::ParameterList>(
                    this->parameter_list)));

              // Initialize
              prec->initialize(false);

              // THIS ACTUALLY COMPUTES THE PRECONDITIONER
              prec->compute();

              // convert the FROSch preconditioner into a Xpetra::Operator
              Teuchos::RCP<XpetraTypes::LinearOperator<Number, MemorySpace>>
                xpetra_prec = Teuchos::rcp_dynamic_cast<
                  XpetraTypes::LinearOperator<Number, MemorySpace>>(prec);

              // convert the FROSch preconditioner into a Tpetra::Operator
              // (The OneLevelOperator is derived from the Xpetra::Operator)
              this->preconditioner =
                internal::XpetraToTpetra<Number, MemorySpace>(xpetra_prec);

              break;
            }
          default:
            DEAL_II_NOT_IMPLEMENTED();
        }
    }



    template <typename Number, typename MemorySpace>
    template <int dim, int spacedim>
    void
    PreconditionFROSch<Number, MemorySpace>::initialize(
      const SparseMatrix<Number, MemorySpace> &A,
      const DoFHandler<dim, spacedim>         &dof_handler,
      const AdditionalData                    &additional_data)
    {
      // If the user did not provide a parameter list, use the AdditionalData
      if (!user_provided_parameter_list)
        {
          Teuchos::RCP<Teuchos::ParameterList> parameters =
            internal::AdditionalData_to_ParameterList<Number, MemorySpace>(
              additional_data);

          // store the parameter list
          this->parameter_list.setParameters(*parameters);
        }

      // get the overlap:
      const int overlap = this->parameter_list.get("Overlap", 1);

      // get the mpi communicatior from the dof handler
      const MPI_Comm mpi_communicator = dof_handler.get_mpi_communicator();

      const Teuchos::RCP<XpetraTypes::MatrixType<Number, MemorySpace>>
        xpetra_matrix = internal::SparseMatrix_to_CRSMatrix(A);

      // Repeated Map
      const Teuchos::RCP<XpetraTypes::MapType<MemorySpace>> uniqueMap =
        Teuchos::rcp(new XpetraTypes::TpetraMapType<MemorySpace>(
          dof_handler.locally_owned_dofs().make_tpetra_map_rcp(mpi_communicator,
                                                               true)));

      // Extract the locally_active_dofs and turn the IndexSet into an
      // Teuchos::Array
      Teuchos::Array<types::signed_global_dof_index> cell_dofs;
      {
        IndexSet locally_active_dofs =
          DoFTools::extract_locally_active_dofs(dof_handler);
        const std::vector<types::global_dof_index> indices =
          locally_active_dofs.get_index_vector();
        std::vector<types::signed_global_dof_index> int_indices(indices.size());
        std::copy(indices.begin(), indices.end(), int_indices.begin());
        cell_dofs = int_indices;
      }

      const Teuchos::RCP<XpetraTypes::MapType<MemorySpace>> repeated_map =
        XpetraTypes::MapFactoryType<MemorySpace>::Build(
          Xpetra::UseTpetra,
          uniqueMap->getGlobalNumElements(),
          cell_dofs,
          0,
          uniqueMap->getComm());

      switch (precondition_type)
        {
          case OneLevel:
            {
              // The one-level Schwarz preconditioner object
              Teuchos::RCP<XpetraTypes::FROSchOneLevelType<Number, MemorySpace>>
                prec(new XpetraTypes::FROSchOneLevelType<Number, MemorySpace>(
                  xpetra_matrix,
                  Teuchos::rcpFromRef<Teuchos::ParameterList>(
                    this->parameter_list)));

              // Initialize
              prec->initialize(overlap, repeated_map);

              // THIS ACTUALLY COMPUTES THE PRECONDITIONER
              prec->compute();

              // convert the FROSch preconditioner into a Xpetra::Operator
              Teuchos::RCP<XpetraTypes::LinearOperator<Number, MemorySpace>>
                xpetra_prec = Teuchos::rcp_dynamic_cast<
                  XpetraTypes::LinearOperator<Number, MemorySpace>>(prec);

              // convert the FROSch preconditioner into a Tpetra::Operator
              // (The OneLevelOperator is derived from the Xpetra::Operator)
              this->preconditioner =
                internal::XpetraToTpetra<Number, MemorySpace>(xpetra_prec);

              break;
            }
          case TwoLevel:
            {
              // The two-level Schwarz preconditioner object
              Teuchos::RCP<XpetraTypes::FROSchTwoLevelType<Number, MemorySpace>>
                prec(new XpetraTypes::FROSchTwoLevelType<Number, MemorySpace>(
                  xpetra_matrix,
                  Teuchos::rcpFromRef<Teuchos::ParameterList>(
                    this->parameter_list)));

              // Initialize
              const unsigned int dimension =
                this->parameter_list.get("Dimension", 2);
              AssertThrow(
                dof_handler.get_fe().n_dofs_per_cell() ==
                  (GeometryInfo<dim>::vertices_per_cell *
                   dof_handler.get_fe().n_dofs_per_vertex()),
                ExcMessage(
                  "Currently, the deal.II - FROSch interface for two-level "
                  "preconditioners is only implemented for the case where "
                  "the DoFs are located only on vertices."));
              prec->initialize(dimension,
                               dof_handler.get_fe().n_dofs_per_vertex(),
                               overlap,
                               FROSch::null,
                               FROSch::null,
                               FROSch::NodeWise,
                               repeated_map);

              // THIS ACTUALLY COMPUTES THE PRECONDITIONER
              prec->compute();

              // convert the FROSch preconditioner into a Xpetra::Operator
              Teuchos::RCP<XpetraTypes::LinearOperator<Number, MemorySpace>>
                xpetra_prec = Teuchos::rcp_dynamic_cast<
                  XpetraTypes::LinearOperator<Number, MemorySpace>>(prec);

              // convert the FROSch preconditioner into a Tpetra::Operator
              // (The OneLevelOperator is derived from the Xpetra::Operator)
              this->preconditioner =
                internal::XpetraToTpetra<Number, MemorySpace>(xpetra_prec);

              break;
            }
          default:
            DEAL_II_NOT_IMPLEMENTED();
        }
    }

  } // namespace TpetraWrappers
} // namespace LinearAlgebra

#endif // DEAL_II_TRILINOS_WITH_TPETRA && DEAL_II_TRILINOS_WITH_SHYLU_DDFROSCH

DEAL_II_NAMESPACE_CLOSE

#endif // dealii_trilinos_tpetra_precondition_frosch_templates_h
