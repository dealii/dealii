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

#include <deal.II/base/geometry_info.h>
#include <deal.II/base/index_set.h>

#include <deal.II/lac/trilinos_tpetra_types.h>
#include <deal.II/lac/trilinos_xpetra_types.h>

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>

#include <string>

#ifdef DEAL_II_TRILINOS_WITH_TPETRA
#  include <deal.II/lac/trilinos_tpetra_precondition.h>
#  include <deal.II/lac/trilinos_tpetra_precondition.templates.h>

#  ifdef DEAL_II_TRILINOS_WITH_SHYLU_DDFROSCH
#    include <FROSch_OneLevelPreconditioner_def.hpp>
#    include <FROSch_SchwarzPreconditioners_fwd.hpp>
#    include <FROSch_Tools_def.hpp>
#    include <ShyLU_DDFROSch_config.h>
#  endif // DEAL_II_TRILINOS_WITH_SHYLU_DDFROSCH

DEAL_II_NAMESPACE_OPEN

namespace LinearAlgebra
{
  namespace TpetraWrappers
  {

#  ifdef DEAL_II_TRILINOS_WITH_SHYLU_DDFROSCH

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
       * @brief Create a Teuchos::ParameterList and fill it with the information from
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
        parameter_list->set("Dimension", additional_data.dimension);
        parameter_list->set("Overlap", additional_data.overlap);
        parameter_list->set("CoarseOperator Type",
                            additional_data.coarse_operator_type);

        { /* sublist: AlgebraicOverlappingOperator */
          Teuchos::RCP<Teuchos::ParameterList> operator_parameter_list =
            Teuchos::rcp(
              new Teuchos::ParameterList("AlgebraicOverlappingOperator"));

          operator_parameter_list->set(
            "Combine Values in Overlap",
            additional_data.combine_values_in_overlap);
          operator_parameter_list->set("Adding Layers Strategy", "CrsMatrix");

          { /* sub-sublist: Solver */
            Teuchos::RCP<Teuchos::ParameterList> solver_parameter_list =
              Teuchos::rcp(new Teuchos::ParameterList("Solver"));
            solver_parameter_list->set("SolverType", "Amesos2");
            solver_parameter_list->set("Solver",
                                       additional_data.subdomain_solver);

            operator_parameter_list->set("Solver", *solver_parameter_list);
          }

          parameter_list->set("AlgebraicOverlappingOperator",
                              *operator_parameter_list);
        }

        { /* sublist: <coarse_operator_type> */
          Teuchos::RCP<Teuchos::ParameterList>
            coarse_operator_type_parameter_list = Teuchos::rcp(
              new Teuchos::ParameterList(additional_data.coarse_operator_type));

          { /* sub-sublist: ExtensionSolver */
            Teuchos::RCP<Teuchos::ParameterList>
              extension_solver_parameter_list =
                Teuchos::rcp(new Teuchos::ParameterList("ExtensionSolver"));
            extension_solver_parameter_list->set("SolverType", "Amesos2");
            extension_solver_parameter_list->set(
              "Solver", additional_data.extension_solver);

            coarse_operator_type_parameter_list->set(
              "ExtensionSolver", *extension_solver_parameter_list);
          }

          { /* sub-sublist: CoarseSolver */
            Teuchos::RCP<Teuchos::ParameterList> coarse_solver_parameter_list =
              Teuchos::rcp(new Teuchos::ParameterList("CoarseSolver"));
            coarse_solver_parameter_list->set("SolverType", "Amesos2");
            coarse_solver_parameter_list->set("Solver",
                                              additional_data.coarse_solver);

            coarse_operator_type_parameter_list->set(
              "CoarseSolver", *coarse_solver_parameter_list);
          }

          parameter_list->set(additional_data.coarse_operator_type,
                              *coarse_operator_type_parameter_list);
        }

        return parameter_list;
      }
    } // namespace internal



    template <typename Number, typename MemorySpace>
    PreconditionFROSch<Number, MemorySpace>::PreconditionFROSch(
      std::string precondition_type)
      : precondition_type(precondition_type)
      , user_provided_parameter_list(false)
    {
      // Check if the precondition type does exist
        AssertThrow(((precondition_type == "One Level") ||
            (precondition_type == "Two Level"))), dealii::ExcNotImplemented());
    }



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
      const int          overlap,
      const unsigned int dimension,
      const std::string  combine_values_in_overlap,
      const std::string  subdomain_solver,
      const std::string  coarse_operator_type,
      const std::string  extension_solver,
      const std::string  coarse_solver)
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
          Teuchos::RCP<Teuchos::ParameterList> parameters =
            internal::AdditionalData_to_ParameterList<Number, MemorySpace>(
              additional_data);

          // store the parameter list
          this->parameter_list.setParameters(*parameters);
        }

      // create a Xpetra::CrsMatrix Object,
      // which is used to create the Xpetra::Matrix
      Teuchos::RCP<XpetraTypes::CrsMatrixType<Number, MemorySpace>>
        xpetra_crsmatrix = Teuchos::rcp(
          new XpetraTypes::TpetraCrsMatrixType<Number, MemorySpace>(
            Teuchos::rcp_const_cast<
              TpetraTypes::MatrixType<Number, MemorySpace>>(A.trilinos_rcp())));

      // Create from the above defined Xpetra::CrsMatrix
      // an Xpetra::Matrix
      Teuchos::RCP<XpetraTypes::MatrixType<Number, MemorySpace>> xpetra_matrix =
        Teuchos::rcp(new XpetraTypes::CrsMatrixWrapType<Number, MemorySpace>(
          xpetra_crsmatrix));

      if (precondition_type == "One Level")
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
        }
      else if (precondition_type == "Two Level")
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
        }
      else
        AssertThrow(false, dealii::ExcNotImplemented());
    }



    template <typename Number, typename MemorySpace>
    template <int dim, int spacedim>
    void
    PreconditionFROSch<Number, MemorySpace>::initialize(
      const SparseMatrix<Number, MemorySpace> &A,
      DoFHandler<dim, spacedim>               &dof_handler,
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
      MPI_Comm mpi_communicator = dof_handler.get_mpi_communicator();

      // create a Xpetra::CrsMatrix Object,
      // which is used to create the Xpetra::Matrix
      Teuchos::RCP<XpetraTypes::CrsMatrixType<Number, MemorySpace>>
        xpetra_crsmatrix = Teuchos::rcp(
          new XpetraTypes::TpetraCrsMatrixType<Number, MemorySpace>(
            Teuchos::rcp_const_cast<
              TpetraTypes::MatrixType<Number, MemorySpace>>(A.trilinos_rcp())));

      // Create from the above defined Xpetra::CrsMatrix
      // an Xpetra::Matrix
      Teuchos::RCP<XpetraTypes::MatrixType<Number, MemorySpace>> xpetra_matrix =
        Teuchos::rcp(new XpetraTypes::CrsMatrixWrapType<Number, MemorySpace>(
          xpetra_crsmatrix));

      // Repeated Map
      Teuchos::RCP<XpetraTypes::MapType<MemorySpace>> uniqueMap =
        Teuchos::rcp(new XpetraTypes::TpetraMapType<MemorySpace>(
          dof_handler.locally_owned_dofs().make_tpetra_map_rcp(mpi_communicator,
                                                               true)));

      // The non-distributed Triangulation class does not provide the function
      // n_locally_owned_active_cells(), so we try to cast the Triangulation
      // either on parallel::distributed::Triangulation or on
      // parallel::shared::Triangulation
      unsigned int n_locally_owned_cells;
      try
        {
          auto &tria = dynamic_cast<
            const parallel::distributed::Triangulation<dim, spacedim> &>(
            dof_handler.get_triangulation());
          n_locally_owned_cells = tria.n_locally_owned_active_cells();
        }
      catch (const std::bad_cast &)
        {
          auto &tria = dynamic_cast<
            const parallel::shared::Triangulation<dim, spacedim> &>(
            dof_handler.get_triangulation());
          n_locally_owned_cells = tria.n_locally_owned_active_cells();
        }

      const unsigned int n_dofs_per_cell =
        dof_handler.get_fe().n_dofs_per_cell();

      Teuchos::Array<types::signed_global_dof_index> cell_dofs(
        n_locally_owned_cells * n_dofs_per_cell);
      std::vector<types::global_dof_index> local_dof_indices(n_dofs_per_cell);

      unsigned int cell_counter = 0;
      for (const auto &cell : dof_handler.active_cell_iterators())
        {
          if (!cell->is_locally_owned())
            continue;

          cell->get_dof_indices(local_dof_indices);
          for (unsigned int i = 0; i < n_dofs_per_cell; ++i)
            cell_dofs[(cell_counter * n_dofs_per_cell) + i] =
              local_dof_indices[i];
          ++cell_counter;
        }

      FROSch::sortunique(cell_dofs);

      Teuchos::RCP<XpetraTypes::MapType<MemorySpace>> repeated_map =
        XpetraTypes::MapFactoryType<MemorySpace>::Build(
          Xpetra::UseTpetra,
          uniqueMap->getGlobalNumElements(),
          cell_dofs,
          0,
          uniqueMap->getComm());

      if (precondition_type == "One Level")
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
        }
      else if (precondition_type == "Two Level")
        {
          // The two-level Schwarz preconditioner object
          Teuchos::RCP<XpetraTypes::FROSchTwoLevelType<Number, MemorySpace>>
            prec(new XpetraTypes::FROSchTwoLevelType<Number, MemorySpace>(
              xpetra_matrix,
              Teuchos::rcpFromRef<Teuchos::ParameterList>(
                this->parameter_list)));

          // Initialize
          const unsigned int dimension = this->parameter_list.get("Dimension", 2);
          AssertThrow(
            n_dofs_per_cell == (GeometryInfo<dim>::vertices_per_cell *
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
        }
      else
        AssertThrow(false, dealii::ExcNotImplemented());
    }

#  endif // DEAL_II_TRILINOS_WITH_SHYLU_DDFROSCH

  } // namespace TpetraWrappers
} // namespace LinearAlgebra

DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_TRILINOS_WITH_TPETRA

#endif // dealii_trilinos_tpetra_precondition_frosch_templates_h
