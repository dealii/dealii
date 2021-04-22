// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2020 by the deal.II authors
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

#ifndef dealii_vector_tools_project_templates_h
#define dealii_vector_tools_project_templates_h


#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/la_parallel_block_vector.h>
#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/lac/la_vector.h>
#include <deal.II/lac/petsc_block_vector.h>
#include <deal.II/lac/petsc_vector.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_control.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/trilinos_parallel_block_vector.h>
#include <deal.II/lac/trilinos_vector.h>

#include <deal.II/matrix_free/operators.h>

#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/vector_tools_boundary.h>
#include <deal.II/numerics/vector_tools_project.h>
#include <deal.II/numerics/vector_tools_rhs.h>

DEAL_II_NAMESPACE_OPEN

namespace VectorTools
{
  namespace internal
  {
    /**
     * Interpolate zero boundary values. We don't need to worry about a
     * mapping here because the function we evaluate for the DoFs is zero in
     * the mapped locations as well as in the original, unmapped locations
     */
    template <int dim, int spacedim, typename number>
    void
    interpolate_zero_boundary_values(
      const DoFHandler<dim, spacedim> &          dof_handler,
      std::map<types::global_dof_index, number> &boundary_values)
    {
      // loop over all boundary faces
      // to get all dof indices of
      // dofs on the boundary. note
      // that in 3d there are cases
      // where a face is not at the
      // boundary, yet one of its
      // lines is, and we should
      // consider the degrees of
      // freedom on it as boundary
      // nodes. likewise, in 2d and
      // 3d there are cases where a
      // cell is only at the boundary
      // by one vertex. nevertheless,
      // since we do not support
      // boundaries with dimension
      // less or equal to dim-2, each
      // such boundary dof is also
      // found from some other face
      // that is actually wholly on
      // the boundary, not only by
      // one line or one vertex
      typename DoFHandler<dim, spacedim>::active_cell_iterator
        cell = dof_handler.begin_active(),
        endc = dof_handler.end();
      std::vector<types::global_dof_index> face_dof_indices;
      for (; cell != endc; ++cell)
        for (auto f : GeometryInfo<dim>::face_indices())
          if (cell->at_boundary(f))
            {
              face_dof_indices.resize(cell->get_fe().n_dofs_per_face(f));
              cell->face(f)->get_dof_indices(face_dof_indices,
                                             cell->active_fe_index());
              for (unsigned int i = 0; i < cell->get_fe().n_dofs_per_face(f);
                   ++i)
                // enter zero boundary values
                // for all boundary nodes
                //
                // we need not care about
                // vector valued elements here,
                // since we set all components
                boundary_values[face_dof_indices[i]] = 0.;
            }
    }

    /**
     * Compute the boundary values to be used in the project() functions.
     */
    template <int dim,
              int spacedim,
              template <int, int> class M_or_MC,
              template <int> class Q_or_QC,
              typename number>
    void
    project_compute_b_v(
      const M_or_MC<dim, spacedim> &             mapping,
      const DoFHandler<dim, spacedim> &          dof,
      const Function<spacedim, number> &         function,
      const bool                                 enforce_zero_boundary,
      const Q_or_QC<dim - 1> &                   q_boundary,
      const bool                                 project_to_boundary_first,
      std::map<types::global_dof_index, number> &boundary_values)
    {
      if (enforce_zero_boundary == true)
        // no need to project boundary
        // values, but enforce
        // homogeneous boundary values
        // anyway
        interpolate_zero_boundary_values(dof, boundary_values);

      else
        // no homogeneous boundary values
        if (project_to_boundary_first == true)
        // boundary projection required
        {
          // set up a list of boundary
          // functions for the
          // different boundary
          // parts. We want the
          // function to hold on
          // all parts of the boundary
          const std::vector<types::boundary_id> used_boundary_ids =
            dof.get_triangulation().get_boundary_ids();

          std::map<types::boundary_id, const Function<spacedim, number> *>
            boundary_functions;
          for (const auto used_boundary_id : used_boundary_ids)
            boundary_functions[used_boundary_id] = &function;
          project_boundary_values(
            mapping, dof, boundary_functions, q_boundary, boundary_values);
        }
    }


    /*
     * MatrixFree implementation of project() for an arbitrary number of
     * components of the FiniteElement.
     */
    template <int components, int dim, typename Number, int spacedim>
    void
    project_matrix_free(
      const Mapping<dim, spacedim> &   mapping,
      const DoFHandler<dim, spacedim> &dof,
      const AffineConstraints<Number> &constraints,
      const Quadrature<dim> &          quadrature,
      const Function<
        spacedim,
        typename LinearAlgebra::distributed::Vector<Number>::value_type>
        &                                         function,
      LinearAlgebra::distributed::Vector<Number> &work_result,
      const bool                                  enforce_zero_boundary,
      const Quadrature<dim - 1> &                 q_boundary,
      const bool                                  project_to_boundary_first)
    {
      Assert(project_to_boundary_first == false, ExcNotImplemented());
      Assert(enforce_zero_boundary == false, ExcNotImplemented());
      (void)enforce_zero_boundary;
      (void)project_to_boundary_first;
      (void)q_boundary;

      AssertDimension(dof.get_fe_collection().size(), 1);

      Assert(dof.get_fe(0).n_components() == function.n_components,
             ExcDimensionMismatch(dof.get_fe(0).n_components(),
                                  function.n_components));
      Assert(dof.get_fe(0).n_components() == components,
             ExcDimensionMismatch(components, dof.get_fe(0).n_components()));

      Quadrature<dim> quadrature_mf;

      if (dof.get_fe(0).reference_cell() ==
          ReferenceCells::get_hypercube<dim>())
        quadrature_mf = QGauss<dim>(dof.get_fe().degree + 2);
      else
        // TODO: since we have currently only implemented a handful quadrature
        // rules for non-hypercube objects, we do not construct a new
        // quadrature rule with degree + 2 here but  use the user-provided
        // quadrature rule (which is guaranteed to be tabulated).
        quadrature_mf = quadrature;

      // set up mass matrix and right hand side
      typename MatrixFree<dim, Number>::AdditionalData additional_data;
      additional_data.tasks_parallel_scheme =
        MatrixFree<dim, Number>::AdditionalData::partition_color;
      additional_data.mapping_update_flags =
        (update_values | update_JxW_values);
      std::shared_ptr<MatrixFree<dim, Number>> matrix_free(
        new MatrixFree<dim, Number>());
      matrix_free->reinit(
        mapping, dof, constraints, quadrature_mf, additional_data);
      using MatrixType = MatrixFreeOperators::MassOperator<
        dim,
        -1,
        0,
        components,
        LinearAlgebra::distributed::Vector<Number>>;
      MatrixType mass_matrix;
      mass_matrix.initialize(matrix_free);
      mass_matrix.compute_diagonal();

      LinearAlgebra::distributed::Vector<Number> rhs, inhomogeneities;
      matrix_free->initialize_dof_vector(work_result);
      matrix_free->initialize_dof_vector(rhs);
      matrix_free->initialize_dof_vector(inhomogeneities);
      constraints.distribute(inhomogeneities);
      inhomogeneities *= -1.;

      {
        create_right_hand_side(
          mapping, dof, quadrature, function, rhs, constraints);

        // account for inhomogeneous constraints
        inhomogeneities.update_ghost_values();
        FEEvaluation<dim, -1, 0, components, Number> phi(*matrix_free);
        for (unsigned int cell = 0; cell < matrix_free->n_cell_batches();
             ++cell)
          {
            phi.reinit(cell);
            phi.read_dof_values_plain(inhomogeneities);
            phi.evaluate(EvaluationFlags::values);
            for (unsigned int q = 0; q < phi.n_q_points; ++q)
              phi.submit_value(phi.get_value(q), q);

            phi.integrate(EvaluationFlags::values);
            phi.distribute_local_to_global(rhs);
          }
        rhs.compress(VectorOperation::add);
      }

      // now invert the matrix
      // Allow for a maximum of 6*n steps to reduce the residual by 10^-12. n
      // steps may not be sufficient, since roundoff errors may accumulate for
      // badly conditioned matrices. This behavior can be observed, e.g. for
      // FE_Q_Hierarchical for degree higher than three.
      ReductionControl control(6 * rhs.size(), 0., 1e-12, false, false);
      SolverCG<LinearAlgebra::distributed::Vector<Number>> cg(control);
      PreconditionJacobi<MatrixType>                       preconditioner;
      preconditioner.initialize(mass_matrix, 1.);
      cg.solve(mass_matrix, work_result, rhs, preconditioner);
      work_result += inhomogeneities;

      constraints.distribute(work_result);
    }



    // Helper interface for the matrix-free implementation of project().
    // Used to determine the number of components.
    template <int dim, typename Number, int spacedim>
    void
    project_matrix_free_component(
      const Mapping<dim, spacedim> &   mapping,
      const DoFHandler<dim, spacedim> &dof,
      const AffineConstraints<Number> &constraints,
      const Quadrature<dim> &          quadrature,
      const Function<
        spacedim,
        typename LinearAlgebra::distributed::Vector<Number>::value_type>
        &                                         function,
      LinearAlgebra::distributed::Vector<Number> &work_result,
      const bool                                  enforce_zero_boundary,
      const Quadrature<dim - 1> &                 q_boundary,
      const bool                                  project_to_boundary_first)
    {
      switch (dof.get_fe(0).n_components())
        {
          case 1:
            project_matrix_free<1>(mapping,
                                   dof,
                                   constraints,
                                   quadrature,
                                   function,
                                   work_result,
                                   enforce_zero_boundary,
                                   q_boundary,
                                   project_to_boundary_first);
            break;

          case 2:
            project_matrix_free<2>(mapping,
                                   dof,
                                   constraints,
                                   quadrature,
                                   function,
                                   work_result,
                                   enforce_zero_boundary,
                                   q_boundary,
                                   project_to_boundary_first);
            break;

          case 3:
            project_matrix_free<3>(mapping,
                                   dof,
                                   constraints,
                                   quadrature,
                                   function,
                                   work_result,
                                   enforce_zero_boundary,
                                   q_boundary,
                                   project_to_boundary_first);
            break;

          case 4:
            project_matrix_free<4>(mapping,
                                   dof,
                                   constraints,
                                   quadrature,
                                   function,
                                   work_result,
                                   enforce_zero_boundary,
                                   q_boundary,
                                   project_to_boundary_first);
            break;

          default:
            Assert(false, ExcInternalError());
        }
    }



    /**
     * Helper interface for the matrix-free implementation of project(): avoid
     * instantiating the other helper functions for more than one VectorType
     * by copying from a LinearAlgebra::distributed::Vector.
     */
    template <int dim, typename VectorType, int spacedim>
    void
    project_matrix_free_copy_vector(
      const Mapping<dim, spacedim> &                             mapping,
      const DoFHandler<dim, spacedim> &                          dof,
      const AffineConstraints<typename VectorType::value_type> & constraints,
      const Quadrature<dim> &                                    quadrature,
      const Function<spacedim, typename VectorType::value_type> &function,
      VectorType &                                               vec_result,
      const bool                 enforce_zero_boundary,
      const Quadrature<dim - 1> &q_boundary,
      const bool                 project_to_boundary_first)
    {
      Assert(vec_result.size() == dof.n_dofs(),
             ExcDimensionMismatch(vec_result.size(), dof.n_dofs()));

      LinearAlgebra::distributed::Vector<typename VectorType::value_type>
        work_result;
      project_matrix_free_component(mapping,
                                    dof,
                                    constraints,
                                    quadrature,
                                    function,
                                    work_result,
                                    enforce_zero_boundary,
                                    q_boundary,
                                    project_to_boundary_first);

      const IndexSet &          locally_owned_dofs = dof.locally_owned_dofs();
      IndexSet::ElementIterator it                 = locally_owned_dofs.begin();
      for (; it != locally_owned_dofs.end(); ++it)
        ::dealii::internal::ElementAccess<VectorType>::set(work_result(*it),
                                                           *it,
                                                           vec_result);
      vec_result.compress(VectorOperation::insert);
    }

    /**
     * Return whether the boundary values try to constrain a degree of freedom
     * that is already constrained to something else
     */
    template <typename number>
    bool
    constraints_and_b_v_are_compatible(
      const AffineConstraints<number> &          constraints,
      std::map<types::global_dof_index, number> &boundary_values)
    {
      for (const auto &boundary_value : boundary_values)
        if (constraints.is_constrained(boundary_value.first))
          // TODO: This looks wrong -- shouldn't it be ==0 in the first
          // condition and && ?
          if (!(constraints.get_constraint_entries(boundary_value.first)
                    ->size() > 0 ||
                (constraints.get_inhomogeneity(boundary_value.first) ==
                 boundary_value.second)))
            return false;

      return true;
    }



    /**
     * Generic implementation of the project() function
     */
    template <int dim,
              int spacedim,
              typename VectorType,
              template <int, int> class M_or_MC,
              template <int> class Q_or_QC>
    void
    do_project(
      const M_or_MC<dim, spacedim> &                             mapping,
      const DoFHandler<dim, spacedim> &                          dof,
      const AffineConstraints<typename VectorType::value_type> & constraints,
      const Q_or_QC<dim> &                                       quadrature,
      const Function<spacedim, typename VectorType::value_type> &function,
      VectorType &                                               vec_result,
      const bool              enforce_zero_boundary,
      const Q_or_QC<dim - 1> &q_boundary,
      const bool              project_to_boundary_first)
    {
      using number = typename VectorType::value_type;
      Assert(dof.get_fe(0).n_components() == function.n_components,
             ExcDimensionMismatch(dof.get_fe(0).n_components(),
                                  function.n_components));
      Assert(vec_result.size() == dof.n_dofs(),
             ExcDimensionMismatch(vec_result.size(), dof.n_dofs()));

      // make up boundary values
      std::map<types::global_dof_index, number> boundary_values;
      project_compute_b_v(mapping,
                          dof,
                          function,
                          enforce_zero_boundary,
                          q_boundary,
                          project_to_boundary_first,
                          boundary_values);

      // check if constraints are compatible (see below)
      const bool constraints_are_compatible =
        constraints_and_b_v_are_compatible<number>(constraints,
                                                   boundary_values);

      // set up mass matrix and right hand side
      Vector<number>  vec(dof.n_dofs());
      SparsityPattern sparsity;
      {
        DynamicSparsityPattern dsp(dof.n_dofs(), dof.n_dofs());
        DoFTools::make_sparsity_pattern(dof,
                                        dsp,
                                        constraints,
                                        !constraints_are_compatible);

        sparsity.copy_from(dsp);
      }
      SparseMatrix<number> mass_matrix(sparsity);
      Vector<number>       tmp(mass_matrix.n());

      // If the constraints object does not conflict with the given boundary
      // values (i.e., it either does not contain boundary values or it contains
      // the same as boundary_values), we can let it call
      // distribute_local_to_global straight away, otherwise we need to first
      // interpolate the boundary values and then condense the matrix and vector
      if (constraints_are_compatible)
        {
          const Function<spacedim, number> *dummy = nullptr;
          MatrixCreator::create_mass_matrix(mapping,
                                            dof,
                                            quadrature,
                                            mass_matrix,
                                            function,
                                            tmp,
                                            dummy,
                                            constraints);
          if (boundary_values.size() > 0)
            MatrixTools::apply_boundary_values(
              boundary_values, mass_matrix, vec, tmp, true);
        }
      else
        {
          // create mass matrix and rhs at once, which is faster.
          MatrixCreator::create_mass_matrix(
            mapping, dof, quadrature, mass_matrix, function, tmp);
          MatrixTools::apply_boundary_values(
            boundary_values, mass_matrix, vec, tmp, true);
          constraints.condense(mass_matrix, tmp);
        }

      // Allow for a maximum of 5*n steps to reduce the residual by 10^-12. n
      // steps may not be sufficient, since roundoff errors may accumulate for
      // badly conditioned matrices
      ReductionControl control(5 * tmp.size(), 0., 1e-12, false, false);
      GrowingVectorMemory<Vector<number>> memory;
      SolverCG<Vector<number>>            cg(control, memory);

      PreconditionSSOR<SparseMatrix<number>> prec;
      prec.initialize(mass_matrix, 1.2);

      cg.solve(mass_matrix, vec, tmp, prec);
      constraints.distribute(vec);

      // copy vec into vec_result. we can't use vec_result itself above, since
      // it may be of another type than Vector<double> and that wouldn't
      // necessarily go together with the matrix and other functions
      for (unsigned int i = 0; i < vec.size(); ++i)
        ::dealii::internal::ElementAccess<VectorType>::set(vec(i),
                                                           i,
                                                           vec_result);
    }

    template <int dim, typename VectorType, int spacedim>
    void
    project_parallel(
      const Mapping<dim, spacedim> &                            mapping,
      const DoFHandler<dim, spacedim> &                         dof,
      const AffineConstraints<typename VectorType::value_type> &constraints,
      const Quadrature<dim> &                                   quadrature,
      const std::function<typename VectorType::value_type(
        const typename DoFHandler<dim, spacedim>::active_cell_iterator &,
        const unsigned int)> &                                  func,
      VectorType &                                              vec_result)
    {
      using Number = typename VectorType::value_type;
      Assert(dof.get_fe(0).n_components() == 1,
             ExcDimensionMismatch(dof.get_fe(0).n_components(), 1));
      Assert(vec_result.size() == dof.n_dofs(),
             ExcDimensionMismatch(vec_result.size(), dof.n_dofs()));

      // set up mass matrix and right hand side
      typename MatrixFree<dim, Number>::AdditionalData additional_data;
      additional_data.tasks_parallel_scheme =
        MatrixFree<dim, Number>::AdditionalData::partition_color;
      additional_data.mapping_update_flags =
        (update_values | update_JxW_values);
      std::shared_ptr<MatrixFree<dim, Number>> matrix_free(
        new MatrixFree<dim, Number>());
      matrix_free->reinit(mapping,
                          dof,
                          constraints,
                          QGauss<1>(dof.get_fe().degree + 2),
                          additional_data);
      using MatrixType = MatrixFreeOperators::
        MassOperator<dim, -1, 0, 1, LinearAlgebra::distributed::Vector<Number>>;
      MatrixType mass_matrix;
      mass_matrix.initialize(matrix_free);
      mass_matrix.compute_diagonal();

      using LocalVectorType = LinearAlgebra::distributed::Vector<Number>;
      LocalVectorType vec, rhs, inhomogeneities;
      matrix_free->initialize_dof_vector(vec);
      matrix_free->initialize_dof_vector(rhs);
      matrix_free->initialize_dof_vector(inhomogeneities);
      constraints.distribute(inhomogeneities);
      inhomogeneities *= -1.;

      // assemble right hand side:
      {
        FEValues<dim> fe_values(mapping,
                                dof.get_fe(),
                                quadrature,
                                update_values | update_JxW_values);

        const unsigned int dofs_per_cell = dof.get_fe().n_dofs_per_cell();
        const unsigned int n_q_points    = quadrature.size();
        Vector<Number>     cell_rhs(dofs_per_cell);
        std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
        typename DoFHandler<dim, spacedim>::active_cell_iterator
          cell = dof.begin_active(),
          endc = dof.end();
        for (; cell != endc; ++cell)
          if (cell->is_locally_owned())
            {
              cell_rhs = 0;
              fe_values.reinit(cell);
              for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
                {
                  const double val_q = func(cell, q_point);
                  for (unsigned int i = 0; i < dofs_per_cell; ++i)
                    cell_rhs(i) += (fe_values.shape_value(i, q_point) * val_q *
                                    fe_values.JxW(q_point));
                }

              cell->get_dof_indices(local_dof_indices);
              constraints.distribute_local_to_global(cell_rhs,
                                                     local_dof_indices,
                                                     rhs);
            }
        rhs.compress(VectorOperation::add);
      }

      mass_matrix.vmult_add(rhs, inhomogeneities);

      // now invert the matrix
      // Allow for a maximum of 5*n steps to reduce the residual by 10^-12. n
      // steps may not be sufficient, since roundoff errors may accumulate for
      // badly conditioned matrices. This behavior can be observed, e.g. for
      // FE_Q_Hierarchical for degree higher than three.
      ReductionControl control(5 * rhs.size(), 0., 1e-12, false, false);
      SolverCG<LinearAlgebra::distributed::Vector<Number>>    cg(control);
      typename PreconditionJacobi<MatrixType>::AdditionalData data(0.8);
      PreconditionJacobi<MatrixType>                          preconditioner;
      preconditioner.initialize(mass_matrix, data);
      cg.solve(mass_matrix, vec, rhs, preconditioner);
      vec += inhomogeneities;

      constraints.distribute(vec);

      const IndexSet &          locally_owned_dofs = dof.locally_owned_dofs();
      IndexSet::ElementIterator it                 = locally_owned_dofs.begin();
      for (; it != locally_owned_dofs.end(); ++it)
        ::dealii::internal::ElementAccess<VectorType>::set(vec(*it),
                                                           *it,
                                                           vec_result);
      vec_result.compress(VectorOperation::insert);
    }



    template <int dim, typename VectorType, int spacedim>
    void
    project_parallel(
      std::shared_ptr<const MatrixFree<dim, typename VectorType::value_type>>
                                                                matrix_free,
      const AffineConstraints<typename VectorType::value_type> &constraints,
      const std::function<VectorizedArray<typename VectorType::value_type>(
        const unsigned int,
        const unsigned int)> &                                  func,
      VectorType &                                              vec_result,
      const unsigned int                                        fe_component)
    {
      const DoFHandler<dim, spacedim> &dof =
        matrix_free->get_dof_handler(fe_component);

      using Number = typename VectorType::value_type;
      Assert(dof.get_fe(0).n_components() == 1,
             ExcDimensionMismatch(dof.get_fe(0).n_components(), 1));
      Assert(vec_result.size() == dof.n_dofs(),
             ExcDimensionMismatch(vec_result.size(), dof.n_dofs()));

      using MatrixType = MatrixFreeOperators::
        MassOperator<dim, -1, 0, 1, LinearAlgebra::distributed::Vector<Number>>;
      MatrixType mass_matrix;
      mass_matrix.initialize(matrix_free, {fe_component});
      mass_matrix.compute_diagonal();

      using LocalVectorType = LinearAlgebra::distributed::Vector<Number>;
      LocalVectorType vec, rhs, inhomogeneities;
      matrix_free->initialize_dof_vector(vec, fe_component);
      matrix_free->initialize_dof_vector(rhs, fe_component);
      matrix_free->initialize_dof_vector(inhomogeneities, fe_component);
      constraints.distribute(inhomogeneities);
      inhomogeneities *= -1.;

      // assemble right hand side:
      {
        FEEvaluation<dim, -1, 0, 1, Number> fe_eval(*matrix_free, fe_component);
        const unsigned int n_cells    = matrix_free->n_cell_batches();
        const unsigned int n_q_points = fe_eval.n_q_points;

        for (unsigned int cell = 0; cell < n_cells; ++cell)
          {
            fe_eval.reinit(cell);
            for (unsigned int q = 0; q < n_q_points; ++q)
              fe_eval.submit_value(func(cell, q), q);

            fe_eval.integrate(EvaluationFlags::values);
            fe_eval.distribute_local_to_global(rhs);
          }
        rhs.compress(VectorOperation::add);
      }

      mass_matrix.vmult_add(rhs, inhomogeneities);

      // now invert the matrix
      // Allow for a maximum of 5*n steps to reduce the residual by 10^-12. n
      // steps may not be sufficient, since roundoff errors may accumulate for
      // badly conditioned matrices. This behavior can be observed, e.g. for
      // FE_Q_Hierarchical for degree higher than three.
      ReductionControl control(5 * rhs.size(), 0., 1e-12, false, false);
      SolverCG<LinearAlgebra::distributed::Vector<Number>>    cg(control);
      typename PreconditionJacobi<MatrixType>::AdditionalData data(0.8);
      PreconditionJacobi<MatrixType>                          preconditioner;
      preconditioner.initialize(mass_matrix, data);
      cg.solve(mass_matrix, vec, rhs, preconditioner);
      vec += inhomogeneities;

      constraints.distribute(vec);

      const IndexSet &          locally_owned_dofs = dof.locally_owned_dofs();
      IndexSet::ElementIterator it                 = locally_owned_dofs.begin();
      for (; it != locally_owned_dofs.end(); ++it)
        ::dealii::internal::ElementAccess<VectorType>::set(vec(*it),
                                                           *it,
                                                           vec_result);
      vec_result.compress(VectorOperation::insert);
    }

    /**
     * Specialization of project() for the case dim==spacedim.
     * Check if we can use the MatrixFree implementation or need
     * to use the matrix based one.
     */
    template <typename VectorType, int dim>
    void
    project(
      const Mapping<dim> &                                      mapping,
      const DoFHandler<dim> &                                   dof,
      const AffineConstraints<typename VectorType::value_type> &constraints,
      const Quadrature<dim> &                                   quadrature,
      const Function<dim, typename VectorType::value_type> &    function,
      VectorType &                                              vec_result,
      const bool                 enforce_zero_boundary,
      const Quadrature<dim - 1> &q_boundary,
      const bool                 project_to_boundary_first)
    {
      // If we can, use the matrix-free implementation
      bool use_matrix_free =
        MatrixFree<dim, typename VectorType::value_type>::is_supported(
          dof.get_fe()) &&
        dof.get_fe().n_base_elements() == 1;

      // enforce_zero_boundary and project_to_boundary_first
      // are not yet supported.
      // We have explicit instantiations only if
      // the number of components is not too high.
      if (enforce_zero_boundary || project_to_boundary_first ||
          dof.get_fe(0).n_components() > 4)
        use_matrix_free = false;

      if (use_matrix_free)
        project_matrix_free_copy_vector(mapping,
                                        dof,
                                        constraints,
                                        quadrature,
                                        function,
                                        vec_result,
                                        enforce_zero_boundary,
                                        q_boundary,
                                        project_to_boundary_first);
      else
        {
          Assert((dynamic_cast<const parallel::TriangulationBase<dim> *>(
                    &(dof.get_triangulation())) == nullptr),
                 ExcNotImplemented());
          do_project(mapping,
                     dof,
                     constraints,
                     quadrature,
                     function,
                     vec_result,
                     enforce_zero_boundary,
                     q_boundary,
                     project_to_boundary_first);
        }
    }
  } // namespace internal

  template <int dim, typename VectorType, int spacedim>
  void
  project(const Mapping<dim, spacedim> &                            mapping,
          const DoFHandler<dim, spacedim> &                         dof,
          const AffineConstraints<typename VectorType::value_type> &constraints,
          const Quadrature<dim> &                                   quadrature,
          const std::function<typename VectorType::value_type(
            const typename DoFHandler<dim, spacedim>::active_cell_iterator &,
            const unsigned int)> &                                  func,
          VectorType &                                              vec_result)
  {
    internal::project_parallel<dim, VectorType, spacedim>(
      mapping, dof, constraints, quadrature, func, vec_result);
  }



  template <int dim, typename VectorType>
  void
  project(std::shared_ptr<const MatrixFree<
            dim,
            typename VectorType::value_type,
            VectorizedArray<typename VectorType::value_type>>>      matrix_free,
          const AffineConstraints<typename VectorType::value_type> &constraints,
          const unsigned int,
          const std::function<VectorizedArray<typename VectorType::value_type>(
            const unsigned int,
            const unsigned int)> &func,
          VectorType &            vec_result,
          const unsigned int      fe_component)
  {
    internal::project_parallel<dim, VectorType, dim>(
      matrix_free, constraints, func, vec_result, fe_component);
  }



  template <int dim, typename VectorType>
  void
  project(std::shared_ptr<const MatrixFree<
            dim,
            typename VectorType::value_type,
            VectorizedArray<typename VectorType::value_type>>>      matrix_free,
          const AffineConstraints<typename VectorType::value_type> &constraints,
          const std::function<VectorizedArray<typename VectorType::value_type>(
            const unsigned int,
            const unsigned int)> &                                  func,
          VectorType &                                              vec_result,
          const unsigned int fe_component)
  {
    project(matrix_free,
            constraints,
            matrix_free->get_dof_handler(fe_component).get_fe().degree + 1,
            func,
            vec_result,
            fe_component);
  }



  template <int dim, typename VectorType, int spacedim>
  void
  project(const Mapping<dim, spacedim> &                            mapping,
          const DoFHandler<dim, spacedim> &                         dof,
          const AffineConstraints<typename VectorType::value_type> &constraints,
          const Quadrature<dim> &                                   quadrature,
          const Function<spacedim, typename VectorType::value_type> &function,
          VectorType &                                               vec_result,
          const bool                 enforce_zero_boundary,
          const Quadrature<dim - 1> &q_boundary,
          const bool                 project_to_boundary_first)
  {
    if (dim == spacedim)
      {
        const Mapping<dim> *const mapping_ptr =
          dynamic_cast<const Mapping<dim> *>(&mapping);
        const DoFHandler<dim> *const dof_ptr =
          dynamic_cast<const DoFHandler<dim> *>(&dof);
        const Function<dim,
                       typename VectorType::value_type> *const function_ptr =
          dynamic_cast<const Function<dim, typename VectorType::value_type> *>(
            &function);
        Assert(mapping_ptr != nullptr, ExcInternalError());
        Assert(dof_ptr != nullptr, ExcInternalError());
        internal::project<VectorType, dim>(*mapping_ptr,
                                           *dof_ptr,
                                           constraints,
                                           quadrature,
                                           *function_ptr,
                                           vec_result,
                                           enforce_zero_boundary,
                                           q_boundary,
                                           project_to_boundary_first);
      }
    else
      {
        Assert(
          (dynamic_cast<const parallel::TriangulationBase<dim, spacedim> *>(
             &(dof.get_triangulation())) == nullptr),
          ExcNotImplemented());
        internal::do_project(mapping,
                             dof,
                             constraints,
                             quadrature,
                             function,
                             vec_result,
                             enforce_zero_boundary,
                             q_boundary,
                             project_to_boundary_first);
      }
  }



  template <int dim, typename VectorType, int spacedim>
  void
  project(const DoFHandler<dim, spacedim> &                         dof,
          const AffineConstraints<typename VectorType::value_type> &constraints,
          const Quadrature<dim> &                                   quadrature,
          const Function<spacedim, typename VectorType::value_type> &function,
          VectorType &                                               vec,
          const bool                 enforce_zero_boundary,
          const Quadrature<dim - 1> &q_boundary,
          const bool                 project_to_boundary_first)
  {
#ifdef _MSC_VER
    Assert(false,
           ExcMessage("Please specify the mapping explicitly "
                      "when building with MSVC!"));
#else
    project(get_default_linear_mapping(dof.get_triangulation()),
            dof,
            constraints,
            quadrature,
            function,
            vec,
            enforce_zero_boundary,
            q_boundary,
            project_to_boundary_first);
#endif
  }



  template <int dim, typename VectorType, int spacedim>
  void
  project(const hp::MappingCollection<dim, spacedim> &              mapping,
          const DoFHandler<dim, spacedim> &                         dof,
          const AffineConstraints<typename VectorType::value_type> &constraints,
          const hp::QCollection<dim> &                              quadrature,
          const Function<spacedim, typename VectorType::value_type> &function,
          VectorType &                                               vec_result,
          const bool                      enforce_zero_boundary,
          const hp::QCollection<dim - 1> &q_boundary,
          const bool                      project_to_boundary_first)
  {
    Assert((dynamic_cast<const parallel::TriangulationBase<dim, spacedim> *>(
              &(dof.get_triangulation())) == nullptr),
           ExcNotImplemented());

    internal::do_project(mapping,
                         dof,
                         constraints,
                         quadrature,
                         function,
                         vec_result,
                         enforce_zero_boundary,
                         q_boundary,
                         project_to_boundary_first);
  }


  template <int dim, typename VectorType, int spacedim>
  void
  project(const DoFHandler<dim, spacedim> &                         dof,
          const AffineConstraints<typename VectorType::value_type> &constraints,
          const hp::QCollection<dim> &                              quadrature,
          const Function<spacedim, typename VectorType::value_type> &function,
          VectorType &                                               vec,
          const bool                      enforce_zero_boundary,
          const hp::QCollection<dim - 1> &q_boundary,
          const bool                      project_to_boundary_first)
  {
    project(hp::StaticMappingQ1<dim, spacedim>::mapping_collection,
            dof,
            constraints,
            quadrature,
            function,
            vec,
            enforce_zero_boundary,
            q_boundary,
            project_to_boundary_first);
  }
} // namespace VectorTools

DEAL_II_NAMESPACE_CLOSE

#endif // dealii_vector_tools_project_templates_h
