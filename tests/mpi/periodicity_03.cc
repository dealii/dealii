/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2008 - 2019 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * The deal.II library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE.md at
 * the top level directory of deal.II.
 *
 * ---------------------------------------------------------------------
 */

//
// check solution for periodicity. The used test case is based on step-22.
// We consider a 3D cube and require periodicity with respect to
// the left and bottom boundary.
// We refine two times adaptively using the Kelly error estimator.
//

#include <deal.II/base/conditional_ostream.h>

#include <deal.II/distributed/grid_refinement.h>

#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/manifold_lib.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/block_sparsity_pattern.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/trilinos_block_sparse_matrix.h>
#include <deal.II/lac/trilinos_parallel_block_vector.h>
#include <deal.II/lac/trilinos_precondition.h>
#include <deal.II/lac/trilinos_solver.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"

namespace Step22
{
  template <int dim>
  class StokesProblem
  {
  public:
    StokesProblem(const unsigned int degree);
    void
    run();

  private:
    void
    setup_dofs();
    void
    assemble_system();
    void
    solve();
    void
    get_point_value(const Point<dim> point,
                    const int        proc,
                    Vector<double> & value) const;
    void
    check_periodicity() const;
    void
    output_results(const unsigned int refinement_cycle) const;
    void
    refine_mesh();

    const unsigned int degree;

    MPI_Comm mpi_communicator;

    parallel::distributed::Triangulation<dim> triangulation;
    FESystem<dim>                             fe;
    DoFHandler<dim>                           dof_handler;

    AffineConstraints<double> constraints;
    std::vector<IndexSet>     owned_partitioning;
    std::vector<IndexSet>     relevant_partitioning;

    TrilinosWrappers::BlockSparseMatrix system_matrix;

    TrilinosWrappers::MPI::BlockVector solution;
    TrilinosWrappers::MPI::BlockVector system_rhs;

    ConditionalOStream pcout;

    FullMatrix<double> rot_matrix;
    Tensor<1, dim>     offset;
  };


  template <class Matrix, class Preconditioner>
  class InverseMatrix : public Preconditioner
  {
  public:
    InverseMatrix(const Matrix &        m,
                  const Preconditioner &preconditioner,
                  const IndexSet &      locally_owned,
                  const MPI_Comm &      mpi_communicator);

    void
    vmult(TrilinosWrappers::MPI::Vector &      dst,
          const TrilinosWrappers::MPI::Vector &src) const;

  private:
    const SmartPointer<const Matrix>         matrix;
    const SmartPointer<const Preconditioner> preconditioner;

    const MPI_Comm *                      mpi_communicator;
    mutable TrilinosWrappers::MPI::Vector tmp;
  };


  template <class Matrix, class Preconditioner>
  InverseMatrix<Matrix, Preconditioner>::InverseMatrix(
    const Matrix &        m,
    const Preconditioner &preconditioner,
    const IndexSet &      locally_owned,
    const MPI_Comm &      mpi_communicator)
    : matrix(&m)
    , preconditioner(&preconditioner)
    , mpi_communicator(&mpi_communicator)
    , tmp(locally_owned, mpi_communicator)
  {}



  template <class Matrix, class Preconditioner>
  void
  InverseMatrix<Matrix, Preconditioner>::vmult(
    TrilinosWrappers::MPI::Vector &      dst,
    const TrilinosWrappers::MPI::Vector &src) const
  {
    SolverControl                           solver_control(src.size(),
                                 1e-6 * src.l2_norm(),
                                 false,
                                 false);
    SolverCG<TrilinosWrappers::MPI::Vector> cg(solver_control);

    tmp = 0.;
    cg.solve(*matrix, tmp, src, *preconditioner);
    dst = tmp;
  }



  template <class Preconditioner>
  class SchurComplement : public TrilinosWrappers::SparseMatrix
  {
  public:
    SchurComplement(const TrilinosWrappers::BlockSparseMatrix &system_matrix,
                    const InverseMatrix<TrilinosWrappers::SparseMatrix,
                                        Preconditioner> &      A_inverse,
                    const IndexSet &                           owned_pres,
                    const MPI_Comm &mpi_communicator);

    void
    vmult(TrilinosWrappers::MPI::Vector &      dst,
          const TrilinosWrappers::MPI::Vector &src) const;

  private:
    const SmartPointer<const TrilinosWrappers::BlockSparseMatrix> system_matrix;
    const SmartPointer<
      const InverseMatrix<TrilinosWrappers::SparseMatrix, Preconditioner>>
                                          A_inverse;
    mutable TrilinosWrappers::MPI::Vector tmp1, tmp2;
  };



  template <class Preconditioner>
  SchurComplement<Preconditioner>::SchurComplement(
    const TrilinosWrappers::BlockSparseMatrix &system_matrix,
    const InverseMatrix<TrilinosWrappers::SparseMatrix, Preconditioner>
      &             A_inverse,
    const IndexSet &owned_vel,
    const MPI_Comm &mpi_communicator)
    : system_matrix(&system_matrix)
    , A_inverse(&A_inverse)
    , tmp1(owned_vel, mpi_communicator)
    , tmp2(tmp1)
  {}


  template <class Preconditioner>
  void
  SchurComplement<Preconditioner>::vmult(
    TrilinosWrappers::MPI::Vector &      dst,
    const TrilinosWrappers::MPI::Vector &src) const
  {
    system_matrix->block(0, 1).vmult(tmp1, src);
    A_inverse->vmult(tmp2, tmp1);
    system_matrix->block(1, 0).vmult(dst, tmp2);
  }



  template <int dim>
  StokesProblem<dim>::StokesProblem (const unsigned int degree)
    :
    degree (degree),
    mpi_communicator (MPI_COMM_WORLD),
    triangulation (mpi_communicator/*,
                   Triangulation<dim>::maximum_smoothing*/),
    fe (FE_Q<dim>(degree+1), dim,
        FE_Q<dim>(degree), 1),
    dof_handler (triangulation),
    pcout (Utilities::MPI::this_mpi_process(mpi_communicator)
           == 0
           ?
           deallog.get_file_stream()
           :
           std::cout,
           (Utilities::MPI::this_mpi_process(mpi_communicator)
            == 0)),
    rot_matrix(dim)
  {
    rot_matrix[0][dim - 1] = -1.;
    rot_matrix[dim - 1][0] = 1.;
    if (dim == 3)
      rot_matrix[1][1] = 1.;
    offset[0] = 0;
    offset[1] = 0;
  }



  template <int dim>
  void
  StokesProblem<dim>::setup_dofs()
  {
    dof_handler.distribute_dofs(fe);

    std::vector<unsigned int> block_component(dim + 1, 0);
    block_component[dim] = 1;
    DoFRenumbering::component_wise(dof_handler, block_component);

    const std::vector<types::global_dof_index> dofs_per_block =
      DoFTools::count_dofs_per_fe_block(dof_handler, block_component);
    const unsigned int n_u = dofs_per_block[0], n_p = dofs_per_block[1];

    {
      owned_partitioning.clear();
      IndexSet locally_owned_dofs = dof_handler.locally_owned_dofs();
      owned_partitioning.push_back(locally_owned_dofs.get_view(0, n_u));
      owned_partitioning.push_back(locally_owned_dofs.get_view(n_u, n_u + n_p));

      relevant_partitioning.clear();
      IndexSet locally_relevant_dofs;
      DoFTools::extract_locally_relevant_dofs(dof_handler,
                                              locally_relevant_dofs);
      relevant_partitioning.push_back(locally_relevant_dofs.get_view(0, n_u));
      relevant_partitioning.push_back(
        locally_relevant_dofs.get_view(n_u, n_u + n_p));

      constraints.clear();
      constraints.reinit(locally_relevant_dofs);

      FEValuesExtractors::Vector velocities(0);
      FEValuesExtractors::Scalar pressure(dim);

      DoFTools::make_hanging_node_constraints(dof_handler, constraints);

      std::vector<
        GridTools::PeriodicFacePair<typename DoFHandler<dim>::cell_iterator>>
        periodicity_vector;

      std::vector<unsigned int> first_vector_components;
      first_vector_components.push_back(0);

      GridTools::collect_periodic_faces(
        dof_handler, 4, 0, 1, periodicity_vector, offset, rot_matrix);

      DoFTools::make_periodicity_constraints<DoFHandler<dim>>(
        periodicity_vector,
        constraints,
        fe.component_mask(velocities),
        first_vector_components);

      VectorTools::interpolate_boundary_values(dof_handler,
                                               1,
                                               Functions::ZeroFunction<dim>(
                                                 dim + 1),
                                               constraints,
                                               fe.component_mask(velocities));

      VectorTools::interpolate_boundary_values(dof_handler,
                                               2,
                                               Functions::ZeroFunction<dim>(
                                                 dim + 1),
                                               constraints,
                                               fe.component_mask(velocities));

      VectorTools::interpolate_boundary_values(dof_handler,
                                               3,
                                               Functions::ZeroFunction<dim>(
                                                 dim + 1),
                                               constraints,
                                               fe.component_mask(velocities));

      VectorTools::interpolate_boundary_values(dof_handler,
                                               5,
                                               Functions::ZeroFunction<dim>(
                                                 dim + 1),
                                               constraints,
                                               fe.component_mask(velocities));
    }
    constraints.close();

    const std::vector<IndexSet> &locally_owned_dofs =
      dof_handler.compute_locally_owned_dofs_per_processor();
    IndexSet locally_active_dofs;
    DoFTools::extract_locally_active_dofs(dof_handler, locally_active_dofs);
    AssertThrow(constraints.is_consistent_in_parallel(locally_owned_dofs,
                                                      locally_active_dofs,
                                                      mpi_communicator,
                                                      /*verbose*/ true),
                ExcInternalError());

    {
      TrilinosWrappers::BlockSparsityPattern bsp(owned_partitioning,
                                                 owned_partitioning,
                                                 relevant_partitioning,
                                                 mpi_communicator);

      DoFTools::make_sparsity_pattern(dof_handler,
                                      bsp,
                                      constraints,
                                      false,
                                      Utilities::MPI::this_mpi_process(
                                        mpi_communicator));

      bsp.compress();

      system_matrix.reinit(bsp);
    }

    system_rhs.reinit(owned_partitioning, mpi_communicator);
    solution.reinit(owned_partitioning,
                    relevant_partitioning,
                    mpi_communicator);
  }



  template <int dim>
  void
  StokesProblem<dim>::assemble_system()
  {
    system_matrix = 0.;
    system_rhs    = 0.;

    QGauss<dim> quadrature_formula(degree + 2);

    FEValues<dim> fe_values(fe,
                            quadrature_formula,
                            update_values | update_quadrature_points |
                              update_JxW_values | update_gradients);

    const unsigned int dofs_per_cell = fe.dofs_per_cell;

    const unsigned int n_q_points = quadrature_formula.size();

    FullMatrix<double> local_matrix(dofs_per_cell, dofs_per_cell);
    Vector<double>     local_rhs(dofs_per_cell);

    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

    const FEValuesExtractors::Vector velocities(0);
    const FEValuesExtractors::Scalar pressure(dim);

    std::vector<SymmetricTensor<2, dim>> symgrad_phi_u(dofs_per_cell);
    std::vector<double>                  div_phi_u(dofs_per_cell);
    std::vector<double>                  phi_p(dofs_per_cell);

    typename DoFHandler<dim>::active_cell_iterator cell =
                                                     dof_handler.begin_active(),
                                                   endc = dof_handler.end();
    for (; cell != endc; ++cell)
      if (cell->is_locally_owned())
        {
          fe_values.reinit(cell);
          local_matrix = 0;
          local_rhs    = 0;

          for (unsigned int q = 0; q < n_q_points; ++q)
            {
              for (unsigned int k = 0; k < dofs_per_cell; ++k)
                {
                  symgrad_phi_u[k] =
                    fe_values[velocities].symmetric_gradient(k, q);
                  div_phi_u[k] = fe_values[velocities].divergence(k, q);
                  phi_p[k]     = fe_values[pressure].value(k, q);
                }

              for (unsigned int i = 0; i < dofs_per_cell; ++i)
                {
                  for (unsigned int j = 0; j <= i; ++j)
                    {
                      local_matrix(i, j) +=
                        (2 * (symgrad_phi_u[i] * symgrad_phi_u[j]) -
                         div_phi_u[i] * phi_p[j] - phi_p[i] * div_phi_u[j] +
                         phi_p[i] * phi_p[j]) *
                        fe_values.JxW(q);
                    }

                  const unsigned int component_i =
                    fe.system_to_component_index(i).first;
                  if (component_i == 0)
                    local_rhs(i) +=
                      fe_values.shape_value(i, q) * 1 * fe_values.JxW(q);
                }
            }


          for (unsigned int i = 0; i < dofs_per_cell; ++i)
            for (unsigned int j = i + 1; j < dofs_per_cell; ++j)
              local_matrix(i, j) = local_matrix(j, i);

          cell->get_dof_indices(local_dof_indices);
          constraints.distribute_local_to_global(local_matrix,
                                                 local_rhs,
                                                 local_dof_indices,
                                                 system_matrix,
                                                 system_rhs);
        }

    system_matrix.compress(VectorOperation::add);
    system_rhs.compress(VectorOperation::add);
  }



  template <int dim>
  void
  StokesProblem<dim>::solve()
  {
    TrilinosWrappers::PreconditionJacobi A_preconditioner;
    A_preconditioner.initialize(system_matrix.block(0, 0));

    const InverseMatrix<TrilinosWrappers::SparseMatrix,
                        TrilinosWrappers::PreconditionJacobi>
      A_inverse(system_matrix.block(0, 0),
                A_preconditioner,
                owned_partitioning[0],
                mpi_communicator);

    TrilinosWrappers::MPI::BlockVector tmp(owned_partitioning,
                                           mpi_communicator);

    {
      TrilinosWrappers::MPI::Vector schur_rhs(owned_partitioning[1],
                                              mpi_communicator);
      A_inverse.vmult(tmp.block(0), system_rhs.block(0));
      system_matrix.block(1, 0).vmult(schur_rhs, tmp.block(0));
      schur_rhs -= system_rhs.block(1);

      SchurComplement<TrilinosWrappers::PreconditionJacobi> schur_complement(
        system_matrix, A_inverse, owned_partitioning[0], mpi_communicator);

      SolverControl solver_control(solution.block(1).size(),
                                   1e-6 * schur_rhs.l2_norm(),
                                   false,
                                   false);
      SolverCG<TrilinosWrappers::MPI::Vector> cg(solver_control);

      TrilinosWrappers::PreconditionAMG preconditioner;
      preconditioner.initialize(system_matrix.block(1, 1));


      cg.solve(schur_complement, tmp.block(1), schur_rhs, preconditioner);

      constraints.distribute(tmp);
      solution.block(1) = tmp.block(1);
    }

    {
      system_matrix.block(0, 1).vmult(tmp.block(0), tmp.block(1));
      tmp.block(0) *= -1;
      tmp.block(0) += system_rhs.block(0);

      A_inverse.vmult(tmp.block(0), tmp.block(0));

      constraints.distribute(tmp);
      solution.block(0) = tmp.block(0);
    }
  }

  template <int dim>
  void
  StokesProblem<dim>::get_point_value(const Point<dim> point,
                                      const int        proc,
                                      Vector<double> & value) const
  {
    try
      {
        const typename DoFHandler<dim>::active_cell_iterator cell =
          GridTools::find_active_cell_around_point(dof_handler, point);

        if (cell->is_locally_owned())
          VectorTools::point_value(dof_handler, solution, point, value);
      }
    catch (GridTools::ExcPointNotFound<dim> &p)
      {
        pcout << "Point: " << point << " is not inside a cell!" << std::endl;
      }


    std::vector<double> tmp(value.size());
    for (unsigned int i = 0; i < value.size(); ++i)
      tmp[i] = value[i];

    std::vector<double> tmp2(value.size());
    MPI_Reduce(&(tmp[0]),
               &(tmp2[0]),
               value.size(),
               MPI_DOUBLE,
               MPI_SUM,
               proc,
               mpi_communicator);

    for (unsigned int i = 0; i < value.size(); ++i)
      value[i] = tmp2[i];
  }

  template <int dim>
  void
  StokesProblem<dim>::check_periodicity() const
  {
    AssertThrow(false, ExcNotImplemented());
  }

  template <>
  void
  StokesProblem<3>::check_periodicity() const
  {
    const unsigned int  dim = 3;
    Quadrature<dim - 1> q_dummy(
      fe.base_element(0).get_unit_face_support_points());
    FEFaceValues<dim> fe_face_values(fe, q_dummy, update_quadrature_points);

    DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active();
    DoFHandler<dim>::active_cell_iterator endc = dof_handler.end();

    std::vector<double> local_quad_points_first, local_quad_points_second;
    std::vector<types::global_dof_index> gdi(fe.dofs_per_cell);

    // first collect all points locally
    for (; cell != endc; ++cell)
      if (cell->is_locally_owned() && cell->at_boundary())
        for (const unsigned int face_no : GeometryInfo<dim>::face_indices())
          {
            const typename DoFHandler<dim>::face_iterator face =
              cell->face(face_no);
            if (face->at_boundary())
              {
                if (face->boundary_id() == 0)
                  {
                    fe_face_values.reinit(cell, face_no);
                    const std::vector<Point<dim>> tmp_points =
                      fe_face_values.get_quadrature_points();
                    for (unsigned int i = 0; i < tmp_points.size(); ++i)
                      for (unsigned int c = 0; c < dim; ++c)
                        local_quad_points_first.push_back(tmp_points[i](c));
                  }
                else if (face->boundary_id() == 4)
                  {
                    fe_face_values.reinit(cell, face_no);
                    const std::vector<Point<dim>> tmp_points =
                      fe_face_values.get_quadrature_points();
                    for (unsigned int i = 0; i < tmp_points.size(); ++i)
                      for (unsigned int c = 0; c < dim; ++c)
                        local_quad_points_second.push_back(tmp_points[i](c));
                  }
              }
          }
    // next exchange them
    std::vector<double> global_quad_points_first;
    {
      // how many elements are sent from which process?
      unsigned int       n_my_elements = local_quad_points_first.size();
      const unsigned int n_processes =
        Utilities::MPI::n_mpi_processes(mpi_communicator);
      std::vector<int> n_elements(n_processes);
      MPI_Allgather(&n_my_elements,
                    1,
                    MPI_INT,
                    &n_elements[0],
                    1,
                    MPI_INT,
                    mpi_communicator);
      std::vector<int> displacements(n_processes + 1);
      displacements[0] = 0;
      for (unsigned int i = 1; i <= n_processes; ++i)
        displacements[i] = displacements[i - 1] + n_elements[i - 1];

      global_quad_points_first.resize(displacements[n_processes]);
      MPI_Allgatherv(&local_quad_points_first[0],
                     n_my_elements,
                     MPI_DOUBLE,
                     &global_quad_points_first[0],
                     &n_elements[0],
                     &displacements[0],
                     MPI_DOUBLE,
                     mpi_communicator);
    }
    std::vector<double> global_quad_points_second;
    {
      // how many elements are sent from which process?
      unsigned int       n_my_elements = local_quad_points_second.size();
      const unsigned int n_processes =
        Utilities::MPI::n_mpi_processes(mpi_communicator);
      std::vector<int> n_elements(n_processes);
      MPI_Allgather(&n_my_elements,
                    1,
                    MPI_INT,
                    &n_elements[0],
                    1,
                    MPI_INT,
                    mpi_communicator);
      std::vector<int> displacements(n_processes + 1);
      displacements[0] = 0;
      for (unsigned int i = 1; i <= n_processes; ++i)
        displacements[i] = displacements[i - 1] + n_elements[i - 1];

      global_quad_points_second.resize(displacements[n_processes]);
      MPI_Allgatherv(&local_quad_points_second[0],
                     n_my_elements,
                     MPI_DOUBLE,
                     &global_quad_points_second[0],
                     &n_elements[0],
                     &displacements[0],
                     MPI_DOUBLE,
                     mpi_communicator);
    }

    // consider points on the boundary with the first indicator
    // get the values of the solution at these points and compare to the values
    // at the corresponding points on the boundary with the second indicator
    {
      for (unsigned int i = 0; i < global_quad_points_first.size(); i += dim)
        {
          Vector<double> value_1(dim + 1);
          Vector<double> value_2(dim + 1);

          Point<dim>     point_1, point_2;
          Vector<double> vector_point_1(dim);
          for (unsigned int c = 0; c < dim; ++c)
            {
              vector_point_1(c) = global_quad_points_first[i + c];
              point_1(c)        = vector_point_1(c);
            }
          Vector<double> vector_point_2(dim);
          rot_matrix.Tvmult(vector_point_2, vector_point_1);
          for (unsigned int c = 0; c < dim; ++c)
            point_2(c) = vector_point_2(c) + offset[c];

          get_point_value(point_1, 0, value_1);
          get_point_value(point_2, 0, value_2);

          if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
            {
              Vector<double> vel_value_1(dim);
              Vector<double> vel_value_2(dim);
              Vector<double> expected_vel(dim);
              for (unsigned int c = 0; c < dim; ++c)
                {
                  vel_value_1(c) = value_1(c);
                  vel_value_2(c) = value_2(c);
                }
              rot_matrix.Tvmult(expected_vel, vel_value_1);
              pcout << "Point 1: " << point_1 << "\t Point 2: " << point_2
                    << std::endl;
              Vector<double> error = expected_vel;
              error -= vel_value_2;
              if (std::abs(error.l2_norm()) > 1e-8)
                {
                  pcout << "first first_rotated second" << std::endl;
                  for (unsigned int c = 0; c < dim; ++c)
                    pcout << vel_value_1(c) << "\t" << expected_vel(c) << "\t"
                          << vel_value_2(c) << std::endl;
                  pcout << std::endl;
                  Assert(false, ExcInternalError());
                }
            }
        }
    }

    // consider points on the boundary with the second indicator
    // get the values of the solution at these points and compare to the values
    // at the corresponding points on the boundary with the first indicator
    {
      for (unsigned int i = 0; i < global_quad_points_second.size(); i += dim)
        {
          Vector<double> value_1(dim + 1);
          Vector<double> value_2(dim + 1);

          Point<dim>     point_1, point_2;
          Vector<double> vector_point_1(dim);
          for (unsigned int c = 0; c < dim; ++c)
            {
              vector_point_1(c) = global_quad_points_second[i + c];
              point_1(c)        = vector_point_1(c);
            }
          Vector<double> vector_point_2(dim);
          for (unsigned int c = 0; c < dim; ++c)
            vector_point_1(c) -= offset[c];
          rot_matrix.vmult(vector_point_2, vector_point_1);
          for (unsigned int c = 0; c < dim; ++c)
            point_2(c) = vector_point_2(c);

          get_point_value(point_1, 0, value_1);
          get_point_value(point_2, 0, value_2);

          if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
            {
              Vector<double> vel_value_1(dim);
              Vector<double> vel_value_2(dim);
              Vector<double> expected_vel(dim);
              for (unsigned int c = 0; c < dim; ++c)
                {
                  vel_value_1(c) = value_1(c);
                  vel_value_2(c) = value_2(c);
                }
              rot_matrix.vmult(expected_vel, vel_value_1);
              pcout << "Point 1: " << point_1 << "\t Point 2: " << point_2
                    << std::endl;
              Vector<double> error = expected_vel;
              error -= vel_value_2;
              if (std::abs(error.l2_norm()) > 1e-8)
                {
                  pcout << "first first_rotated second" << std::endl;
                  for (unsigned int c = 0; c < dim; ++c)
                    pcout << vel_value_1(c) << "\t" << expected_vel(c) << "\t"
                          << vel_value_2(c) << std::endl;
                  pcout << std::endl;
                  Assert(false, ExcInternalError());
                }
            }
        }
    }
  }


  template <int dim>
  void
  StokesProblem<dim>::output_results(const unsigned int refinement_cycle) const
  {
    std::vector<std::string> solution_names(dim, "velocity");
    solution_names.push_back("pressure");

    std::vector<DataComponentInterpretation::DataComponentInterpretation>
      data_component_interpretation(
        dim, DataComponentInterpretation::component_is_part_of_vector);
    data_component_interpretation.push_back(
      DataComponentInterpretation::component_is_scalar);

    DataOut<dim> data_out;
    data_out.attach_dof_handler(dof_handler);
    data_out.add_data_vector(solution,
                             solution_names,
                             DataOut<dim>::type_dof_data,
                             data_component_interpretation);
    Vector<float> subdomain(triangulation.n_active_cells());
    for (unsigned int i = 0; i < subdomain.size(); ++i)
      subdomain(i) = triangulation.locally_owned_subdomain();
    data_out.add_data_vector(subdomain, "subdomain");
    data_out.build_patches(degree + 1);

    std::ostringstream filename;
    filename
      << "solution-" << Utilities::int_to_string(refinement_cycle, 2) << "."
      << Utilities::int_to_string(triangulation.locally_owned_subdomain(), 2)
      << ".vtu";

    std::ofstream output(filename.str().c_str());
    data_out.write_vtu(output);

    if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
      {
        std::vector<std::string> filenames;
        for (unsigned int i = 0;
             i < Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);
             ++i)
          filenames.push_back(std::string("solution-") +
                              Utilities::int_to_string(refinement_cycle, 2) +
                              "." + Utilities::int_to_string(i, 2) + ".vtu");
        const std::string pvtu_master_filename =
          ("solution-" + Utilities::int_to_string(refinement_cycle, 2) +
           ".pvtu");
        std::ofstream pvtu_master(pvtu_master_filename.c_str());
        data_out.write_pvtu_record(pvtu_master, filenames);
      }
  }



  template <int dim>
  void
  StokesProblem<dim>::refine_mesh()
  {
    Vector<float> estimated_error_per_cell(triangulation.n_active_cells());

    FEValuesExtractors::Scalar pressure(dim);
    KellyErrorEstimator<dim>::estimate(
      dof_handler,
      QGauss<dim - 1>(degree + 1),
      std::map<types::boundary_id, const Function<dim> *>(),
      solution,
      estimated_error_per_cell,
      fe.component_mask(pressure));

    parallel::distributed::GridRefinement::refine_and_coarsen_fixed_number(
      triangulation, estimated_error_per_cell, 0.3, 0.0);
    triangulation.execute_coarsening_and_refinement();
  }



  template <int dim>
  void
  StokesProblem<dim>::run()
  {
    GridGenerator::hyper_cube(triangulation, 0., 1., true);

    std::vector<GridTools::PeriodicFacePair<
      typename parallel::distributed::Triangulation<dim>::cell_iterator>>
      periodicity_vector;

    GridTools::collect_periodic_faces(
      triangulation, 4, 0, 1, periodicity_vector, offset, rot_matrix);
    triangulation.add_periodicity(periodicity_vector);

    triangulation.refine_global(1);

    for (unsigned int refinement_cycle = 0; refinement_cycle < 3;
         ++refinement_cycle)
      {
        pcout << "Refinement cycle " << refinement_cycle << std::endl;

        if (refinement_cycle > 0)
          refine_mesh();

        setup_dofs();

        pcout << "   Assembling..." << std::endl << std::flush;
        assemble_system();

        pcout << "   Solving..." << std::endl << std::flush;
        solve();

        //        output_results (refinement_cycle);
        check_periodicity();

        pcout << std::endl;
      }
  }
} // namespace Step22



int
main(int argc, char *argv[])
{
  try
    {
      using namespace Step22;

      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

      if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
        {
          std::ofstream logfile("output");
          deallog.attach(logfile, false);
          {
            StokesProblem<3> flow_problem(1);
            flow_problem.run();
          }
        }
      else
        {
          StokesProblem<3> flow_problem(1);
          flow_problem.run();
        }
    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Exception on processing: " << std::endl
                << exc.what() << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;

      return 1;
    }
  catch (...)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Unknown exception!" << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }

  return 0;
}
