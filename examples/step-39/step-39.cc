/* ------------------------------------------------------------------------
 *
 * SPDX-License-Identifier: LGPL-2.1-or-later
 * Copyright (C) 2010 - 2024 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * Part of the source code is dual licensed under Apache-2.0 WITH
 * LLVM-exception OR LGPL-2.1-or-later. Detailed license information
 * governing the source code and code contributions can be found in
 * LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
 *
 * ------------------------------------------------------------------------
 *
 * Author: Guido Kanschat, Texas A&M University, 2009
 */


// The include files for the linear algebra: A regular SparseMatrix, which in
// turn will include the necessary files for SparsityPattern and Vector
// classes.
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/precondition_block.h>
#include <deal.II/lac/block_vector.h>

// Include files for setting up the mesh
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>

// Include files for FiniteElement classes and DoFHandler.
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_dgp.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/mapping_q1.h>

// The include files for using the MeshWorker framework
#include <deal.II/meshworker/dof_info.h>
#include <deal.II/meshworker/integration_info.h>
#include <deal.II/meshworker/assembler.h>
#include <deal.II/meshworker/loop.h>

// The include file for local integrators associated with the Laplacian
#include <deal.II/integrators/laplace.h>

// Support for multigrid methods
#include <deal.II/multigrid/mg_tools.h>
#include <deal.II/multigrid/multigrid.h>
#include <deal.II/multigrid/mg_matrix.h>
#include <deal.II/multigrid/mg_transfer.h>
#include <deal.II/multigrid/mg_coarse.h>
#include <deal.II/multigrid/mg_smoother.h>

// Finally, we take our exact solution from the library as well as quadrature
// and additional tools.
#include <deal.II/base/function_lib.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/data_out.h>

#include <iostream>
#include <fstream>

// All classes of the deal.II library are in the namespace dealii. In order to
// save typing, we tell the compiler to search names in there as well.
namespace Step39
{
  using namespace dealii;

  // This is the function we use to set the boundary values and also the exact
  // solution we compare to.
  Functions::SlitSingularityFunction<2> exact_solution;

  // @sect3{The local integrators}

  // The MeshWorker::loop() function separates what needs to be done for
  // local integration, from the loops over cells and
  // faces. It does this by calling functions that integrate over a cell,
  // a boundary face, or an interior face, and letting them create the
  // local contributions and then in a separate step calling a function
  // that moves these local contributions into the global objects.
  // We will use this approach for computing the
  // matrices, the right hand side, the error estimator, and the actual
  // error computation in the functions below. For each of these operations,
  // we provide a namespace that contains a set of functions for cell, boundary,
  // and interior face contributions.
  //
  // All the information needed for these local integration is provided by
  // MeshWorker::DoFInfo<dim> and MeshWorker::IntegrationInfo<dim>. In each
  // case, the functions' signatures is fixed: MeshWorker::loop() wants to call
  // functions with a specific set of arguments, so the signature of the
  // functions cannot be changed.

  // The first namespace defining local integrators is responsible for
  // assembling the global matrix as well as the level matrices.
  // On each cell, we integrate the Dirichlet form as well as the
  // Nitsche boundary conditions and the interior penalty fluxes between
  // cells.
  //
  // The boundary and flux terms need a penalty parameter, which should be
  // adjusted to the cell size and the polynomial degree. We compute it
  // in two steps: First, we compute on each cell
  // $K_i$ the value $P_i = p_i(p_i+1)/h_i$, where
  // $p_i$ is the polynomial degree on cell $K_i$ and $h_i$ is the length of
  // $K_i$ orthogonal to the current face. Second, if exactly one of the two
  // cells adjacent to the face has children, its penalty is multiplied
  // by two (to account for the fact that the mesh size $h_i$ there is
  // only half that previously computed); it is possible that both adjacent
  // cells are refined, in which case we are integrating over a non-active
  // face and no adjustment is necessary. Finally, we return the average
  // of the two penalty values.
  namespace MatrixIntegrator
  {
    template <int dim>
    double ip_penalty_factor(const MeshWorker::DoFInfo<dim> &dinfo1,
                             const MeshWorker::DoFInfo<dim> &dinfo2,
                             unsigned int                    deg1,
                             unsigned int                    deg2)
    {
      const unsigned int normal1 =
        GeometryInfo<dim>::unit_normal_direction[dinfo1.face_number];
      const unsigned int normal2 =
        GeometryInfo<dim>::unit_normal_direction[dinfo2.face_number];
      const unsigned int deg1sq = (deg1 == 0) ? 1 : deg1 * (deg1 + 1);
      const unsigned int deg2sq = (deg2 == 0) ? 1 : deg2 * (deg2 + 1);

      double penalty1 = deg1sq / dinfo1.cell->extent_in_direction(normal1);
      double penalty2 = deg2sq / dinfo2.cell->extent_in_direction(normal2);
      if (dinfo1.cell->has_children() && !dinfo2.cell->has_children())
        penalty1 *= 2;
      else if (!dinfo1.cell->has_children() && dinfo2.cell->has_children())
        penalty2 *= 2;

      const double penalty = 0.5 * (penalty1 + penalty2);
      return penalty;
    }


    template <int dim>
    void cell(MeshWorker::DoFInfo<dim>         &dinfo,
              MeshWorker::IntegrationInfo<dim> &info)
    {
      FullMatrix<double> &M = dinfo.matrix(0, false).matrix;

      for (unsigned int k = 0; k < info.fe_values().n_quadrature_points; ++k)
        {
          const double dx = info.fe_values().JxW(k);

          for (unsigned int i = 0; i < info.fe_values().dofs_per_cell; ++i)
            {
              const double Mii = (info.fe_values().shape_grad(i, k) *
                                  info.fe_values().shape_grad(i, k) * dx);

              M(i, i) += Mii;

              for (unsigned int j = i + 1; j < info.fe_values().dofs_per_cell;
                   ++j)
                {
                  const double Mij = info.fe_values().shape_grad(j, k) *
                                     info.fe_values().shape_grad(i, k) * dx;

                  M(i, j) += Mij;
                  M(j, i) += Mij;
                }
            }
        }
    }


    // Boundary faces use the Nitsche method to impose boundary values:
    template <int dim>
    void boundary(MeshWorker::DoFInfo<dim>         &dinfo,
                  MeshWorker::IntegrationInfo<dim> &info)
    {
      const FEValuesBase<dim> &fe_face_values = info.fe_values(0);

      FullMatrix<double> &M = dinfo.matrix(0, false).matrix;
      AssertDimension(M.n(), fe_face_values.dofs_per_cell);
      AssertDimension(M.m(), fe_face_values.dofs_per_cell);

      const unsigned int polynomial_degree =
        info.fe_values(0).get_fe().tensor_degree();

      const double ip_penalty =
        ip_penalty_factor(dinfo, dinfo, polynomial_degree, polynomial_degree);

      for (unsigned int k = 0; k < fe_face_values.n_quadrature_points; ++k)
        {
          const double          dx = fe_face_values.JxW(k);
          const Tensor<1, dim> &n  = fe_face_values.normal_vector(k);

          for (unsigned int i = 0; i < fe_face_values.dofs_per_cell; ++i)
            for (unsigned int j = 0; j < fe_face_values.dofs_per_cell; ++j)
              M(i, j) += (2. * fe_face_values.shape_value(i, k) * ip_penalty *
                            fe_face_values.shape_value(j, k) -
                          (n * fe_face_values.shape_grad(i, k)) *
                            fe_face_values.shape_value(j, k) -
                          (n * fe_face_values.shape_grad(j, k)) *
                            fe_face_values.shape_value(i, k)) *
                         dx;
        }
    }

    // Interior faces use the interior penalty method:
    template <int dim>
    void face(MeshWorker::DoFInfo<dim>         &dinfo1,
              MeshWorker::DoFInfo<dim>         &dinfo2,
              MeshWorker::IntegrationInfo<dim> &info1,
              MeshWorker::IntegrationInfo<dim> &info2)
    {
      const FEValuesBase<dim> &fe_face_values_1 = info1.fe_values(0);
      const FEValuesBase<dim> &fe_face_values_2 = info2.fe_values(0);

      FullMatrix<double> &M11 = dinfo1.matrix(0, false).matrix;
      FullMatrix<double> &M12 = dinfo1.matrix(0, true).matrix;
      FullMatrix<double> &M21 = dinfo2.matrix(0, true).matrix;
      FullMatrix<double> &M22 = dinfo2.matrix(0, false).matrix;

      AssertDimension(M11.n(), fe_face_values_1.dofs_per_cell);
      AssertDimension(M11.m(), fe_face_values_1.dofs_per_cell);
      AssertDimension(M12.n(), fe_face_values_1.dofs_per_cell);
      AssertDimension(M12.m(), fe_face_values_1.dofs_per_cell);
      AssertDimension(M21.n(), fe_face_values_1.dofs_per_cell);
      AssertDimension(M21.m(), fe_face_values_1.dofs_per_cell);
      AssertDimension(M22.n(), fe_face_values_1.dofs_per_cell);
      AssertDimension(M22.m(), fe_face_values_1.dofs_per_cell);

      const unsigned int polynomial_degree =
        info1.fe_values(0).get_fe().tensor_degree();
      const double ip_penalty =
        ip_penalty_factor(dinfo1, dinfo2, polynomial_degree, polynomial_degree);

      const double nui = 1.;
      const double nue = 1.;
      const double nu  = .5 * (nui + nue);

      for (unsigned int k = 0; k < fe_face_values_1.n_quadrature_points; ++k)
        {
          const double          dx = fe_face_values_1.JxW(k);
          const Tensor<1, dim> &n  = fe_face_values_1.normal_vector(k);

          for (unsigned int i = 0; i < fe_face_values_1.dofs_per_cell; ++i)
            {
              for (unsigned int j = 0; j < fe_face_values_1.dofs_per_cell; ++j)
                {
                  const double vi   = fe_face_values_1.shape_value(i, k);
                  const double dnvi = n * fe_face_values_1.shape_grad(i, k);
                  const double ve   = fe_face_values_2.shape_value(i, k);
                  const double dnve = n * fe_face_values_2.shape_grad(i, k);
                  const double ui   = fe_face_values_1.shape_value(j, k);
                  const double dnui = n * fe_face_values_1.shape_grad(j, k);
                  const double ue   = fe_face_values_2.shape_value(j, k);
                  const double dnue = n * fe_face_values_2.shape_grad(j, k);

                  M11(i, j) += (-.5 * nui * dnvi * ui - .5 * nui * dnui * vi +
                                nu * ip_penalty * ui * vi) *
                               dx;
                  M12(i, j) += (.5 * nui * dnvi * ue - .5 * nue * dnue * vi -
                                nu * ip_penalty * vi * ue) *
                               dx;
                  M21(i, j) += (-.5 * nue * dnve * ui + .5 * nui * dnui * ve -
                                nu * ip_penalty * ui * ve) *
                               dx;
                  M22(i, j) += (.5 * nue * dnve * ue + .5 * nue * dnue * ve +
                                nu * ip_penalty * ue * ve) *
                               dx;
                }
            }
        }
    }
  } // namespace MatrixIntegrator

  // The second set of local integrators builds the right hand side. In our
  // example, the right hand side function is zero, such that only the boundary
  // condition is set here in weak form.
  namespace RHSIntegrator
  {
    template <int dim>
    void cell(MeshWorker::DoFInfo<dim> &, MeshWorker::IntegrationInfo<dim> &)
    {}


    template <int dim>
    void boundary(MeshWorker::DoFInfo<dim>         &dinfo,
                  MeshWorker::IntegrationInfo<dim> &info)
    {
      const FEValuesBase<dim> &fe           = info.fe_values();
      Vector<double>          &local_vector = dinfo.vector(0).block(0);

      std::vector<double> boundary_values(fe.n_quadrature_points);
      exact_solution.value_list(fe.get_quadrature_points(), boundary_values);

      const unsigned int degree  = fe.get_fe().tensor_degree();
      const double       penalty = 2. * degree * (degree + 1) *
                             dinfo.face->measure() / dinfo.cell->measure();

      for (unsigned k = 0; k < fe.n_quadrature_points; ++k)
        for (unsigned int i = 0; i < fe.dofs_per_cell; ++i)
          local_vector(i) +=
            (-penalty * fe.shape_value(i, k) // (-sigma * v_i(x_k)
             +
             fe.normal_vector(k) * fe.shape_grad(i, k)) // + n * grad v_i(x_k))
            * boundary_values[k] * fe.JxW(k);           // u^D(x_k) * dx
    }


    template <int dim>
    void face(MeshWorker::DoFInfo<dim> &,
              MeshWorker::DoFInfo<dim> &,
              MeshWorker::IntegrationInfo<dim> &,
              MeshWorker::IntegrationInfo<dim> &)
    {}
  } // namespace RHSIntegrator

  // The third local integrator is responsible for the contributions to the
  // error estimate. This is the standard energy estimator due to Karakashian
  // and Pascal (2003).
  // The cell contribution is the Laplacian of the discrete solution, since
  // the right hand side is zero.
  namespace Estimator
  {
    template <int dim>
    void cell(MeshWorker::DoFInfo<dim>         &dinfo,
              MeshWorker::IntegrationInfo<dim> &info)
    {
      const FEValuesBase<dim> &fe = info.fe_values();

      const std::vector<Tensor<2, dim>> &DDuh = info.hessians[0][0];
      for (unsigned k = 0; k < fe.n_quadrature_points; ++k)
        {
          const double t = dinfo.cell->diameter() * trace(DDuh[k]);
          dinfo.value(0) += t * t * fe.JxW(k);
        }
      dinfo.value(0) = std::sqrt(dinfo.value(0));
    }

    // At the boundary, we use simply a weighted form of the boundary residual,
    // namely the norm of the difference between the finite element solution and
    // the correct boundary condition.
    template <int dim>
    void boundary(MeshWorker::DoFInfo<dim>         &dinfo,
                  MeshWorker::IntegrationInfo<dim> &info)
    {
      const FEValuesBase<dim> &fe = info.fe_values();

      std::vector<double> boundary_values(fe.n_quadrature_points);
      exact_solution.value_list(fe.get_quadrature_points(), boundary_values);

      const std::vector<double> &uh = info.values[0][0];

      const unsigned int degree  = fe.get_fe().tensor_degree();
      const double       penalty = 2. * degree * (degree + 1) *
                             dinfo.face->measure() / dinfo.cell->measure();

      for (unsigned k = 0; k < fe.n_quadrature_points; ++k)
        {
          const double diff = boundary_values[k] - uh[k];
          dinfo.value(0) += penalty * diff * diff * fe.JxW(k);
        }
      dinfo.value(0) = std::sqrt(dinfo.value(0));
    }


    // Finally, on interior faces, the estimator consists of the jumps of the
    // solution and its normal derivative, weighted appropriately.
    template <int dim>
    void face(MeshWorker::DoFInfo<dim>         &dinfo1,
              MeshWorker::DoFInfo<dim>         &dinfo2,
              MeshWorker::IntegrationInfo<dim> &info1,
              MeshWorker::IntegrationInfo<dim> &info2)
    {
      const FEValuesBase<dim>           &fe   = info1.fe_values();
      const std::vector<double>         &uh1  = info1.values[0][0];
      const std::vector<double>         &uh2  = info2.values[0][0];
      const std::vector<Tensor<1, dim>> &Duh1 = info1.gradients[0][0];
      const std::vector<Tensor<1, dim>> &Duh2 = info2.gradients[0][0];

      const unsigned int degree = fe.get_fe().tensor_degree();
      const double       penalty1 =
        degree * (degree + 1) * dinfo1.face->measure() / dinfo1.cell->measure();
      const double penalty2 =
        degree * (degree + 1) * dinfo2.face->measure() / dinfo2.cell->measure();
      const double penalty = penalty1 + penalty2;
      const double h       = dinfo1.face->measure();

      for (unsigned k = 0; k < fe.n_quadrature_points; ++k)
        {
          const double diff1 = uh1[k] - uh2[k];
          const double diff2 =
            fe.normal_vector(k) * Duh1[k] - fe.normal_vector(k) * Duh2[k];
          dinfo1.value(0) +=
            (penalty * diff1 * diff1 + h * diff2 * diff2) * fe.JxW(k);
        }
      dinfo1.value(0) = std::sqrt(dinfo1.value(0));
      dinfo2.value(0) = dinfo1.value(0);
    }
  } // namespace Estimator

  // Finally we have an integrator for the error. Since the energy norm for
  // discontinuous Galerkin problems not only involves the difference of the
  // gradient inside the cells, but also the jump terms across faces and at
  // the boundary, we cannot just use VectorTools::integrate_difference().
  // Instead, we use the MeshWorker interface to compute the error ourselves.

  // There are several different ways to define this energy norm, but all of
  // them are equivalent to each other uniformly with mesh size (some not
  // uniformly with polynomial degree). Here, we choose @f[ \|u\|_{1,h} =
  // \sum_{K\in \mathbb T_h} \|\nabla u\|_K^2 + \sum_{F \in F_h^i}
  // 4\sigma_F\|\average{ u \mathbf n}\|^2_F + \sum_{F \in F_h^b}
  // 2\sigma_F\|u\|^2_F @f]
  //
  // Below, the first function is, as always, the integration
  // on cells. There is currently no good
  // interface in MeshWorker that would allow us to access values of regular
  // functions in the quadrature points. Thus, we have to create the vectors
  // for the exact function's values and gradients inside the cell
  // integrator. After that, everything is as before and we just add up the
  // squares of the differences.

  // Additionally to computing the error in the energy norm, we use the
  // capability of the mesh worker to compute two functionals at the same time
  // and compute the <i>L<sup>2</sup></i>-error in the same loop. Obviously,
  // this one does not have any jump terms and only appears in the integration
  // on cells.

  namespace ErrorIntegrator
  {
    template <int dim>
    void cell(MeshWorker::DoFInfo<dim>         &dinfo,
              MeshWorker::IntegrationInfo<dim> &info)
    {
      const FEValuesBase<dim>    &fe = info.fe_values();
      std::vector<Tensor<1, dim>> exact_gradients(fe.n_quadrature_points);
      std::vector<double>         exact_values(fe.n_quadrature_points);

      exact_solution.gradient_list(fe.get_quadrature_points(), exact_gradients);
      exact_solution.value_list(fe.get_quadrature_points(), exact_values);

      const std::vector<Tensor<1, dim>> &Duh = info.gradients[0][0];
      const std::vector<double>         &uh  = info.values[0][0];

      for (unsigned k = 0; k < fe.n_quadrature_points; ++k)
        {
          double sum = 0;
          for (unsigned int d = 0; d < dim; ++d)
            {
              const double diff = exact_gradients[k][d] - Duh[k][d];
              sum += diff * diff;
            }
          const double diff = exact_values[k] - uh[k];
          dinfo.value(0) += sum * fe.JxW(k);
          dinfo.value(1) += diff * diff * fe.JxW(k);
        }
      dinfo.value(0) = std::sqrt(dinfo.value(0));
      dinfo.value(1) = std::sqrt(dinfo.value(1));
    }


    template <int dim>
    void boundary(MeshWorker::DoFInfo<dim>         &dinfo,
                  MeshWorker::IntegrationInfo<dim> &info)
    {
      const FEValuesBase<dim> &fe = info.fe_values();

      std::vector<double> exact_values(fe.n_quadrature_points);
      exact_solution.value_list(fe.get_quadrature_points(), exact_values);

      const std::vector<double> &uh = info.values[0][0];

      const unsigned int degree  = fe.get_fe().tensor_degree();
      const double       penalty = 2. * degree * (degree + 1) *
                             dinfo.face->measure() / dinfo.cell->measure();

      for (unsigned k = 0; k < fe.n_quadrature_points; ++k)
        {
          const double diff = exact_values[k] - uh[k];
          dinfo.value(0) += penalty * diff * diff * fe.JxW(k);
        }
      dinfo.value(0) = std::sqrt(dinfo.value(0));
    }


    template <int dim>
    void face(MeshWorker::DoFInfo<dim>         &dinfo1,
              MeshWorker::DoFInfo<dim>         &dinfo2,
              MeshWorker::IntegrationInfo<dim> &info1,
              MeshWorker::IntegrationInfo<dim> &info2)
    {
      const FEValuesBase<dim>   &fe  = info1.fe_values();
      const std::vector<double> &uh1 = info1.values[0][0];
      const std::vector<double> &uh2 = info2.values[0][0];

      const unsigned int degree = fe.get_fe().tensor_degree();
      const double       penalty1 =
        degree * (degree + 1) * dinfo1.face->measure() / dinfo1.cell->measure();
      const double penalty2 =
        degree * (degree + 1) * dinfo2.face->measure() / dinfo2.cell->measure();
      const double penalty = penalty1 + penalty2;

      for (unsigned k = 0; k < fe.n_quadrature_points; ++k)
        {
          const double diff = uh1[k] - uh2[k];
          dinfo1.value(0) += (penalty * diff * diff) * fe.JxW(k);
        }
      dinfo1.value(0) = std::sqrt(dinfo1.value(0));
      dinfo2.value(0) = dinfo1.value(0);
    }
  } // namespace ErrorIntegrator


  // @sect3{The main class}

  // This class does the main job, like in previous examples. For a
  // description of the functions declared here, please refer to the
  // implementation below.
  template <int dim>
  class InteriorPenaltyProblem
  {
  public:
    using CellInfo = MeshWorker::IntegrationInfo<dim>;

    InteriorPenaltyProblem();

    void run(unsigned int n_steps);

  private:
    void   setup_system();
    void   assemble_matrix();
    void   assemble_mg_matrix();
    void   assemble_right_hand_side();
    void   error();
    double estimate();
    void   solve();
    void   output_results(const unsigned int cycle) const;

    // The member objects related to the discretization are here.
    Triangulation<dim>   triangulation;
    const MappingQ1<dim> mapping;
    const FE_DGQ<2>      fe;
    DoFHandler<dim>      dof_handler;

    // Then, we have the matrices and vectors related to the global discrete
    // system.
    SparsityPattern      sparsity;
    SparseMatrix<double> matrix;
    Vector<double>       solution;
    Vector<double>       right_hand_side;
    BlockVector<double>  estimates;

    // Finally, we have a group of sparsity patterns and sparse matrices
    // related to the multilevel preconditioner.  First, we have a level
    // matrix and its sparsity pattern.
    MGLevelObject<SparsityPattern>      mg_sparsity;
    MGLevelObject<SparseMatrix<double>> mg_matrix;

    // When we perform multigrid with local smoothing on locally refined
    // meshes, additional matrices are required; see Kanschat (2004). Here is
    // the sparsity pattern for these edge matrices. We only need one, because
    // the pattern of the up matrix is the transpose of that of the down
    // matrix. Actually, we do not care too much about these details, since
    // the MeshWorker is filling these matrices.
    MGLevelObject<SparsityPattern> mg_sparsity_dg_interface;
    // The flux matrix at the refinement edge, coupling fine level degrees of
    // freedom to coarse level.
    MGLevelObject<SparseMatrix<double>> mg_matrix_dg_down;
    // The transpose of the flux matrix at the refinement edge, coupling
    // coarse level degrees of freedom to fine level.
    MGLevelObject<SparseMatrix<double>> mg_matrix_dg_up;
  };


  // The constructor simply sets up the coarse grid and the DoFHandler.
  template <int dim>
  InteriorPenaltyProblem<dim>::InteriorPenaltyProblem()
    : triangulation(Triangulation<dim>::limit_level_difference_at_vertices)
    , mapping()
    , fe(3)
    , dof_handler(triangulation)
    , estimates(1)
  {
    GridGenerator::hyper_cube_slit(triangulation, -1, 1);
  }


  // In this function, we set up the dimension of the linear system and the
  // sparsity patterns for the global matrix as well as the level matrices.
  template <int dim>
  void InteriorPenaltyProblem<dim>::setup_system()
  {
    // First, we use the finite element to distribute degrees of freedom over
    // the mesh and number them.
    dof_handler.distribute_dofs(fe);
    dof_handler.distribute_mg_dofs();
    unsigned int n_dofs = dof_handler.n_dofs();
    // Then, we already know the size of the vectors representing finite
    // element functions.
    solution.reinit(n_dofs);
    right_hand_side.reinit(n_dofs);

    // Next, we set up the sparsity pattern for the global matrix. Since we do
    // not know the row sizes in advance, we first fill a temporary
    // DynamicSparsityPattern object and copy it to the regular
    // SparsityPattern once it is complete.
    DynamicSparsityPattern dsp(n_dofs);
    DoFTools::make_flux_sparsity_pattern(dof_handler, dsp);
    sparsity.copy_from(dsp);
    matrix.reinit(sparsity);

    const unsigned int n_levels = triangulation.n_levels();
    // The global system is set up, now we attend to the level matrices. We
    // resize all matrix objects to hold one matrix per level.
    mg_matrix.resize(0, n_levels - 1);
    mg_matrix.clear_elements();
    mg_matrix_dg_up.resize(0, n_levels - 1);
    mg_matrix_dg_up.clear_elements();
    mg_matrix_dg_down.resize(0, n_levels - 1);
    mg_matrix_dg_down.clear_elements();
    // It is important to update the sparsity patterns after <tt>clear()</tt>
    // was called for the level matrices, since the matrices lock the sparsity
    // pattern through the SmartPointer and Subscriptor mechanism.
    mg_sparsity.resize(0, n_levels - 1);
    mg_sparsity_dg_interface.resize(0, n_levels - 1);

    // Now all objects are prepared to hold one sparsity pattern or matrix per
    // level. What's left is setting up the sparsity patterns on each level.
    for (unsigned int level = mg_sparsity.min_level();
         level <= mg_sparsity.max_level();
         ++level)
      {
        // These are roughly the same lines as above for the global matrix,
        // now for each level.
        DynamicSparsityPattern dsp(dof_handler.n_dofs(level));
        MGTools::make_flux_sparsity_pattern(dof_handler, dsp, level);
        mg_sparsity[level].copy_from(dsp);
        mg_matrix[level].reinit(mg_sparsity[level]);

        // Additionally, we need to initialize the transfer matrices at the
        // refinement edge between levels. They are stored at the index
        // referring to the finer of the two indices, thus there is no such
        // object on level 0.
        if (level > 0)
          {
            DynamicSparsityPattern dsp;
            dsp.reinit(dof_handler.n_dofs(level - 1),
                       dof_handler.n_dofs(level));
            MGTools::make_flux_sparsity_pattern_edge(dof_handler, dsp, level);
            mg_sparsity_dg_interface[level].copy_from(dsp);
            mg_matrix_dg_up[level].reinit(mg_sparsity_dg_interface[level]);
            mg_matrix_dg_down[level].reinit(mg_sparsity_dg_interface[level]);
          }
      }
  }


  // In this function, we assemble the global system matrix, where by global
  // we indicate that this is the matrix of the discrete system we solve and
  // it is covering the whole mesh.
  template <int dim>
  void InteriorPenaltyProblem<dim>::assemble_matrix()
  {
    // First, we need t set up the object providing the values we
    // integrate. This object contains all FEValues and FEFaceValues objects
    // needed and also maintains them automatically such that they always
    // point to the current cell. To this end, we need to tell it first, where
    // and what to compute. Since we are not doing anything fancy, we can rely
    // on their standard choice for quadrature rules.
    //
    // Since their default update flags are minimal, we add what we need
    // additionally, namely the values and gradients of shape functions on all
    // objects (cells, boundary and interior faces). Afterwards, we are ready
    // to initialize the container, which will create all necessary
    // FEValuesBase objects for integration.
    MeshWorker::IntegrationInfoBox<dim> info_box;
    UpdateFlags update_flags = update_values | update_gradients;
    info_box.add_update_flags_all(update_flags);
    info_box.initialize(fe, mapping);

    // This is the object into which we integrate local data. It is filled by
    // the local integration routines in `MatrixIntegrator` and then used by the
    // assembler to distribute the information into the global matrix.
    MeshWorker::DoFInfo<dim> dof_info(dof_handler);

    // Furthermore, we need an object that assembles the local matrix into the
    // global matrix. These assembler objects have all the knowledge
    // of the structures of the target object, in this case a
    // SparseMatrix, possible constraints and the mesh structure.
    MeshWorker::Assembler::MatrixSimple<SparseMatrix<double>> assembler;
    assembler.initialize(matrix);

    // Now, we throw everything into a MeshWorker::loop<dim, dim>(), which here
    // traverses all active cells of the mesh, computes cell and face matrices
    // and assembles them into the global matrix. We use the variable
    // <tt>dof_handler</tt> here in order to use the global numbering of
    // degrees of freedom.
    MeshWorker::loop<dim, dim>(dof_handler.begin_active(),
                               dof_handler.end(),
                               dof_info,
                               info_box,
                               &MatrixIntegrator::cell<dim>,
                               &MatrixIntegrator::boundary<dim>,
                               &MatrixIntegrator::face<dim>,
                               assembler);
  }


  // Now, we do the same for the level matrices. Not too surprisingly, this
  // function looks like a twin of the previous one. Indeed, there are only
  // two minor differences.
  template <int dim>
  void InteriorPenaltyProblem<dim>::assemble_mg_matrix()
  {
    MeshWorker::IntegrationInfoBox<dim> info_box;
    UpdateFlags update_flags = update_values | update_gradients;
    info_box.add_update_flags_all(update_flags);
    info_box.initialize(fe, mapping);

    MeshWorker::DoFInfo<dim> dof_info(dof_handler);

    // Obviously, the assembler needs to be replaced by one filling level
    // matrices. Note that it automatically fills the edge matrices as well.
    MeshWorker::Assembler::MGMatrixSimple<SparseMatrix<double>> assembler;
    assembler.initialize(mg_matrix);
    assembler.initialize_fluxes(mg_matrix_dg_up, mg_matrix_dg_down);

    // Here is the other difference to the previous function: we run
    // over all cells, not only the active ones. And we use functions
    // ending on <code>_mg</code> since we need the degrees of freedom
    // on each level, not the global numbering.
    MeshWorker::loop<dim, dim>(dof_handler.begin_mg(),
                               dof_handler.end_mg(),
                               dof_info,
                               info_box,
                               &MatrixIntegrator::cell<dim>,
                               &MatrixIntegrator::boundary<dim>,
                               &MatrixIntegrator::face<dim>,
                               assembler);
  }


  // Here we have another clone of the assemble function. The difference to
  // assembling the system matrix consists in that we assemble a vector here.
  template <int dim>
  void InteriorPenaltyProblem<dim>::assemble_right_hand_side()
  {
    MeshWorker::IntegrationInfoBox<dim> info_box;
    UpdateFlags                         update_flags =
      update_quadrature_points | update_values | update_gradients;
    info_box.add_update_flags_all(update_flags);
    info_box.initialize(fe, mapping);

    MeshWorker::DoFInfo<dim> dof_info(dof_handler);

    // Since this assembler allows us to fill several vectors, the interface is
    // a little more complicated as above. The pointers to the vectors have to
    // be stored in an AnyData object. While this seems to cause two extra
    // lines of code here, it actually comes handy in more complex
    // applications.
    MeshWorker::Assembler::ResidualSimple<Vector<double>> assembler;
    AnyData                                               data;
    data.add<Vector<double> *>(&right_hand_side, "RHS");
    assembler.initialize(data);

    MeshWorker::loop<dim, dim>(dof_handler.begin_active(),
                               dof_handler.end(),
                               dof_info,
                               info_box,
                               &RHSIntegrator::cell<dim>,
                               &RHSIntegrator::boundary<dim>,
                               &RHSIntegrator::face<dim>,
                               assembler);

    right_hand_side *= -1.;
  }


  // Now that we have coded all functions building the discrete linear system,
  // it is about time that we actually solve it.
  template <int dim>
  void InteriorPenaltyProblem<dim>::solve()
  {
    // The solver of choice is conjugate gradient.
    SolverControl            control(1000, 1.e-12, false, true);
    SolverCG<Vector<double>> solver(control);

    // Now we are setting up the components of the multilevel
    // preconditioner. First, we need transfer between grid levels. The object
    // we are using here generates sparse matrices for these transfers.
    MGTransferPrebuilt<Vector<double>> mg_transfer;
    mg_transfer.build(dof_handler);

    // Then, we need an exact solver for the matrix on the coarsest level.
    FullMatrix<double> coarse_matrix;
    coarse_matrix.copy_from(mg_matrix[0]);
    MGCoarseGridHouseholder<double, Vector<double>> mg_coarse;
    mg_coarse.initialize(coarse_matrix);

    // While transfer and coarse grid solver are pretty much generic, more
    // flexibility is offered for the smoother. First, we choose Gauss-Seidel
    // as our smoothing method.
    GrowingVectorMemory<Vector<double>> mem;
    using RELAXATION = PreconditionSOR<SparseMatrix<double>>;
    mg::SmootherRelaxation<RELAXATION, Vector<double>> mg_smoother;
    RELAXATION::AdditionalData                         smoother_data(1.);
    mg_smoother.initialize(mg_matrix, smoother_data);

    // Do two smoothing steps on each level.
    mg_smoother.set_steps(2);
    // Since the SOR method is not symmetric, but we use conjugate gradient
    // iteration below, here is a trick to make the multilevel preconditioner
    // a symmetric operator even for nonsymmetric smoothers.
    mg_smoother.set_symmetric(true);
    // The smoother class optionally implements the variable V-cycle, which we
    // do not want here.
    mg_smoother.set_variable(false);

    // Finally, we must wrap our matrices in an object having the required
    // multiplication functions.
    mg::Matrix<Vector<double>> mgmatrix(mg_matrix);
    mg::Matrix<Vector<double>> mgdown(mg_matrix_dg_down);
    mg::Matrix<Vector<double>> mgup(mg_matrix_dg_up);

    // Now, we are ready to set up the V-cycle operator and the multilevel
    // preconditioner.
    Multigrid<Vector<double>> mg(
      mgmatrix, mg_coarse, mg_transfer, mg_smoother, mg_smoother);
    // Let us not forget the edge matrices needed because of the adaptive
    // refinement.
    mg.set_edge_flux_matrices(mgdown, mgup);

    // After all preparations, wrap the Multigrid object into another object,
    // which can be used as a regular preconditioner,
    PreconditionMG<dim, Vector<double>, MGTransferPrebuilt<Vector<double>>>
      preconditioner(dof_handler, mg, mg_transfer);
    // and use it to solve the system.
    solver.solve(matrix, solution, right_hand_side, preconditioner);
  }


  // Another clone of the assemble function. The big difference to the
  // previous ones is here that we also have an input vector.
  template <int dim>
  double InteriorPenaltyProblem<dim>::estimate()
  {
    // The results of the estimator are stored in a vector with one entry per
    // cell. Since cells in deal.II are not numbered, we have to create our
    // own numbering in order to use this vector. For the assembler used below
    // the information in which component of a vector the result is stored is
    // transmitted by the user_index variable for each cell. We need to set this
    // numbering up here.
    //
    // On the other hand, somebody might have used the user indices
    // already. So, let's be good citizens and save them before tampering with
    // them.
    std::vector<unsigned int> old_user_indices;
    triangulation.save_user_indices(old_user_indices);

    estimates.block(0).reinit(triangulation.n_active_cells());
    unsigned int i = 0;
    for (const auto &cell : triangulation.active_cell_iterators())
      cell->set_user_index(i++);

    // This starts like before,
    MeshWorker::IntegrationInfoBox<dim> info_box;
    const unsigned int                  n_gauss_points =
      dof_handler.get_fe().tensor_degree() + 1;
    info_box.initialize_gauss_quadrature(n_gauss_points,
                                         n_gauss_points + 1,
                                         n_gauss_points);

    // but now we need to notify the info box of the finite element function we
    // want to evaluate in the quadrature points. First, we create an AnyData
    // object with this vector, which is the solution we just computed.
    AnyData solution_data;
    solution_data.add<const Vector<double> *>(&solution, "solution");

    // Then, we tell the Meshworker::VectorSelector for cells, that we need
    // the second derivatives of this solution (to compute the
    // Laplacian). Therefore, the Boolean arguments selecting function values
    // and first derivatives a false, only the last one selecting second
    // derivatives is true.
    info_box.cell_selector.add("solution", false, false, true);
    // On interior and boundary faces, we need the function values and the
    // first derivatives, but not second derivatives.
    info_box.boundary_selector.add("solution", true, true, false);
    info_box.face_selector.add("solution", true, true, false);

    // And we continue as before, with the exception that the default update
    // flags are already adjusted to the values and derivatives we requested
    // above.
    info_box.add_update_flags_boundary(update_quadrature_points);
    info_box.initialize(fe, mapping, solution_data, solution);

    MeshWorker::DoFInfo<dim> dof_info(dof_handler);

    // The assembler stores one number per cell, but else this is the same as
    // in the computation of the right hand side.
    MeshWorker::Assembler::CellsAndFaces<double> assembler;
    AnyData                                      out_data;
    out_data.add<BlockVector<double> *>(&estimates, "cells");
    assembler.initialize(out_data, false);

    MeshWorker::loop<dim, dim>(dof_handler.begin_active(),
                               dof_handler.end(),
                               dof_info,
                               info_box,
                               &Estimator::cell<dim>,
                               &Estimator::boundary<dim>,
                               &Estimator::face<dim>,
                               assembler);

    // Right before we return the result of the error estimate, we restore the
    // old user indices.
    triangulation.load_user_indices(old_user_indices);
    return estimates.block(0).l2_norm();
  }

  // Here we compare our finite element solution with the (known) exact
  // solution and compute the mean quadratic error of the gradient and the
  // function itself. This function is a clone of the estimation function
  // right above.

  // Since we compute the error in the energy and the
  // <i>L<sup>2</sup></i>-norm, respectively, our block vector needs two
  // blocks here.
  template <int dim>
  void InteriorPenaltyProblem<dim>::error()
  {
    BlockVector<double> errors(2);
    errors.block(0).reinit(triangulation.n_active_cells());
    errors.block(1).reinit(triangulation.n_active_cells());

    std::vector<unsigned int> old_user_indices;
    triangulation.save_user_indices(old_user_indices);
    unsigned int i = 0;
    for (const auto &cell : triangulation.active_cell_iterators())
      cell->set_user_index(i++);

    MeshWorker::IntegrationInfoBox<dim> info_box;
    const unsigned int                  n_gauss_points =
      dof_handler.get_fe().tensor_degree() + 1;
    info_box.initialize_gauss_quadrature(n_gauss_points,
                                         n_gauss_points + 1,
                                         n_gauss_points);

    AnyData solution_data;
    solution_data.add<Vector<double> *>(&solution, "solution");

    info_box.cell_selector.add("solution", true, true, false);
    info_box.boundary_selector.add("solution", true, false, false);
    info_box.face_selector.add("solution", true, false, false);

    info_box.add_update_flags_cell(update_quadrature_points);
    info_box.add_update_flags_boundary(update_quadrature_points);
    info_box.initialize(fe, mapping, solution_data, solution);

    MeshWorker::DoFInfo<dim> dof_info(dof_handler);

    MeshWorker::Assembler::CellsAndFaces<double> assembler;
    AnyData                                      out_data;
    out_data.add<BlockVector<double> *>(&errors, "cells");
    assembler.initialize(out_data, false);

    MeshWorker::loop<dim, dim>(dof_handler.begin_active(),
                               dof_handler.end(),
                               dof_info,
                               info_box,
                               &ErrorIntegrator::cell<dim>,
                               &ErrorIntegrator::boundary<dim>,
                               &ErrorIntegrator::face<dim>,
                               assembler);
    triangulation.load_user_indices(old_user_indices);

    deallog << "energy-error: " << errors.block(0).l2_norm() << std::endl;
    deallog << "L2-error:     " << errors.block(1).l2_norm() << std::endl;
  }


  // Create graphical output. We produce the filename by collating the
  // name from its various components, including the refinement cycle
  // that we output with two digits.
  template <int dim>
  void
  InteriorPenaltyProblem<dim>::output_results(const unsigned int cycle) const
  {
    const std::string filename =
      "sol-" + Utilities::int_to_string(cycle, 2) + ".gnuplot";

    deallog << "Writing solution to <" << filename << ">..." << std::endl
            << std::endl;
    std::ofstream gnuplot_output(filename);

    DataOut<dim> data_out;
    data_out.attach_dof_handler(dof_handler);
    data_out.add_data_vector(solution, "u");
    data_out.add_data_vector(estimates.block(0), "est");

    data_out.build_patches();

    data_out.write_gnuplot(gnuplot_output);
  }

  // And finally the adaptive loop, more or less like in previous examples.
  template <int dim>
  void InteriorPenaltyProblem<dim>::run(unsigned int n_steps)
  {
    deallog << "Element: " << fe.get_name() << std::endl;
    for (unsigned int s = 0; s < n_steps; ++s)
      {
        deallog << "Step " << s << std::endl;
        if (estimates.block(0).size() == 0)
          triangulation.refine_global(1);
        else
          {
            GridRefinement::refine_and_coarsen_fixed_fraction(
              triangulation, estimates.block(0), 0.5, 0.0);
            triangulation.execute_coarsening_and_refinement();
          }

        deallog << "Triangulation " << triangulation.n_active_cells()
                << " cells, " << triangulation.n_levels() << " levels"
                << std::endl;

        setup_system();
        deallog << "DoFHandler " << dof_handler.n_dofs() << " dofs, level dofs";
        for (unsigned int l = 0; l < triangulation.n_levels(); ++l)
          deallog << ' ' << dof_handler.n_dofs(l);
        deallog << std::endl;

        deallog << "Assemble matrix" << std::endl;
        assemble_matrix();
        deallog << "Assemble multilevel matrix" << std::endl;
        assemble_mg_matrix();
        deallog << "Assemble right hand side" << std::endl;
        assemble_right_hand_side();
        deallog << "Solve" << std::endl;
        solve();
        error();
        deallog << "Estimate " << estimate() << std::endl;
        output_results(s);
      }
  }
} // namespace Step39



int main()
{
  try
    {
      using namespace dealii;
      using namespace Step39;

      deallog.depth_console(2);
      std::ofstream logfile("deallog");
      deallog.attach(logfile);

      InteriorPenaltyProblem<2> test1;
      test1.run(12);
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
