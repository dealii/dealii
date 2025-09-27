/* ------------------------------------------------------------------------
 *
 * SPDX-License-Identifier: LGPL-2.1-or-later
 * Copyright (C) 2021 - 2024 by the deal.II authors
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
 * Authors: Andrea Bonito and Diane Guignard, 2021.
 */

// @sect3{Include files}

// All the include files have already been discussed in previous tutorials.
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_system.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>

#include <deal.II/numerics/data_out.h>

#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
// The following three header files are for the solvers.
// The linear system is solved using a direct method
// while the conjugate gradient method is used to solve
// the local problems for the lifting terms.
#include <deal.II/lac/sparse_direct.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>

#include <fstream>
#include <iostream>


namespace Step82
{
  using namespace dealii;

  // @sect3{The <code>BiLaplacianLDGLift</code> class template}

  // The main class of this program is similar to that of step-3
  // or step-20, as well as many other tutorial programs. The key
  // function here is <code>compute_discrete_hessians()</code> which, as its
  // name suggests, computes the discrete Hessians needed for the assembly of
  // the matrix $A$.
  template <int dim>
  class BiLaplacianLDGLift
  {
  public:
    BiLaplacianLDGLift(const unsigned int n_refinements,
                       const unsigned int fe_degree,
                       const double       penalty_jump_grad,
                       const double       penalty_jump_val);

    void run();

  private:
    void make_grid();
    void setup_system();
    void assemble_system();
    void assemble_matrix();
    void assemble_rhs();

    void solve();

    void compute_errors();
    void output_results() const;

    // As indicated by its name, the function
    // <code>assemble_local_matrix()</code> is used for the assembly of the
    // (local) @ref GlossMassMatrix "mass matrix" used to compute the two lifting terms (see the matrix
    // $\boldsymbol{M}_c$ introduced in the introduction when describing the
    // computation of $b_e$). The function
    // <code>compute_discrete_hessians()</code> computes the required discrete
    // Hessians: the discrete Hessians of the basis functions with support on
    // the current <code>cell</code> (stored in the output variable
    // <code>discrete_hessians</code>) and the basis functions with support on a
    // neighbor of the current <code>cell</code> (stored in the output variable
    // <code>discrete_hessians_neigh</code>). More precisely,
    // <code>discrete_hessians[i][q_point]</code> stores $H_h(\varphi_i)(x_q)$,
    // where $\varphi_i$ is a basis function with support on cell, while
    // <code>discrete_hessians_neigh[face_no][i][q_point]</code> stores
    // $H_h(\varphi_i)(x_q)$, where $\varphi_i$ is a basis function of the
    // neighboring cell adjacent to the face
    // <code>face=cell->face(face_no)</code>.
    void assemble_local_matrix(const FEValues<dim> &fe_values_lift,
                               const unsigned int   n_q_points,
                               FullMatrix<double>  &local_matrix);

    void compute_discrete_hessians(
      const typename DoFHandler<dim>::active_cell_iterator &cell,
      std::vector<std::vector<Tensor<2, dim>>>             &discrete_hessians,
      std::vector<std::vector<std::vector<Tensor<2, dim>>>>
        &discrete_hessians_neigh);

    Triangulation<dim> triangulation;

    const unsigned int n_refinements;

    const FE_DGQ<dim> fe;
    DoFHandler<dim>   dof_handler;

    // We also need a variable that describes the finite element space
    // $[\mathbb{V}_h]^{d\times d}$ used for the two lifting
    // operators. The other member variables below are as in most of the other
    // tutorial programs.
    const FESystem<dim> fe_lift;

    SparsityPattern      sparsity_pattern;
    SparseMatrix<double> matrix;
    Vector<double>       rhs;
    Vector<double>       solution;

    // Finally, the last two variables correspond to the penalty coefficients
    // $\gamma_1$ and $\gamma_0$ for the jump of $\nabla_hu_h$ and $u_h$,
    // respectively.
    const double penalty_jump_grad;
    const double penalty_jump_val;
  };



  // @sect3{Equation data}

  // This class implement the right-hand side $f=\Delta^2 u$ corresponding to
  // the manufactured solution $u$ defined in the introduction.
  template <int dim>
  class RightHandSide : public Function<dim>
  {
  public:
    RightHandSide()
      : Function<dim>()
    {}

    virtual double value(const Point<dim>  &p,
                         const unsigned int component = 0) const override;
  };



  template <int dim>
  double RightHandSide<dim>::value(const Point<dim> &p,
                                   const unsigned int /*component*/) const
  {
    double return_value = 0.0;

    if (dim == 2)
      {
        return_value = 24.0 * Utilities::fixed_power<2>(p[1] * (1.0 - p[1])) +
                       +24.0 * Utilities::fixed_power<2>(p[0] * (1.0 - p[0])) +
                       2.0 * (2.0 - 12.0 * p[0] + 12.0 * p[0] * p[0]) *
                         (2.0 - 12.0 * p[1] + 12.0 * p[1] * p[1]);
      }
    else if (dim == 3)
      {
        return_value = 24.0 * Utilities::fixed_power<2>(p[1] * (1.0 - p[1]) *
                                                        p[2] * (1.0 - p[2])) +
                       24.0 * Utilities::fixed_power<2>(p[0] * (1.0 - p[0]) *
                                                        p[2] * (1.0 - p[2])) +
                       24.0 * Utilities::fixed_power<2>(p[0] * (1.0 - p[0]) *
                                                        p[1] * (1.0 - p[1])) +
                       2.0 * (2.0 - 12.0 * p[0] + 12.0 * p[0] * p[0]) *
                         (2.0 - 12.0 * p[1] + 12.0 * p[1] * p[1]) *
                         Utilities::fixed_power<2>(p[2] * (1.0 - p[2])) +
                       2.0 * (2.0 - 12.0 * p[0] + 12.0 * p[0] * p[0]) *
                         (2.0 - 12.0 * p[2] + 12.0 * p[2] * p[2]) *
                         Utilities::fixed_power<2>(p[1] * (1.0 - p[1])) +
                       2.0 * (2.0 - 12.0 * p[1] + 12.0 * p[1] * p[1]) *
                         (2.0 - 12.0 * p[2] + 12.0 * p[2] * p[2]) *
                         Utilities::fixed_power<2>(p[0] * (1.0 - p[0]));
      }
    else
      DEAL_II_NOT_IMPLEMENTED();

    return return_value;
  }



  // This class implement the manufactured (exact) solution $u$. To compute the
  // errors, we need the value of $u$ as well as its gradient and its Hessian.
  template <int dim>
  class ExactSolution : public Function<dim>
  {
  public:
    ExactSolution()
      : Function<dim>()
    {}

    virtual double value(const Point<dim>  &p,
                         const unsigned int component = 0) const override;

    virtual Tensor<1, dim>
    gradient(const Point<dim>  &p,
             const unsigned int component = 0) const override;

    virtual SymmetricTensor<2, dim>
    hessian(const Point<dim>  &p,
            const unsigned int component = 0) const override;
  };



  template <int dim>
  double ExactSolution<dim>::value(const Point<dim> &p,
                                   const unsigned int /*component*/) const
  {
    double return_value = 0.0;

    if (dim == 2)
      {
        return_value =
          Utilities::fixed_power<2>(p[0] * (1.0 - p[0]) * p[1] * (1.0 - p[1]));
      }
    else if (dim == 3)
      {
        return_value = Utilities::fixed_power<2>(
          p[0] * (1.0 - p[0]) * p[1] * (1.0 - p[1]) * p[2] * (1.0 - p[2]));
      }
    else
      DEAL_II_NOT_IMPLEMENTED();

    return return_value;
  }



  template <int dim>
  Tensor<1, dim>
  ExactSolution<dim>::gradient(const Point<dim> &p,
                               const unsigned int /*component*/) const
  {
    Tensor<1, dim> return_gradient;

    if (dim == 2)
      {
        return_gradient[0] =
          (2.0 * p[0] - 6.0 * Utilities::fixed_power<2>(p[0]) +
           4.0 * Utilities::fixed_power<3>(p[0])) *
          Utilities::fixed_power<2>(p[1] * (1.0 - p[1]));
        return_gradient[1] =
          (2.0 * p[1] - 6.0 * Utilities::fixed_power<2>(p[1]) +
           4.0 * Utilities::fixed_power<3>(p[1])) *
          Utilities::fixed_power<2>(p[0] * (1.0 - p[0]));
      }
    else if (dim == 3)
      {
        return_gradient[0] =
          (2.0 * p[0] - 6.0 * Utilities::fixed_power<2>(p[0]) +
           4.0 * Utilities::fixed_power<3>(p[0])) *
          Utilities::fixed_power<2>(p[1] * (1.0 - p[1]) * p[2] * (1.0 - p[2]));
        return_gradient[1] =
          (2.0 * p[1] - 6.0 * Utilities::fixed_power<2>(p[1]) +
           4.0 * Utilities::fixed_power<3>(p[1])) *
          Utilities::fixed_power<2>(p[0] * (1.0 - p[0]) * p[2] * (1.0 - p[2]));
        return_gradient[2] =
          (2.0 * p[2] - 6.0 * Utilities::fixed_power<2>(p[2]) +
           4.0 * Utilities::fixed_power<3>(p[2])) *
          Utilities::fixed_power<2>(p[0] * (1.0 - p[0]) * p[1] * (1.0 - p[1]));
      }
    else
      DEAL_II_NOT_IMPLEMENTED();

    return return_gradient;
  }



  template <int dim>
  SymmetricTensor<2, dim>
  ExactSolution<dim>::hessian(const Point<dim> &p,
                              const unsigned int /*component*/) const
  {
    SymmetricTensor<2, dim> return_hessian;

    if (dim == 2)
      {
        return_hessian[0][0] = (2.0 - 12.0 * p[0] + 12.0 * p[0] * p[0]) *
                               Utilities::fixed_power<2>(p[1] * (1.0 - p[1]));
        return_hessian[0][1] =
          (2.0 * p[0] - 6.0 * Utilities::fixed_power<2>(p[0]) +
           4.0 * Utilities::fixed_power<3>(p[0])) *
          (2.0 * p[1] - 6.0 * Utilities::fixed_power<2>(p[1]) +
           4.0 * Utilities::fixed_power<3>(p[1]));
        return_hessian[1][1] = (2.0 - 12.0 * p[1] + 12.0 * p[1] * p[1]) *
                               Utilities::fixed_power<2>(p[0] * (1.0 - p[0]));
      }
    else if (dim == 3)
      {
        return_hessian[0][0] =
          (2.0 - 12.0 * p[0] + 12.0 * p[0] * p[0]) *
          Utilities::fixed_power<2>(p[1] * (1.0 - p[1]) * p[2] * (1.0 - p[2]));
        return_hessian[0][1] =
          (2.0 * p[0] - 6.0 * Utilities::fixed_power<2>(p[0]) +
           4.0 * Utilities::fixed_power<3>(p[0])) *
          (2.0 * p[1] - 6.0 * Utilities::fixed_power<2>(p[1]) +
           4.0 * Utilities::fixed_power<3>(p[1])) *
          Utilities::fixed_power<2>(p[2] * (1.0 - p[2]));
        return_hessian[0][2] =
          (2.0 * p[0] - 6.0 * Utilities::fixed_power<2>(p[0]) +
           4.0 * Utilities::fixed_power<3>(p[0])) *
          (2.0 * p[2] - 6.0 * Utilities::fixed_power<2>(p[2]) +
           4.0 * Utilities::fixed_power<3>(p[2])) *
          Utilities::fixed_power<2>(p[1] * (1.0 - p[1]));
        return_hessian[1][1] =
          (2.0 - 12.0 * p[1] + 12.0 * p[1] * p[1]) *
          Utilities::fixed_power<2>(p[0] * (1.0 - p[0]) * p[2] * (1.0 - p[2]));
        return_hessian[1][2] =
          (2.0 * p[1] - 6.0 * Utilities::fixed_power<2>(p[1]) +
           4.0 * Utilities::fixed_power<3>(p[1])) *
          (2.0 * p[2] - 6.0 * Utilities::fixed_power<2>(p[2]) +
           4.0 * Utilities::fixed_power<3>(p[2])) *
          Utilities::fixed_power<2>(p[0] * (1.0 - p[0]));
        return_hessian[2][2] =
          (2.0 - 12.0 * p[2] + 12.0 * p[2] * p[2]) *
          Utilities::fixed_power<2>(p[0] * (1.0 - p[0]) * p[1] * (1.0 - p[1]));
      }
    else
      DEAL_II_NOT_IMPLEMENTED();

    return return_hessian;
  }



  // @sect3{Implementation of the <code>BiLaplacianLDGLift</code> class}

  // @sect4{BiLaplacianLDGLift::BiLaplacianLDGLift}

  // In the constructor, we set the polynomial degree of the two finite element
  // spaces, we associate the corresponding DoF handlers to the triangulation,
  // and we set the two penalty coefficients.
  template <int dim>
  BiLaplacianLDGLift<dim>::BiLaplacianLDGLift(const unsigned int n_refinements,
                                              const unsigned int fe_degree,
                                              const double penalty_jump_grad,
                                              const double penalty_jump_val)
    : n_refinements(n_refinements)
    , fe(fe_degree)
    , dof_handler(triangulation)
    , fe_lift(FE_DGQ<dim>(fe_degree), dim * dim)
    , penalty_jump_grad(penalty_jump_grad)
    , penalty_jump_val(penalty_jump_val)
  {}



  // @sect4{BiLaplacianLDGLift::make_grid}

  // To build a mesh for $\Omega=(0,1)^d$, we simply call the function
  // <code>GridGenerator::hyper_cube</code> and then refine it using
  // <code>refine_global</code>. The number of refinements is hard-coded
  // here.
  template <int dim>
  void BiLaplacianLDGLift<dim>::make_grid()
  {
    std::cout << "Building the mesh............." << std::endl;

    GridGenerator::hyper_cube(triangulation, 0.0, 1.0);

    triangulation.refine_global(n_refinements);

    std::cout << "Number of active cells: " << triangulation.n_active_cells()
              << std::endl;
  }



  // @sect4{BiLaplacianLDGLift::setup_system}

  // In the following function, we set up the degrees of freedom, the sparsity
  // pattern, the size of the matrix $A$, and the size of the solution and
  // right-hand side vectors
  // $\boldsymbol{U}$ and $\boldsymbol{F}$. For the sparsity pattern, we cannot
  // directly use the function <code>DoFTools::make_flux_sparsity_pattern</code>
  // (as we would do for instance for the SIPG method) because we need to take
  // into account the interactions of a neighboring cell with another
  // neighboring cell as described in the introduction. The extended sparsity
  // pattern is built by iterating over all the active cells. For the current
  // cell, we collect all its degrees of freedom as well as the degrees of
  // freedom of all its neighboring cells, and then couple everything with
  // everything.
  template <int dim>
  void BiLaplacianLDGLift<dim>::setup_system()
  {
    dof_handler.distribute_dofs(fe);

    std::cout << "Number of degrees of freedom: " << dof_handler.n_dofs()
              << std::endl;

    DynamicSparsityPattern dsp(dof_handler.n_dofs(), dof_handler.n_dofs());

    const auto dofs_per_cell = fe.dofs_per_cell;

    for (const auto &cell : dof_handler.active_cell_iterators())
      {
        std::vector<types::global_dof_index> dofs(dofs_per_cell);
        cell->get_dof_indices(dofs);

        for (unsigned int f = 0; f < cell->n_faces(); ++f)
          if (!cell->face(f)->at_boundary())
            {
              const auto neighbor_cell = cell->neighbor(f);

              std::vector<types::global_dof_index> tmp(dofs_per_cell);
              neighbor_cell->get_dof_indices(tmp);

              dofs.insert(std::end(dofs), std::begin(tmp), std::end(tmp));
            }

        for (const auto i : dofs)
          for (const auto j : dofs)
            {
              dsp.add(i, j);
              dsp.add(j, i);
            }
      }

    sparsity_pattern.copy_from(dsp);


    matrix.reinit(sparsity_pattern);
    rhs.reinit(dof_handler.n_dofs());

    solution.reinit(dof_handler.n_dofs());

    // At the end of the function, we output this sparsity pattern as
    // a scalable vector graphic. You can visualize it by loading this
    // file in most web browsers:
    std::ofstream out("sparsity-pattern.svg");
    sparsity_pattern.print_svg(out);
  }



  // @sect4{BiLaplacianLDGLift::assemble_system}

  // This function simply calls the two functions responsible
  // for the assembly of the matrix and the right-hand side.
  template <int dim>
  void BiLaplacianLDGLift<dim>::assemble_system()
  {
    std::cout << "Assembling the system............." << std::endl;

    assemble_matrix();
    assemble_rhs();

    std::cout << "Done. " << std::endl;
  }



  // @sect4{BiLaplacianLDGLift::assemble_matrix}

  // This function assembles the matrix $A$ whose entries are defined
  // by $A_{ij}=A_h(\varphi_j,\varphi_i)$ which involves the product of
  // discrete Hessians and the penalty terms.
  template <int dim>
  void BiLaplacianLDGLift<dim>::assemble_matrix()
  {
    matrix = 0;

    const QGauss<dim>     quad(fe.degree + 1);
    const QGauss<dim - 1> quad_face(fe.degree + 1);

    const unsigned int n_q_points      = quad.size();
    const unsigned int n_q_points_face = quad_face.size();

    FEValues<dim> fe_values(fe, quad, update_hessians | update_JxW_values);

    FEFaceValues<dim> fe_face(
      fe, quad_face, update_values | update_gradients | update_normal_vectors);

    FEFaceValues<dim> fe_face_neighbor(
      fe, quad_face, update_values | update_gradients | update_normal_vectors);

    const unsigned int n_dofs = fe_values.dofs_per_cell;

    std::vector<types::global_dof_index> local_dof_indices(n_dofs);
    std::vector<types::global_dof_index> local_dof_indices_neighbor(n_dofs);
    std::vector<types::global_dof_index> local_dof_indices_neighbor_2(n_dofs);

    // As indicated in the introduction, the following matrices are used for
    // the contributions of the products of the discrete Hessians.
    FullMatrix<double> stiffness_matrix_cc(n_dofs,
                                           n_dofs); // interactions cell / cell
    FullMatrix<double> stiffness_matrix_cn(
      n_dofs, n_dofs); // interactions cell / neighbor
    FullMatrix<double> stiffness_matrix_nc(
      n_dofs, n_dofs); // interactions neighbor / cell
    FullMatrix<double> stiffness_matrix_nn(
      n_dofs, n_dofs); // interactions neighbor / neighbor
    FullMatrix<double> stiffness_matrix_n1n2(
      n_dofs, n_dofs); // interactions neighbor1 / neighbor2
    FullMatrix<double> stiffness_matrix_n2n1(
      n_dofs, n_dofs); // interactions neighbor2 / neighbor1

    // The following matrices are used for the contributions of the two
    // penalty terms:
    FullMatrix<double> ip_matrix_cc(n_dofs, n_dofs); // interactions cell / cell
    FullMatrix<double> ip_matrix_cn(n_dofs,
                                    n_dofs); // interactions cell / neighbor
    FullMatrix<double> ip_matrix_nc(n_dofs,
                                    n_dofs); // interactions neighbor / cell
    FullMatrix<double> ip_matrix_nn(n_dofs,
                                    n_dofs); // interactions neighbor / neighbor

    std::vector<std::vector<Tensor<2, dim>>> discrete_hessians(
      n_dofs, std::vector<Tensor<2, dim>>(n_q_points));
    std::vector<std::vector<std::vector<Tensor<2, dim>>>>
      discrete_hessians_neigh(GeometryInfo<dim>::faces_per_cell,
                              discrete_hessians);

    for (const auto &cell : dof_handler.active_cell_iterators())
      {
        fe_values.reinit(cell);
        cell->get_dof_indices(local_dof_indices);

        // We now compute all the discrete Hessians that are not vanishing
        // on the current cell, i.e., the discrete Hessian of all the basis
        // functions with support on the current cell or on one of its
        // neighbors.
        compute_discrete_hessians(cell,
                                  discrete_hessians,
                                  discrete_hessians_neigh);

        // First, we compute and add the interactions of the degrees of freedom
        // of the current cell.
        stiffness_matrix_cc = 0;
        for (unsigned int q = 0; q < n_q_points; ++q)
          {
            const double dx = fe_values.JxW(q);

            for (unsigned int i = 0; i < n_dofs; ++i)
              for (unsigned int j = 0; j < n_dofs; ++j)
                {
                  const Tensor<2, dim> &H_i = discrete_hessians[i][q];
                  const Tensor<2, dim> &H_j = discrete_hessians[j][q];

                  stiffness_matrix_cc(i, j) += scalar_product(H_j, H_i) * dx;
                }
          }

        for (unsigned int i = 0; i < n_dofs; ++i)
          for (unsigned int j = 0; j < n_dofs; ++j)
            {
              matrix(local_dof_indices[i], local_dof_indices[j]) +=
                stiffness_matrix_cc(i, j);
            }

        // Next, we compute and add the interactions of the degrees of freedom
        // of the current cell with those of its neighbors. Note that the
        // interactions of the degrees of freedom of a neighbor with those of
        // the same neighbor are included here.
        for (unsigned int face_no = 0; face_no < cell->n_faces(); ++face_no)
          {
            const typename DoFHandler<dim>::face_iterator face =
              cell->face(face_no);

            const bool at_boundary = face->at_boundary();
            if (!at_boundary)
              {
                // There is nothing to be done if boundary face (the liftings of
                // the Dirichlet BCs are accounted for in the assembly of the
                // RHS; in fact, nothing to be done in this program since we
                // prescribe homogeneous BCs).

                const typename DoFHandler<dim>::active_cell_iterator
                  neighbor_cell = cell->neighbor(face_no);
                neighbor_cell->get_dof_indices(local_dof_indices_neighbor);

                stiffness_matrix_cn = 0;
                stiffness_matrix_nc = 0;
                stiffness_matrix_nn = 0;
                for (unsigned int q = 0; q < n_q_points; ++q)
                  {
                    const double dx = fe_values.JxW(q);

                    for (unsigned int i = 0; i < n_dofs; ++i)
                      {
                        for (unsigned int j = 0; j < n_dofs; ++j)
                          {
                            const Tensor<2, dim> &H_i = discrete_hessians[i][q];
                            const Tensor<2, dim> &H_j = discrete_hessians[j][q];

                            const Tensor<2, dim> &H_i_neigh =
                              discrete_hessians_neigh[face_no][i][q];
                            const Tensor<2, dim> &H_j_neigh =
                              discrete_hessians_neigh[face_no][j][q];

                            stiffness_matrix_cn(i, j) +=
                              scalar_product(H_j_neigh, H_i) * dx;
                            stiffness_matrix_nc(i, j) +=
                              scalar_product(H_j, H_i_neigh) * dx;
                            stiffness_matrix_nn(i, j) +=
                              scalar_product(H_j_neigh, H_i_neigh) * dx;
                          }
                      }
                  }

                for (unsigned int i = 0; i < n_dofs; ++i)
                  {
                    for (unsigned int j = 0; j < n_dofs; ++j)
                      {
                        matrix(local_dof_indices[i],
                               local_dof_indices_neighbor[j]) +=
                          stiffness_matrix_cn(i, j);
                        matrix(local_dof_indices_neighbor[i],
                               local_dof_indices[j]) +=
                          stiffness_matrix_nc(i, j);
                        matrix(local_dof_indices_neighbor[i],
                               local_dof_indices_neighbor[j]) +=
                          stiffness_matrix_nn(i, j);
                      }
                  }

              } // boundary check
          }     // for face

        // We now compute and add the interactions of the degrees of freedom of
        // a neighboring cells with those of another neighboring cell (this is
        // where we need the extended sparsity pattern).
        for (unsigned int face_no = 0; face_no < cell->n_faces() - 1; ++face_no)
          {
            const typename DoFHandler<dim>::face_iterator face =
              cell->face(face_no);

            const bool at_boundary = face->at_boundary();
            if (!at_boundary)
              { // nothing to be done if boundary face (the liftings of the
                // Dirichlet BCs are accounted for in the assembly of the RHS;
                // in fact, nothing to be done in this program since we
                // prescribe homogeneous BCs)


                for (unsigned int face_no_2 = face_no + 1;
                     face_no_2 < cell->n_faces();
                     ++face_no_2)
                  {
                    const typename DoFHandler<dim>::face_iterator face_2 =
                      cell->face(face_no_2);

                    const bool at_boundary_2 = face_2->at_boundary();
                    if (!at_boundary_2)
                      {
                        const typename DoFHandler<dim>::active_cell_iterator
                          neighbor_cell = cell->neighbor(face_no);
                        neighbor_cell->get_dof_indices(
                          local_dof_indices_neighbor);
                        const typename DoFHandler<dim>::active_cell_iterator
                          neighbor_cell_2 = cell->neighbor(face_no_2);
                        neighbor_cell_2->get_dof_indices(
                          local_dof_indices_neighbor_2);

                        stiffness_matrix_n1n2 = 0;
                        stiffness_matrix_n2n1 = 0;

                        for (unsigned int q = 0; q < n_q_points; ++q)
                          {
                            const double dx = fe_values.JxW(q);

                            for (unsigned int i = 0; i < n_dofs; ++i)
                              for (unsigned int j = 0; j < n_dofs; ++j)
                                {
                                  const Tensor<2, dim> &H_i_neigh =
                                    discrete_hessians_neigh[face_no][i][q];
                                  const Tensor<2, dim> &H_j_neigh =
                                    discrete_hessians_neigh[face_no][j][q];

                                  const Tensor<2, dim> &H_i_neigh2 =
                                    discrete_hessians_neigh[face_no_2][i][q];
                                  const Tensor<2, dim> &H_j_neigh2 =
                                    discrete_hessians_neigh[face_no_2][j][q];

                                  stiffness_matrix_n1n2(i, j) +=
                                    scalar_product(H_j_neigh2, H_i_neigh) * dx;
                                  stiffness_matrix_n2n1(i, j) +=
                                    scalar_product(H_j_neigh, H_i_neigh2) * dx;
                                }
                          }

                        for (unsigned int i = 0; i < n_dofs; ++i)
                          for (unsigned int j = 0; j < n_dofs; ++j)
                            {
                              matrix(local_dof_indices_neighbor[i],
                                     local_dof_indices_neighbor_2[j]) +=
                                stiffness_matrix_n1n2(i, j);
                              matrix(local_dof_indices_neighbor_2[i],
                                     local_dof_indices_neighbor[j]) +=
                                stiffness_matrix_n2n1(i, j);
                            }
                      } // boundary check face_2
                  }     // for face_2
              }         // boundary check face_1
          }             // for face_1


        // Finally, we compute and add the two penalty terms.
        for (unsigned int face_no = 0; face_no < cell->n_faces(); ++face_no)
          {
            const typename DoFHandler<dim>::face_iterator face =
              cell->face(face_no);

            const double mesh_inv = 1.0 / face->diameter(); // h_e^{-1}
            const double mesh3_inv =
              1.0 / Utilities::fixed_power<3>(face->diameter()); // h_e^{-3}

            fe_face.reinit(cell, face_no);

            ip_matrix_cc = 0; // filled in any case (boundary or interior face)

            const bool at_boundary = face->at_boundary();
            if (at_boundary)
              {
                for (unsigned int q = 0; q < n_q_points_face; ++q)
                  {
                    const double dx = fe_face.JxW(q);

                    for (unsigned int i = 0; i < n_dofs; ++i)
                      for (unsigned int j = 0; j < n_dofs; ++j)
                        {
                          ip_matrix_cc(i, j) += penalty_jump_grad * mesh_inv *
                                                fe_face.shape_grad(j, q) *
                                                fe_face.shape_grad(i, q) * dx;
                          ip_matrix_cc(i, j) += penalty_jump_val * mesh3_inv *
                                                fe_face.shape_value(j, q) *
                                                fe_face.shape_value(i, q) * dx;
                        }
                  }
              }
            else
              { // interior face

                const typename DoFHandler<dim>::active_cell_iterator
                                   neighbor_cell = cell->neighbor(face_no);
                const unsigned int face_no_neighbor =
                  cell->neighbor_of_neighbor(face_no);

                // In the next step, we need to have a global way to compare the
                // cells in order to not calculate the same jump term twice:
                if (neighbor_cell->id() < cell->id())
                  continue; // skip this face (already considered)
                else
                  {
                    fe_face_neighbor.reinit(neighbor_cell, face_no_neighbor);
                    neighbor_cell->get_dof_indices(local_dof_indices_neighbor);

                    ip_matrix_cn = 0;
                    ip_matrix_nc = 0;
                    ip_matrix_nn = 0;

                    for (unsigned int q = 0; q < n_q_points_face; ++q)
                      {
                        const double dx = fe_face.JxW(q);

                        for (unsigned int i = 0; i < n_dofs; ++i)
                          {
                            for (unsigned int j = 0; j < n_dofs; ++j)
                              {
                                ip_matrix_cc(i, j) +=
                                  penalty_jump_grad * mesh_inv *
                                  fe_face.shape_grad(j, q) *
                                  fe_face.shape_grad(i, q) * dx;
                                ip_matrix_cc(i, j) +=
                                  penalty_jump_val * mesh3_inv *
                                  fe_face.shape_value(j, q) *
                                  fe_face.shape_value(i, q) * dx;

                                ip_matrix_cn(i, j) -=
                                  penalty_jump_grad * mesh_inv *
                                  fe_face_neighbor.shape_grad(j, q) *
                                  fe_face.shape_grad(i, q) * dx;
                                ip_matrix_cn(i, j) -=
                                  penalty_jump_val * mesh3_inv *
                                  fe_face_neighbor.shape_value(j, q) *
                                  fe_face.shape_value(i, q) * dx;

                                ip_matrix_nc(i, j) -=
                                  penalty_jump_grad * mesh_inv *
                                  fe_face.shape_grad(j, q) *
                                  fe_face_neighbor.shape_grad(i, q) * dx;
                                ip_matrix_nc(i, j) -=
                                  penalty_jump_val * mesh3_inv *
                                  fe_face.shape_value(j, q) *
                                  fe_face_neighbor.shape_value(i, q) * dx;

                                ip_matrix_nn(i, j) +=
                                  penalty_jump_grad * mesh_inv *
                                  fe_face_neighbor.shape_grad(j, q) *
                                  fe_face_neighbor.shape_grad(i, q) * dx;
                                ip_matrix_nn(i, j) +=
                                  penalty_jump_val * mesh3_inv *
                                  fe_face_neighbor.shape_value(j, q) *
                                  fe_face_neighbor.shape_value(i, q) * dx;
                              }
                          }
                      }
                  } // face not visited yet

              } // boundary check

            for (unsigned int i = 0; i < n_dofs; ++i)
              {
                for (unsigned int j = 0; j < n_dofs; ++j)
                  {
                    matrix(local_dof_indices[i], local_dof_indices[j]) +=
                      ip_matrix_cc(i, j);
                  }
              }

            if (!at_boundary)
              {
                for (unsigned int i = 0; i < n_dofs; ++i)
                  {
                    for (unsigned int j = 0; j < n_dofs; ++j)
                      {
                        matrix(local_dof_indices[i],
                               local_dof_indices_neighbor[j]) +=
                          ip_matrix_cn(i, j);
                        matrix(local_dof_indices_neighbor[i],
                               local_dof_indices[j]) += ip_matrix_nc(i, j);
                        matrix(local_dof_indices_neighbor[i],
                               local_dof_indices_neighbor[j]) +=
                          ip_matrix_nn(i, j);
                      }
                  }
              }

          } // for face
      }     // for cell
  }



  // @sect4{BiLaplacianLDGLift::assemble_rhs}

  // This function assemble the right-hand side of the system. Since we consider
  // homogeneous Dirichlet boundary conditions, imposed weakly in the bilinear
  // form using the Nitsche approach, it only involves the contribution of the
  // forcing term $\int_{\Omega}fv_h$.
  template <int dim>
  void BiLaplacianLDGLift<dim>::assemble_rhs()
  {
    rhs = 0;

    const QGauss<dim> quad(fe.degree + 1);
    FEValues<dim>     fe_values(
      fe, quad, update_values | update_quadrature_points | update_JxW_values);

    const unsigned int n_dofs     = fe_values.dofs_per_cell;
    const unsigned int n_quad_pts = quad.size();

    const RightHandSide<dim> right_hand_side;

    Vector<double>                       local_rhs(n_dofs);
    std::vector<types::global_dof_index> local_dof_indices(n_dofs);

    for (const auto &cell : dof_handler.active_cell_iterators())
      {
        fe_values.reinit(cell);
        cell->get_dof_indices(local_dof_indices);

        local_rhs = 0;
        for (unsigned int q = 0; q < n_quad_pts; ++q)
          {
            const double dx = fe_values.JxW(q);

            for (unsigned int i = 0; i < n_dofs; ++i)
              {
                local_rhs(i) +=
                  right_hand_side.value(fe_values.quadrature_point(q)) *
                  fe_values.shape_value(i, q) * dx;
              }
          }

        for (unsigned int i = 0; i < n_dofs; ++i)
          rhs(local_dof_indices[i]) += local_rhs(i);
      }
  }



  // @sect4{BiLaplacianLDGLift::solve}

  // To solve the linear system $A\boldsymbol{U}=\boldsymbol{F}$,
  // we proceed as in step-74 and use a direct method. We could
  // as well use an iterative method, for instance the conjugate
  // gradient method as in step-3.
  template <int dim>
  void BiLaplacianLDGLift<dim>::solve()
  {
    SparseDirectUMFPACK A_direct;
    A_direct.initialize(matrix);
    A_direct.vmult(solution, rhs);
  }



  // @sect4{BiLaplacianLDGLift::compute_errors}

  // This function computes the discrete $H^2$, $H^1$ and $L^2$ norms of
  // the error $u-u_h$, where $u$ is the exact solution and $u_h$ is
  // the approximate solution. See the introduction for the definition
  // of the norms.
  template <int dim>
  void BiLaplacianLDGLift<dim>::compute_errors()
  {
    double error_H2 = 0;
    double error_H1 = 0;
    double error_L2 = 0;

    const QGauss<dim>     quad(fe.degree + 1);
    const QGauss<dim - 1> quad_face(fe.degree + 1);

    FEValues<dim> fe_values(fe,
                            quad,
                            update_values | update_gradients | update_hessians |
                              update_quadrature_points | update_JxW_values);

    FEFaceValues<dim> fe_face(fe,
                              quad_face,
                              update_values | update_gradients |
                                update_quadrature_points | update_JxW_values);

    FEFaceValues<dim> fe_face_neighbor(fe,
                                       quad_face,
                                       update_values | update_gradients);

    const unsigned int n_q_points      = quad.size();
    const unsigned int n_q_points_face = quad_face.size();

    // We introduce some variables for the exact solution...
    const ExactSolution<dim> u_exact;

    // ...and for the approximate solution:
    std::vector<double>         solution_values_cell(n_q_points);
    std::vector<Tensor<1, dim>> solution_gradients_cell(n_q_points);
    std::vector<Tensor<2, dim>> solution_hessians_cell(n_q_points);

    std::vector<double>         solution_values(n_q_points_face);
    std::vector<double>         solution_values_neigh(n_q_points_face);
    std::vector<Tensor<1, dim>> solution_gradients(n_q_points_face);
    std::vector<Tensor<1, dim>> solution_gradients_neigh(n_q_points_face);

    for (const auto &cell : dof_handler.active_cell_iterators())
      {
        fe_values.reinit(cell);

        fe_values.get_function_values(solution, solution_values_cell);
        fe_values.get_function_gradients(solution, solution_gradients_cell);
        fe_values.get_function_hessians(solution, solution_hessians_cell);

        // We first add the <i>bulk</i> terms.
        for (unsigned int q = 0; q < n_q_points; ++q)
          {
            const double dx = fe_values.JxW(q);

            error_H2 += (u_exact.hessian(fe_values.quadrature_point(q)) -
                         solution_hessians_cell[q])
                          .norm_square() *
                        dx;
            error_H1 += (u_exact.gradient(fe_values.quadrature_point(q)) -
                         solution_gradients_cell[q])
                          .norm_square() *
                        dx;
            error_L2 += Utilities::fixed_power<2>(
                          u_exact.value(fe_values.quadrature_point(q)) -
                          solution_values_cell[q]) *
                        dx;
          } // for quadrature points

        // We then add the face contributions.
        for (unsigned int face_no = 0; face_no < cell->n_faces(); ++face_no)
          {
            const typename DoFHandler<dim>::face_iterator face =
              cell->face(face_no);

            const double mesh_inv = 1.0 / face->diameter(); // h^{-1}
            const double mesh3_inv =
              1.0 / Utilities::fixed_power<3>(face->diameter()); // h^{-3}

            fe_face.reinit(cell, face_no);

            fe_face.get_function_values(solution, solution_values);
            fe_face.get_function_gradients(solution, solution_gradients);

            const bool at_boundary = face->at_boundary();
            if (at_boundary)
              {
                for (unsigned int q = 0; q < n_q_points_face; ++q)
                  {
                    const double dx = fe_face.JxW(q);
                    const double u_exact_q =
                      u_exact.value(fe_face.quadrature_point(q));
                    const Tensor<1, dim> u_exact_grad_q =
                      u_exact.gradient(fe_face.quadrature_point(q));

                    error_H2 +=
                      mesh_inv *
                      (u_exact_grad_q - solution_gradients[q]).norm_square() *
                      dx;
                    error_H2 += mesh3_inv *
                                Utilities::fixed_power<2>(u_exact_q -
                                                          solution_values[q]) *
                                dx;
                    error_H1 += mesh_inv *
                                Utilities::fixed_power<2>(u_exact_q -
                                                          solution_values[q]) *
                                dx;
                  }
              }
            else
              { // interior face

                const typename DoFHandler<dim>::active_cell_iterator
                                   neighbor_cell = cell->neighbor(face_no);
                const unsigned int face_no_neighbor =
                  cell->neighbor_of_neighbor(face_no);

                // In the next step, we need to have a global way to compare the
                // cells in order to not calculate the same jump term twice:
                if (neighbor_cell->id() < cell->id())
                  continue; // skip this face (already considered)
                else
                  {
                    fe_face_neighbor.reinit(neighbor_cell, face_no_neighbor);

                    fe_face.get_function_values(solution, solution_values);
                    fe_face_neighbor.get_function_values(solution,
                                                         solution_values_neigh);
                    fe_face.get_function_gradients(solution,
                                                   solution_gradients);
                    fe_face_neighbor.get_function_gradients(
                      solution, solution_gradients_neigh);

                    for (unsigned int q = 0; q < n_q_points_face; ++q)
                      {
                        const double dx = fe_face.JxW(q);

                        // To compute the jump term, we use the fact that
                        // $\jump{u}=0$ and
                        // $\jump{\nabla u}=\mathbf{0}$ since $u\in
                        // H^2(\Omega)$.
                        error_H2 +=
                          mesh_inv *
                          (solution_gradients_neigh[q] - solution_gradients[q])
                            .norm_square() *
                          dx;
                        error_H2 +=
                          mesh3_inv *
                          Utilities::fixed_power<2>(solution_values_neigh[q] -
                                                    solution_values[q]) *
                          dx;
                        error_H1 +=
                          mesh_inv *
                          Utilities::fixed_power<2>(solution_values_neigh[q] -
                                                    solution_values[q]) *
                          dx;
                      }
                  } // face not visited yet

              } // boundary check

          } // for face

      } // for cell

    error_H2 = std::sqrt(error_H2);
    error_H1 = std::sqrt(error_H1);
    error_L2 = std::sqrt(error_L2);

    std::cout << "DG H2 norm of the error: " << error_H2 << std::endl;
    std::cout << "DG H1 norm of the error: " << error_H1 << std::endl;
    std::cout << "   L2 norm of the error: " << error_L2 << std::endl;
  }



  // @sect4{BiLaplacianLDGLift::output_results}

  // This function, which writes the solution to a vtk file,
  // is copied from step-3.
  template <int dim>
  void BiLaplacianLDGLift<dim>::output_results() const
  {
    DataOut<dim> data_out;
    data_out.attach_dof_handler(dof_handler);
    data_out.add_data_vector(solution, "solution");
    data_out.build_patches();

    std::ofstream output("solution.vtk");
    data_out.write_vtk(output);
  }



  // @sect4{BiLaplacianLDGLift::assemble_local_matrix}

  // As already mentioned above, this function is used to assemble
  // the (local) mass matrices needed for the computations of the
  // lifting terms. We reiterate that only the basis functions with
  // support on the current cell are considered.
  template <int dim>
  void BiLaplacianLDGLift<dim>::assemble_local_matrix(
    const FEValues<dim> &fe_values_lift,
    const unsigned int   n_q_points,
    FullMatrix<double>  &local_matrix)
  {
    const FEValuesExtractors::Tensor<2> tau_ext(0);

    const unsigned int n_dofs = fe_values_lift.dofs_per_cell;

    local_matrix = 0;
    for (unsigned int q = 0; q < n_q_points; ++q)
      {
        const double dx = fe_values_lift.JxW(q);

        for (unsigned int m = 0; m < n_dofs; ++m)
          for (unsigned int n = 0; n < n_dofs; ++n)
            {
              local_matrix(m, n) +=
                scalar_product(fe_values_lift[tau_ext].value(n, q),
                               fe_values_lift[tau_ext].value(m, q)) *
                dx;
            }
      }
  }



  // @sect4{BiLaplacianLDGLift::compute_discrete_hessians}

  // This function is the main novelty of this program. It computes
  // the discrete Hessian $H_h(\varphi)$ for all the basis functions
  // $\varphi$ of $\mathbb{V}_h$ supported on the current cell and
  // those supported on a neighboring cell. The first argument
  // indicates the current cell (referring to the global DoFHandler
  // object), while the other two arguments are output variables that
  // are filled by this function.
  //
  // In the following, we need to evaluate finite element shape
  // functions for the `fe_lift` finite element on the current
  // cell. Like for example in step-61, this "lift" space is defined
  // on every cell individually; as a consequence, there is no global
  // DoFHandler associated with this because we simply have no need
  // for such a DoFHandler. That leaves the question of what we should
  // initialize the FEValues and FEFaceValues objects with when we ask
  // them to evaluate shape functions of `fe_lift` on a concrete
  // cell. If we simply provide the first argument to this function,
  // `cell`, to FEValues::reinit(), we will receive an error message
  // that the given `cell` belongs to a DoFHandler that has a
  // different finite element associated with it than the `fe_lift`
  // object we want to evaluate. Fortunately, there is a relatively
  // easy solution: We can call FEValues::reinit() with a cell that
  // points into a triangulation -- the same cell, but not associated
  // with a DoFHandler, and consequently no finite element space. In
  // that case, FEValues::reinit() will skip the check that would
  // otherwise lead to an error message. All we have to do is to convert
  // the DoFHandler cell iterator into a Triangulation cell iterator;
  // see the first couple of lines of the function below to see how
  // this can be done.
  template <int dim>
  void BiLaplacianLDGLift<dim>::compute_discrete_hessians(
    const typename DoFHandler<dim>::active_cell_iterator &cell,
    std::vector<std::vector<Tensor<2, dim>>>             &discrete_hessians,
    std::vector<std::vector<std::vector<Tensor<2, dim>>>>
      &discrete_hessians_neigh)
  {
    const typename Triangulation<dim>::cell_iterator cell_lift =
      static_cast<typename Triangulation<dim>::cell_iterator>(cell);

    const QGauss<dim>     quad(fe.degree + 1);
    const QGauss<dim - 1> quad_face(fe.degree + 1);

    const unsigned int n_q_points      = quad.size();
    const unsigned int n_q_points_face = quad_face.size();

    // The information we need from the basis functions of
    // $\mathbb{V}_h$: <code>fe_values</code> is needed to add
    // the broken Hessian part of the discrete Hessian, while
    // <code>fe_face</code> and <code>fe_face_neighbor</code>
    // are used to compute the right-hand sides for the local
    // problems.
    FEValues<dim> fe_values(fe, quad, update_hessians | update_JxW_values);

    FEFaceValues<dim> fe_face(
      fe, quad_face, update_values | update_gradients | update_normal_vectors);

    FEFaceValues<dim> fe_face_neighbor(
      fe, quad_face, update_values | update_gradients | update_normal_vectors);

    const unsigned int n_dofs = fe_values.dofs_per_cell;

    // The information needed from the basis functions
    // of the finite element space for the lifting terms:
    // <code>fe_values_lift</code> is used for the (local)
    // mass matrix (see $\boldsymbol{M}_c$ in the introduction),
    // while <code>fe_face_lift</code> is used to compute the
    // right-hand sides (see $\boldsymbol{G}_c$ for $b_e$).
    FEValues<dim> fe_values_lift(fe_lift,
                                 quad,
                                 update_values | update_JxW_values);

    FEFaceValues<dim> fe_face_lift(
      fe_lift, quad_face, update_values | update_gradients | update_JxW_values);

    const FEValuesExtractors::Tensor<2> tau_ext(0);

    const unsigned int n_dofs_lift = fe_values_lift.dofs_per_cell;
    FullMatrix<double> local_matrix_lift(n_dofs_lift, n_dofs_lift);

    Vector<double> local_rhs_re(n_dofs_lift), local_rhs_be(n_dofs_lift),
      coeffs_re(n_dofs_lift), coeffs_be(n_dofs_lift), coeffs_tmp(n_dofs_lift);

    SolverControl            solver_control(1000, 1e-12);
    SolverCG<Vector<double>> solver(solver_control);

    double factor_avg; // 0.5 for interior faces, 1.0 for boundary faces

    fe_values.reinit(cell);
    fe_values_lift.reinit(cell_lift);

    // We start by assembling the (local) mass matrix used for the computation
    // of the lifting terms $r_e$ and $b_e$.
    assemble_local_matrix(fe_values_lift, n_q_points, local_matrix_lift);

    for (unsigned int i = 0; i < n_dofs; ++i)
      for (unsigned int q = 0; q < n_q_points; ++q)
        {
          discrete_hessians[i][q] = 0;

          for (unsigned int face_no = 0;
               face_no < discrete_hessians_neigh.size();
               ++face_no)
            {
              discrete_hessians_neigh[face_no][i][q] = 0;
            }
        }

    // In this loop, we compute the discrete Hessian at each quadrature point
    // $x_q$ of <code>cell</code> for each basis function supported on
    // <code>cell</code>, namely we fill-in the variable
    // <code>discrete_hessians[i][q]</code>. For the lifting terms, we need to
    // add the contribution of all the faces of <code>cell</code>.
    for (unsigned int i = 0; i < n_dofs; ++i)
      {
        coeffs_re = 0;
        coeffs_be = 0;

        for (unsigned int face_no = 0; face_no < cell->n_faces(); ++face_no)
          {
            const typename DoFHandler<dim>::face_iterator face =
              cell->face(face_no);

            const bool at_boundary = face->at_boundary();

            // Recall that by convention, the average of a function across a
            // boundary face $e$ reduces to the trace of the function on the
            // only element adjacent to $e$, namely there is no factor
            // $\frac{1}{2}$. We distinguish between the two cases (the current
            // face lies in the interior or on the boundary of the domain) using
            // the variable <code>factor_avg</code>.
            factor_avg = 0.5;
            if (at_boundary)
              {
                factor_avg = 1.0;
              }

            fe_face.reinit(cell, face_no);
            fe_face_lift.reinit(cell_lift, face_no);

            local_rhs_re = 0;
            for (unsigned int q = 0; q < n_q_points_face; ++q)
              {
                const double         dx     = fe_face_lift.JxW(q);
                const Tensor<1, dim> normal = fe_face.normal_vector(
                  q); // same as fe_face_lift.normal_vector(q)

                for (unsigned int m = 0; m < n_dofs_lift; ++m)
                  {
                    local_rhs_re(m) +=
                      factor_avg *
                      (fe_face_lift[tau_ext].value(m, q) * normal) *
                      fe_face.shape_grad(i, q) * dx;
                  }
              }

            // Here, <code>local_rhs_be(m)</code> corresponds to $G_m$
            // introduced in the comments about the implementation of the
            // lifting $b_e$ in the case
            // $\varphi=\varphi^c$.
            local_rhs_be = 0;
            for (unsigned int q = 0; q < n_q_points_face; ++q)
              {
                const double         dx     = fe_face_lift.JxW(q);
                const Tensor<1, dim> normal = fe_face.normal_vector(
                  q); // same as fe_face_lift.normal_vector(q)

                for (unsigned int m = 0; m < n_dofs_lift; ++m)
                  {
                    local_rhs_be(m) += factor_avg *
                                       fe_face_lift[tau_ext].divergence(m, q) *
                                       normal * fe_face.shape_value(i, q) * dx;
                  }
              }

            coeffs_tmp = 0;
            solver.solve(local_matrix_lift,
                         coeffs_tmp,
                         local_rhs_re,
                         PreconditionIdentity());
            coeffs_re += coeffs_tmp;

            coeffs_tmp = 0;
            solver.solve(local_matrix_lift,
                         coeffs_tmp,
                         local_rhs_be,
                         PreconditionIdentity());
            coeffs_be += coeffs_tmp;

          } // for face

        for (unsigned int q = 0; q < n_q_points; ++q)
          {
            discrete_hessians[i][q] += fe_values.shape_hessian(i, q);

            for (unsigned int m = 0; m < n_dofs_lift; ++m)
              {
                discrete_hessians[i][q] -=
                  coeffs_re[m] * fe_values_lift[tau_ext].value(m, q);
              }

            for (unsigned int m = 0; m < n_dofs_lift; ++m)
              {
                discrete_hessians[i][q] +=
                  coeffs_be[m] * fe_values_lift[tau_ext].value(m, q);
              }
          }
      } // for dof i



    // In this loop, we compute the discrete Hessian at each quadrature point
    // $x_q$ of <code>cell</code> for each basis function supported on a
    // neighboring <code>neighbor_cell</code> of <code>cell</code>, namely we
    // fill-in the variable <code>discrete_hessians_neigh[face_no][i][q]</code>.
    // For the lifting terms, we only need to add the contribution of the
    // face adjacent to <code>cell</code> and <code>neighbor_cell</code>.
    for (unsigned int face_no = 0; face_no < cell->n_faces(); ++face_no)
      {
        const typename DoFHandler<dim>::face_iterator face =
          cell->face(face_no);

        const bool at_boundary = face->at_boundary();

        if (!at_boundary)
          {
            // For non-homogeneous Dirichlet BCs, we would need to
            // compute the lifting of the prescribed BC (see the
            // "Possible Extensions" section for more details).

            const typename DoFHandler<2, dim>::active_cell_iterator
                               neighbor_cell = cell->neighbor(face_no);
            const unsigned int face_no_neighbor =
              cell->neighbor_of_neighbor(face_no);
            fe_face_neighbor.reinit(neighbor_cell, face_no_neighbor);

            for (unsigned int i = 0; i < n_dofs; ++i)
              {
                coeffs_re = 0;
                coeffs_be = 0;

                fe_face_lift.reinit(cell_lift, face_no);

                local_rhs_re = 0;
                for (unsigned int q = 0; q < n_q_points_face; ++q)
                  {
                    const double         dx = fe_face_lift.JxW(q);
                    const Tensor<1, dim> normal =
                      fe_face_neighbor.normal_vector(q);

                    for (unsigned int m = 0; m < n_dofs_lift; ++m)
                      {
                        local_rhs_re(m) +=
                          0.5 * (fe_face_lift[tau_ext].value(m, q) * normal) *
                          fe_face_neighbor.shape_grad(i, q) * dx;
                      }
                  }

                // Here, <code>local_rhs_be(m)</code> corresponds to $G_m$
                // introduced in the comments about the implementation of the
                // lifting $b_e$ in the case
                // $\varphi=\varphi^n$.
                local_rhs_be = 0;
                for (unsigned int q = 0; q < n_q_points_face; ++q)
                  {
                    const double         dx = fe_face_lift.JxW(q);
                    const Tensor<1, dim> normal =
                      fe_face_neighbor.normal_vector(q);

                    for (unsigned int m = 0; m < n_dofs_lift; ++m)
                      {
                        local_rhs_be(m) +=
                          0.5 * fe_face_lift[tau_ext].divergence(m, q) *
                          normal * fe_face_neighbor.shape_value(i, q) * dx;
                      }
                  }

                solver.solve(local_matrix_lift,
                             coeffs_re,
                             local_rhs_re,
                             PreconditionIdentity());
                solver.solve(local_matrix_lift,
                             coeffs_be,
                             local_rhs_be,
                             PreconditionIdentity());

                for (unsigned int q = 0; q < n_q_points; ++q)
                  {
                    for (unsigned int m = 0; m < n_dofs_lift; ++m)
                      {
                        discrete_hessians_neigh[face_no][i][q] -=
                          coeffs_re[m] * fe_values_lift[tau_ext].value(m, q);
                      }

                    for (unsigned int m = 0; m < n_dofs_lift; ++m)
                      {
                        discrete_hessians_neigh[face_no][i][q] +=
                          coeffs_be[m] * fe_values_lift[tau_ext].value(m, q);
                      }
                  }

              } // for dof i
          }     // boundary check
      }         // for face
  }



  // @sect4{BiLaplacianLDGLift::run}
  template <int dim>
  void BiLaplacianLDGLift<dim>::run()
  {
    make_grid();

    setup_system();
    assemble_system();

    solve();

    compute_errors();
    output_results();
  }

} // namespace Step82



// @sect3{The <code>main</code> function}

// This is the <code>main</code> function. We define here the number of mesh
// refinements, the polynomial degree for the two finite element spaces
// (for the solution and the two liftings) and the two penalty coefficients.
// We can also change the dimension to run the code in 3d.
int main()
{
  try
    {
      const unsigned int n_ref = 3; // number of mesh refinements

      const unsigned int degree =
        2; // FE degree for u_h and the two lifting terms

      const double penalty_grad =
        1.0; // penalty coefficient for the jump of the gradients
      const double penalty_val =
        1.0; // penalty coefficient for the jump of the values

      Step82::BiLaplacianLDGLift<2> problem(n_ref,
                                            degree,
                                            penalty_grad,
                                            penalty_val);

      problem.run();
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
