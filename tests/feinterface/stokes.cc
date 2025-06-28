/* ------------------------------------------------------------------------
 *
 * SPDX-License-Identifier: LGPL-2.1-or-later
 * Copyright (C) 2021 - 2025 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * Part of the source code is dual licensed under Apache-2.0 WITH
 * LLVM-exception OR LGPL-2.1-or-later. Detailed license information
 * governing the source code and code contributions can be found in
 * LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
 *
 * ------------------------------------------------------------------------

 * Author:
 *         Timo Heister, Clemson University
 */

// test FEInterfaceValues for a DG Stokes problem.

#include <deal.II/base/flow_function.h>
#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/utilities.h>

#include <deal.II/dofs/block_info.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_interface_values.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/block_sparse_matrix.h>
#include <deal.II/lac/block_sparsity_pattern.h>
#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_direct.h>

#include <deal.II/meshworker/mesh_loop.h>
#include <deal.II/meshworker/scratch_data.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/vector_tools.h>

#include <fstream>
#include <sstream>

#include "../tests.h"

namespace StokesTests
{

  struct CopyDataFace
  {
    FullMatrix<double>                   cell_matrix;
    std::vector<types::global_dof_index> joint_dof_indices;
    std::array<double, 2>                values;
    std::array<unsigned int, 2>          cell_indices;
  };



  struct CopyData
  {
    FullMatrix<double>                   cell_matrix;
    Vector<double>                       cell_rhs;
    std::vector<types::global_dof_index> local_dof_indices;
    std::vector<CopyDataFace>            face_data;
    double                               value;
    unsigned int                         cell_index;



    template <class Iterator>
    void
    reinit(const Iterator &cell, const unsigned int dofs_per_cell)
    {
      cell_matrix.reinit(dofs_per_cell, dofs_per_cell);
      cell_rhs.reinit(dofs_per_cell);
      local_dof_indices.resize(dofs_per_cell);
      cell->get_dof_indices(local_dof_indices);
    }
  };



  template <class MatrixType, class VectorType>
  inline void
  copy(const CopyData                  &c,
       const AffineConstraints<double> &constraints,
       MatrixType                      &system_matrix,
       VectorType                      &system_rhs)
  {
    constraints.distribute_local_to_global(c.cell_matrix,
                                           c.cell_rhs,
                                           c.local_dof_indices,
                                           system_matrix,
                                           system_rhs);
    for (auto &cdf : c.face_data)
      {
        const unsigned int dofs_per_cell = cdf.joint_dof_indices.size();
        for (unsigned int i = 0; i < dofs_per_cell; ++i)
          for (unsigned int k = 0; k < dofs_per_cell; ++k)
            system_matrix.add(cdf.joint_dof_indices[i],
                              cdf.joint_dof_indices[k],
                              cdf.cell_matrix(i, k));
      }
  }

  // @sect3{Functions for Solution and Righthand side}
  //
  // The class Solution is used to define the boundary conditions and to
  // compute errors of the numerical solution. Note that we need to define the
  // values and gradients in order to compute L2 and H1 errors. Here we
  // decided to separate the implementations for 2d and 3d using template
  // specialization.
  //
  // Note that the first dim components are the velocity components
  // and the last is the pressure.
  template <int dim>
  class Solution : public Function<dim>
  {
  public:
    Solution()
      : Function<dim>(dim + 1)
    {}
    virtual double
    value(const Point<dim> &p, const unsigned int component = 0) const;
    virtual Tensor<1, dim>
    gradient(const Point<dim> &p, const unsigned int component = 0) const;
  };

  template <>
  double
  Solution<2>::value(const Point<2> &p, const unsigned int component) const
  {
    using numbers::PI;
    const double x = p[0];
    const double y = p[1];
    // zero on BD's
    if (component == 0)
      return PI * sin(PI * x) * sin(PI * x) * sin(2.0 * PI * y);
    if (component == 1)
      return -PI * sin(PI * y) * sin(PI * y) * sin(2.0 * PI * x);
    if (component == 2)
      return cos(PI * x) * sin(PI * y);

    return 0;
  }

  template <>
  double
  Solution<3>::value(const Point<3> &p, const unsigned int component) const
  {
    Assert(component <= 3 + 1, ExcIndexRange(component, 0, 3 + 1));

    using numbers::PI;
    const double x = p[0];
    const double y = p[1];
    const double z = p[2];

    if (component == 0)
      return 2. * PI * sin(PI * x) * sin(PI * x) * sin(2.0 * PI * y) *
             sin(2.0 * PI * z);
    if (component == 1)
      return -PI * sin(PI * y) * sin(PI * y) * sin(2.0 * PI * x) *
             sin(2.0 * PI * z);
    if (component == 2)
      return -PI * sin(PI * z) * sin(PI * z) * sin(2.0 * PI * x) *
             sin(2.0 * PI * y);
    if (component == 3)
      return sin(PI * x) * cos(PI * y) * sin(PI * z);

    return 0;
  }

  // Note that for the gradient we need to return a Tensor<1,dim>
  template <>
  Tensor<1, 2>
  Solution<2>::gradient(const Point<2> &p, const unsigned int component) const
  {
    Assert(component <= 2, ExcIndexRange(component, 0, 2 + 1));

    using numbers::PI;
    const double x = p[0];
    const double y = p[1];

    Tensor<1, 2> return_value;
    if (component == 0)
      {
        return_value[0] = PI * PI * sin(2.0 * PI * y) * sin(2.0 * PI * x);
        return_value[1] =
          2.0 * PI * PI * sin(PI * x) * sin(PI * x) * cos(2.0 * PI * y);
      }
    else if (component == 1)
      {
        return_value[0] =
          -2.0 * PI * PI * sin(PI * y) * sin(PI * y) * cos(2.0 * PI * x);
        return_value[1] = -PI * PI * sin(2.0 * PI * y) * sin(2.0 * PI * x);
      }
    else if (component == 2)
      {
        return_value[0] = PI * cos(PI * x) * cos(PI * y);
        return_value[1] = -PI * sin(PI * x) * sin(PI * y);
      }

    return return_value;
  }

  template <>
  Tensor<1, 3>
  Solution<3>::gradient(const Point<3> &p, const unsigned int component) const
  {
    Assert(component <= 3, ExcIndexRange(component, 0, 3 + 1));

    using numbers::PI;
    const double x = p[0];
    const double y = p[1];
    const double z = p[2];

    Tensor<1, 3> return_value;
    if (component == 0)
      {
        return_value[0] =
          2 * PI * PI * sin(2 * PI * x) * sin(2 * PI * y) * sin(2 * PI * z);
        return_value[1] = 4 * PI * PI * sin(PI * x) * sin(PI * x) *
                          cos(2 * PI * y) * sin(2 * PI * z);
        return_value[2] = 4 * PI * PI * sin(PI * x) * sin(PI * x) *
                          cos(2 * PI * z) * sin(2 * PI * y);
      }
    else if (component == 1)
      {
        return_value[0] = -2 * PI * PI * sin(PI * y) * sin(PI * y) *
                          cos(2 * PI * x) * sin(2 * PI * z);
        return_value[1] =
          -PI * PI * sin(2 * PI * x) * sin(2 * PI * y) * sin(2 * PI * z);
        return_value[2] = -2 * PI * PI * sin(PI * y) * sin(PI * y) *
                          cos(2 * PI * z) * sin(2 * PI * x);
      }
    else if (component == 2)
      {
        return_value[0] = -2 * PI * PI * sin(PI * z) * sin(PI * z) *
                          cos(2 * PI * x) * sin(2 * PI * y);
        return_value[1] = -2 * PI * PI * sin(PI * z) * sin(PI * z) *
                          cos(2 * PI * y) * sin(2 * PI * x);
        return_value[2] =
          -PI * PI * sin(2 * PI * x) * sin(2 * PI * y) * sin(2 * PI * z);
      }
    else if (component == 3)
      {
        return_value[0] = PI * cos(PI * x) * cos(PI * y) * sin(PI * z);
        return_value[1] = -PI * sin(PI * x) * sin(PI * y) * sin(PI * z);
        return_value[2] = PI * sin(PI * x) * cos(PI * y) * cos(PI * z);
      }

    return return_value;
  }



  // Implementation of $f$. See the introduction for more information.
  template <int dim>
  class RightHandSide : public Function<dim>
  {
  public:
    RightHandSide()
      : Function<dim>(dim + 1)
    {}

    virtual double
    value(const Point<dim> &p, const unsigned int component = 0) const;
  };

  template <>
  double
  RightHandSide<2>::value(const Point<2> &p, const unsigned int component) const
  {
    Assert(component <= 2, ExcIndexRange(component, 0, 2 + 1));

    using numbers::PI;
    double x  = p[0];
    double y  = p[1];
    double nu = 1.0;

    // RHS for 0 BD's
    if (component == 0)
      return -nu * 2.0 * PI * PI * PI *
               (-2.0 * sin(PI * x) * sin(PI * x) + cos(2. * PI * x)) *
               sin(2.0 * PI * y) -
             PI * sin(PI * x) * sin(PI * y);
    if (component == 1)
      return nu * 2.0 * PI * PI * PI * (2.0 * cos(2.0 * PI * y) - 1) *
               sin(2.0 * PI * x) +
             PI * cos(PI * x) * cos(PI * y);
    if (component == 2)
      return 0.0;

    return 0.0;
  }

  template <>
  double
  RightHandSide<3>::value(const Point<3> &p, const unsigned int component) const
  {
    Assert(component <= 3, ExcIndexRange(component, 0, 3 + 1));

    using numbers::PI;
    double x = p[0];
    double y = p[1];
    double z = p[2];

    if (component == 0)
      return 4. * PI * PI * PI *
               (4. * sin(PI * x) * sin(PI * x) - cos(2. * PI * x)) *
               sin(2. * PI * y) * sin(2. * PI * z) +
             PI * cos(PI * x) * cos(PI * y) * sin(PI * z);
    if (component == 1)
      return -2. * PI * PI * PI *
               (4. * sin(PI * y) * sin(PI * y) - cos(2. * PI * y)) *
               sin(2. * PI * x) * sin(2. * PI * z) +
             PI * (-1) * sin(PI * y) * sin(PI * x) * sin(PI * z);
    if (component == 2)
      return -2. * PI * PI * PI *
               (4. * sin(PI * z) * sin(PI * z) - cos(2. * PI * z)) *
               sin(2. * PI * x) * sin(2. * PI * y) +
             PI * cos(PI * z) * sin(PI * x) * cos(PI * y);
    if (component == 3)
      return 0.0;

    return 0.0;
  }



  // @sect3{The StokesProblem class}
  //
  // This is the main class of the problem.
  template <int dim>
  class StokesProblem
  {
  public:
    StokesProblem(FiniteElement<dim> &fe, const unsigned int pressure_degree);
    void
    run();

  private:
    void
    setup_dofs();
    void
    assemble_system_mesh_loop();
    void
    solve();
    void
    compute_errors(unsigned int k);

    const unsigned int pressure_degree;

    Triangulation<dim>  triangulation;
    FiniteElement<dim> &fe;
    DoFHandler<dim>     dof_handler;

    AffineConstraints<> constraints;

    BlockSparsityPattern      sparsity_pattern;
    BlockSparseMatrix<double> system_matrix;
    SparseMatrix<double>      pressure_mass_matrix;

    BlockVector<double> solution;
    BlockVector<double> system_rhs;

    double last_l2_error;
    double last_H1_error;
    double last_Hdiv_error1;
    double last_Hdiv_error2;
  };



  template <int dim>
  StokesProblem<dim>::StokesProblem(FiniteElement<dim> &fe,
                                    const unsigned int  pressure_degree)
    : pressure_degree(pressure_degree)
    , triangulation(Triangulation<dim>::maximum_smoothing)
    , fe(fe)
    , dof_handler(triangulation)
  {}

  // @sect4{StokesProblem::setup_dofs}

  // This function sets up the DoFHandler, matrices, vectors, and Multigrid
  // structures (if needed).
  template <int dim>
  void
  StokesProblem<dim>::setup_dofs()
  {
    system_matrix.clear();
    pressure_mass_matrix.clear();

    // The main DoFHandler only needs active DoFs, so we are not calling
    // distribute_mg_dofs() here
    dof_handler.distribute_dofs(fe);

    // This block structure separates the dim velocity components from
    // the pressure component (used for reordering). Note that we have
    // 2 instead of dim+1 blocks like in step-22, because our FESystem
    // is nested and the dim velocity components appear as one block.
    std::vector<unsigned int> block_component(2);
    block_component[0] = 0;
    block_component[1] = 1;

    // Velocities start at component 0:
    const FEValuesExtractors::Vector velocities(0);

    // ILU behaves better if we apply a reordering to reduce filling. There
    // is no advantage in doing this for the other solvers.
    DoFRenumbering::Cuthill_McKee(dof_handler);


    // This ensures that all velocities DoFs are enumerated before the
    // pressure unknowns. This allows us to use blocks for vectors and
    // matrices and allows us to get the same DoF numbering for
    // dof_handler and velocity_dof_handler.
    DoFRenumbering::block_wise(dof_handler);

    std::vector<types::global_dof_index> dofs_per_block =
      DoFTools::count_dofs_per_fe_block(dof_handler, block_component);
    const unsigned int n_u = dofs_per_block[0], n_p = dofs_per_block[1];

    deallog << "\tNumber of active cells: " << triangulation.n_active_cells()
            << std::endl
            << "\tNumber of degrees of freedom: " << dof_handler.n_dofs()
            << " (" << n_u << '+' << n_p << ')' << std::endl;

    {
      constraints.reinit();
      DoFTools::make_hanging_node_constraints(dof_handler, constraints);
      VectorTools::interpolate_boundary_values(dof_handler,
                                               0,
                                               Solution<dim>(),
                                               constraints,
                                               fe.component_mask(velocities));
      constraints.close();
    }

    {
      BlockDynamicSparsityPattern csp(dofs_per_block, dofs_per_block);
      DoFTools::make_flux_sparsity_pattern(dof_handler, csp, constraints);
      sparsity_pattern.copy_from(csp);
    }
    system_matrix.reinit(sparsity_pattern);

    solution.reinit(dofs_per_block);
    system_rhs.reinit(dofs_per_block);
  }



  template <int dim>
  void
  StokesProblem<dim>::assemble_system_mesh_loop()
  {
    system_matrix = 0;
    system_rhs    = 0;

    typedef decltype(dof_handler.begin_active()) Iterator;
    const RightHandSide<dim>                     rhs_function;
    const Solution<dim>                          boundary_function;
    const FEValuesExtractors::Vector             velocities(0);
    const FEValuesExtractors::Scalar             pressure(dim);

    const double nu = 1.0;

    auto penalty_parameter = [](const double degree,
                                const double extent1,
                                const double extent2) -> double {
      return 4.0 * degree * (degree + 1.0) * 0.5 *
             (1.0 / extent1 + 1.0 / extent2);
    };

    auto cell_worker = [&](const Iterator               &cell,
                           MeshWorker::ScratchData<dim> &scratch_data,
                           CopyData                     &copy_data) {
      const FEValues<dim> &fe_v = scratch_data.reinit(cell);

      const unsigned int dofs_per_cell = fe_v.dofs_per_cell;
      const unsigned int n_q_points    = fe_v.get_quadrature().size();

      copy_data.reinit(cell, dofs_per_cell);

      const std::vector<double> &JxW = fe_v.get_JxW_values();

      const double                nu = 1.0;
      std::vector<Vector<double>> rhs_values(n_q_points,
                                             Vector<double>(dim + 1));
      rhs_function.vector_value_list(fe_v.get_quadrature_points(), rhs_values);
      Tensor<1, dim> force_f;

      for (unsigned int point = 0; point < fe_v.n_quadrature_points; ++point)
        {
          for (unsigned int d = 0; d < dim; ++d)
            force_f[d] = rhs_values[point](d);
          for (unsigned int i = 0; i < fe_v.dofs_per_cell; ++i)
            {
              for (unsigned int j = 0; j < fe_v.dofs_per_cell; ++j)
                copy_data.cell_matrix(i, j) +=
                  (
                    // nu \nabla v : \nabla u
                    nu * scalar_product(fe_v[velocities].gradient(i, point),
                                        fe_v[velocities].gradient(j, point))
                    // -q, div u
                    - fe_v[pressure].value(i, point) *
                        fe_v[velocities].divergence(j, point)
                    // -p, div v
                    - fe_v[pressure].value(j, point) *
                        fe_v[velocities].divergence(i, point)
                    // p,q
                    + fe_v[pressure].value(j, point) *
                        fe_v[pressure].value(i, point)) *
                  JxW[point];

              copy_data.cell_rhs(i) +=
                // f,v
                (force_f * fe_v[velocities].value(i, point)) * JxW[point];
            }
        }
    };

    auto boundary_worker = [&](const Iterator               &cell,
                               const unsigned int           &face_no,
                               MeshWorker::ScratchData<dim> &scratch_data,
                               CopyData                     &copy_data) {
      const FEFaceValuesBase<dim> &fe_fv = scratch_data.reinit(cell, face_no);

      const auto &q_points = fe_fv.get_quadrature_points();

      const std::vector<double>         &JxW     = fe_fv.get_JxW_values();
      const std::vector<Tensor<1, dim>> &normals = fe_fv.get_normal_vectors();

      std::vector<Vector<double>> g_values(q_points.size(),
                                           Vector<double>(dim + 1));
      boundary_function.vector_value_list(q_points, g_values);
      Tensor<1, dim> g;

      const double degree =
        std::max(1.0, static_cast<double>(fe_fv.get_fe().degree));
      const double extent1 = cell->extent_in_direction(
        GeometryInfo<dim>::unit_normal_direction[face_no]);
      const double penalty = penalty_parameter(degree, extent1, extent1);

      for (unsigned int point = 0; point < q_points.size(); ++point)
        {
          for (unsigned int d = 0; d < dim; ++d)
            g[d] = g_values[point](d);

          for (unsigned int i = 0; i < fe_fv.dofs_per_cell; ++i)
            for (unsigned int j = 0; j < fe_fv.dofs_per_cell; ++j)
              copy_data.cell_matrix(i, j) +=
                (
                  // - nu (\nabla u n) . v
                  -nu *
                    ((fe_fv[velocities].gradient(j, point) * normals[point]) *
                     fe_fv[velocities].value(i, point))

                  // - nu u . (\nabla v n)  // NIPG: use +
                  -
                  nu * (fe_fv[velocities].value(j, point) *
                        (fe_fv[velocities].gradient(i, point) * normals[point]))

                  // + nu * penalty u . v
                  + nu * penalty *
                      (fe_fv[velocities].value(j, point) *
                       fe_fv[velocities].value(i, point))

                  // p (v.n)
                  + fe_fv[pressure].value(j, point) *
                      scalar_product(fe_fv[velocities].value(i, point),
                                     normals[point])

                  // q (u.n)
                  + fe_fv[pressure].value(i, point) *
                      scalar_product(fe_fv[velocities].value(j, point),
                                     normals[point])

                    ) *
                JxW[point];

          for (unsigned int i = 0; i < fe_fv.dofs_per_cell; ++i)
            copy_data.cell_rhs(i) +=
              (
                // -nu g . (\nabla v n) // NIPG: use +
                -nu * scalar_product(g,
                                     (fe_fv[velocities].gradient(i, point) *
                                      normals[point]))

                // +nu penalty g . v
                + nu * penalty *
                    scalar_product(g, fe_fv[velocities].value(i, point))

                // q (g.n) (weak normal component of boundary condition)
                + fe_fv[pressure].value(i, point) *
                    scalar_product(g, normals[point])) *
              JxW[point];
        }
    };

    auto face_worker = [&](const Iterator               &cell,
                           const unsigned int           &f,
                           const unsigned int           &sf,
                           const Iterator               &ncell,
                           const unsigned int           &nf,
                           const unsigned int           &nsf,
                           MeshWorker::ScratchData<dim> &scratch_data,
                           CopyData                     &copy_data) {
      const FEInterfaceValues<dim> &fe_fv =
        scratch_data.reinit(cell, f, sf, ncell, nf, nsf);
      FEInterfaceViews::Scalar<dim> interface_scalar(fe_fv, pressure.component);
      FEInterfaceViews::Vector<dim> interface_vector(
        fe_fv, velocities.first_vector_component);

      copy_data.face_data.emplace_back();
      CopyDataFace      &copy_data_face = copy_data.face_data.back();
      const unsigned int dofs_per_cell  = fe_fv.n_current_interface_dofs();

      copy_data_face.joint_dof_indices = fe_fv.get_interface_dof_indices();
      copy_data_face.cell_matrix.reinit(dofs_per_cell, dofs_per_cell);

      const std::vector<double>         &JxW     = fe_fv.get_JxW_values();
      const std::vector<Tensor<1, dim>> &normals = fe_fv.get_normal_vectors();
      const auto &q_points = fe_fv.get_quadrature_points();

      double       nu = 1.0;
      const double degree =
        std::max(1.0, static_cast<double>(fe_fv.get_fe().degree));
      const double extent1 = cell->measure() / cell->face(f)->measure();
      const double extent2 = ncell->measure() / ncell->face(nf)->measure();
      const double penalty = penalty_parameter(degree, extent1, extent2);

      for (unsigned int point = 0; point < q_points.size(); ++point)
        {
          for (unsigned int i = 0; i < dofs_per_cell; ++i)
            for (unsigned int j = 0; j < dofs_per_cell; ++j)
              copy_data_face.cell_matrix(i, j) +=
                (
                  // - nu {\nabla u}n . [v] (consistency)
                  -nu *
                    (fe_fv[velocities].average_of_gradients(j, point) *
                     normals[point]) *
                    fe_fv[velocities].jump_in_values(i, point)

                  // - nu [u] . {\nabla v}n  (symmetry) // NIPG: use +
                  - nu * fe_fv[velocities].jump_in_values(j, point) *
                      (fe_fv[velocities].average_of_gradients(i, point) *
                       normals[point])

                  // nu sigma [u].[v] (penalty)
                  + nu * penalty *
                      scalar_product(fe_fv[velocities].jump_in_values(j, point),
                                     fe_fv[velocities].jump_in_values(i, point))

                  // {p} ([v].n)
                  + fe_fv[pressure].average_of_values(j, point) *
                      scalar_product(fe_fv[velocities].jump_in_values(i, point),
                                     normals[point])

                  // {q} ([u].n)
                  + fe_fv[pressure].average_of_values(i, point) *
                      scalar_product(fe_fv[velocities].jump_in_values(j, point),
                                     normals[point])) *
                JxW[point];
        }
    };

    auto copier = [&](const CopyData &c) {
      copy(c, constraints, system_matrix, system_rhs);
    };

    const unsigned int n_gauss_points = pressure_degree + 2;
    const UpdateFlags  cell_flags     = update_values | update_gradients |
                                   update_quadrature_points | update_JxW_values;
    const UpdateFlags face_flags = update_values | update_gradients |
                                   update_quadrature_points |
                                   update_normal_vectors | update_JxW_values;

    const QGauss<dim>     quadrature(n_gauss_points);
    const QGauss<dim - 1> face_quadrature(n_gauss_points);

    static MappingQ1<dim>        mapping;
    MeshWorker::ScratchData<dim> scratch_data(
      mapping, fe, quadrature, cell_flags, face_quadrature, face_flags);
    CopyData cd;
    MeshWorker::mesh_loop(dof_handler.begin_active(),
                          dof_handler.end(),
                          cell_worker,
                          copier,
                          scratch_data,
                          cd,
                          MeshWorker::assemble_own_cells |
                            MeshWorker::assemble_boundary_faces |
                            MeshWorker::assemble_own_interior_faces_once,
                          boundary_worker,
                          face_worker);

    pressure_mass_matrix.reinit(sparsity_pattern.block(1, 1));
    pressure_mass_matrix.copy_from(system_matrix.block(1, 1));
    system_matrix.block(1, 1) = 0;
  }



  template <int dim>
  void
  StokesProblem<dim>::solve()
  {
    solution = 0;
    SparseDirectUMFPACK A_direct;
    A_direct.initialize(system_matrix);
    A_direct.vmult(solution, system_rhs);

    constraints.distribute(solution);
  }



  // @sect4{StokesProblem::process_solution}

  // This function computes the L2 and H1 errors of the solution. For this,
  // we need to make sure the pressure has mean zero.
  template <int dim>
  void
  StokesProblem<dim>::compute_errors(unsigned int k)
  {
    const ComponentSelectFunction<dim> velocity_mask(std::make_pair(0, dim),
                                                     dim + 1);
    const ComponentSelectFunction<dim> pressure_mask(dim, 1.0, dim + 1);

    Vector<float> difference_per_cell(triangulation.n_active_cells());
    VectorTools::integrate_difference(dof_handler,
                                      solution,
                                      Solution<dim>(),
                                      difference_per_cell,
                                      QGauss<dim>(pressure_degree + 3),
                                      VectorTools::L2_norm,
                                      &velocity_mask);

    const double Velocity_L2_error = difference_per_cell.l2_norm();

    VectorTools::integrate_difference(dof_handler,
                                      solution,
                                      Solution<dim>(),
                                      difference_per_cell,
                                      QGauss<dim>(pressure_degree + 3),
                                      VectorTools::H1_norm,
                                      &velocity_mask);

    const double Velocity_H1_error = difference_per_cell.l2_norm();

    VectorTools::integrate_difference(dof_handler,
                                      solution,
                                      Functions::ZeroFunction<dim>(dim + 1),
                                      difference_per_cell,
                                      QGauss<dim>(pressure_degree + 3),
                                      VectorTools::Hdiv_seminorm,
                                      &velocity_mask);

    const double Velocity_Hdiv_error1 = difference_per_cell.l2_norm();

    static double last_Pressure_L2_error = 0;

    const double mean_pressure = VectorTools::compute_mean_value(
      dof_handler, QGauss<dim>(pressure_degree + 3), solution, dim);
    solution.block(1).add(-mean_pressure);

    VectorTools::integrate_difference(dof_handler,
                                      solution,
                                      Solution<dim>(),
                                      difference_per_cell,
                                      QGauss<dim>(pressure_degree + 3),
                                      VectorTools::L2_norm,
                                      &pressure_mask);
    const double Pressure_L2_error =
      VectorTools::compute_global_error(triangulation,
                                        difference_per_cell,
                                        VectorTools::L2_norm);

    deallog << " At " << k + 1 << "th mesh" << std::endl
            << " L2 error:     " << std::setw(12) << Velocity_L2_error
            << std::setw(0) << " L2_Conv_rate: " << std::setw(6)
            << (k == 0 ? 0 : last_l2_error / Velocity_L2_error) << std::endl
            << " H1 error:     " << std::setw(12) << Velocity_H1_error
            << std::setw(0) << " H1_Conv_rate: " << std::setw(6)
            << (k == 0 ? 0 : last_H1_error / Velocity_H1_error) << std::endl
            << " Hdiv error1:  " << std::setw(12) << Velocity_Hdiv_error1
            << std::setw(0) << " Hdiv_Conv_rate1: " << std::setw(6)
            << (k == 0 ? 0 : last_Hdiv_error1 / Velocity_Hdiv_error1)
            << std::endl
            << " L2 pressure:  " << std::setw(12) << Pressure_L2_error
            << std::setw(0) << " rate: " << std::setw(6)
            << (k == 0 ? 0 : last_Pressure_L2_error / Pressure_L2_error)
            << std::endl
            << std::setw(0) << std::endl;
    last_l2_error          = Velocity_L2_error;
    last_H1_error          = Velocity_H1_error;
    last_Hdiv_error1       = Velocity_Hdiv_error1;
    last_Pressure_L2_error = Pressure_L2_error;
  }



  // @sect4{StokesProblem::run}

  // The last step in the Stokes class is, as usual, the function that
  // generates the initial grid and calls the other functions in the
  // respective order.
  template <int dim>
  void
  StokesProblem<dim>::run()
  {
    GridGenerator::hyper_cube(triangulation);

    triangulation.refine_global(1);

    deallog << "  Now running with " << fe.get_name() << std::endl;

    for (unsigned int refinement_cycle = 0; refinement_cycle < 5;
         ++refinement_cycle)
      {
        if (refinement_cycle > 0)
          triangulation.refine_global();

        setup_dofs();

        int assemble_type = 2;

        assemble_system_mesh_loop();

        solve();

        compute_errors(refinement_cycle);
      }
  }
} // namespace StokesTests



int
main()
{
  try
    {
      using namespace StokesTests;
      const int dim = 2;

      initlog();

      std::unique_ptr<FiniteElement<dim>> fe;
      const int                           degree = 2;

      Assert(degree >= 1, ExcMessage("invalid degree!"));
      fe = std::make_unique<FESystem<dim>>(
        FESystem<dim>(FE_DGQ<dim>(degree), dim), 1, FE_DGQ<dim>(degree - 1), 1);

      deallog << fe->get_name() << ": degree=" << fe->degree
              << " tensor_degree=" << fe->tensor_degree() << std::endl;
      StokesProblem<dim> flow_problem(*fe.get(), degree);

      flow_problem.run();
    }
  catch (const std::exception &exc)
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
