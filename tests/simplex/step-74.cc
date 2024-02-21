// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2021 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Step-74 on a simplex mesh. Following incompatible modifications had to be
// made:
//  - The functions GridGenerator::hyper_cube() and GridGenerator::hyper_cube()
//    are only working for hypercube meshes. Use these functions to create
//    temporary meshes and convert them to simplex mesh.
//  - Local refinement is replaced by global refinement.


#include <deal.II/base/convergence_table.h>
#include <deal.II/base/function.h>
#include <deal.II/base/function_lib.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_interface_values.h>
#include <deal.II/fe/fe_pyramid_p.h>
#include <deal.II/fe/fe_simplex_p.h>
#include <deal.II/fe/fe_simplex_p_bubbles.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_wedge_p.h>
#include <deal.II/fe/mapping_fe.h>
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/sparse_direct.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>

#include <deal.II/meshworker/copy_data.h>
#include <deal.II/meshworker/mesh_loop.h>
#include <deal.II/meshworker/scratch_data.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/derivative_approximation.h>
#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"

#define USE_SIMPLEX

namespace Step74
{

  enum class Test_Case
  {
    convergence_rate,
    l_singularity
  };

  template <int dim>
  class SmoothSolution : public Function<dim>
  {
  public:
    SmoothSolution()
      : Function<dim>()
    {}
    virtual void
    value_list(const std::vector<Point<dim>> &points,
               std::vector<double>           &values,
               const unsigned int             component = 0) const override;
    virtual Tensor<1, dim>
    gradient(const Point<dim>  &point,
             const unsigned int component = 0) const override;
  };

  template <int dim>
  void
  SmoothSolution<dim>::value_list(const std::vector<Point<dim>> &points,
                                  std::vector<double>           &values,
                                  const unsigned int /*component*/) const
  {
    using numbers::PI;
    for (unsigned int i = 0; i < values.size(); ++i)
      values[i] =
        std::sin(2. * PI * points[i][0]) * std::sin(2. * PI * points[i][1]);
  }

  template <int dim>
  Tensor<1, dim>
  SmoothSolution<dim>::gradient(const Point<dim> &point,
                                const unsigned int /*component*/) const
  {
    Tensor<1, dim> return_value;
    using numbers::PI;
    return_value[0] =
      2. * PI * std::cos(2. * PI * point[0]) * std::sin(2. * PI * point[1]);
    return_value[1] =
      2. * PI * std::sin(2. * PI * point[0]) * std::cos(2. * PI * point[1]);
    return return_value;
  }

  template <int dim>
  class SmoothRightHandSide : public Function<dim>
  {
  public:
    SmoothRightHandSide()
      : Function<dim>()
    {}
    virtual void
    value_list(const std::vector<Point<dim>> &points,
               std::vector<double>           &values,
               const unsigned int /*component*/) const override;
  };

  template <int dim>
  void
  SmoothRightHandSide<dim>::value_list(const std::vector<Point<dim>> &points,
                                       std::vector<double>           &values,
                                       const unsigned int /*component*/) const
  {
    using numbers::PI;
    for (unsigned int i = 0; i < values.size(); ++i)
      values[i] = 8. * PI * PI * std::sin(2. * PI * points[i][0]) *
                  std::sin(2. * PI * points[i][1]);
  }

  template <int dim>
  class SingularRightHandSide : public Function<dim>
  {
  public:
    SingularRightHandSide()
      : Function<dim>()
    {}
    virtual void
    value_list(const std::vector<Point<dim>> &points,
               std::vector<double>           &values,
               const unsigned int /*component*/) const override;

  private:
    Functions::LSingularityFunction ref;
  };

  template <int dim>
  void
  SingularRightHandSide<dim>::value_list(const std::vector<Point<dim>> &points,
                                         std::vector<double>           &values,
                                         const unsigned int /*component*/) const
  {
    for (unsigned int i = 0; i < values.size(); ++i)
      values[i] = -ref.laplacian(points[i]);
  }

  double
  compute_penalty(const unsigned int fe_degree,
                  const double       cell_extend_left,
                  const double       cell_extend_right)
  {
    const double degree = std::max(1., static_cast<double>(fe_degree));
    return degree * (degree + 1.) * 0.5 *
           (1. / cell_extend_left + 1. / cell_extend_right);
  }


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
    reinit(const Iterator &cell, unsigned int dofs_per_cell)
    {
      cell_matrix.reinit(dofs_per_cell, dofs_per_cell);
      cell_rhs.reinit(dofs_per_cell);
      local_dof_indices.resize(dofs_per_cell);
      cell->get_dof_indices(local_dof_indices);
    }
  };

  template <int dim>
  class SIPGLaplace
  {
  public:
    SIPGLaplace(const Test_Case &test_case);
    void
    run();

  private:
    void
    setup_system();
    void
    assemble_system();
    void
    solve();
    void
    refine_grid();
    void
    output_results(const unsigned int cycle) const;

    void
    compute_errors();
    void
    compute_error_estimate();
    double
    compute_energy_norm();

    Triangulation<dim> triangulation;
    const unsigned     degree;

#ifdef USE_SIMPLEX
    const QGaussSimplex<dim>     quadrature;
    const QGaussSimplex<dim - 1> face_quadrature;
    const QGaussSimplex<dim>     quadrature_2;
    const QGaussSimplex<dim - 1> face_quadrature_2;
    const MappingFE<dim>         mapping;
    const FE_SimplexDGP<dim>     fe;
#else
    const QGauss<dim>     quadrature;
    const QGauss<dim - 1> face_quadrature;
    const QGauss<dim>     quadrature_2;
    const QGauss<dim - 1> face_quadrature_2;
    const MappingQ1<dim>  mapping;
    const FE_DGQ<dim>     fe;
#endif

    using ScratchData = MeshWorker::ScratchData<dim>;

    DoFHandler<dim> dof_handler;

    SparsityPattern      sparsity_pattern;
    SparseMatrix<double> system_matrix;
    Vector<double>       solution;
    Vector<double>       system_rhs;

    Vector<double> estimated_error_square_per_cell;
    Vector<double> energy_norm_square_per_cell;

    ConvergenceTable convergence_table;

    const double diffusion_coefficient = 1.;

    const Test_Case test_case;

    std::unique_ptr<Function<dim>> exact_solution;
    std::unique_ptr<Function<dim>> rhs_function;
  };

  template <int dim>
  SIPGLaplace<dim>::SIPGLaplace(const Test_Case &test_case)
    : degree(2)
    , quadrature(degree + 1)
    , face_quadrature(degree + 1)
    , quadrature_2(degree + 1)
    , face_quadrature_2(degree + 1)
#ifdef USE_SIMPLEX
    , mapping(FE_SimplexP<dim>(1))
#else
    , mapping()
#endif
    , fe(degree)
    , dof_handler(triangulation)
    , test_case(test_case)
  {
    if (test_case == Test_Case::convergence_rate)
      {
        exact_solution = std::make_unique<SmoothSolution<dim>>();
        rhs_function   = std::make_unique<SmoothRightHandSide<dim>>();
      }

    else if (test_case == Test_Case::l_singularity)
      {
        exact_solution = std::make_unique<Functions::LSingularityFunction>();
        rhs_function   = std::make_unique<SingularRightHandSide<dim>>();
      }
    else
      AssertThrow(false, ExcNotImplemented());
  }

  template <int dim>
  void
  SIPGLaplace<dim>::setup_system()
  {
    dof_handler.distribute_dofs(fe);
    DynamicSparsityPattern dsp(dof_handler.n_dofs());
    DoFTools::make_flux_sparsity_pattern(dof_handler, dsp);
    sparsity_pattern.copy_from(dsp);

    system_matrix.reinit(sparsity_pattern);
    solution.reinit(dof_handler.n_dofs());
    system_rhs.reinit(dof_handler.n_dofs());
  }

  template <int dim>
  void
  SIPGLaplace<dim>::assemble_system()
  {
    const auto cell_worker =
      [&](const auto &cell, auto &scratch_data, auto &copy_data) {
        const FEValues<dim> &fe_v          = scratch_data.reinit(cell);
        const unsigned int   dofs_per_cell = fe_v.dofs_per_cell;
        copy_data.reinit(cell, dofs_per_cell);

        const auto        &q_points    = scratch_data.get_quadrature_points();
        const unsigned int n_q_points  = q_points.size();
        const std::vector<double> &JxW = scratch_data.get_JxW_values();

        std::vector<double> rhs(n_q_points);
        rhs_function->value_list(q_points, rhs);

        for (unsigned int point = 0; point < n_q_points; ++point)
          for (unsigned int i = 0; i < fe_v.dofs_per_cell; ++i)
            {
              for (unsigned int j = 0; j < fe_v.dofs_per_cell; ++j)
                copy_data.cell_matrix(i, j) +=
                  diffusion_coefficient *     // nu
                  fe_v.shape_grad(i, point) * // grad v_h
                  fe_v.shape_grad(j, point) * // grad u_h
                  JxW[point];                 // dx

              copy_data.cell_rhs(i) += fe_v.shape_value(i, point) * // v_h
                                       rhs[point] *                 // f
                                       JxW[point];                  // dx
            }
      };

    const auto boundary_worker = [&](const auto         &cell,
                                     const unsigned int &face_no,
                                     auto               &scratch_data,
                                     auto               &copy_data) {
      const FEFaceValuesBase<dim> &fe_fv = scratch_data.reinit(cell, face_no);

      Assert(fe_fv.get_cell() == cell, ExcInternalError());
      Assert(fe_fv.get_face_number() == face_no, ExcInternalError());

      const auto        &q_points      = scratch_data.get_quadrature_points();
      const unsigned int n_q_points    = q_points.size();
      const unsigned int dofs_per_cell = fe_fv.dofs_per_cell;

      const std::vector<double>         &JxW = scratch_data.get_JxW_values();
      const std::vector<Tensor<1, dim>> &normals =
        scratch_data.get_normal_vectors();

      std::vector<double> g(n_q_points);
      exact_solution->value_list(q_points, g);

      const double extent1 = cell->measure() / cell->face(face_no)->measure();
      const double penalty = compute_penalty(degree, extent1, extent1);

      for (unsigned int point = 0; point < n_q_points; ++point)
        {
          for (unsigned int i = 0; i < dofs_per_cell; ++i)
            for (unsigned int j = 0; j < dofs_per_cell; ++j)
              copy_data.cell_matrix(i, j) +=
                (-diffusion_coefficient *        // - nu
                   fe_fv.shape_value(i, point) * // v_h
                   (fe_fv.shape_grad(j, point) * // (grad u_h .
                    normals[point])              //  n)

                 - diffusion_coefficient *         // - nu
                     (fe_fv.shape_grad(i, point) * // (grad v_h .
                      normals[point]) *            //  n)
                     fe_fv.shape_value(j, point)   // u_h

                 + diffusion_coefficient * penalty * // + nu sigma
                     fe_fv.shape_value(i, point) *   // v_h
                     fe_fv.shape_value(j, point)     // u_h

                 ) *
                JxW[point]; // dx

          for (unsigned int i = 0; i < dofs_per_cell; ++i)
            copy_data.cell_rhs(i) +=
              (-diffusion_coefficient *        // - nu
                 (fe_fv.shape_grad(i, point) * // (grad v_h .
                  normals[point]) *            //  n)
                 g[point]                      // g


               + diffusion_coefficient * penalty *        // + nu sigma
                   fe_fv.shape_value(i, point) * g[point] // v_h g

               ) *
              JxW[point]; // dx
        }
    };

    const auto face_worker = [&](const auto         &cell,
                                 const unsigned int &f,
                                 const unsigned int &sf,
                                 const auto         &ncell,
                                 const unsigned int &nf,
                                 const unsigned int &nsf,
                                 auto               &scratch_data,
                                 auto               &copy_data) {
      const FEInterfaceValues<dim> &fe_iv =
        scratch_data.reinit(cell, f, sf, ncell, nf, nsf);

      Assert(fe_iv.get_cell(0) == cell, ExcInternalError());
      Assert(fe_iv.get_face_number(0) == f, ExcInternalError());
      Assert(fe_iv.get_cell(1) == ncell, ExcInternalError());
      Assert(fe_iv.get_face_number(1) == nf, ExcInternalError());

      copy_data.face_data.emplace_back();
      CopyDataFace      &copy_data_face = copy_data.face_data.back();
      const unsigned int n_dofs_face    = fe_iv.n_current_interface_dofs();
      copy_data_face.joint_dof_indices  = fe_iv.get_interface_dof_indices();
      copy_data_face.cell_matrix.reinit(n_dofs_face, n_dofs_face);

      const std::vector<double>         &JxW     = fe_iv.get_JxW_values();
      const std::vector<Tensor<1, dim>> &normals = fe_iv.get_normal_vectors();

      const double extent1 = cell->measure() / cell->face(f)->measure();
      const double extent2 = ncell->measure() / ncell->face(nf)->measure();
      const double penalty = compute_penalty(degree, extent1, extent2);

      for (const unsigned int point : fe_iv.quadrature_point_indices())
        {
          for (const unsigned int i : fe_iv.dof_indices())
            for (const unsigned int j : fe_iv.dof_indices())
              copy_data_face.cell_matrix(i, j) +=
                (-diffusion_coefficient *                 // - nu
                   fe_iv.jump_in_shape_values(i, point) * // [v_h]
                   (fe_iv.average_of_shape_gradients(j,
                                                     point) * // ({grad u_h} .
                    normals[point])                           //  n)

                 -
                 diffusion_coefficient *                         // - nu
                   (fe_iv.average_of_shape_gradients(i, point) * // (grad v_h .
                    normals[point]) *                            //  n)
                   fe_iv.jump_in_shape_values(j, point)          // [u_h]

                 + diffusion_coefficient * penalty *        // + nu sigma
                     fe_iv.jump_in_shape_values(i, point) * // [v_h]
                     fe_iv.jump_in_shape_values(j, point)   // [u_h]

                 ) *
                JxW[point]; // dx
        }
    };

    AffineConstraints<double> constraints;
    constraints.close();
    const auto copier = [&](const auto &c) {
      constraints.distribute_local_to_global(c.cell_matrix,
                                             c.cell_rhs,
                                             c.local_dof_indices,
                                             system_matrix,
                                             system_rhs);

      for (auto &cdf : c.face_data)
        {
          constraints.distribute_local_to_global(cdf.cell_matrix,
                                                 cdf.joint_dof_indices,
                                                 system_matrix);
        }
    };

    UpdateFlags cell_flags = update_values | update_gradients |
                             update_quadrature_points | update_JxW_values;
    UpdateFlags face_flags = update_values | update_gradients |
                             update_quadrature_points | update_normal_vectors |
                             update_JxW_values;

    ScratchData scratch_data(
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
  }

  template <int dim>
  void
  SIPGLaplace<dim>::solve()
  {
    SparseDirectUMFPACK A_direct;
    A_direct.initialize(system_matrix);
    A_direct.vmult(solution, system_rhs);
  }

  template <int dim>
  void
  SIPGLaplace<dim>::output_results(const unsigned int cycle) const
  {
    std::string filename = "sol_Q" + Utilities::int_to_string(degree, 1) + "-" +
                           Utilities::int_to_string(cycle, 2) + ".vtu";
    std::ofstream output(filename);

    DataOut<dim> data_out;
    data_out.attach_dof_handler(dof_handler);
    data_out.add_data_vector(solution, "u", DataOut<dim>::type_dof_data);
    data_out.build_patches(mapping);
    data_out.write_vtu(output);
  }

  template <int dim>
  void
  SIPGLaplace<dim>::compute_error_estimate()
  {
    estimated_error_square_per_cell.reinit(triangulation.n_active_cells());

    const auto cell_worker =
      [&](const auto &cell, auto &scratch_data, auto &copy_data) {
        const FEValues<dim> &fe_v = scratch_data.reinit(cell);

        copy_data.cell_index = cell->active_cell_index();

        const auto                &q_points   = fe_v.get_quadrature_points();
        const unsigned int         n_q_points = q_points.size();
        const std::vector<double> &JxW        = fe_v.get_JxW_values();

        std::vector<Tensor<2, dim>> hessians(n_q_points);
        fe_v.get_function_hessians(solution, hessians);

        std::vector<double> rhs(n_q_points);
        rhs_function->value_list(q_points, rhs);

        const double hk                   = cell->diameter();
        double       residual_norm_square = 0;

        for (unsigned int point = 0; point < n_q_points; ++point)
          {
            const double residual =
              rhs[point] + diffusion_coefficient * trace(hessians[point]);
            residual_norm_square += residual * residual * JxW[point];
          }
        copy_data.value = hk * hk * residual_norm_square;
      };

    const auto boundary_worker = [&](const auto         &cell,
                                     const unsigned int &face_no,
                                     auto               &scratch_data,
                                     auto               &copy_data) {
      const FEFaceValuesBase<dim> &fe_fv = scratch_data.reinit(cell, face_no);

      const auto    &q_points   = fe_fv.get_quadrature_points();
      const unsigned n_q_points = q_points.size();

      const std::vector<double> &JxW = fe_fv.get_JxW_values();

      std::vector<double> g(n_q_points);
      exact_solution->value_list(q_points, g);

      std::vector<double> sol_u(n_q_points);
      fe_fv.get_function_values(solution, sol_u);

      const double extent1 = cell->measure() / cell->face(face_no)->measure();
      const double penalty = compute_penalty(degree, extent1, extent1);

      double difference_norm_square = 0.;
      for (unsigned int point = 0; point < q_points.size(); ++point)
        {
          const double diff = (g[point] - sol_u[point]);
          difference_norm_square += diff * diff * JxW[point];
        }
      copy_data.value += penalty * difference_norm_square;
    };

    const auto face_worker = [&](const auto         &cell,
                                 const unsigned int &f,
                                 const unsigned int &sf,
                                 const auto         &ncell,
                                 const unsigned int &nf,
                                 const unsigned int &nsf,
                                 auto               &scratch_data,
                                 auto               &copy_data) {
      const FEInterfaceValues<dim> &fe_iv =
        scratch_data.reinit(cell, f, sf, ncell, nf, nsf);

      copy_data.face_data.emplace_back();
      CopyDataFace &copy_data_face = copy_data.face_data.back();

      copy_data_face.cell_indices[0] = cell->active_cell_index();
      copy_data_face.cell_indices[1] = ncell->active_cell_index();

      const std::vector<double>         &JxW     = fe_iv.get_JxW_values();
      const std::vector<Tensor<1, dim>> &normals = fe_iv.get_normal_vectors();

      const auto        &q_points   = fe_iv.get_quadrature_points();
      const unsigned int n_q_points = q_points.size();

      std::vector<double> jump(n_q_points);
      fe_iv.get_jump_in_function_values(solution, jump);

      std::vector<Tensor<1, dim>> grad_jump(n_q_points);
      fe_iv.get_jump_in_function_gradients(solution, grad_jump);

      const double h = cell->face(f)->diameter();

      const double extent1 = cell->measure() / cell->face(f)->measure();
      const double extent2 = ncell->measure() / ncell->face(nf)->measure();
      const double penalty = compute_penalty(degree, extent1, extent2);

      double flux_jump_square = 0;
      double u_jump_square    = 0;
      for (unsigned int point = 0; point < n_q_points; ++point)
        {
          u_jump_square += jump[point] * jump[point] * JxW[point];
          const double flux_jump = grad_jump[point] * normals[point];
          flux_jump_square +=
            diffusion_coefficient * flux_jump * flux_jump * JxW[point];
        }
      copy_data_face.values[0] =
        0.5 * h * (flux_jump_square + penalty * u_jump_square);
      copy_data_face.values[1] = copy_data_face.values[0];
    };

    const auto copier = [&](const auto &copy_data) {
      if (copy_data.cell_index != numbers::invalid_unsigned_int)
        estimated_error_square_per_cell[copy_data.cell_index] +=
          copy_data.value;
      for (auto &cdf : copy_data.face_data)
        for (unsigned int j = 0; j < 2; ++j)
          estimated_error_square_per_cell[cdf.cell_indices[j]] += cdf.values[j];
    };

    UpdateFlags cell_flags =
      update_hessians | update_quadrature_points | update_JxW_values;
    UpdateFlags face_flags = update_values | update_gradients |
                             update_quadrature_points | update_JxW_values |
                             update_normal_vectors;

    ScratchData scratch_data(
      mapping, fe, quadrature, cell_flags, face_quadrature, face_flags);

    CopyData cd;
    MeshWorker::mesh_loop(dof_handler.begin_active(),
                          dof_handler.end(),
                          cell_worker,
                          copier,
                          scratch_data,
                          cd,
                          MeshWorker::assemble_own_cells |
                            MeshWorker::assemble_own_interior_faces_once |
                            MeshWorker::assemble_boundary_faces,
                          boundary_worker,
                          face_worker);
  }

  template <int dim>
  double
  SIPGLaplace<dim>::compute_energy_norm()
  {
    energy_norm_square_per_cell.reinit(triangulation.n_active_cells());

    const auto cell_worker =
      [&](const auto &cell, auto &scratch_data, auto &copy_data) {
        const FEValues<dim> &fe_v = scratch_data.reinit(cell);

        copy_data.cell_index = cell->active_cell_index();

        const auto                &q_points   = fe_v.get_quadrature_points();
        const unsigned int         n_q_points = q_points.size();
        const std::vector<double> &JxW        = fe_v.get_JxW_values();

        std::vector<Tensor<1, dim>> grad_u(n_q_points);
        fe_v.get_function_gradients(solution, grad_u);

        std::vector<Tensor<1, dim>> grad_exact(n_q_points);
        exact_solution->gradient_list(q_points, grad_exact);

        double norm_square = 0;
        for (unsigned int point = 0; point < n_q_points; ++point)
          {
            norm_square +=
              (grad_u[point] - grad_exact[point]).norm_square() * JxW[point];
          }
        copy_data.value = norm_square;
      };

    const auto boundary_worker = [&](const auto         &cell,
                                     const unsigned int &face_no,
                                     auto               &scratch_data,
                                     auto               &copy_data) {
      const FEFaceValuesBase<dim> &fe_fv = scratch_data.reinit(cell, face_no);

      const auto    &q_points   = fe_fv.get_quadrature_points();
      const unsigned n_q_points = q_points.size();

      const std::vector<double> &JxW = fe_fv.get_JxW_values();

      std::vector<double> g(n_q_points);
      exact_solution->value_list(q_points, g);

      std::vector<double> sol_u(n_q_points);
      fe_fv.get_function_values(solution, sol_u);

      const double extent1 = cell->measure() / cell->face(face_no)->measure();
      const double penalty = compute_penalty(degree, extent1, extent1);

      double difference_norm_square = 0.;
      for (unsigned int point = 0; point < q_points.size(); ++point)
        {
          const double diff = (g[point] - sol_u[point]);
          difference_norm_square += diff * diff * JxW[point];
        }
      copy_data.value += penalty * difference_norm_square;
    };

    const auto face_worker = [&](const auto         &cell,
                                 const unsigned int &f,
                                 const unsigned int &sf,
                                 const auto         &ncell,
                                 const unsigned int &nf,
                                 const unsigned int &nsf,
                                 auto               &scratch_data,
                                 auto               &copy_data) {
      const FEInterfaceValues<dim> &fe_iv =
        scratch_data.reinit(cell, f, sf, ncell, nf, nsf);

      copy_data.face_data.emplace_back();
      CopyDataFace &copy_data_face = copy_data.face_data.back();

      copy_data_face.cell_indices[0] = cell->active_cell_index();
      copy_data_face.cell_indices[1] = ncell->active_cell_index();

      const std::vector<double> &JxW = fe_iv.get_JxW_values();

      const auto        &q_points   = fe_iv.get_quadrature_points();
      const unsigned int n_q_points = q_points.size();

      std::vector<double> jump(n_q_points);
      fe_iv.get_jump_in_function_values(solution, jump);

      const double extent1 = cell->measure() / cell->face(f)->measure();
      const double extent2 = ncell->measure() / ncell->face(nf)->measure();
      const double penalty = compute_penalty(degree, extent1, extent2);

      double u_jump_square = 0;
      for (unsigned int point = 0; point < n_q_points; ++point)
        {
          u_jump_square += jump[point] * jump[point] * JxW[point];
        }
      copy_data_face.values[0] = 0.5 * penalty * u_jump_square;
      copy_data_face.values[1] = copy_data_face.values[0];
    };

    const auto copier = [&](const auto &copy_data) {
      if (copy_data.cell_index != numbers::invalid_unsigned_int)
        energy_norm_square_per_cell[copy_data.cell_index] += copy_data.value;
      for (auto &cdf : copy_data.face_data)
        for (unsigned int j = 0; j < 2; ++j)
          energy_norm_square_per_cell[cdf.cell_indices[j]] += cdf.values[j];
    };

    UpdateFlags cell_flags =
      update_gradients | update_quadrature_points | update_JxW_values;
    UpdateFlags face_flags =
      update_values | update_quadrature_points | update_JxW_values;

    ScratchData scratch_data(
      mapping, fe, quadrature_2, cell_flags, face_quadrature_2, face_flags);

    CopyData cd;
    MeshWorker::mesh_loop(dof_handler.begin_active(),
                          dof_handler.end(),
                          cell_worker,
                          copier,
                          scratch_data,
                          cd,
                          MeshWorker::assemble_own_cells |
                            MeshWorker::assemble_own_interior_faces_once |
                            MeshWorker::assemble_boundary_faces,
                          boundary_worker,
                          face_worker);
    const double energy_error =
      std::sqrt(energy_norm_square_per_cell.l1_norm());
    return energy_error;
  }

  template <int dim>
  void
  SIPGLaplace<dim>::refine_grid()
  {
    const double refinement_fraction = 0.1;

    GridRefinement::refine_and_coarsen_fixed_number(
      triangulation, estimated_error_square_per_cell, refinement_fraction, 0.);

    triangulation.execute_coarsening_and_refinement();
  }

  template <int dim>
  void
  SIPGLaplace<dim>::compute_errors()
  {
    double L2_error, H1_error;

    {
      Vector<float> difference_per_cell(triangulation.n_active_cells());
      VectorTools::integrate_difference(mapping,
                                        dof_handler,
                                        solution,
                                        *(exact_solution.get()),
                                        difference_per_cell,
                                        quadrature_2,
                                        VectorTools::L2_norm);

      L2_error = VectorTools::compute_global_error(triangulation,
                                                   difference_per_cell,
                                                   VectorTools::L2_norm);
    }

    {
      Vector<float> difference_per_cell(triangulation.n_active_cells());
      VectorTools::integrate_difference(mapping,
                                        dof_handler,
                                        solution,
                                        *(exact_solution.get()),
                                        difference_per_cell,
                                        quadrature_2,
                                        VectorTools::H1_seminorm);

      H1_error = VectorTools::compute_global_error(triangulation,
                                                   difference_per_cell,
                                                   VectorTools::H1_seminorm);
    }

    convergence_table.add_value("L2", L2_error);
    convergence_table.add_value("H1", H1_error);
    const double energy_error = compute_energy_norm();
    convergence_table.add_value("Energy", energy_error);

    deallog << "  Error in the L2 norm         : " << L2_error << std::endl
            << "  Error in the H1 seminorm     : " << H1_error << std::endl
            << "  Error in the energy norm     : " << energy_error << std::endl;
  }

  template <int dim>
  void
  SIPGLaplace<dim>::run()
  {
    unsigned int max_cycle = 3;
    for (unsigned int cycle = 0; cycle < max_cycle; ++cycle)
      {
        deallog << "Cycle " << cycle << std::endl;

        switch (test_case)
          {
            case Test_Case::convergence_rate:
              {
                if (cycle == 0)
                  {
#ifdef USE_SIMPLEX
                    Triangulation<dim> temp;
                    GridGenerator::hyper_cube(temp);
                    GridGenerator::convert_hypercube_to_simplex_mesh(
                      temp, triangulation);
#else
                    GridGenerator::hyper_cube(triangulation);
#endif
                    triangulation.refine_global(2);
                  }
                else
                  {
                    triangulation.refine_global();
                  }
                break;
              }
            case Test_Case::l_singularity:
              {
                if (cycle == 0)
                  {
#ifdef USE_SIMPLEX
                    Triangulation<dim> temp;
                    GridGenerator::hyper_L(temp);
                    GridGenerator::convert_hypercube_to_simplex_mesh(
                      temp, triangulation);
#else
                    GridGenerator::hyper_L(triangulation);
#endif
                    triangulation.refine_global(3);
                  }
                else
                  {
                    triangulation.refine_global();
                  }
                break;
              }
            default:
              {
                DEAL_II_NOT_IMPLEMENTED();
              }
          }
        deallog << "  Number of active cells       : "
                << triangulation.n_active_cells() << std::endl;
        setup_system();

        deallog << "  Number of degrees of freedom : " << dof_handler.n_dofs()
                << std::endl;

        assemble_system();
        solve();
        output_results(cycle);
        {
          convergence_table.add_value("cycle", cycle);
          convergence_table.add_value("cells", triangulation.n_active_cells());
          convergence_table.add_value("dofs", dof_handler.n_dofs());
        }
        compute_errors();

        if (test_case == Test_Case::l_singularity)
          {
            compute_error_estimate();
            deallog << "  Estimated error              : "
                    << std::sqrt(estimated_error_square_per_cell.l1_norm())
                    << std::endl;

            convergence_table.add_value(
              "Estimator",
              std::sqrt(estimated_error_square_per_cell.l1_norm()));
          }
        deallog << std::endl;
      }
    {
      convergence_table.set_precision("L2", 3);
      convergence_table.set_precision("H1", 3);
      convergence_table.set_precision("Energy", 3);

      convergence_table.set_scientific("L2", true);
      convergence_table.set_scientific("H1", true);
      convergence_table.set_scientific("Energy", true);

      if (test_case == Test_Case::l_singularity)
        {
          convergence_table.set_precision("Estimator", 3);
          convergence_table.set_scientific("Estimator", true);
        }
      if (test_case == Test_Case::convergence_rate)
        {
          convergence_table.evaluate_convergence_rates(
            "L2", ConvergenceTable::reduction_rate_log2);
          convergence_table.evaluate_convergence_rates(
            "H1", ConvergenceTable::reduction_rate_log2);
        }
      deallog << "degree = " << degree << std::endl;
      convergence_table.write_text(
        deallog.get_file_stream(),
        TableHandler::TextOutputFormat::org_mode_table);
    }
  }
} // namespace Step74


int
main()
{
  initlog();

  try
    {
      using namespace Step74;
      Test_Case      test_case = Test_Case::l_singularity;
      SIPGLaplace<2> problem(test_case);
      problem.run();
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
    };

  return 0;
}
