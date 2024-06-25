/* ------------------------------------------------------------------------
 *
 * SPDX-License-Identifier: LGPL-2.1-or-later
 * Copyright (C) 2020 - 2024 by the deal.II authors
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
 * Authors: Timo Heister and Jiaqi Zhang, Clemson University, 2020
 */

// The first few files have already been covered in previous examples and will
// thus not be further commented on:
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/function_lib.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/sparse_direct.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q1.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/numerics/data_out.h>
// Here the discontinuous finite elements and FEInterfaceValues are defined.
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_interface_values.h>

#include <deal.II/numerics/derivative_approximation.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/base/convergence_table.h>

#include <deal.II/meshworker/copy_data.h>
#include <deal.II/meshworker/mesh_loop.h>
#include <deal.II/meshworker/scratch_data.h>

namespace Step74
{
  using namespace dealii;

  // @sect3{Equation data}
  // Here we define two test cases: convergence_rate for a smooth function
  // and l_singularity for the Functions::LSingularityFunction.
  enum class TestCase
  {
    convergence_rate,
    l_singularity
  };



  // A smooth solution for the convergence test:
  template <int dim>
  class SmoothSolution : public Function<dim>
  {
  public:
    SmoothSolution()
      : Function<dim>()
    {}

    virtual void value_list(const std::vector<Point<dim>> &points,
                            std::vector<double>           &values,
                            const unsigned int component = 0) const override;

    virtual Tensor<1, dim>
    gradient(const Point<dim>  &point,
             const unsigned int component = 0) const override;
  };



  template <int dim>
  void SmoothSolution<dim>::value_list(const std::vector<Point<dim>> &points,
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



  // The corresponding right-hand side of the smooth function:
  template <int dim>
  class SmoothRightHandSide : public Function<dim>
  {
  public:
    SmoothRightHandSide()
      : Function<dim>()
    {}

    virtual void value_list(const std::vector<Point<dim>> &points,
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



  // The right-hand side that corresponds to the function
  // Functions::LSingularityFunction, where we
  // assume that the diffusion coefficient $\nu = 1$:
  template <int dim>
  class SingularRightHandSide : public Function<dim>
  {
  public:
    SingularRightHandSide()
      : Function<dim>()
    {}

    virtual void value_list(const std::vector<Point<dim>> &points,
                            std::vector<double>           &values,
                            const unsigned int /*component*/) const override;

  private:
    const Functions::LSingularityFunction ref;
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



  // @sect3{Auxiliary functions}
  // This function computes the penalty $\sigma$.
  double get_penalty_factor(const unsigned int fe_degree,
                            const double       cell_extent_left,
                            const double       cell_extent_right)
  {
    const unsigned int degree = std::max(1U, fe_degree);
    return degree * (degree + 1.) * 0.5 *
           (1. / cell_extent_left + 1. / cell_extent_right);
  }


  // @sect3{The CopyData}
  // In the following, we define "Copy" objects for the MeshWorker::mesh_loop(),
  // which is essentially the same as step-12. Note that the
  // "Scratch" object is not defined here because we use
  // MeshWorker::ScratchData<dim> instead. (The use of "Copy" and "Scratch"
  // objects is extensively explained in the WorkStream namespace documentation.
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
    void reinit(const Iterator &cell, const unsigned int dofs_per_cell)
    {
      cell_matrix.reinit(dofs_per_cell, dofs_per_cell);
      cell_rhs.reinit(dofs_per_cell);
      local_dof_indices.resize(dofs_per_cell);
      cell->get_dof_indices(local_dof_indices);
    }
  };



  // @sect3{The SIPGLaplace class}
  // After these preparations, we proceed with the main class of this program,
  // called `SIPGLaplace`. The overall structure of the class is as in many
  // of the other tutorial programs. Major differences will only come up in the
  // implementation of the assemble functions, since we use FEInterfaceValues to
  // assemble face terms.
  template <int dim>
  class SIPGLaplace
  {
  public:
    SIPGLaplace(const TestCase &test_case);
    void run();

  private:
    void setup_system();
    void assemble_system();
    void solve();
    void refine_grid();
    void output_results(const unsigned int cycle) const;

    void   compute_errors();
    void   compute_error_estimate();
    double compute_energy_norm_error();

    Triangulation<dim>    triangulation;
    const unsigned int    degree;
    const QGauss<dim>     quadrature;
    const QGauss<dim - 1> face_quadrature;
    const QGauss<dim>     quadrature_overintegration;
    const QGauss<dim - 1> face_quadrature_overintegration;
    const MappingQ1<dim>  mapping;

    using ScratchData = MeshWorker::ScratchData<dim>;

    const FE_DGQ<dim> fe;
    DoFHandler<dim>   dof_handler;

    SparsityPattern      sparsity_pattern;
    SparseMatrix<double> system_matrix;
    Vector<double>       solution;
    Vector<double>       system_rhs;

    // The remainder of the class's members are used for the following:
    // - Vectors to store error estimator square and energy norm square per
    // cell.
    // - Print convergence rate and errors on the screen.
    // - The fiffusion coefficient $\nu$ is set to 1.
    // - Members that store information about the test case to be computed.
    Vector<double> estimated_error_square_per_cell;
    Vector<double> energy_norm_square_per_cell;

    ConvergenceTable convergence_table;

    const double diffusion_coefficient = 1.;

    const TestCase                       test_case;
    std::unique_ptr<const Function<dim>> exact_solution;
    std::unique_ptr<const Function<dim>> rhs_function;
  };

  // The constructor here takes the test case as input and then
  // determines the correct solution and right-hand side classes. The
  // remaining member variables are initialized in the obvious way.
  template <int dim>
  SIPGLaplace<dim>::SIPGLaplace(const TestCase &test_case)
    : degree(3)
    , quadrature(degree + 1)
    , face_quadrature(degree + 1)
    , quadrature_overintegration(degree + 2)
    , face_quadrature_overintegration(degree + 2)
    , mapping()
    , fe(degree)
    , dof_handler(triangulation)
    , test_case(test_case)
  {
    if (test_case == TestCase::convergence_rate)
      {
        exact_solution = std::make_unique<const SmoothSolution<dim>>();
        rhs_function   = std::make_unique<const SmoothRightHandSide<dim>>();
      }

    else if (test_case == TestCase::l_singularity)
      {
        exact_solution =
          std::make_unique<const Functions::LSingularityFunction>();
        rhs_function = std::make_unique<const SingularRightHandSide<dim>>();
      }
    else
      AssertThrow(false, ExcNotImplemented());
  }



  template <int dim>
  void SIPGLaplace<dim>::setup_system()
  {
    dof_handler.distribute_dofs(fe);
    DynamicSparsityPattern dsp(dof_handler.n_dofs());
    DoFTools::make_flux_sparsity_pattern(dof_handler, dsp);
    sparsity_pattern.copy_from(dsp);

    system_matrix.reinit(sparsity_pattern);
    solution.reinit(dof_handler.n_dofs());
    system_rhs.reinit(dof_handler.n_dofs());
  }



  // @sect3{The assemble_system function}
  // The assemble function here is similar to that in step-12 and step-47.
  // Different from assembling by hand, we just need to focus
  // on assembling on each cell, each boundary face, and each
  // interior face. The loops over cells and faces are handled
  // automatically by MeshWorker::mesh_loop().
  //
  // The function starts by defining a local (lambda) function that is
  // used to integrate the cell terms:
  template <int dim>
  void SIPGLaplace<dim>::assemble_system()
  {
    const auto cell_worker =
      [&](const typename DoFHandler<dim>::active_cell_iterator &cell,
          ScratchData                                          &scratch_data,
          CopyData                                             &copy_data) {
        const FEValues<dim> &fe_v          = scratch_data.reinit(cell);
        const unsigned int   dofs_per_cell = fe_v.dofs_per_cell;
        copy_data.reinit(cell, dofs_per_cell);

        const std::vector<Point<dim>> &q_points =
          scratch_data.get_quadrature_points();
        const unsigned int         n_q_points = q_points.size();
        const std::vector<double> &JxW        = scratch_data.get_JxW_values();

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

    // Next, we need a function that assembles face integrals on the boundary:
    const auto boundary_worker =
      [&](const typename DoFHandler<dim>::active_cell_iterator &cell,
          const unsigned int                                   &face_no,
          ScratchData                                          &scratch_data,
          CopyData                                             &copy_data) {
        const FEFaceValuesBase<dim> &fe_fv = scratch_data.reinit(cell, face_no);

        const std::vector<Point<dim>> &q_points =
          scratch_data.get_quadrature_points();
        const unsigned int n_q_points    = q_points.size();
        const unsigned int dofs_per_cell = fe_fv.dofs_per_cell;

        const std::vector<double>         &JxW = scratch_data.get_JxW_values();
        const std::vector<Tensor<1, dim>> &normals =
          scratch_data.get_normal_vectors();

        std::vector<double> g(n_q_points);
        exact_solution->value_list(q_points, g);

        const double extent1 = cell->measure() / cell->face(face_no)->measure();
        const double penalty = get_penalty_factor(degree, extent1, extent1);

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

    // Finally, a function that assembles face integrals on interior
    // faces. To reinitialize FEInterfaceValues, we need to pass
    // cells, face and subface indices (for adaptive refinement) to
    // the reinit() function of FEInterfaceValues:
    const auto face_worker =
      [&](const typename DoFHandler<dim>::cell_iterator &cell,
          const unsigned int                            &f,
          const unsigned int                            &sf,
          const typename DoFHandler<dim>::cell_iterator &ncell,
          const unsigned int                            &nf,
          const unsigned int                            &nsf,
          ScratchData                                   &scratch_data,
          CopyData                                      &copy_data) {
        const FEInterfaceValues<dim> &fe_iv =
          scratch_data.reinit(cell, f, sf, ncell, nf, nsf);

        copy_data.face_data.emplace_back();
        CopyDataFace      &copy_data_face = copy_data.face_data.back();
        const unsigned int n_dofs_face    = fe_iv.n_current_interface_dofs();
        copy_data_face.joint_dof_indices  = fe_iv.get_interface_dof_indices();
        copy_data_face.cell_matrix.reinit(n_dofs_face, n_dofs_face);

        const std::vector<double>         &JxW     = fe_iv.get_JxW_values();
        const std::vector<Tensor<1, dim>> &normals = fe_iv.get_normal_vectors();

        const double extent1 = cell->measure() / cell->face(f)->measure();
        const double extent2 = ncell->measure() / ncell->face(nf)->measure();
        const double penalty = get_penalty_factor(degree, extent1, extent2);

        for (const unsigned int point : fe_iv.quadrature_point_indices())
          {
            for (const unsigned int i : fe_iv.dof_indices())
              for (const unsigned int j : fe_iv.dof_indices())
                copy_data_face.cell_matrix(i, j) +=
                  (-diffusion_coefficient *                     // - nu
                     fe_iv.jump_in_shape_values(i, point) *     // [v_h]
                     (fe_iv.average_of_shape_gradients(j,       //
                                                       point) * // ({grad u_h} .
                      normals[point])                           //  n)

                   - diffusion_coefficient * // - nu
                       (fe_iv.average_of_shape_gradients(i,
                                                         point) * // (grad v_h .
                        normals[point]) *                         //  n)
                       fe_iv.jump_in_shape_values(j, point)       // [u_h]

                   + diffusion_coefficient * penalty *        // + nu sigma
                       fe_iv.jump_in_shape_values(i, point) * // [v_h]
                       fe_iv.jump_in_shape_values(j, point)   // [u_h]

                   ) *
                  JxW[point]; // dx
          }
      };

    // The following lambda function will then copy data into the
    // global matrix and right-hand side.  Though there are no hanging
    // node constraints in DG discretization, we define an empty
    // AffineConstraints object that allows us to use the
    // AffineConstraints::distribute_local_to_global() functionality.
    AffineConstraints<double> constraints;
    constraints.close();
    const auto copier = [&](const CopyData &c) {
      constraints.distribute_local_to_global(c.cell_matrix,
                                             c.cell_rhs,
                                             c.local_dof_indices,
                                             system_matrix,
                                             system_rhs);

      // Copy data from interior face assembly to the global matrix.
      for (const CopyDataFace &cdf : c.face_data)
        {
          constraints.distribute_local_to_global(cdf.cell_matrix,
                                                 cdf.joint_dof_indices,
                                                 system_matrix);
        }
    };


    // With the assembly functions defined, we can now create
    // ScratchData and CopyData objects, and pass them together with
    // the lambda functions above to MeshWorker::mesh_loop(). In
    // addition, we need to specify that we want to assemble on
    // interior faces exactly once.
    const UpdateFlags cell_flags = update_values | update_gradients |
                                   update_quadrature_points | update_JxW_values;
    const UpdateFlags face_flags = update_values | update_gradients |
                                   update_quadrature_points |
                                   update_normal_vectors | update_JxW_values;

    ScratchData scratch_data(
      mapping, fe, quadrature, cell_flags, face_quadrature, face_flags);
    CopyData copy_data;

    MeshWorker::mesh_loop(dof_handler.begin_active(),
                          dof_handler.end(),
                          cell_worker,
                          copier,
                          scratch_data,
                          copy_data,
                          MeshWorker::assemble_own_cells |
                            MeshWorker::assemble_boundary_faces |
                            MeshWorker::assemble_own_interior_faces_once,
                          boundary_worker,
                          face_worker);
  }



  // @sect3{The solve() and output_results() function}
  // The following two functions are entirely standard and without difficulty.
  template <int dim>
  void SIPGLaplace<dim>::solve()
  {
    SparseDirectUMFPACK A_direct;
    A_direct.initialize(system_matrix);
    A_direct.vmult(solution, system_rhs);
  }



  template <int dim>
  void SIPGLaplace<dim>::output_results(const unsigned int cycle) const
  {
    const std::string filename = "sol_Q" + Utilities::int_to_string(degree, 1) +
                                 "-" + Utilities::int_to_string(cycle, 2) +
                                 ".vtu";
    std::ofstream output(filename);

    DataOut<dim> data_out;
    data_out.attach_dof_handler(dof_handler);
    data_out.add_data_vector(solution, "u", DataOut<dim>::type_dof_data);
    data_out.build_patches(mapping);
    data_out.write_vtu(output);
  }


  // @sect3{The compute_error_estimate() function}
  // The assembly of the error estimator here is quite similar to
  // that of the global matrix and right-had side and can be handled
  // by the MeshWorker::mesh_loop() framework. To understand what
  // each of the local (lambda) functions is doing, recall first that
  // the local cell residual is defined as
  // $h_K^2 \left\| f + \nu \Delta u_h \right\|_K^2$:
  template <int dim>
  void SIPGLaplace<dim>::compute_error_estimate()
  {
    const auto cell_worker =
      [&](const typename DoFHandler<dim>::active_cell_iterator &cell,
          ScratchData                                          &scratch_data,
          CopyData                                             &copy_data) {
        const FEValues<dim> &fe_v = scratch_data.reinit(cell);

        copy_data.cell_index = cell->active_cell_index();

        const std::vector<Point<dim>> &q_points = fe_v.get_quadrature_points();
        const unsigned int             n_q_points = q_points.size();
        const std::vector<double>     &JxW        = fe_v.get_JxW_values();

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

    // Next compute boundary terms $\sum_{f\in \partial K \cap \partial \Omega}
    // \sigma \left\| [  u_h-g_D ]  \right\|_f^2  $:
    const auto boundary_worker =
      [&](const typename DoFHandler<dim>::active_cell_iterator &cell,
          const unsigned int                                   &face_no,
          ScratchData                                          &scratch_data,
          CopyData                                             &copy_data) {
        const FEFaceValuesBase<dim> &fe_fv = scratch_data.reinit(cell, face_no);

        const std::vector<Point<dim>> &q_points = fe_fv.get_quadrature_points();
        const unsigned                 n_q_points = q_points.size();

        const std::vector<double> &JxW = fe_fv.get_JxW_values();

        std::vector<double> g(n_q_points);
        exact_solution->value_list(q_points, g);

        std::vector<double> sol_u(n_q_points);
        fe_fv.get_function_values(solution, sol_u);

        const double extent1 = cell->measure() / cell->face(face_no)->measure();
        const double penalty = get_penalty_factor(degree, extent1, extent1);

        double difference_norm_square = 0.;
        for (unsigned int point = 0; point < q_points.size(); ++point)
          {
            const double diff = (g[point] - sol_u[point]);
            difference_norm_square += diff * diff * JxW[point];
          }
        copy_data.value += penalty * difference_norm_square;
      };

    // And finally interior face terms $\sum_{f\in \partial K}\lbrace \sigma
    // \left\| [u_h]  \right\|_f^2   +  h_f \left\|  [\nu \nabla u_h \cdot
    // \mathbf n ] \right\|_f^2 \rbrace$:
    const auto face_worker =
      [&](const typename DoFHandler<dim>::cell_iterator &cell,
          const unsigned int                            &f,
          const unsigned int                            &sf,
          const typename DoFHandler<dim>::cell_iterator &ncell,
          const unsigned int                            &nf,
          const unsigned int                            &nsf,
          ScratchData                                   &scratch_data,
          CopyData                                      &copy_data) {
        const FEInterfaceValues<dim> &fe_iv =
          scratch_data.reinit(cell, f, sf, ncell, nf, nsf);

        copy_data.face_data.emplace_back();
        CopyDataFace &copy_data_face = copy_data.face_data.back();

        copy_data_face.cell_indices[0] = cell->active_cell_index();
        copy_data_face.cell_indices[1] = ncell->active_cell_index();

        const std::vector<double>         &JxW     = fe_iv.get_JxW_values();
        const std::vector<Tensor<1, dim>> &normals = fe_iv.get_normal_vectors();

        const std::vector<Point<dim>> &q_points = fe_iv.get_quadrature_points();
        const unsigned int             n_q_points = q_points.size();

        std::vector<double> jump(n_q_points);
        fe_iv.get_jump_in_function_values(solution, jump);

        std::vector<Tensor<1, dim>> grad_jump(n_q_points);
        fe_iv.get_jump_in_function_gradients(solution, grad_jump);

        const double h = cell->face(f)->diameter();

        const double extent1 = cell->measure() / cell->face(f)->measure();
        const double extent2 = ncell->measure() / ncell->face(nf)->measure();
        const double penalty = get_penalty_factor(degree, extent1, extent2);

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

    // Having computed local contributions for each cell, we still
    // need a way to copy these into the global vector that will hold
    // the error estimators for all cells:
    const auto copier = [&](const CopyData &copy_data) {
      if (copy_data.cell_index != numbers::invalid_unsigned_int)
        estimated_error_square_per_cell[copy_data.cell_index] +=
          copy_data.value;
      for (const CopyDataFace &cdf : copy_data.face_data)
        for (unsigned int j = 0; j < 2; ++j)
          estimated_error_square_per_cell[cdf.cell_indices[j]] += cdf.values[j];
    };

    // After all of this set-up, let's do the actual work: We resize
    // the vector into which the results will be written, and then
    // drive the whole process using the MeshWorker::mesh_loop()
    // function.
    estimated_error_square_per_cell.reinit(triangulation.n_active_cells());

    const UpdateFlags cell_flags =
      update_hessians | update_quadrature_points | update_JxW_values;
    const UpdateFlags face_flags = update_values | update_gradients |
                                   update_quadrature_points |
                                   update_JxW_values | update_normal_vectors;

    ScratchData scratch_data(
      mapping, fe, quadrature, cell_flags, face_quadrature, face_flags);

    CopyData copy_data;
    MeshWorker::mesh_loop(dof_handler.begin_active(),
                          dof_handler.end(),
                          cell_worker,
                          copier,
                          scratch_data,
                          copy_data,
                          MeshWorker::assemble_own_cells |
                            MeshWorker::assemble_own_interior_faces_once |
                            MeshWorker::assemble_boundary_faces,
                          boundary_worker,
                          face_worker);
  }

  // @sect3{The compute_energy_norm_error() function}
  // Next, we evaluate the accuracy in terms of the energy norm.
  // This function is similar to the assembling of the error estimator above.
  // Here we compute the square of the energy norm defined by
  // @f[
  //   \|u \|_{1,h}^2 = \sum_{K \in \Gamma_h} \nu\|\nabla u \|_K^2 +
  //   \sum_{f \in F_i} \sigma \| [ u ] \|_f^2 +
  //   \sum_{f \in F_b} \sigma  \|u\|_f^2.
  // @f]
  // Therefore the corresponding error is
  // @f[
  //   \|u -u_h \|_{1,h}^2 = \sum_{K \in \Gamma_h} \nu\|\nabla (u_h - u)  \|_K^2
  //   + \sum_{f \in F_i} \sigma  \|[ u_h ] \|_f^2 + \sum_{f \in F_b}\sigma
  //   \|u_h-g_D\|_f^2.
  // @f]
  template <int dim>
  double SIPGLaplace<dim>::compute_energy_norm_error()
  {
    energy_norm_square_per_cell.reinit(triangulation.n_active_cells());

    // Assemble $\sum_{K \in \Gamma_h} \nu\|\nabla (u_h - u)  \|_K^2 $.
    const auto cell_worker =
      [&](const typename DoFHandler<dim>::active_cell_iterator &cell,
          ScratchData                                          &scratch_data,
          CopyData                                             &copy_data) {
        const FEValues<dim> &fe_v = scratch_data.reinit(cell);

        copy_data.cell_index = cell->active_cell_index();

        const std::vector<Point<dim>> &q_points = fe_v.get_quadrature_points();
        const unsigned int             n_q_points = q_points.size();
        const std::vector<double>     &JxW        = fe_v.get_JxW_values();

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
        copy_data.value = diffusion_coefficient * norm_square;
      };

    // Assemble $\sum_{f \in F_b}\sigma  \|u_h-g_D\|_f^2$.
    const auto boundary_worker =
      [&](const typename DoFHandler<dim>::active_cell_iterator &cell,
          const unsigned int                                   &face_no,
          ScratchData                                          &scratch_data,
          CopyData                                             &copy_data) {
        const FEFaceValuesBase<dim> &fe_fv = scratch_data.reinit(cell, face_no);

        const std::vector<Point<dim>> &q_points = fe_fv.get_quadrature_points();
        const unsigned                 n_q_points = q_points.size();

        const std::vector<double> &JxW = fe_fv.get_JxW_values();

        std::vector<double> g(n_q_points);
        exact_solution->value_list(q_points, g);

        std::vector<double> sol_u(n_q_points);
        fe_fv.get_function_values(solution, sol_u);

        const double extent1 = cell->measure() / cell->face(face_no)->measure();
        const double penalty = get_penalty_factor(degree, extent1, extent1);

        double difference_norm_square = 0.;
        for (unsigned int point = 0; point < q_points.size(); ++point)
          {
            const double diff = (g[point] - sol_u[point]);
            difference_norm_square += diff * diff * JxW[point];
          }
        copy_data.value += penalty * difference_norm_square;
      };

    // Assemble $\sum_{f \in F_i} \sigma  \| [ u_h ] \|_f^2$.
    const auto face_worker =
      [&](const typename DoFHandler<dim>::cell_iterator &cell,
          const unsigned int                            &f,
          const unsigned int                            &sf,
          const typename DoFHandler<dim>::cell_iterator &ncell,
          const unsigned int                            &nf,
          const unsigned int                            &nsf,
          ScratchData                                   &scratch_data,
          CopyData                                      &copy_data) {
        const FEInterfaceValues<dim> &fe_iv =
          scratch_data.reinit(cell, f, sf, ncell, nf, nsf);

        copy_data.face_data.emplace_back();
        CopyDataFace &copy_data_face = copy_data.face_data.back();

        copy_data_face.cell_indices[0] = cell->active_cell_index();
        copy_data_face.cell_indices[1] = ncell->active_cell_index();

        const std::vector<double> &JxW = fe_iv.get_JxW_values();

        const std::vector<Point<dim>> &q_points = fe_iv.get_quadrature_points();
        const unsigned int             n_q_points = q_points.size();

        std::vector<double> jump(n_q_points);
        fe_iv.get_jump_in_function_values(solution, jump);

        const double extent1 = cell->measure() / cell->face(f)->measure();
        const double extent2 = ncell->measure() / ncell->face(nf)->measure();
        const double penalty = get_penalty_factor(degree, extent1, extent2);

        double u_jump_square = 0;
        for (unsigned int point = 0; point < n_q_points; ++point)
          {
            u_jump_square += jump[point] * jump[point] * JxW[point];
          }
        copy_data_face.values[0] = 0.5 * penalty * u_jump_square;
        copy_data_face.values[1] = copy_data_face.values[0];
      };

    const auto copier = [&](const CopyData &copy_data) {
      if (copy_data.cell_index != numbers::invalid_unsigned_int)
        energy_norm_square_per_cell[copy_data.cell_index] += copy_data.value;
      for (const CopyDataFace &cdf : copy_data.face_data)
        for (unsigned int j = 0; j < 2; ++j)
          energy_norm_square_per_cell[cdf.cell_indices[j]] += cdf.values[j];
    };

    const UpdateFlags cell_flags =
      update_gradients | update_quadrature_points | update_JxW_values;
    UpdateFlags face_flags =
      update_values | update_quadrature_points | update_JxW_values;

    const ScratchData scratch_data(mapping,
                                   fe,
                                   quadrature_overintegration,
                                   cell_flags,
                                   face_quadrature_overintegration,
                                   face_flags);

    CopyData copy_data;
    MeshWorker::mesh_loop(dof_handler.begin_active(),
                          dof_handler.end(),
                          cell_worker,
                          copier,
                          scratch_data,
                          copy_data,
                          MeshWorker::assemble_own_cells |
                            MeshWorker::assemble_own_interior_faces_once |
                            MeshWorker::assemble_boundary_faces,
                          boundary_worker,
                          face_worker);
    const double energy_error =
      std::sqrt(energy_norm_square_per_cell.l1_norm());
    return energy_error;
  }



  // @sect3{The refine_grid() function}
  template <int dim>
  void SIPGLaplace<dim>::refine_grid()
  {
    const double refinement_fraction = 0.1;

    GridRefinement::refine_and_coarsen_fixed_number(
      triangulation, estimated_error_square_per_cell, refinement_fraction, 0.);

    triangulation.execute_coarsening_and_refinement();
  }



  // @sect3{The compute_errors() function}
  // We compute three errors in the $L_2$ norm, $H_1$ seminorm, and
  // the energy norm, respectively. These are then printed to screen,
  // but also stored in a table that records how these errors decay
  // with mesh refinement and which can be output in one step at the
  // end of the program.
  template <int dim>
  void SIPGLaplace<dim>::compute_errors()
  {
    double L2_error, H1_error, energy_error;

    {
      Vector<float> difference_per_cell(triangulation.n_active_cells());
      VectorTools::integrate_difference(mapping,
                                        dof_handler,
                                        solution,
                                        *(exact_solution.get()),
                                        difference_per_cell,
                                        quadrature_overintegration,
                                        VectorTools::L2_norm);

      L2_error = VectorTools::compute_global_error(triangulation,
                                                   difference_per_cell,
                                                   VectorTools::L2_norm);
      convergence_table.add_value("L2", L2_error);
    }

    {
      Vector<float> difference_per_cell(triangulation.n_active_cells());
      VectorTools::integrate_difference(mapping,
                                        dof_handler,
                                        solution,
                                        *(exact_solution.get()),
                                        difference_per_cell,
                                        quadrature_overintegration,
                                        VectorTools::H1_seminorm);

      H1_error = VectorTools::compute_global_error(triangulation,
                                                   difference_per_cell,
                                                   VectorTools::H1_seminorm);
      convergence_table.add_value("H1", H1_error);
    }

    {
      energy_error = compute_energy_norm_error();
      convergence_table.add_value("Energy", energy_error);
    }

    std::cout << "  Error in the L2 norm         : " << L2_error << std::endl
              << "  Error in the H1 seminorm     : " << H1_error << std::endl
              << "  Error in the energy norm     : " << energy_error
              << std::endl;
  }



  // @sect3{The run() function}
  template <int dim>
  void SIPGLaplace<dim>::run()
  {
    const unsigned int max_cycle =
      (test_case == TestCase::convergence_rate ? 6 : 20);
    for (unsigned int cycle = 0; cycle < max_cycle; ++cycle)
      {
        std::cout << "Cycle " << cycle << std::endl;

        switch (test_case)
          {
            case TestCase::convergence_rate:
              {
                if (cycle == 0)
                  {
                    GridGenerator::hyper_cube(triangulation);

                    triangulation.refine_global(2);
                  }
                else
                  {
                    triangulation.refine_global(1);
                  }
                break;
              }

            case TestCase::l_singularity:
              {
                if (cycle == 0)
                  {
                    GridGenerator::hyper_L(triangulation);
                    triangulation.refine_global(3);
                  }
                else
                  {
                    refine_grid();
                  }
                break;
              }

            default:
              {
                DEAL_II_NOT_IMPLEMENTED();
              }
          }

        std::cout << "  Number of active cells       : "
                  << triangulation.n_active_cells() << std::endl;
        setup_system();

        std::cout << "  Number of degrees of freedom : " << dof_handler.n_dofs()
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

        if (test_case == TestCase::l_singularity)
          {
            compute_error_estimate();
            std::cout << "  Estimated error              : "
                      << std::sqrt(estimated_error_square_per_cell.l1_norm())
                      << std::endl;

            convergence_table.add_value(
              "Estimator",
              std::sqrt(estimated_error_square_per_cell.l1_norm()));
          }
        std::cout << std::endl;
      }

    // Having run all of our computations, let us tell the convergence
    // table how to format its data and output it to screen:
    convergence_table.set_precision("L2", 3);
    convergence_table.set_precision("H1", 3);
    convergence_table.set_precision("Energy", 3);

    convergence_table.set_scientific("L2", true);
    convergence_table.set_scientific("H1", true);
    convergence_table.set_scientific("Energy", true);

    if (test_case == TestCase::convergence_rate)
      {
        convergence_table.evaluate_convergence_rates(
          "L2", ConvergenceTable::reduction_rate_log2);
        convergence_table.evaluate_convergence_rates(
          "H1", ConvergenceTable::reduction_rate_log2);
      }
    if (test_case == TestCase::l_singularity)
      {
        convergence_table.set_precision("Estimator", 3);
        convergence_table.set_scientific("Estimator", true);
      }

    std::cout << "degree = " << degree << std::endl;
    convergence_table.write_text(
      std::cout, TableHandler::TextOutputFormat::org_mode_table);
  }
} // namespace Step74



// @sect3{The main() function}
// The following <code>main</code> function is similar to previous examples as
// well, and need not be commented on.
int main()
{
  try
    {
      using namespace dealii;
      using namespace Step74;

      const TestCase test_case = TestCase::l_singularity;

      SIPGLaplace<2> problem(test_case);
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
    };

  return 0;
}
