/* ------------------------------------------------------------------------
 *
 * SPDX-License-Identifier: LGPL-2.1-or-later
 * Copyright (C) 2024 - 2025 by the deal.II authors
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
 * This program was contributed by Vladimir Yushutin and Timo Heister, Clemson
 * University, 2023.
 */

#include <deal.II/base/convergence_table.h>
#include <deal.II/base/function.h>
#include <deal.II/base/numbers.h>
#include <deal.II/base/point.h>
#include <deal.II/base/quadrature.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/timer.h>
#include <deal.II/distributed/grid_refinement.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_nothing.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_update_flags.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q1.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>
#include <deal.II/hp/fe_collection.h>
#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_control.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/sparse_direct.h>
#include <deal.II/lac/sparsity_pattern.h>
#include <deal.II/lac/sparsity_tools.h>
#include <deal.II/lac/trilinos_precondition.h>
#include <deal.II/lac/trilinos_solver.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/trilinos_vector.h>
#include <deal.II/lac/vector.h>
#include <deal.II/meshworker/mesh_loop.h>
#include <deal.II/meshworker/scratch_data.h>
#include <deal.II/non_matching/fe_immersed_values.h>
#include <deal.II/non_matching/fe_values.h>
#include <deal.II/non_matching/mesh_classifier.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/vector_tools.h>

namespace Step90
{
  using namespace dealii;
  using VectorType = TrilinosWrappers::MPI::Vector;
  using MatrixType = TrilinosWrappers::SparseMatrix;

  // The parallelization in this tutorial relies on the Trilinos library. We
  // will grant to some cells empty finite element spaces FE_Nothing as done in
  // step-85, but this time active DoFs will be only assigned to cell which are
  // intersected by the surface approximation.
  enum class ActiveFEIndex : types::fe_index
  {
    lagrange = 0,
    nothing  = 1
  };

  // @sect3{Exact surface}
  // The following class defines the surface using the implicit level set
  // representation. The exact surface normal uses the Cartesian gradient of the
  // level set function. The exact Hessian is needed for the construction of the
  // test case only.
  template <int dim>
  class TamarindShape : public Function<dim>
  {
  public:
    TamarindShape()
      : Function<dim>()
    {}
    double value(const Point<dim>  &point,
                 const unsigned int component = 0) const override
    {
      AssertIndexRange(component, this->n_components);
      (void)component;
      Assert(dim == 3, ExcNotImplemented());

      return 0.25 * Utilities::pow(point[0], 2) + Utilities::pow(point[1], 2) +
             4.0 * Utilities::pow(point[2], 2) *
               std::pow(1.0 + 0.5 * std::sin(numbers::PI * point[0]), -2) -
             1.0;
    }

    Tensor<1, dim> gradient(const Point<dim>  &point,
                            const unsigned int component = 0) const override
    {
      AssertIndexRange(component, this->n_components);
      (void)component;
      Assert(dim == 3, ExcNotImplemented());

      Tensor<1, dim> grad;
      grad[0] = 0.5 * point[0] +
                (-2.0) * 4.0 * Utilities::pow(point[2], 2) *
                  std::pow(1.0 + 0.5 * std::sin(numbers::PI * point[0]), -3) *
                  (0.5 * numbers::PI * std::cos(numbers::PI * point[0]));
      grad[1] = 2.0 * point[1];
      grad[2] = (2.0) * 4.0 * point[2] *
                std::pow(1.0 + 0.5 * std::sin(numbers::PI * point[0]), -2);

      return grad;
    }

    SymmetricTensor<2, dim>
    hessian(const Point<dim>  &point,
            const unsigned int component = 0) const override
    {
      AssertIndexRange(component, this->n_components);
      (void)component;
      Assert(dim == 3, ExcNotImplemented());

      SymmetricTensor<2, dim> hessian;

      hessian[0][0] =
        0.5 +
        8.0 * Utilities::pow(point[2], 2) *
          (3.0 * std::pow(1.0 + 0.5 * std::sin(numbers::PI * point[0]), -4) *
             Utilities::pow(0.5 * numbers::PI *
                              std::cos(numbers::PI * point[0]),
                            2) +
           std::pow(1.0 + 0.5 * std::sin(numbers::PI * point[0]), -3) * 0.5 *
             numbers::PI * numbers::PI * std::sin(numbers::PI * point[0]));
      hessian[0][1] = 0.0;
      hessian[0][2] =
        (-8.0) * point[2] *
        std::pow(1.0 + 0.5 * std::sin(numbers::PI * point[0]), -3) *
        numbers::PI * std::cos(numbers::PI * point[0]);

      hessian[1][1] = 2.0;
      hessian[1][2] = 0.0;

      hessian[2][2] =
        8.0 * std::pow(1.0 + 0.5 * std::sin(numbers::PI * point[0]), -2);

      return hessian;
    }
  };

  // @sect3{Exact solution}
  // The following class defines the chosen exact solution and its surface
  // gradient. The exact solution we try to reproduce is $u=xy$ and it may be
  // evaluated away from
  // $\Gamma$ as any other function of Cartesian points. Also note that the
  // gradient() method returns the surface gradient $\nabla_\Gamma u$ of the
  // exact solution.
  template <int dim>
  class AnalyticalSolution : public Function<dim>
  {
  private:
    const TamarindShape<dim> tamarind;

  public:
    AnalyticalSolution()
      : Function<dim>()
    {}
    double value(const Point<dim>  &point,
                 const unsigned int component = 0) const override;

    Tensor<1, dim> gradient(const Point<dim>  &point,
                            const unsigned int component = 0) const override;
  };

  template <int dim>
  double AnalyticalSolution<dim>::value(const Point<dim>  &point,
                                        const unsigned int component) const
  {
    AssertIndexRange(component, this->n_components);
    (void)component;
    return point[0] * point[1];
  }

  template <int dim>
  Tensor<1, dim>
  AnalyticalSolution<dim>::gradient(const Point<dim>  &point,
                                    const unsigned int component) const
  {
    AssertIndexRange(component, this->n_components);
    (void)component;

    const Tensor<1, dim> grad   = tamarind.gradient(point, component);
    const Tensor<1, dim> normal = grad / grad.norm();

    Tensor<1, dim> projector_first_column = -normal[0] * normal;
    projector_first_column[0] += 1.0;

    Tensor<1, dim> projector_second_column = -normal[1] * normal;
    projector_second_column[1] += 1.0;

    Tensor<1, dim> surface_gradient =
      point[1] * projector_first_column + point[0] * projector_second_column;

    return surface_gradient;
  }

  // @sect3{Exact forcing}
  // We choose the right hand side equal to the evaluation of the surface
  // Laplacian for a manufactured solution $u$.
  // This corresponds to the exact forcing $f=-\Delta_\Gamma u+u$:
  template <int dim>
  class RightHandSide : public Function<dim>
  {
    const TamarindShape<dim> tamarind;

  public:
    RightHandSide()
      : Function<dim>()
    {}

    virtual double value(const Point<dim>  &p,
                         const unsigned int component = 0) const override;
  };

  template <int dim>
  double RightHandSide<dim>::value(const Point<dim>  &point,
                                   const unsigned int component) const
  {
    AssertIndexRange(component, this->n_components);
    (void)component;
    Assert(dim == 3, ExcNotImplemented());

    const Tensor<1, dim>          grad    = tamarind.gradient(point, component);
    const Tensor<1, dim>          normal  = grad / grad.norm();
    const SymmetricTensor<2, dim> hessian = tamarind.hessian(point, component);

    double mean_curv = 0.0;
    for (int j = 0; j < 3; j++)
      for (int k = 0; k < 3; k++)
        mean_curv += ((j == k ? 1 : 0) - normal[j] * normal[k]) * hessian[j][k];
    mean_curv /= grad.norm();

    return point[0] * point[1] + 2.0 * normal[0] * normal[1] +
           mean_curv * (point[1] * normal[0] + point[0] * normal[1]);
  }

  // @sect3{Scratch and Copy objects for TraceFEM}
  // Since the assembly procedure will be performed via MeshWorker, we need a
  // Scratch object that handles the Non-Matching FEValues effectively.
  // The input arguments of its constructor are discussed in the solver class
  // below.
  template <int dim>
  struct ScratchData
  {
    ScratchData(const Mapping<dim>                     &mapping,
                const hp::FECollection<dim>            &fe_collection,
                const NonMatching::MeshClassifier<dim> &mesh_classifier,
                const DoFHandler<dim>                  &level_set_dof_handler,
                const VectorType                       &level_set,
                const NonMatching::RegionUpdateFlags nonmatching_update_flags,
                const Quadrature<dim>               &quadrature,
                const Quadrature<1>                 &quadrature_edge,
                const UpdateFlags cell_update_flags = update_values |
                                                      update_gradients |
                                                      update_quadrature_points |
                                                      update_JxW_values)
      : fe_values(
          mapping,
          fe_collection[static_cast<types::fe_index>(ActiveFEIndex::lagrange)],
          quadrature,
          cell_update_flags)
      , region_update_flags(nonmatching_update_flags)
      , quadrature_1D(quadrature_edge)
      , fe_collection(fe_collection)
      , mesh_classifier(mesh_classifier)
      , level_set_dof_handler(level_set_dof_handler)
      , level_set(level_set)
      , level_set_fe_values(mapping,
                            level_set_dof_handler.get_fe(),
                            quadrature,
                            cell_update_flags)
      , non_matching_fe_values(fe_collection,
                               quadrature_edge,
                               nonmatching_update_flags,
                               mesh_classifier,
                               level_set_dof_handler,
                               level_set)
    {}

    ScratchData(const ScratchData<dim> &scratch_data)
      : fe_values(scratch_data.fe_values.get_mapping(),
                  scratch_data.fe_values.get_fe(),
                  scratch_data.fe_values.get_quadrature(),
                  scratch_data.fe_values.get_update_flags())
      , region_update_flags(scratch_data.region_update_flags)
      , quadrature_1D(scratch_data.quadrature_1D)
      , fe_collection(scratch_data.fe_collection)
      , mesh_classifier(scratch_data.mesh_classifier)
      , level_set_dof_handler(scratch_data.level_set_dof_handler)
      , level_set(scratch_data.level_set)
      , level_set_fe_values(scratch_data.level_set_fe_values.get_mapping(),
                            scratch_data.level_set_fe_values.get_fe(),
                            scratch_data.level_set_fe_values.get_quadrature(),
                            scratch_data.level_set_fe_values.get_update_flags())
      , non_matching_fe_values(fe_collection,
                               quadrature_1D,
                               region_update_flags,
                               mesh_classifier,
                               level_set_dof_handler,
                               level_set)
    {}

    // The following FEValues object is used for the standard quadrature on
    // cells involving the FE space of the solution. In TraceFEM, we need this
    // quadrature due to the stabilization term.  In addition, a cell quadrature
    // for the FE space of the level set is defined.
    FEValues<dim>                           fe_values;
    const NonMatching::RegionUpdateFlags    region_update_flags;
    const Quadrature<1>                    &quadrature_1D;
    const hp::FECollection<dim>            &fe_collection;
    const NonMatching::MeshClassifier<dim> &mesh_classifier;
    const DoFHandler<dim>                  &level_set_dof_handler;
    const VectorType                       &level_set;
    FEValues<dim>                           level_set_fe_values;
    NonMatching::FEValues<dim>              non_matching_fe_values;
  };

  // The MeshWorker framework also requires a "copy" data structure that is
  // filled by the worker function working on a cell or face, and whose contents
  // are then later copied into global matrices and vectors. This CopyData
  // object is customized for TraceFEM. In particular, the implementation of the
  // normal-gradient volume stabilization relies on it.
  template <int dim>
  struct CopyData
  {
    FullMatrix<double>                   cell_matrix;
    Vector<double>                       cell_rhs;
    std::vector<types::global_dof_index> local_dof_indices;

    void reinit(const typename DoFHandler<dim>::active_cell_iterator &cell)
    {
      const unsigned int dofs_per_cell = cell->get_fe().n_dofs_per_cell();
      cell_matrix.reinit(dofs_per_cell, dofs_per_cell);
      cell_rhs.reinit(dofs_per_cell);
      local_dof_indices.resize(dofs_per_cell);
      cell->get_dof_indices(local_dof_indices);
    }
  };

  template <int dim>
  struct CopyDataError
  {
    unsigned int cell_index;
    double       cell_L2_error_sqr;
    double       cell_H1_error_sqr;
    double       cell_stab_sqr;

    void reinit(const typename DoFHandler<dim>::active_cell_iterator &cell)
    {
      cell_index        = cell->active_cell_index();
      cell_L2_error_sqr = 0.0;
      cell_H1_error_sqr = 0.0;
      cell_stab_sqr     = 0.0;
    }
  };

  // @sect3{Normal-gradient stabilization form of TraceFEM}
  // The following class corresponds to the stabilization form,
  // its contribution to the global matrix and to the error. More specifically,
  // the method needs_cell_worker() indicates
  // whether the bilinear form of the stabilization, unlike the main bilinear
  // form of Laplace-Beltrami operator, needs the bulk cell quadratures. The
  // cell worker which is useful in an accumulation by MeshWorkers is provided
  // by the assemble_cell_worker() method. The remaining method
  // evaluate_cell_worker() computes the stabilization error for the solution
  // $u_h$, i.e $s_h(u_h,u_h)$. Also note that the method needs_cell_worker()
  // indicates that the assembly and the evaluation of the form does require a
  // bulk cell quadrature. This methodology may be utilized in the MeshWorker.
  // The stabilization scaling is specified by
  // $\mathrm{stabilization\_parameter}\cdot
  // h^\mathrm{stabilization\_exponent}$. For elliptic problems with smooth
  // solutions we can choose any
  // $-1\leq \mathrm{stabilization\_exponent} \leq 1$ and
  // a sufficiently large $\mathrm{stabilization\_parameter}$ that depends of
  // $\Gamma$.
  template <int dim>
  class NormalGradientVolumeStabilization
  {
  public:
    NormalGradientVolumeStabilization()
      : stabilization_parameter(1.0)
      , stabilization_exponent(-1.0)
    {}

    bool needs_cell_worker() const
    {
      return true;
    }

    // We define the stabilization form here assuming that ScratchData and
    // CopyData arguments are initialized properly. The local contribution of
    // the stabilization from this cell to the global matrix is given in
    // assemble_cell_worker() and, later in evaluate_cell_worker(), the
    // local bilinear form of the stabilization is evaluated on the solution.
    // Note the gradients of the discrete level set are computed
    // in the bulk cell quadrature points, which, upon normalization, give the
    // discrete normal vector in a bulk cell.
    void assemble_cell_worker(
      VectorType                                           &level_set,
      const typename DoFHandler<dim>::active_cell_iterator &cell,
      ScratchData<dim>                                     &scratch_data,
      CopyData<dim>                                        &copy_data) const
    {
      const FEValues<dim> &fe_values = scratch_data.fe_values;
      const FEValues<dim> &level_set_fe_values =
        scratch_data.level_set_fe_values;

      const std::vector<double> &JxW_cell = fe_values.get_JxW_values();

      std::vector<Tensor<1, dim>> grad_level_set(
        level_set_fe_values.get_quadrature().size());
      level_set_fe_values.get_function_gradients(level_set, grad_level_set);

      const double factor =
        stabilization_parameter *
        std::pow(cell->minimum_vertex_distance(), stabilization_exponent);
      for (const unsigned int q : fe_values.quadrature_point_indices())
        {
          const Tensor<1, dim> &normal =
            grad_level_set[q] / grad_level_set[q].norm();
          for (const unsigned int i : fe_values.dof_indices())
            for (const unsigned int j : fe_values.dof_indices())
              copy_data.cell_matrix(i, j) +=
                factor * (normal * fe_values.shape_grad(i, q)) *
                (normal * fe_values.shape_grad(j, q)) * JxW_cell[q];
        }
    }

    void evaluate_cell_worker(
      VectorType                                           &solution,
      VectorType                                           &level_set,
      const typename DoFHandler<dim>::active_cell_iterator &cell,
      ScratchData<dim>                                     &scratch_data,
      CopyDataError<dim>                                   &copy_data) const
    {
      double                     cell_stab_sqr = 0.0;
      const FEValues<dim>       &fe_values     = scratch_data.fe_values;
      const std::vector<double> &JxW_cell      = fe_values.get_JxW_values();
      const unsigned int n_q_points = fe_values.get_quadrature_points().size();
      const FEValues<dim> &level_set_fe_values =
        scratch_data.level_set_fe_values;

      std::vector<Tensor<1, dim>> level_set_grad(n_q_points);
      level_set_fe_values.get_function_gradients(level_set, level_set_grad);

      std::vector<Tensor<1, dim>> sol_grad(n_q_points);
      fe_values.get_function_gradients(solution, sol_grad);

      const double factor =
        stabilization_parameter *
        std::pow(cell->minimum_vertex_distance(), stabilization_exponent);

      for (const unsigned int q : fe_values.quadrature_point_indices())
        {
          const Tensor<1, dim> normal =
            level_set_grad[q] / level_set_grad[q].norm();

          const double stabilization_at_point = normal * sol_grad[q];
          cell_stab_sqr +=
            factor * Utilities::pow(stabilization_at_point, 2) * JxW_cell[q];
        }
      copy_data.cell_stab_sqr = cell_stab_sqr;
    }

  private:
    const double stabilization_parameter;
    const double stabilization_exponent;
  };

  // @sect3{Laplace--Beltrami solver}
  // The main class whose method run() performs the computation.
  // One may adjust main parameters of TraceFEM in the constructor.
  // The other methods are discussed below.
  template <int dim>
  class LaplaceBeltramiSolver
  {
  public:
    LaplaceBeltramiSolver();
    void run();

  private:
    void make_grid();

    void localize_surface();

    void setup_discrete_level_set();

    void distribute_dofs();

    void initialize_matrices();

    void assemble_system();

    void solve();

    void mark_intersected();

    void refine_grid();

    void compute_errors();

    void output_level_set(unsigned int);

    void output_solution();

    MPI_Comm mpi_communicator;

    // The surface of interest corresponds to the zero contour of the following
    // exact level set function:
    const TamarindShape<dim> tamarind;

    // The manufactured solution to the Laplace--Beltrami problem and the
    // corresponding right-hand side:
    const AnalyticalSolution<dim> analytical_solution;
    const RightHandSide<dim>      right_hand_side;

    // There is a single triangulation which is shared by the discretizations of
    // the solution and of the level set.
    parallel::distributed::Triangulation<dim, dim> triangulation;
    ConditionalOStream                             pcout;
    TimerOutput                                    computing_timer;

    // We need two separate FE spaces.
    // The first manages the TraceFEM space which is active on intersected
    // elements. The second manages the discrete
    // level set function that describes the geometry of the surface.
    // Also, the degrees of the FE spaces and the corresponding DoFHandler
    // objects are given in the following:
    const unsigned int    fe_degree;
    hp::FECollection<dim> fe_collection;
    DoFHandler<dim>       dof_handler;

    const unsigned int level_set_fe_degree;
    const FE_Q<dim>    level_set_fe;
    DoFHandler<dim>    level_set_dof_handler;

    const MappingQ1<dim> mapping;

    // Since we will adaptively refine the bulk triangulation, two constraints
    // are needed: one for the solution space and another for the level set
    // space.
    AffineConstraints<double> constraints;
    AffineConstraints<double> level_set_constraints;

    // Discrete vectors initialized with dof_handler and level_set_dof_handler.
    VectorType    completely_distributed_solution;
    VectorType    locally_relevant_solution;
    VectorType    locally_relevant_exact;
    VectorType    level_set;
    Vector<float> active_fe_indicator;

    // The following NonMatching::MeshClassifier object is used to
    // separate intersected elements and non-intersected ones.
    // We will then use different finite elements from an hp::FECollection for
    // these two categories:
    NonMatching::MeshClassifier<dim> mesh_classifier;

    // The first bulk quadrature is required for the
    // for TraceFEM stabilization, while the integration over implicit surface
    // is based on the last, one-dimensional rule.
    const QGauss<dim> cell_quadrature;
    const QGauss<1>   quadrature_1D;

    // Any TraceFEM needs a stabilization, and we choose the normal-gradient,
    // volume stabilization.
    const NormalGradientVolumeStabilization<dim> stabilization_scheme;

    // Discrete right-hand side and the final matrix corresponding to
    // dof_handler.
    VectorType      global_rhs;
    MatrixType      global_matrix;
    SparsityPattern sparsity_pattern;
    IndexSet        locally_owned_dofs;
    IndexSet        locally_relevant_dofs;

    // Depending on the type of the quadrature, surface, face or volume, we need
    // to define different update flags.
    NonMatching::RegionUpdateFlags surface_update_flags;

    // The following variables are used to display the results of the
    // convergence test:
    ConvergenceTable convergence_table;
  };

  template <int dim>
  LaplaceBeltramiSolver<dim>::LaplaceBeltramiSolver()
    : mpi_communicator(MPI_COMM_WORLD)
    , tamarind()
    , analytical_solution()
    , right_hand_side()
    , triangulation(mpi_communicator)
    , pcout(std::cout,
            (Utilities::MPI::this_mpi_process(mpi_communicator) == 0))
    , computing_timer(mpi_communicator,
                      pcout,
                      TimerOutput::never,
                      TimerOutput::wall_times)
    , fe_degree(1)
    , fe_collection(FE_Q<dim>(fe_degree), FE_Nothing<dim>())
    , dof_handler(triangulation)
    , level_set_fe_degree(1)
    , level_set_fe(level_set_fe_degree)
    , level_set_dof_handler(triangulation)
    , mapping()
    , mesh_classifier(level_set_dof_handler, level_set)
    , cell_quadrature(fe_degree + 1)
    , quadrature_1D(fe_degree + 1)
    , stabilization_scheme()
  {
    surface_update_flags.surface =
      update_values | update_gradients | update_JxW_values |
      update_quadrature_points | update_normal_vectors;
  }

  // @sect3{Geometric approximation}
  // Let us start with a function that creates the background mesh, using a
  // domain size chosen to avoid situations in which level set function vanishes
  // at mesh vertices. The initial refinement helps the level set to approximate
  // the surface meaningfully.
  //
  // In following next method we construct the discrete level set and determine
  // which cells are intersected. Note that all cells, intersected and
  // non-intersected, have a corresponding active_fe_indicator.
  // Similarly, the exact level set function is approximated on the whole
  // triangulation and postprocessed afterward  resulting on a surface
  // approximation with no gaps.
  template <int dim>
  void LaplaceBeltramiSolver<dim>::make_grid()
  {
    pcout << "Creating background mesh..."
          << "\n"
          << std::flush;
    const double cube_side = 2.008901281;
    GridGenerator::hyper_cube(triangulation, -cube_side, cube_side);
    triangulation.refine_global(3);
  }

  template <int dim>
  void LaplaceBeltramiSolver<dim>::setup_discrete_level_set()
  {
    pcout
      << "Setting up discrete level set function and reclassifying cells... "
      << "\n"
      << std::flush;
    TimerOutput::Scope t(computing_timer, "setup_level_set");

    active_fe_indicator.reinit(triangulation.n_active_cells());
    level_set_dof_handler.distribute_dofs(level_set_fe);
    level_set_constraints.clear();
    const IndexSet level_set_locally_relevant_dofs =
      DoFTools::extract_locally_relevant_dofs(level_set_dof_handler);
    level_set_constraints.reinit(level_set_dof_handler.locally_owned_dofs(),
                                 level_set_locally_relevant_dofs);
    DoFTools::make_hanging_node_constraints(level_set_dof_handler,
                                            level_set_constraints);
    level_set_constraints.close();

    // Here is where the geometric information enters the code. Next, using the
    // discrete level set, we mark the cell which are intersected by its zero
    // contour. Finally, once the triangulation's cells are classified, we
    // determine which cells are active.
    VectorType tmp_sol(level_set_dof_handler.locally_owned_dofs(),
                       mpi_communicator);
    VectorTools::interpolate(level_set_dof_handler, tamarind, tmp_sol);
    level_set_constraints.distribute(tmp_sol);

    level_set.reinit(level_set_locally_relevant_dofs,
                     level_set_dof_handler.locally_owned_dofs(),
                     mpi_communicator);
    level_set = tmp_sol;

    mesh_classifier.reclassify();

    for (const auto &cell : dof_handler.active_cell_iterators() |
                              IteratorFilters::LocallyOwnedCell())
      {
        if (mesh_classifier.location_to_level_set(cell) ==
            NonMatching::LocationToLevelSet::intersected)
          cell->set_active_fe_index(
            static_cast<types::fe_index>(ActiveFEIndex::lagrange));
        else
          cell->set_active_fe_index(
            static_cast<types::fe_index>(ActiveFEIndex::nothing));
      }
  }

  // The method fills in the indicator telling which cells are intersected. It
  // is used in the adaptive refinement near the surface.
  template <int dim>
  void LaplaceBeltramiSolver<dim>::mark_intersected()
  {
    pcout << "Determining cells with active FE index..."
          << "\n"
          << std::flush;
    for (const auto &cell : dof_handler.active_cell_iterators() |
                              IteratorFilters::LocallyOwnedCell())
      {
        if (mesh_classifier.location_to_level_set(cell) ==
            NonMatching::LocationToLevelSet::intersected)
          active_fe_indicator[cell->active_cell_index()] = 1.0;
      }
  }


  // We refine only intersected cells with active_fe_indicator=1. We are calling
  // GridRefinement::refine_and_coarsen_fixed_fraction() instead of the
  // GridRefinement::refine_and_coarsen_fixed_number() function called in most
  // other tutorials because the number of non-intersected cells also grows
  // interfering with the number of active, intersected cells.
  template <int dim>
  void LaplaceBeltramiSolver<dim>::refine_grid()
  {
    TimerOutput::Scope t(computing_timer, "refine");
    pcout << "Refining near surface..."
          << "\n"
          << std::flush;
    parallel::distributed::GridRefinement::refine_and_coarsen_fixed_fraction(
      triangulation, active_fe_indicator, 1.0, 0.0);

    triangulation.execute_coarsening_and_refinement();
  }

  // As the surface is properly approximated by several adaptive steps, we may
  // now distribute the degrees of
  // freedom on   cells which are intersected by the discrete approximation.
  // Next, we initialize matrices for active DoFs and apply the constraints for
  // the solution.
  template <int dim>
  void LaplaceBeltramiSolver<dim>::distribute_dofs()
  {
    pcout << "Distributing degrees of freedom... "
          << "\n"
          << std::flush;
    dof_handler.distribute_dofs(fe_collection);
    locally_owned_dofs = dof_handler.locally_owned_dofs();
    locally_relevant_dofs =
      DoFTools::extract_locally_relevant_dofs(dof_handler);
    completely_distributed_solution.reinit(dof_handler.locally_owned_dofs(),
                                           mpi_communicator);
    locally_relevant_solution.reinit(locally_owned_dofs,
                                     locally_relevant_dofs,
                                     mpi_communicator);
    global_rhs.reinit(locally_owned_dofs, mpi_communicator);

    const unsigned int dof_handler_size = dof_handler.n_dofs();
    const unsigned int level_set_dof_handler_size =
      level_set_dof_handler.n_dofs();

    convergence_table.add_value("LevelSet dofs", level_set_dof_handler_size);
    convergence_table.evaluate_convergence_rates(
      "LevelSet dofs", ConvergenceTable::reduction_rate_log2);

    convergence_table.add_value("Active dofs", dof_handler_size);
    convergence_table.evaluate_convergence_rates(
      "Active dofs", ConvergenceTable::reduction_rate_log2);
  }

  template <int dim>
  void LaplaceBeltramiSolver<dim>::initialize_matrices()
  {
    pcout << "Initializing the matrix... "
          << "\n"
          << std::flush;

    DynamicSparsityPattern dsp(dof_handler.n_dofs(),
                               dof_handler.n_dofs(),
                               locally_relevant_dofs);
    constraints.reinit(locally_owned_dofs, locally_relevant_dofs);

    DoFTools::make_hanging_node_constraints(dof_handler, constraints);
    constraints.close();
    DoFTools::make_sparsity_pattern(dof_handler, dsp, constraints);

    SparsityTools::distribute_sparsity_pattern(dsp,
                                               locally_owned_dofs,
                                               mpi_communicator,
                                               locally_relevant_dofs);
    global_matrix.reinit(locally_owned_dofs,
                         locally_owned_dofs,
                         dsp,
                         mpi_communicator);
  }

  // @sect3{Assembly and surface accumulation}
  // We use a MeshWorker to assemble the linear problem efficiently.
  // This cell worker does not do anything for non-intersected cells.
  template <int dim>
  void LaplaceBeltramiSolver<dim>::assemble_system()
  {
    pcout << "Assembling... "
          << "\n"
          << std::flush;
    TimerOutput::Scope t(computing_timer, "assembly");

    const auto cell_worker =
      [&](const typename DoFHandler<dim>::active_cell_iterator &cell,
          ScratchData<dim>                                     &scratch_data,
          CopyData<dim>                                        &copy_data) {
        if (mesh_classifier.location_to_level_set(cell) ==
              NonMatching::LocationToLevelSet::intersected &&
            cell->is_locally_owned())
          {
            // Once we know that the cell is intersected, we construct the
            // unfitted quadratures for the solutions FE space on the cell.
            scratch_data.non_matching_fe_values.reinit(cell);
            copy_data.reinit(cell);
            copy_data.cell_matrix = 0;
            copy_data.cell_rhs    = 0;
            const std::optional<NonMatching::FEImmersedSurfaceValues<dim>>
              &surface_fe_values =
                scratch_data.non_matching_fe_values.get_surface_fe_values();

            if (surface_fe_values)
              {
                const std::vector<double> &JxW_surface =
                  surface_fe_values->get_JxW_values();

                // The accumulation of the surface integrals, including the
                // forcing, is performed here.
                for (unsigned int q :
                     surface_fe_values->quadrature_point_indices())
                  {
                    const Tensor<1, dim> &normal =
                      surface_fe_values->normal_vector(q);

                    for (const unsigned int i :
                         surface_fe_values->dof_indices())
                      {
                        copy_data.cell_rhs(i) +=
                          surface_fe_values->shape_value(i, q) *
                          right_hand_side.value(
                            surface_fe_values->quadrature_point(q)) *
                          JxW_surface[q];

                        for (const unsigned int j :
                             surface_fe_values->dof_indices())
                          {
                            copy_data.cell_matrix(i, j) +=
                              (surface_fe_values->shape_value(i, q) *
                               surface_fe_values->shape_value(j, q)) *
                              JxW_surface[q];
                            copy_data.cell_matrix(i, j) +=
                              (surface_fe_values->shape_grad(i, q) -
                               (normal * surface_fe_values->shape_grad(i, q)) *
                                 normal) *
                              (surface_fe_values->shape_grad(j, q) -
                               (normal * surface_fe_values->shape_grad(j, q)) *
                                 normal) *
                              JxW_surface[q];
                          }
                      }
                  }

                // The normal-gradient volume stabilization form needs a bulk
                // cell integration while other types of stabilization may need
                // face quadratures, for example. So we check it first. The cell
                // was provided by the solution's DoFHandler,
                //  so we recast it as a level set's DoFHandler cell.
                //  However, it is the same geometric entity of the common
                //  triangulation.
                if (stabilization_scheme.needs_cell_worker())
                  {
                    typename DoFHandler<dim>::active_cell_iterator
                      level_set_cell =
                        cell->as_dof_handler_iterator(level_set_dof_handler);
                    scratch_data.fe_values.reinit(cell);
                    scratch_data.level_set_fe_values.reinit(level_set_cell);
                    stabilization_scheme.assemble_cell_worker(level_set,
                                                              cell,
                                                              scratch_data,
                                                              copy_data);
                  }
              }
          }
      };

    // Next, the copier worker distributes the local contributions from
    // the CopyData taking into account the constraints. Finally, the
    // MeshWorker goes over all cells provided by the solutions'
    // DoFHandler. Note that this includes non-intersected cells as
    // well, but the cell worker does nothing on them.
    const auto copier = [&](const CopyData<dim> &c) {
      constraints.distribute_local_to_global(c.cell_matrix,
                                             c.cell_rhs,
                                             c.local_dof_indices,
                                             global_matrix,
                                             global_rhs);
    };

    ScratchData<dim> scratch_data(mapping,
                                  fe_collection,
                                  mesh_classifier,
                                  level_set_dof_handler,
                                  level_set,
                                  surface_update_flags,
                                  cell_quadrature,
                                  quadrature_1D);

    CopyData<dim> copy_data;

    MeshWorker::mesh_loop(dof_handler.begin_active(),
                          dof_handler.end(),
                          cell_worker,
                          copier,
                          scratch_data,
                          copy_data,
                          MeshWorker::assemble_own_cells);

    global_matrix.compress(VectorOperation::add);
    global_rhs.compress(VectorOperation::add);
  }

  // In the following, we solve the resulting linear system of equations. We
  // either use a direct solver or AMG.
  template <int dim>
  void LaplaceBeltramiSolver<dim>::solve()
  {
    TimerOutput::Scope t(computing_timer, "solve");
    bool               apply_direct_solver = false;
    const double       relative_error      = 1e-9 * global_rhs.l2_norm();
    unsigned int       n_iterations        = 0;
    if (apply_direct_solver)
      {
        pcout << "Solving directly... " << '\n' << std::flush;
        SolverControl solver_control(100, relative_error);
        TrilinosWrappers::SolverDirect::AdditionalData data;
        TrilinosWrappers::SolverDirect trilinos(solver_control, data);
        trilinos.solve(global_matrix,
                       completely_distributed_solution,
                       global_rhs);
      }
    else
      {
        Timer timer;
        pcout << "Solving with AMG... "
              << "\n"
              << std::flush;
        const unsigned int max_iterations = 500;
        SolverControl      solver_control(max_iterations, relative_error);
        const std::vector<std::vector<bool>> constant_modes =
          DoFTools::extract_constant_modes(dof_handler);
        TrilinosWrappers::PreconditionAMG preconditioner_stiffness;
        TrilinosWrappers::PreconditionAMG::AdditionalData Amg_data;
        Amg_data.constant_modes        = constant_modes;
        Amg_data.elliptic              = true;
        Amg_data.higher_order_elements = false;
        Amg_data.smoother_sweeps       = 2;
        Amg_data.aggregation_threshold = 0.02;
        Amg_data.output_details        = true;
        preconditioner_stiffness.initialize(global_matrix);

        SolverCG<VectorType> cg(solver_control);
        cg.solve(global_matrix,
                 completely_distributed_solution,
                 global_rhs,
                 preconditioner_stiffness);
        n_iterations = solver_control.last_step();
      }
    constraints.distribute(completely_distributed_solution);
    locally_relevant_solution = completely_distributed_solution;

    convergence_table.add_value("Iterations", n_iterations);
  }

  // Similarly to what we do in the assembly() function,
  // a MeshWorker is used to accumulate errors
  // including the stabilization term. At the end,  we collect the results,
  // and print them out.
  template <int dim>
  void LaplaceBeltramiSolver<dim>::compute_errors()
  {
    pcout << "Evaluating errors on the surface..."
          << "\n"
          << std::flush;
    TimerOutput::Scope t(computing_timer, "eval_errors");
    double             error_L2_sqr   = 0.0;
    double             error_H1_sqr   = 0.0;
    double             error_stab_sqr = 0.0;
    const auto         cell_worker    = [&](const auto &cell,
                                 auto       &scratch_data,
                                 auto       &copy_data) {
      if (mesh_classifier.location_to_level_set(cell) ==
            NonMatching::LocationToLevelSet::intersected &&
          cell->is_locally_owned())
        {
          double cell_L2_error_sqr = 0.0;
          double cell_H1_error_sqr = 0.0;

          copy_data.reinit(cell);
          scratch_data.non_matching_fe_values.reinit(cell);

          const std::optional<NonMatching::FEImmersedSurfaceValues<dim>>
            &surface_fe_values =
              scratch_data.non_matching_fe_values.get_surface_fe_values();

          if (surface_fe_values)
            {
              const std::vector<double> &JxW_surface =
                surface_fe_values->get_JxW_values();
              const unsigned int n_q_points =
                surface_fe_values->n_quadrature_points;

              std::vector<double> sol(n_q_points);
              surface_fe_values->get_function_values(locally_relevant_solution,
                                                     sol);

              std::vector<Tensor<1, dim>> sol_grad(n_q_points);
              surface_fe_values->get_function_gradients(
                locally_relevant_solution, sol_grad);

              for (const unsigned int q :
                   surface_fe_values->quadrature_point_indices())
                {
                  const Point<dim> &point =
                    surface_fe_values->quadrature_point(q);
                  const Tensor<1, dim> &normal =
                    surface_fe_values->normal_vector(q);
                  const double error_at_point =
                    sol.at(q) - analytical_solution.value(point);
                  const Tensor<1, dim> grad_error_at_point =
                    (sol_grad.at(q) - (normal * sol_grad.at(q)) * normal -
                     analytical_solution.gradient(point));

                  cell_L2_error_sqr +=
                    Utilities::pow(error_at_point, 2) * JxW_surface[q];
                  cell_H1_error_sqr +=
                    grad_error_at_point * grad_error_at_point * JxW_surface[q];
                }
              copy_data.cell_L2_error_sqr = cell_L2_error_sqr;
              copy_data.cell_H1_error_sqr = cell_H1_error_sqr;

              if (stabilization_scheme.needs_cell_worker())
                {
                  typename DoFHandler<dim>::active_cell_iterator
                    level_set_cell =
                      cell->as_dof_handler_iterator(level_set_dof_handler);
                  scratch_data.fe_values.reinit(cell);
                  scratch_data.level_set_fe_values.reinit(level_set_cell);
                  stabilization_scheme.evaluate_cell_worker(
                    locally_relevant_solution,
                    level_set,
                    cell,
                    scratch_data,
                    copy_data);
                }
            }
        }
    };

    const auto copier = [&](const auto &copy_data) {
      if (copy_data.cell_index < active_fe_indicator.size())
        {
          error_L2_sqr += copy_data.cell_L2_error_sqr;
          error_H1_sqr += copy_data.cell_H1_error_sqr;
          error_stab_sqr += copy_data.cell_stab_sqr;
        }
    };

    ScratchData<dim> scratch_data(mapping,
                                  fe_collection,
                                  mesh_classifier,
                                  level_set_dof_handler,
                                  level_set,
                                  surface_update_flags,
                                  cell_quadrature,
                                  quadrature_1D);

    CopyDataError<dim> copy_data;

    MeshWorker::mesh_loop(dof_handler.begin_active(),
                          dof_handler.end(),
                          cell_worker,
                          copier,
                          scratch_data,
                          copy_data,
                          MeshWorker::assemble_own_cells);

    const double error_L2 =
      std::sqrt(Utilities::MPI::sum(error_L2_sqr, mpi_communicator));
    const double error_semiH1 =
      std::sqrt(Utilities::MPI::sum(error_H1_sqr, mpi_communicator));
    const double error_stab =
      std::sqrt(Utilities::MPI::sum(error_stab_sqr, mpi_communicator));

    convergence_table.add_value("L2 Error", error_L2);
    convergence_table.evaluate_convergence_rates(
      "L2 Error", ConvergenceTable::reduction_rate_log2);
    convergence_table.set_scientific("L2 Error", true);

    convergence_table.add_value("H1 error", error_semiH1);
    convergence_table.evaluate_convergence_rates(
      "H1 error", ConvergenceTable::reduction_rate_log2);
    convergence_table.set_scientific("H1 error", true);

    convergence_table.add_value("Stab norm", error_stab);
    convergence_table.evaluate_convergence_rates(
      "Stab norm", ConvergenceTable::reduction_rate_log2);
    convergence_table.set_scientific("Stab norm", true);
  }

  // The following two methods perform VTK output of the preliminary mesh
  // refinements for the geometry approximation and of the TraceFEM solution.
  // The important difference between the two is that the non-intersected cells
  // are excluded from the output saving considerable amount of time and
  // storage.
  template <int dim>
  void LaplaceBeltramiSolver<dim>::output_level_set(const unsigned int cycle)
  {
    pcout << "Writing vtu file for surface... " << '\n' << std::flush;
    TimerOutput::Scope t(computing_timer, "output_level_set");
    DataOut<dim>       data_out;
    data_out.add_data_vector(level_set_dof_handler, level_set, "level_set");
    data_out.add_data_vector(active_fe_indicator, "ref_indicator");
    data_out.build_patches();

    data_out.write_vtu_in_parallel("surface_" + std::to_string(cycle) + ".vtu",
                                   mpi_communicator);
  }

  template <int dim>
  void LaplaceBeltramiSolver<dim>::output_solution()
  {
    pcout << "Writing vtu file... " << std::flush;
    TimerOutput::Scope t(computing_timer, "output_solution");
    Vector<double>     exact(dof_handler.locally_owned_dofs().size());

    VectorTools::interpolate(dof_handler, analytical_solution, exact);
    DataOut<dim> data_out;
    data_out.add_data_vector(dof_handler,
                             locally_relevant_solution,
                             "solution");
    data_out.add_data_vector(dof_handler, exact, "exact");
    data_out.add_data_vector(level_set_dof_handler, level_set, "level_set");

    data_out.set_cell_selection(
      [this](const typename Triangulation<dim>::cell_iterator &cell) {
        return cell->is_active() && cell->is_locally_owned() &&
               mesh_classifier.location_to_level_set(cell) ==
                 NonMatching::LocationToLevelSet::intersected;
      });
    data_out.build_patches();

    data_out.write_vtu_in_parallel("solution.vtu", mpi_communicator);
  }

  // The method localize_surface() generates iteratively a surface approximation
  // as described above. Once the surface approximation is constructed, the main
  // logic of the solver is executed as presented in the method run().
  template <int dim>
  void LaplaceBeltramiSolver<dim>::localize_surface()
  {
    unsigned int preliminary_levels = 3;
    for (unsigned int localization_cycle = 0;
         localization_cycle < preliminary_levels;
         ++localization_cycle)
      {
        pcout << std::endl
              << "Preliminary refinement #" << localization_cycle << std::endl;
        setup_discrete_level_set();
        mark_intersected();
        output_level_set(localization_cycle);
        refine_grid();
      }
    computing_timer.reset();
  }

  template <int dim>
  void LaplaceBeltramiSolver<dim>::run()
  {
    make_grid();
    localize_surface();
    const unsigned int convergence_levels = 3;
    for (unsigned int cycle = 0; cycle < convergence_levels; ++cycle)
      {
        pcout << std::endl << "Convergence refinement #" << cycle << std::endl;
        setup_discrete_level_set();
        distribute_dofs();
        initialize_matrices();
        assemble_system();
        solve();
        compute_errors();
        if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
          convergence_table.write_text(pcout.get_stream());

        computing_timer.print_summary();
        computing_timer.reset();
        if (cycle < convergence_levels - 1)
          {
            mark_intersected();
            refine_grid();
          }
        else
          output_solution();

        computing_timer.print_summary();
        computing_timer.reset();
      }
  }
} // namespace Step90

int main(int argc, char *argv[])
{
  try
    {
      using namespace dealii;
      using namespace Step90;
      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
      LaplaceBeltramiSolver<3>         LB_solver;
      LB_solver.run();
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
