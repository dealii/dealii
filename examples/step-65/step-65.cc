/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2019 by the deal.II authors
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

 * This tutorial program was contributed by Martin Kronbichler
 */

// @sect3{Include files}

// The include files for this tutorial are essentially the same as in
// step-6. Importantly, the TransfiniteInterpolationManifold class we
// will be using is provided by `deal.II/grid/manifold_lib.h`.

#include <deal.II/base/timer.h>

#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/manifold_lib.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q_generic.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/vector_tools.h>

#include <fstream>

// The only new include file is the one for the MappingQCache class.
#include <deal.II/fe/mapping_q_cache.h>


namespace Step65
{
  using namespace dealii;


  // @sect3{Analytical solution and coefficient}

  // In this tutorial program, we want to solve the Poisson equation
  // with a coefficient that jumps along a sphere of radius 0.5, and
  // using a constant right hand side of value $f(\mathbf{x}) = -3$. (This
  // setup is similar to step-5 and step-6, but the concrete values
  // for the coefficient and the right hand side are different.)
  // Due to the jump in the
  // coefficient, the analytical solution must have a kink where the
  // coefficient switches from one value to the other. To keep things simple,
  // we select an analytical solution that is quadratic in all components,
  // i.e., $u(x,y,z) = x^2 + y^2 + z^2$ in the ball of radius 0.5 and
  // $u(x,y,z) = 0.1(x^2 + y^2 + z^2) + 0.25-0.025$ in the outer part of the
  // domain. This analytical solution is compatible with the right hand side
  // in case the coefficient is 0.5 in the inner ball and 5 outside. It is
  // also continuous along the circle of radius 0.5.
  template <int dim>
  class ExactSolution : public Function<dim>
  {
  public:
    virtual double value(const Point<dim> &p,
                         const unsigned int /*component*/ = 0) const override
    {
      if (p.norm_square() < 0.25)
        return p.norm_square();
      else
        return 0.1 * p.norm_square() + (0.25 - 0.025);
    }

    virtual Tensor<1, dim>
    gradient(const Point<dim> &p,
             const unsigned int /*component*/ = 0) const override
    {
      if (p.norm_square() < 0.25)
        return 2. * p;
      else
        return 0.2 * p;
    }
  };


  template <int dim>
  double coefficient(const Point<dim> &p)
  {
    if (p.norm_square() < 0.25)
      return 0.5;
    else
      return 5.0;
  }



  // @sect3{The PoissonProblem class}
  //
  // The implementation of the Poisson problem is very similar to what
  // we used in the step-5 tutorial program. The two main differences
  // are that we pass a mapping object to the various steps in the
  // program in order to switch between two mapping representations as
  // explained in the introduction, and the `timer` object (of type
  // TimerOutput) that will be used for measuring the run times in the
  // various cases. (The concept of mapping objects was first
  // introduced in step-10 and step-11, in case you want to look up
  // the use of these classes.)
  template <int dim>
  class PoissonProblem
  {
  public:
    PoissonProblem();
    void run();

  private:
    void create_grid();
    void setup_system(const Mapping<dim> &mapping);
    void assemble_system(const Mapping<dim> &mapping);
    void solve();
    void postprocess(const Mapping<dim> &mapping);

    Triangulation<dim> triangulation;
    FE_Q<dim>          fe;
    DoFHandler<dim>    dof_handler;

    AffineConstraints<double> constraints;
    SparsityPattern           sparsity_pattern;
    SparseMatrix<double>      system_matrix;
    Vector<double>            solution;
    Vector<double>            system_rhs;

    TimerOutput timer;
  };



  // In the constructor, we set up the timer object to record wall times but
  // be quiet during the normal execution. We will query it for timing details
  // in the `PoissonProblem::run()` function. Furthermore, we select a
  // relatively high polynomial degree of three for the finite element in use.
  template <int dim>
  PoissonProblem<dim>::PoissonProblem()
    : fe(3)
    , dof_handler(triangulation)
    , timer(std::cout, TimerOutput::never, TimerOutput::wall_times)
  {}



  // @sect3{Grid creation and initialization of the manifolds}
  //
  // The next function presents the typical usage of
  // TransfiniteInterpolationManifold. The first step is to create the desired
  // grid, which can be done by composition of two grids from
  // GridGenerator. The inner ball mesh is simple enough: We run
  // GridGenerator::hyper_cube() centered at the origin with radius 0.5 (third
  // function argument). The second mesh is more interesting and constructed
  // as follows: We want to have a mesh that is spherical in the interior but
  // flat on the outer surface. Furthermore, the mesh topology of the inner
  // ball should be compatible with the outer grid in the sense that their
  // vertices coincide so as to allow the two grid to be merged. The grid coming
  // out of GridGenerator::hyper_shell fulfills the requirements on the inner
  // side in case it is created with $2d$ coarse cells (6 coarse cells in 3D
  // which we are going to use) &ndash; this is the same number of cells as
  // there are boundary faces for the ball. For the outer surface, we use the
  // fact that the 6 faces on the surface of the shell without a manifold
  // attached would degenerate to the surface of a cube. What we are still
  // missing is the radius of the outer shell boundary. Since we desire a cube
  // of extent
  // $[-1, 1]$ and the 6-cell shell puts its 8 outer vertices at the 8
  // opposing diagonals, we must translate the points $(\pm 1, \pm 1, \pm 1)$
  // into a radius: Clearly, the radius must be $\sqrt{d}$ in $d$ dimensions,
  // i.e., $\sqrt{3}$ for the three-dimensional case we want to consider.
  //
  // Thus, we have a plan: After creating the inner triangulation for
  // the ball and the one for the outer shell, we merge those two
  // grids but remove all manifolds that the functions in
  // GridGenerator may have set from the resulting triangulation, to
  // ensure that we have full control over manifolds. In particular,
  // we want additional points added on the boundary during refinement
  // to follow a flat manifold description. To start the process of
  // adding more appropriate manifold ids, we assign the manifold id 0
  // to all mesh entities (cells, faces, lines), which will later be
  // associated with the TransfiniteInterpolationManifold. Then, we
  // must identify the faces and lines that are along the sphere of
  // radius 0.5 and mark them with a different manifold id, so as to then
  // assign a SphericalManifold to those. We will choose the manifold
  // id of 1. Since we have thrown away all manifolds that pre-existed
  // after calling GridGenerator::hyper_ball(), we manually go through
  // the cells of the mesh and all their faces. We have found a face
  // on the sphere if all four vertices have a radius of 0.5, or, as
  // we write in the program, have $r^2-0.25 \approx 0$. Note that we call
  // `cell->face(f)->set_all_manifold_ids(1)` to set the manifold id
  // both on the faces and the surrounding lines. Furthermore, we want
  // to distinguish the cells inside the ball and outside the ball by
  // a material id for visualization, corresponding to the picture in the
  // introduction.
  template <int dim>
  void PoissonProblem<dim>::create_grid()
  {
    Triangulation<dim> tria_inner;
    GridGenerator::hyper_ball(tria_inner, Point<dim>(), 0.5);

    Triangulation<dim> tria_outer;
    GridGenerator::hyper_shell(
      tria_outer, Point<dim>(), 0.5, std::sqrt(dim), 2 * dim);

    GridGenerator::merge_triangulations(tria_inner, tria_outer, triangulation);

    triangulation.reset_all_manifolds();
    triangulation.set_all_manifold_ids(0);

    for (const auto &cell : triangulation.cell_iterators())
      {
        for (const auto &face : cell->face_iterators())
          {
            bool face_at_sphere_boundary = true;
            for (unsigned int v = 0;
                 v < GeometryInfo<dim - 1>::vertices_per_cell;
                 ++v)
              {
                if (std::abs(face->vertex(v).norm_square() - 0.25) > 1e-12)
                  face_at_sphere_boundary = false;
              }
            if (face_at_sphere_boundary)
              face->set_all_manifold_ids(1);
          }
        if (cell->center().norm_square() < 0.25)
          cell->set_material_id(1);
        else
          cell->set_material_id(0);
      }

    // With all cells, faces and lines marked appropriately, we can
    // attach the Manifold objects to those numbers. The entities with
    // manifold id 1 will get a spherical manifold, whereas the other
    // entities, which have the manifold id 0, will be assigned the
    // TransfiniteInterpolationManifold. As mentioned in the
    // introduction, we must explicitly initialize the manifold with
    // the current mesh using a call to
    // TransfiniteInterpolationManifold::initialize() in order to pick
    // up the coarse mesh cells and the manifolds attached to the
    // boundaries of those cells. We also note that the manifold
    // objects we create locally in this function are allowed to go
    // out of scope (as they do at the end of the function scope),
    // because the Triangulation object internally copies them.
    //
    // With all manifolds attached, we will finally go about and refine the
    // mesh a few times to create a sufficiently large test case.
    triangulation.set_manifold(1, SphericalManifold<dim>());

    TransfiniteInterpolationManifold<dim> transfinite_manifold;
    transfinite_manifold.initialize(triangulation);
    triangulation.set_manifold(0, transfinite_manifold);

    triangulation.refine_global(9 - 2 * dim);
  }



  // @sect3{Setup of data structures}
  //
  // The following function is well-known from other tutorials in that
  // it enumerates the degrees of freedom, creates a constraint object
  // and sets up a sparse matrix for the linear system. The only thing
  // worth mentioning is the fact that the function receives a
  // reference to a mapping object that we then pass to the
  // VectorTools::interpolate_boundary_values() function to ensure
  // that our boundary values are evaluated on the high-order mesh
  // used for assembly. In the present example, it does not really
  // matter because the outer surfaces are flat, but for curved outer
  // cells this leads to more accurate approximation of the boundary
  // values.
  template <int dim>
  void PoissonProblem<dim>::setup_system(const Mapping<dim> &mapping)
  {
    dof_handler.distribute_dofs(fe);
    std::cout << "   Number of active cells:       "
              << triangulation.n_global_active_cells() << std::endl;
    std::cout << "   Number of degrees of freedom: " << dof_handler.n_dofs()
              << std::endl;

    {
      TimerOutput::Scope scope(timer, "Compute constraints");

      constraints.clear();

      DoFTools::make_hanging_node_constraints(dof_handler, constraints);
      VectorTools::interpolate_boundary_values(
        mapping, dof_handler, 0, ExactSolution<dim>(), constraints);

      constraints.close();
    }

    DynamicSparsityPattern dsp(dof_handler.n_dofs());
    DoFTools::make_sparsity_pattern(dof_handler, dsp, constraints, false);

    sparsity_pattern.copy_from(dsp);
    system_matrix.reinit(sparsity_pattern);

    solution.reinit(dof_handler.n_dofs());
    system_rhs.reinit(dof_handler.n_dofs());
  }


  // @sect3{Assembly of the system matrix and right hand side}
  //
  // The function that assembles the linear system is also well known
  // from the previous tutorial programs. One thing to note is that we
  // set the number of quadrature points to the polynomial degree plus
  // two, not the degree plus one as in most other tutorials. This is
  // because we expect some extra accuracy as the mapping also
  // involves a degree one more than the polynomials for the solution.
  //
  // The only somewhat unusual code in the assembly is the way we compute the
  // cell matrix. Rather than using three nested loop over the quadrature
  // point index, the row, and the column of the matrix, we first collect the
  // derivatives of the shape function, multiplied by the square root of the
  // product of the coefficient and the integration factor `JxW` in a separate
  // matrix `partial_matrix`. To compute the cell matrix, we then execute
  // `cell_matrix = partial_matrix * transpose(partial_matrix)` in the line
  // `partial_matrix.mTmult(cell_matrix, partial_matrix);`. To understand why
  // this works, we realize that the matrix-matrix multiplication performs a
  // summation over the columns of `partial_matrix`. If we denote the
  // coefficient by $a(\mathbf{x}_q)$, the entries in the temporary matrix are
  // $\sqrt{\text{det}(J) w_q a(x)} \frac{\partial \varphi_i(\boldsymbol
  // \xi_q)}{\partial x_k}$. If we take the product of the <i>i</i>th row with
  // the <i>j</i>th column of that matrix, we compute a nested sum involving
  // $\sum_q \sum_{k=1}^d \sqrt{\text{det}(J) w_q a(x)} \frac{\partial
  // \varphi_i(\boldsymbol \xi_q)}{\partial x_k} \sqrt{\text{det}(J) w_q a(x)}
  // \frac{\partial \varphi_j(\boldsymbol \xi_q)}{\partial x_k} = \sum_q
  // \sum_{k=1}^d\text{det}(J) w_q a(x)\frac{\partial \varphi_i(\boldsymbol
  // \xi_q)}{\partial x_k} \frac{\partial \varphi_j(\boldsymbol
  // \xi_q)}{\partial x_k}$, which is exactly the terms needed for the
  // bilinear form of the Laplace equation.
  //
  // The reason for choosing this somewhat unusual scheme is due to the heavy
  // work involved in computing the cell matrix for a relatively high
  // polynomial degree in 3D. As we want to highlight the cost of the mapping
  // in this tutorial program, we better do the assembly in an optimized way
  // in order to not chase bottlenecks that have been solved by the community
  // already. Matrix-matrix multiplication is one of the best optimized
  // kernels in the HPC context, and the FullMatrix::mTmult() function will
  // call into those optimized BLAS functions. If the user has provided a good
  // BLAS library when configuring deal.II (like OpenBLAS or Intel's MKL), the
  // computation of the cell matrix will execute close to the processor's peak
  // arithmetic performance. As a side note, we mention that despite an
  // optimized matrix-matrix multiplication, the current strategy is
  // sub-optimal in terms of complexity as the work to be done is proportional
  // to $(p+1)^9$ operations for degree $p$ (this also applies to the usual
  // evaluation with FEValues). One could compute the cell matrix with
  // $\mathcal O((p+1)^7)$ operations by utilizing the tensor product
  // structure of the shape functions, as is done by the matrix-free framework
  // in deal.II. We refer to step-37 and the documentation of the
  // tensor-product-aware evaluators FEEvaluation for details on how an even
  // more efficient cell matrix computation could be realized.
  template <int dim>
  void PoissonProblem<dim>::assemble_system(const Mapping<dim> &mapping)
  {
    TimerOutput::Scope scope(timer, "Assemble linear system");

    const QGauss<dim> quadrature_formula(fe.degree + 2);
    FEValues<dim>     fe_values(mapping,
                            fe,
                            quadrature_formula,
                            update_values | update_gradients |
                              update_quadrature_points | update_JxW_values);

    const unsigned int dofs_per_cell = fe.dofs_per_cell;
    const unsigned int n_q_points    = quadrature_formula.size();

    FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
    Vector<double>     cell_rhs(dofs_per_cell);
    FullMatrix<double> partial_matrix(dofs_per_cell, dim * n_q_points);
    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

    for (const auto &cell : dof_handler.active_cell_iterators())
      {
        cell_rhs = 0.;
        fe_values.reinit(cell);

        for (unsigned int q_index = 0; q_index < n_q_points; ++q_index)
          {
            const double current_coefficient =
              coefficient(fe_values.quadrature_point(q_index));
            for (unsigned int i = 0; i < dofs_per_cell; ++i)
              {
                for (unsigned int d = 0; d < dim; ++d)
                  partial_matrix(i, q_index * dim + d) =
                    std::sqrt(fe_values.JxW(q_index) * current_coefficient) *
                    fe_values.shape_grad(i, q_index)[d];
                cell_rhs(i) +=
                  (fe_values.shape_value(i, q_index) * // phi_i(x_q)
                   (-dim) *                            // f(x_q)
                   fe_values.JxW(q_index));            // dx
              }
          }

        partial_matrix.mTmult(cell_matrix, partial_matrix);

        cell->get_dof_indices(local_dof_indices);
        constraints.distribute_local_to_global(
          cell_matrix, cell_rhs, local_dof_indices, system_matrix, system_rhs);
      }
  }



  // @sect3{Solution of the linear system}
  //
  // For solving the linear system, we pick a simple Jacobi-preconditioned
  // conjugate gradient solver, similar to the settings in the early tutorials.
  template <int dim>
  void PoissonProblem<dim>::solve()
  {
    TimerOutput::Scope scope(timer, "Solve linear system");

    SolverControl            solver_control(1000, 1e-12);
    SolverCG<Vector<double>> solver(solver_control);

    PreconditionJacobi<SparseMatrix<double>> preconditioner;
    preconditioner.initialize(system_matrix);

    solver.solve(system_matrix, solution, system_rhs, preconditioner);
    constraints.distribute(solution);

    std::cout << "   Number of solver iterations:  "
              << solver_control.last_step() << std::endl;
  }



  // @sect3{Output of the solution and computation of errors}
  //
  // In the next function we do various post-processing steps with the
  // solution, all of which involve the mapping in one way or the other.
  //
  // The first operation we do is to write the solution as well as the
  // material ids to a VTU file. This is similar to what was done in many
  // other tutorial programs. The new ingredient presented in this tutorial
  // program is that we want to ensure that the data written to the file
  // used for visualization is actually a faithful representation of what
  // is used internally by deal.II. That is because most of the visualization
  // data formats only represent cells by their vertex coordinates, but
  // have no way of representing the curved boundaries that are used
  // in deal.II when using higher order mappings -- in other words, what
  // you see in the visualization tool is not actually what you are computing
  // on. (The same, incidentally, is true when using higher order shape
  // functions: Most visualization tools only render bilinear/trilinear
  // representations. This is discussed in detail in DataOut::build_patches().)
  //
  // So we need to ensure that a high-order representation is written
  // to the file. We need to consider two particular topics. Firstly, we tell
  // the DataOut object via the DataOutBase::VtkFlags that we intend to
  // interpret the subdivisions of the elements as a high-order Lagrange
  // polynomial rather than a collection of bilinear patches.
  // Recent visualization programs, like ParaView version 5.5
  // or newer, can then render a high-order solution (see a <a
  // href="https://github.com/dealii/dealii/wiki/Notes-on-visualizing-high-order-output">wiki
  // page</a> for more details).
  // Secondly, we need to make sure that the mapping is passed to the
  // DataOut::build_patches() method. Finally, the DataOut class only prints
  // curved faces for <i>boundary</i> cells by default, so we need to ensure
  // that also inner cells are printed in a curved representation via the
  // mapping.
  template <int dim>
  void PoissonProblem<dim>::postprocess(const Mapping<dim> &mapping)
  {
    {
      TimerOutput::Scope scope(timer, "Write output");

      DataOut<dim> data_out;

      DataOutBase::VtkFlags flags;
      flags.write_higher_order_cells = true;
      data_out.set_flags(flags);

      data_out.attach_dof_handler(dof_handler);
      data_out.add_data_vector(solution, "solution");

      Vector<double> material_ids(triangulation.n_active_cells());
      for (const auto &cell : triangulation.active_cell_iterators())
        material_ids[cell->active_cell_index()] = cell->material_id();
      data_out.add_data_vector(material_ids, "material_ids");

      data_out.build_patches(mapping,
                             fe.degree,
                             DataOut<dim>::curved_inner_cells);

      std::ofstream file(
        ("solution-" +
         std::to_string(triangulation.n_global_levels() - 10 + 2 * dim) +
         ".vtu")
          .c_str());

      data_out.write_vtu(file);
    }

    // The next operation in the postprocessing function is to compute the $L_2$
    // and $H^1$ errors against the analytical solution. As the analytical
    // solution is a quadratic polynomial, we expect a very accurate result at
    // this point. If we were solving on a simple mesh with planar faces and a
    // coefficient whose jumps are aligned with the faces between cells, then
    // we would expect the numerical result to coincide with the
    // analytical solution up to roundoff accuracy. However, since we are using
    // deformed cells following a sphere, which are only tracked by
    // polynomials of degree 4 (one more than the degree for the finite
    // elements), we will see that there is an error around $10^{-7}$. We could
    // get more accuracy by increasing the polynomial degree or refining the
    // mesh.
    {
      TimerOutput::Scope scope(timer, "Compute error norms");

      Vector<double> norm_per_cell_p(triangulation.n_active_cells());

      VectorTools::integrate_difference(mapping,
                                        dof_handler,
                                        solution,
                                        ExactSolution<dim>(),
                                        norm_per_cell_p,
                                        QGauss<dim>(fe.degree + 2),
                                        VectorTools::L2_norm);
      std::cout << "   L2 error vs exact solution:   "
                << norm_per_cell_p.l2_norm() << std::endl;

      VectorTools::integrate_difference(mapping,
                                        dof_handler,
                                        solution,
                                        ExactSolution<dim>(),
                                        norm_per_cell_p,
                                        QGauss<dim>(fe.degree + 2),
                                        VectorTools::H1_norm);
      std::cout << "   H1 error vs exact solution:   "
                << norm_per_cell_p.l2_norm() << std::endl;
    }

    // The final post-processing operation we do here is to compute an error
    // estimate with the KellyErrorEstimator. We use the exact same settings
    // as in the step-6 tutorial program, except for the fact that we also
    // hand in the mapping to ensure that errors are evaluated along the
    // curved element, consistent with the remainder of the program. However,
    // we do not really use the result here to drive a mesh adaptation step
    // (that would refine the mesh around the material interface along the
    // sphere), as the focus here is on the cost of this operation.
    {
      TimerOutput::Scope scope(timer, "Compute error estimator");

      Vector<float> estimated_error_per_cell(triangulation.n_active_cells());
      KellyErrorEstimator<dim>::estimate(
        mapping,
        dof_handler,
        QGauss<dim - 1>(fe.degree + 1),
        std::map<types::boundary_id, const Function<dim> *>(),
        solution,
        estimated_error_per_cell);
      std::cout << "   Max cell-wise error estimate: "
                << estimated_error_per_cell.linfty_norm() << std::endl;
    }
  }



  // @sect3{The PoissonProblem::run() function}
  //
  // Finally, we define the `run()` function that controls how we want to
  // execute this program (which is called by the main() function in the usual
  // way). We start by calling the `create_grid()` function that sets up our
  // geometry with the appropriate manifolds. We then run two instances of a
  // solver chain, starting from the setup of the equations, the assembly of
  // the linear system, its solution with a simple iterative solver, and the
  // postprocessing discussed above. The two instances differ in the way they
  // use the mapping. The first uses a conventional MappingQGeneric mapping
  // object which we initialize to a degree one more than we use for the
  // finite element &ndash; after all, we expect the geometry representation
  // to be the bottleneck as the analytic solution is only a quadratic
  // polynomial. (In reality, things are interlinked to quite some extent
  // because the evaluation of the polynomials in real coordinates involves
  // the mapping of a higher-degree polynomials, which represent some smooth
  // rational functions. As a consequence, higher-degree polynomials still pay
  // off, so it does not make sense to increase the degree of the mapping
  // further.) Once the first pass is completed, we let the timer print a
  // summary of the compute times of the individual stages.
  template <int dim>
  void PoissonProblem<dim>::run()
  {
    create_grid();

    {
      std::cout << std::endl
                << "====== Running with the basic MappingQGeneric class ====== "
                << std::endl
                << std::endl;

      MappingQGeneric<dim> mapping(fe.degree + 1);
      setup_system(mapping);
      assemble_system(mapping);
      solve();
      postprocess(mapping);

      timer.print_summary();
      timer.reset();
    }

    // For the second instance, we instead set up the MappingQCache class. Its
    // use is very simple: After constructing it (with the degree, given that
    // we want it to show the correct degree functionality in other contexts),
    // we fill the cache via the MappingQCache::initialize() function. At this
    // stage, we specify which mapping we want to use (obviously, the same
    // MappingQGeneric as previously in order to repeat the same computations)
    // for the cache, and then run through the same functions again, now
    // handing in the modified mapping. In the end, we again print the
    // accumulated wall times since the reset to see how the times compare to
    // the original setting.
    {
      std::cout
        << "====== Running with the optimized MappingQCache class ====== "
        << std::endl
        << std::endl;

      MappingQCache<dim> mapping(fe.degree + 1);
      {
        TimerOutput::Scope scope(timer, "Initialize mapping cache");
        mapping.initialize(triangulation, MappingQGeneric<dim>(fe.degree + 1));
      }
      std::cout << "   Memory consumption cache:     "
                << 1e-6 * mapping.memory_consumption() << " MB" << std::endl;

      setup_system(mapping);
      assemble_system(mapping);
      solve();
      postprocess(mapping);

      timer.print_summary();
    }
  }
} // namespace Step65



int main()
{
  Step65::PoissonProblem<3> test_program;
  test_program.run();
  return 0;
}
