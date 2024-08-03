/* ------------------------------------------------------------------------
 *
 * SPDX-License-Identifier: LGPL-2.1-or-later
 * Copyright (C) 2006 - 2024 by the deal.II authors
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
 * Authors: Wolfgang Bangerth, Texas A&M University, 2006, 2007;
 *          Denis Davydov, University of Erlangen-Nuremberg, 2016;
 *          Marc Fehling, Colorado State University, 2020.
 */


// @sect3{Include files}

// The first few files have already been covered in previous examples and will
// thus not be further commented on.
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/utilities.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/affine_constraints.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/error_estimator.h>

// These are the new files we need. The first and second provide the
// FECollection and the <i>hp</i> version of the FEValues class as described in
// the introduction of this program. The next one provides the functionality
// for automatic $hp$-adaptation, for which we will use the estimation
// algorithms based on decaying series expansion coefficients that are part of
// the last two files.
#include <deal.II/hp/fe_collection.h>
#include <deal.II/hp/fe_values.h>
#include <deal.II/hp/refinement.h>
#include <deal.II/fe/fe_series.h>
#include <deal.II/numerics/smoothness_estimator.h>

// The last set of include files are standard C++ headers.
#include <fstream>
#include <iostream>


// Finally, this is as in previous programs:
namespace Step27
{
  using namespace dealii;


  // @sect3{The main class}

  // The main class of this program looks very much like the one already used
  // in the first few tutorial programs, for example the one in step-6. The
  // main difference is that we have merged the refine_grid and output_results
  // functions into one since we will also want to output some of the
  // quantities used in deciding how to refine the mesh (in particular the
  // estimated smoothness of the solution).
  //
  // As far as member variables are concerned, we use the same structure as
  // already used in step-6, but we need collections instead of
  // individual finite element, quadrature, and face quadrature objects. We
  // will fill these collections in the constructor of the class. The last
  // variable, <code>max_degree</code>, indicates the maximal polynomial
  // degree of shape functions used.
  template <int dim>
  class LaplaceProblem
  {
  public:
    LaplaceProblem();
    ~LaplaceProblem();

    void run();

  private:
    void setup_system();
    void assemble_system();
    void solve();
    void create_coarse_grid();
    void postprocess(const unsigned int cycle);

    Triangulation<dim> triangulation;

    DoFHandler<dim>          dof_handler;
    hp::FECollection<dim>    fe_collection;
    hp::QCollection<dim>     quadrature_collection;
    hp::QCollection<dim - 1> face_quadrature_collection;

    AffineConstraints<double> constraints;

    SparsityPattern      sparsity_pattern;
    SparseMatrix<double> system_matrix;

    Vector<double> solution;
    Vector<double> system_rhs;

    const unsigned int max_degree;
  };



  // @sect3{Equation data}
  //
  // Next, let us define the right hand side function for this problem. It is
  // $x+1$ in 1d, $(x+1)(y+1)$ in 2d, and so on.
  template <int dim>
  class RightHandSide : public Function<dim>
  {
  public:
    virtual double value(const Point<dim>  &p,
                         const unsigned int component) const override;
  };


  template <int dim>
  double RightHandSide<dim>::value(const Point<dim> &p,
                                   const unsigned int /*component*/) const
  {
    double product = 1;
    for (unsigned int d = 0; d < dim; ++d)
      product *= (p[d] + 1);
    return product;
  }



  // @sect3{Implementation of the main class}

  // @sect4{LaplaceProblem::LaplaceProblem constructor}

  // The constructor of this class is fairly straightforward. It associates
  // the DoFHandler object with the triangulation, and then sets the
  // maximal polynomial degree to 7 (in 1d and 2d) or 5 (in 3d and higher). We
  // do so because using higher order polynomial degrees becomes prohibitively
  // expensive, especially in higher space dimensions.
  //
  // Following this, we fill the collections of finite element, and cell and
  // face quadrature objects. We start with quadratic elements, and each
  // quadrature formula is chosen so that it is appropriate for the matching
  // finite element in the hp::FECollection object.
  template <int dim>
  LaplaceProblem<dim>::LaplaceProblem()
    : dof_handler(triangulation)
    , max_degree(dim <= 2 ? 7 : 5)
  {
    for (unsigned int degree = 2; degree <= max_degree; ++degree)
      {
        fe_collection.push_back(FE_Q<dim>(degree));
        quadrature_collection.push_back(QGauss<dim>(degree + 1));
        face_quadrature_collection.push_back(QGauss<dim - 1>(degree + 1));
      }
  }


  // @sect4{LaplaceProblem::~LaplaceProblem destructor}

  // The destructor is unchanged from what we already did in step-6:
  template <int dim>
  LaplaceProblem<dim>::~LaplaceProblem()
  {
    dof_handler.clear();
  }


  // @sect4{LaplaceProblem::setup_system}
  //
  // This function is again a verbatim copy of what we already did in
  // step-6. Despite function calls with exactly the same names and arguments,
  // the algorithms used internally are different in some aspect since the
  // dof_handler variable here is in $hp$-mode.
  template <int dim>
  void LaplaceProblem<dim>::setup_system()
  {
    dof_handler.distribute_dofs(fe_collection);

    solution.reinit(dof_handler.n_dofs());
    system_rhs.reinit(dof_handler.n_dofs());

    constraints.clear();
    DoFTools::make_hanging_node_constraints(dof_handler, constraints);
    VectorTools::interpolate_boundary_values(dof_handler,
                                             0,
                                             Functions::ZeroFunction<dim>(),
                                             constraints);
    constraints.close();

    DynamicSparsityPattern dsp(dof_handler.n_dofs(), dof_handler.n_dofs());
    DoFTools::make_sparsity_pattern(dof_handler, dsp, constraints, false);
    sparsity_pattern.copy_from(dsp);

    system_matrix.reinit(sparsity_pattern);
  }



  // @sect4{LaplaceProblem::assemble_system}

  // This is the function that assembles the global matrix and right hand side
  // vector from the local contributions of each cell. Its main working is as
  // has been described in many of the tutorial programs before. The
  // significant deviations are the ones necessary for <i>hp</i> finite
  // element methods. In particular, that we need to use a collection of
  // FEValues object (implemented through the hp::FEValues class), and that we
  // have to eliminate constrained degrees of freedom already when copying
  // local contributions into global objects. Both of these are explained in
  // detail in the introduction of this program.
  //
  // One other slight complication is the fact that because we use different
  // polynomial degrees on different cells, the matrices and vectors holding
  // local contributions do not have the same size on all cells. At the
  // beginning of the loop over all cells, we therefore each time have to
  // resize them to the correct size (given by <code>dofs_per_cell</code>).
  // Because these classes are implemented in such a way that reducing the size
  // of a matrix or vector does not release the currently allocated memory
  // (unless the new size is zero), the process of resizing at the beginning of
  // the loop will only require re-allocation of memory during the first few
  // iterations. Once we have found in a cell with the maximal finite element
  // degree, no more re-allocations will happen because all subsequent
  // <code>reinit</code> calls will only set the size to something that fits the
  // currently allocated memory. This is important since allocating memory is
  // expensive, and doing so every time we visit a new cell would take
  // significant compute time.
  template <int dim>
  void LaplaceProblem<dim>::assemble_system()
  {
    hp::FEValues<dim> hp_fe_values(fe_collection,
                                   quadrature_collection,
                                   update_values | update_gradients |
                                     update_quadrature_points |
                                     update_JxW_values);

    RightHandSide<dim> rhs_function;

    FullMatrix<double> cell_matrix;
    Vector<double>     cell_rhs;

    std::vector<types::global_dof_index> local_dof_indices;

    for (const auto &cell : dof_handler.active_cell_iterators())
      {
        const unsigned int dofs_per_cell = cell->get_fe().n_dofs_per_cell();

        cell_matrix.reinit(dofs_per_cell, dofs_per_cell);
        cell_matrix = 0;

        cell_rhs.reinit(dofs_per_cell);
        cell_rhs = 0;

        hp_fe_values.reinit(cell);

        const FEValues<dim> &fe_values = hp_fe_values.get_present_fe_values();

        std::vector<double> rhs_values(fe_values.n_quadrature_points);
        rhs_function.value_list(fe_values.get_quadrature_points(), rhs_values);

        for (unsigned int q_point = 0; q_point < fe_values.n_quadrature_points;
             ++q_point)
          for (unsigned int i = 0; i < dofs_per_cell; ++i)
            {
              for (unsigned int j = 0; j < dofs_per_cell; ++j)
                cell_matrix(i, j) +=
                  (fe_values.shape_grad(i, q_point) * // grad phi_i(x_q)
                   fe_values.shape_grad(j, q_point) * // grad phi_j(x_q)
                   fe_values.JxW(q_point));           // dx

              cell_rhs(i) += (fe_values.shape_value(i, q_point) * // phi_i(x_q)
                              rhs_values[q_point] *               // f(x_q)
                              fe_values.JxW(q_point));            // dx
            }

        local_dof_indices.resize(dofs_per_cell);
        cell->get_dof_indices(local_dof_indices);

        constraints.distribute_local_to_global(
          cell_matrix, cell_rhs, local_dof_indices, system_matrix, system_rhs);
      }
  }



  // @sect4{LaplaceProblem::solve}

  // The function solving the linear system is entirely unchanged from
  // previous examples. We simply try to reduce the initial residual (which
  // equals the $l_2$ norm of the right hand side) by a certain factor:
  template <int dim>
  void LaplaceProblem<dim>::solve()
  {
    SolverControl            solver_control(system_rhs.size(),
                                 1e-12 * system_rhs.l2_norm());
    SolverCG<Vector<double>> cg(solver_control);

    PreconditionSSOR<SparseMatrix<double>> preconditioner;
    preconditioner.initialize(system_matrix, 1.2);

    cg.solve(system_matrix, solution, system_rhs, preconditioner);

    constraints.distribute(solution);
  }



  // @sect4{LaplaceProblem::postprocess}

  // After solving the linear system, we will want to postprocess the
  // solution. Here, all we do is to estimate the error, estimate the local
  // smoothness of the solution as described in the introduction, then write
  // graphical output, and finally refine the mesh in both $h$ and $p$
  // according to the indicators computed before. We do all this in the same
  // function because we want the estimated error and smoothness indicators
  // not only for refinement, but also include them in the graphical output.
  template <int dim>
  void LaplaceProblem<dim>::postprocess(const unsigned int cycle)
  {
    // Let us start with computing estimated error and smoothness indicators,
    // which each are one number for each active cell of our
    // triangulation. For the error indicator, we use the KellyErrorEstimator
    // class as always.
    Vector<float> estimated_error_per_cell(triangulation.n_active_cells());
    KellyErrorEstimator<dim>::estimate(
      dof_handler,
      face_quadrature_collection,
      std::map<types::boundary_id, const Function<dim> *>(),
      solution,
      estimated_error_per_cell);

    // Estimating the smoothness is performed with the method of decaying
    // expansion coefficients as outlined in the introduction. We will first
    // need to create an object capable of transforming the finite element
    // solution on every single cell into a sequence of Fourier series
    // coefficients. The SmoothnessEstimator namespace offers a factory function
    // for such a FESeries::Fourier object that is optimized for the process of
    // estimating smoothness. The actual determination of the decay of Fourier
    // coefficients on every individual cell then happens in the last function.
    Vector<float> smoothness_indicators(triangulation.n_active_cells());
    FESeries::Fourier<dim> fourier =
      SmoothnessEstimator::Fourier::default_fe_series(fe_collection);
    SmoothnessEstimator::Fourier::coefficient_decay(fourier,
                                                    dof_handler,
                                                    solution,
                                                    smoothness_indicators);

    // Next we want to generate graphical output. In addition to the two
    // estimated quantities derived above, we would also like to output the
    // polynomial degree of the finite elements used on each of the elements
    // on the mesh.
    //
    // The way to do that requires that we loop over all cells and poll the
    // active finite element index of them using
    // <code>cell-@>active_fe_index()</code>. We then use the result of this
    // operation and query the finite element collection for the finite
    // element with that index, and finally determine the polynomial degree of
    // that element. The result we put into a vector with one element per
    // cell. The DataOut class requires this to be a vector of
    // <code>float</code> or <code>double</code>, even though our values are
    // all integers, so that is what we use:
    {
      Vector<float> fe_degrees(triangulation.n_active_cells());
      for (const auto &cell : dof_handler.active_cell_iterators())
        fe_degrees(cell->active_cell_index()) =
          fe_collection[cell->active_fe_index()].degree;

      // With now all data vectors available -- solution, estimated errors and
      // smoothness indicators, and finite element degrees --, we create a
      // DataOut object for graphical output and attach all data:
      DataOut<dim> data_out;

      data_out.attach_dof_handler(dof_handler);
      data_out.add_data_vector(solution, "solution");
      data_out.add_data_vector(estimated_error_per_cell, "error");
      data_out.add_data_vector(smoothness_indicators, "smoothness");
      data_out.add_data_vector(fe_degrees, "fe_degree");
      data_out.build_patches();

      // The final step in generating output is to determine a file name, open
      // the file, and write the data into it (here, we use VTK format):
      const std::string filename =
        "solution-" + Utilities::int_to_string(cycle, 2) + ".vtk";
      std::ofstream output(filename);
      data_out.write_vtk(output);
    }

    // After this, we would like to actually refine the mesh, in both $h$ and
    // $p$. The way we are going to do this is as follows: first, we use the
    // estimated error to flag those cells for refinement that have the
    // largest error. This is what we have always done:
    {
      GridRefinement::refine_and_coarsen_fixed_number(triangulation,
                                                      estimated_error_per_cell,
                                                      0.3,
                                                      0.03);

      // Next we would like to figure out which of the cells that have been
      // flagged for refinement should actually have $p$ increased instead of
      // $h$ decreased. The strategy we choose here is that we look at the
      // smoothness indicators of those cells that are flagged for refinement,
      // and increase $p$ for those with a smoothness larger than a certain
      // relative threshold. In other words, for every cell for which (i) the
      // refinement flag is set, (ii) the smoothness indicator is larger than
      // the threshold, and (iii) we still have a finite element with a
      // polynomial degree higher than the current one in the finite element
      // collection, we will assign a future FE index that corresponds to a
      // polynomial with degree one higher than it currently is. The following
      // function is capable of doing exactly this. Absent any better
      // strategies, we will set the threshold via interpolation between the
      // minimal and maximal smoothness indicators on cells flagged for
      // refinement. Since the corner singularities are strongly localized, we
      // will favor $p$- over $h$-refinement quantitatively. We achieve this
      // with a low threshold by setting a small interpolation factor of 0.2. In
      // the same way, we deal with cells that are going to be coarsened and
      // decrease their polynomial degree when their smoothness indicator is
      // below the corresponding threshold determined on cells to be coarsened.
      hp::Refinement::p_adaptivity_from_relative_threshold(
        dof_handler, smoothness_indicators, 0.2, 0.2);

      // The above function only determines whether the polynomial degree will
      // change via future FE indices, but does not manipulate the
      // $h$-refinement flags. So for cells that are flagged for both refinement
      // categories, we prefer $p$- over $h$-refinement. The following function
      // call ensures that only one of $p$- or $h$-refinement is imposed, and
      // not both at once.
      hp::Refinement::choose_p_over_h(dof_handler);

      // For grid adaptive refinement, we ensure a 2:1 mesh balance by limiting
      // the difference of refinement levels of neighboring cells to one by
      // calling Triangulation::prepare_coarsening_and_refinement(). We would
      // like to achieve something similar for the p-levels of neighboring
      // cells: levels of future finite elements are not allowed to differ by
      // more than a specified difference. With its default parameters, a call
      // of hp::Refinement::limit_p_level_difference() ensures that their level
      // difference is limited to one. This will not necessarily decrease the
      // number of hanging nodes in the domain, but makes sure that high order
      // polynomials are not constrained to much lower polynomials on faces,
      // e.g., fifth order to second order polynomials.
      triangulation.prepare_coarsening_and_refinement();
      hp::Refinement::limit_p_level_difference(dof_handler);

      // At the end of this procedure, we then refine the mesh. During this
      // process, children of cells undergoing bisection inherit their mother
      // cell's finite element index. Further, future finite element indices
      // will turn into active ones, so that the new finite elements will be
      // assigned to cells after the next call of DoFHandler::distribute_dofs().
      triangulation.execute_coarsening_and_refinement();
    }
  }


  // @sect4{LaplaceProblem::create_coarse_grid}

  // The following function is used when creating the initial grid. The grid we
  // would like to create is actually similar to the one from step-14, i.e., the
  // square domain with the square hole in the middle. It can be generated by
  // exactly the same function. However, since its implementation is only a
  // specialization of the 2d case, we will present a different way of creating
  // this domain which is dimension independent.
  //
  // We first create a hypercube triangulation with enough cells so that it
  // already holds our desired domain $[-1,1]^d$, subdivided into $4^d$ cells.
  // We then remove those cells in the center of the domain by testing the
  // coordinate values of the vertices on each cell. In the end, we refine the
  // so created grid globally as usual.
  template <int dim>
  void LaplaceProblem<dim>::create_coarse_grid()
  {
    Triangulation<dim> cube;
    GridGenerator::subdivided_hyper_cube(cube, 4, -1., 1.);

    std::set<typename Triangulation<dim>::active_cell_iterator> cells_to_remove;
    for (const auto &cell : cube.active_cell_iterators())
      for (unsigned int v = 0; v < GeometryInfo<dim>::vertices_per_cell; ++v)
        if (cell->vertex(v).square() < .1)
          cells_to_remove.insert(cell);

    GridGenerator::create_triangulation_with_removed_cells(cube,
                                                           cells_to_remove,
                                                           triangulation);

    triangulation.refine_global(3);
  }



  // @sect4{LaplaceProblem::run}

  // This function implements the logic of the program, as did the respective
  // function in most of the previous programs already, see for example step-6.
  //
  // Basically, it contains the adaptive loop: in the first iteration create a
  // coarse grid, and then set up the linear system, assemble it, solve, and
  // postprocess the solution including mesh refinement. Then start over
  // again. In the meantime, also output some information for those staring at
  // the screen trying to figure out what the program does:
  template <int dim>
  void LaplaceProblem<dim>::run()
  {
    for (unsigned int cycle = 0; cycle < 6; ++cycle)
      {
        std::cout << "Cycle " << cycle << ':' << std::endl;

        if (cycle == 0)
          create_coarse_grid();

        setup_system();

        std::cout << "   Number of active cells      : "
                  << triangulation.n_active_cells() << std::endl
                  << "   Number of degrees of freedom: " << dof_handler.n_dofs()
                  << std::endl
                  << "   Number of constraints       : "
                  << constraints.n_constraints() << std::endl;

        assemble_system();
        solve();
        postprocess(cycle);
      }
  }
} // namespace Step27


// @sect3{The main function}

// The main function is again verbatim what we had before: wrap creating and
// running an object of the main class into a <code>try</code> block and catch
// whatever exceptions are thrown, thereby producing meaningful output if
// anything should go wrong:
int main()
{
  try
    {
      using namespace Step27;

      LaplaceProblem<2> laplace_problem;
      laplace_problem.run();
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
