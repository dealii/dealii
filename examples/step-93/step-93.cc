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
 * Authors : Sam Scheuerman, Wolfgang Bangerth, Colorado State University, 2024.
 */


// @sect3{Include files and other top matter}

#include <deal.II/base/function.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/timer.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_enriched.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/hp/fe_collection.h>
#include <deal.II/hp/fe_values.h>

#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_minres.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>

#include <fstream>
#include <iostream>


namespace Step93
{
  using namespace dealii;


  // @sect3{Target and indicator functions}

  // We start by defining a function class for $\bar u$; it represents the heat
  // profile we want to match.
  //
  // This class has two member variables:
  //
  //   - `center`: A Point object representing the center of the function
  //   - `radius`: A double representing the radius of the circular step, or the
  //   standard deviation of a Gaussian
  template <int dim>
  class TargetFunction : public Function<dim>
  {
  public:
    TargetFunction()
      : Function<dim>(1)
      , center()
      , radius(0.1)
    {}

    // Next, we define an overloaded constructor for TargetFunction, so we can
    // choose to set the center and radius
    //
    // Parameters:
    //
    //   - `center`: A constant Point pointer used to set the member variable
    //   center
    //   - `radius`: A constant double used to set the member variable radius
    TargetFunction(const unsigned int n_components,
                   const Point<dim>  &center,
                   const double       radius = .3)
      : Function<dim>(n_components)
      , center(center)
      , radius(radius)
    {}

    virtual double value(const Point<dim>  &p,
                         const unsigned int component = 0) const override;

  private:
    const Point<dim> center;
    const double     radius;
  };

  // The value() function returns the value of the function at point `p` and
  // component index `component`. In this case, if the component corresponds to
  // the solution $u$, then the function returns a value based on either a step
  // function or a Gaussian. If the component corresponds to any of the other
  // variables, the function always returns 0.
  //
  // Note: In the documentation we discuss using a Gaussian target function, but
  // this is not included in the code for readability. To run the code with a
  // Gaussian target function, simply replace the code in the `component == 0`
  // braces by
  //
  // `return std::exp(-((p-center).norm()*(p-center).norm())/(radius*radius));`
  template <int dim>
  double TargetFunction<dim>::value(const Point<dim>  &p,
                                    const unsigned int component) const
  {
    if (component == 0)
      {
        if ((p - center).norm() <= radius)
          return 1;
        else
          return 0;
      }
    else
      return 0;
  }

  // The next class we define is one for a circular indicator function. Unlike
  // the target function, this does not need a component argument because we
  // have to manually address where it gets used in the code. Objects of this
  // function type correspond to the nonlocal dofs.
  //
  //  Parameters:
  //
  //   - `center`: A constant Point object giving the center of the indicator
  //   region.
  //   - `radius`: The radius of the region.
  template <int dim>
  class CircularIndicatorFunction : public Function<dim>
  {
  public:
    CircularIndicatorFunction();

    CircularIndicatorFunction(const Point<dim> &center, const double radius);

    virtual double value(const Point<dim>  &p,
                         const unsigned int component = 0) const override;

  private:
    const Point<dim> center;
    const double     radius;
  };

  template <int dim>
  CircularIndicatorFunction<dim>::CircularIndicatorFunction()
    : center(Point<dim>())
    , radius(20)
  {}


  template <int dim>
  CircularIndicatorFunction<dim>::CircularIndicatorFunction(
    const Point<dim> &center,
    const double      radius)
    : center(center)
    , radius(radius)
  {}

  template <int dim>
  double CircularIndicatorFunction<dim>::value(
    const Point<dim> &p,
    const unsigned int /* component */) const
  {
    if ((center - p).norm() <= radius)
      return 1;
    else
      return 0;
  }


  // @sect3{The principal class}

  // The main class is very similar to step-4 in structure, given that this is
  // a relatively simple program that does not use adaptive mesh refinement.
  // However, there are three new member variables:
  //
  //   - `nonlocal_dofs`: A `std::vector` of dof indices that stores the
  //   dof index for the nonlocal dofs.
  //
  //   - `heat_functions`: A `std::vector` of CircularIndicatorFunction objects;
  //   these are the heat sources.
  //
  //   - `target_function`: This is the function we want to match. We store it
  //   as a class variable because it is used both in assemble_system() and
  //   output_results().
  template <int dim>
  class Step93
  {
  public:
    Step93();

    void run();

  private:
    void make_grid();
    void setup_system();
    void assemble_system();
    void solve();
    void output_results() const;

    Triangulation<dim> triangulation;
    DoFHandler<dim>    dof_handler;

    hp::FECollection<dim> fe_collection;
    hp::QCollection<dim>  quadrature_collection;

    AffineConstraints<double> constraints;

    SparsityPattern      sparsity_pattern;
    SparseMatrix<double> system_matrix;

    Vector<double> solution;
    Vector<double> system_rhs;

    std::vector<types::global_dof_index> nonlocal_dofs;

    std::vector<CircularIndicatorFunction<dim>> heat_functions;

    const TargetFunction<dim> target_function;
  };


  // @sect4{The Step93 constructor}
  // In this constructor, we set up several of the fundamental data
  // structures this program needs. Specifically, we generate the
  // finite element collection, which is basically a list of all the
  // possible finite elements we could use on each cell. This
  // collection has two elements: one FESystem that has two degree 2
  // FE_Q elements and one FE_Nothing element (to be used on all cells
  // that are not "special"), and one FESystem that has two degree 2
  // FE_Q elements and one degree 0 FE_DGQ element (to be used on
  // those "special" cells we use to anchor the non-local degrees of
  // freedom, as discussed in the introduction).
  //
  // Where we have a collection of elements, we then also need a
  // collection of quadratures -- which here has only a single element
  // because we can use the same quadrature on all cells.

  // The default constructor also  initializes
  // dof_handler and target_function,
  // and it generates a vector of CircularIndicatorFunctions objects which
  // represent the heat functions.
  template <int dim>
  Step93<dim>::Step93()
    : dof_handler(triangulation)
    , target_function(3,
                      (dim == 1 ? Point<dim>(0.5) :
                       dim == 2 ? Point<dim>(0.5, 0.5) :
                                  Point<dim>(0.5, 0.5, 0.5)))
  {
    // Here, we generate the finite element collection, which is basically a
    // list of all the possible finite elements we could use on each cell. This
    // collection has two elements: one FESystem that has two degree 2 FE_Q
    // elements and one FE_Nothing element, and one FESystem that has two
    // degree 2 FE_Q elements and one degree 0 FE_DGQ element.
    fe_collection.push_back(
      FESystem<dim>(FE_Q<dim>(2), 2, FE_Nothing<dim>(), 1));
    fe_collection.push_back(FESystem<dim>(FE_Q<dim>(2), 2, FE_DGQ<dim>(0), 1));

    // The quadrature collection is just one degree 3 QGauss element.
    quadrature_collection.push_back(QGauss<dim>(3));

    // Here, we choose the center points for the heat functions to be the
    // vertices of a lattice, in this case the corners of a hypercube
    // centered at 0, with side length 1. The block of code below
    // then creates the CircularIndicatorFunction
    // objects for the main class. To make these, we check the dimension of the
    // problem and create 2, 4, or 8 function objects.
    switch (dim)
      {
        case (1):
          heat_functions.emplace_back(Point<dim>({0.5}), 0.2);
          heat_functions.emplace_back(Point<dim>({-0.5}), 0.2);
          break;
        case (2):
          heat_functions.emplace_back(Point<dim>({0.5, 0.5}), 0.2);
          heat_functions.emplace_back(Point<dim>({0.5, -0.5}), 0.2);
          heat_functions.emplace_back(Point<dim>({-0.5, 0.5}), 0.2);
          heat_functions.emplace_back(Point<dim>({-0.5, -0.5}), 0.2);
          break;
        case (3):
          heat_functions.emplace_back(Point<dim>({0.5, 0.5, 0.5}), 0.2);
          heat_functions.emplace_back(Point<dim>({0.5, 0.5, -0.5}), 0.2);
          heat_functions.emplace_back(Point<dim>({0.5, -0.5, 0.5}), 0.2);
          heat_functions.emplace_back(Point<dim>({0.5, -0.5, -0.5}), 0.2);
          heat_functions.emplace_back(Point<dim>({-0.5, 0.5, 0.5}), 0.2);
          heat_functions.emplace_back(Point<dim>({-0.5, 0.5, -0.5}), 0.2);
          heat_functions.emplace_back(Point<dim>({-0.5, -0.5, 0.5}), 0.2);
          heat_functions.emplace_back(Point<dim>({-0.5, -0.5, -0.5}), 0.2);
          break;
        default:
          DEAL_II_ASSERT_UNREACHABLE();
          break;
      }
  }

  // @sect4{Step93::make_grid()}

  // The make_grid() function makes a hypercube grid, see step-4:
  template <int dim>
  void Step93<dim>::make_grid()
  {
    GridGenerator::hyper_cube(triangulation, -1, 1);
    triangulation.refine_global(7);

    std::cout << "Number of active cells: " << triangulation.n_active_cells()
              << std::endl;
  }


  // @sect4{Step93::setup_system()}

  // The `setup_system()` function is similar to step-4, except we have to add a
  // few steps to prepare for the nonlocal dofs.
  template <int dim>
  void Step93<dim>::setup_system()
  {
    // Here, we loop over the cells and set the FESystem index to 1, which
    // corresponds to the system with 2 FE_Q elements and one FE_DGQ element. We
    // do this until we have enough dofs for each heat function. Note that we
    // use the global active cell index to measure when to stop the
    // loop. This allows the loop to be run in parallel with no alteration.
    // Then, we call DoFHandler::distribute_dofs() to actually enumerate all
    // degrees of freedom.
    for (const auto &cell : dof_handler.active_cell_iterators())
      {
        if (cell->global_active_cell_index() < heat_functions.size())
          {
            cell->set_active_fe_index(1);
          }
        else
          {
            break;
          }
      }
    dof_handler.distribute_dofs(fe_collection);

    // Once we've assigned dofs, the code block below counts the
    // number of dofs in the system, and outputs to the console. In
    // other contexts, we might want to use *block* matrices (see, for
    // example, step-20 or step-22) to build more efficient linear
    // solvers; here, we will just put everything into one big matrix
    // and so knowing the number of unknowns for each of the variables
    // is purely for informational purposes.
    const std::vector<types::global_dof_index> dofs_per_component =
      DoFTools::count_dofs_per_fe_component(dof_handler);
    const unsigned int n_dofs_u = dofs_per_component[0],
                       n_dofs_l = dofs_per_component[1],
                       n_dofs_c = dofs_per_component[2];
    std::cout << "Number of degrees of freedom: " << n_dofs_u << "+" << n_dofs_l
              << "+" << n_dofs_c << " = " << n_dofs_u + n_dofs_l + n_dofs_c
              << std::endl;

    // Finally, we need to extract the indices of the finite elements which
    // correspond to the non-local dofs.
    //
    // First, we make a component mask, which is `false` except for
    // the third component. This will extract only the dofs from the
    // third component of the FE system. Next, we actually extract the
    // dofs, and store them in an IndexSet variable.  Finally, we add
    // each extracted index to the member array `nonlocal_dofs`.
    const ComponentMask component_mask_c({false, false, true});
    const IndexSet      indices_c =
      DoFTools::extract_dofs(dof_handler, component_mask_c);

    for (const types::global_dof_index non_local_index : indices_c)
      nonlocal_dofs.push_back(non_local_index);

    std::cout << "Number of nonlocal dofs: " << nonlocal_dofs.size()
              << std::endl;


    // The mesh we created above is not locally refined, and so there
    // are no hanging node constraints to keep track of. But it does
    // not hurt to just use the same setup we have used starting in
    // step-6 of building a constraints object that contains hanging
    // node constraints (anticipating that perhaps we'd want to do
    // adaptive mesh refinement in a later step) into which we then
    // also put the constraints for boundary values on both the $u$
    // and $\lambda$ variables.
    //
    // Because the nonlocal degrees of freedom use discontinuous
    // elements, they do not contribute to boundary values
    // (discontinuous elements do not have degrees of freedom
    // logically located on the boundary that could be interpolated)
    // and we do not need to exclude these solution components
    // explicitly when calling
    // VectorTools::interpolate_boundary_values()
    constraints.clear();
    DoFTools::make_hanging_node_constraints(dof_handler, constraints);
    VectorTools::interpolate_boundary_values(dof_handler,
                                             0,
                                             Functions::ZeroFunction<dim>(3),
                                             constraints);
    constraints.close();

    // The remainder of the function deals with building the sparsity
    // pattern. It consists of two parts: The entries that result from
    // the usual integration of products of shape functions, and then
    // the entries that result from integrals that contain nonlocal
    // degrees of freedom (which one can think of as associated with
    // shape functions that are constant across the entire
    // domain). The first part is easily built using standard tools:
    DynamicSparsityPattern dsp(dof_handler.n_dofs(), dof_handler.n_dofs());
    DoFTools::make_sparsity_pattern(dof_handler, dsp, constraints, false);

    // The other part is more awkward. We have matrix entries that
    // result from terms such as $\int_{\Omega}\varphi_j f_k$ where
    // each $\varphi_j$ is a shape function associated to $\lambda$
    // and each $f_k$ is a
    // characteristic function of a part of the domain. As a
    // consequence, we end up with a nonzero matrix entry $A_{jk}$
    // (along with its transpose $A_{kj}$) if there is overlap between
    // $\varphi_j$ and $f_k$. In practice, because we will use
    // quadrature, this means that we end up with a quadrature point
    // on a cell on which $\varphi_j$ lives and at which $f_k$ is not
    // zero. (We will implicitly assume that a shape function that
    // lives on the current cell is nonzero at all quadrature points
    // -- an assumption that is generally true unless one chooses
    // specific quadrature formulas.) Determining which sparsity
    // pattern entries we need to add then essentially comes down to
    // "simulating" what would happen if we actually computed
    // integrals and which matrix entries would end up being non-zero
    // in the process. As a consequence, the following code's general
    // structure looks very similar to what we will do for the
    // nonlocal contributions in the `assemble_system()` function
    // below.
    //
    // To get this started, we create an `hp_fe_values` that we
    // will only use to query the quadrature point locations.  The
    // non-local dofs will need to interact with the second component
    // of the fe system (namely, $\lambda$), so we also declare a
    // variable that will help us extract this scalar field for use
    // below.
    hp::FEValues<dim> hp_fe_values(fe_collection,
                                   quadrature_collection,
                                   update_quadrature_points);

    // Then, we loop over the cells, then over the quadrature points, and
    // finally over the indices, as if we were constructing a mass matrix.
    // However, what we instead do here is check two things. First, we check if
    // the quadrature point is within the radius of a circular indicator
    // function that represents our non-local dof. If so
    // then we add an entry to the sparse matrix at the
    // (nonlocal dof index, lambda dof index) entry and the (lambda dof index,
    // nonlocal dof index) entry for all lambda degrees of freedom. (Because the
    // matrix we solve with has both the lambda-nonlocal interacting block and
    // its transpose, we need to add two entries each time.)
    std::vector<types::global_dof_index> local_dof_indices;
    for (const auto &cell : dof_handler.active_cell_iterators())
      {
        hp_fe_values.reinit(cell);

        const FEValues<dim> &fe_values = hp_fe_values.get_present_fe_values();

        local_dof_indices.resize(fe_values.dofs_per_cell);
        cell->get_dof_indices(local_dof_indices);

        for (const unsigned int q_index : fe_values.quadrature_point_indices())
          {
            const Point<dim> q_point = fe_values.quadrature_point(q_index);
            for (const unsigned int i : fe_values.dof_indices())
              {
                if (fe_values.get_fe().system_to_component_index(i).first ==
                    1) // 'i' is a lambda shape function
                  {
                    for (unsigned int j = 0; j < heat_functions.size(); ++j)
                      if (heat_functions[j].value(q_point) != 0)
                        {
                          dsp.add(local_dof_indices[i], nonlocal_dofs[j]);
                          dsp.add(nonlocal_dofs[j], local_dof_indices[i]);
                        }
                  }
              }
          }
      }

    // The rest (below) is standard setup code, see step-4:
    sparsity_pattern.copy_from(dsp);

    system_matrix.reinit(sparsity_pattern);

    solution.reinit(dof_handler.n_dofs());
    system_rhs.reinit(dof_handler.n_dofs());
  }


  // @sect4{Step93::assemble_system()}

  // The `assemble_system()` function works very similar to how is
  // does in other tutorial programs (cf. step-4, step-6, step-8, and
  // for the vector-valued case see step-22).  However, there is an
  // additional component to constructing the system matrix, because
  // we need to handle the nonlocal dofs manually.
  template <int dim>
  void Step93<dim>::assemble_system()
  {
    // First, we do a standard loop setup for constructing the system matrix.
    hp::FEValues<dim> hp_fe_values(fe_collection,
                                   quadrature_collection,
                                   update_values | update_gradients |
                                     update_quadrature_points |
                                     update_JxW_values);

    FullMatrix<double> cell_matrix;
    Vector<double>     cell_rhs;

    std::vector<types::global_dof_index> local_dof_indices;

    const FEValuesExtractors::Scalar u(0);
    const FEValuesExtractors::Scalar lambda(1);

    for (const auto &cell : dof_handler.active_cell_iterators())
      {
        const unsigned int dofs_per_cell = cell->get_fe().n_dofs_per_cell();

        cell_matrix.reinit(dofs_per_cell, dofs_per_cell);
        cell_rhs.reinit(dofs_per_cell);
        hp_fe_values.reinit(cell);

        cell_matrix = 0;
        cell_rhs    = 0;

        const FEValues<dim> &fe_values = hp_fe_values.get_present_fe_values();

        local_dof_indices.resize(fe_values.dofs_per_cell);

        cell->get_dof_indices(local_dof_indices);


        // In the loop over quadrature points, we start by building
        // all of the usual terms that are bilinear in shape functions
        // corresponding to the $u$ and $\lambda$ variables:
        for (const unsigned int q_index : fe_values.quadrature_point_indices())
          {
            const double JxW = fe_values.JxW(q_index);
            for (const unsigned int i : fe_values.dof_indices())
              {
                const double phi_i_u = fe_values[u].value(i, q_index),
                             phi_i_l = fe_values[lambda].value(i, q_index);

                const Tensor<1, dim> grad_i_u =
                                       fe_values[u].gradient(i, q_index),
                                     grad_i_l =
                                       fe_values[lambda].gradient(i, q_index);

                for (const unsigned int j : fe_values.dof_indices())
                  {
                    const double phi_j_u = fe_values[u].value(j, q_index);

                    const Tensor<1, dim> grad_j_u =
                                           fe_values[u].gradient(j, q_index),
                                         grad_j_l =
                                           fe_values[lambda].gradient(j,
                                                                      q_index);

                    cell_matrix(i, j) += phi_i_u * phi_j_u * JxW;
                    cell_matrix(i, j) += -grad_i_u * grad_j_l * JxW;
                    cell_matrix(i, j) += -grad_i_l * grad_j_u * JxW;
                  }

                const Point<dim> q_point = fe_values.quadrature_point(q_index);
                cell_rhs(i) += (phi_i_u * target_function.value(q_point) * JxW);


                // For the integrals that involve the nonlocal dofs,
                // we make use of the quadrature point again. To
                // compute the integrals, we loop over the
                // heat functions, adding the numeric integral of each
                // heat equation with each $\lambda$ shape function,
                // at the appropriate indices (which we found in
                // `setup_system()`). Note that if we try to add 0 to
                // a matrix entry we have not previously indicated should
                // be nonzero, there will not be a problem; but if we
                // try to add a nonzero value to an entry not
                // previously added to the sparsity pattern, we will
                // get an error. In other words, the following lines
                // of the code check that we adjusted the sparsity
                // pattern correctly in the previous function.
                for (unsigned int j = 0; j < heat_functions.size(); ++j)
                  {
                    system_matrix.add(local_dof_indices[i],
                                      nonlocal_dofs[j],
                                      heat_functions[j].value(q_point) *
                                        phi_i_l * JxW);
                    system_matrix.add(nonlocal_dofs[j],
                                      local_dof_indices[i],
                                      heat_functions[j].value(q_point) *
                                        phi_i_l * JxW);
                  }
              }
          }

        // Finally, we copy the local contributions to the linear
        // system into the global matrix and right hand side vector,
        // taking into account hanging node and boundary values
        // constraints:
        constraints.distribute_local_to_global(
          cell_matrix, cell_rhs, local_dof_indices, system_matrix, system_rhs);
      }
  }


  // @sect4{Step93::solve()}

  // The solve() function works similar to how it is done in step-6
  // and step-8, except we need to use a different solver because the
  // linear problem we are trying to solve is a saddle point problem
  // for which the Conjugate Gradient algorithm is not
  // applicable. But, because the matrix is symmetric, we can use
  // SolverMinRes, an iterative solver specialized for symmetric
  // indefinite problems. This solver could be improved with the use
  // of preconditioners, but we don't do that here for simplicity (see
  // the Possibilities for Extensions section below).
  //
  // As you will see in the output, given that we are not using a
  // preconditioner, we need a *lot* of iterations to solve this
  // linear system. We set the maximum to one million, more than we
  // need of course, but an indication that this is not an efficient
  // solver. For smaller problems, one can also use a direct solver
  // (see step-29) for which you would just replace the main part of
  // this function by the following three lines of code:
  // @code
  //   SparseDirectUMFPACK direct_solver;
  //   direct_solver.initialize(system_matrix);
  //   direct_solver.vmult(solution, system_rhs);
  // @endcode
  template <int dim>
  void Step93<dim>::solve()
  {
    Timer timer;
    timer.start();

    std::cout << "Beginning solve..." << std::endl;

    SolverControl solver_control(1'000'000, 1e-6 * system_rhs.l2_norm());
    SolverMinRes<Vector<double>> solver(solver_control);

    solver.solve(system_matrix, solution, system_rhs, PreconditionIdentity());

    timer.stop();

    std::cout << "Wall time: " << timer.wall_time() << "s" << std::endl;
    std::cout << "Solved in " << solver_control.last_step()
              << " MINRES iterations." << std::endl;
  }


  // @sect4{Step93::output_results()}

  // The `output_results()` function is a bit more robust for this program than
  // is typical. This is because, in order to visualize the heat sources we have
  // optimized, we need to do extra work and interpolate them onto a mesh. We do
  // this by instantiating a new DoFHandler object and then using the helper
  // function VectorTools::interpolate().
  //
  // The top of the function is as always when using vector-valued
  // elements (see, for example, step-22) and simply outputs all of
  // the solution variables on the mesh cells they are defined on:
  template <int dim>
  void Step93<dim>::output_results() const
  {
    const std::vector<std::string> solution_names = {"u", "lambda", "c"};

    const std::vector<DataComponentInterpretation::DataComponentInterpretation>
      interpretation = {DataComponentInterpretation::component_is_scalar,
                        DataComponentInterpretation::component_is_scalar,
                        DataComponentInterpretation::component_is_scalar};

    DataOut<dim> data_out;
    data_out.add_data_vector(dof_handler,
                             solution,
                             solution_names,
                             interpretation);

    // The non-local degrees of freedom are of course defined on the
    // first several cells, but at least logically are considered to
    // live everywhere (or, if you prefer, nowhere at all since they
    // do not represent spatial functions). But, conceptually, we use
    // them as multipliers for the heat sources, and so while the
    // coefficients $C^k$ are non-local, the heat source $\sum_k C^k
    // f_k(\mathbf x)$ *is* a spatially variable function. It would be
    // nice if we could visualize that as well. The same is true for
    // the target heat distribution $\bar u$ we are trying to match.
    //
    // To do so, we
    // create a new dof handler to output the target function
    // and heat plate values, and associate it with
    // a finite element with a degree that matches what we used to solve for $u$
    // and $\lambda$, although in reality this is an arbitrary choice:
    DoFHandler<dim> new_dof_handler(triangulation);

    const FE_Q<dim> new_fe(2);
    new_dof_handler.distribute_dofs(new_fe);

    // To get started with the visualization, we need a vector which
    // stores the interpolated target function. We create the vector,
    // interpolate the target function $\bar u$ onto the mesh, then
    // add the data to our data_out object.
    Vector<double> target(new_dof_handler.n_dofs());
    VectorTools::interpolate(new_dof_handler,
                             ScalarFunctionFromFunctionObject<dim, double>(
                               [&](const Point<dim> &x) {
                                 return target_function.value(x, 0);
                               }),
                             target);
    data_out.add_data_vector(new_dof_handler, target, "u_bar");

    // In order to visualize the sum of the heat sources $\sum_k C^k
    // f_k(\mathbf x)$, we create a vector which will store the
    // interpolated values of this function.  Then, we loop through
    // the heat functions, create a vector to store the interpolated
    // data, call the VectorTools::interpolate() function to fill the
    // vector, multiply the interpolated data by the nonlocal dof
    // value $C^k$ (so that the heat plate is set to the correct
    // temperature), and then add this data to the sum of heat
    // sources. Because we can, we also add the vector for each source
    // individually to the DataOut object, so that they can be
    // visualized individually.
    Vector<double> full_heat_profile(new_dof_handler.n_dofs());

    for (unsigned int i = 0; i < heat_functions.size(); ++i)
      {
        Vector<double> hot_plate_i(new_dof_handler.n_dofs());

        VectorTools::interpolate(new_dof_handler,
                                 heat_functions[i],
                                 hot_plate_i);

        hot_plate_i *= solution[nonlocal_dofs[i]];
        full_heat_profile += hot_plate_i;

        const std::string data_name =
          "Heat_Source_" + Utilities::int_to_string(i);
        data_out.add_data_vector(new_dof_handler, hot_plate_i, data_name);
      }

    // Once all the heat functions have been combined, we add them to the
    // data_out object, and output everything into a file:
    data_out.add_data_vector(new_dof_handler,
                             full_heat_profile,
                             "Full_Heat_Profile");

    data_out.build_patches();

    std::ofstream output("solution.vtu");
    data_out.write_vtu(output);

    // Finally, we output the nonlocal coefficient values to the console:
    std::cout << "The c coefficients are " << std::endl;
    for (long unsigned int i = 0; i < nonlocal_dofs.size(); ++i)
      {
        std::cout << "\tc" << i + 1 << ": " << solution[nonlocal_dofs[i]]
                  << std::endl;
      }
  }


  // @sect4{Step93::run()}

  // The run() function runs through each step of the program, nothing new here:
  template <int dim>
  void Step93<dim>::run()
  {
    make_grid();
    setup_system();
    assemble_system();
    solve();
    output_results();
  }
} // namespace Step93


// @sect3{The main() function}

// The `main()` function looks essentially like that of most other tutorial
// programs.
int main()
{
  try
    {
      Step93::Step93<2> heat_optimization_problem;
      heat_optimization_problem.run();
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
