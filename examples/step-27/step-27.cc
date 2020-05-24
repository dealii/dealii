/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2006 - 2020 by the deal.II authors
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

 *
 * Authors: Wolfgang Bangerth, Texas A&M University, 2006, 2007;
 *          Denis Davydov, University of Erlangen-Nuremberg, 2016.
 */


// @sect3{Include files}

// The first few files have already been covered in previous examples and will
// thus not be further commented on.
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
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
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/error_estimator.h>

// These are the new files we need. The first and second provide <i>hp</i>
// versions of the DoFHandler and FEValues classes as described in the
// introduction of this program. The last one provides Fourier transformation
// class on the unit cell.
#include <deal.II/hp/dof_handler.h>
#include <deal.II/hp/fe_values.h>
#include <deal.II/fe/fe_series.h>

// The last set of include files are standard C++ headers. We need support for
// complex numbers when we compute the Fourier transform.
#include <fstream>
#include <iostream>
#include <complex>


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
  // estimated smoothness of the solution). There is also a function that
  // computes this estimated smoothness, as discussed in the introduction.
  //
  // As far as member variables are concerned, we use the same structure as
  // already used in step-6, but instead of a regular DoFHandler we use an
  // object of type hp::DoFHandler, and we need collections instead of
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
    void estimate_smoothness(Vector<float> &smoothness_indicators);
    void postprocess(const unsigned int cycle);
    std::pair<bool, unsigned int> predicate(const TableIndices<dim> &indices);

    Triangulation<dim> triangulation;

    hp::DoFHandler<dim>      dof_handler;
    hp::FECollection<dim>    fe_collection;
    hp::QCollection<dim>     quadrature_collection;
    hp::QCollection<dim - 1> face_quadrature_collection;

    hp::QCollection<dim>                    fourier_q_collection;
    std::unique_ptr<FESeries::Fourier<dim>> fourier;
    std::vector<double>                     ln_k;
    Table<dim, std::complex<double>>        fourier_coefficients;

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
    virtual double value(const Point<dim> & p,
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
  // the hp::DoFHandler object with the triangulation, and then sets the
  // maximal polynomial degree to 7 (in 1d and 2d) or 5 (in 3d and higher). We
  // do so because using higher order polynomial degrees becomes prohibitively
  // expensive, especially in higher space dimensions.
  //
  // Following this, we fill the collections of finite element, and cell and
  // face quadrature objects. We start with quadratic elements, and each
  // quadrature formula is chosen so that it is appropriate for the matching
  // finite element in the hp::FECollection object.
  //
  // Finally, we initialize FESeries::Fourier object which will be used to
  // calculate coefficient in Fourier series as described in the introduction.
  // In addition to the hp::FECollection, we need to provide quadrature rules
  // hp::QCollection for integration on the reference cell.
  //
  // In order to resize fourier_coefficients Table, we use the following
  // auxiliary function
  template <int dim, typename T>
  void resize(Table<dim, T> &coeff, const unsigned int N)
  {
    TableIndices<dim> size;
    for (unsigned int d = 0; d < dim; d++)
      size[d] = N;
    coeff.reinit(size);
  }

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

    // As described in the introduction, we define the Fourier vectors ${\bf
    // k}$ for which we want to compute Fourier coefficients of the solution
    // on each cell as follows. In 2d, we will need coefficients corresponding
    // to vectors ${\bf k}=(2 \pi i, 2\pi j)^T$ for which $\sqrt{i^2+j^2}\le N$,
    // with $i,j$ integers and $N$ being the maximal polynomial degree we use
    // for the finite elements in this program. The FESeries::Fourier class'
    // constructor first parameter $N$ defines the number of coefficients in 1D
    // with the total number of coefficients being $N^{dim}$. Although we will
    // not use coefficients corresponding to
    // $\sqrt{i^2+j^2}> N$ and $i+j==0$, the overhead of their calculation is
    // minimal. The transformation matrices for each FiniteElement will be
    // calculated only once the first time they are required in the course of
    // hp-adaptive refinement. Because we work on the unit cell, we can do all
    // this work without a mapping from reference to real cell and consequently
    // can precalculate these matrices. The calculation of expansion
    // coefficients for a particular set of local degrees of freedom on a given
    // cell then follows as a simple matrix-vector product.
    // The 3d case is handled analogously.
    const unsigned int N = max_degree;

    // We will need to assemble the matrices that do the Fourier transforms
    // for each of the finite elements we deal with, i.e. the matrices ${\cal
    // F}_{{\bf k},j}$ defined in the introduction. We have to do that for
    // each of the finite elements in use. To that end we need a quadrature
    // rule. In this example we use the same quadrature formula for each
    // finite element, namely that is obtained by iterating a
    // 2-point Gauss formula as many times as the maximal exponent we use for
    // the term $e^{i{\bf k}\cdot{\bf x}}$:
    QGauss<1>      base_quadrature(2);
    QIterated<dim> quadrature(base_quadrature, N);
    for (unsigned int i = 0; i < fe_collection.size(); i++)
      fourier_q_collection.push_back(quadrature);

    // Now we are ready to set-up the FESeries::Fourier object
    const std::vector<unsigned int> n_coefficients_per_direction(
      fe_collection.size(), N);
    fourier =
      std::make_unique<FESeries::Fourier<dim>>(n_coefficients_per_direction,
                                               fe_collection,
                                               fourier_q_collection);

    // We need to resize the matrix of fourier coefficients according to the
    // number of modes N.
    resize(fourier_coefficients, N);
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
  // dof_handler variable here is an hp object.
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
  // resize them to the correct size (given by
  // <code>dofs_per_cell</code>). Because these classes are implement in such
  // a way that reducing the size of a matrix or vector does not release the
  // currently allocated memory (unless the new size is zero), the process of
  // resizing at the beginning of the loop will only require re-allocation of
  // memory during the first few iterations. Once we have found in a cell with
  // the maximal finite element degree, no more re-allocations will happen
  // because all subsequent <code>reinit</code> calls will only set the size
  // to something that fits the currently allocated memory. This is important
  // since allocating memory is expensive, and doing so every time we visit a
  // new cell would take significant compute time.
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
        const unsigned int dofs_per_cell = cell->get_fe().dofs_per_cell;

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
    // class as always. Estimating the smoothness is done in the respective
    // function of this class; that function is discussed further down below:
    Vector<float> estimated_error_per_cell(triangulation.n_active_cells());
    KellyErrorEstimator<dim>::estimate(
      dof_handler,
      face_quadrature_collection,
      std::map<types::boundary_id, const Function<dim> *>(),
      solution,
      estimated_error_per_cell);


    Vector<float> smoothness_indicators(triangulation.n_active_cells());
    estimate_smoothness(smoothness_indicators);

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
    // all integers, so that it what we use:
    {
      Vector<float> fe_degrees(triangulation.n_active_cells());
      for (const auto &cell : dof_handler.active_cell_iterators())
        fe_degrees(cell->active_cell_index()) =
          fe_collection[cell->active_fe_index()].degree;

      // With now all data vectors available -- solution, estimated errors and
      // smoothness indicators, and finite element degrees --, we create a
      // DataOut object for graphical output and attach all data. Note that
      // the DataOut class has a second template argument (which defaults to
      // DoFHandler@<dim@>, which is why we have never seen it in previous
      // tutorial programs) that indicates the type of DoF handler to be
      // used. Here, we have to use the hp::DoFHandler class:
      DataOut<dim, hp::DoFHandler<dim>> data_out;

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
      // threshold. For this, we first have to determine the maximal and
      // minimal values of the smoothness indicators of all flagged cells,
      // which we do using a loop over all cells and comparing current minimal
      // and maximal values. (We start with the minimal and maximal values of
      // <i>all</i> cells, a range within which the minimal and maximal values
      // on cells flagged for refinement must surely lie.) Absent any better
      // strategies, we will then set the threshold above which will increase
      // $p$ instead of reducing $h$ as the mean value between minimal and
      // maximal smoothness indicators on cells flagged for refinement:
      float max_smoothness = *std::min_element(smoothness_indicators.begin(),
                                               smoothness_indicators.end()),
            min_smoothness = *std::max_element(smoothness_indicators.begin(),
                                               smoothness_indicators.end());
      for (const auto &cell : dof_handler.active_cell_iterators())
        if (cell->refine_flag_set())
          {
            max_smoothness =
              std::max(max_smoothness,
                       smoothness_indicators(cell->active_cell_index()));
            min_smoothness =
              std::min(min_smoothness,
                       smoothness_indicators(cell->active_cell_index()));
          }
      const float threshold_smoothness = (max_smoothness + min_smoothness) / 2;

      // With this, we can go back, loop over all cells again, and for those
      // cells for which (i) the refinement flag is set, (ii) the smoothness
      // indicator is larger than the threshold, and (iii) we still have a
      // finite element with a polynomial degree higher than the current one
      // in the finite element collection, we then increase the polynomial
      // degree and in return remove the flag indicating that the cell should
      // undergo bisection. For all other cells, the refinement flags remain
      // untouched:
      for (const auto &cell : dof_handler.active_cell_iterators())
        if (cell->refine_flag_set() &&
            (smoothness_indicators(cell->active_cell_index()) >
             threshold_smoothness) &&
            (cell->active_fe_index() + 1 < fe_collection.size()))
          {
            cell->clear_refine_flag();
            cell->set_active_fe_index(cell->active_fe_index() + 1);
          }

      // At the end of this procedure, we then refine the mesh. During this
      // process, children of cells undergoing bisection inherit their mother
      // cell's finite element index:
      triangulation.execute_coarsening_and_refinement();
    }
  }


  // @sect4{LaplaceProblem::create_coarse_grid}

  // The following function is used when creating the initial grid. It is a
  // specialization for the 2d case, i.e. a corresponding function needs to be
  // implemented if the program is run in anything other then 2d. The function
  // is actually stolen from step-14 and generates the same mesh used already
  // there, i.e. the square domain with the square hole in the middle. The
  // meaning of the different parts of this function are explained in the
  // documentation of step-14:
  template <>
  void LaplaceProblem<2>::create_coarse_grid()
  {
    const unsigned int dim = 2;

    const std::vector<Point<2>> vertices = {
      {-1.0, -1.0}, {-0.5, -1.0}, {+0.0, -1.0}, {+0.5, -1.0}, {+1.0, -1.0}, //
      {-1.0, -0.5}, {-0.5, -0.5}, {+0.0, -0.5}, {+0.5, -0.5}, {+1.0, -0.5}, //
      {-1.0, +0.0}, {-0.5, +0.0}, {+0.5, +0.0}, {+1.0, +0.0},               //
      {-1.0, +0.5}, {-0.5, +0.5}, {+0.0, +0.5}, {+0.5, +0.5}, {+1.0, +0.5}, //
      {-1.0, +1.0}, {-0.5, +1.0}, {+0.0, +1.0}, {+0.5, +1.0}, {+1.0, +1.0}};

    const std::vector<std::array<int, GeometryInfo<dim>::vertices_per_cell>>
      cell_vertices = {{{0, 1, 5, 6}},
                       {{1, 2, 6, 7}},
                       {{2, 3, 7, 8}},
                       {{3, 4, 8, 9}},
                       {{5, 6, 10, 11}},
                       {{8, 9, 12, 13}},
                       {{10, 11, 14, 15}},
                       {{12, 13, 17, 18}},
                       {{14, 15, 19, 20}},
                       {{15, 16, 20, 21}},
                       {{16, 17, 21, 22}},
                       {{17, 18, 22, 23}}};

    const unsigned int n_cells = cell_vertices.size();

    std::vector<CellData<dim>> cells(n_cells, CellData<dim>());
    for (unsigned int i = 0; i < n_cells; ++i)
      {
        for (unsigned int j = 0; j < GeometryInfo<dim>::vertices_per_cell; ++j)
          cells[i].vertices[j] = cell_vertices[i][j];
        cells[i].material_id = 0;
      }

    triangulation.create_triangulation(vertices, cells, SubCellData());
    triangulation.refine_global(3);
  }



  // @sect4{LaplaceProblem::run}

  // This function implements the logic of the program, as did the respective
  // function in most of the previous programs already, see for example
  // step-6.
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


  // @sect4{LaplaceProblem::estimate_smoothness}

  // As described in the introduction, we will need to take the maximum
  // absolute value of fourier coefficients which correspond to $k$-vector
  // $|{\bf k}|= const$. To filter the coefficients Table we
  // will use the FESeries::process_coefficients() which requires a predicate
  // to be specified. The predicate should operate on TableIndices and return
  // a pair of <code>bool</code> and <code>unsigned int</code>. The latter
  // is the value of the map from TableIndicies to unsigned int.  It is
  // used to define subsets of coefficients from which we search for the one
  // with highest absolute value, i.e. $l^\infty$-norm. The <code>bool</code>
  // parameter defines which indices should be used in processing. In the
  // current case we are interested in coefficients which correspond to
  // $0 < i*i+j*j < N*N$ and $0 < i*i+j*j+k*k < N*N$ in 2D and 3D, respectively.
  template <int dim>
  std::pair<bool, unsigned int>
  LaplaceProblem<dim>::predicate(const TableIndices<dim> &ind)
  {
    unsigned int v = 0;
    for (unsigned int i = 0; i < dim; i++)
      v += ind[i] * ind[i];
    if (v > 0 && v < max_degree * max_degree)
      return std::make_pair(true, v);
    else
      return std::make_pair(false, v);
  }

  // This last function of significance implements the algorithm to estimate
  // the smoothness exponent using the algorithms explained in detail in the
  // introduction. We will therefore only comment on those points that are of
  // implementational importance.
  template <int dim>
  void
  LaplaceProblem<dim>::estimate_smoothness(Vector<float> &smoothness_indicators)
  {
    // Since most of the hard work is done for us in FESeries::Fourier and
    // we set up the object of this class in the constructor, what we are left
    // to do here is apply this class to calculate coefficients and then
    // perform linear regression to fit their decay slope.


    // First thing to do is to loop over all cells and do our work there, i.e.
    // to locally do the Fourier transform and estimate the decay coefficient.
    // We will use the following array as a scratch array in the loop to store
    // local DoF values:
    Vector<double> local_dof_values;

    // Then here is the loop:
    for (const auto &cell : dof_handler.active_cell_iterators())
      {
        // Inside the loop, we first need to get the values of the local
        // degrees of freedom (which we put into the
        // <code>local_dof_values</code> array after setting it to the right
        // size) and then need to compute the Fourier transform by multiplying
        // this vector with the matrix ${\cal F}$ corresponding to this finite
        // element. This is done by calling FESeries::Fourier::calculate(),
        // that has to be provided with the <code>local_dof_values</code>,
        // <code>cell->active_fe_index()</code> and a Table to store
        // coefficients.
        local_dof_values.reinit(cell->get_fe().dofs_per_cell);
        cell->get_dof_values(solution, local_dof_values);

        fourier->calculate(local_dof_values,
                           cell->active_fe_index(),
                           fourier_coefficients);

        // The next thing, as explained in the introduction, is that we wanted
        // to only fit our exponential decay of Fourier coefficients to the
        // largest coefficients for each possible value of $|{\bf k}|$. To
        // this end, we use FESeries::process_coefficients() to rework
        // coefficients into the desired format. We'll only take those Fourier
        // coefficients with the largest magnitude for a given value of $|{\bf
        // k}|$ and thereby need to use VectorTools::Linfty_norm:
        std::pair<std::vector<unsigned int>, std::vector<double>> res =
          FESeries::process_coefficients<dim>(
            fourier_coefficients,
            [this](const TableIndices<dim> &indices) {
              return this->predicate(indices);
            },
            VectorTools::Linfty_norm);

        Assert(res.first.size() == res.second.size(), ExcInternalError());

        // The first vector in the <code>std::pair</code> will store values of
        // the predicate, that is $i*i+j*j= const$ or $i*i+j*j+k*k = const$ in
        // 2D or 3D respectively. This vector will be the same for all the cells
        // so we can calculate logarithms of the corresponding Fourier vectors
        // $|{\bf k}|$ only once in the whole hp-refinement cycle:
        if (ln_k.size() == 0)
          {
            ln_k.resize(res.first.size(), 0);
            for (unsigned int f = 0; f < ln_k.size(); f++)
              ln_k[f] =
                std::log(2.0 * numbers::PI * std::sqrt(1. * res.first[f]));
          }

        // We have to calculate the logarithms of absolute values of
        // coefficients and use it in a linear regression fit to obtain $\mu$.
        for (double &residual_element : res.second)
          residual_element = std::log(residual_element);

        std::pair<double, double> fit =
          FESeries::linear_regression(ln_k, res.second);

        // The final step is to compute the Sobolev index $s=\mu-\frac d2$ and
        // store it in the vector of estimated values for each cell:
        smoothness_indicators(cell->active_cell_index()) =
          -fit.first - 1. * dim / 2;
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
