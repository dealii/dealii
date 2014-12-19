/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2006 - 2013 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * The deal.II library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE at
 * the top level of the deal.II distribution.
 *
 * ---------------------------------------------------------------------

 *
 * Author: Wolfgang Bangerth, Texas A&M University, 2006, 2007
 */


// @sect3{Include files}

// The first few files have already been covered in previous examples and will
// thus not be further commented on.
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/utilities.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/constraint_matrix.h>
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

// These are the new files we need. The first one provides an alternative to
// the usual SparsityPattern class and the CompressedSparsityPattern class
// already discussed in step-11 and step-18. The last two provide <i>hp</i>
// versions of the DoFHandler and FEValues classes as described in the
// introduction of this program.
#include <deal.II/lac/compressed_set_sparsity_pattern.h>
#include <deal.II/hp/dof_handler.h>
#include <deal.II/hp/fe_values.h>

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
    LaplaceProblem ();
    ~LaplaceProblem ();

    void run ();

  private:
    void setup_system ();
    void assemble_system ();
    void solve ();
    void create_coarse_grid ();
    void estimate_smoothness (Vector<float> &smoothness_indicators) const;
    void postprocess (const unsigned int cycle);

    Triangulation<dim>   triangulation;

    hp::DoFHandler<dim>      dof_handler;
    hp::FECollection<dim>    fe_collection;
    hp::QCollection<dim>     quadrature_collection;
    hp::QCollection<dim-1>   face_quadrature_collection;

    ConstraintMatrix     constraints;

    SparsityPattern      sparsity_pattern;
    SparseMatrix<double> system_matrix;

    Vector<double>       solution;
    Vector<double>       system_rhs;

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
    RightHandSide () : Function<dim> () {}

    virtual double value (const Point<dim>   &p,
                          const unsigned int  component) const;
  };


  template <int dim>
  double
  RightHandSide<dim>::value (const Point<dim>   &p,
                             const unsigned int  /*component*/) const
  {
    double product = 1;
    for (unsigned int d=0; d<dim; ++d)
      product *= (p[d]+1);
    return product;
  }




  // @sect3{Implementation of the main class}

  // @sect4{LaplaceProblem::LaplaceProblem}

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
  template <int dim>
  LaplaceProblem<dim>::LaplaceProblem ()
    :
    dof_handler (triangulation),
    max_degree (dim <= 2 ? 7 : 5)
  {
    for (unsigned int degree=2; degree<=max_degree; ++degree)
      {
        fe_collection.push_back (FE_Q<dim>(degree));
        quadrature_collection.push_back (QGauss<dim>(degree+1));
        face_quadrature_collection.push_back (QGauss<dim-1>(degree+1));
      }
  }


  // @sect4{LaplaceProblem::~LaplaceProblem}

  // The destructor is unchanged from what we already did in step-6:
  template <int dim>
  LaplaceProblem<dim>::~LaplaceProblem ()
  {
    dof_handler.clear ();
  }


  // @sect4{LaplaceProblem::setup_system}
  //
  // This function is again a verbatim copy of what we already did in
  // step-6. Despite function calls with exactly the same names and arguments,
  // the algorithms used internally are different in some aspect since the
  // dof_handler variable here is an hp object.
  template <int dim>
  void LaplaceProblem<dim>::setup_system ()
  {
    dof_handler.distribute_dofs (fe_collection);

    solution.reinit (dof_handler.n_dofs());
    system_rhs.reinit (dof_handler.n_dofs());

    constraints.clear ();
    DoFTools::make_hanging_node_constraints (dof_handler,
                                             constraints);
    VectorTools::interpolate_boundary_values (dof_handler,
                                              0,
                                              ZeroFunction<dim>(),
                                              constraints);
    constraints.close ();

    CompressedSetSparsityPattern csp (dof_handler.n_dofs(),
                                      dof_handler.n_dofs());
    DoFTools::make_sparsity_pattern (dof_handler, csp, constraints, false);
    sparsity_pattern.copy_from (csp);

    system_matrix.reinit (sparsity_pattern);
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
  void LaplaceProblem<dim>::assemble_system ()
  {
    hp::FEValues<dim> hp_fe_values (fe_collection,
                                    quadrature_collection,
                                    update_values    |  update_gradients |
                                    update_quadrature_points  |  update_JxW_values);

    const RightHandSide<dim> rhs_function;

    FullMatrix<double>   cell_matrix;
    Vector<double>       cell_rhs;

    std::vector<types::global_dof_index> local_dof_indices;

    typename hp::DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();
    for (; cell!=endc; ++cell)
      {
        const unsigned int   dofs_per_cell = cell->get_fe().dofs_per_cell;

        cell_matrix.reinit (dofs_per_cell, dofs_per_cell);
        cell_matrix = 0;

        cell_rhs.reinit (dofs_per_cell);
        cell_rhs = 0;

        hp_fe_values.reinit (cell);

        const FEValues<dim> &fe_values = hp_fe_values.get_present_fe_values ();

        std::vector<double>  rhs_values (fe_values.n_quadrature_points);
        rhs_function.value_list (fe_values.get_quadrature_points(),
                                 rhs_values);

        for (unsigned int q_point=0;
             q_point<fe_values.n_quadrature_points;
             ++q_point)
          for (unsigned int i=0; i<dofs_per_cell; ++i)
            {
              for (unsigned int j=0; j<dofs_per_cell; ++j)
                cell_matrix(i,j) += (fe_values.shape_grad(i,q_point) *
                                     fe_values.shape_grad(j,q_point) *
                                     fe_values.JxW(q_point));

              cell_rhs(i) += (fe_values.shape_value(i,q_point) *
                              rhs_values[q_point] *
                              fe_values.JxW(q_point));
            }

        local_dof_indices.resize (dofs_per_cell);
        cell->get_dof_indices (local_dof_indices);

        constraints.distribute_local_to_global (cell_matrix, cell_rhs,
                                                local_dof_indices,
                                                system_matrix, system_rhs);
      }
  }



  // @sect4{LaplaceProblem::solve}

  // The function solving the linear system is entirely unchanged from
  // previous examples. We simply try to reduce the initial residual (which
  // equals the $l_2$ norm of the right hand side) by a certain factor:
  template <int dim>
  void LaplaceProblem<dim>::solve ()
  {
    SolverControl           solver_control (system_rhs.size(),
                                            1e-8*system_rhs.l2_norm());
    SolverCG<>              cg (solver_control);

    PreconditionSSOR<> preconditioner;
    preconditioner.initialize(system_matrix, 1.2);

    cg.solve (system_matrix, solution, system_rhs,
              preconditioner);

    constraints.distribute (solution);
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
  void LaplaceProblem<dim>::postprocess (const unsigned int cycle)
  {
    // Let us start with computing estimated error and smoothness indicators,
    // which each are one number for each active cell of our
    // triangulation. For the error indicator, we use the KellyErrorEstimator
    // class as always. Estimating the smoothness is done in the respective
    // function of this class; that function is discussed further down below:
    Vector<float> estimated_error_per_cell (triangulation.n_active_cells());
    KellyErrorEstimator<dim>::estimate (dof_handler,
                                        face_quadrature_collection,
                                        typename FunctionMap<dim>::type(),
                                        solution,
                                        estimated_error_per_cell);


    Vector<float> smoothness_indicators (triangulation.n_active_cells());
    estimate_smoothness (smoothness_indicators);

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
      Vector<float> fe_degrees (triangulation.n_active_cells());
      {
        typename hp::DoFHandler<dim>::active_cell_iterator
        cell = dof_handler.begin_active(),
        endc = dof_handler.end();
        for (unsigned int index=0; cell!=endc; ++cell, ++index)
          fe_degrees(index)
            = fe_collection[cell->active_fe_index()].degree;
      }

      // With now all data vectors available -- solution, estimated errors and
      // smoothness indicators, and finite element degrees --, we create a
      // DataOut object for graphical output and attach all data. Note that
      // the DataOut class has a second template argument (which defaults to
      // DoFHandler@<dim@>, which is why we have never seen it in previous
      // tutorial programs) that indicates the type of DoF handler to be
      // used. Here, we have to use the hp::DoFHandler class:
      DataOut<dim,hp::DoFHandler<dim> > data_out;

      data_out.attach_dof_handler (dof_handler);
      data_out.add_data_vector (solution, "solution");
      data_out.add_data_vector (estimated_error_per_cell, "error");
      data_out.add_data_vector (smoothness_indicators, "smoothness");
      data_out.add_data_vector (fe_degrees, "fe_degree");
      data_out.build_patches ();

      // The final step in generating output is to determine a file name, open
      // the file, and write the data into it (here, we use VTK format):
      const std::string filename = "solution-" +
                                   Utilities::int_to_string (cycle, 2) +
                                   ".vtk";
      std::ofstream output (filename.c_str());
      data_out.write_vtk (output);
    }

    // After this, we would like to actually refine the mesh, in both $h$ and
    // $p$. The way we are going to do this is as follows: first, we use the
    // estimated error to flag those cells for refinement that have the
    // largest error. This is what we have always done:
    {
      GridRefinement::refine_and_coarsen_fixed_number (triangulation,
                                                       estimated_error_per_cell,
                                                       0.3, 0.03);

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
      float max_smoothness = *std::min_element (smoothness_indicators.begin(),
                                                smoothness_indicators.end()),
                             min_smoothness = *std::max_element (smoothness_indicators.begin(),
                                                                 smoothness_indicators.end());
      {
        typename hp::DoFHandler<dim>::active_cell_iterator
        cell = dof_handler.begin_active(),
        endc = dof_handler.end();
        for (unsigned int index=0; cell!=endc; ++cell, ++index)
          if (cell->refine_flag_set())
            {
              max_smoothness = std::max (max_smoothness,
                                         smoothness_indicators(index));
              min_smoothness = std::min (min_smoothness,
                                         smoothness_indicators(index));
            }
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
      {
        typename hp::DoFHandler<dim>::active_cell_iterator
        cell = dof_handler.begin_active(),
        endc = dof_handler.end();
        for (unsigned int index=0; cell!=endc; ++cell, ++index)
          if (cell->refine_flag_set()
              &&
              (smoothness_indicators(index) > threshold_smoothness)
              &&
              (cell->active_fe_index()+1 < fe_collection.size()))
            {
              cell->clear_refine_flag();
              cell->set_active_fe_index (cell->active_fe_index() + 1);
            }
      }

      // At the end of this procedure, we then refine the mesh. During this
      // process, children of cells undergoing bisection inherit their mother
      // cell's finite element index:
      triangulation.execute_coarsening_and_refinement ();
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
  void LaplaceProblem<2>::create_coarse_grid ()
  {
    const unsigned int dim = 2;

    static const Point<2> vertices_1[]
      = {  Point<2> (-1.,   -1.),
           Point<2> (-1./2, -1.),
           Point<2> (0.,    -1.),
           Point<2> (+1./2, -1.),
           Point<2> (+1,    -1.),

           Point<2> (-1.,   -1./2.),
           Point<2> (-1./2, -1./2.),
           Point<2> (0.,    -1./2.),
           Point<2> (+1./2, -1./2.),
           Point<2> (+1,    -1./2.),

           Point<2> (-1.,   0.),
           Point<2> (-1./2, 0.),
           Point<2> (+1./2, 0.),
           Point<2> (+1,    0.),

           Point<2> (-1.,   1./2.),
           Point<2> (-1./2, 1./2.),
           Point<2> (0.,    1./2.),
           Point<2> (+1./2, 1./2.),
           Point<2> (+1,    1./2.),

           Point<2> (-1.,   1.),
           Point<2> (-1./2, 1.),
           Point<2> (0.,    1.),
           Point<2> (+1./2, 1.),
           Point<2> (+1,    1.)
        };
    const unsigned int
    n_vertices = sizeof(vertices_1) / sizeof(vertices_1[0]);
    const std::vector<Point<dim> > vertices (&vertices_1[0],
                                             &vertices_1[n_vertices]);
    static const int cell_vertices[][GeometryInfo<dim>::vertices_per_cell]
    = {{0, 1, 5, 6},
      {1, 2, 6, 7},
      {2, 3, 7, 8},
      {3, 4, 8, 9},
      {5, 6, 10, 11},
      {8, 9, 12, 13},
      {10, 11, 14, 15},
      {12, 13, 17, 18},
      {14, 15, 19, 20},
      {15, 16, 20, 21},
      {16, 17, 21, 22},
      {17, 18, 22, 23}
    };
    const unsigned int
    n_cells = sizeof(cell_vertices) / sizeof(cell_vertices[0]);

    std::vector<CellData<dim> > cells (n_cells, CellData<dim>());
    for (unsigned int i=0; i<n_cells; ++i)
      {
        for (unsigned int j=0;
             j<GeometryInfo<dim>::vertices_per_cell;
             ++j)
          cells[i].vertices[j] = cell_vertices[i][j];
        cells[i].material_id = 0;
      }

    triangulation.create_triangulation (vertices,
                                        cells,
                                        SubCellData());
    triangulation.refine_global (3);
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
  void LaplaceProblem<dim>::run ()
  {
    for (unsigned int cycle=0; cycle<6; ++cycle)
      {
        std::cout << "Cycle " << cycle << ':' << std::endl;

        if (cycle == 0)
          create_coarse_grid ();

        setup_system ();

        std::cout << "   Number of active cells:       "
                  << triangulation.n_active_cells()
                  << std::endl
                  << "   Number of degrees of freedom: "
                  << dof_handler.n_dofs()
                  << std::endl
                  << "   Number of constraints       : "
                  << constraints.n_constraints()
                  << std::endl;

        assemble_system ();
        solve ();
        postprocess (cycle);
      }
  }


  // @sect4{LaplaceProblem::estimate_smoothness}

  // This last function of significance implements the algorithm to estimate
  // the smoothness exponent using the algorithms explained in detail in the
  // introduction. We will therefore only comment on those points that are of
  // implementational importance.
  template <int dim>
  void
  LaplaceProblem<dim>::
  estimate_smoothness (Vector<float> &smoothness_indicators) const
  {
    // The first thing we need to do is to define the Fourier vectors ${\bf
    // k}$ for which we want to compute Fourier coefficients of the solution
    // on each cell. In 2d, we pick those vectors ${\bf k}=(\pi i, \pi j)^T$
    // for which $\sqrt{i^2+j^2}\le N$, with $i,j$ integers and $N$ being the
    // maximal polynomial degree we use for the finite elements in this
    // program. The 3d case is handled analogously. 1d and dimensions higher
    // than 3 are not implemented, and we guard our implementation by making
    // sure that we receive an exception in case someone tries to compile the
    // program for any of these dimensions.
    //
    // We exclude ${\bf k}=0$ to avoid problems computing $|{\bf k}|^{-mu}$
    // and $\ln |{\bf k}|$. The other vectors are stored in the field
    // <code>k_vectors</code>. In addition, we store the square of the
    // magnitude of each of these vectors (up to a factor $\pi^2$) in the
    // <code>k_vectors_magnitude</code> array -- we will need that when we
    // attempt to find out which of those Fourier coefficients corresponding
    // to Fourier vectors of the same magnitude is the largest:
    const unsigned int N = max_degree;

    std::vector<Tensor<1,dim> > k_vectors;
    std::vector<unsigned int>   k_vectors_magnitude;
    switch (dim)
      {
      case 2:
      {
        for (unsigned int i=0; i<N; ++i)
          for (unsigned int j=0; j<N; ++j)
            if (!((i==0) && (j==0))
                &&
                (i*i + j*j < N*N))
              {
                k_vectors.push_back (Point<dim>(numbers::PI * i,
                                                numbers::PI * j));
                k_vectors_magnitude.push_back (i*i+j*j);
              }

        break;
      }

      case 3:
      {
        for (unsigned int i=0; i<N; ++i)
          for (unsigned int j=0; j<N; ++j)
            for (unsigned int k=0; k<N; ++k)
              if (!((i==0) && (j==0) && (k==0))
                  &&
                  (i*i + j*j + k*k < N*N))
                {
                  k_vectors.push_back (Point<dim>(numbers::PI * i,
                                                  numbers::PI * j,
                                                  numbers::PI * k));
                  k_vectors_magnitude.push_back (i*i+j*j+k*k);
                }

        break;
      }

      default:
        Assert (false, ExcNotImplemented());
      }

    // After we have set up the Fourier vectors, we also store their total
    // number for simplicity, and compute the logarithm of the magnitude of
    // each of these vectors since we will need it many times over further
    // down below:
    const unsigned n_fourier_modes = k_vectors.size();
    std::vector<double> ln_k (n_fourier_modes);
    for (unsigned int i=0; i<n_fourier_modes; ++i)
      ln_k[i] = std::log (k_vectors[i].norm());


    // Next, we need to assemble the matrices that do the Fourier transforms
    // for each of the finite elements we deal with, i.e. the matrices ${\cal
    // F}_{{\bf k},j}$ defined in the introduction. We have to do that for
    // each of the finite elements in use. Note that these matrices are
    // complex-valued, so we can't use the FullMatrix class. Instead, we use
    // the Table class template.
    std::vector<Table<2,std::complex<double> > >
    fourier_transform_matrices (fe_collection.size());

    // In order to compute them, we of course can't perform the Fourier
    // transform analytically, but have to approximate it using quadrature. To
    // this end, we use a quadrature formula that is obtained by iterating a
    // 2-point Gauss formula as many times as the maximal exponent we use for
    // the term $e^{i{\bf k}\cdot{\bf x}}$:
    QGauss<1>      base_quadrature (2);
    QIterated<dim> quadrature (base_quadrature, N);

    // With this, we then loop over all finite elements in use, reinitialize
    // the respective matrix ${\cal F}$ to the right size, and integrate each
    // entry of the matrix numerically as ${\cal F}_{{\bf k},j}=\sum_q
    // e^{i{\bf k}\cdot {\bf x}}\varphi_j({\bf x}_q) w_q$, where $x_q$ are the
    // quadrature points and $w_q$ are the quadrature weights. Note that the
    // imaginary unit $i=\sqrt{-1}$ is obtained from the standard C++ classes
    // using <code>std::complex@<double@>(0,1)</code>.

    // Because we work on the unit cell, we can do all this work without a
    // mapping from reference to real cell and consequently do not need the
    // FEValues class.
    for (unsigned int fe=0; fe<fe_collection.size(); ++fe)
      {
        fourier_transform_matrices[fe].reinit (n_fourier_modes,
                                               fe_collection[fe].dofs_per_cell);

        for (unsigned int k=0; k<n_fourier_modes; ++k)
          for (unsigned int j=0; j<fe_collection[fe].dofs_per_cell; ++j)
            {
              std::complex<double> sum = 0;
              for (unsigned int q=0; q<quadrature.size(); ++q)
                {
                  const Point<dim> x_q = quadrature.point(q);
                  sum += std::exp(std::complex<double>(0,1) *
                                  (k_vectors[k] * x_q)) *
                         fe_collection[fe].shape_value(j,x_q) *
                         quadrature.weight(q);
                }
              fourier_transform_matrices[fe](k,j)
                = sum / std::pow(2*numbers::PI, 1.*dim/2);
            }
      }

    // The next thing is to loop over all cells and do our work there, i.e. to
    // locally do the Fourier transform and estimate the decay coefficient. We
    // will use the following two arrays as scratch arrays in the loop and
    // allocate them here to avoid repeated memory allocations:
    std::vector<std::complex<double> > fourier_coefficients (n_fourier_modes);
    Vector<double>                     local_dof_values;

    // Then here is the loop:
    typename hp::DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();
    for (unsigned int index=0; cell!=endc; ++cell, ++index)
      {
        // Inside the loop, we first need to get the values of the local
        // degrees of freedom (which we put into the
        // <code>local_dof_values</code> array after setting it to the right
        // size) and then need to compute the Fourier transform by multiplying
        // this vector with the matrix ${\cal F}$ corresponding to this finite
        // element. We need to write out the multiplication by hand because
        // the objects holding the data do not have <code>vmult</code>-like
        // functions declared:
        local_dof_values.reinit (cell->get_fe().dofs_per_cell);
        cell->get_dof_values (solution, local_dof_values);

        for (unsigned int f=0; f<n_fourier_modes; ++f)
          {
            fourier_coefficients[f] = 0;

            for (unsigned int i=0; i<cell->get_fe().dofs_per_cell; ++i)
              fourier_coefficients[f] +=
                fourier_transform_matrices[cell->active_fe_index()](f,i)
                *
                local_dof_values(i);
          }

        // The next thing, as explained in the introduction, is that we wanted
        // to only fit our exponential decay of Fourier coefficients to the
        // largest coefficients for each possible value of $|{\bf k}|$. To
        // this end, we create a map that for each magnitude $|{\bf k}|$
        // stores the largest $|\hat U_{{\bf k}}|$ found so far, i.e. we
        // overwrite the existing value (or add it to the map) if no value for
        // the current $|{\bf k}|$ exists yet, or if the current value is
        // larger than the previously stored one:
        std::map<unsigned int, double> k_to_max_U_map;
        for (unsigned int f=0; f<n_fourier_modes; ++f)
          if ((k_to_max_U_map.find (k_vectors_magnitude[f]) ==
               k_to_max_U_map.end())
              ||
              (k_to_max_U_map[k_vectors_magnitude[f]] <
               std::abs (fourier_coefficients[f])))
            k_to_max_U_map[k_vectors_magnitude[f]]
              = std::abs (fourier_coefficients[f]);
        // Note that it comes in handy here that we have stored the magnitudes
        // of vectors as integers, since this way we do not have to deal with
        // round-off-sized differences between different values of $|{\bf
        // k}|$.

        // As the final task, we have to calculate the various contributions
        // to the formula for $\mu$. We'll only take those Fourier
        // coefficients with the largest magnitude for a given value of $|{\bf
        // k}|$ as explained above:
        double  sum_1           = 0,
                sum_ln_k        = 0,
                sum_ln_k_square = 0,
                sum_ln_U        = 0,
                sum_ln_U_ln_k   = 0;
        for (unsigned int f=0; f<n_fourier_modes; ++f)
          if (k_to_max_U_map[k_vectors_magnitude[f]] ==
              std::abs (fourier_coefficients[f]))
            {
              sum_1 += 1;
              sum_ln_k += ln_k[f];
              sum_ln_k_square += ln_k[f]*ln_k[f];
              sum_ln_U += std::log (std::abs (fourier_coefficients[f]));
              sum_ln_U_ln_k += std::log (std::abs (fourier_coefficients[f])) *
                               ln_k[f];
            }

        // With these so-computed sums, we can now evaluate the formula for
        // $\mu$ derived in the introduction:
        const double mu
          = (1./(sum_1*sum_ln_k_square - sum_ln_k*sum_ln_k)
             *
             (sum_ln_k*sum_ln_U - sum_1*sum_ln_U_ln_k));

        // The final step is to compute the Sobolev index $s=\mu-\frac d2$ and
        // store it in the vector of estimated values for each cell:
        smoothness_indicators(index) = mu - 1.*dim/2;
      }
  }
}


// @sect3{The main function}

// The main function is again verbatim what we had before: wrap creating and
// running an object of the main class into a <code>try</code> block and catch
// whatever exceptions are thrown, thereby producing meaningful output if
// anything should go wrong:
int main ()
{
  try
    {
      using namespace dealii;
      using namespace Step27;

      deallog.depth_console (0);

      LaplaceProblem<2> laplace_problem;
      laplace_problem.run ();
    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl << std::endl
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
      std::cerr << std::endl << std::endl
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
