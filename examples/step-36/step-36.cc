/* ------------------------------------------------------------------------
 *
 * SPDX-License-Identifier: LGPL-2.1-or-later
 * Copyright (C) 2009 - 2024 by the deal.II authors
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
 * Authors: Toby D. Young, Polish Academy of Sciences,
 *          Wolfgang Bangerth, Texas A&M University
 */

// @sect3{Include files}

// As mentioned in the introduction, this program is essentially only a
// slightly revised version of step-4. As a consequence, most of the following
// include files are as used there, or at least as used already in previous
// tutorial programs:
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/function_parser.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/utilities.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/full_matrix.h>

// IndexSet is used to set the size of each PETScWrappers::MPI::Vector:
#include <deal.II/base/index_set.h>

// PETSc appears here because SLEPc depends on this library:
#include <deal.II/lac/petsc_sparse_matrix.h>
#include <deal.II/lac/petsc_vector.h>

// And then we need to actually import the interfaces for solvers that SLEPc
// provides:
#include <deal.II/lac/slepc_solver.h>

// We also need some standard C++:
#include <fstream>
#include <iostream>

// Finally, as in previous programs, we import all the deal.II class and
// function names into the namespace into which everything in this program
// will go:
namespace Step36
{
  using namespace dealii;

  // @sect3{The <code>EigenvalueProblem</code> class template}

  // Following is the class declaration for the main class template. It looks
  // pretty much exactly like what has already been shown in step-4:
  template <int dim>
  class EigenvalueProblem
  {
  public:
    EigenvalueProblem(const std::string &prm_file);
    void run();

  private:
    void         make_grid_and_dofs();
    void         assemble_system();
    unsigned int solve();
    void         output_results() const;

    Triangulation<dim> triangulation;
    const FE_Q<dim>    fe;
    DoFHandler<dim>    dof_handler;

    // With these exceptions: For our eigenvalue problem, we need both a
    // @ref GlossStiffnessMatrix "stiffness matrix" for the left hand side as well as a @ref GlossMassMatrix "mass matrix" for
    // the right hand side. We also need not just one solution function, but a
    // whole set of these for the eigenfunctions we want to compute, along
    // with the corresponding eigenvalues:
    PETScWrappers::SparseMatrix             stiffness_matrix, mass_matrix;
    std::vector<PETScWrappers::MPI::Vector> eigenfunctions;
    std::vector<double>                     eigenvalues;

    // And then we need an object that will store several run-time parameters
    // that we will specify in an input file:
    ParameterHandler parameters;

    // Finally, we will have an object that contains "constraints" on our
    // degrees of freedom. This could include hanging node constraints if we
    // had adaptively refined meshes (which we don't have in the current
    // program). Here, we will store the constraints for boundary nodes
    // $U_i=0$.
    AffineConstraints<double> constraints;
  };

  // @sect3{Implementation of the <code>EigenvalueProblem</code> class}

  // @sect4{EigenvalueProblem::EigenvalueProblem}

  // First up, the constructor. The main new part is handling the run-time
  // input parameters. We need to declare their existence first, and then read
  // their values from the input file whose name is specified as an argument
  // to this function:
  template <int dim>
  EigenvalueProblem<dim>::EigenvalueProblem(const std::string &prm_file)
    : fe(1)
    , dof_handler(triangulation)
  {
    // TODO investigate why the minimum number of refinement steps required to
    // obtain the correct eigenvalue degeneracies is 6
    parameters.declare_entry(
      "Global mesh refinement steps",
      "5",
      Patterns::Integer(0, 20),
      "The number of times the 1-cell coarse mesh should "
      "be refined globally for our computations.");
    parameters.declare_entry("Number of eigenvalues/eigenfunctions",
                             "5",
                             Patterns::Integer(0, 100),
                             "The number of eigenvalues/eigenfunctions "
                             "to be computed.");
    parameters.declare_entry("Potential",
                             "0",
                             Patterns::Anything(),
                             "A functional description of the potential.");

    parameters.parse_input(prm_file);
  }


  // @sect4{EigenvalueProblem::make_grid_and_dofs}

  // The next function creates a mesh on the domain $[-1,1]^d$, refines it as
  // many times as the input file calls for, and then attaches a DoFHandler to
  // it and initializes the matrices and vectors to their correct sizes. We
  // also build the constraints that correspond to the boundary values
  // $u|_{\partial\Omega}=0$.
  //
  // For the matrices, we use the PETSc wrappers. These have the ability to
  // allocate memory as necessary as non-zero entries are added. This seems
  // inefficient: we could as well first compute the sparsity pattern,
  // initialize the matrices with it, and as we then insert entries we can be
  // sure that we do not need to re-allocate memory and free the one used
  // previously. One way to do that would be to use code like this:
  // @code
  //   DynamicSparsityPattern
  //      dsp (dof_handler.n_dofs(),
  //           dof_handler.n_dofs());
  //   DoFTools::make_sparsity_pattern (dof_handler, dsp);
  //   dsp.compress ();
  //   stiffness_matrix.reinit (dsp);
  //   mass_matrix.reinit (dsp);
  // @endcode
  // instead of the two <code>reinit()</code> calls for the
  // stiffness and mass matrices below.
  //
  // This doesn't quite work, unfortunately. The code above may lead to a few
  // entries in the non-zero pattern to which we only ever write zero entries;
  // most notably, this holds true for off-diagonal entries for those rows and
  // columns that belong to boundary nodes. This shouldn't be a problem, but
  // for whatever reason, PETSc's ILU preconditioner, which we use to solve
  // linear systems in the eigenvalue solver, doesn't like these extra entries
  // and aborts with an error message.
  //
  // In the absence of any obvious way to avoid this, we simply settle for the
  // second best option, which is have PETSc allocate memory as
  // necessary. That said, since this is not a time critical part, this whole
  // affair is of no further importance.
  template <int dim>
  void EigenvalueProblem<dim>::make_grid_and_dofs()
  {
    GridGenerator::hyper_cube(triangulation, -1, 1);
    triangulation.refine_global(
      parameters.get_integer("Global mesh refinement steps"));
    dof_handler.distribute_dofs(fe);

    DoFTools::make_zero_boundary_constraints(dof_handler, constraints);
    constraints.close();

    stiffness_matrix.reinit(dof_handler.n_dofs(),
                            dof_handler.n_dofs(),
                            dof_handler.max_couplings_between_dofs());
    mass_matrix.reinit(dof_handler.n_dofs(),
                       dof_handler.n_dofs(),
                       dof_handler.max_couplings_between_dofs());

    // The next step is to take care of the eigenspectrum. In this case, the
    // outputs are eigenvalues and eigenfunctions, so we set the size of the
    // list of eigenfunctions and eigenvalues to be as large as we asked for
    // in the input file. When using a PETScWrappers::MPI::Vector, the Vector
    // is initialized using an IndexSet. IndexSet is used not only to resize the
    // PETScWrappers::MPI::Vector but it also associates an index in the
    // PETScWrappers::MPI::Vector with a degree of freedom (see step-40 for a
    // more detailed explanation). The function complete_index_set() creates
    // an IndexSet where every valid index is part of the set. Note that this
    // program can only be run sequentially and will throw an exception if used
    // in parallel.
    IndexSet eigenfunction_index_set = dof_handler.locally_owned_dofs();
    eigenfunctions.resize(
      parameters.get_integer("Number of eigenvalues/eigenfunctions"));
    for (auto &eigenfunction : eigenfunctions)
      eigenfunction.reinit(eigenfunction_index_set, MPI_COMM_WORLD);

    eigenvalues.resize(eigenfunctions.size());
  }


  // @sect4{EigenvalueProblem::assemble_system}

  // Here, we assemble the global stiffness and mass matrices from local
  // contributions $A^K_{ij} = \int_K \nabla\varphi_i(\mathbf x) \cdot
  // \nabla\varphi_j(\mathbf x) + V(\mathbf x)\varphi_i(\mathbf
  // x)\varphi_j(\mathbf x)$ and $M^K_{ij} = \int_K \varphi_i(\mathbf
  // x)\varphi_j(\mathbf x)$ respectively. This function should be immediately
  // familiar if you've seen previous tutorial programs. The only thing new
  // would be setting up an object that described the potential $V(\mathbf x)$
  // using the expression that we got from the input file. We then need to
  // evaluate this object at the quadrature points on each cell. If you've
  // seen how to evaluate function objects (see, for example the coefficient
  // in step-5), the code here will also look rather familiar.
  template <int dim>
  void EigenvalueProblem<dim>::assemble_system()
  {
    const QGauss<dim> quadrature_formula(fe.degree + 1);

    FEValues<dim> fe_values(fe,
                            quadrature_formula,
                            update_values | update_gradients |
                              update_quadrature_points | update_JxW_values);

    const unsigned int dofs_per_cell = fe.n_dofs_per_cell();
    const unsigned int n_q_points    = quadrature_formula.size();

    FullMatrix<double> cell_stiffness_matrix(dofs_per_cell, dofs_per_cell);
    FullMatrix<double> cell_mass_matrix(dofs_per_cell, dofs_per_cell);

    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

    FunctionParser<dim> potential;
    potential.initialize(FunctionParser<dim>::default_variable_names(),
                         parameters.get("Potential"),
                         typename FunctionParser<dim>::ConstMap());

    std::vector<double> potential_values(n_q_points);
    for (const auto &cell : dof_handler.active_cell_iterators())
      {
        fe_values.reinit(cell);
        cell_stiffness_matrix = 0;
        cell_mass_matrix      = 0;

        potential.value_list(fe_values.get_quadrature_points(),
                             potential_values);

        for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
          for (unsigned int i = 0; i < dofs_per_cell; ++i)
            for (unsigned int j = 0; j < dofs_per_cell; ++j)
              {
                cell_stiffness_matrix(i, j) +=           //
                  (fe_values.shape_grad(i, q_point) *    //
                     fe_values.shape_grad(j, q_point)    //
                   +                                     //
                   potential_values[q_point] *           //
                     fe_values.shape_value(i, q_point) * //
                     fe_values.shape_value(j, q_point)   //
                   ) *                                   //
                  fe_values.JxW(q_point);                //

                cell_mass_matrix(i, j) +=              //
                  (fe_values.shape_value(i, q_point) * //
                   fe_values.shape_value(j, q_point)   //
                   ) *                                 //
                  fe_values.JxW(q_point);              //
              }

        // Now that we have the local matrix contributions, we transfer them
        // into the global objects and take care of zero boundary constraints:
        cell->get_dof_indices(local_dof_indices);

        constraints.distribute_local_to_global(cell_stiffness_matrix,
                                               local_dof_indices,
                                               stiffness_matrix);
        constraints.distribute_local_to_global(cell_mass_matrix,
                                               local_dof_indices,
                                               mass_matrix);
      }

    // At the end of the function, we tell PETSc that the matrices have now
    // been fully assembled and that the sparse matrix representation can now
    // be compressed as no more entries will be added:
    stiffness_matrix.compress(VectorOperation::add);
    mass_matrix.compress(VectorOperation::add);


    // Before leaving the function, we calculate spurious eigenvalues,
    // introduced to the system by zero Dirichlet constraints. As
    // discussed in the introduction, the use of Dirichlet boundary
    // conditions coupled with the fact that the degrees of freedom
    // located at the boundary of the domain remain part of the linear
    // system we solve, introduces a number of spurious eigenvalues.
    // Below, we output the interval within which they all lie to
    // ensure that we can ignore them should they show up in our
    // computations.
    double min_spurious_eigenvalue = std::numeric_limits<double>::max(),
           max_spurious_eigenvalue = std::numeric_limits<double>::lowest();

    for (unsigned int i = 0; i < dof_handler.n_dofs(); ++i)
      if (constraints.is_constrained(i))
        {
          const double ev         = stiffness_matrix(i, i) / mass_matrix(i, i);
          min_spurious_eigenvalue = std::min(min_spurious_eigenvalue, ev);
          max_spurious_eigenvalue = std::max(max_spurious_eigenvalue, ev);
        }

    std::cout << "   Spurious eigenvalues are all in the interval " << '['
              << min_spurious_eigenvalue << ',' << max_spurious_eigenvalue
              << ']' << std::endl;
  }


  // @sect4{EigenvalueProblem::solve}

  // This is the key new functionality of the program. Now that the system is
  // set up, here is a good time to actually solve the problem: As with other
  // examples this is done using a "solve" routine. Essentially, it works as
  // in other programs: you set up a SolverControl object that describes the
  // accuracy to which we want to solve the linear systems, and then we select
  // the kind of solver we want. Here we choose the Krylov-Schur solver of
  // SLEPc, a pretty fast and robust choice for this kind of problem:
  template <int dim>
  unsigned int EigenvalueProblem<dim>::solve()
  {
    // We start here, as we normally do, by assigning convergence control we
    // want:
    SolverControl                    solver_control(dof_handler.n_dofs(), 1e-9);
    SLEPcWrappers::SolverKrylovSchur eigensolver(solver_control);

    // Before we actually solve for the eigenfunctions and -values, we have to
    // also select which set of eigenvalues to solve for. Lets select those
    // eigenvalues and corresponding eigenfunctions with the smallest real
    // part (in fact, the problem we solve here is symmetric and so the
    // eigenvalues are purely real). After that, we can actually let SLEPc do
    // its work:
    eigensolver.set_which_eigenpairs(EPS_SMALLEST_REAL);

    eigensolver.set_problem_type(EPS_GHEP);

    eigensolver.solve(stiffness_matrix,
                      mass_matrix,
                      eigenvalues,
                      eigenfunctions,
                      eigenfunctions.size());

    // The output of the call above is a set of vectors and values. In
    // eigenvalue problems, the eigenfunctions are only determined up to a
    // constant that can be fixed pretty arbitrarily. Knowing nothing about
    // the origin of the eigenvalue problem, SLEPc has no other choice than to
    // normalize the eigenvectors to one in the $l_2$ (vector)
    // norm. Unfortunately this norm has little to do with any norm we may be
    // interested from a eigenfunction perspective: the $L_2(\Omega)$ norm, or
    // maybe the $L_\infty(\Omega)$ norm.
    //
    // Let us choose the latter and rescale eigenfunctions so that they have
    // $\|\phi_i(\mathbf x)\|_{L^\infty(\Omega)}=1$ instead of
    // $\|\Phi\|_{l_2}=1$ (where $\phi_i$ is the $i$th eigen<i>function</i>
    // and $\Phi_i$ the corresponding vector of nodal values). For the $Q_1$
    // elements chosen here, we know that the maximum of the function
    // $\phi_i(\mathbf x)$ is attained at one of the nodes, so $\max_{\mathbf
    // x}\phi_i(\mathbf x)=\max_j (\Phi_i)_j$, making the normalization in the
    // $L_\infty$ norm trivial. Note that this doesn't work as easily if we
    // had chosen $Q_k$ elements with $k>1$: there, the maximum of a function
    // does not necessarily have to be attained at a node, and so
    // $\max_{\mathbf x}\phi_i(\mathbf x)\ge\max_j (\Phi_i)_j$ (although the
    // equality is usually nearly true).
    for (auto &eigenfunction : eigenfunctions)
      eigenfunction /= eigenfunction.linfty_norm();

    // Finally return the number of iterations it took to converge:
    return solver_control.last_step();
  }


  // @sect4{EigenvalueProblem::output_results}

  // This is the last significant function of this program. It uses the
  // DataOut class to generate graphical output from the eigenfunctions for
  // later visualization. It works as in many of the other tutorial programs.
  //
  // The whole collection of functions is then output as a single VTK file.
  template <int dim>
  void EigenvalueProblem<dim>::output_results() const
  {
    DataOut<dim> data_out;

    data_out.attach_dof_handler(dof_handler);

    for (unsigned int i = 0; i < eigenfunctions.size(); ++i)
      data_out.add_data_vector(eigenfunctions[i],
                               "eigenfunction_" + Utilities::int_to_string(i));

    // The only thing worth discussing may be that because the potential is
    // specified as a function expression in the input file, it would be nice
    // to also have it as a graphical representation along with the
    // eigenfunctions. The process to achieve this is relatively
    // straightforward: we build an object that represents $V(\mathbf x)$ and
    // then we interpolate this continuous function onto the finite element
    // space. The result we also attach to the DataOut object for
    // visualization.
    Vector<double> projected_potential(dof_handler.n_dofs());
    {
      FunctionParser<dim> potential;
      potential.initialize(FunctionParser<dim>::default_variable_names(),
                           parameters.get("Potential"),
                           typename FunctionParser<dim>::ConstMap());
      VectorTools::interpolate(dof_handler, potential, projected_potential);
    }
    data_out.add_data_vector(projected_potential, "interpolated_potential");

    data_out.build_patches();

    std::ofstream output("eigenvectors.vtk");
    data_out.write_vtk(output);
  }


  // @sect4{EigenvalueProblem::run}

  // This is the function which has the top-level control over everything. It
  // is almost exactly the same as in step-4:
  template <int dim>
  void EigenvalueProblem<dim>::run()
  {
    make_grid_and_dofs();

    std::cout << "   Number of active cells:       "
              << triangulation.n_active_cells() << std::endl
              << "   Number of degrees of freedom: " << dof_handler.n_dofs()
              << std::endl;

    assemble_system();

    const unsigned int n_iterations = solve();
    std::cout << "   Solver converged in " << n_iterations << " iterations."
              << std::endl;

    output_results();

    std::cout << std::endl;
    for (unsigned int i = 0; i < eigenvalues.size(); ++i)
      std::cout << "      Eigenvalue " << i << " : " << eigenvalues[i]
                << std::endl;
  }
} // namespace Step36

// @sect3{The <code>main</code> function}
int main(int argc, char **argv)
{
  try
    {
      using namespace dealii;
      using namespace Step36;

      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);


      // This program can only be run in serial. Otherwise, throw an exception.
      AssertThrow(Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD) == 1,
                  ExcMessage(
                    "This program can only be run in serial, use ./step-36"));

      EigenvalueProblem<2> problem("step-36.prm");
      problem.run();
    }

  // All the while, we are watching out if any exceptions should have been
  // generated. If that is so, we panic...
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

  // If no exceptions are thrown, then we tell the program to stop monkeying
  // around and exit nicely:
  std::cout << std::endl << "   Job done." << std::endl;

  return 0;
}
