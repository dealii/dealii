/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2001 - 2013 by the deal.II authors
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
 * Author: Wolfgang Bangerth, University of Heidelberg, 2001
 */


// As usual, the program starts with a rather long list of include files which
// you are probably already used to by now:
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/table_handler.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>

// Just this one is new: it declares a class
// <code>CompressedSparsityPattern</code>, which we will use and explain
// further down below.
#include <deal.II/lac/compressed_sparsity_pattern.h>

// We will make use of the std::find algorithm of the C++ standard library, so
// we have to include the following file for its declaration:
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <cmath>

// The last step is as in all previous programs:
namespace Step11
{
  using namespace dealii;

  // Then we declare a class which represents the solution of a Laplace
  // problem. As this example program is based on step-5, the class looks
  // rather the same, with the sole structural difference that the functions
  // <code>assemble_system</code> now calls <code>solve</code> itself, and is
  // thus called <code>assemble_and_solve</code>, and that the output function
  // was dropped since the solution function is so boring that it is not worth
  // being viewed.
  //
  // The only other noteworthy change is that the constructor takes a value
  // representing the polynomial degree of the mapping to be used later on,
  // and that it has another member variable representing exactly this
  // mapping. In general, this variable will occur in real applications at the
  // same places where the finite element is declared or used.
  template <int dim>
  class LaplaceProblem
  {
  public:
    LaplaceProblem (const unsigned int mapping_degree);
    void run ();

  private:
    void setup_system ();
    void assemble_and_solve ();
    void solve ();

    Triangulation<dim>   triangulation;
    FE_Q<dim>            fe;
    DoFHandler<dim>      dof_handler;
    MappingQ<dim>        mapping;

    SparsityPattern      sparsity_pattern;
    SparseMatrix<double> system_matrix;
    ConstraintMatrix     mean_value_constraints;

    Vector<double>       solution;
    Vector<double>       system_rhs;

    TableHandler         output_table;
  };



  // Construct such an object, by initializing the variables. Here, we use
  // linear finite elements (the argument to the <code>fe</code> variable
  // denotes the polynomial degree), and mappings of given order. Print to
  // screen what we are about to do.
  template <int dim>
  LaplaceProblem<dim>::LaplaceProblem (const unsigned int mapping_degree) :
    fe (1),
    dof_handler (triangulation),
    mapping (mapping_degree)
  {
    std::cout << "Using mapping with degree " << mapping_degree << ":"
              << std::endl
              << "============================"
              << std::endl;
  }



  // The first task is to set up the variables for this problem. This includes
  // generating a valid <code>DoFHandler</code> object, as well as the
  // sparsity patterns for the matrix, and the object representing the
  // constraints that the mean value of the degrees of freedom on the boundary
  // be zero.
  template <int dim>
  void LaplaceProblem<dim>::setup_system ()
  {
    // The first task is trivial: generate an enumeration of the degrees of
    // freedom, and initialize solution and right hand side vector to their
    // correct sizes:
    dof_handler.distribute_dofs (fe);
    solution.reinit (dof_handler.n_dofs());
    system_rhs.reinit (dof_handler.n_dofs());

    // Next task is to construct the object representing the constraint that
    // the mean value of the degrees of freedom on the boundary shall be
    // zero. For this, we first want a list of those nodes which are actually
    // at the boundary. The <code>DoFTools</code> class has a function that
    // returns an array of Boolean values where <code>true</code> indicates
    // that the node is at the boundary. The second argument denotes a mask
    // selecting which components of vector valued finite elements we want to
    // be considered. This sort of information is encoded using the
    // ComponentMask class (see also @ref GlossComponentMask). Since we have a
    // scalar finite element anyway, this mask in reality should have only one
    // entry with a <code>true</code> value. However, the ComponentMask class
    // has semantics that allow it to represents a mask of indefinite size
    // whose every element equals <code>true</code> when one just default
    // constructs such an object, so this is what we'll do here.
    std::vector<bool> boundary_dofs (dof_handler.n_dofs(), false);
    DoFTools::extract_boundary_dofs (dof_handler,
                                     ComponentMask(),
                                     boundary_dofs);

    // Now first for the generation of the constraints: as mentioned in the
    // introduction, we constrain one of the nodes on the boundary by the
    // values of all other DoFs on the boundary. So, let us first pick out the
    // first boundary node from this list. We do that by searching for the
    // first <code>true</code> value in the array (note that
    // <code>std::find</code> returns an iterator to this element), and
    // computing its distance to the overall first element in the array to get
    // its index:
    const unsigned int first_boundary_dof
      = std::distance (boundary_dofs.begin(),
                       std::find (boundary_dofs.begin(),
                                  boundary_dofs.end(),
                                  true));

    // Then generate a constraints object with just this one constraint. First
    // clear all previous content (which might reside there from the previous
    // computation on a once coarser grid), then add this one line
    // constraining the <code>first_boundary_dof</code> to the sum of other
    // boundary DoFs each with weight -1. Finally, close the constraints
    // object, i.e. do some internal bookkeeping on it for faster processing
    // of what is to come later:
    mean_value_constraints.clear ();
    mean_value_constraints.add_line (first_boundary_dof);
    for (unsigned int i=first_boundary_dof+1; i<dof_handler.n_dofs(); ++i)
      if (boundary_dofs[i] == true)
        mean_value_constraints.add_entry (first_boundary_dof,
                                          i, -1);
    mean_value_constraints.close ();

    // Next task is to generate a sparsity pattern. This is indeed a tricky
    // task here. Usually, we just call
    // <code>DoFTools::make_sparsity_pattern</code> and condense the result
    // using the hanging node constraints. We have no hanging node constraints
    // here (since we only refine globally in this example), but we have this
    // global constraint on the boundary. This poses one severe problem in
    // this context: the <code>SparsityPattern</code> class wants us to state
    // beforehand the maximal number of entries per row, either for all rows
    // or for each row separately. There are functions in the library which
    // can tell you this number in case you just have hanging node constraints
    // (namely <code>DoFHandler::max_coupling_between_dofs</code>), but how is
    // this for the present case? The difficulty arises because the
    // elimination of the constrained degree of freedom requires a number of
    // additional entries in the matrix at places that are not so simple to
    // determine. We would therefore have a problem had we to give a maximal
    // number of entries per row here.
    //
    // Since this can be so difficult that no reasonable answer can be given
    // that allows allocation of only a reasonable amount of memory, there is
    // a class <code>CompressedSparsityPattern</code>, that can help us out
    // here. It does not require that we know in advance how many entries rows
    // could have, but allows just about any length. It is thus significantly
    // more flexible in case you do not have good estimates of row lengths,
    // however at the price that building up such a pattern is also
    // significantly more expensive than building up a pattern for which you
    // had information in advance. Nevertheless, as we have no other choice
    // here, we'll just build such an object by initializing it with the
    // dimensions of the matrix and calling another function
    // <code>DoFTools::make_sparsity_pattern</code> to get the sparsity
    // pattern due to the differential operator, then condense it with the
    // constraints object which adds those positions in the sparsity pattern
    // that are required for the elimination of the constraint.
    CompressedSparsityPattern csp (dof_handler.n_dofs(),
                                   dof_handler.n_dofs());
    DoFTools::make_sparsity_pattern (dof_handler, csp);
    mean_value_constraints.condense (csp);

    // Finally, once we have the full pattern, we can initialize an object of
    // type <code>SparsityPattern</code> from it and in turn initialize the
    // matrix with it. Note that this is actually necessary, since the
    // <code>CompressedSparsityPattern</code> is so inefficient compared to
    // the <code>SparsityPattern</code> class due to the more flexible data
    // structures it has to use, that we can impossibly base the sparse matrix
    // class on it, but rather need an object of type
    // <code>SparsityPattern</code>, which we generate by copying from the
    // intermediate object.
    //
    // As a further sidenote, you will notice that we do not explicitly have
    // to <code>compress</code> the sparsity pattern here. This, of course, is
    // due to the fact that the <code>copy_from</code> function generates a
    // compressed object right from the start, to which you cannot add new
    // entries anymore. The <code>compress</code> call is therefore implicit
    // in the <code>copy_from</code> call.
    sparsity_pattern.copy_from (csp);
    system_matrix.reinit (sparsity_pattern);
  }



  // The next function then assembles the linear system of equations, solves
  // it, and evaluates the solution. This then makes three actions, and we
  // will put them into eight true statements (excluding declaration of
  // variables, and handling of temporary vectors). Thus, this function is
  // something for the very lazy. Nevertheless, the functions called are
  // rather powerful, and through them this function uses a good deal of the
  // whole library. But let's look at each of the steps.
  template <int dim>
  void LaplaceProblem<dim>::assemble_and_solve ()
  {

    // First, we have to assemble the matrix and the right hand side. In all
    // previous examples, we have investigated various ways how to do this
    // manually. However, since the Laplace matrix and simple right hand sides
    // appear so frequently in applications, the library provides functions
    // for actually doing this for you, i.e. they perform the loop over all
    // cells, setting up the local matrices and vectors, and putting them
    // together for the end result.
    //
    // The following are the two most commonly used ones: creation of the
    // Laplace matrix and creation of a right hand side vector from body or
    // boundary forces. They take the mapping object, the
    // <code>DoFHandler</code> object representing the degrees of freedom and
    // the finite element in use, a quadrature formula to be used, and the
    // output object. The function that creates a right hand side vector also
    // has to take a function object describing the (continuous) right hand
    // side function.
    //
    // Let us look at the way the matrix and body forces are integrated:
    const unsigned int gauss_degree
      = std::max (static_cast<unsigned int>(std::ceil(1.*(mapping.get_degree()+1)/2)),
                  2U);
    MatrixTools::create_laplace_matrix (mapping, dof_handler,
                                        QGauss<dim>(gauss_degree),
                                        system_matrix);
    VectorTools::create_right_hand_side (mapping, dof_handler,
                                         QGauss<dim>(gauss_degree),
                                         ConstantFunction<dim>(-2),
                                         system_rhs);
    // That's quite simple, right?
    //
    // Two remarks are in order, though: First, these functions are used in a
    // lot of contexts. Maybe you want to create a Laplace or mass matrix for
    // a vector values finite element; or you want to use the default Q1
    // mapping; or you want to assembled the matrix with a coefficient in the
    // Laplace operator. For this reason, there are quite a large number of
    // variants of these functions in the <code>MatrixCreator</code> and
    // <code>MatrixTools</code> classes. Whenever you need a slightly
    // different version of these functions than the ones called above, it is
    // certainly worthwhile to take a look at the documentation and to check
    // whether something fits your needs.
    //
    // The second remark concerns the quadrature formula we use: we want to
    // integrate over bilinear shape functions, so we know that we have to use
    // at least a Gauss2 quadrature formula. On the other hand, we want to
    // have the quadrature rule to have at least the order of the boundary
    // approximation. Since the order of Gauss-r is 2r, and the order of the
    // boundary approximation using polynomials of degree p is p+1, we know
    // that 2r@>=p+1. Since r has to be an integer and (as mentioned above)
    // has to be at least 2, this makes up for the formula above computing
    // <code>gauss_degree</code>.
    //
    // Since the generation of the body force contributions to the right hand
    // side vector was so simple, we do that all over again for the boundary
    // forces as well: allocate a vector of the right size and call the right
    // function. The boundary function has constant values, so we can generate
    // an object from the library on the fly, and we use the same quadrature
    // formula as above, but this time of lower dimension since we integrate
    // over faces now instead of cells:
    Vector<double> tmp (system_rhs.size());
    VectorTools::create_boundary_right_hand_side (mapping, dof_handler,
                                                  QGauss<dim-1>(gauss_degree),
                                                  ConstantFunction<dim>(1),
                                                  tmp);
    // Then add the contributions from the boundary to those from the interior
    // of the domain:
    system_rhs += tmp;
    // For assembling the right hand side, we had to use two different vector
    // objects, and later add them together. The reason we had to do so is
    // that the <code>VectorTools::create_right_hand_side</code> and
    // <code>VectorTools::create_boundary_right_hand_side</code> functions
    // first clear the output vector, rather than adding up their results to
    // previous contents. This can reasonably be called a design flaw in the
    // library made in its infancy, but unfortunately things are as they are
    // for some time now and it is difficult to change such things that
    // silently break existing code, so we have to live with that.

    // Now, the linear system is set up, so we can eliminate the one degree of
    // freedom which we constrained to the other DoFs on the boundary for the
    // mean value constraint from matrix and right hand side vector, and solve
    // the system. After that, distribute the constraints again, which in this
    // case means setting the constrained degree of freedom to its proper
    // value
    mean_value_constraints.condense (system_matrix);
    mean_value_constraints.condense (system_rhs);

    solve ();
    mean_value_constraints.distribute (solution);

    // Finally, evaluate what we got as solution. As stated in the
    // introduction, we are interested in the H1 semi-norm of the
    // solution. Here, as well, we have a function in the library that does
    // this, although in a slightly non-obvious way: the
    // <code>VectorTools::integrate_difference</code> function integrates the
    // norm of the difference between a finite element function and a
    // continuous function. If we therefore want the norm of a finite element
    // field, we just put the continuous function to zero. Note that this
    // function, just as so many other ones in the library as well, has at
    // least two versions, one which takes a mapping as argument (which we
    // make us of here), and the one which we have used in previous examples
    // which implicitly uses <code>MappingQ1</code>.  Also note that we take a
    // quadrature formula of one degree higher, in order to avoid
    // superconvergence effects where the solution happens to be especially
    // close to the exact solution at certain points (we don't know whether
    // this might be the case here, but there are cases known of this, and we
    // just want to make sure):
    Vector<float> norm_per_cell (triangulation.n_active_cells());
    VectorTools::integrate_difference (mapping, dof_handler,
                                       solution,
                                       ZeroFunction<dim>(),
                                       norm_per_cell,
                                       QGauss<dim>(gauss_degree+1),
                                       VectorTools::H1_seminorm);
    // Then, the function just called returns its results as a vector of
    // values each of which denotes the norm on one cell. To get the global
    // norm, a simple computation shows that we have to take the l2 norm of
    // the vector:
    const double norm = norm_per_cell.l2_norm();

    // Last task -- generate output:
    output_table.add_value ("cells", triangulation.n_active_cells());
    output_table.add_value ("|u|_1", norm);
    output_table.add_value ("error", std::fabs(norm-std::sqrt(3.14159265358/2)));
  }



  // The following function solving the linear system of equations is copied
  // from step-5 and is explained there in some detail:
  template <int dim>
  void LaplaceProblem<dim>::solve ()
  {
    SolverControl           solver_control (1000, 1e-12);
    SolverCG<>              cg (solver_control);

    PreconditionSSOR<> preconditioner;
    preconditioner.initialize(system_matrix, 1.2);

    cg.solve (system_matrix, solution, system_rhs,
              preconditioner);
  }



  // Finally the main function controlling the different steps to be
  // performed. Its content is rather straightforward, generating a
  // triangulation of a circle, associating a boundary to it, and then doing
  // several cycles on subsequently finer grids. Note again that we have put
  // mesh refinement into the loop header; this may be something for a test
  // program, but for real applications you should consider that this implies
  // that the mesh is refined after the loop is executed the last time since
  // the increment clause (the last part of the three-parted loop header) is
  // executed before the comparison part (the second one), which may be rather
  // costly if the mesh is already quite refined. In that case, you should
  // arrange code such that the mesh is not further refined after the last
  // loop run (or you should do it at the beginning of each run except for the
  // first one).
  template <int dim>
  void LaplaceProblem<dim>::run ()
  {
    GridGenerator::hyper_ball (triangulation);
    static const HyperBallBoundary<dim> boundary;
    triangulation.set_boundary (0, boundary);

    for (unsigned int cycle=0; cycle<6; ++cycle, triangulation.refine_global(1))
      {
        setup_system ();
        assemble_and_solve ();
      };

    // After all the data is generated, write a table of results to the
    // screen:
    output_table.set_precision("|u|_1", 6);
    output_table.set_precision("error", 6);
    output_table.write_text (std::cout);
    std::cout << std::endl;
  }
}



// Finally the main function. It's structure is the same as that used in
// several of the previous examples, so probably needs no more explanation.
int main ()
{
  try
    {
      dealii::deallog.depth_console (0);
      std::cout.precision(5);

      // This is the main loop, doing the computations with mappings of linear
      // through cubic mappings. Note that since we need the object of type
      // <code>LaplaceProblem@<2@></code> only once, we do not even name it,
      // but create an unnamed such object and call the <code>run</code>
      // function of it, subsequent to which it is immediately destroyed
      // again.
      for (unsigned int mapping_degree=1; mapping_degree<=3; ++mapping_degree)
        Step11::LaplaceProblem<2>(mapping_degree).run ();
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
    };

  return 0;
}
