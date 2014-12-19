/* ---------------------------------------------------------------------
 *
 * Copyright (C) 1999 - 2014 by the deal.II authors
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
 * Authors: Wolfgang Bangerth, 1999,
 *          Guido Kanschat, 2011
 */


// @sect3{Many new include files}

// These include files are already known to you. They declare the classes
// which handle triangulations and enumeration of degrees of freedom:
#include <deal.II/grid/tria.h>
#include <deal.II/dofs/dof_handler.h>
// And this is the file in which the functions are declared that create grids:
#include <deal.II/grid/grid_generator.h>

// The next three files contain classes which are needed for loops over all
// cells and to get the information from the cell objects. The first two have
// been used before to get geometric information from cells; the last one is
// new and provides information about the degrees of freedom local to a cell:
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/dofs/dof_accessor.h>

// In this file contains the description of the Lagrange interpolation finite
// element:
#include <deal.II/fe/fe_q.h>

// And this file is needed for the creation of sparsity patterns of sparse
// matrices, as shown in previous examples:
#include <deal.II/dofs/dof_tools.h>

// The next two file are needed for assembling the matrix using quadrature on
// each cell. The classes declared in them will be explained below:
#include <deal.II/fe/fe_values.h>
#include <deal.II/base/quadrature_lib.h>

// The following three include files we need for the treatment of boundary
// values:
#include <deal.II/base/function.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>

// We're now almost to the end. The second to last group of include files is
// for the linear algebra which we employ to solve the system of equations
// arising from the finite element discretization of the Laplace equation. We
// will use vectors and full matrices for assembling the system of equations
// locally on each cell, and transfer the results into a sparse matrix. We
// will then use a Conjugate Gradient solver to solve the problem, for which
// we need a preconditioner (in this program, we use the identity
// preconditioner which does nothing, but we need to include the file anyway):
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/compressed_sparsity_pattern.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>

// Finally, this is for output to a file and to the console:
#include <deal.II/numerics/data_out.h>
#include <fstream>
#include <iostream>

// ...and this is to import the deal.II namespace into the global scope:
using namespace dealii;

// @sect3{The <code>Step3</code> class}

// Instead of the procedural programming of previous examples, we encapsulate
// everything into a class for this program. The class consists of functions
// which each perform certain aspects of a finite element program, a `main'
// function which controls what is done first and what is done next, and a
// list of member variables.

// The public part of the class is rather short: it has a constructor and a
// function `run' that is called from the outside and acts as something like
// the `main' function: it coordinates which operations of this class shall be
// run in which order. Everything else in the class, i.e. all the functions
// that actually do anything, are in the private section of the class:
class Step3
{
public:
  Step3 ();

  void run ();

  // Then there are the member functions that mostly do what their names
  // suggest and whose have been discussed in the introduction already. Since
  // they do not need to be called from outside, they are made private to this
  // class.

private:
  void make_grid ();
  void setup_system ();
  void assemble_system ();
  void solve ();
  void output_results () const;

  // And finally we have some member variables. There are variables describing
  // the triangulation and the global numbering of the degrees of freedom (we
  // will specify the exact polynomial degree of the finite element in the
  // constructor of this class)...
  Triangulation<2>     triangulation;
  FE_Q<2>              fe;
  DoFHandler<2>        dof_handler;

  // ...variables for the sparsity pattern and values of the system matrix
  // resulting from the discretization of the Laplace equation...
  SparsityPattern      sparsity_pattern;
  SparseMatrix<double> system_matrix;

  // ...and variables which will hold the right hand side and solution
  // vectors.
  Vector<double>       solution;
  Vector<double>       system_rhs;
};

// @sect4{Step3::Step3}

// Here comes the constructor. It does not much more than first to specify
// that we want bi-linear elements (denoted by the parameter to the finite
// element object, which indicates the polynomial degree), and to associate
// the dof_handler variable to the triangulation we use. (Note that the
// triangulation isn't set up with a mesh at all at the present time, but the
// DoFHandler doesn't care: it only wants to know which triangulation it will
// be associated with, and it only starts to care about an actual mesh once
// you try to distribute degree of freedom on the mesh using the
// distribute_dofs() function.) All the other member variables of the Step3
// class have a default constructor which does all we want.
Step3::Step3 ()
  :
  fe (1),
  dof_handler (triangulation)
{}


// @sect4{Step3::make_grid}

// Now, the first thing we've got to do is to generate the triangulation on
// which we would like to do our computation and number each vertex with a
// degree of freedom. We have seen these two steps in step-1 and step-2
// before, respectively.
//
// This function does the first part, creating the mesh.  We create the grid
// and refine all cells five times. Since the initial grid (which is the
// square $[-1,1] \times [-1,1]$) consists of only one cell, the final grid
// has 32 times 32 cells, for a total of 1024.
//
// Unsure that 1024 is the correct number? We can check that by outputting the
// number of cells using the <code>n_active_cells()</code> function on the
// triangulation.
void Step3::make_grid ()
{
  GridGenerator::hyper_cube (triangulation, -1, 1);
  triangulation.refine_global (5);

  std::cout << "Number of active cells: "
            << triangulation.n_active_cells()
            << std::endl;
}

// @note We call the Triangulation::n_active_cells() function, rather than
// Triangulation::n_cells(). Here, <i>active</i> means the cells that aren't
// refined any further. We stress the adjective "active" since there are more
// cells, namely the parent cells of the finest cells, their parents, etc, up
// to the one cell which made up the initial grid. Of course, on the next
// coarser level, the number of cells is one quarter that of the cells on the
// finest level, i.e. 256, then 64, 16, 4, and 1. If you called
// <code>triangulation.n_cells()</code> instead in the code above, you would
// consequently get a value of 1365 instead. On the other hand, the number of
// cells (as opposed to the number of active cells) is not typically of much
// interest, so there is no good reason to print it.


// @sect4{Step3::setup_system}

// Next we enumerate all the degrees of freedom and set up matrix and vector
// objects to hold the system data. Enumerating is done by using
// DoFHandler::distribute_dofs(), as we have seen in the step-2 example. Since
// we use the FE_Q class and have set the polynomial degree to 1 in the
// constructor, i.e. bilinear elements, this associates one degree of freedom
// with each vertex. While we're at generating output, let us also take a look
// at how many degrees of freedom are generated:
void Step3::setup_system ()
{
  dof_handler.distribute_dofs (fe);
  std::cout << "Number of degrees of freedom: "
            << dof_handler.n_dofs()
            << std::endl;
  // There should be one DoF for each vertex. Since we have a 32 times 32
  // grid, the number of DoFs should be 33 times 33, or 1089.

  // As we have seen in the previous example, we set up a sparsity pattern by
  // first creating a temporary structure, tagging those entries that might be
  // nonzero, and then copying the data over to the SparsityPattern object
  // that can then be used by the system matrix.
  CompressedSparsityPattern c_sparsity(dof_handler.n_dofs());
  DoFTools::make_sparsity_pattern (dof_handler, c_sparsity);
  sparsity_pattern.copy_from(c_sparsity);

  // Note that the SparsityPattern object does not hold the values of the
  // matrix, it only stores the places where entries are. The entries
  // themselves are stored in objects of type SparseMatrix, of which our
  // variable system_matrix is one.
  //
  // The distinction between sparsity pattern and matrix was made to allow
  // several matrices to use the same sparsity pattern. This may not seem
  // relevant here, but when you consider the size which matrices can have,
  // and that it may take some time to build the sparsity pattern, this
  // becomes important in large-scale problems if you have to store several
  // matrices in your program.
  system_matrix.reinit (sparsity_pattern);

  // The last thing to do in this function is to set the sizes of the right
  // hand side vector and the solution vector to the right values:
  solution.reinit (dof_handler.n_dofs());
  system_rhs.reinit (dof_handler.n_dofs());
}

// @sect4{Step3::assemble_system}


// The next step is to compute the entries of the matrix and right hand side
// that form the linear system from which we compute the solution. This is the
// central function of each finite element program and we have discussed the
// primary steps in the introduction already.
//
// The general approach to assemble matrices and vectors is to loop over all
// cells, and on each cell compute the contribution of that cell to the global
// matrix and right hand side by quadrature. The point to realize now is that
// we need the values of the shape functions at the locations of quadrature
// points on the real cell. However, both the finite element shape functions
// as well as the quadrature points are only defined on the reference
// cell. They are therefore of little help to us, and we will in fact hardly
// ever query information about finite element shape functions or quadrature
// points from these objects directly.
//
// Rather, what is required is a way to map this data from the reference cell
// to the real cell. Classes that can do that are derived from the Mapping
// class, though one again often does not have to deal with them directly:
// many functions in the library can take a mapping object as argument, but
// when it is omitted they simply resort to the standard bilinear Q1
// mapping. We will go this route, and not bother with it for the moment (we
// come back to this in step-10, step-11, and step-12).
//
// So what we now have is a collection of three classes to deal with: finite
// element, quadrature, and mapping objects. That's too much, so there is one
// type of class that orchestrates information exchange between these three:
// the FEValues class. If given one instance of each three of these objects
// (or two, and an implicit linear mapping), it will be able to provide you
// with information about values and gradients of shape functions at
// quadrature points on a real cell.
//
// Using all this, we will assemble the linear system for this problem in the
// following function:
void Step3::assemble_system ()
{
  // Ok, let's start: we need a quadrature formula for the evaluation of the
  // integrals on each cell. Let's take a Gauss formula with two quadrature
  // points in each direction, i.e. a total of four points since we are in
  // 2D. This quadrature formula integrates polynomials of degrees up to three
  // exactly (in 1D). It is easy to check that this is sufficient for the
  // present problem:
  QGauss<2>  quadrature_formula(2);
  // And we initialize the object which we have briefly talked about above. It
  // needs to be told which finite element we want to use, and the quadrature
  // points and their weights (jointly described by a Quadrature object). As
  // mentioned, we use the implied Q1 mapping, rather than specifying one
  // ourselves explicitly. Finally, we have to tell it what we want it to
  // compute on each cell: we need the values of the shape functions at the
  // quadrature points (for the right hand side $(\varphi_i,f)$), their
  // gradients (for the matrix entries $(\nabla \varphi_i, \nabla
  // \varphi_j)$), and also the weights of the quadrature points and the
  // determinants of the Jacobian transformations from the reference cell to
  // the real cells.
  //
  // This list of what kind of information we actually need is given as a
  // collection of flags as the third argument to the constructor of
  // FEValues. Since these values have to be recomputed, or updated, every
  // time we go to a new cell, all of these flags start with the prefix
  // <code>update_</code> and then indicate what it actually is that we want
  // updated. The flag to give if we want the values of the shape functions
  // computed is #update_values; for the gradients it is
  // #update_gradients. The determinants of the Jacobians and the quadrature
  // weights are always used together, so only the products (Jacobians times
  // weights, or short <code>JxW</code>) are computed; since we need them, we
  // have to list #update_JxW_values as well:
  FEValues<2> fe_values (fe, quadrature_formula,
                         update_values | update_gradients | update_JxW_values);
  // The advantage of this approach is that we can specify what kind of
  // information we actually need on each cell. It is easily understandable
  // that this approach can significantly speed up finite element computations,
  // compared to approaches where everything, including second derivatives,
  // normal vectors to cells, etc are computed on each cell, regardless of
  // whether they are needed or not.

  // For use further down below, we define two shortcuts for values that will
  // be used very frequently. First, an abbreviation for the number of degrees
  // of freedom on each cell (since we are in 2D and degrees of freedom are
  // associated with vertices only, this number is four, but we rather want to
  // write the definition of this variable in a way that does not preclude us
  // from later choosing a different finite element that has a different
  // number of degrees of freedom per cell, or work in a different space
  // dimension).
  //
  // Secondly, we also define an abbreviation for the number of quadrature
  // points (here that should be four). In general, it is a good idea to use
  // their symbolic names instead of hard-coding these numbers even if you know
  // them, since you may want to change the quadrature formula and/or finite
  // element at some time; the program will just work with these changes,
  // without the need to change anything in this function.
  //
  // The shortcuts, finally, are only defined to make the following loops a
  // bit more readable. You will see them in many places in larger programs,
  // and `dofs_per_cell' and `n_q_points' are more or less by convention the
  // standard names for these purposes:
  const unsigned int   dofs_per_cell = fe.dofs_per_cell;
  const unsigned int   n_q_points    = quadrature_formula.size();

  // Now, we said that we wanted to assemble the global matrix and vector
  // cell-by-cell. We could write the results directly into the global matrix,
  // but this is not very efficient since access to the elements of a sparse
  // matrix is slow. Rather, we first compute the contribution of each cell in
  // a small matrix with the degrees of freedom on the present cell, and only
  // transfer them to the global matrix when the computations are finished for
  // this cell. We do the same for the right hand side vector. So let's first
  // allocate these objects (these being local objects, all degrees of freedom
  // are coupling with all others, and we should use a full matrix object
  // rather than a sparse one for the local operations; everything will be
  // transferred to a global sparse matrix later on):
  FullMatrix<double>   cell_matrix (dofs_per_cell, dofs_per_cell);
  Vector<double>       cell_rhs (dofs_per_cell);

  // When assembling the contributions of each cell, we do this with the local
  // numbering of the degrees of freedom (i.e. the number running from zero
  // through dofs_per_cell-1). However, when we transfer the result into the
  // global matrix, we have to know the global numbers of the degrees of
  // freedom. When we query them, we need a scratch (temporary) array for
  // these numbers (see the discussion at the end of the introduction for
  // the type, types::global_dof_index, used here):
  std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);

  // Now for the loop over all cells. We have seen before how this works, so
  // the following code should be familiar including the conventional names
  // for these variables:
  DoFHandler<2>::active_cell_iterator
  cell = dof_handler.begin_active(),
  endc = dof_handler.end();
  for (; cell!=endc; ++cell)
    {
      // @note As already mentioned in step-1, there is a more convenient way
      // of writing such loops if your compiler supports the C++11
      // standard. See @ref CPP11 "the deal.II C++11 page" to see
      // how this works.
      //
      // We are now sitting on one cell, and we would like the values and
      // gradients of the shape functions be computed, as well as the
      // determinants of the Jacobian matrices of the mapping between
      // reference cell and true cell, at the quadrature points. Since all
      // these values depend on the geometry of the cell, we have to have the
      // FEValues object re-compute them on each cell:
      fe_values.reinit (cell);

      // Next, reset the local cell's contributions to global matrix and
      // global right hand side to zero, before we fill them:
      cell_matrix = 0;
      cell_rhs = 0;

      // Now it is time to start integration over the cell, which we
      // do by looping over all quadrature points, which we will
      // number by q_index.
      for (unsigned int q_index=0; q_index<n_q_points; ++q_index)
        {
          // First assemble the matrix: For the Laplace problem, the
          // matrix on each cell is the integral over the gradients of
          // shape function i and j. Since we do not integrate, but
          // rather use quadrature, this is the sum over all
          // quadrature points of the integrands times the determinant
          // of the Jacobian matrix at the quadrature point times the
          // weight of this quadrature point. You can get the gradient
          // of shape function $i$ at quadrature point with number q_index by
          // using <code>fe_values.shape_grad(i,q_index)</code>; this
          // gradient is a 2-dimensional vector (in fact it is of type
          // Tensor@<1,dim@>, with here dim=2) and the product of two
          // such vectors is the scalar product, i.e. the product of
          // the two shape_grad function calls is the dot
          // product. This is in turn multiplied by the Jacobian
          // determinant and the quadrature point weight (that one
          // gets together by the call to FEValues::JxW() ). Finally,
          // this is repeated for all shape functions $i$ and $j$:
          for (unsigned int i=0; i<dofs_per_cell; ++i)
            for (unsigned int j=0; j<dofs_per_cell; ++j)
              cell_matrix(i,j) += (fe_values.shape_grad (i, q_index) *
                                   fe_values.shape_grad (j, q_index) *
                                   fe_values.JxW (q_index));

          // We then do the same thing for the right hand side. Here,
          // the integral is over the shape function i times the right
          // hand side function, which we choose to be the function
          // with constant value one (more interesting examples will
          // be considered in the following programs).
          for (unsigned int i=0; i<dofs_per_cell; ++i)
            cell_rhs(i) += (fe_values.shape_value (i, q_index) *
                            1 *
                            fe_values.JxW (q_index));
        }
      // Now that we have the contribution of this cell, we have to transfer
      // it to the global matrix and right hand side. To this end, we first
      // have to find out which global numbers the degrees of freedom on this
      // cell have. Let's simply ask the cell for that information:
      cell->get_dof_indices (local_dof_indices);

      // Then again loop over all shape functions i and j and transfer the
      // local elements to the global matrix. The global numbers can be
      // obtained using local_dof_indices[i]:
      for (unsigned int i=0; i<dofs_per_cell; ++i)
        for (unsigned int j=0; j<dofs_per_cell; ++j)
          system_matrix.add (local_dof_indices[i],
                             local_dof_indices[j],
                             cell_matrix(i,j));

      // And again, we do the same thing for the right hand side vector.
      for (unsigned int i=0; i<dofs_per_cell; ++i)
        system_rhs(local_dof_indices[i]) += cell_rhs(i);
    }


  // Now almost everything is set up for the solution of the discrete
  // system. However, we have not yet taken care of boundary values (in fact,
  // Laplace's equation without Dirichlet boundary values is not even uniquely
  // solvable, since you can add an arbitrary constant to the discrete
  // solution). We therefore have to do something about the situation.
  //
  // For this, we first obtain a list of the degrees of freedom on the
  // boundary and the value the shape function shall have there. For
  // simplicity, we only interpolate the boundary value function, rather than
  // projecting it onto the boundary. There is a function in the library which
  // does exactly this: VectorTools::interpolate_boundary_values(). Its
  // parameters are (omitting parameters for which default values exist and
  // that we don't care about): the DoFHandler object to get the global
  // numbers of the degrees of freedom on the boundary; the component of the
  // boundary where the boundary values shall be interpolated; the boundary
  // value function itself; and the output object.
  //
  // The component of the boundary is meant as follows: in many cases, you may
  // want to impose certain boundary values only on parts of the boundary. For
  // example, you may have inflow and outflow boundaries in fluid dynamics, or
  // clamped and free parts of bodies in deformation computations of
  // bodies. Then you will want to denote these different parts of the
  // boundary by different numbers and tell the interpolate_boundary_values
  // function to only compute the boundary values on a certain part of the
  // boundary (e.g. the clamped part, or the inflow boundary). By default, all
  // boundaries have the number `0', and since we have not changed that, this
  // is still so; therefore, if we give `0' as the desired portion of the
  // boundary, this means we get the whole boundary. If you have boundaries
  // with kinds of boundaries, you have to number them differently. The
  // function call below will then only determine boundary values for parts of
  // the boundary.
  //
  // The function describing the boundary values is an object of type Function
  // or of a derived class. One of the derived classes is ZeroFunction, which
  // describes (not unexpectedly) a function which is zero everywhere. We
  // create such an object in-place and pass it to the
  // VectorTools::interpolate_boundary_values() function.
  //
  // Finally, the output object is a list of pairs of global degree of freedom
  // numbers (i.e. the number of the degrees of freedom on the boundary) and
  // their boundary values (which are zero here for all entries). This mapping
  // of DoF numbers to boundary values is done by the <code>std::map</code>
  // class.
  std::map<types::global_dof_index,double> boundary_values;
  VectorTools::interpolate_boundary_values (dof_handler,
                                            0,
                                            ZeroFunction<2>(),
                                            boundary_values);
  // Now that we got the list of boundary DoFs and their respective boundary
  // values, let's use them to modify the system of equations
  // accordingly. This is done by the following function call:
  MatrixTools::apply_boundary_values (boundary_values,
                                      system_matrix,
                                      solution,
                                      system_rhs);
}


// @sect4{Step3::solve}

// The following function simply solves the discretized equation. As the
// system is quite a large one for direct solvers such as Gauss elimination or
// LU decomposition, we use a Conjugate Gradient algorithm. You should
// remember that the number of variables here (only 1089) is a very small
// number for finite element computations, where 100.000 is a more usual
// number.  For this number of variables, direct methods are no longer usable
// and you are forced to use methods like CG.
void Step3::solve ()
{
  // First, we need to have an object that knows how to tell the CG algorithm
  // when to stop. This is done by using a SolverControl object, and as
  // stopping criterion we say: stop after a maximum of 1000 iterations (which
  // is far more than is needed for 1089 variables; see the results section to
  // find out how many were really used), and stop if the norm of the residual
  // is below $10^{-12}$. In practice, the latter criterion will be the one
  // which stops the iteration:
  SolverControl           solver_control (1000, 1e-12);
  // Then we need the solver itself. The template parameters to the SolverCG
  // class are the matrix type and the type of the vectors, but the empty
  // angle brackets indicate that we simply take the default arguments (which
  // are <code>SparseMatrix@<double@></code> and
  // <code>Vector@<double@></code>):
  SolverCG<>              solver (solver_control);

  // Now solve the system of equations. The CG solver takes a preconditioner
  // as its fourth argument. We don't feel ready to delve into this yet, so we
  // tell it to use the identity operation as preconditioner:
  solver.solve (system_matrix, solution, system_rhs,
                PreconditionIdentity());
  // Now that the solver has done its job, the solution variable contains the
  // nodal values of the solution function.
}


// @sect4{Step3::output_results}

// The last part of a typical finite element program is to output the results
// and maybe do some postprocessing (for example compute the maximal stress
// values at the boundary, or the average flux across the outflow, etc). We
// have no such postprocessing here, but we would like to write the solution
// to a file.
void Step3::output_results () const
{
  // To write the output to a file, we need an object which knows about output
  // formats and the like. This is the DataOut class, and we need an object of
  // that type:
  DataOut<2> data_out;
  // Now we have to tell it where to take the values from which it shall
  // write. We tell it which DoFHandler object to use, and the solution vector
  // (and the name by which the solution variable shall appear in the output
  // file). If we had more than one vector which we would like to look at in
  // the output (for example right hand sides, errors per cell, etc) we would
  // add them as well:
  data_out.attach_dof_handler (dof_handler);
  data_out.add_data_vector (solution, "solution");
  // After the DataOut object knows which data it is to work on, we have to
  // tell it to process them into something the back ends can handle. The
  // reason is that we have separated the frontend (which knows about how to
  // treat DoFHandler objects and data vectors) from the back end (which knows
  // many different output formats) and use an intermediate data format to
  // transfer data from the front- to the backend. The data is transformed
  // into this intermediate format by the following function:
  data_out.build_patches ();

  // Now we have everything in place for the actual output. Just open a file
  // and write the data into it, using GNUPLOT format (there are other
  // functions which write their data in postscript, AVS, GMV, or some other
  // format):
  std::ofstream output ("solution.gpl");
  data_out.write_gnuplot (output);
}


// @sect4{Step3::run}

// Finally, the last function of this class is the main function which calls
// all the other functions of the <code>Step3</code> class. The order in which
// this is done resembles the order in which most finite element programs
// work. Since the names are mostly self-explanatory, there is not much to
// comment about:
void Step3::run ()
{
  make_grid ();
  setup_system ();
  assemble_system ();
  solve ();
  output_results ();
}


// @sect3{The <code>main</code> function}

// This is the main function of the program. Since the concept of a main
// function is mostly a remnant from the pre-object era in C/C++ programming,
// it often does not much more than creating an object of the top-level class
// and calling its principle function. This is what is done here as well:
int main ()
{
  Step3 laplace_problem;
  laplace_problem.run ();

  return 0;
}
